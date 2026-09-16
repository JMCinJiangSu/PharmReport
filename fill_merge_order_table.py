#!/usr/bin/env python3
"""用 Jinja2/docxtpl 将 JSON 的 merge_order 写入原始 Word 模板。"""

from __future__ import annotations

import argparse
import json
import tempfile
import zipfile
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET

import jinja2
from docxtpl import DocxTemplate

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
META_FIELDS = {"gene_symbol", "hgvs_c", "hgvs_p", "group_key", "count"}
A4_WIDTH_TWIPS, A4_HEIGHT_TWIPS = "11906", "16838"
ET.register_namespace("w", W)


def w(tag: str) -> str:
    return f"{{{W}}}{tag}"


def sample_columns(records: list[dict[str, Any]]) -> list[str]:
    """按字段名确定频率列，忽略仅作冗余展示的 freq_N 字段。"""
    columns = [key for key in records[0] if key not in META_FIELDS and not key.startswith("freq_")]
    if not columns:
        raise ValueError("var.merge_order 中未找到样本频率列。")
    return columns


def validate_column_count(records: list[dict[str, Any]], columns: list[str]) -> None:
    """count 决定动态生成的频率列数，且必须与样本字段数相同。"""
    counts = {record.get("count") for record in records}
    if len(counts) != 1 or not isinstance(next(iter(counts)), int):
        raise ValueError("merge_order 的 count 必须是所有记录一致的整数。")
    count = next(iter(counts))
    if count < 1 or count != len(columns):
        raise ValueError(f"count 为 {count}，但检测到 {len(columns)} 个样本频率列。")


def build_context(payload: dict[str, Any]) -> dict[str, Any]:
    """构造 docxtpl 上下文；每条记录的 VAF 列数等于 count。"""
    records = payload.get("var", {}).get("merge_order", [])
    if not isinstance(records, list) or not records or not all(isinstance(row, dict) for row in records):
        raise ValueError("JSON 的 var.merge_order 必须是非空对象数组。")

    columns = sample_columns(records)
    validate_column_count(records, columns)
    rows = []
    for record in records:
        rows.append({
            "gene_symbol": str(record.get("gene_symbol", "")),
            "hgvs_c": str(record.get("hgvs_c", "")),
            "hgvs_p": str(record.get("hgvs_p", "")),
            # 同一 sample 字段恒定写入同一列，不依赖 freq_1/freq_2 的位置。
            "vaf_columns": [str(record.get(column) or "-") for column in columns],
        })
    return {"table": {"sample_columns": columns, "rows": rows}}


def set_cell_text(cell: ET.Element, value: str) -> None:
    """用单个 Jinja 标签替换单元格显示文本，保留该单元格原有样式。"""
    text_nodes = cell.findall(f".//{w('t')}")
    if not text_nodes:
        paragraph = cell.find(w("p"))
        if paragraph is None:
            paragraph = ET.SubElement(cell, w("p"))
        run = ET.SubElement(paragraph, w("r"))
        text_nodes = [ET.SubElement(run, w("t"))]
    text_nodes[0].text = value
    for node in text_nodes[1:]:
        node.text = ""


def remove_vertical_merge(cell: ET.Element) -> None:
    cell_properties = cell.find(w("tcPr"))
    if cell_properties is not None:
        merge = cell_properties.find(w("vMerge"))
        if merge is not None:
            cell_properties.remove(merge)


def make_jinja_template(source_path: Path, template_path: Path) -> None:
    """在临时副本中加入 docxtpl 的动态列标签，不改动原始 col_fill.docx。"""
    with zipfile.ZipFile(source_path) as source:
        root = ET.fromstring(source.read("word/document.xml"))
        table = root.find(f".//{w('tbl')}")
        if table is None:
            raise ValueError("Word 模板中未找到表格。")
        table_rows = table.findall(w("tr"))
        if len(table_rows) < 5:
            raise ValueError("Word 模板表格至少需要标题行、表头行和 3 个数据模板行。")

        # 删除原合并标题行；将第二行改为使用 {%tc %} 的动态频率列表头。
        title_row, header_row, for_row, body_row, end_row = table_rows[:5]
        header_cells = header_row.findall(w("tc"))
        body_cells = body_row.findall(w("tc"))
        if len(header_cells) < 6 or len(body_cells) < 6:
            raise ValueError("Word 模板必须包含基因、变异和至少 4 个示例频率列。")
        table.remove(title_row)
        for cell in header_cells[:2]:
            remove_vertical_merge(cell)
        for cell in header_cells[5:]:
            header_row.remove(cell)
        for cell, value in zip(header_row.findall(w("tc")), [
            "基因", "变异", "{%tc for sample in table.sample_columns %}",
            "基因变异丰度（VAF）\n{{ sample }}", "{%tc endfor %}",
        ]):
            set_cell_text(cell, value)

        for cell in body_cells[5:]:
            body_row.remove(cell)
        for cell, value in zip(body_row.findall(w("tc")), [
            "{{ item.gene_symbol }}", "{{ item.hgvs_c }}（{{ item.hgvs_p }}）",
            "{%tc for vaf in item.vaf_columns %}", "{{ vaf }}", "{%tc endfor %}",
        ]):
            set_cell_text(cell, value)
        set_cell_text(for_row.findall(w("tc"))[0], "{%tr for item in table.rows %}")
        set_cell_text(end_row.findall(w("tc"))[0], "{%tr endfor %}")

        # 删除原有示例数据行；保留 Jinja 表头、tr 循环起止行及其正文行。
        retained_rows = {header_row, for_row, body_row, end_row}
        for row in list(table.findall(w("tr"))):
            if row not in retained_rows:
                table.remove(row)
        table_properties = table.find(w("tblPr"))
        if table_properties is not None:
            layout = table_properties.find(w("tblLayout"))
            if layout is None:
                layout = ET.SubElement(table_properties, w("tblLayout"))
            layout.set(w("type"), "autofit")
        section = root.find(f".//{w('sectPr')}")
        if section is not None:
            page_size = section.find(w("pgSz"))
            if page_size is not None:
                page_size.set(w("w"), A4_WIDTH_TWIPS)
                page_size.set(w("h"), A4_HEIGHT_TWIPS)
                page_size.attrib.pop(w("orient"), None)

        with zipfile.ZipFile(template_path, "w", zipfile.ZIP_DEFLATED) as target:
            for item in source.infolist():
                content = ET.tostring(root, encoding="utf-8", xml_declaration=True) if item.filename == "word/document.xml" else source.read(item.filename)
                target.writestr(item, content)


def fill_word(json_path: Path, template_path: Path, output_path: Path) -> None:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    context = build_context(payload)
    with tempfile.TemporaryDirectory() as temporary_directory:
        jinja_template = Path(temporary_directory) / "merge_order_template.docx"
        make_jinja_template(template_path, jinja_template)
        document = DocxTemplate(jinja_template)
        document.render(context, jinja_env=jinja2.Environment(autoescape=False))
        document.save(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="用 Jinja2 Word 模板填充动态 merge_order 频率表")
    parser.add_argument("--json", type=Path, default=Path("decb550cb19a11f1814b0242c0a87008_to_word.json"))
    parser.add_argument("--template", type=Path, default=Path("col_fill.docx"))
    parser.add_argument("--output", type=Path, default=Path("col_fill_filled.docx"))
    args = parser.parse_args()
    fill_word(args.json, args.template, args.output)


if __name__ == "__main__":
    main()
