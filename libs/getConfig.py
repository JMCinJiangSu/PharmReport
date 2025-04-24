#-*- coding:gbk -*-
import os
import json
import xlrd
import re

'''
Discription
	
	��ȡconfig�ļ����е������ļ�������Ϊ��������������õĸ�ʽ��	

'''

# 1. ��ȡjson�����ļ������Ϣ
def getconfigfile(config):
	return json.loads(open(os.path.join(config, "config.json"), encoding="utf-8").read())

def clinicalNumStran(config): return getconfigfile(config)["clinical_num_stran"]
def functionNumStran(config): return getconfigfile(config)["function_num_stran"]
def senseStran(config): return getconfigfile(config)["sense_stran"]
def typeStran(config): return getconfigfile(config)["type_stran"]
def evidenceTypeStran(config): return getconfigfile(config)["evidence_type_stran"]
def CP40Gene(config): return getconfigfile(config)["cp40_gene"]
def CP40_RECOMGene(config): return getconfigfile(config)["rpt_guideline_recom"]
def CP40_SplitGene(config): return getconfigfile(config)["cp40_split_gene"]
def PAN116_RECOMGene(config): return getconfigfile(config)["rpt_guideline_recom_116"]
def PAN116_other_gene(config): return getconfigfile(config)["gene_116_other"]

# 2. ��ȡxlsx�����ļ������Ϣ
def getconfigxlsx(config):
	file_path = os.path.join(config, "config.xlsx")
	xls = xlrd.open_workbook(file_path)
	# 2.1. ������ɽ���ñ�������ؼ�����������ܡ��Ŵ����պͷ���ֵ
	fdzs = xls.sheet_by_name("rpt_HRR_FDZS-v4")
	fdzs_dict = {}
	for num in range(1, fdzs.nrows):
		rows = fdzs.row_values(num)
		fdzs_dict.setdefault(rows[0], {})
		fdzs_dict[rows[0]] = {
				"disease" : rows[1],
				"gene_info" : rows[3],
				"disease_info" : rows[4],
				"risk_info" : rows[5]
			}
	# 2.2. ��������LC10���ñ��������򡢰��ֺͱ��췢����
	fjzl = xls.sheet_by_name("gene_freq_inter-v3")
	fjzl_dict = {}
	for num in range(1, fjzl.nrows):
		rows = fjzl.row_values(num)
		fjzl_dict.setdefault((rows[0], rows[1]), "")
		fjzl_dict[(rows[0], rows[1])] = rows[2]
	# 2.3. io, �����㽭����
	io = xls.sheet_by_name("rpt_io")
	io_data = {}
	for num in range(1, io.nrows):
		rows = io.row_values(num)
		io_data[rows[1]] = rows[2]

	return fdzs_dict, fjzl_dict, io_data

# ��ʱ��һ�����ã��ȼ�Ԫ�������غ��ٽ��ã�������ɸѡ���찴�°��Ի��ǰ��²��Խ��з���-2023.05.23
def get_gene_class(config):
	database_for_gene_class = os.path.join(config, "gene_class.xlsx")
	rpt_xls = xlrd.open_workbook(database_for_gene_class)
	rpt_sheet_name = rpt_xls.sheet_names()
	result = {}
	for each in rpt_sheet_name:
		sheet = rpt_xls.sheet_by_name(each)
		for num in range(1, sheet.nrows):
			rows = sheet.row_values(num)
			if rows[1] not in result.keys():
				result.setdefault(rows[1], [])
			result[rows[1]].append(rows[0])
	return result

# ��������ת�����ӱ�����������ȡ
# 1. �ñ���Ϊ������ͻ�䡱�����¡���
# 2. �ñ���Ϊ���ں���ͻ�䡱��
def typeStran_from_inter(var):
	# Ϊ���±�������ʱ��������type�ֶ����ݣ����ص�type�����д���
	type_dict = {
		"�ں�����ͻ��" : "Intronic",
		"5'UTR��ͻ��" : "5'UTR",
		"3'UTR��ͻ��" : "3'UTR",
		"��������ͻ��" : "FlankingRegion5",
		"��������ͻ��" : "FlankingRegion3"
	}

	type_from_inter =  re.split("��", var["variant_desc_cn"][0:-1])[0].replace("�ñ���Ϊ","")
	type_cn = type_from_inter if type_from_inter not in type_dict.keys() else "--"
	type_en = var["type"] if type_cn != "--" else type_dict.get(type_from_inter, var["type"])
	return type_cn, type_en

# ��ȡxlsx�����ļ����shrj_cp_gene_info���Ϻ��ʼ�CP40������ܣ�-2023.09.18
def get_shrj_geneinfo(config):
	file_path = os.path.join(config, "config.xlsx")
	xls = xlrd.open_workbook(file_path)
	gene_info = xls.sheet_by_name("shrj_cp_gene_info")
	gene_info_dic = {}
	for num in range(1, gene_info.nrows):
		rows = gene_info.row_values(num)
		gene_info_dic[rows[0]] = rows[1]

	return gene_info_dic