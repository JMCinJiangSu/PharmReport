#-*- coding:gbk -*-
import re
from libs.getConfig import clinicalNumStran, functionNumStran, typeStran, typeStran_from_inter
from libs.getEvi import varRegimen
import copy
from libs.specialRequest import varInfo_XAJDY, varInter_FJZL, varInfo_SYX
from libs.hgvsp_To_hgvsp_abbr import splitAA
from libs.specialRequest import varInfo_FDZS, varInfo_FJZL
from libs.rule import S_function
from libs.rule import decimal_float, decimal_percen
from libs.getConfig import get_gene_class
import datetime
from decimal import Decimal

'''
Discription
	
	����snvindel��ʽ�� 

'''

def process_snvindel(jsonDict, config):
	snvindel = copy.deepcopy(jsonDict["snvindel"])
	#print  ("######### 3.0.1 step0", datetime.datetime.now())
	function = get_gene_class(config)["function"]
	clinical = get_gene_class(config)["clinical"]
	clinicalNumStran_dict = clinicalNumStran(config)
	functionNumStran_dict = functionNumStran(config)
	for var in snvindel:
		# �²���/�°��Զ��У�����ϵƥ���²��ԣ���ϸ��ƥ���°��ԣ�
		# �²���/�°��Խ���һ������ϵ/��ϸ����ƥ�������ݵ��Ǹ���
		# �²���/�°��Ծ��ޣ���ϵ/��ϸ����Ĭ��Ϊ3��
		if var["clinical_significance"] != "-" or var["function_classification"] != "-":
			var["clinic_num_g"] = clinicalNumStran_dict.get(var["clinical_significance"], 3) if var["clinical_significance"] != "-" else \
								  clinicalNumStran_dict.get(var["function_classification"], 3) 
			var["clinic_num_s"] = functionNumStran_dict.get(var["function_classification"], 3) if var["function_classification"] != "-" else \
								  functionNumStran_dict.get(var["clinical_significance"], 3)
		else:
			var["clinic_num_g"] = 3
			var["clinic_num_s"] = 3
		# AD3101 ����һ���ֶΣ�ֻ���²��Խ���ķ����²��Եȼ�������������°��Խ���ȼ����,2025��3��6��
		var['clinic_ad3101_g'] = var['clinic_num_g'] if var["clinical_significance"] != "-" and var['function_classification'] == '-' else 3
		var['clinic_ad3101_s'] = var['clinic_num_s'] if var['function_classification'] != '-' else 3

		# freq ���ݰٷ�����ʽ����ת��ΪС����ʽ-����v4-��쿷ң�2023.11.09
		# freq �����ַ�����ʽ����ת��ΪС����ʽ-���γ���2024.03.20
		freq_key = ["freq", "freq_sc", "freq_rc", "freq_ctrl", "freq_case"]
		for key in freq_key:
			if key in var.keys():
				var[key] = float(var[key].replace("%", ""))/100 if re.search("%", str(var[key])) else float(var[key]) if isinstance(var[key], str) and var[key] else var[key]
		# ���ݽ���-2023.11.09

		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["clinic_num_s"], var["top_level"] = S_function(var)
		var["freq_str"] = decimal_percen(var["freq"]) if var["freq"] else ""
		var["freq_sc_str"] = decimal_percen(var["freq_sc"]) if var["freq_sc"] else ""
		# ����freq_rcδ���ص����
		var["freq_rc"] = var["freq_rc"] if "freq_rc" in var.keys() and var["freq_rc"] else var["freq_ctrl"] if "freq_ctrl" in var.keys() and var["freq_ctrl"] else ""
		var["freq_rc_str"] = decimal_percen(var["freq_rc"]) if "freq_rc" in var.keys() and var["freq_rc"] else decimal_percen(var["freq_ctrl"]) if var["freq_ctrl"] else ""
		var["freq_ctrl_str"] = decimal_percen(var["freq_ctrl"]) if var["freq_ctrl"] else ""
		var["freq_ss_str"] = decimal_percen(var["freq_ss"]) if var["freq_ss"] else ""
		var["freq_case_str"] = decimal_percen(var["freq_case"]) if var["freq_case"] else ""
		
		# ���⴦��
		# 1. hgvs_p������-Ĭ�ϸ�ʽ�� �������ݣ���ֹ��������δ��hgvs_p��
		var["hgvs_p"] = var["hgvs_p"] if var["hgvs_p"] else "p.?"
		if not re.search("=", var["hgvs_p"]) and var["hgvs_p"] != "p.?" and not re.search("\(", var["hgvs_p"]):
			var["hgvs_p"] = var["hgvs_p"].replace("p.", "p.(")+")"
		# 2. hgvs_p_abbr������-Ĭ�ϸ�ʽ
		var["hgvs_p_abbr"] = var["hgvs_p_abbr"] if var["hgvs_p_abbr"] else ""
		if var["hgvs_p_abbr"] and not re.search("=", var["hgvs_p_abbr"]) and var["hgvs_p_abbr"] != "p.?" and not re.search("\(", var["hgvs_p_abbr"]):
			var["hgvs_p_abbr"] = var["hgvs_p_abbr"].replace("p.", "p.(")+")"
		# 3. hgvs_p����ͬ��ͻ�䣬���಻������-��������
		var["hgvs_p_ZJZL"] = var["hgvs_p"].replace("(", "").replace(")", "") if not re.search("=", var["hgvs_p"]) else var["hgvs_p"]
		# 4. hgvs_pȫ����������-������CP40������һ
		var["hgvs_p_FJFY"] = var["hgvs_p"].replace("(", "").replace(")", "")
		# 5. ����С�汾�ŵ�ת¼��		
		var["transcript_primary_simple"] = re.split("\.", str(var["transcript_primary"]))[0] if var["transcript_primary"] else ""
		# 6. ��������ת��Ϊ���ģ���������ʱ�������Ϳ��ܻ�Ϊ�գ��ᵼ�³��򱨴�������¼���
		#var["type_"] = var["type"] if var["type"] else ""
		var["type_cn"], var["type"] = typeStran_from_inter(var)
		var["type_cn"] = var["type_cn"] if var["type_cn"] else ""
		# ��ʱ��һ��AltDepth����֣��һ���ȼ�Ԫ���غ�����ϵͳ�Դ��ģ�2023.07.03
		if var["depth"] and var["freq"]:
			var["altdepth"] = str(Decimal(float(var["depth"]) * float(var["freq"])).quantize(Decimal("0"), rounding="ROUND_HALF_UP"))
		else:
			var["altdepth"] = "��Ȼ�Ƶ��Ϊ�գ����˹���д��"

		# ������������
		# 1. ����������������
		var["varInter_FJZL"] = varInter_FJZL(var["hgvs_p"], var["type"], var["gene_region"], var["hgvs_c"])
		# 2. ��������һ����X��������/�ں���XXXͻ��
		var["varInfo_XAJDY"] = varInfo_XAJDY(var["gene_region"], var["type_cn"]) if var["type_cn"] != "--" else varInfo_XAJDY(var["gene_region"], "ͻ��")
		# 3. �����ɣ�XX��������/�ں���XX hgvs_p/hgvs_c XXͻ��
		var["varInfo_SYX"] = varInfo_SYX(var["gene_region"], var["type_cn"], var["type"], var["hgvs_c"], var["hgvs_p_ZJZL"])
		# 4. ������ĸ������Ϊ�գ���ʹ�ñ���ű�ת���Ľ��-2022.09.07
		var["hgvs_p_abbr"] = splitAA(var["hgvs_p"]) if not var["hgvs_p_abbr"] else var["hgvs_p_abbr"]
		# 5. ������ɽ����������������
		var["var_info_forFDZS"] = varInfo_FDZS(var["gene_symbol"], var["gene_region"], var["variant_desc_cn"], config)
		# 6. �㽭��������������ͳһ�¸�ʽ����Ҫ�ӡ������γɹ������˻�ʧ��ĵ��ס������Լ��ڱ������Ϊ�˱���֪ʶ���и�����Ϣ��ʽ��ͳһ�����пո���ɱ���չʾ���ң����Ԥ����һ��
		var["variant_desc_cn_ZJZL"] = var["variant_desc_cn"].strip()[0:-1] if var["variant_desc_cn"] and var["variant_desc_cn"].strip()[-1] == "��" else var["variant_desc_cn"]
		# 7. �������������ر���Ƶ�������Ϣ����Դ���ñ������Ʒ�������
		var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
		# 8. ����ҽԺ��XX��������/�ں���XX hgvs_p/hgvs_c XXͻ�� ͬ�����ɣ�����hgvs_p��Ҫ��������-2023.08.03
		var["varInfo_SZYY"] = varInfo_SYX(var["gene_region"], var["type_cn"], var["type"], var["hgvs_c"], var["hgvs_p"])
		# 9. XW0413 EGFR���ӱ�ǩ Ex19del Ex20ins  ���γ�
		# EGFR_tag����Ҫ���ż� "p." 2024��6��17�� ���γ�
		#var['EGFR_tag'] = 'Ex19 del' if 'EGFR Exon19 del' in var['var_category_names'] else 'Ex20 ins' if 'EGFR Exon20 ins' in var['var_category_names'] \
		#	else var['hgvs_p'].replace('p.', '') if not re.search('=', var['hgvs_p']) and var['hgvs_p'] != 'p.?' else var['hgvs_p']
		# json���ظĳ�null�ˣ��Ӹ��ж� 2024��10��29�� jmc
		if var['var_category_names']:
			var['EGFR_tag'] = 'Ex19 del' if 'EGFR Exon19 del' in var['var_category_names'] else 'Ex20 ins' if 'EGFR Exon20 ins' in var['var_category_names'] else var['hgvs_p_ZJZL'].replace('p.', '') if var['hgvs_p_ZJZL'] != 'p.?' else var['hgvs_p_ZJZL']
		#var['EGFR_tag'] = 'Ex19 del' if 'EGFR Exon19 del' in var['var_category_names'] else 'Ex20 ins' if 'EGFR Exon20 ins' in var['var_category_names'] \
		#	else var['hgvs_p_ZJZL'].replace('p.', '') if var['hgvs_p_ZJZL'] != 'p.?' else var['hgvs_p_ZJZL']
	#print  ("######### 3.0.2 step2", datetime.datetime.now())
	# ��Ƶ��������
	snvindel = sorted(snvindel, key=lambda i:i["freq"], reverse=True)
	#print  ("######### 3.0.3 step3", datetime.datetime.now())
	return snvindel