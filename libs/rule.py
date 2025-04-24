#-*- coding:gbk -*-
import re
import copy
from decimal import Decimal

'''
	ͳһ����ģ���еĸ��ֹ���
	��ϸ����
	# I�ࣺ����/�������/Ԥ��֤����ߵȼ�ΪA/B
	# II�ࣺ����/�������/Ԥ��֤����ߵȼ�ΪC/D
	# ����������չ��أ�����/�������/Ԥ��δƥ�䵽֤�ݣ��ұ���Ϊ�²�/�����²���
	# III�ࣺ����/�������/Ԥ��δƥ�䵽֤�ݣ���Ϊ���岻��λ��
'''

## ����ʹ�õ��ж������������
# 1. ���ر�����Դ����ϵ/��ϸ����δ֪��Դ�İ���ϸ������
def origin(var):
	if var["var_origin"] == "germline":
		return "G"
	else:
		return "S"

# 2. �ж��Ƿ�ƥ�䵽����/�������/Ԥ���֤��
def judgeRegimen(var):
	if var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()):
		return 1
	else:
		return 0

# 3. ��������Ƿ���Ŀ������б��У���������Ŀ���ת���б��ʽ���д���
def matchGene(var, genelist):
	return set(re.split(",", var["gene_symbol"])) & set(genelist)

# 4. �ж������Ƿ�Ϊ��ϵ[ָ���ȼ�]����ϵ[ָ���ȼ�]
def judge_var(var, s_list, g_list):
	if (var["var_origin"] == "germline" and var["clinic_num_g"] in g_list) or (var["var_origin"] != "germline" and var["clinic_num_s"] in s_list):
		return 1
	else:
		return 0

# 5. ���ر���������Snvindel����hgvs_p��hgvs_c
def get_varsimpleinfo(var):
	if var["bio_category"] == "Snvindel":
		return var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
	elif var["bio_category"] == "Cnv":
		return "����"
	elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
		return "{0}-{1}�ں�".format(var["five_prime_gene"], var["three_prime_gene"])
## ��ʹ���ж�����-����

# ��ϸ�������ΪI/II/����������չ���/III��-���ֻ��� �����ڴ󲿷�ģ��
def s_var_rule(var_data):
	result = {}
	result["level_I"] = [var for var in var_data if var["clinic_num_s"] in [5] and origin(var)=="S" and judgeRegimen(var)]
	result["level_II"] = [var for var in var_data if var["clinic_num_s"] in [4] and origin(var)=="S" and judgeRegimen(var)]
	result["level_onco_nodrug"] = [var for var in var_data if var["clinic_num_s"] in [5, 4] and origin(var)=="S" and not judgeRegimen(var)]
	result["level_III"] = [var for var in var_data if var["clinic_num_s"] == 3 and origin(var)=="S"]
	# IV�ࣺ ������Ҫչʾ�����б��HRR-20220915����
	result["level_IV"] =  sorted([var for var in var_data if var["clinic_num_s"] in [1,2] and origin(var)=="S"], key=lambda i:i["clinic_num_s"], reverse=True)
	# III�಻����ͬ��ͻ��-���ø���ʡ��CP40-20221012 "Synonymous_Substitution"
	result["level_III_without_Syn"] = list(filter(lambda i : "type" in i.keys() and i["type"] != "Synonymous_Substitution", result["level_III"]))
	return result

# ��ϸ�������ΪI/II/����������չ���/III��-�ֵ������� ������PTM��BPTM�Ͳ��ֶ���ģ��
def s_var_rule_gene(var_data, gene):
	result = {}
	result["level_I"] = [var for var in var_data if matchGene(var, [gene]) and var["clinic_num_s"] == 5 and origin(var)=="S" and judgeRegimen(var)]
	result["level_II"] = [var for var in var_data if matchGene(var, [gene]) and var["clinic_num_s"] == 4 and origin(var)=="S" and judgeRegimen(var)]
	result["level_onco_nodrug"] = [var for var in var_data if matchGene(var, [gene]) and var["clinic_num_s"] in [5, 4] and origin(var)=="S" and not judgeRegimen(var)]
	result["level_III"] = [var for var in var_data if matchGene(var, [gene]) and var["clinic_num_s"] == 3 and origin(var)=="S"]
	return result

# ��ϸ�������ΪI/II/����������չ���/III��-�ֶ������ �������Ϻ�ʮԺHRR�Ͳ��ֶ���ģ��
def s_var_rule_genelist(var_data, genelist):
	result = {}
	result["level_I"] = [var for var in var_data if matchGene(var, genelist) and var["clinic_num_s"] == 5 and origin(var) == "S" and judgeRegimen(var)]
	result["level_II"] = [var for var in var_data if matchGene(var, genelist) and var["clinic_num_s"] == 4 and origin(var) == "S" and judgeRegimen(var)]
	result["level_onco_nodrug"] = [var for var in var_data if matchGene(var, genelist) and var["clinic_num_s"] in [5, 4] and origin(var)=="S" and not judgeRegimen(var)]
	result["level_III"] = [var for var in var_data if matchGene(var, genelist) and var["clinic_num_s"] == 3 and origin(var)=="S"]
	return result

# ��ϵ��Ϊ5/4/3/2/1��
def g_var_rule(var_data):
	result = {}
	for i in range(1, 6):
		result["level_"+str(i)] = [var for var in var_data if var["clinic_num_g"] == i and origin(var) == "G"]
	return result

# ��ϵ��Ϊ5/4/3/2/1�࣬��ƥ������б�
def g_var_rule_genelist(var_data, genelist):
	result = {}
	for i in range(1, 6):
		result["level_"+str(i)] = [var for var in var_data if var["clinic_num_g"] == i and origin(var) == "G" and matchGene(var, genelist)]
	return result

# ��ϵ�����������/�������/Ԥ��֤����ȡI/II��
def g_var_regimen_rule(var_data):
	result = {}
	result["regimen_level_I"] = [var for var in var_data if var["clinic_num_s"] == 5 and origin(var)=="G" and judgeRegimen(var)]
	result["regimen_level_II"] = [var for var in var_data if var["clinic_num_s"] == 4 and origin(var)=="G" and judgeRegimen(var)]
	return result

# ��ϵ+��ϸ�������������/�������/Ԥ��֤����ȡI/II��
def var_regimen_rule(var_data):
	result = {}
	result["level_I"] = [var for var in var_data if var["clinic_num_s"] == 5 and judgeRegimen(var)]
	result["level_II"] = [var for var in var_data if var["clinic_num_s"] == 4 and judgeRegimen(var)]
	return result

# �㽭������ϵ+��ϸ������ҩ����/Ԥ��/������ϱ������������н�չʾ��ҩ������ҩ����Ԥ��/������ϵ�д�����ް�����ҩ��ʾ����ģ������н���ʵ�֣�������ټ��ֶΣ�
def process_result_regimen_ZJZL(var_list_zjzl):
	for var in var_list_zjzl:
		if "Predictive" not in var["evi_sum"]["evi_split"]:
			var["evi_sum"]["evi_split"]["Predictive"] = {"note" : "without_regimen"}
	return var_list_zjzl

# PTM/BPTM��Ҫչʾ����������������ΪTP53��POLE��BRCA������֤�ݣ�����������������չ��ر��죬�еĻ�����������߸ġ�
def var_bptm_rule(var_data, gene):
	result = {}
	# ��ʽ1��I/II��III��
	result[gene+"_level12"] = [var for var in var_data if var["clinic_num_s"] in [5, 4] and origin(var)=="S" and judgeRegimen(var) and matchGene(var, [gene])]
	result[gene+"_level3"] = [var for var in var_data if var["clinic_num_s"] == 3 and origin(var)=="S" and matchGene(var, [gene])]
	# ����������չ��ر���ŵ�III����-2023.09.12
	result[gene+"_level3"] = [var for var in var_data if (var["clinic_num_s"] == 3 or (var["clinic_num_s"] in [5, 4]) and not judgeRegimen(var)) and origin(var)=="S" and matchGene(var, [gene])]
	# 2023.09.12-�������

	# ��ʽ2�����ʽ1�����������δ������죬��Ҫ���
	result[gene+"_level12_withECtype"] = copy.deepcopy(result[gene+"_level12"]) if result[gene+"_level12"] else [{"gene_symbol" : gene, "result" : "nofound"}]
	result[gene+"_level3_withECtype"] = copy.deepcopy(result[gene+"_level3"]) if result[gene+"_level3"] else [{"gene_symbol" : gene, "result" : "nofound"}]
	# ��ʽ3��I/II/III�����һ��
	level_data = s_var_rule_gene(var_data, gene)
	result[gene+"_withEC_type"] = level_data["level_I"] + level_data["level_II"] + level_data["level_onco_nodrug"] + level_data["level_III"]
	result[gene+"_withEC_type"] = result[gene+"_withEC_type"] if result[gene+"_withEC_type"] else [{"gene_symbol" : gene, "result" : "nofound"}]
	return result

# ��ϵ-�����ۺ����������������Ҫ�ֿ�չʾ��EPCAM��MLH1��MSH2��MSH6��PMS2����չʾ���������3��4��5�����
def var_lyn5_rule(var_data, gene_list):
	result = {}
	for gene in gene_list:
		if gene not in result.keys():
			result.setdefault(gene, [])
		result[gene] = [var for var in var_data if var["clinic_num_g"] in [5, 4, 3] and origin(var)=="G" and matchGene(var, [gene])]
	return result

# ��ȡδ��⵽I/II����ϸ������Ļ���
def nofoundPath_genelist(var_data, gene_list):
	return sorted(list(set(gene_list) - set([var["gene_symbol"] for var in var_data if var["clinic_num_s"] in [4,5] and origin(var)=="S" and judgeRegimen(var)])))

# ��ȡδ��⵽��ϸ������Ļ���(����I/II/III/����������չ���)
def nofound_genelist(var_data, gene_list):
	return sorted(list(set(gene_list) - set([var["gene_symbol"] for var in var_data if var["clinic_num_s"] in [3,4,5] and origin(var) == "S"])))

# �²�/�����²�������ҩ�ı��죬�ڽ�����ܱ��й�ΪIII��
def S_level(var):
	s_level = var["clinic_num_s"] if judgeRegimen(var) else 3
	return s_level

# ��ϸ���ȼ�����ߵȼ�ABΪI�࣬CDΪII�࣬��֤���򰴽���ķ���
# ��ȡ֤����ߵȼ�
def S_function(var):
	regimen_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive", "Prognostic", "Diagnostic"]] \
					 if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 else []

	clinic_num_s =  5 if set(["A", "B"]) & set(regimen_level) else \
					4 if set(["C", "D"]) & set(regimen_level) else \
					var["clinic_num_s"]
					
	top_level = "A" if "A" in regimen_level else \
				"B" if "B" in regimen_level else \
				"C" if "C" in regimen_level else \
				"D" if "D" in regimen_level else \
				"N"
	return clinic_num_s, top_level

def decimal_float(a):
	'''
	�������룬��������
	'''
	#return str(Decimal(str(a)).quantize(Decimal("0.01"), rounding="ROUND_HALF_UP"))

	# ����Ϊ<5��ȥ��>=5��1����excel����һ�£����п��ܺ�ϵͳչʾ�Ĳ�ͬ-2023.07.04
	return "{:.2f}".format(float(int(float(a) * 100 + 0.5) / 100))

def decimal_percen(a):
	'''
	�������룬����ٷ���
	'''
	#return str(Decimal(str(float(a)*100)).quantize(Decimal("0.01"), rounding="ROUND_HALF_UP"))+"%"

	# ����Ϊ<5��ȥ��>=5��1����excel����һ�£����п��ܺ�ϵͳչʾ�Ĳ�ͬ-2023.07.04
	return "{:.2%}".format(float(int(float(a) * 10000 + 0.5) / 10000))

# ʶ���Ƿ�Ϊ��ֵ
def is_number(i):
	try:
		float(i)
		return True
	except:
		pass 
	if i.isnumeric():
		return True
	return False