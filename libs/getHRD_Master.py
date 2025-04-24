#-*- coding:gbk -*-
import copy
import re
from libs import listResultToDict
from libs.getEvi import varRegimen

'''
Discription
	
	�ýű�������ȡHRD�����Master��HRD C�Ĺ���ͬ��HRD C�ķ�bin�ļ����ˣ�������ȷ�libs�ļ��а�
	��ȡ���ݣ�
	1. GSS��ֵ
	2. �ж��ʿ��Ƿ�ϸ��ⲿ��ģ�����ж� �ϸ��׼��cellularity >= 0.3 and baf_noise <= 0.055 and depth_noise <= 0.35��
	3. BRCA�����
	4. HRRͨ·�������� 

'''
# ���ر�����
def getVar(a):
	if a["bio_category"] == "Snvindel":
		return "{0} {1}".format(a["gene_symbol"], a["hgvs_p"]) if a["hgvs_p"] != "p.?" else "{0} {1}".format(a["gene_symbol"], a["hgvs_c"])
	elif a["bio_category"] == "Cnv":
		return "{0} ����".format(a["gene_symbol"])
	elif a["bio_category"] in ["Sv", "PSeqRnaSv"]:
		return "{0}:{1}-{2}:{3} �ں�".format(a["five_prime_gene"], a["five_prime_cds"], a["three_prime_gene"], a["three_prime_cds"])

# ������������ϵ45����ϵ45 + �����ڼ������б���
def judge_var(var, gene_list):
	if ((var["clinic_num_s"] in [4, 5] and var["var_origin"] != "germline") or (var["clinic_num_g"] in [4, 5] and var["var_origin"] == "germline")) and \
		set(re.split(",", var["gene_symbol"])) & set(gene_list):
		return 1
	else:
		return 0

def getgss(jsonDict, var_data, config):
	result = {}
	result["gss"] = listResultToDict.ListToDict(copy.deepcopy(jsonDict["gss"])) if "gss" in jsonDict.keys() and jsonDict["gss"] else {}
	# ��ȡBRCA����
	for gene in ["BRCA1", "BRCA2"]:
		result[gene] = [var for var in var_data if judge_var(var, [gene])]
	# ��ȡHRRͨ·��������������ع̶�����
	#gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDH1","CDK12","CHEK1","CHEK2","FANCA",\
	#			 "FANCL","HDAC2","PALB2","PPP2R2A","PTEN","RAD51B","RAD51C","RAD51D","RAD54L"]
	# �����б����-2023.04.17
	gene_list = ["ATM", "ATR", "BARD1", "BRCA1","BRCA2", "BRIP1", "CDK12","CHEK1","CHEK2", "FANCA",\
				"FANCL", "HDAC2", "MRE11", "NBN", "PALB2","PPP2R2A", "RAD51B","RAD51C","RAD51D","RAD54L"]
	result["summary"] = ", ".join([getVar(var) for var in var_data if judge_var(var, gene_list)])

	# �������-2023.01.03
	if "evi_sum" in result["gss"].keys() and result["gss"]["evi_sum"]:
		# hrd���Ʒ���ת��
		result["gss"]["evi_sum"] = varRegimen(jsonDict, result["gss"]["evi_sum"], config, result["gss"])
		# hrd���Ʒ�������
		regimen_sum = [
			{
				"regimen_name":i["regimen_name"], 
				"evidence_type":i["evidence_type"], 
				"clinical_significance_cn":i["clinical_significance_cn"], 
				"evi_conclusion_simple":i["evi_conclusion_simple"],
				"regimen_name_py":i["regimen_name_py"]
				} 
			for i in result["gss"]["evi_sum"]["regimen_evi_sum"]
			]
		# BRCA���Ʒ�������
		BRCA_data = [var for var in var_data if judge_var(var, ["BRCA1", "BRCA2"])]
		for var in BRCA_data:
			regimen_sum += [
				{
					"regimen_name":i["regimen_name"], 
					"evidence_type":i["evidence_type"], 
					"clinical_significance_cn":i["clinical_significance_cn"], 
					"evi_conclusion_simple":i["evi_conclusion_simple"],
					"regimen_name_py":i["regimen_name_py"]
					} 
				for i in var["evi_sum"]["regimen_evi_sum"]
				]
		# hrd + BRCA���Ʒ���ȥ�ء�����
		regimen_sum_redup = []
		for i in regimen_sum:
			if i not in regimen_sum_redup:
				regimen_sum_redup.append(i)	
		regimen_sum_redup = sorted(regimen_sum_redup, key=lambda i:(i["evi_conclusion_simple"], i["clinical_significance_cn"], i["regimen_name_py"]))

		# HRD�ȼ��жϣ��ж����ݰ���hrd��BRCA�Ľ���� ע��ȼ��ж�ֻ�������ơ�������Ϻ�Ԥ��
		regimen_level = [i["evi_conclusion_simple"] for i in regimen_sum_redup if i["evidence_type"] in ["Predictive", "Prognostic", "Diagnostic"]]
		result["gss"]["level_num"] = 5 if set(["A", "B"]) & set(regimen_level) else 4 if set(["C", "D"]) & set(regimen_level) else 3

		# ���Ʒ�������չʾ
		result["gss"]["regimen"] = {}
		for regimen in regimen_sum_redup:
			if regimen["evidence_type"] not in result["gss"]["regimen"]:
				result["gss"]["regimen"].setdefault(regimen["evidence_type"], [])
			result["gss"]["regimen"][regimen["evidence_type"]].append(regimen)

	return result