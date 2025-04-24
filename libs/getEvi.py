#-*- coding:gbk -*-
import re
from libs.getConfig import senseStran
from libs.getPinyin import topinyin
from libs.getInterRef import getRef_from_inter

'''
Discription	

	�������Ӧ�����Ʒ�����֤��������ת��Ϊ�ʺ����ĸ�ʽ�� 

'''

# ��֤ͬ�������ĺϲ����Ʒ���չʾ
def merge_Predictive_evi(datainfo):
	tmp_dict = {}
	for evi in datainfo:
		tmp_dict.setdefault(evi["evi_interpretation"], [])
		tmp_dict[evi["evi_interpretation"]].append(
			{
				"regimen_name" : evi["regimen_name"],
				"evi_conclusion_simple" : evi["evi_conclusion_simple"],
				"clinical_significance_cn" : evi["clinical_significance_cn"],
				"regimen_name_py" : evi["regimen_name_py"],
				"evi_conclusion" : evi["evi_conclusion"]
			}
		)

	merge_result = []
	for k, v  in tmp_dict.items():
		merge_result.append(
			{
				"regimen_name" : "��".join([i["regimen_name"] for i in v]),
				"evi_conclusion_simple" : "/".join([i["evi_conclusion_simple"] for i in v]),
				"clinical_significance_cn" : "/".join([i["clinical_significance_cn"] for i in v]),
				"regimen_name_py" : "/".join([i["regimen_name_py"] for i in v]),
				"evi_interpretation" : k,
				"evi_conclusion" : "/".join([i["evi_conclusion"] for i in v])
			}
		)

	return merge_result

def varRegimen(jsonDict, evi_sum, config, var):
	data = {}
	data["refer_evi"] = []
	data["evi_split"] = {}
	data["refer_evi_risk"] = []

	# 1. �����ֶΡ�֤�������������������
	for evi in evi_sum:
		## �����ֶ�-���ڱ���չʾ-1.�ٴ�����ת��Ϊ����
		evi["clinical_significance_cn"] = senseStran(config).get(evi["clinical_significance"], evi["clinical_significance"])
		## �����ֶ�-���ڱ���չʾ-2.�ȼ���A0��A1��ת��ΪA
		evi["evi_conclusion_simple"] = evi["evi_conclusion"][0] if evi["evi_conclusion"] else ""
		## �����ֶ�-���ڱ���չʾ-3.�������ף����Ʒ�������������������FDA��MNPA��EMA
		evi["regimen_refer_agency_ZJZL"] = "/".join(list(set(["FDA","NMPA","EMA"]) & set(re.split(",", evi["regimen_refer_agency"])))) \
										   if "regimen_refer_agency" in evi.keys() and evi["regimen_refer_agency"] \
										   else "-"
		
		## �����ֶ�-��������-1. ���Ʒ���ת��Ϊƴ�������ں�������
		evi["regimen_name_py"] = topinyin(evi["regimen_name"]) if evi["regimen_name"] else "0"
		## �����ֶ�-��������-2. ���ڶ����к���ҩ������������ǰ��ҩ��
		evi["sense_rule"] = "0" if re.search("Sensitive", evi["clinical_significance"]) else \
							"1" if re.search("Resistant", evi["clinical_significance"]) else \
							evi["clinical_significance"]
		
		# ��������-1.֤������ȥ��ĩβ�ո�
		evi["evi_interpretation"] = evi["evi_interpretation"].strip() if evi["evi_interpretation"] else ""
		
	# ��֤�ݽ������򣬰����Ʒ����ȼ�������/��ҩ�����Ʒ���ƴ������ĸͳһ��д��������������
	#evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	# RET �ںϣ���ҩ�������-��֮ǰ�Ŀ������ᡢ�������ᡢ�����������Ϊ�������ᡢ�������ᡢ��������-2023.05.11
	drug_list_rule = {
		"��������" : 0, 
		"��������" : 1,
		"��������" : 2
		}
	if "bio_category" in var.keys() and var["bio_category"] in ["Sv", "PSeqRnaSv"] and (var["five_prime_gene"] == "RET" or var["three_prime_gene"] == "RET"):
		evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], drug_list_rule.get(i["regimen_name"], 3) , i["regimen_name_py"].upper()))
	else:
		evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	# ����������-2023.05.11

	# 2. ��������֤�ݽ��вο�������ȡ��֤�ݷ��ࡢ����֤�ݺϲ���
	for evi in evi_sum:
		## ��ȡ�ο�����-֤������
		data["refer_evi"].extend(getRef_from_inter(jsonDict, evi["evi_interpretation"]))

		## ��ȡ�ο�����-�Ŵ�����
		if evi["evidence_type"] == "Predisposing":
			data["refer_evi_risk"].extend(getRef_from_inter(jsonDict, evi["evi_interpretation"]))
		
		## ����evi_split����֤�ݽ��з��࣬���ơ�������ϡ�Ԥ���Ŵ�����
		if evi["evidence_type"] not in data["evi_split"].keys():
			data["evi_split"].setdefault(evi["evidence_type"], [])
		data["evi_split"][evi["evidence_type"]].append(evi)

		## Ŀǰ����Predictiveʱ��Ҫ�ϲ���֤ͬ�ݣ��������á�Predictive_merge���ֶΣ����ݱ����������ѡ��
		if "Predictive" in data["evi_split"].keys():
			data["evi_split"]["Predictive_merge"] = merge_Predictive_evi(data["evi_split"]["Predictive"])

	# һЩ��������
	## 1. ��������֤��
	data["regimen_evi_sum"] = evi_sum
	## 2. ����A���������Ʒ���
	data["regimen_FDA_S"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"]
			} 
		for var in evi_sum if re.search("Sensitive",var["clinical_significance"]) and var["evi_conclusion_simple"] == "A"
		]
	## 3. ���ܷ�A���������Ʒ���
	data["regimen_noFDA_S"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"]
			} 
		for var in evi_sum if re.search("Sensitive",var["clinical_significance"]) and var["evi_conclusion_simple"] != "A"
		]

	## 4. �����������Ʒ���-����֤������-2023.08.15
	data["regimen_S"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"],
			"evi_interpretation" : var["evi_interpretation"]
			} 
		for var in evi_sum if re.search("Sensitive",var["clinical_significance"])
		]

	## 5. ������ҩ���Ʒ���-����֤������-2023.08.15
	data["regimen_R"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"],
			"evi_interpretation" : var["evi_interpretation"]
			} 
		for var in evi_sum if re.search("Resistant",var["clinical_significance"])
		]

	## 6. �������Ʒ����б�ת��Ϊ�ַ���
	data["regimen_S_str"] = "��".join([i["regimen_name"] for i in data["regimen_S"]])

	## 7. ��ҩ���Ʒ����б�ת��Ϊ�ַ���
	data["regimen_R_str"] = "��".join([i["regimen_name"] for i in data["regimen_R"]])

	# �����ɣ�ҩ����
	# I�� FDA/NMPA/NCCNҩ�A0/A1/A2/C3 
	# I�� �ٴ�����ҩ�B1/B2/B3/C1/C2/C4/D1/D2/D3/D4/D5/D6
	# II�� FDA/NMPA��׼���������ֵ�ҩ� C3
	# II�� �ٴ�����ҩ� C1/C2/C4/D1/D2/D3/D4/D5/D6
	# ���ܽ�Ϊ�����ֶΣ��������ٴ�����ҩ�I/II��ɸ��ݱ���ȼ����ж���Ӧ�ò������
	SYX_regimen_appr = [evi for evi in evi_sum if evi["evi_conclusion"] in ["A0","A1","A2","C3"]]
	SYX_regimen_clinic = [evi for evi in evi_sum if evi["evi_conclusion"] in ["B1","B2","B3","C1","C2","C4","D1","D2","D3","D4","D5","D6"]]
	#�������Ʒ�����Ԥ�󡢸�����ϡ����յȲ��
	data["SYX_regiman_appr"] = {}
	for evi in SYX_regimen_appr:
		if evi["evidence_type"] not in  data["SYX_regiman_appr"].keys():
			data["SYX_regiman_appr"].setdefault(evi["evidence_type"], [])
		data["SYX_regiman_appr"][evi["evidence_type"]].append(evi)
		
	data["SYX_regimen_clinic"] = {}
	for evi in SYX_regimen_clinic:
		if evi["evidence_type"] not in  data["SYX_regimen_clinic"].keys():
			data["SYX_regimen_clinic"].setdefault(evi["evidence_type"], [])
		data["SYX_regimen_clinic"][evi["evidence_type"]].append(evi)

	return data