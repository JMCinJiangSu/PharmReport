#-*- coding:gbk -*-
import copy
import re
from libs import listResultToDict
'''
Discription
	
	ҩ���Ⱥ�˳���ƶ�������������Ʒ��������ж���ҩ����ִ��� 

'''

def sort_rule(jsonDict):
	rule = []
	#gss_copy = copy.deepcopy(jsonDict["gss"]) if "gss" in jsonDict.keys() and jsonDict["gss"] else {}
	# gssӦ�÷��ظ��ֵ䣬�°汾������б���߼Ӹ�����-2023.09.26
	gss_copy = listResultToDict.ListToDict(copy.deepcopy(jsonDict["gss"])) if "gss" in jsonDict.keys() and jsonDict["gss"] else {}
	# �������-2023.09.26
	var_list = copy.deepcopy(jsonDict["snvindel"]+\
							 jsonDict["cnv"]+\
							 jsonDict["sv"]+\
							 jsonDict["rna_sv"]+\
							 jsonDict["knb"]+\
							 jsonDict["msi"]+\
							 jsonDict["tmb"]+\
							 jsonDict["pdl1"]+\
							 jsonDict["mlpa"]+\
							 jsonDict["hrd"]+\
							 [gss_copy])
	drug_list = []
	for var in var_list:
		if "evi_sum" in var.keys():
			drug_list += [evi["regimen_name"] for evi in var["evi_sum"] if evi["regimen_name"] and re.search("Sensitive", evi["clinical_significance"])]
	for drug in drug_list:
		for i in re.split("\+", drug):
			if i.strip() not in rule:
				rule.append(i.strip())
	return rule