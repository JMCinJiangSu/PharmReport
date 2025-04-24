#-*- coding:gbk -*-
import copy
import re
from libs import getDrugSortRule
from functools import reduce
import itertools
'''
Discription
	
	�ýű�������ȡ����ҩ�������Ϣ��
	1. ��Ҫ�������򣨰�ҩ������Ⱥ�˳��
	2. �漰ҩ����ܵĲ���: ���졢MSI��TMB��PD-L1��HRD
	* BRCAֻ�漰������(snvindel��MLPA) 

'''

def getDrug(jsonDict):
	rule = getDrugSortRule.sort_rule(jsonDict)
	drug_list = copy.deepcopy(jsonDict["drug"])
	drug_result = []
	# ɾ��TMB��������
	# �ټ�MP IVD���в���TMB-H��ҩ�����
	# ֻ��MP��448��TMB�����ݣ�������Ŀ��Ӱ��
	# ��TMB-H��[var]�У���ɾ���÷��ӱ�־��
	tmb_biomarker = {"biomarker_type" : "TMB-H"}
	for drug in drug_list:
		drug["name"] = drug["general_name_cn"] if drug["general_name_cn"] else drug["general_name_en"]
		#var �е�hgvs_p������
		if drug["var"]:
			for a in drug["var"]:
				if a and "hgvs_p" in a.keys() and a["hgvs_p"]:
					a["hgvs_p"] = a["hgvs_p"].replace("p.", "p.(")+")" if not re.search("=", a["hgvs_p"]) and a["hgvs_p"] != "p.?" and not re.search("\(", a["hgvs_p"]) else a["hgvs_p"]
				# hgvs_p_abbr�������-20220923
				if a and "hgvs_p_abbr" in a.keys() and a["hgvs_p_abbr"]:
					a["hgvs_p_abbr"] = a["hgvs_p_abbr"].replace("p.", "p.(")+")" if not re.search("=", a["hgvs_p_abbr"]) and a["hgvs_p_abbr"] != "p.?" else a["hgvs_p_abbr"]
				if a == tmb_biomarker:
					drug["var"].remove(a)
				# ��MLPA��varͳһΪgene exon del/gene exon dup(dup ��ʱ��ҩ)-20220902
				if "biomarker_type" in a.keys() and a["biomarker_type"] and re.search("BRCA", a["biomarker_type"]) and re.search("Loss|Gain", a["biomarker_type"]):
					biomarker_type = re.split(":", a["biomarker_type"])
					a["biomarker_type"] = biomarker_type[1]+" "+biomarker_type[2]+" del" if "Loss" in biomarker_type else biomarker_type[1]+" "+biomarker_type[2]+" dup"

		# ��ҩ����������޶���FDA��NMPA��
		drug["approval_organization"] = list(set(drug["approval_organization"]) & set(["FDA", "NMPA"]))
		# ��Ӧ֢ȥ��-2022.08.10
		#drug["adaptation_disease_cn"] = reduce(lambda x, y:x if y in x else x + [y], [[],]+drug["adaptation_disease_cn"]) if drug["adaptation_disease_cn"] else []
		# ��Ӧ֢���ݡ�\n�����в�֣���ȥ��-2022.08.18
		adaptation = list(itertools.chain(*[re.split("\n", i.strip()) for i in drug["adaptation_disease_cn"]])) if drug["adaptation_disease_cn"] else []
		drug["adaptation_disease_cn"] = reduce(lambda x, y:x if y in x else x + [y], [[],]+adaptation)		

		if str(drug["name"]) in rule and drug["approval_organization"] and drug["var"]:
			drug_result.append(drug)

	# ��߲�ȷ���᲻����drug���ݳ���ʵ��ѡ��ҩ��ķ�Χ
	drug_result = sorted(drug_result, key=lambda i:rule.index(i["name"]))

	return drug_result