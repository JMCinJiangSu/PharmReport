#-*- coding:gbk -*-
import copy
import re
from libs import getRegimenSortRule
from functools import reduce
import itertools
'''
Discription
	
	�ýű�������ȡ���Ʒ���������Ϣ��
	1. ��Ҫ�������򣨰�ҩ������Ⱥ�˳��
	2. �漰ҩ����ܵĲ���: ���졢MSI��TMB��PD-L1��HRD
	* BRCAֻ�漰������(snvindel��MLPA) 

'''

def getRegimen(jsonDict):
	rule = getRegimenSortRule.sort_rule(jsonDict)
	regimen_list = copy.deepcopy(jsonDict["therapeutic_regimen"]) if "therapeutic_regimen" in jsonDict.keys() else []
	regimen_result = []
	# ɾ��TMB��������
	# �ټ�MP IVD���в���TMB-H��ҩ�����
	# ֻ��MP��448��TMB�����ݣ�������Ŀ��Ӱ��
	# ��TMB-H��[var]�У���ɾ���÷��ӱ�־��
	tmb_biomarker = {"biomarker_type" : "TMB-H"}
	for regimen in regimen_list:
		# var�е�hgvs_p������
		if "var" in regimen.keys() and regimen["var"]:
			# var�Ӹ�ȥ��
			regimen["var"] = reduce(lambda x, y:x if y in x else x + [y], [[],]+regimen["var"])
			for a in regimen["var"]:
				if a and "hgvs_p" in a.keys() and a["hgvs_p"]:
					a["hgvs_p"] = a["hgvs_p"].replace("p.", "p.(")+")" if not re.search("=", a["hgvs_p"]) and a["hgvs_p"] != "p.?" and not re.search("\(", a["hgvs_p"]) else a["hgvs_p"]
				# hgvs_p_abbr�������-20220923
				if a and "hgvs_p_abbr" in a.keys() and a["hgvs_p_abbr"]:
					a["hgvs_p_abbr"] = a["hgvs_p_abbr"].replace("p.", "p.(")+")" if not re.search("=", a["hgvs_p_abbr"]) and a["hgvs_p_abbr"] != "p.?" else a["hgvs_p_abbr"]
				if a == tmb_biomarker:
					regimen["var"].remove(a)
				# ��MLPA��varͳһΪgene exon del/gene exon dup(dup ��ʱ��ҩ)-20220902
				if "biomarker_type" in a.keys() and a["biomarker_type"] and re.search("BRCA", a["biomarker_type"]) and re.search("Loss|Gain", a["biomarker_type"]):
					biomarker_type = re.split(":", a["biomarker_type"])
					a["biomarker_type"] = biomarker_type[1]+" "+biomarker_type[2]+" del" if "Loss" in biomarker_type else biomarker_type[1]+" "+biomarker_type[2]+" dup"
		# ��Ӣ�����Ʒ�����ҩ������ܻ���ͬ����������ж�����ͬʱӢ�����ſ�-2022.09.28
		if regimen["regimen_cn"] and regimen["regimen_en"] and regimen["regimen_cn"] == regimen["regimen_en"]:
			regimen["regimen_en"] = ""

		regimen["name"] = regimen["regimen_cn"] if regimen["regimen_cn"] else regimen["regimen_en"]
		# ��Ӧ֢���ݡ�\n�����в�֣���ȥ��-2022.08.18
		adaptation = list(itertools.chain(*[re.split("\n", i.strip()) for i in regimen["adaptation_disease_cn"]])) if regimen["adaptation_disease_cn"] else []
		regimen["adaptation_disease_cn"] = reduce(lambda x, y:x if y in x else x + [y], [[],]+adaptation)	
		if str(regimen["name"]) in rule and regimen["approval_organization"] and ("FDA" in regimen["approval_organization"] or "NMPA" in regimen["approval_organization"]) and "var" in regimen.keys() and regimen["var"]:
			#print(regimen["name"], regimen["approval_organization"], regimen["var"])
			regimen_result.append(regimen)

	# ��߲�ȷ���᲻����drug���ݳ���ʵ��ѡ��ҩ��ķ�Χ
	regimen_result = sorted(regimen_result, key=lambda i:rule.index(i["name"]))
		
	# �ںϸ�ʽ����Ϊgene1:region1-gene2:region2 -2022.08.11
	for regimen in regimen_result:
		if "var" in regimen.keys() and regimen["var"]:
			for var in regimen["var"]:
				if var and "hgvs" in var.keys() and var["hgvs"] and re.search("-", var["hgvs"]):
					five_prime_cds = "-".join(re.split("-", (re.split(":", var["hgvs"])[2]))[:-1]) if not re.search("--", var["hgvs"]) else re.split("_", (re.split("--", var["hgvs"])[0]))[-1]
					three_prime_cds = re.split(":", var["hgvs"])[-1] if not re.search("--", var["hgvs"]) else re.split("_", (re.split("--", var["hgvs"])[1]))[-1]
					# ��һ������-2023.10.19
					# var_hgvs�¸�ʽ��gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, �ɵ�Ϊgene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds����xxx:exon1��xxx:exon2
					five_prime_cds = re.split(":", five_prime_cds)[-1] if re.search(":", five_prime_cds) else five_prime_cds
					three_prime_cds = re.split(":", three_prime_cds)[-1] if re.search(":", three_prime_cds) else three_prime_cds
					# �������-2023.10.19
					five_prime_gene = re.split(":", var["hgvs"])[0] if not re.search("--", var["hgvs"]) else re.split(":", (re.split("--", var["hgvs"])[0]))[0]
					three_prime_gene = re.split(":", re.split("-", var["hgvs"])[-1])[0] if not re.search("--", var["hgvs"]) else re.split(":", (re.split("--", var["hgvs"])[1]))[0]
					var["hgvs2"] = five_prime_gene+":"+five_prime_cds+"-"+three_prime_gene+":"+three_prime_cds+" �ں�"

	# ���£���ҩ����Ӧ֢����drug��ģ���ҩ���������Ʒ����е�-2022.08.15
	# ��ȡdrug�е���Ӧ֢��Ϣ
	drug_adapation = {}
	drug_info = copy.deepcopy(jsonDict["drug"]) if "drug" in jsonDict.keys() else []
	for drug in drug_info:
		drug_name = drug["general_name_cn"].strip() if drug["general_name_cn"] else drug["general_name_en"].strip() if drug["general_name_en"] else ""
		# ��Ӧ֢���ݡ�\n�����в�֣���ȥ��-2022.08.18
		adaptation = list(itertools.chain(*[re.split("\n", i.strip()) for i in drug["adaptation_disease_cn"]])) if drug["adaptation_disease_cn"] else []
		drug_adapation[drug_name] = {
			"adaptation_disease_cn" : reduce(lambda x, y:x if y in x else x + [y], [[],]+adaptation),
			"approval_organization" : drug["approval_organization"] if drug["approval_organization"] else [],
			"trade_name_cn" : drug["trade_name_cn"] if drug["trade_name_cn"] else "",
			"trade_name_en" : drug["trade_name_en"] if drug["trade_name_en"] else "",
			}
	
	# ���Ʒ�����ҩ������Ӧ֢��Ϣ����drug���У���ʹ��drug�����ݣ�û�еĻ�������regimen����
	for regimen in regimen_result:
		if len(regimen["drug_details"]) == 1:
			regimen_name = regimen["regimen_cn"].strip() if regimen["regimen_cn"] else regimen["regimen_en"].strip() if regimen["regimen_en"] else ""
			regimen["adaptation_disease_cn"] = drug_adapation[regimen_name]["adaptation_disease_cn"] if regimen_name in drug_adapation.keys() else regimen["adaptation_disease_cn"] 
			regimen["approval_organization"] = drug_adapation[regimen_name]["approval_organization"] if regimen_name in drug_adapation.keys() else regimen["approval_organization"]
			regimen["trade_name_cn"] = ""
			regimen["trade_name_en"] = ""
			regimen["trade_name_cn"] = drug_adapation[regimen_name]["trade_name_cn"] if regimen_name in drug_adapation.keys() else regimen["trade_name_cn"]
			regimen["trade_name_en"] = drug_adapation[regimen_name]["trade_name_en"] if regimen_name in drug_adapation.keys() else regimen["trade_name_en"]

	return regimen_result