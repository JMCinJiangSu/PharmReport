#-*- coding:gbk -*-
import copy
import re

'''
Discription
	
	�ýű�������ȡ���ܻ����ٴ�������Ϣ��
	����
	1. �й��ٴ�����3��+����ٴ�����3��
	2. ����ѡ��������ļ�ġ��ٴ����ڿ���� 

'''
# ����ȥ�����Ʒ����е�"Drug: /Biological: /Procedure: "���ַ�
def stran_druginfo(druginfo):
	regimen_dict = {
		"Drug" : "_Drug",
		"Biological" : "_Biological",
		"Procedure" : "_Procedure",
		"Other" : "_Other",
		"Device" : "_Device",
		"Behavioral" : "_Behavioral",
		"Radiation" : "_Radiation",
		"Diagnostic Test" : "_Diagnostic Test",
		"Genetic" : "_Genetic"
	}
	for k, v in regimen_dict.items():
		if k in druginfo:
			druginfo = druginfo.replace(k, v)
	return druginfo

def getClinic(jsonDict):
	# json���ص������м�⵽�Ļ��������ٴ�������Ϣ�������ɸѡ
	clinic_list = copy.deepcopy(jsonDict["clinic_trial"])
	gene_list = []
	result = []
	for i in clinic_list:
		if i["gene_symbol"] not in gene_list:
			gene_list.append(i["gene_symbol"])

	for gene in gene_list:
		clinic_cn = [i for i in clinic_list if re.search("CTR", i["clinicaltrial_number"]) and i["gene_symbol"] == gene]
		clinic_en = [i for i in clinic_list if re.search("NCT", i["clinicaltrial_number"]) and "Drug" in i["interventions"] and i["gene_symbol"] == gene]
		
		# �й��ٴ�������Ϣ������ȡ�ٴ����ڿ��������
		if clinic_cn:
			clinic_cn = sorted(clinic_cn, key=lambda i:i["phase"])
			for i in clinic_cn:
				i["interventions"] = [i["interventions"]]
			result += clinic_cn if len(clinic_cn) <= 3 else clinic_cn[0:3]

		# ����ٴ�������Ϣ������ȡ�ٴ����ڿ�����������ҩ�������
		# ���Ʒ�����ʽ��
		# �ٴ����ڸ�Ϊ���ĸ�ʽ
		phase_dict = {
			"Phase 1" : "I��",
			"Phase 2" : "II��",
			"Phase 3" : "III��",
			"Phase 4" : "IV��",
			"Phase 1 Phase 2" : "I��/II��",
			"Phase 2 Phase 3" : "II��/III��",
			"Early Phase 1" : "Early I��",
			"Phase 1;Phase 2" : "I��/II��"
		}
		if clinic_en:
			clinic_en = sorted(clinic_en, key=lambda i:i["phase"], reverse = True)
			for i in clinic_en:
				druginfo = stran_druginfo(i["interventions"])
				druglist = []
				for j in re.split("_", druginfo):
					if re.search("Drug", j):
						j = j.replace("Drug: ","")
						druglist.append(j)
				i["interventions"] = druglist
				i["phase"] = phase_dict.get(i["phase"], i["phase"])
			result += clinic_en if len(clinic_en) <= 3 else clinic_en[0:3]
			
	return result