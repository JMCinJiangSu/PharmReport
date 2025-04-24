#-*- coding:gbk -*-
#import os
#import xlrd
import re
from libs.getConfig import getconfigxlsx
from libs.rule import judge_var

'''
Discription
	
	����������ؼ������ 
	���ã�
	1. v3�ټ�ͨ��MP IVD
	2. v3�ټ�ͨ��116
	3. v3����Master ���ټ�ĸ�ʽ��һ����Ϊ�˷�ֹ���������ٴ���һ��
	
'''
def get_io_detect(var_data, config):
	io = {}
	io["result"], io["io_p_summary"], io["io_n_summary"] = io_detect(var_data)
	io["result_116"], io["io_p_summary_116"], io["io_n_summary_116"] = io_detect_for_116(var_data)
	io["io_p_summary_ZJZL"], io["io_n_summary_ZJZL"] = io_detect_for_ZJZL(var_data, config)
	io["result_new_cnvlist"], io["io_p_summary_new_cnvlist"], io["io_n_summary_new_cnvlist"] = io_detect_new_cnvlist(var_data)
	return io

def io_detect(var_data):
	# ���ؽ���е�io_result�������IO��", ".join(io_p_list), ", ".join(io_n_list)�����������С��
	io_result = {}
	# ������ϸ��I/II/����������չ��ر���+��ϵ�²�/�����²��Ա���
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF4","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF4", "FGF19"]

	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]

	for var in level_12_var:
		# ��չʾ�����Ļ���
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("����")
		# ��չʾ�ںϵĻ���
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# �ںϿ��ܻ����exon��ͬ�ϵ㲻ͬ�ı��죬�����л��ظ�����߼Ӹ�ȥ��-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"�ں�" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"�ں�")
		# �������չʾSnvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summaryչʾ
	both_cnv_list = ["CCND1","FGF3","FGF4","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("�ں�", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	#io_p_list = []
	#io_n_list = []
	#for k, v in io_result.items():
	#	io_p_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_P])
	#	io_n_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_N])


	# ����CNV��ͻ��
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF4" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF4/FGF19����")
		#print (", ".join(io_n_list))
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_new_cnvlist(var_data):
	# ���ؽ���е�io_result�������IO��", ".join(io_p_list), ", ".join(io_n_list)�����������С��
	io_result = {}
	# ������ϸ��I/II/����������չ��ر���+��ϵ�²�/�����²��Ա���
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	#io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
	#			 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF4","FGF19"]
	#cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF4", "FGF19"]
	# ɾ��FGF4����-2022307.20
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]

	for var in level_12_var:
		# ��չʾ�����Ļ���
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("����")
		# ��չʾ�ںϵĻ���
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# �ںϿ��ܻ����exon��ͬ�ϵ㲻ͬ�ı��죬�����л��ظ�����߼Ӹ�ȥ��-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"�ں�" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"�ں�")
		# �������չʾSnvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summaryչʾ
	#both_cnv_list = ["CCND1","FGF3","FGF4","FGF19"]
	# ɾ��FGF4����-2022307.20
	both_cnv_list = ["CCND1","FGF3","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("�ں�", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	#io_p_list = []
	#io_n_list = []
	#for k, v in io_result.items():
	#	io_p_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_P])
	#	io_n_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_N])


	# ����CNV��ͻ��
	#if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF4" in io_result.keys() and "FGF19" in io_result.keys():
	#	io_n_list.append("CCND1/FGF3/FGF4/FGF19����")
	# ɾ��FGF4����-2022307.20
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19����")
	#print (", ".join(io_n_list))
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_for_116(var_data):
	# ���ؽ���е�io_result�������IO��", ".join(io_p_list), ", ".join(io_n_list)�����������С��
	io_result = {}
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","FANCA","MRE11","PALB2","MLH1","MSH2","MSH6",\
				 "PMS2","POLE","POLD1","TP53","KRAS","CD274","ARID1A","TERT","CDK12"]
	#io_gene_N = ["EGFR","ALK","CDKN2A","CDKN2B","STK11","JAK1","JAK2","APC","CTNNB1","PTEN"]
	# ����CCND1/FGF3/FGF19������-2023.07.12
	io_gene_N = ["EGFR","ALK","CDKN2A","CDKN2B","STK11","JAK1","JAK2","APC","CTNNB1","PTEN", "CCND1", "FGF3", "FGF19"]
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	for var in level_12_var:
		# ��չʾ�����Ļ���
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in ["CD274", "CCND1", "FGF3", "FGF19"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("����")
		# ��չʾ�ںϵĻ���
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+"�ں�")
			# ��һ�������ɵ� Gene1-ALK(G8:A12)�������� -2023.08.02
			if "ALK_syx" not in io_result.keys():
				io_result.setdefault("ALK_syx", [])
			io_result["ALK_syx"].append("{0}-{1}({2}{3}:{4}{5})��������".format(var["five_prime_gene"], var["three_prime_gene"], var["five_prime_gene"][0], \
								var["five_prime_cds"].replace("exon", ""), var["three_prime_gene"][0], var["three_prime_cds"].replace("exon", "")))
			# ��������-2023.08.02
		# �������չʾSnvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summaryչʾ
	both_cnv_list = ["CCND1", "FGF3", "FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("�ں�", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	# ����CNV��ͻ��
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19����")
	
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_for_ZJZL(var_data, config):
	fdzs_dict, fjzl_database, Data = getconfigxlsx(config)
	# io_result�������IO��
	io_result = {}

	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF4","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF4", "FGF19"]
	# չʾI/II�������������չ��ر���
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	for var in level_12_var:
		# ��չʾ�����Ļ���
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(
				{
					"var_type" : "cnv", 
					"var_info" : "����"
					}
				)
		# ��չʾ�ںϵĻ���
		# ���ܴ��������ںϻ����ڼ�ⷶΧ��������ϵͳjson���ص�gene_symbolΪgene1,gene2���������ʶ��
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(
				{
					"var_type" : "sv", 
					"three_prime_gene" : var["three_prime_gene"], 
					"three_prime_cds" : var["three_prime_cds"], 
					"five_prime_gene" : var["five_prime_gene"], 
					"five_prime_cds" : var["five_prime_cds"]
					}
				)
		# �������չʾSnvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(
				{
					"var_type" : "snvindel", 
					"hgvs_p" : var["hgvs_p"], 
					"hgvs_c" : var["hgvs_c"], 
					"var_origin" : var["var_origin"], 
					"gene_region" : var["gene_region"], 
					"transcript_primary" : var["transcript_primary"]
					}
				)
	# summaryչʾ
	io_inter = []
	for k, v in io_result.items():
		if k not in ["CCND1","FGF3","FGF4","FGF19"]:
			if k in Data.keys():
				io_inter.append({
					"gene_symbol" : k,
					"var_info" : v,
					"inter" : Data[k]
				})
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF4" in io_result.keys() and "FGF19" in io_result.keys():
		io_inter.append({
			"gene_symbol" : "CCND1/FGF3/FGF4/FGF19",
			"var_info" : io_result["CCND1"],
			"inter" : Data["CCND1"] 
		})
	
	io_p_list = [i for i in io_inter if i["gene_symbol"] in io_gene_P]
	io_n_list = [i for i in io_inter if i["gene_symbol"] in io_gene_N or i["gene_symbol"] == "CCND1/FGF3/FGF4/FGF19"]

	return io_p_list, io_n_list