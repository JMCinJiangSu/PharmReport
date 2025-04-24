#-*- coding:gbk -*-
import jinja2
import re
from docxtpl import RichText
import jinja2.filters
from libs.rule import nofound_genelist
from libs.rule import decimal_float

'''
�Զ���jinja2������
'''

### ������ ###
# 1. ������ɽ����ҽԺCP40��ҩ�����.join()չʾΪһ��
def regimen_sum(a):
	result = []
	evi_sum = a["evi_sum"]
	if "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Predictive"]:
				result.append("{0}��{1}��{2}����".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Prognostic"]:
				result.append("{0}��{1}��{2}����".format("Ԥ��"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Diagnostic"]:
				result.append("{0}��{1}��{2}����".format("�������", " / ", i["evi_conclusion_simple"]))
	return "��".join(result) if result else "-"
jinja2.filters.FILTERS["regimen_sum_filter"] = regimen_sum 

# 2. ������������HRR-��BRCA����-������
# gene_symbol gene_region hgvs_c hgvs_p, clinic_num_g.strans
def var_brca_sum(var_list):
	result = []
	clinic_trans = {5 : "�²��Ա���", 4 : "�����²��Ա���"}
	for var in var_list:
		if var["type"] == "Loss":
			result.append("{0} {1} del, {2}".format(var["gene_symbol"], var["value"], "�����²��Ա���"))

		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}, {4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], clinic_trans[var["clinic_num_g"]]))
			else:
				result.append("{0} {1} {2}, {3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], clinic_trans[var["clinic_num_g"]]))
	return "��".join(result)
jinja2.filters.FILTERS["var_brca_sum"] = var_brca_sum

# 3. ������������HRR-��BRCA����-����������
def inter_brca_sum(var_list):
	result_dict = {}
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			var["clinic_num_g"] = 4
		elif var["type"] == "Gain":
			var["clinic_num_g"] = 3

		if var["clinic_num_g"] not in result_dict.keys():
			result_dict.setdefault(var["clinic_num_g"], [])
		result_dict[var["clinic_num_g"]].append(var)
	clinic_trans = {5:"�²��Ա���", 4:"�����²��Ա���", 3:"���岻��ȷ����", 2:"�������Ա���"}
	for i in [5,4,3,2]:
		if i in result_dict.keys():
			result.append(str(len(result_dict[i]))+"��"+clinic_trans[i])
	return "��".join(result)
jinja2.filters.FILTERS["inter_brca_sum"] = inter_brca_sum

# ������116-���С��
def summary_SYX_116(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			result.append(var["gene_symbol"]+" ����")
		elif var["bio_category"] == "Sv":
			result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" �ں�")
	return ", ".join(result)
jinja2.filters.FILTERS["summary_SYX_116"] = summary_SYX_116

# ��ɽ����HRDC�����С��-2023.06.02
def summary_ZSRM_hrd(info):
	# ���ԣ����ڱ����ͼ����������BRCA1����BRCA1 gene_region hgvs_c hgvs_p I����죬HRD״̬�����Լ�TP53����TP53 gene_region hgvs_c hgvs_p II����졣
	# ���ԣ����ڱ����ͼ�������δ����²��Ի�����²��Ա��졣
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I��", 4 : "II��"}
	
	brca_result = []
	for var in var_dict["ec_type"]["BRCA1_level12"] + var_dict["ec_type"]["BRCA2_level12"]:
		if var["hgvs_p"] != "p.?":
			brca_result.append("{0}����{0} {1} {2} {3} {4}����".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			brca_result.append("{0}����{0} {1} {2} {3}����".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	brca_result_str = "��".join(brca_result)
	
	hrd_result = "HRD״̬����" if hrd["var_id"] == "HRD+" else "HRD״̬����"

	other_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["gene_symbol"] not in ["BRCA1", "BRCA2"]:
			if var["hgvs_p"] != "p.?":
				other_result.append("{0}����{0} {1} {2} {3} {4}����".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
			else:
				other_result.append("{0}����{0} {1} {2} {3}����".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	other_result_str = ",".join(other_result)
	hrd_result_str = "�Լ�".join([hrd_result, other_result_str])

	result = []
	if brca_result_str:
		result.append(brca_result_str)
	if hrd_result_str:
		result.append(hrd_result_str)

	# ��װ-����/����
	if var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		return "���ڱ����ͼ����������{0}��".format("��".join(result))
	else:
		return "���ڱ����ͼ�������δ����²��Ի������²��Ա��졣"
jinja2.filters.FILTERS["summary_ZSRM_hrd"] = summary_ZSRM_hrd

##############################################################################################################

### ͨ���� ###

### �����б��� ###
# 1. ��ϸ������������������չ��ر����б�
def return_somatic_onconodrug_list(info):
	var = info[0]
	sample = info[1]
	if re.search("Master", sample["prod_names"]):
		if not sample["control_sample_id"]:
			return var["var_somatic"]["level_onco_nodrug"] + var["var_germline_nodrug"]
		else:
			return var["var_somatic"]["level_onco_nodrug"]
	else:
		return var["var_somatic"]["level_onco_nodrug"]
jinja2.filters.FILTERS["somatic_onconodrug_list"] = return_somatic_onconodrug_list

# 2. ��ϸ����������ϸ��/δ֪��Դ���죬III������б�
def return_somatic_level3_list(info):
	var = info[0]
	sample = info[1]
	if sample["prod_names"] == "BPTM��5����":
		return var["ec_type"]["POLE_level3"]+var["ec_type"]["TP53_level3"]+var["ec_type"]["BRCA1_level3"]+var["ec_type"]["BRCA2_level3"]
	elif sample["prod_names"] == "PTM��3����":
		return var["ec_type"]["POLE_level3"]+var["ec_type"]["TP53_level3"]
	else:
		return var["var_somatic"]["level_III"]
jinja2.filters.FILTERS["somatic_level3_list"] = return_somatic_level3_list

# 3. ������������б�δ��⵽����Ļ���ҲҪչʾ
def return_gLS5_list(gLS5):
	result = []
	gene_list = ["EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	for gene in gene_list:
		if gene in gLS5.keys() and gLS5[gene]:
			result.extend(gLS5[gene])
		else:
			result.append({"gene_symbol" : gene})
	return result
jinja2.filters.FILTERS["gLS5_list"] = return_gLS5_list
	
# 4. ����BRCA����������δ��⵽����Ļ���ҲҪչʾ
def return_brca_list(brca):
	result = []
	brca1_list = brca["snv_s"]["B1_L5"] + brca["snv_s"]["B1_L4"]+brca["mlpa"]["B1_Loss"]+brca["mlpa"]["B1_Gain"]+brca["snv_s"]["B1_L3"]
	brca2_list = brca["snv_s"]["B2_L5"] + brca["snv_s"]["B2_L4"]+brca["mlpa"]["B2_Loss"]+brca["mlpa"]["B2_Gain"]+brca["snv_s"]["B2_L3"]
	if brca1_list:
		result.extend(brca1_list)
	else:
		result.append({"gene_symbol" : "BRCA1"})
	
	if brca2_list:
		result.extend(brca2_list)
	else:
		result.append({"gene_symbol" : "BRCA2"})
	return result
jinja2.filters.FILTERS["brca_list"] = return_brca_list
	

#----------------------------------------------------------------------------------------------------------------

### �������� ###
# 1. ��ϸ����������ϸ��/δ֪��Դ���� I��II���������б�
def return_somatic_level12_inter(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	if re.search("BRCA", sample["prod_names"]):
		# var_origin != "germline"
		brca_result = []
		for var in var_brca["snv_s"]["B1_L5"] + var_brca["snv_s"]["B2_L5"] + var_brca["snv_s"]["B1_L4"] + var_brca["snv_s"]["B2_L4"]:
			if var["var_origin"] != "germline":
				brca_result.append(var)
		return brca_result
	elif re.search("BPTM��5����", sample["prod_names"]):
		return var["ec_type"]["POLE_level12"] + var["ec_type"]["TP53_level12"] + var["ec_type"]["BRCA1_level12"] + var["ec_type"]["BRCA2_level12"]
	elif re.search("PTM��3����", sample["prod_names"]):
		return var["ec_type"]["POLE_level12"] + var["ec_type"]["TP53_level12"]
	# �����MP�������ʱ��var_originΪgermline��������ҩ��ҲҪչʾ
	elif re.search("Master Panel����֯��", sample["prod_names"]) and not sample["control_sample_id"]:
		return var["var_somatic"]["level_I"] + var["var_germline"]["regimen_level_I"] + var["var_somatic"]["level_II"] + var["var_germline"]["regimen_level_II"]
	else:
		return var["var_somatic"]["level_I"] + var["var_somatic"]["level_II"]
jinja2.filters.FILTERS["somatic_level12_inter"] = return_somatic_level12_inter

# 2. ��ϵ��������ϵ����4/5���б����ڽ��
def return_germline_level45_inter(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	if re.search("BRCA", sample["prod_names"]):
		return var_brca["snv_m"]["B1_G_L5"]+var_brca["snv_m"]["B2_G_L5"]+var_brca["snv_m"]["B1_G_L4"]+var_brca["snv_m"]["B2_G_L4"]+\
			   var_brca["mlpa"]["B1_Loss"] + var_brca["mlpa"]["B2_Loss"]
	elif re.search("HRR", sample["prod_names"]):
		return var["var_germline"]["level_5"]+var["var_germline"]["level_4"]+var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"]
	else:
		return var["var_germline"]["level_5"]+var["var_germline"]["level_4"]
jinja2.filters.FILTERS["germline_level45_inter"] = return_germline_level45_inter

#3. ��ϵ��������ϵ����3���б����ڽ��
def return_germline_level3_inter(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	if re.search("BRCA", sample["prod_names"]):
		return var_brca["mlpa"]["B1_Gain"]+var_brca["mlpa"]["B2_Gain"]+var_brca["snv_m"]["B1_G_L3"]+var_brca["snv_m"]["B2_G_L3"]
	elif re.search("HRR", sample["prod_names"]):
		return var_brca["mlpa"]["B1_Gain"]+var_brca["mlpa"]["B2_Gain"]+var["var_germline"]["level_3"]
	else:
		return var["var_germline"]["level_3"]
jinja2.filters.FILTERS["germline_level3_inter"] = return_germline_level3_inter

#----------------------------------------------------------------------------------------------------------------

### �������е����� ###
# 1. ���ػ��򣬵�������
def return_gene_symbol(var):
	result = []
	if "bio_category" in var.keys() and var["bio_category"]:
		if var["bio_category"] in ["Snvindel", "Cnv"]:
			result.append(var["gene_symbol"])
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if "," in var["gene_symbol"]:
				if var["five_prime_gene"] != var["three_prime_gene"]:
					result.append(var["five_prime_gene"])
					result.append(var["three_prime_gene"])
				else:
					result.append(var["five_prime_gene"])
			else:
				result.append(var["gene_symbol"])
		else:
			result.append(var["gene_symbol"])
	else:
		# gene_symbol ���ܲ����ڣ���߼Ӹ�����-2023.09.22
		if "gene_symbol" in var.keys():
			result.append(var["gene_symbol"])
		# �������-2023.09.22
	rt = RichText()
	rt.add("\n".join(result), italic=True)
	return rt
jinja2.filters.FILTERS["gene_symbol"] = return_gene_symbol

# 2. ��������
def var_info(var):
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_region"]+" "+var["hgvs_c"])
			result.append(var["transcript_primary"])
		elif var["bio_category"] == "Cnv":
			if var['cnv_type'] == 'loss':
				result.append("ȱʧ")
			else:
				result.append("����")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 ��Ծ")
				result.append(var["five_prime_transcript"])
			else:
				result.append("{0}:{1}-{2}:{3}�ں�".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"]))
				result.append("{0}/{1}".format(var["five_prime_transcript"], var["three_prime_transcript"]))
		elif var["bio_category"] == "PMLPA":
			if var["type"] == "Loss":
				result.append(var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["value"]+" dup")
	rt = RichText()
	rt.add("\n".join(result))
	return rt
	#return result
jinja2.filters.FILTERS["var_info"] = var_info


# 3. ���ر����ȣ�����BRCA�ȣ���ϵչʾ����/�Ӻϣ���ϸ��չʾƵ�ʣ�MLPA��-
def freq_stran(info):
	'''
	��ϸ��snvindel��freq_str
	��ϵsnvindel��������Ŀfreq�� �Ϻ���Ŀfreq_rc��ע��MP��Ҫ�����Ƿ��ж����������޶�����������ԴΪgermline��ֱ��չʾfreq_str
	sv��CP40 copies������freq_str
	PSeqRnaSv��freq
	'''
	var = info[0]
	sample = info[1]
	if "var_origin" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["var_origin"] == "germline":
				if "Master" in sample["prod_names"]:
					if sample["control_sample_id"]:
						return "����" if float(var["freq_rc"]) >= 0.85 else "�Ӻ�"
					else:
						return var["freq_str"]
					pass 
				elif re.search("116|76|25|21|18", sample["prod_names"]):
					print (var["freq_rc"])
					return "����" if var["freq_rc"] and float(var["freq_rc"]) >= 0.85 else "�Ӻ�" if var["freq_rc"] and float(var["freq_rc"]) < 0.85 else "δ��ȡ��freq_rc��"
				else:
					return "����" if float(var["freq"]) >= 0.85 else "�Ӻ�"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "Cnv":
			return var["cn_mean"]
		elif var["bio_category"] == "Sv":
			# ����OncoPro����֯�� XW1402
			if re.search("Classic|CRC12|OncoPro����֯��", sample["prod_names"]):
				return str(var["copies"])+" copies"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "PSeqRnaSv":
			return str(var["freq"])+" copies"
	else:
		if var["type"] in ["Loss", "Gain"]:
			return "/"
jinja2.filters.FILTERS["freq_stran"] = freq_stran

# 4. �ٴ�����-ҩ��
def significance_regimen(var):
	result = []
	#rt = RichText()
	evi_sum = var["evi_sum"]
	if evi_sum and "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Predictive"]:
				result.append("{0}��{1}��{2}����".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Prognostic"]:
				result.append("{0}��{1}��{2}����".format("Ԥ��"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Diagnostic"]:
				result.append("{0}��{1}��{2}����".format("�������", " / ", i["evi_conclusion_simple"]))
		#rt.add("\n".join(result)) 
	if not result:
		result = ["-"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
	#return result
jinja2.filters.FILTERS["significance_regimen"] = significance_regimen

# 5. �����ٴ����壬����BRCA��
def clinic_stran(info):
	var = info[0]
	s_type = info[1]
	g_dict = {
		5 : "�²���",
		4 : "�����²���",
		3 : "���岻��ȷ",
		2 : "��������",
		1 : "����"
	}
	s_dict1 = {
		5 : "I��-ǿ�ٴ�����",
		4 : "II��-Ǳ���ٴ�����",
		3 : "III��-�ٴ����岻��",
		2 : "IV��-����/��������",
		1 : "IV��-����/��������"
	}
	s_dict2 = {
		5 : "I��",
		4 : "II��",
		3 : "III��",
		2 : "IV��",
		1 : "IV��"
	}
	# ҩ���������
	s_dict3 = {
		5 : "I�����",
		4 : "II�����",
		3 : "III�����",
		2 : "IV�����",
		1 : "IV�����"
	}
	s_dict = s_dict1 if s_type == "s1" else s_dict2 if s_type == "s2" else s_dict3 if s_type == "s3" else {}
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			return g_dict.get(var["clinic_num_g"], "")
		else:
			# ����������չ��ع�ΪIII��
			if var["clinic_num_s"] in [4, 5] and (not var["evi_sum"]["evi_split"] or \
				(var["evi_sum"]["evi_split"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()))) :
				return s_dict.get(var[3])
			else:
				return s_dict.get(var["clinic_num_s"], "")
	else:
		if var["type"] == "Loss":
			return g_dict.get(4)
		elif var["type"] == "Gain":
			return g_dict.get(3)
jinja2.filters.FILTERS["clinic_stran"] = clinic_stran

# 5.1 �����ٴ����壬����MP
def clinic_stran_MP(info):
	var = info[0]
	s_type = info[1]
	sample = info[2]
	g_dict = {
		5 : "�²���",
		4 : "�����²���",
		3 : "���岻��ȷ",
		2 : "��������",
		1 : "����"
	}
	s_dict1 = {
		5 : "I��-ǿ�ٴ�����",
		4 : "II��-Ǳ���ٴ�����",
		3 : "III��-�ٴ����岻��",
		2 : "IV��-����/��������",
		1 : "IV��-����/��������"
	}
	s_dict2 = {
		5 : "I��",
		4 : "II��",
		3 : "III��",
		2 : "IV��",
		1 : "IV��"
	}
	s_dict = s_dict1 if s_type == "s1" else s_dict2 if s_type == "s2" else {}
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			if sample["control_sample_id"]:
				return g_dict.get(var["clinic_num_g"], "")
			else:
				return s_dict.get(var["clinic_num_s"], "")
		else:
			# ����������չ��ع�ΪIII��
			if var["clinic_num_s"] in [4, 5] and (not var["evi_sum"]["evi_split"] or \
				(var["evi_sum"]["evi_split"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()))) :
				return s_dict.get(var[3])
			else:
				return s_dict.get(var["clinic_num_s"], "")
	else:
		if var["type"] == "Loss":
			return g_dict.get(4)
		elif var["type"] == "Gain":
			return g_dict.get(3)
jinja2.filters.FILTERS["clinic_stran_MP"] = clinic_stran_MP


# 6. ������Դת��
def var_origin_stran(var):
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			return "��ϵ"
		elif var["var_origin"] == "somatic":
			return "��ϸ��"
		else:
			return "����"
	else:
		if var["type"] in ["Loss", "Gain"]:
			return "��ϵ"
jinja2.filters.FILTERS["var_origin_stran"] = var_origin_stran

# 7. ���С�����ͷ-���-������BRCA������
def title_freq_brca(prod_names):
	if prod_names in ["BRCA1/BRCA2��ȫѪ��", "BRCA1/2�������ӣ�", "HRR��ȫѪ��", "�����ۺ���"]:
		return  "������"
	elif prod_names in ["BRCA1/BRCA2����֯��"]:
		return "���"
	elif prod_names in ["BRCA1/BRCA2����֯ ȫѪ��"]:
		return "������/���"
jinja2.filters.FILTERS["title_freq_brca"] = title_freq_brca

# 8. III ������ͷ-���-��ϸ��
def title_freq_III(sample):
	if sample["prod_names"] in ["BRCA1/BRCA2��ȫѪ��", "BRCA1/2�������ӣ�", "HRR��ȫѪ��", "�����ۺ���"]:
		return  "������"
	elif sample["prod_names"] in ["BRCA1/BRCA2����֯��","BPTM��5����", "PTM��3����", "HRR����֯��", "HRD Complete����֯��", "10����ѪҺ��", \
								  "Pan116��ѪҺ��", "LC76��ѪҺ��", "CRC25��ѪҺ��", "TC21��ѪҺ��", "GA18��ѪҺ��"]:
		return "���"
	elif sample["prod_names"] in ["BRCA1/BRCA2����֯ ȫѪ��", "HRR����֯ ȫѪ��"]:
		return "���"
	elif sample["prod_names"] in ["10������֯��", "Classic Panel", "CRC12-MSI",\
								  "Pan116����֯��", "LC76����֯��", "CRC25����֯��", "TC21����֯��", "GA18����֯��",\
								  "Master Panel����֯��", "Master Panel��ѪҺ��"]:
		return "���/������"
jinja2.filters.FILTERS["title_freq_III"] = title_freq_III

# 9. ����������չ��ر���-���
# Ŀǰ�漰����Ŀ��HRR��֯/HRR���/HRD/LC10��֯/LC10ѪҺ/PAN116��֯/PAN116ѪҺ/LC76��֯/LC76ѪҺ/CRC25��֯/CRC25ѪҺ/TC21��֯/TC21ѪҺ/GA18��֯/GA18ѪҺ
# Master��֯/MasterѪҺ/CP40/CRC12
# ���/������
def title_freq_onconodrug(sample):
	if sample["prod_names"] in ["HRR����֯��", "HRR����֯ ȫѪ��", "HRD Complete����֯��", "10����ѪҺ��", \
							   "Pan116��ѪҺ��", "LC76��ѪҺ��", "CRC25��ѪҺ��", "TC21��ѪҺ��", "GA18��ѪҺ��"]:
		return "���"
	elif sample["prod_names"] in ["Pan116����֯��", "LC76����֯��", "CRC25����֯��", "TC21����֯��", "GA18����֯��",\
								 "Master Panel����֯��", "Master Panel��ѪҺ��", "Classic Panel", "CRC12-MSI", "10������֯��"]:
		return "���/������"
jinja2.filters.FILTERS["title_freq_onconodrug"] = title_freq_onconodrug

# 10. �������Ʊ�-���
# ������MP/116/LC10/CP40/CRC12
def title_freq_targetRegimen(sample):
	freq = ["���"]
	if re.search("��֯", sample["prod_names"]):
		freq.append("������")
	if sample["prod_names"] in ["Classic Panel", "CRC12-MSI"]:
		freq.append("������")
	if sample["prod_names"] == "Master Panel��ѪҺ��":
		freq.append("������")
	if sample["control_sample_id"]:
		freq.append("������")
	return "/".join(freq)
jinja2.filters.FILTERS["title_freq_targetRegimen"] = title_freq_targetRegimen

# 11 PARP������ܱ�-��ȱ�ͷ
def title_parp(sample):
	if sample["prod_names"] in ["BRCA1/BRCA2��ȫѪ��", "HRR��ȫѪ��"]:
		return "������"
	elif sample["prod_names"] in ["HRR����֯��", "HRD Complete����֯��"]:
		return "���"
	elif sample["prod_names"] in ["HRR����֯ ȫѪ��"]:
		return "������/���"
jinja2.filters.FILTERS["title_parp"] = title_parp
#----------------------------------------------------------------------------------------------------------------

### ���С���� ###
# 1. PD-L1
def pdl1_summary(pdl1):
	return "���ԡ�" if pdl1["result"] == "����" else "���ԣ�{0}Ϊ{1}��".format(pdl1["type"], pdl1["value"])
jinja2.filters.FILTERS["pdl1_summary"] = pdl1_summary

# 2. MSI
def msi_summary(msi):
	return "΢�����ȶ���MSS����" if msi["var_id"] == "MSS" else "΢���ǲ��ȶ���MSI-H����" if msi["var_id"] == "MSI-H" else "δ����MSI�����"
jinja2.filters.FILTERS["msi_summary"] = msi_summary

# 3. TMB
def tmb_summary(tmb):
	TMB_result = "��" if tmb["var_id"] == "TMB-L" else "��"
	return "{0} Muts/Mb, ����ͻ�为�ɽ�{1}��{2}����".format(tmb["TMB_value"], TMB_result, "TMB-H" if tmb["var_id"] == "TMB-H" else "TMB-L")
jinja2.filters.FILTERS["tmb_summary"] = tmb_summary

# 4. GEP
def gep_summary(gep):
	return "GEP��ֵΪ{0}��".format(gep["gep_score"])
jinja2.filters.FILTERS["gep_summary"] = gep_summary

# 5. TME
def tme_summary(tme):
	tme_dict = {
		"IE/F" : "���߸���/��ά������(IE/F)",
		"IE" : "���߸���/����ά������(IE)",
		"F" : "��ά������(F)",
		"D" : "���߻�Į��(D)"
	}
	return "TME����Ϊ{0}".format(tme_dict.get(tme["tme_type"], tme["tme_type"]))
jinja2.filters.FILTERS["tme_summary"] = tme_summary

# 6. ���߼������Ƽ���Ч��ػ���
def io_summary(io):
	result = []
	if io["io_p_summary"]:
		result.append(io["io_p_summary"]+"����Ч����أ�")
	if io["io_n_summary"]:
		result.append(io["io_n_summary"]+"����Ч����أ�")
	if not result:
		result = ["δ�����ر���"]
	rt = RichText()
	rt.add("��\n".join(result)+"��")
	return rt
jinja2.filters.FILTERS["io_summary"] = io_summary

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

# 7. HRD�����
def hrd_summary(info):
	gss = info[0]
	sample = info[1]
	result = []
	hrd_result = "HRD����" if gss["BRCA1"] or gss["BRCA2"] or float(gss["gss"]["gsscore"]) >= 50 else "HRD����"
	# �ӽ�������ο���һЩ���-2023.08.29
	# δ��⵽BRCA�²�/�����²�����+��������һ������
	# 1��baf_noise > 0.055
	# 2) depth_noise > 0.35
	# 3) ������ϸ������
	# 4��������ϸ���������Ǹ�ʽ����ȷ��
	# 5��������ϸ����������ʽ��ȷ��������ֵ<30
	note = "����������ο���" if not gss["BRCA1"] and not gss["BRCA2"] and \
							(float(gss["gss"]["baf_noise"]) > 0.055 or \
							 float(gss["gss"]["depth_noise"]) > 0.35 or \
							 not sample["tumor_content"] or \
							 (sample["tumor_content"] and not is_number(sample["tumor_content_num"])) or \
							 (sample["tumor_content"] and is_number(sample["tumor_content_num"]) and float(sample["tumor_content_num"]) < 30)) \
							else ""
	result.append(hrd_result + note)
	if gss["summary"]:
		result.append("HRRͨ·��ػ���ͻ�䣺{0}".format(gss["summary"]))
	rt = RichText()
	rt.add("��\n".join(result)+"��")
	return rt
jinja2.filters.FILTERS["hrd_summary"] = hrd_summary

# 8. GA�����
def ga_summary(info):
	ga = info[0]
	msi = info[1]
	result = []
	if ga["ebv_type"]["ebv_type"] == "P" and msi["var_id"] == "MSI-H":
		if ga["ebv_sum"]:
			result.append("EB������Ⱦ�ͣ�EBV����{0}".fromat(ga["ebv_sum"]))
		else:
			result.append("EB������Ⱦ�ͣ�EBV��")
		result.append("΢���ǲ��ȶ��ͣ�MSI��")
	elif ga["ebv_type"]["ebv_type"] == "P":
		result.append("EB������Ⱦ�ͣ�EBV��")
	elif msi["var_id"] == "MSI-H":
		result.append("΢���ǲ��ȶ��ͣ�MSI��")
	else:
		if ga["gs_sum"]:
			result.append("�������ȶ��ͣ�GS��")
		if ga["cin_sum"]:
			result.append("Ⱦɫ�岻�ȶ��ͣ�CIN��")
	if not result:
		result = ["-"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["ga_summary"] = ga_summary

# 9. �ӹ���Ĥ�����ӷ��ͽ��
def ec_summary(ec_type):
	ec_dict = {
		"POLE-ultramutated type EC" : "POLEͻ���ͣ�POLE mutation��POLE mut��",
		"MSI-H type EC" : "�����޸�����ȱ�ݣ�Mismatch repair deficiency��MMRd��",
		"CNH type EC" : "TP53����ͻ�䣨p53 abnormality��p53 abn��",
		"CNL type EC" : "�������Է����ף�Non-specific molecular profile��NSMP��"
	}
	return ec_dict.get(ec_type, ec_type)
jinja2.filters.FILTERS["ec_summary"] = ec_summary
	
# 10. ��ϵ���죨��Ҫչʾ������죩
def var_g_summary(var_list):
	g_var_list = var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"]
	result = []
	if g_var_list:
		for var in g_var_list:
			hgvs = var["hgvs_c"] if var["hgvs_p"] == "p.?" else var["hgvs_p"]
			clinic_g_cn = "�²��Ա���" if var["clinic_num_g"] == 5 else "�����²��Ա���"
			result.append("���{0} {1}��Ϊ{2}".format(var["gene_symbol"], hgvs, clinic_g_cn))
	else:
		result = ["�ڼ�ⷶΧ�ڣ�δ����²���/�����²��Ա���"]
	return "��".join(result)+"��"
jinja2.filters.FILTERS["var_g_summary"] = var_g_summary

# 10. ��ϸ��/��Դ��������
def var_s_summary(info):
	var_list = info[0]
	sample = info[1]
	result = []
	def sum_var(vlist):
		v_result = []
		for var in vlist:
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					v_result.append(var["gene_symbol"]+" "+var["hgvs_p"])
				else:
					v_result.append(var["gene_symbol"]+" "+var["hgvs_c"])
			elif var["bio_category"] == "Cnv":
				v_result.append(var["gene_symbol"]+" ����")
			elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					v_result.append("MET exon14 ��Ծ")
				else:
					if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" �ں�" not in v_result:
						v_result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" �ں�")
		return ", ".join(v_result)
	# �ж�������ʱ
	c_var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
				 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"])
	c_var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	c_var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"])
	c_var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	# �޶�������ʱ
	var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
				 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"] +\
				 var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"])
	var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
							var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_germline_nodrug"])
	var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
								  var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	if sample["control_sample_id"]:
		str_1 = "���{0}����ϸ�����죬���о����ٴ�����ı�����{1}��������������չ��ر�����{2}����".format(str(c_var_all_num), str(c_var_onco_drug_num), str(c_var_onco_nodrug_num))
		str_2 = "�����ٴ�����ı�����{0}��".format(c_var_onco_drug_str)
		if c_var_all_num != 0:
			if c_var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "�ڼ�ⷶΧ�ڣ�δ�����ϸ�����졣"
	else:
		str_1 = "���{0}�����죬���о����ٴ�����ı�����{1}��������������չ��ر�����{2}����".format(str(var_all_num), str(var_onco_drug_num), str(var_onco_nodrug_num))
		str_2 = "�����ٴ�����ı�����{0}��".format(var_onco_drug_str)
		if var_all_num != 0:
			if var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "�ڼ�ⷶΧ�ڣ�δ������졣"
jinja2.filters.FILTERS["var_s_summary"] = var_s_summary

#----------------------------------------------------------------------------------------------------------------

### HRD ###
# 1. BRCA�����
def hrd_brca_result(var_list):
	result = []
	if var_list:
		for var in var_list:
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_region"]+" "+var["hgvs_c"])
	else:
		result = ["δ����²��Ի������²��Ա���"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["hrd_brca_result"] = hrd_brca_result

### ���Ʋ���-��� ###
def get_evi_sum(evi_sum):
	result = []
	if evi_sum and "evi_split" in evi_sum.keys():
		if "Predictive_merge" in evi_sum["evi_split"].keys() and evi_sum["evi_split"]["Predictive_merge"]:
			result.extend([{"regimen" : a["regimen_name"], "inter" : a["evi_interpretation"]} for a in evi_sum["evi_split"]["Predictive_merge"]])
		if "Prognostic" in evi_sum["evi_split"].keys() and evi_sum["evi_split"]["Prognostic"]:
			result.extend([{"regimen" : "Ԥ�����", "inter" : a["evi_interpretation"]} for a in evi_sum["evi_split"]["Prognostic"]])
		if "Diagnostic" in evi_sum["evi_split"].keys() and evi_sum["evi_split"]["Diagnostic"]:
			result.extend([{"regimen" : "����������", "inter" : a["evi_interpretation"]} for a in evi_sum["evi_split"]["Diagnostic"]])
	if not result:
		result = [{"regimen" : "", "inter": "Ŀǰ���ڸñ�����ٴ�����ʵ���в���ȷ��"}]
	return result
jinja2.filters.FILTERS["evi_sum"] = get_evi_sum

### ������ ###
# 1. ���Ʒ�������-�����־��
def approval_regimen_biomarker(info):
	biomaker_list = info[0]
	judge_mergeMET = info[1]
	result = []
	for i in biomaker_list:
		if "hgvs_c" in i.keys() and i["hgvs_c"]:
			if i["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"], i["hgvs_p"]))
			else:
				result.append("{0} {1} {2}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"]))
		elif "cnv_type" in i.keys() and i["cnv_type"]:
			result.append("{0} ����".format(i["gene_symbol"]))
		elif "five_prime_gene" in i.keys() and i["five_prime_gene"]:
			if i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
				result.append("MET exon14 ��Ծ")
			else:
				# ���²��hgvs��CP40��region����ins
				if "hgvs" in i.keys() and i["hgvs"]:
					i["five_prime_region"] = "-".join(re.split("-", (re.split(":", i["hgvs"])[2]))[:-1]) \
						  					 if not re.search("--", i["hgvs"]) \
											 else re.split("_", (re.split("--", i["hgvs"])[0]))[-1]
					i["three_prime_region"] = re.split(":", i["hgvs"])[-1] \
											  if not re.search("--", i["hgvs"]) \
											  else re.split("_", (re.split("--", i["hgvs"])[1]))[-1]
					# ��һ������-2023.10.19
					# var_hgvs�¸�ʽ��gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, �ɵ�Ϊgene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds����xxx:exon1��xxx:exon2
					i["five_prime_region"] = re.split(":", i["five_prime_region"])[-1] if re.search(":", i["five_prime_region"]) else i["five_prime_region"]
					i["three_prime_region"] = re.split(":", i["three_prime_region"])[-1] if re.search(":", i["three_prime_region"]) else i["three_prime_region"]
					# �������-2023.10.19
				result.append("{0}:{1}:{2}-{3}:{4}:{5}".format(i["five_prime_gene"], i["five_prime_transcript"], i["five_prime_region"],\
															   i["three_prime_gene"], i["three_prime_transcript"], i["three_prime_region"]))
		elif "biomarker_type" in i.keys() and i["biomarker_type"]:
			if i["biomarker_type"] == "KRAS/NRAS/BRAF WT":
				result.append("KRAS/NRAS/BRAF Ұ����")
			elif i["biomarker_type"] == "HRD-":
				result.append("HRD����")
			elif i["biomarker_type"] == "HRD+":
				result.append("HRD����")
			else:
				result.append(i["biomarker_type"])
		else:
			result.append("�޷��ֱ�ķ��ӱ�־�")
	# ������MET 14��ԾDNA/RNA����Ļ�����ɾ��RNA��Ľ����������DNA
	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)
	#print (result_redup)
	#print (judge_mergeMET)
	if judge_mergeMET:
		if "MET exon14 ��Ծ" in result_redup:
			result_redup.remove("MET exon14 ��Ծ")
	if not result_redup:
		result_redup = ["-"]
	rt = RichText()
	rt.add("\n".join(result_redup))
	return rt
	#return result_redup
jinja2.filters.FILTERS["approval_regimen_biomarker"] = approval_regimen_biomarker

# 2. PARP������ܱ����ݲ�Ʒ��HRR��HRD��ѡ���Ӧչʾ�б�
def choose_parp_list(info):
	hrr_list = info[0]["cdx"]["format5_forHRR"]
	hrd_list = info[0]["cdx"]["format4_forHRDC"]
	sample = info[1]
	if sample["prod_names"] in ["HRR��ȫѪ��", "HRR����֯��", "HRR����֯ ȫѪ��"]:
		return hrr_list
	elif sample["prod_names"] in ["HRD Complete����֯��"]:
		return hrd_list
jinja2.filters.FILTERS["choose_parp_list"] = choose_parp_list

# ��ϸ���ȼ�ת��
def somatic_class_stran(clinic_num_s):
	if clinic_num_s == 5:
		return "I��"
	elif clinic_num_s == 4:
		return "II��"
jinja2.filters.FILTERS["somatic_class_stran"] = somatic_class_stran


# 5. �����ϸ���-�����
def detect_result(var):
	result = []
	if var["bio_category"] == "Snvindel":
		if var["hgvs_p"] != "p.?":
			result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
		else:
			result.append(var["gene_region"]+" "+var["hgvs_c"])
		result.append(var["transcript_primary"])
	elif var["bio_category"] == "Cnv":
		result.append("����")
	elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
		if var["five_prime_gene"] == var["three_prime_gene"] and var["three_prime_gene"] == "MET":
			result.append("MET exon14 ��Ծ")
			result.append(var["three_prime_transcript"])
		else:
			result.append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"�ں�")
			result.append(var["five_prime_transcript"]+"/"+var["three_prime_transcript"])
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["detect_result"] = detect_result

# 6. �����ϸ���-Ƶ��-��ϸ��
def detect_freq_s(var):
	freq = ""
	if var["bio_category"] == "Snvindel":
		freq = var["freq_str"]
	elif var["bio_category"] == "Cnv":
		freq = var["cn_mean"]
	elif var["bio_category"] == "Sv":
		freq = var["freq_str"]
	elif var["bio_category"] == "PSeqRnaSv":
		freq = str(var["freq"])+" copies"
	return freq
jinja2.filters.FILTERS["freq_s"] = detect_freq_s

# 7. �����ϸ���-Ƶ��-��ϵ
# �����Ϻ���Ŀ���ж������������
def detect_freq_g_SH(var):
	return "δ��ȡ�����ݣ�" if not var["freq_sc"] else "����" if float(var["freq_sc"]) >= 0.85 else "�Ӻ�"
jinja2.filters.FILTERS["freq_g_SH"] = detect_freq_g_SH

# ioչʾ����Ϊ����
def io_stran(io_list):
	if not io_list:
		io_list = ["-"]
	rt = RichText()
	rt.add("\n".join(io_list))
	return rt
jinja2.filters.FILTERS["io_stran"] = io_stran

# NCCNָ���Ƽ����������������չʾ-MP
def cdx_type1(info):
	cdx_list = info[0]
	sample = info[1]
	rt = RichText()
	result = []
	if cdx_list:
		for var in cdx_list:
			if var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					result.append("MET exon14 ��Ծ��{0}".format(var["freq"]))
				else:
					if "var_info_M" in var.keys() and var["var_info_M"]:
						result.append(var["var_info_M"]+"��"+var["freq"])
					else:
						result.append(var["var_info"]+"��"+var["freq"])
			else:
				if var["gene_symbol"] in ["BRCA1", "BRCA2"]:
					if sample["control_sample_id"]:
						if var["var_origin"] == "germline":
							result.append("{0}��{1}".format(var["var_info"], var["freq_rc_str"]))
							result.append("��ϵ����")
						else:
							result.append("{0}��{1}".format(var["var_info"], var["freq"]))
							result.append("��ϸ������")
					else:
						result.append("{0}��{1}".format(var["var_info"], var["freq"]))
				else:
					result.append("{0}��{1}".format(var["var_info"], var["freq"]))
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type1"] = cdx_type1

# NCCNָ���Ƽ����������������չʾ-116��������
def cdx_type2_var_info(var):
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				return var["hgvs_p"]
			else:
				return var["hgvs_c"]
		elif var["bio_category"] == "Cnv":
			return "����"
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				return "MET exon14 ��Ծ"
			else:
				return "{0}-{1} �ں�".format(var["five_prime_gene"], var["three_prime_gene"])
	else:
		return "δ��⵽"
jinja2.filters.FILTERS["cdx_type2_var_info"] = cdx_type2_var_info

# NCCNָ���Ƽ����������������չʾ-116���
def cdx_type2_freq(info):
	var = info[0]
	sample = info[1]
	rt = RichText()
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["var_origin"] == "germline":
				result.append("δ��ȡ��freq_rc��" if not var["freq_rc"] else "����" if float(var["freq_rc"]) >= 0.85 else "�Ӻ�")
				#rt.add("δ��ȡ��freq_rc��" if not var["freq_rc"] else "����" if float(var["freq_rc"]) >= 0.85 else "�Ӻ�", \
				#	   size = 17, color = "#000000", font = "Source Han Sans Normal")
			else:
				result.append(var["freq_str"])
				#rt.add(var["freq_str"], size = 17, color = "#000000", font = "Source Han Sans Normal")
			if var["gene_symbol"] in ["BRCA1", "BRCA2"] and sample["control_sample_id"]:
				result.append("��ϵ����" if var["var_origin"] == "germline" else "��ϸ������")
				#rt.add("\n��ϵ����" if var["var_origin"] == "germline" else "\n��ϸ������", size = 17, color = "#000000", font = "Source Han Sans Normal")
		elif var["bio_category"] == "Cnv":
			result.append(var["cn_mean"])
			#rt.add(var["cn_mean"], size = 17, color = "#000000", font = "Source Han Sans Normal")
		elif var["bio_category"] == "Sv":
			result.append(var["freq_str"])
			#rt.add(var["freq_str"], size = 17, color = "#000000", font = "Source Han Sans Normal")
	else:
		result.append("-")
		#rt.add("-", size = 17, color = "#000000", font = "Source Han Sans Normal")
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type2_freq"] = cdx_type2_freq

# NCCNָ���Ƽ����������������չʾ-116�ٴ�����
def cdx_type2_class(var):
	stran = {
		"germline" : {
			5 : "�²���",
			4 : "�����²���"
		},
		"somatic" : {
			5 : "I��",
			4 : "II��",
			3 : "III��"
		}
	}
	if "bio_category" in var.keys():
		if var["var_origin"] == "germline":
			return stran["germline"].get(var["clinic_num_g"])
		else:
			return stran["somatic"].get(var["clinic_num_s"])
	else:
		return "-"
jinja2.filters.FILTERS["cdx_type2_class"] = cdx_type2_class

# NCCNָ���Ƽ���������-CP40
def cdx_type3_var_info(var):
	result = []
	if "var_info" in var.keys() and var["var_info"]:
		if "MET-MET" in var["var_info"]:
			result.append("MET exon14 skipping")
		elif var["merge_sv_list"]:
			for a in var["merge_sv_list"]:
				result.append(a+"�ں�")
		else:
			result.append(var["var_info"])
		if var["note"]:
			result.append("(MET exon14 skipping)")
	else:
		result = ["δ��⵽"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type3_var_info"] = cdx_type3_var_info

# NCCNָ���Ƽ���������-���
def cdx_type3_freq(var):
	result = []
	if "freq" in var.keys() and var["freq"]:
		if "�ں�" in var["var_info"]:
			result.append(str(var["freq"])+" copies")
		else:
			result.append(var["freq"])
		if var["note_freq"]:
			result.append("("+var["note_freq"]+")")
	else:
		result.append("-")
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type3_freq"] = cdx_type3_freq

# NCCNָ���Ƽ���������-�ٴ�����
def cdx_type3_class(var):
	if "level" in var.keys() and var["level"]:
		if var["level"] == "5":
			return "I��"
		elif var["level"] == "4":
			return "II��"
		else:
			return "III��"
	else:
		return "-"
jinja2.filters.FILTERS["cdx_type3_class"] = cdx_type3_class

# ѡ��CRC25��TC21��GA18���������
def choose_116_list(info):
	CRC25_list = info[0]["CRC25"]
	GA18_list = info[0]["GA18"]
	TC21_list = info[0]["TC21"]
	prod_names = info[1]
	if prod_names in ["CRC25����֯��", "CRC25��ѪҺ��"]:
		return CRC25_list
	elif prod_names in ["GA18����֯��", "GA18��ѪҺ��"]:
		return GA18_list
	elif prod_names in ["TC21����֯��", "TC21��ѪҺ��"]:
		return TC21_list
jinja2.filters.FILTERS["choose_116_list"] = choose_116_list

### ��Ʒ���� ###
# 1. ���
def product_state_index(info):
	var_brca = info[0]
	sample = info[1]
	if sample["prod_names"] in ["BRCA1/BRCA2��ȫѪ��",  "BRCA1/2�������ӣ�"]:
		if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"] or\
																	   var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
			if re.search("����", sample["tumor_names_cn"]):
				return "5"
			else:
				return "7"
		else:
			if re.search("����", sample["tumor_names_cn"]):
				return "4"
			else:
				return "6"
	elif sample["prod_names"] in ["BRCA1/BRCA2����֯��"]:
		if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"]:
			return "7"
		else:
			return "6"
	elif sample["prod_names"] in ["BRCA1/BRCA2����֯ ȫѪ��"]:
		if var_brca["snv_m"]["B1_G_L5"] or var_brca["snv_m"]["B2_G_L5"] or var_brca["snv_m"]["B1_G_L4"] or var_brca["snv_m"]["B2_G_L4"] or\
																	   var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
			return "7"
		else:
			return "6"
	elif sample["prod_names"] in ["�����ۺ���"]:
		return "4"
	# ��һ��HRRȫѪ��������4������6
	elif sample["prod_names"] in ["HRR��ȫѪ��"]:
		if re.search("����", sample["tumor_names_cn"]):
			return "4"
		else:
			return "6"
	else:
		return "6"
jinja2.filters.FILTERS["product_state_index"] = product_state_index

### sample ###
# 1. ��������
# 116ѪҺ��Ŀ���������ͻὫ�����Ͷ��յ�дһ��
def sample_type(sample):
	result = ""
	if sample["prod_names"] in ["Pan116��ѪҺ��", "TC21��ѪҺ��", "GA18��ѪҺ��", "LC76��ѪҺ��", "CRC25��ѪҺ��", "Master Panel��ѪҺ��"] and sample["control_sample_id"]:
		if sample["sample_type"] and sample["control_sample_type"]:
			if sample["sample_type"] != sample["control_sample_type"]:
				result = sample["sample_type"]+","+sample["control_sample_type"]
			else:
				result = sample["sample_type"]
		elif sample["sample_type"] and not sample["control_sample_type"]:
			result = sample["sample_type"]
		elif not sample["sample_type"] and sample["control_sample_type"]:
			result = sample["control_sample_type"]
	else:
		result = sample["sample_type"]
	if not result:
		result = "-"
	return result
jinja2.filters.FILTERS["sample_type"] = sample_type

# 2. ��������
# 116ѪҺ��Ŀ�����������Ὣ�����Ͷ��յ�дһ��
def sample_amount(sample):
	result = ""
	if sample["prod_names"] in ["Pan116��ѪҺ��", "TC21��ѪҺ��", "GA18��ѪҺ��", "LC76��ѪҺ��", "CRC25��ѪҺ��", "Master Panel��ѪҺ��"] and sample["control_sample_id"]:
		if sample["sample_type"] and sample["control_sample_type"]:
			if sample["sample_type"] != sample["control_sample_type"]:
				result = sample["sample_amount"]+","+sample["control_sample_amount"]
			else:
				result = sample["sample_amount"]
		elif sample["sample_type"] and not sample["control_sample_type"]:
			result = sample["sample_amount"]
		elif not sample["sample_type"] and sample["control_sample_type"]:
			result = sample["control_sample_amount"]
	else:
		result = sample["sample_amount"]
	if not result:
		result = "-"
	return result
jinja2.filters.FILTERS["sample_amount"] = sample_amount

# 3. �����ɼ�����
# ѪҺ��Ŀ�ɼ����ڴ�blood_collection_date����ȡ��������gather_data����ȡ
def gather_data(sample):
	result = ""
	if sample["prod_names"] in ["BRCA1/BRCA2��ȫѪ��", "HRR��ȫѪ��", "10����ѪҺ��", "61�Ŵ�����", "Pan116��ѪҺ��", "�����ۺ���",\
							    "TC21��ѪҺ��", "GA18��ѪҺ��", "LC76��ѪҺ��", "CRC25��ѪҺ��", "BPTM��ȫѪ��", "Master Panel��ѪҺ��", "BRCA1/2�������ӣ�"]:
		result = sample["blood_collection_date"]
	else:
		result = sample["gather_data"]
	if not result:
		result = "-"
	return result
jinja2.filters.FILTERS["gather_data"] = gather_data


### ���Ʒ�������-���� ###
# MP��֯�����ѳ��������ٰ���ǰ���ٰ�ʱչʾHRD��������������������Ʒ������֣���Ҫ���˵�HRD+��HRD-������
def approval_regimen_filter_hrd(info):
	regimen_list = info[0]
	sample = info[1]
	hrd_p = {"biomarker_type" : "HRD+"}
	hrd_n = {"biomarker_type" : "HRD-"}
	result = []
	if "Master" in sample["prod_names"] and not set(["���ٰ�", "�ѳ���", "ǰ���ٰ�"])&set(sample["tumor_list"]) and "ʵ����" not in sample["tumor_names_cn"]:
		for regimen in regimen_list:
			for a in regimen["var"]:
				if a == hrd_p or a == hrd_n:
					regimen["var"].remove(a)
			if regimen["var"]:
				result.append(regimen)
	else:
		result = regimen_list
	return result
jinja2.filters.FILTERS["approval_regimen_filter_hrd"] = approval_regimen_filter_hrd


### QC ###
# 1. HRRȫѪ
def qc_gHRR(qc, lib_quality_control):
	ngs_qc = qc["dna_data_qc"] if "dna_data_qc" in qc.keys() else {}
	lib_qc = lib_quality_control["lib_dna_qc"] if "lib_dna_qc" in lib_quality_control.keys() else {}
	if float(ngs_qc["cleandata_q30_num"]) >= 0.75 and float(ngs_qc["depth_ssbc_num"]) >= 100:
		if lib_qc and "dna_qty" in lib_qc.keys() and lib_qc["dna_qty"] and "library_qty" in lib_qc.keys() and lib_qc["library_qty"]:
			if float(lib_qc["dna_qty"]) >= 20 and float(lib_qc["library_qty"]) >= 200:
				return "�ϸ�"
			else:
				return "����"
		else:
			return "�ϸ��ʿ�����ȱʧ���벹�����ݺ�����������"
	else:
		return "����"
jinja2.filters.FILTERS["qc_gHRR"] = qc_gHRR

# 2. HRR��֯

qc_stand = {
	"gHRR" : {},
	"tHRR" : {
		"qc" : {
			"dna_data_qc" : {
				"cleandata_q30_num" : 0.75,
				"depth_ssbc_num" : 300
			}
		},
		"lib_quality_control" : {
			"lib_dna_qc" : {
				"dna_qty" : 30,
				"library_qty" : 200
			}
		},
		"tumor_content" : 20
	}
}

def qc_tHRR(qc, lib_quality_control, sample, prod):
	stand = qc_stand.get(prod)
	# qc���
	qc_result = []
	if "dna_data_qc" in qc.keys() and qc["dna_data_qc"]:
		for item in stand["qc"]["dna_data_qc"].keys():
			if item in qc["dna_data_qc"].keys() and float(qc["dna_data_qc"][item]) >= stand["dna_data_qc"][item]:
				qc_result.append("�ϸ�")
			else:
				qc_result.append("����")
	# ʪʵ����
	if "lib_dna_qc" in lib_quality_control.keys() and lib_quality_control["lib_dna_qc"]:
		for item in stand["lib_quality_control"]["dna_data_qc"].keys():
			pass

# v4�ټ�ͨ��ģ��ʹ��
# MPС��չʾ�������-��ϸ������
def var_sum_s(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			result.append(var["gene_symbol"]+" ����")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			# �ںϿ��ܻ����ظ���rna exon��ͬ���ϵ㲻ͬ�������
			if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" �ں�" not in result:
				result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" �ں�")
	return ", ".join(result)
jinja2.filters.FILTERS["var_sum_s_filter"] = var_sum_s

# MPС��չʾ�������-��ϵ����
def var_sum_g(var_list):
	result = []
	for var in var_list:
		hgvs = var["hgvs_c"] if var["hgvs_p"] == "p.?" else var["hgvs_p"]
		clinic_g_cn = "�²�����" if var["clinic_num_g"] == 5 else "�����²�����"
		result.append("���{0} {1}��Ϊ{2}".format(var["gene_symbol"], hgvs, clinic_g_cn))
	return "��".join(result)
jinja2.filters.FILTERS["var_sum_g_filter"] = var_sum_g


######################################## ����ģ�������#################################################
# ��������CP40���������Ҫ�г����о�֤�ݡ��ο�����
# �������£�
# 1. A��C3�ȼ������ƥ�䵽NCCN�ģ�չʾ��NCCN�ٴ�ʵ��ָ�ϡ�
# 2. ��PMID��չʾPMID�ţ�û�е�չʾNCT�ţ���û�е�չʾб�ܣ����е�չʾPMID��
# 3. NCCN�ͣ�PMID or NCT�����ܻ�ͬʱ���֣�NCCNҪչʾ����PMID or NCT��������1��
def get_refer_CQXN(evi_sum):
	# ��ȡPMID
	def getPMID_from_inter(inter):
		pmid_list = []
		mat = re.compile(r"PMID.\s?\d+")
		for i in mat.findall(str(inter)):
			if re.search(":|: |��|�� ", i):
				pmid = (re.split(":|: |��|�� ", i))[1].replace(" ", "")
			else:
				pmid = (re.split("PMID", i))[1]
			pmid_list.append(pmid)
		return pmid_list
	# ��ȡNCT
	def getNCT_from_inter(inter):
		nct_list = []
		mat = re.compile(r"NCT\d+")
		for i in mat.findall(str(inter)):
			nct_list.append(i)
		return nct_list
	# �����о�֤�ݣ��������ơ�������Ϻ�Ԥ��
	evi_list = []
	for i in ["Predictive", "Prognostic", "Diagnostic"]:
		if i in evi_sum["evi_split"].keys() and evi_sum["evi_split"][i]:
			evi_list.extend(evi_sum["evi_split"][i])
	# ��ȡ�ο�����
	refer = []
	for evi in evi_list:
		# �ж��Ƿ���nccn
		if "A" in evi["evi_conclusion"] or "C3" in evi["evi_conclusion"]:
			if re.search("NCCN", evi["evi_interpretation"]):
				refer.append("NCCN�ٴ�ʵ��ָ��")
		# �ж�pmid��NCT
		pmid_list = getPMID_from_inter(evi["evi_interpretation"])
		nct_list = getNCT_from_inter(evi["evi_interpretation"])
		if pmid_list:
			for i in pmid_list:
				refer.append("PMID: "+str(i))
		else:
			if nct_list:
				refer.extend(nct_list)
	# ȥ��
	refer_redup = []
	for i in refer:
		if i not in refer_redup:
			refer_redup.append(i)
	if not refer_redup:
		refer_redup = ["/"]
	rt = RichText()
	rt.add("\n".join(refer_redup))
	return rt
jinja2.filters.FILTERS["refer_CQXN"] = get_refer_CQXN

def ZDY_HRR_germline_summary(var_list):
	var_dict = {}
	for var in var_list:
		if var["type"] == "Loss":
			var_info = "{0} {1} del".format(var["gene_symbol"], var["value"])
			var["clinic_num_g"] = 4
		else:
			if var["hgvs_p"] != "p.?":
				var_info = "{0} {1}:{2}:{3}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"], var["hgvs_p"])
			else:
				var_info = "{0} {1}:{2}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"])
		if var["clinic_num_g"] not in var_dict.keys():
			var_dict.setdefault(var["clinic_num_g"], [])
		var_dict[var["clinic_num_g"]].append(var_info)
	result = []
	if 5 in var_dict.keys():
		result.append("���{0}���²��Ա��죬{1}".format(str(len(var_dict[5])), ", ".join(var_dict[5])))
	if 4 in var_dict.keys():
		result.append("���{0}�������²��Ա��죬{1}".format(str(len(var_dict[4])), ", ".join(var_dict[4])))
	if not var_dict:
		result.append("δ����²��������²�����ϵ����")
	return "��".join(result)
jinja2.filters.FILTERS["ZDY_HRR_germline_summary"] = ZDY_HRR_germline_summary

# ��ɽ����HRDC�����С��-����-2023.07.12
def summary_ZSRM_hrd_v2(info):
	# ���ԣ����ڱ����ͼ����������HRD״̬������/���ԡ���BRCA gene_region hgvs_c hgvs_p I����졢TP53 gene_region hgvs_c hgvs_p��
	# ���ԣ����ڱ����ͼ����������HRD������/���ԡ���δ���I�ࡢII����졣
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I��", 4 : "II��"}

	var_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["hgvs_p"] != "p.?":
			var_result.append("{0} {1} {2} {3} {4}����".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			var_result.append("{0} {1} {2} {3}����".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	
	hrd_result = "HRD״̬����" if hrd["var_id"] == "HRD+" else "HRD״̬����"

	if var_result:
		return "���ڱ����ͼ����������"+hrd_result+"��"+"��".join(var_result)+"��"
	else:
		return "���ڱ����ͼ����������"+hrd_result+"��δ���I�ࡢII����졣"
jinja2.filters.FILTERS["summary_ZSRM_hrd_v2"] = summary_ZSRM_hrd_v2

# ��������-�������ɾ�����Ϊ3�ķǱ��������죨���ӱ�����⣩
def filter_III_var_CQXN(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] in ["Sv", "Cnv"]:
			result.append(var)
		elif var["bio_category"] == "Snvindel":
			if var["clinic_num_s"] == 3 and var["type"] in ["Intronic", "3'UTR", "5'UTR", "FlankingRegion3", "FlankingRegion5"]:
				pass
			else:
				result.append(var)
	return result
jinja2.filters.FILTERS["filter_III_var_CQXN"] = filter_III_var_CQXN

# ���ڶ�ԺHRR�����ܻ����ٴ������չʾBRCA����-2023.07.17
def filter_clinicalTrial_SZEY(clinical_list):
	return [a for a in clinical_list if a["gene_symbol"] in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["filter_clinicalTrial_SZEY"] = filter_clinicalTrial_SZEY

# ���ڶ�ԺHRR��ҩ����ܽ�չʾBRCA������ص�-2023.07.17
def filter_drug_SZEY(drug_list):
	for drug in drug_list:
		drug["var_filter"] = [a for a in drug["var"] if ("gene_symbol" in a.keys() and a["gene_symbol"] and a["gene_symbol"] in ["BRCA1", "BRCA2"]) or\
							  "biomarker_type" in a.keys() and "BRCA1" in a["biomarker_type"] or \
							  "biomarker_type" in a.keys() and "BRCA2" in a["biomarker_type"]]
	return [a for a in drug_list if a["var_filter"]]
jinja2.filters.FILTERS["filter_drug_SZEY"] = filter_drug_SZEY

# �ӱ�ʡ����CP40��I������չʾABҩ�II������չʾCDҩ��
# info��ʽΪ[msi, knb, var_level_I, var_level_II, regimen_approval_list]
def filter_drug_HNRM(info):
	# ��������б�����ָ���б�ָ���ȼ�ҩ��Ľ������ʽΪ[("��������", var1), ("��������", var2), ("��������", var1)]
	def get_list(var_list, level_list):
		result = []
		for var in var_list:
			bio = []
			if "bio_category" in var.keys():
				if var["bio_category"] == "Snvindel":
					bio.append({
						"bio_category" : var["bio_category"],
						"gene_symbol" : var["gene_symbol"],
						"hgvs_c" : var["hgvs_c"],
						"hgvs_p" : var["hgvs_p"]
					})
				elif var["bio_category"] == "Cnv":
					bio.append({"bio_category" : var["bio_category"], "gene_symbol" : var["gene_symbol"]})
				elif var["bio_category"] == "Sv":
					if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
						bio.append({"bio_category" : var["bio_category"], "hgvs" : "MET exon14 skipping"})
					else:
						for i in var["merge_sv_list"]:
							bio.append({"bio_category" : var["bio_category"], "hgvs" : i+"�ں�"})
			elif "var_id" in var.keys():
				if var["var_id"] == "KRAS/NRAS/BRAF WT":
					bio.append({"var_id" : "KRAS/NRAS/BRAF Ұ����"})
				elif var["var_id"] == "MSI-H":
					bio.append({"var_id" : "MSI-H"})
			if "evi_sum" in var.keys() and var["evi_sum"] and "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				for evi in var["evi_sum"]["evi_split"]["Predictive"]:
					if evi["evi_conclusion_simple"] in level_list:
						result.append((evi["regimen_name"], bio))
		return result
	
	# 1. ����MSI-H AB��KNB AB��I��AB��II��CDҩ��
	regimen_list_filter = []
	regimen_list_filter.extend(get_list([info[0]], ["A", "B"]))
	regimen_list_filter.extend(get_list([info[1]], ["A", "B"]))
	regimen_list_filter.extend(get_list(info[2], ["A", "B"]))
	regimen_list_filter.extend(get_list(info[3], ["C", "D"]))

	# 2. ��ʽת��{"��������" : [var1, var2], "��������" : [var1]}
	regimen_dict_filter = {}
	for i in regimen_list_filter:
		regimen_dict_filter.setdefault(i[0], [])
		for j in i[1]:
			regimen_dict_filter[i[0]].append(j)

	# �����Ʒ��������б���й��ˣ����������˺�ı����б�var_hnrm
	regimen_approval_list = info[4]
	regimen_approval_filter = []
	for regimen in regimen_approval_list:
		regimen_name = regimen["regimen_cn"] if regimen["regimen_cn"] else regimen["regimen_en"]
		if regimen_name in regimen_dict_filter.keys():
			regimen["var_hnrm"] = regimen_dict_filter[regimen_name]
			regimen_approval_filter.append(regimen)

	return regimen_approval_filter
jinja2.filters.FILTERS["filter_drug_HNRM"] = filter_drug_HNRM

# ��������ΰ���Ҫ�ж�EGFR���������
def judge_EGFR(info):
	var_list = info[0]
	sample = info[1]
	gene_list = [var["gene_symbol"] for var in var_list]
	if "EGFR" not in gene_list and "�ΰ�" in sample["tumor_list"] and "ʵ����" not in sample["tumor_names_cn"]:
		return "EGFRnone"
jinja2.filters.FILTERS["judge_EGFR"] = judge_EGFR

# ��������ΰ���Ҫ�ж�EGFR T790Mͻ�����
def sort_var_forhz_egfr(info):
	var_list = info[0]
	sample = info[1]
	for var in var_list:
		var["egfr_sort"] = 0
	t790m = {"gene_symbol" : "EGFR", "t790m_var_info" : "δ���T790M����", "egfr_sort" : 1}
	egfr_hgvs_p_list = [var["hgvs_p"].replace("p.","").replace("(","").replace(")","") for var in var_list if var["gene_symbol"] == "EGFR" and var["bio_category"] == "Snvindel"]
	
	for var in var_list:
		if var["bio_category"] == "Cnv" and var["gene_symbol"] == "EGFR":
			egfr_hgvs_p_list.append("����")

	if "T790M" not in egfr_hgvs_p_list and egfr_hgvs_p_list and "�ΰ�" in sample["tumor_list"] and "ʵ����" not in sample["tumor_names_cn"]:
		# ����EGFR���һ����������index
		egfr_last_index = 0
		for var in var_list:
			if var["gene_symbol"] == "EGFR":
				egfr_last_index = var_list.index(var)
		# ����T790M�����
		var_list.insert(egfr_last_index+1, t790m)

	return var_list
jinja2.filters.FILTERS["sort_var_forhz_egfr"] = sort_var_forhz_egfr

# ������ɽCP40-С�Ჿ��չʾ����λ��������ʽ�ϻ�������Ҫ�ϲ���Ԫ��-2023.07.31
def chemo_stran_fdzs(chemo_list):
	result = {}
	for i in chemo_list:
		result.setdefault(i["gene_symbol"], [])
		result[i["gene_symbol"]].append({
			"dbsnp" : i["dbsnp"],
			"genotype" : i["genotype"]
		})
	chemo_result = []
	for k, v in result.items():
		v_sort = sorted(v, key=lambda i:i["dbsnp"])
		chemo_result.append({
			"gene_symbol" : k,
			"info" : v_sort
		})
	return sorted(chemo_result, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["chemo_stran_fdzs"] = chemo_stran_fdzs

# ������һ�����ο�����-2023.07.31
def refer_nccn_fjfy(sample):
	refer = {
		"��Сϸ���ΰ�" : [
			"�й���Сϸ���ΰ�����EGFR+T790M����ͻ����ר�ҹ�ʶ��2018��棩",
			"������������NSCLC�е��ٴ�Ӧ���й�ר�ҹ�ʶ��2020�棩",
			"��Сϸ���ΰ����Ӳ������ٴ�ʵ��ָ�ϣ�2021�棩",
			"�й���Сϸ���ΰ�RET�����ں��ٴ����ר�ҹ�ʶ��2021��棩",
			"��Сϸ���ΰ�MET�ٴ�����й�ר�ҹ�ʶ��2022��棩",
			"��Сϸ���ΰ�ϸ�봩��ϸ��ѧ�걾������ר�ҹ�ʶ��2022��棩",
			"��Сϸ���ΰ����Խ�Ĥǻ��Һ���Ӳ������й�ר�ҹ�ʶ��2022��棩",
			"��Сϸ���ΰ��ںϻ������ٴ�ʵ���й�ר�ҹ�ʶ��2023��棩",
			"�й����ڷ�Сϸ���ΰ�BRAFͻ������ר�ҹ�ʶ��2023��棩"
		],
		"��ֱ����" : [
			"��ֱ�������������ʵ����΢���ǲ��ȶ��Լ���й�ר�ҹ�ʶ��2019��棩",
			"��ֱ�������ӱ�־���ٴ�����й�ר�ҹ�ʶ��2021��棩",
			"��ֱ�������Ӽ���ͨ�������й�ר�ҹ�ʶ��2021��棩"
		],
		"�ѳ���" : [
			"�ѳ���Ƥ�԰�BRCA��������й�ר�����ۣ�2017��棩",
			"������һ����������BRCA1_2������ָ��(2019���)",
			"��Ƥ���ѳ���PARP���Ƽ���������־������й�ר�ҹ�ʶ��2020��棩",
			"BRCA1_2���ݽ���й�ר�ҹ�ʶ(2021���)"
		],
		"ǰ���ٰ�" : [
			"�й�ǰ���ٰ����߻�����ר�ҹ�ʶ��2020��棩",
			"ǰ���ٰ�ͬԴ�����޸������⼰������ר�ҹ�ʶ��2022��棩"
		],
		"���ٰ�" : [
			"�й����ٰ�����BRCA���������ٴ�Ӧ��ר�ҹ�ʶ��2018��棩",
			"������һ����������BRCA1_2������ָ��(2019���)",
			"����/ת�������ٰ���־���ٴ�Ӧ��ר�ҹ�ʶ��2019��棩",
			"�������ٰ��������ȵ������й�ר�ҹ�ʶ��2021��棩",
			"���ڰб�ָ�����ٰ���׼���Ʊ�־���ٴ�Ӧ��ר�ҹ�ʶ��2022��棩"
		],
		"θ��" : [
			"θ��HER2���ָ�ϣ�2016��棩",
			"HER2��������θ�����Ӱ������Ƶ��й�ר�ҹ�ʶ��2016��棩",
			"θ����ͨ�������ٴ�Ӧ���й�ר�ҹ�ʶ��2022��棩"
		]
	}
	result = []
	if "ʵ����" not in sample["tumor_names_cn"]:
		set_list = set(sample["tumor_list"]) & set([i for i in refer.keys()])
		for tumor in set_list:
			result.extend(refer[tumor])
	return result
jinja2.filters.FILTERS["refer_nccn_fjfy"] = refer_nccn_fjfy

# ɾ�������е�BCL2L11����-2023.08.02
def remove_bcl2l11(var_list):
	for var in var_list:
		if var["gene_symbol"] == "BCL2L11" and "hgvs_c" in var.keys() and var["hgvs_c"] and var["hgvs_c"] == "c.394+1479_394+4381del":
			var_list.remove(var)
	return var_list
jinja2.filters.FILTERS["remove_bcl2l11"] = remove_bcl2l11

# ��������һLC10/PAN116/CP40��ɾ��KRAS��NRAS��BRAF��KNB�к��С�CSCO����֤��
# ���迼��ɾ��֤�ݺ�Ա���ȼ���Ӱ��-2023.08.15
def filter_csco(evi_list):
	return [i for i in evi_list if "CSCO" not in i["evi_interpretation"]]
jinja2.filters.FILTERS["filter_csco"] = filter_csco

# ����ʡ��BRCA-���С�Ჿ�֣���Ҫ����б�壬���ұ���֮���á����-2023.08.28
# ��ÿ�����죨�������һ��������Ӹ��ٺţ�������ģ����ʹ��forѭ��չʾ������
def ahsl_brca_sum(info):
	var_list = info[0]
	judge_intumor = info[1]
	var_g = {5 : "�²��Ա���", 4 : "�����²��Ա���"}
	for var in var_list:
		if var["type"] == "Loss":
				var["clinic_num_g_stran"] = "�����²��Ա���"
		else:
			if var["var_origin"] == "germline":
				var["clinic_num_g_stran"] = var_g.get(var["clinic_num_g"])
			else:
				if var["clinic_num_g"] in [4, 5]:
					var["clinic_num_s_stran"] = "I�����" if judge_intumor == "intumor" else "II�����"
				else:
					var["clinic_num_s_stran"] = "III�����"
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			if var["type"] == "Loss" or var["var_origin"] == "germline":
				var["clinic_num_g_stran"] = var["clinic_num_g_stran"]+"��"
			else:
				var["clinic_num_s_stran"] = var["clinic_num_s_stran"]+"��"
	return var_list
jinja2.filters.FILTERS["ahsl_brca_sum"] = ahsl_brca_sum

# �Ϸ�ҽԺCP40�� NCCNָ���Ƽ���־���������ܱ���Ҫ���Ϊһ/��/��������չʾ-2023.09.01
# �ù��������Է���ָ���ȼ����
def split_cdx_nfyy(info):
	var_list = [var for var in info[0] if "var_info" in var.keys()]
	level = info[1]
	return [var for var in var_list if var["level"] == str(int(level))]
jinja2.filters.FILTERS["split_cdx_nfyy"] = split_cdx_nfyy

# �Ϻ��ʼ�CP40��������Ҫչʾ����/Ұ����-2023.09.18
# wt ΪҰ���ͣ�ֻ��һ��ΪrefΪ�Ӻϱ����ͣ���������ΪrefΪ���ϱ�����
def judge_chemo_genotype(info):
	dbsnp = info[0]
	alt = info[1]
	chemo_dict = {
		"rs3918290" : {"ref": "C", "wt": "C/C"},
		"rs55886062" : {"ref": "A", "wt": "A/A"},
		"rs67376798" : {"ref": "T", "wt": "T/T"},
		"rs75017182" : {"ref": "G", "wt": "G/G"},
		"rs10929302" : {"ref": "G", "wt": "G/G"},
		"rs4148323" : {"ref": "G", "wt": "G/G"},
		"rs8175347" : {"ref": "(TA)6", "wt": "(TA)6/(TA)6"}
		}
	split_alt = re.split("/", alt)
	if dbsnp in chemo_dict.keys():
		if split_alt[0] == chemo_dict.get(dbsnp)["ref"] and split_alt[1] == chemo_dict.get(dbsnp)["ref"]:
			return "Ұ����"
		else:
			if len(set(split_alt)) == 1:
				return "���ϱ�����"
			else:
				return "�Ӻϱ�����"
		#elif split_alt[0] != chemo_dict.get(dbsnp)["ref"] and split_alt[1] != chemo_dict.get(dbsnp)["ref"]:
		#	return "���ϱ�����"
		#else:
		#	return "�Ӻϱ�����"
	else:
		return "λ�㲻��������"
jinja2.filters.FILTERS["judge_chemo_genotype"] = judge_chemo_genotype

# ������116���ӹ���Ĥ��ʱ�����С��ͱ�����Ҫ��POLE��TP53������ΪI�࣬�ұ�������һ��Ҫ�����֤�ݡ�
def syx_ec_var_class(info):
	var_list = info[0]
	class_num = info[1]
	level_I = [var for var in var_list if var["clinic_num_s"] == 5 or (var["clinic_num_s"] == 4 and var["gene_symbol"] in ["POLE", "TP53"])]
	level_II = [var for var in var_list if var["clinic_num_s"] == 4 and var["gene_symbol"] not in ["POLE", "TP53"]]
	if class_num == 1:
		return level_I
	elif class_num == 2:
		return level_II
	else:
		return []
jinja2.filters.FILTERS["syx_ec_var_class"] = syx_ec_var_class

# ��ɽ����HRDC�����С��-����-2023.10.30
def summary_ZSRM_hrd_v3(info):
	# ���ԣ����ڱ����ͼ����������HRD״̬������/���ԡ���BRCA gene_region hgvs_c hgvs_p I����졢TP53 gene_region hgvs_c hgvs_p��
	# ���ԣ����ڱ����ͼ����������HRD������/���ԡ���δ���I�ࡢII����졣
	# 2023.10.30�������ݣ���������v1.3.0�汾��hrd�����gss�ֶ��е������
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I��", 4 : "II��"}
	gss = info[2]

	var_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["hgvs_p"] != "p.?":
			var_result.append("{0} {1} {2} {3} {4}����".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			var_result.append("{0} {1} {2} {3}����".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	
	#hrd_result = "HRD״̬����" if hrd["var_id"] == "HRD+" else "HRD״̬����"
	if hrd and "var_id" in hrd.keys():
		hrd_result = "HRD״̬����" if hrd["var_id"] == "HRD+" else "HRD״̬����"
	elif gss and "var_id" in gss.keys():
		hrd_result = "HRD״̬����" if gss["var_id"] == "HRD+" else "HRD״̬����"
	else:
		hrd_result = "δ��ȡ��HRD�����"

	if var_result:
		return "���ڱ����ͼ����������"+hrd_result+"��"+"��".join(var_result)+"��"
	else:
		return "���ڱ����ͼ����������"+hrd_result+"��δ���I�ࡢII����졣"
jinja2.filters.FILTERS["summary_ZSRM_hrd_v3"] = summary_ZSRM_hrd_v3

# ��ȡ�����б��е�Snvindel-2023.11.16
def filter_snvindel(var_list):
	return [var for var in var_list if var["bio_category"] == "Snvindel"]
jinja2.filters.FILTERS["filter_snvindel"] = filter_snvindel

# ��ȡ�����б��е�RNA SV-2023.11.16
def filter_sv(var_list):
	return [var for var in var_list if var["bio_category"] in ["Sv", "PSeqRnaSv"]]
jinja2.filters.FILTERS["filter_sv"] = filter_sv

# ��ȡ�����б��е�EGFR -2023.11.22
def filter_egfr(var_list):
	return [var for var in var_list if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR"]
jinja2.filters.FILTERS["filter_egfr"] = filter_egfr

#  XW5301-���ܻ�������������δ��⵽����Ļ���-2023.11.27
def getSum_forXW5301(var_list, gene_rule, gene_transcript):
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_rule:
				var["transcript_primary_xw5031"] = gene_transcript.get(gene)
				result.append(var)
	for gene in set(gene_rule) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"transcript_primary_xw5031" : gene_transcript.get(gene)
		})
	return sorted(result, key = lambda i:gene_rule.index(i["gene_symbol"]))

# XW5301-snvindel����-2023.11.27
def sort_for5301_snvindel(var_list):
	gene_rule = ["MET", "EGFR", "ALK", "KRAS", "ROS1", "RET", "ERBB2", "BRAF", "NRAS", "PIK3CA"]
	gene_transcript = {
		"MET" : "NM_000245",
		"EGFR" : "NM_005228",
		"ALK" : "NM_004304",
		"KRAS" : "NM_033360",
		"ROS1" : "NM_002944",
		"RET" : "NM_020975",
		"ERBB2" : "NM_004448",
		"BRAF" : "NM_004333",
		"NRAS" : "NM_002524",
		"PIK3CA" : "NM_006218"
	}
	return getSum_forXW5301(var_list, gene_rule, gene_transcript)
jinja2.filters.FILTERS["sort_for5301_snvindel"] = sort_for5301_snvindel

# XW5301-sv����-2023.11.27
def sort_for5301_sv(var_list):
	gene_rule = ["ALK", "ROS1", "RET"]
	gene_transcript = {
		"ALK" : "NM_004304",
		"ROS1" : "NM_002944",
		"RET" : "NM_020975"
	}
	return getSum_forXW5301(var_list, gene_rule, gene_transcript)
jinja2.filters.FILTERS["sort_for5301_sv"] = sort_for5301_sv

# AD4701-tBRCA �����С��-2023.11.30
def ad4701_tbrca_sum(var_list):
	clinic_5_count = len([var for var in var_list if var["clinic_num_g"] == 5])
	clinic_4_count = len([var for var in var_list if var["clinic_num_g"] == 4])
	if clinic_5_count and clinic_4_count:
		return "����ʵ����{0}���²��Ա����{1}�������²��Ա��졣".format(str(clinic_5_count), str(clinic_4_count))
	elif clinic_5_count and not clinic_4_count:
		return "����ʵ����{0}���²��Ա��졣".format(str(clinic_5_count))
	elif not clinic_5_count and clinic_4_count:
		return "����ʵ����{0}�������²��Ա��졣".format(str(clinic_4_count))
	else:
		return "����ʵ��δ����²��Ի������²��Ա��졣"
jinja2.filters.FILTERS["ad4701_tbrca_sum"] = ad4701_tbrca_sum

# XW2402 CP ��ȡG12C��� 2023��12��27��
def filter_G12C(var_list):
	return [var for var in var_list if var['gene_symbol'] == 'KRAS' and var['hgvs_p'] == 'p.(G12C)']
jinja2.filters.FILTERS['filter_G12C'] = filter_G12C

# ��ȡ�����б��е�CNV-2023��12��27��
def filter_cnv(var_list):
	return [var for var in var_list if var["bio_category"] == "Cnv"]
jinja2.filters.FILTERS["filter_cnv"] = filter_cnv

# XW5101 ���ܽ��
def getSum_forXW5101(var_list, gene_rule):
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_rule:
				result.append(var)
	for gene in set(gene_rule) - set(detect_gene):
		result.append({'gene_symbol' : gene})
	return sorted(result, key = lambda i:gene_rule.index(i["gene_symbol"]))

def sort_5101_snv(var_list):
    gene_rule = ['FLT3', 'KMT2A', 'NPM1', 'NUP98', 'TP53']
    return getSum_forXW5101(var_list, gene_rule)
jinja2.filters.FILTERS['sort_5101_snv'] = sort_5101_snv

def sort_5101_sv(var_list):
    gene_rule = ['KMT2A', 'NPM1', 'NUP98']
    return getSum_forXW5101(var_list, gene_rule)
jinja2.filters.FILTERS['sort_5101_sv'] = sort_5101_sv

def getSum_forXW4205(var_list, gene_rule, gene_transcript):
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_rule:
				var["transcript_primary_xw4205"] = gene_transcript.get(gene)
				result.append(var)
	for gene in sorted(set(gene_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"transcript_primary_xw4205" : gene_transcript.get(gene)
		})
	return result

def sort_for4205_snvindel(var_list):
	#gene_rule = ["MET", "EGFR", "ALK", "KRAS", "ROS1", "RET", "ERBB2", "BRAF", "NRAS", "PIK3CA"]
	gene_rule = ['ALK', 'BRAF', 'EGFR', 'ERBB2', 'KRAS', 'MET', 'NRAS', 'PIK3CA', 'RET', 'ROS1']
	gene_transcript = {
		"ALK" : "NM_004304",
		"BRAF" : "NM_004333",
		"EGFR" : "NM_005228",
		"ERBB2" : "NM_004448",
		"KRAS" : "NM_033360",
		"MET" : "NM_000245",
		"NRAS" : "NM_002524",
		"PIK3CA" : "NM_006218",
		"RET" : "NM_020975",
		"ROS1" : "NM_002944"
	}
	return getSum_forXW4205(var_list, gene_rule, gene_transcript)
jinja2.filters.FILTERS["sort_for4205_snvindel"] = sort_for4205_snvindel

def filter_snv_cnv(var_list):
	return [var for var in var_list if var['bio_category'] in ['Snvindel', 'Cnv']]
jinja2.filters.FILTERS['filter_snv_cnv'] = filter_snv_cnv

def filter_hd_xw1402(info):
	qc = info[0]
	hd = info[1]
	result = []
	detect_gene = []
	for var in hd:
		if var['var_auto_result'] == 'T':
			detect_gene.append(var['gene_symbol'])
	if not detect_gene:
		result.append('HD����')
	else:
		if 'MTAP' in detect_gene and 'CDKN2A' in detect_gene:
			result.append('HD���ԣ����MTAP����ȱʧ��CDKN2A����ȱʧ')
		elif 'MTAP' in detect_gene:
			result.append('HD���ԣ����MTAP����ȱʧ��δ���CDKN2A����ȱʧ')
		elif 'CDKN2A' in detect_gene:
			result.append('HD���ԣ����CDKN2A����ȱʧ��δ���MTAP����ȱʧ')

	if qc['dna_data_qc']['snp_cover_ratio_num'] < 0.9:
		result = ['N/A']
	else:
		if not hd:
			result = ['HD����']
	
	sum = 'N/A' if result == ['N/A'] else 'HD����' if result == ['HD����'] else ''.join(result)

	return sum
jinja2.filters.FILTERS['filter_hd_xw1402'] =  filter_hd_xw1402

def filter_xw1402_dna(var_list):
	dna_rule = ['ALK', 'BRAF', 'BRCA1', 'BRCA2', 'CLDN18', 'CDKN2A', 'EGFR', 'ERBB2', 'FGFR2','IDH1', 'KEAP1', 'KRAS', 'MTAP', 'MET', 'MLH1', 'MSH2', 'MSH6', 'NTRK1', 'NTRK2', 'NTRK3', 'PMS2', 'RET', 'ROS1', 'STK11']
	result = []
	detect_gene = []
	for var in var_list:
		if var['gene_symbol'] in dna_rule:
			detect_gene.append(var['gene_symbol'])
			result.append(var)
	for gene in sorted(set(dna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	
	return sorted(result, key = lambda i:dna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1402_dna'] = filter_xw1402_dna
	

def filter_xw1402_rna(var_list):
	rna_rule = ['ALK', 'BRAF', 'CLDN18', 'EGFR', 'FGFR2', 'MET', 'NRG1', 'NTRK1', 'NTRK2', 'NTRK3', 'RET', 'ROS1']
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(',', var['gene_symbol'])):
			detect_gene.append(gene)
			if gene in rna_rule:
				if gene == 'MET':
					var['gene_symbol'] = 'MET'
				result.append(var)
	for gene in sorted(set(rna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	return sorted(result, key = lambda i:rna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1402_rna'] = filter_xw1402_rna

def filter_xw0258_hd(var_list):
	result = []
	if var_list:
		for var in var_list:
			if var['var_auto_result'] == "T" and var['type'] == "homodel":
				result.append(var)
	return result
jinja2.filters.FILTERS['xw0258_hd'] = filter_xw0258_hd

def filter_xw1404_dna(var_list):
	dna_rule = ['ALK', 'BRAF', 'CDKN2A', 'EGFR', 'ERBB2', 'KRAS', 'MTAP', 'MET', 'NTRK1', 'NTRK2', 'NTRK3', 'RET', 'ROS1']
	result = []
	detect_gene = []
	for var in var_list:
		if var['gene_symbol'] in dna_rule:
			detect_gene.append(var['gene_symbol'])
			result.append(var)
	for gene in sorted(set(dna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	
	return sorted(result, key = lambda i:dna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1404_dna'] = filter_xw1404_dna
	

def filter_xw1404_rna(var_list):
	rna_rule = ['ALK', 'BRAF', 'EGFR', 'MET', 'NTRK1', 'NTRK2', 'NTRK3', 'RET', 'ROS1']
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(',', var['gene_symbol'])):
			detect_gene.append(gene)
			if gene in rna_rule:
				if gene == 'MET':
					var['gene_symbol'] = 'MET'
				result.append(var)
	for gene in sorted(set(rna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	return sorted(result, key = lambda i:rna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1404_rna'] = filter_xw1404_rna

def filter_xw1405_dna(var_list):
	dna_rule = ['ALK', 'BRAF', 'BRCA1', 'BRCA2', 'CLDN18', 'CDKN2A', 'EGFR', 'ERBB2', 'FGFR2', 'HRAS', 'IDH1', 'IDH2', 'KEAP1', 'KRAS', 'MTAP', 'MET', 'MLH1', 'MSH2', 'MSH6', 'NRAS', 'NTRK1', 'NTRK2', 'NTRK3', 'PALB2', 'PMS2', 'RET', 'ROS1', 'STK11']
	result = []
	detect_gene = []
	for var in var_list:
		if var['gene_symbol'] in dna_rule:
			detect_gene.append(var['gene_symbol'])
			result.append(var)
	for gene in sorted(set(dna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	
	return sorted(result, key = lambda i:dna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1405_dna'] = filter_xw1405_dna

def filter_xw1405_rna(var_list):
	rna_rule = ['ALK', 'BRAF', 'CLDN18', 'EGFR', 'FGFR2', 'MET', 'NTRK1', 'NTRK2', 'NTRK3', 'NRG1', 'RET', 'ROS1']
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(',', var['gene_symbol'])):
			detect_gene.append(gene)
			if gene in rna_rule:
				if gene == 'MET':
					var['gene_symbol'] = 'MET'
				result.append(var)
	for gene in sorted(set(rna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	return sorted(result, key = lambda i:rna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1405_rna'] = filter_xw1405_rna

def xw1402_hd(hd):
	target = ['MTAP', 'CDKN2A']
	if hd:
		detect_gene = [var['gene_symbol'] for var in hd] if hd else []
		notdetect = list(set(target) - set(detect_gene))
		detect_gene = list(set(detect_gene))
		if len(detect_gene) == 2:
			return hd
		elif len(detect_gene) == 1:
			hd.append({
				'gene_symbol' : ''.join(notdetect),
				'region' : '/',
				'var_auto_result' : 'F'
			})
	hd = sorted(hd, key = lambda i:target.index(i["gene_symbol"]))
	return hd
jinja2.filters.FILTERS['xw1402_hd'] = xw1402_hd

def ad3101_inter(info):
	var = info[0]
	sample = info[1]
	if sample['product_name'] == 'AD3101' and sample['prod_names'] == 'Master Panel����֯��':
		return var["var_somatic"]["level_I"] + var["var_germline"]["level_5"] + var["var_somatic"]["level_II"] + var["var_somatic"]["level_onco_nodrug"] + var["var_germline"]["level_4"]
jinja2.filters.FILTERS['ad3101_inter'] = ad3101_inter