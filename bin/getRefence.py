#-*- coding:gbk -*-
import xlrd
import os
import re
from libs.getInterRef import getRef_from_inter
'''
Discription
	
	�ýű�������ȡ�̶��ο����׺Ͷ�̬�ο����ס� 

'''
#BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#fixed_refer_path = os.path.join(BASE_DIR, "config/reference.xlsx")

# ��ʽ�������ļ�
def stran_xlrd(sheet_name, config):
	fixed_refer_path = os.path.join(config, "report_requirenment.xlsx")
	xls = xlrd.open_workbook(fixed_refer_path)
	refer_sheet = xls.sheet_by_name(sheet_name)
	key = refer_sheet.row_values(0)
	Data = []
	for num in range(1, refer_sheet.nrows):
		rows = refer_sheet.row_values(num)
		tmpdict = {}
		for i in range(len(key)):
			tmpdict[key[i]] = rows[i]
		Data.append(tmpdict)
	return Data	

# �̶����׸�ʽ����
def stran_fixed_refer(i):
	# authors ȡǰ������","���ӣ��޿ո񣬺��et al.
	authors = i['authors'].split(',')
	if len(authors) > 3:
		i['authors'] = ','.join([authors[0], authors[1], authors[2], 'et al.'])
	# date ��ʽΪ(2013)
	if i['date']:
		i['date'] = '(' + re.search(r'\d{4}', str(i["date"]).replace(".0", "")).group(0) + ')'
	# title ����"."
	if i['title']:
		i['title'] = i['title'].strip('.') + '.'
	# pmid ��[]
	if i['PMID']:
		i['PMID'] = ''.join(['[', 'PMID:', str(int(i['PMID'])), ']'])
	result = " ".join([i["authors"], i["date"], i["title"], i["journal"], i["vol"], i["PMID"]])

	return result


# �̶��ο�����
# ����Ҫθ�����ӷ��͡�GEP��TME����
def getfixed_refer(report_name, tumor_list, KNB, MSI, config):
	Data = stran_xlrd("reference-fixed", config)
	# ƥ��
	fixed_refer = []
	refer_type_dict = {}
	for i in [refer for refer in Data if refer["docx_template"] == report_name]:
		refer_type_dict.setdefault(i["refer_type"], [])
		refer_type_dict[i["refer_type"]].append(i)
	# 1. fixed���̶��ο�����ֱ��չʾ
	if "fixed" in refer_type_dict.keys():
		fixed_refer += refer_type_dict["fixed"]
	# 2. tumor: ֻƥ�䰩��
	if "tumor" in refer_type_dict.keys():
		for refer in refer_type_dict["tumor"]:
			if refer["tumor_or_type"] in tumor_list:
				fixed_refer.append(refer)
	# 3. tumor_and_result: ƥ�䰩�ֺͽ����Ŀǰ��4��
	if "tumor_and_result" in refer_type_dict.keys():
		for refer in refer_type_dict["tumor_and_result"]:
			# 3.1. θ�����ӷ���GA_result
			#    ��Ϊθ������tumor_or_typeΪGA_type�����ݷ���ѡ���Ӧ�ο����ף�GA_EB��GA_MSI��GA_stable��
			# �̶����ף��������ַ���
#			if "θ��" in tumor_list and refer["tumor_or_type"] == "GA_type":
#				pass

			# 3.2. KNB
			if refer["tumor_or_type"] == "KNB" and KNB:
				fixed_refer.append(refer)

			# 3.3. GEP
			if "�ΰ�" in tumor_list and refer["tumor_or_type"] == "GEP":
				fixed_refer.append(refer)

			# 3.4. TME
			if "�ΰ�" in tumor_list and refer["tumor_or_type"] == "TME":
				fixed_refer.append(refer)

	# 4. result �� ��ƥ������Ŀǰ��1��
	if "result" in refer_type_dict.keys():
		for refer in refer_type_dict["result"]:
			# 4.1. MSI-H
			if refer["tumor_or_type"] == "MSI-H" and "var_id" in MSI.keys() and MSI["var_id"] and MSI["var_id"] == "MSI-H":
				fixed_refer.append(refer)

	fixed_refer_stran = [stran_fixed_refer(i) for i in fixed_refer]
	
	return fixed_refer_stran

# ��̬�ο�����
# ����̬���׽��в�֣����水��ʹ��
def getdynamic_refer(jsonDict, var_data, msi, hrd, var_brca):
	dynamic_refer = {}
	# Ĭ�ϲο�����
	# 1. s_var12����ȡ������ܡ������������ơ���ϡ�Ԥ���Ŵ����ղο����ף�����I/II��
	dynamic_refer["s_var12"] = []
	for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]:
		dynamic_refer["s_var12"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["s_var12"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["s_var12"].extend(var["evi_sum"]["refer_evi"])
	# 2. s_var_onco_nodrug ��ȡ������ܡ�����������������������չ��ر���
	dynamic_refer["s_var_onco_nodrug"] = []
	for var in var_data["var_somatic"]["level_onco_nodrug"]:
		dynamic_refer["s_var_onco_nodrug"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["s_var_onco_nodrug"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# 3. s_var3����ȡ������ܡ�������������III�����
	dynamic_refer["s_var3"] = []
	for var in var_data["var_somatic"]["level_III"]:
		dynamic_refer["s_var3"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["s_var3"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# 4. g_var12: ��ȡ������ܡ������������ơ���ϡ�Ԥ���Ŵ����ղο�����
	dynamic_refer["g_var45"] = []
	for var in var_data["var_germline"]["level_4"] + var_data["var_germline"]["level_5"]:
		dynamic_refer["g_var45"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var45"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["g_var45"].extend(var["evi_sum"]["refer_evi"])

	### MLPA ����©����-del�ӵ�g_var45�У���쿷�-2023.11.07
	for var in var_brca["mlpa"]["B1_Loss"] + var_brca["mlpa"]["B2_Loss"]:
		dynamic_refer["g_var45"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var45"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["g_var45"].extend(var["evi_sum"]["refer_evi"])
	### MLPA del������-2023.11.07

	# 5. g_var3: ��ȡ������ܡ�������
	dynamic_refer["g_var3"] = []
	for var in var_data["var_germline"]["level_3"]:
		dynamic_refer["g_var3"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var3"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))

	### MLPA ����©����-dup�ӵ�g_var3�У���쿷�-2023.11.07
	for var in var_brca["mlpa"]["B1_Gain"] + var_brca["mlpa"]["B2_Gain"]:
		dynamic_refer["g_var3"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var3"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	### MLPA del������-2023.11.07

	# 5. KNB����ȡKNB����ʱ�����ơ���ϡ�Ԥ���Ŵ����ղο�����
	dynamic_refer["knb"] = var_data["knb"]["evi_sum"]["refer_evi"] if var_data["knb"] else []
	# 6. EC��EC�����ٴ��������ƾ���
	dynamic_refer["ec_type"] = var_data["ec_type"]["evi_sum"]["refer_evi"] if var_data["ec_type"] and "evi_sum" in var_data["ec_type"].keys() else []
	# 7. MSI��MSI-Hʱ�����Ʋ���
	dynamic_refer["msi"] = msi["evi_sum"]["refer_evi"] if msi and msi["var_id"] == "MSI-H" else []
	# 8. HRD�����Ʋ���
	dynamic_refer["hrd"] = hrd["evi_sum"]["refer_evi"] if hrd and hrd["evi_sum"] else []
	# 2023.06.08-hrd�����gss�еĲο�����
	if "gss" in var_data.keys() and var_data["gss"] and "gss" in var_data["gss"].keys() and var_data["gss"]["gss"] and "evi_sum" in var_data["gss"]["gss"].keys()\
		and var_data["gss"]["gss"]["evi_sum"]:
		dynamic_refer["hrd"].extend(var_data["gss"]["gss"]["evi_sum"]["refer_evi"])

	'''
	# 9. ������HRD C����չʾI/II�������
	dynamic_refer["s_var12_syx"] = []
	for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]:
		dynamic_refer["s_var12_syx"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# ������������III�������-2023.05.16
	dynamic_refer["s_var3_syx"] = []
	for var in var_data["var_somatic"]["level_onco_nodrug"] + var_data["var_somatic"]["level_III"]:
		dynamic_refer["s_var3_syx"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# ����������ʡ��116 10 ������III��������+������-2023.05.24
	dynamic_refer["fjsl_116_lc10III"] = []
	for var in var_data["special"]["var_pan116_split"]["gene10_level_III"]:
		if var["bio_category"] == "Sv" and "," in var["gene_symbol"]:
			dynamic_refer["fjsl_116_lc10III"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["fjsl_116_lc10III"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["fjsl_116_lc10III"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["fjsl_116_lc10III"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	'''

	# 10. �ֱ��ȡ��ϵ�������ͻ������еĲο����ף� 4/5����컹����ȡ�Ŵ������еĲο�����-2023.02.23
	dynamic_refer["g_var45_inter"] = []
	dynamic_refer["g_var45_gene"] = []
	dynamic_refer["g_var45_risk"] = []
	# ����MLPA del-2023.11.07
	#for var in var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"]:
	for var in var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"] + var_brca["mlpa"]["B1_Loss"] + var_brca["mlpa"]["B2_Loss"]:
		dynamic_refer["g_var45_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["g_var45_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var45_risk"].extend(var["evi_sum"]["refer_evi_risk"])
	dynamic_refer["g_var3_inter"] = []
	dynamic_refer["g_var3_gene"] = []
	for var in var_data["var_germline"]["level_3"] + var_brca["mlpa"]["B1_Gain"] + var_brca["mlpa"]["B2_Gain"]:
		dynamic_refer["g_var3_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["g_var3_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
	
	# 11. ��ȡI/II�����ı������ͻ�����ܲο�����-2023.09.01
	dynamic_refer["s_var12_gene"] = []
	dynamic_refer["s_var12_inter"] = []
	for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]:
		dynamic_refer["s_var12_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["s_var12_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# ��������-2023.09.01
	# ����ϸ������Ļ�����ܡ����������ٴ�֤�ݲ�ֿ����ˣ�ģ�����������-2023.11.03
	dynamic_refer["s_var12_evi"] = []
	for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]:
		dynamic_refer["s_var12_evi"].extend(var["evi_sum"]["refer_evi"])

	dynamic_refer["s_varonco_gene"] = []
	dynamic_refer["s_varonco_inter"] = []
	for var in var_data["var_somatic"]["level_onco_nodrug"]:
		dynamic_refer["s_varonco_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["s_varonco_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	
	dynamic_refer["s_var3_gene"] = []
	dynamic_refer["s_var3_inter"] = []
	for var in var_data["var_somatic"]["level_III"]:
		dynamic_refer["s_var3_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["s_var3_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
 
	return dynamic_refer	 