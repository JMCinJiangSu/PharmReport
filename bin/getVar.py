#-*- coding:gbk -*-
import re
from libs.processTheSnvindel import process_snvindel
from libs.processTheCNV import process_cnv
from libs.processTheSV import process_sv
from libs.processTheKNB import process_knb
from libs.processTheEC_type import process_ec_type
from libs.getInterRef import getRef_from_inter
from libs.NCCN_Gene_detect import getNCCN_detect
from libs.mergeMET14 import mergeVar
from libs.detectResultSummary import getSummary_detect
import copy
from libs.specialRequest import judge_GA_tumor_KNB
from libs.specialRequest import sv_shfk
#from libs.io_detect import io_detect
#from libs.io_detect import io_detect_for_116
#from libs.io_detect import io_detect_for_ZJZL
from libs.io_detect import get_io_detect
from libs.specialRequest import CP40_split_gene, PAN116_split_gene
from libs.specialRequest import CP40_FJFY_summary
from libs.processTheMLPA import process_mlpa
from libs.processTheGA_type import process_ga_type
from libs.specialRequest import Master_sum
from libs.getHRD_Master import getgss
from libs.rule import s_var_rule, s_var_rule_gene, s_var_rule_genelist, g_var_rule, g_var_regimen_rule, var_regimen_rule, g_var_rule_genelist
from libs.rule import process_result_regimen_ZJZL, var_bptm_rule, var_lyn5_rule, nofoundPath_genelist, nofound_genelist
from libs.specialRequest import HRR_FJFY_summary, varS3_splitfor301, sort_var_for301, sort_var_for301_2, get_summary_GDRM_gHRR
from libs.specialRequest import varInter_FJZL, var_summary_GDRM, var_GZZL, var_summary_SHRJ, var_summary_FJZL, var_SD
from libs.specialRequest import varInfo_XAJDY, var_summary_CQFY, BRCA_FJFY_summary, var_summary_WHXH, var_level123_SDZL, brca_judge_inApprovalTumor, varnum_ZDY_HRR
from libs.specialRequest import shrj_cp40_var
from libs.specialRequest import sort_for_xw0417
from libs.specialRequest import getSum_for_tXW6002
from libs.specialRequest import getSum_for_gXW6002
from libs.specialRequest import getSum_for_XW5902
from libs.specialRequest import var_ad3101_summary, getSum_for_tXW6701, getSum_for_gXW6701, getSum_for_XW6003

'''
Discription 
	
	��ȡjson�ļ��е�snvindel/cnv/sv/rna_sv��mlpa��ת��Ϊ�����ڴ󲿷ֱ���ģ�巽�����ĸ�ʽ����ʽ�ݶ����£����ȼ���ֱ������㲻ͬģ����������
	var : {
		var_somatic_without_rnasv : {  Ϊ�����㱨��ģ����RNA SV��DNA SVһ����ϲ���һ�����չʾ����������ǲ𿪵ģ�����DNA SV��DNA��RNA����SV
			level_I : [],  I����죬������Ʒ����ȼ�ΪA/B
			level_II : [], II����죬������Ʒ����ȼ�ΪC/D
			level_onco_nodrug : [], ��ҩ������������չ����
			level_III : [], III�����
			},
		var_rna_sv : {   Ϊ�����㱨��ģ����RNA SV��DNA SVһ����ϲ���һ�����չʾ����������ǲ𿪵ģ� ����RNA SV��DNA��RNA����SV
			level_I : [],
			level_II : [],
			level_noco_nodrug : [],
			level_III : []
		},
		var_somatic : {   Ϊ�����㱨��ģ����RNA SV��DNA SVһ����ϲ���һ�����չʾ����������Ǻϲ���
			level_I : [],  I����죬������Ʒ����ȼ�ΪA/B
			level_II : [], II����죬������Ʒ����ȼ�ΪC/D
			level_onco_nodrug : [], ��ҩ������������չ����
			level_III : [], III�����	
		},
		var_germline : {
			level_5 : [], �²��Ա���
			level_4 : [], �����²��Ա���
			level_3 : [], ���岻������
			level_2 : [], �������Ա��죨����ģ����Ҫչʾ���ԵĽ������HRR��
			level_1 : [], ���Ա��죨����ģ����Ҫչʾ���ԵĽ������HRR��
			regimen_level_I : [],
			regimen_level_II : []
		},
		var_for_regimen_without_rnasv : { ��ϵ+��ϸ���������Ʒ����ı���, ����DNA SV��DNA��RNA����SV
			level_I : [],
			level_II : []
		},
		var_for_regimen : { ��ϵ+��ϸ���������Ʒ����ı���
			level_I : [],
			level_II : []
		},
		knb : {},
		var_ec_type : { ������CP40��PTM��BPTM��ͬչʾ����ı���ģ��
			POLE_level12 : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ֿ�չʾ�ı��
			POLE_level3 : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ֿ�չʾ�ı��
			TP53_level12 : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ֿ�չʾ�ı��
			TP53_level3 : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ֿ�չʾ�ı��
			BRCA1_level12 : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ֿ�չʾ�ı��
			BRCA1_level3 : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ֿ�չʾ�ı��
			BRCA2_level12 : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ֿ�չʾ�ı��
			BRCA2_level3 : []�� �������ӹ���Ĥ�����ӷ��ͺͱ������ֿ�չʾ�ı��
			POLE_level12_withECtype : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ϲ�չʾ�ı��
			POLE_level3_withECtype : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ϲ�չʾ�ı��
			TP53_level12_withECtype : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ϲ�չʾ�ı��
			TP53_level3_withECtype : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ϲ�չʾ�ı��
			BRCA1_level12_withECtype : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ϲ�չʾ�ı��
			BRCA1_level3_withECtype : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ϲ�չʾ�ı��
			BRCA2_level12_withECtype : [], �������ӹ���Ĥ�����ӷ��ͺͱ������ϲ�չʾ�ı��
			BRCA2_level3_withECtype : []�� �������ӹ���Ĥ�����ӷ��ͺͱ������ϲ�չʾ�ı��
		},
		cdx : {}, ��������Ƽ�����������
		summary_result : {} ���������

'''

def getVar(jsonDict, config, report_name):
	data = {}
	snvindel = copy.deepcopy(process_snvindel(jsonDict, config))
	cnv = copy.deepcopy(process_cnv(jsonDict, config))
	sv_combination, rna_sv, rna_sv_only = copy.deepcopy(process_sv(jsonDict, config))
	mlpa, mlpa_image, mlpa_image_del = process_mlpa(jsonDict, config)

	# ����Ϊ�Ƿΰ���3�������Ҳ��չʾBCL2L11-2022.09.06
	# ��Ժģ��ҪչʾBCL2L11����3�࣬����ټӸ��ټ������-2022.09.14
	# ����ټӸ��ټ�ͨ��ģ�������
	# 2023.06.01-���µ��ټ�116/76ģ�壬��¼��󲻷�BCL2L11�ˣ��ò���ע�͵�����������в���չʾBCL2L11���죬����������˲�����˵���
	#if re.search("Pan116|LC76", jsonDict["sample_info"]["prod_names"]) and jsonDict["sample_info"]["report_module_type"] == "rummage" and re.search("rummage", report_name):
	#	for var in snvindel:
	#		if var["gene_symbol"] == "BCL2L11" and var["hgvs_c"] and var["hgvs_c"] == "c.394+1479_394+4381del":
	#			snvindel.remove(var)

	

	### ��������
	# �����㽭����Master�ı��������ͼ쵥λ���㽭ʡ����ҽԺ����Ʒ��Master��ҵ�����ͣ���Ժ
	# ������򣺰���ϵ/��ϸ����I/II/����������չ���/III�ࡢ��ͬ�ȼ������Ʒ�����ߵȼ���snvindel>sv>cnv��Ƶ�ʽ�
	if re.search("�㽭ʡ����ҽԺ", jsonDict["sample_info"]["company"]) and re.search("Master|3231����", jsonDict["sample_info"]["prod_names"]) and jsonDict["sample_info"]["report_module_type"] == "hospital":
		var_origin_rule = {"somatic" : 0, "germline" : 1}
		top_level_rule = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "N" : 4}
		var_type_rule = {"Snvindel" : 0, "Sv" : 1, "PSeqRnaSv" : 2, "Cnv" : 3}
		var_data_without_rnasv = sorted(snvindel + cnv + sv_combination, key=lambda i : (var_origin_rule.get(i["var_origin"], 0), top_level_rule.get(i["top_level"]), var_type_rule.get(i["bio_category"])))
		var_data_rna_sv = sorted(rna_sv, key=lambda i : (var_origin_rule.get(i["var_origin"], 0), top_level_rule.get(i["top_level"]), var_type_rule.get(i["bio_category"])))
		var_data =  sorted(snvindel + cnv + sv_combination + rna_sv_only, key=lambda i : (var_origin_rule.get(i["var_origin"], 0), top_level_rule.get(i["top_level"]), var_type_rule.get(i["bio_category"])))
	
	# ����ģ�嶼��������������
	else:
		# ���� ��I��II���������͡�Ƶ�ʣ�Ƶ���ڸ����ִ����ʱ���Ѿ������ˣ�
		rule = {"Snvindel" : 0, "Cnv" : 1, "Sv" : 2, "PSeqRnaSv" : 3}
		# ����������򣬰�I/II��ͬ�ȼ�ʱ�����Ʒ�����ߵȼ��������򣩡��������͡�Ƶ�ʣ�Ƶ���ڸ����ִ����ʱ���Ѿ������ˣ�
		# ͬ�ȼ�ʱ�����Ʒ�����ߵȼ��������� A > B > C > D > N
		top_level_rule = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "N" : 4}
		#var_data_without_rnasv = sorted(snvindel + cnv + sv_combination, key=lambda i : (i["clinic_num_s"], top_level_rule.get(i["top_level"]), i["gene_symbol"], rule.get(i["bio_category"])))
		#var_data_rna_sv = sorted(rna_sv, key=lambda i : (i["clinic_num_s"], top_level_rule.get(i["top_level"]), i["gene_symbol"], rule.get(i["bio_category"])))
		#var_data = sorted(snvindel + cnv + sv_combination + rna_sv_only, key=lambda i : (i["clinic_num_s"], top_level_rule.get(i["top_level"]), i["gene_symbol"], rule.get(i["bio_category"])))
		# ���������򣬰�clinic_num_s�ŵ�top_level���棬����clinic_num_s��������ת��-2023.02.01
		clinic_num_s_rule = {5 : 0, 4 : 1, 3 : 2, 2 : 3, 1 : 4}
		var_data_without_rnasv = sorted(snvindel + cnv + sv_combination, key=lambda i : (top_level_rule.get(i["top_level"]), clinic_num_s_rule.get(i["clinic_num_s"]), i["gene_symbol"], rule.get(i["bio_category"])))
		var_data_rna_sv = sorted(rna_sv, key=lambda i : (top_level_rule.get(i["top_level"]), clinic_num_s_rule.get(i["clinic_num_s"]), i["gene_symbol"], rule.get(i["bio_category"])))
		var_data = sorted(snvindel + cnv + sv_combination + rna_sv_only, key=lambda i : (top_level_rule.get(i["top_level"]), clinic_num_s_rule.get(i["clinic_num_s"]), i["gene_symbol"], rule.get(i["bio_category"])))
	### 

	# �ϲ�MET 14��Ծͻ��
	if jsonDict['sample_info']['product_name'] != 'XW1402':
		var_data_without_rnasv, judge_mergeMET_without_rnasv = mergeVar(var_data_without_rnasv)
		var_data, judge_mergeMET = mergeVar(var_data)

	# ��ϸ���������
	data["var_somatic_without_rnasv"] = s_var_rule(var_data_without_rnasv)
	data["var_somatic_rna_sv"] = s_var_rule(var_data_rna_sv)
	data["var_somatic"] = s_var_rule(var_data)

	# ��ϵ�������
	data["var_germline"] = {**g_var_rule(var_data), **g_var_regimen_rule(var_data)}
	
	# ��ϵ+��ϸ������ҩ����/Ԥ��/������ϱ�������
	data["var_for_regimen"] = var_regimen_rule(var_data)
	data["var_for_regimen_without_rnasv"] = var_regimen_rule(var_data_without_rnasv)

	# �㽭������ϵ+��ϸ������ҩ����/Ԥ��/������ϱ������������н�չʾ��ҩ������ҩ����Ԥ��/������ϵ�д�����ް�����ҩ��ʾ����ģ������н���ʵ�֣�������ټ��ֶΣ�
	# �����ڰ�����ҩ��ʾ���֣�������ҩ���ͻ��ǿ�����ͨ�õ�var_for_regimen
	data["var_for_regimen_ZJZL"] = copy.deepcopy(data["var_for_regimen"])
	data["var_for_regimen_ZJZL"]["level_I"] = process_result_regimen_ZJZL(data["var_for_regimen_ZJZL"]["level_I"])
	data["var_for_regimen_ZJZL"]["level_II"] = process_result_regimen_ZJZL(data["var_for_regimen_ZJZL"]["level_II"])
	# ������ϵ�²�/�����²����ް�����ҩ��λ����ͷ���5.1 �ٴ�������ȷ�ı��켰������ҩ������߼��ܹ֣������ǵ�������ߵ�����һ���б�չʾ��
	data["var_germline_nodrug"] = [var for var in var_data if var["clinic_num_g"] in [5, 4] and var["var_origin"] == "germline" and (("evi_sum" in var.keys() and not var["evi_sum"]["evi_split"]) or ("evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys())))]

	# knb
	data["knb"] = copy.deepcopy(process_knb(jsonDict, config))

	# �ӹ���Ĥ��������BRCA1��BRCA2�Ľ������������HRD��
	data["ec_type"] = {}
	if jsonDict["ec_type"]:
		data["ec_type"] = copy.deepcopy(process_ec_type(jsonDict, config))
		if data["ec_type"]["evi_sum"]["evi_split"] and "Prognostic" in data["ec_type"]["evi_sum"]["evi_split"].keys():
			# Ԥ�� Ŀǰֻ��better��poor��������ĳ��������б�
			data["ec_type"]["clinical_significance"] = [i["clinical_significance"] for i in data["ec_type"]["evi_sum"]["evi_split"]["Prognostic"]]
			# ֤���������ϲ�Ϊ�ַ���
			data["ec_type"]["evi_interpretation"] = "".join([i["evi_interpretation"] for i in data["ec_type"]["evi_sum"]["evi_split"]["Prognostic"]])
			# ĿǰCP40�е��ٴ�֤�ݣ�����ʾ�������ʾ�����������
			data["ec_type"]["evi_inter_cp40"] = re.split("�����ʾ��", data["ec_type"]["evi_interpretation"] )[-1] if data["ec_type"]["evi_interpretation"] else ""
	# TP53��POLE��BRCA1��BRCA2������������ΪI/II���III�࣬���ڱ���ͷ��ӷ��ͷֿ�չʾ�����
	ec_gene_list = ["TP53", "POLE", "BRCA1", "BRCA2"]
	for gene in ec_gene_list:
		data["ec_type"].update(var_bptm_rule(var_data, gene))

	# ����BPTM��ȫѪ�������������Ŀvar_originΪgermline, ��������ϸ����-20221027
	# BRCA������BRCA��Ŀ�����ݽṹ����߽�����POLE��TP53
	data["ec_type"]["POLE_g"] = g_var_rule_genelist(var_data, ["POLE"])
	data["ec_type"]["TP53_g"] = g_var_rule_genelist(var_data, ["TP53"])
	# ��������-20221027

	# MLPA���
	#data["mlpa"], data["mlpa_image"], data["mlpa_image_del"] = process_mlpa(jsonDict, config)
	data["mlpa"] = mlpa
	data["mlpa_image"] = mlpa_image
	data["mlpa_image_del"] = mlpa_image_del

	# ������ϼ����
	data["cdx"] = getNCCN_detect(var_data, jsonDict["sample_info"]["tumor_list"], data["mlpa"], config)

	# ������б�
	#data["summary_result"] = getSummary_detect(var_data, data["mlpa"], config)

	# �����������
	#data["io"] = {}
	#data["io"]["result"], data["io"]["io_p_summary"], data["io"]["io_n_summary"] = io_detect(var_data)
	#data["io"]["result_116"], data["io"]["io_p_summary_116"], data["io"]["io_n_summary_116"] = io_detect_for_116(var_data)
	#data["io"]["io_p_summary_ZJZL"], data["io"]["io_n_summary_ZJZL"] = io_detect_for_ZJZL(var_data, config)
	data["io"] = get_io_detect(var_data, config)

	# ��ϵ-�����ۺ����������������Ҫ�ֿ�չʾ��EPCAM��MLH1��MSH2��MSH6��PMS2��
	# չʾ���������3��4��5�����
	gLS5_gene_list = ["EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	data["gLS5"] = var_lyn5_rule(var_data, gLS5_gene_list)

	# CRC12-MSI��Ҫչʾδ��⵽I/II�����Ļ���
	crc12_gene = ["KRAS", "NRAS", "BRAF", "POLE", "PIK3CA", "ERBB2", "ALK", "FGFR3", "NTRK1", "NTRK3", "RET"]
	data["crc12_withoutPathVar_geneList"] = nofoundPath_genelist(var_data, crc12_gene)

	# LC10��չʾδ�쵽I/II�����Ļ���
	lc10_gene = ["BRAF", "EGFR", "ERBB2", "KRAS", "MET", "ALK", "ROS1", "RET", "NRAS", "PIK3CA"]
	data["lc10_withoutPathVar_geneList"] = nofoundPath_genelist(var_data, lc10_gene) 
	# LC10�ټ�ͨ��ģ��չʾδ��⵽I/II/III�����Ļ���
	data["lc10_withoutVar_geneList"] = nofound_genelist(var_data, lc10_gene)

	# ��������
	data["special"] = {}
	# 1. ������ɽCP40���ж�������������KNB������
	#data["special"]["FDZS_GA_KNB"] = judge_GA_tumor_KNB(var_data, jsonDict["sample_info"]["tumor_list"])
	# 1. ������ɽCP40���ж�������������KNB������-����Ϊֻƥ��tumor_names_cn-20220927
	#data["special"]["FDZS_GA_KNB"] = judge_GA_tumor_KNB(var_data, [jsonDict["sample_info"]["tumor_names_cn"]])
	# 2. �Ϻ��ο�CP40�ںϷ���չʾ
	#data["special"]["SHFK_SV"] = sv_shfk(var_data)
	# 3. ������һ���ĳ�CP40 I/II������Ϊ10����+30����չʾ
	'''
	data["special"]["var_cp40_split"] = {}
	data["special"]["var_cp40_split"]["gene10_level_I_II"], data["special"]["var_cp40_split"]["gene30_level_I_II"], \
	data["special"]["var_cp40_split"]["gene10_level_III"], data["special"]["var_cp40_split"]["gene30_level_III"],\
	data["special"]["var_cp40_split"]["gene10_level_nj_onco_nodrug"], data["special"]["var_cp40_split"]["gene30_level_nj_onco_nodrug"],\
	data["special"]["var_cp40_split"]["gene10_level_nj_III"], data["special"]["var_cp40_split"]["gene30_level_nj_III"] = CP40_split_gene(data)
	# 4. ������һCP40 I��II�����С��
	data["special"]["var_cp40_FJFY"] = {}
	#data["special"]["var_cp40_FJFY"]["level_I_sum"], data["special"]["var_cp40_FJFY"]["level_II_sum"] = CP40_FJFY_summary(data["var_somatic"]["level_I"], data["var_somatic"]["level_II"])
	data["special"]["var_cp40_FJFY"] = CP40_FJFY_summary(data["var_somatic"]["level_I"], data["var_somatic"]["level_II"],data["var_somatic"]["level_III"],data["var_somatic"]["level_onco_nodrug"])
	'''
	# 5. ͳ��Master��ϸ����ʽ��DNA/RNA�ֿ�ͳ�ƣ�
	data["special"]["Master_level_I"], data["special"]["Master_level_II"], data["special"]["Master_level_onco_nodrug"], data["special"]["Master_level_III"] = Master_sum(var_data)
	# 6. ����I/I/III�����Ļ����б�
	'''
	data["special"]["SYX_levelI"] = "��".join(sorted(set([var["gene_symbol"] for var in data["var_somatic"]["level_I"]])))
	data["special"]["SYX_levelII"] = "��".join(sorted(set([var["gene_symbol"] for var in data["var_somatic"]["level_II"]])))
	data["special"]["SYX_levelIII"] = "��".join(sorted(set([var["gene_symbol"] for var in data["var_somatic"]["level_III"] + data["var_somatic"]["level_onco_nodrug"]])))
	# 7. ��ϵ�����������б�����BRCA dup
	data["special"]["gene_var3_list"] = "��".join(sorted(set([var["gene_symbol"] for var in (data["var_germline"]["level_3"]+data["mlpa"]["B1_Gain"]+data["mlpa"]["B2_Gain"])])))
	'''
	#print (data["special"]["gene_var3_list"])
	# 8. ����DNA/RNA MET14�����ֶ�
	data["special"]["judge_mergeMET"] = judge_mergeMET if jsonDict['sample_info']['product_name'] != 'XW1402' else ''
	# 9. ������һHRR�ܽ�
	'''
	data["special"]["FJFY_HRR_sum"] = HRR_FJFY_summary(data)
	# 10. 301 MP ������죬һ��չʾ��������
	data["special"]["varS3_301"] = varS3_splitfor301(data["var_somatic"]["level_III"])
	# 10. 301 gHRR ������죬һ��չʾ�������죬���迼��MLPA dup-2023.09.22
	data["special"]["varS3_301_germline"] = varS3_splitfor301(data["mlpa"]["B1_Gain"]+data["mlpa"]["B2_Gain"]+data["var_germline"]["level_3"])
	# 2023.09.22������
	# 11. 301 MP ��������
	data["special"]["var_sort_301"] = {}
		# 11.1 �ж���ʱ�ġ������ٴ�����ı��족
	data["special"]["var_sort_301"]["Match_drug"] = sort_var_for301(data["var_somatic"]["level_I"] + data["var_somatic"]["level_II"])
		# 11.2 �ж���ʱ�ġ�����������չ��ر��족
	data["special"]["var_sort_301"]["Match_withoutdrug"] = sort_var_for301(data["var_somatic"]["level_onco_nodrug"])
		# 11.3 �޶���ʱ�ġ������ٴ�����ı��족
	data["special"]["var_sort_301"]["single_drug"] = sort_var_for301(data["var_germline"]["regimen_level_I"] + data["var_somatic"]["level_I"] +\
																	data["var_germline"]["regimen_level_II"] + data["var_somatic"]["level_II"])
		# 11.4 �޶���ʱ�ġ�����������չ��ر��족
	data["special"]["var_sort_301"]["single_withoutdrug"] = sort_var_for301(data["var_germline_nodrug"] + data["var_somatic"]["level_onco_nodrug"])

	# 12. 301 MP �������2
	data["special"]["var_sort_301_2"] = {}
	# 12.1 �ж���ʱ�ġ���ϸ�����족
	data["special"]["var_sort_301_2"]["somatic"] = sort_var_for301_2(data["var_somatic"]["level_I"] + data["var_somatic"]["level_II"] +\
																	data["var_somatic"]["level_onco_nodrug"])
	# 12.2 �޶���ʱ�ġ�������족
	data["special"]["var_sort_301_2"]["var"] = sort_var_for301_2(data["var_germline"]["regimen_level_I"] + data["var_somatic"]["level_I"] +\
																data["var_germline"]["regimen_level_II"] + data["var_somatic"]["level_II"]+\
																data["var_germline_nodrug"] + data["var_somatic"]["level_onco_nodrug"])
	
	# 13. �㶫ʡ����gHRR���С��
	data["special"]["summary_GDRM_gHRR"], data["special"]["summary_GDRM_tHRR"], data["special"]["summary_GDRM_ptHRR"] = get_summary_GDRM_gHRR([var for var in data["var_germline"]["level_5"] if var["bio_category"] == "Snvindel"]\
																+data["mlpa"]["B1_Loss"]\
																+data["mlpa"]["B2_Loss"]\
																+[var for var in data["var_germline"]["level_4"] if var["bio_category"] == "Snvindel"], \
																[var for var in data["var_somatic"]["level_I"] + data["var_somatic"]["level_II"] if var["bio_category"] == "Snvindel"])
	# 14. ������һ 116 I/II������Ϊ10����+��������չʾ
	data["special"]["var_pan116_split"] = {}
	data["special"]["var_pan116_split"]["gene10_level_I_II"], \
	data["special"]["var_pan116_split"]["other_level_I_II"], \
	data["special"]["var_pan116_split"]["gene10_level_III"], \
	data["special"]["var_pan116_split"]["other_level_III"], \
	data["special"]["var_pan116_split"]["gene10_germline_45"], \
	data["special"]["var_pan116_split"]["other_germline_45"] = PAN116_split_gene(data)
	# 15. ͳ����ϸ��/��Դδ�����Ƕ�̬�Ա�����-����֣��һHRR
	data["special"]["varnum_ZDY_HRR"] = varnum_ZDY_HRR(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_III"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_IV"])
	# 16. �Ϻ��ʼ�CP40-��Ҫ�������������б��˳�򱣳ֲ��䣨I/II/III�����������ͻ������Ϊ��Ҫ����I/II��+��������I/II�࣬���������Ҫȥ��-2023.09.18
	data["special"]["shrj_cp40"] = {}
	data["special"]["shrj_cp40"]["summary_cdx"],\
	data["special"]["shrj_cp40"]["var_inter"],\
	data["special"]["shrj_cp40"]["gene_inter"],\
	data["special"]["shrj_cp40"]["gene_inter2"] = shrj_cp40_var(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"], config)
	#print (data["special"]["shrj_cp40"]["var_inter"])
	'''
	# xw0417���
	data['special']['xw0417'] = sort_for_xw0417(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"])
	data['special']['tXW6002'] = getSum_for_tXW6002(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"])
	data['special']['gXW6002'] = getSum_for_gXW6002(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"])
	data['special']['xw5902'] = getSum_for_XW5902(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"], jsonDict['hd'])
	data['special']['ad3101'] = var_ad3101_summary(data)
	data['special']['tXW6701'] = getSum_for_tXW6701(data)
	data['special']['gXW6701'] = getSum_for_gXW6701(data)
	data['special']['xw6003'] = getSum_for_XW6003(data, jsonDict['hd'])

	# BCL2L11����2���ں�����ϵȱʧ��̬�Լ���� 
	data["BCL2L11"] = ""
	for var in jsonDict["snvindel"]:
		if var["gene_symbol"] == "BCL2L11" and var["hgvs_c"] == "c.394+1479_394+4381del":
			data["BCL2L11"] = "T"
			break

	# θ�����ӷ��ͽ��
	data["GA_type"] = process_ga_type(jsonDict, var_data)

	# Master HD
	data["gss"] = getgss(jsonDict, var_data, config) 

	# ��ϸ��I/II��Snvindel����IGVͼ-20221020
	# ���ø�����һCP40ģ��
	data["igv_I_II"] = [var["igvplot"] for var in data["var_somatic"]["level_I"] + data["var_somatic"]["level_II"] if var["bio_category"] == "Snvindel" and "igvplot" in var.keys() and var["igvplot"]]

	#��ϵ4/5��Snvindel����IGVͼ-2022.11.18
	# ���ø�����һHRRģ��
	data["igv_4_5"] = [var["igvplot"] for var in data["var_germline"]["level_5"] + data["var_germline"]["level_4"] if var["bio_category"] == "Snvindel" and "igvplot" in var.keys() and var["igvplot"]]

	# ��ϸ������������չ��ر���snvindel IGVͼ-2023.06.27
	# �����人Э��HRR
	data["igv_onconodrug"] = [var["igvplot"] for var in data["var_somatic"]["level_onco_nodrug"] if var["bio_category"] == "Snvindel" and "igvplot" in var.keys() and var["igvplot"]]
	

	#### BRCA #####################################################################
	data_BRCA = {}
	data_BRCA["snv_s"] = {}
	data_BRCA["snv_m"] = {}
	#data_BRCA["mlpa"], data_BRCA["mlpa_image"], data_BRCA["mlpa_image_del"] = process_mlpa(jsonDict, config)
	data_BRCA["mlpa"] = mlpa
	data_BRCA["mlpa_image"] = mlpa
	data_BRCA["mlpa_image_del"] = mlpa_image_del
	data_BRCA["igvplot"] = {}
	data_BRCA["refer"] = {}
	
	for gene in ["BRCA1", "BRCA2"]:
		# snvindel
		for i in range(1, 6):
			# ������֯/����ѪҺ���ȼ��жϣ���ϵ��clinic_num_g; ��ϸ����clinic_num_s ���ܴ���L5�д�����ϵ�����²���L4�д���ϵ�²���ע��ʶ�𣡣�
			data_BRCA["snv_s"]["B"+gene[-1]+"_L"+str(i)] = [var for var in snvindel if var["gene_symbol"] == gene and ((var["var_origin"] == "germline" and var["clinic_num_g"] == i) or (var["var_origin"] != "germline" and var["clinic_num_s"] == i))]
			# ����ȡ��ϵ������������ģ������Ҫ���Ϊ���ݵ�ѪҺ����
			data_BRCA["snv_m"]["B"+gene[-1]+"_G_L"+str(i)] = [var for var in snvindel if var["gene_symbol"] == gene and var["var_origin"] == "germline" and var["clinic_num_g"] == i]

		# IGV
		data_BRCA["igvplot"]["all"] = [var["igvplot"] for var in snvindel if ((var["var_origin"] == "germline" and var["clinic_num_g"] in [4, 5]) or (var["var_origin"] != "germline" and var["clinic_num_s"] in [4, 5])) and var["gene_symbol"] in ["BRCA1","BRCA2"]]

		# ����ȡ��ϵ������������ģ������Ҫ���Ϊ���ݵ�ѪҺ����
		data_BRCA["igvplot"]["germline"] = [var["igvplot"] for var in snvindel if var["var_origin"] == "germline" and var["clinic_num_g"] in [4, 5] and var["gene_symbol"] in ["BRCA1","BRCA2"]]

		# ��ȡ�ο����ף���ͬģ��չʾ���ݲ�ͬ
		# Ŀǰjosn�вο������ǻ�����ʽ���޷�ֱ��ʹ�ã��ʽű�����������ȡ
		data_BRCA["refer"] = {}
		# 1. ��������PMID����
		data_BRCA["refer"]["var"] = []
		data_BRCA["refer"]["var_G"] = []
		# 2. �������PMID����
		data_BRCA["refer"]["gene"] = []
		data_BRCA["refer"]["gene_G"] = []
		# 3. ֤������PMID����
		data_BRCA["refer"]["evi"] = []
		data_BRCA["refer"]["evi_G"] = []
		for var in snvindel:
			## a. BRCA ��������Ŀ
			# ��������������ΪBRCA����-2023.05.06
			if ((var["var_origin"] == "germline" and var["clinic_num_g"] in [4, 5]) or (var["var_origin"] != "germline" and var["clinic_num_s"] in [4, 5])) and \
			var["gene_symbol"] in ["BRCA1", "BRCA2"]:
				data_BRCA["refer"]["var"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
				data_BRCA["refer"]["gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
				data_BRCA["refer"]["evi"].extend(var["evi_sum"]["refer_evi"])
			## b. BRCA ���������Ŀ
			if var["var_origin"] == "germline" and var["clinic_num_g"] in [4, 5] and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
				data_BRCA["refer"]["var_G"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
				data_BRCA["refer"]["gene_G"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
				data_BRCA["refer"]["evi_G"].extend(var["evi_sum"]["refer_evi"])
			
	## MLPA����©���ˣ�����һ��-��쿷�-2023.11.07
	for var in mlpa["B1_Loss"] + mlpa["B2_Loss"]:
		data_BRCA["refer"]["var"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		data_BRCA["refer"]["gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		data_BRCA["refer"]["evi"].extend(var["evi_sum"]["refer_evi"])
		data_BRCA["refer"]["var_G"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		data_BRCA["refer"]["gene_G"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		data_BRCA["refer"]["evi_G"].extend(var["evi_sum"]["refer_evi"])
	for var in mlpa["B1_Gain"] + mlpa["B2_Gain"]:
		data_BRCA["refer"]["var"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		data_BRCA["refer"]["gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		data_BRCA["refer"]["var_G"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		data_BRCA["refer"]["gene_G"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
	# MLPA����������-2023.11.07

	# ����ܽ�����������
	data_BRCA["special"] = {}	
	# ����-�㶫����С��
	data_BRCA["special"]["summary_GDRM_G"], data_BRCA["special"]["summary_GDRM_S"] = var_summary_GDRM(data_BRCA)
	# ����-����������4,5��ͻ��
	data_BRCA["special"]["GZZL_level13_G"], data_BRCA["special"]["GZZL_level13_S"] = var_GZZL(data_BRCA)
	# ����-�Ϻ��ʼã��������ƴ��
	data_BRCA["special"]["summary_SHRJ"] = var_summary_SHRJ(data_BRCA)
	# ����-����������BRCA1�����²��Ա��죩
	data_BRCA["special"]["summary_FJZL"] = var_summary_FJZL(data_BRCA)
	# ����-ɽ��ʡ��/ɽ����³��4��5��ͻ��
	data_BRCA["special"]["SD_level1_2_G"], data_BRCA["special"]["SD_level1_2_S"] = var_SD(data_BRCA)
	# ����-�人Э���ܽᣨ1,2,3���4,5�ࣩ
	data_BRCA["special"]["summary_WHXH_level45"], data_BRCA["special"]["summary_WHXH_level123"], data_BRCA["special"]["summary_WHXH_level123_G"] = var_summary_WHXH(data_BRCA)
	# ����-ɽ��������4,5��
	data_BRCA["special"]["summary_SDZL_level123_G"], data_BRCA["special"]["summary_SDZL_level123_S"] = var_level123_SDZL(data_BRCA)
	# ����-�ж���ϸ��������ߵȼ�
	data_BRCA["special"]["judge_intumor"] = brca_judge_inApprovalTumor(data_BRCA)
	# ����-���츽һ-��ԺȫѪ�ܽ�
	data_BRCA["special"]["sum_g_CQFY"] = var_summary_CQFY(data_BRCA)
	#print (data["special"]["sum_g_CQFY"])
	# ����-������һ-2022.11.11
	data_BRCA["special"]["summary_FJFY"] = BRCA_FJFY_summary(data_BRCA)
	

	### HRR for SHSY ######################################
	data_HRR_SHSY = {}
	# ����������򣬰�I/II��ͬ�ȼ�ʱ�����Ʒ�����ߵȼ��������򣩡��������͡�Ƶ�ʣ�Ƶ���ڸ����ִ����ʱ���Ѿ������ˣ�
	# ͬ�ȼ�ʱ�����Ʒ�����ߵȼ��������� A > B > C > D > N
	shsy_top_level_rule = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "N" : 4}
	shsy_var_data = sorted(snvindel, key=lambda i : (i["clinic_num_s"], shsy_top_level_rule.get(i["top_level"]), i["gene_symbol"]))
	brca_list = ["BRCA1", "BRCA2"]
	hrr_list = ["ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","FANCA","FANCL","MRE11A","NBN","PALB2","PPP2R2A","RAD51B","RAD51C","RAD51D","RAD54L"]
	other_list = ["AR","CDH1","ESR1","HDAC2","HOXB13","PTEN","STK11","TP53","BRAF","ERBB2","KRAS","NRAS","PIK3CA"]

	# ��ϸ���������
	data_HRR_SHSY["var_s_brca"] = s_var_rule_genelist(shsy_var_data, brca_list)
	data_HRR_SHSY["var_s_hrr"] = s_var_rule_genelist(shsy_var_data, hrr_list)
	data_HRR_SHSY["var_s_other"] = s_var_rule_genelist(shsy_var_data, other_list)

	# ��ϵ�������
	data_HRR_SHSY["var_g_brca"] = g_var_rule_genelist(shsy_var_data, brca_list)
	data_HRR_SHSY["var_g_hrr"] = g_var_rule_genelist(shsy_var_data, hrr_list)
	data_HRR_SHSY["var_g_other"] = g_var_rule_genelist(shsy_var_data, other_list)
	

	
	return data, data_BRCA, data_HRR_SHSY