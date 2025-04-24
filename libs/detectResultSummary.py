#-*- coding:gbk -*-
import os
import re
from libs.getConfig import CP40Gene
from libs.getConfig import CP40_SplitGene
import copy
from libs.rule import S_level, s_var_rule
from libs.getConfig import PAN116_other_gene

'''
Discription
	
	��������ܣ���ʽ���֡� 
	
'''
def getSummary_detect(var_data, mlpa_data, config):
	summary_result = {}
	summary_result["format1_forCP"] = format1(var_data, config)
	summary_result["format2_forCP_splitGene"] = {}
	summary_result["format2_forCP_splitGene"]["gene10"], \
	summary_result["format2_forCP_splitGene"]["gene30"], \
	summary_result["format2_forCP_splitGene"]["gene_116_other"] = format2(var_data, config)
	summary_result["format3_forHRDC_QL"] = format3(var_data, mlpa_data)
	summary_result["format4_forHRR_QL"] = format4(var_data, mlpa_data)
	summary_result["format5_for116"] = {}
	summary_result["format5_forLC10"] = {}
	summary_result["format5_for116"]["CRC25"], \
	summary_result["format5_for116"]["TC21"], \
	summary_result["format5_for116"]["GA18"], \
	summary_result["format5_forLC10"]["NBBL"] = format5(var_data)
	summary_result["format6_forLC10_YNZL"] = format6(var_data)
	summary_result["format7_forHRR_YNZL"] = format7(var_data)
	summary_result["format8_forHRDC_YNZL"] = format8(var_data)
	summary_result["format9_forLC10_FJZL"] = format9(var_data)
	summary_result["format10_forYNZL"] = {}
	summary_result["format10_forYNZL"]["CP40"], \
	summary_result["format10_forYNZL"]["TC21"], \
	summary_result["format10_forYNZL"]["GA18"], \
	summary_result["format10_forYNZL"]["CRC12"] = format10(var_data)	
	summary_result["format11_forFDZS_CP40"], summary_result["format11_forFDZS_CP40_sum"] = format11(var_data)
	summary_result["format12_forGDRM_gHRR"] = format12(var_data, mlpa_data)
	summary_result["format13_forBRCA_YNZL"] = format13(var_data, mlpa_data)
	summary_result["format14_forBPTM_YNZL"] = format14(var_data)

	return summary_result

# ����1���������򡢼�����ݺͼ�������������Ͱ���Snvindel��CNV��SV������I/II/III�����
# ��Ҫ����CP40��Ʒ
def getSum_forCP40(var_data_copy,gene_dict):
	gene_list = [k for k in gene_dict.keys()]
	result = []
	var_list = [var for var in var_data_copy if (var["var_origin"] == "germline" and var["clinic_num_g"] in [4,5]) or \
				(var["var_origin"] != "germline" and var["clinic_num_s"] in [3,4,5])]
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_dict.keys():
				var_info = ""
				note = var["judge_mergeMET"] if "judge_mergeMET" in var.keys() else ""
				note_freq = var["freq_2"] if "judge_mergeMET" in var.keys() else ""
				freq = ""
				if var["bio_category"] == "Snvindel":
					freq = var["freq_str"]
				elif var["bio_category"] == "Cnv":
					freq = var["cn_mean"]
				else:
					if "freq_str" in var.keys() and var["freq_str"]:
						freq = var["freq_str"]
					elif "copies" in var.keys() and var["copies"]:
						freq = var["copies"]
					elif "reads" in var.keys() and var["reads"]:
						freq = var["reads"]
				level = str(int(S_level(var)))
				result.append({
					"gene_symbol" : gene,
					"var_info" : var["five_prime_gene"]+"-"+var["three_prime_gene"]+"�ں�" if var["bio_category"] in ["Sv", "PSeqRnaSv"] else "����" if var["bio_category"] == "Cnv" else var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"],
					"freq" : freq,
					"region" : gene_dict[gene],
					"level" : level,
					"note" : note,
					"note_freq" : note_freq,
					"var_origin" : var["var_origin"],
					"freq_rc" : var["freq_rc"] if "freq_rc" in var.keys() and var["freq_rc"] else "",
					"clinic_num_g" : var["clinic_num_g"] if "clinic_num_g" in var.keys() and var["clinic_num_g"] else "",
					"var_detail" : var
				})
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"region" : gene_dict[gene]
		})
	return sorted(result, key=lambda i:i["gene_symbol"])

# ��ʽ1��������CP40����ɽ����³/ɽ��ʡ��CP40
# �������򡢼�����ݺͼ������ÿ������һ�У�����I/II/III�����
def format1(var_data, config):
	var_data_copy = copy.deepcopy(var_data)
	cp40_dict = CP40Gene(config)
	cp40_summary = getSum_forCP40(var_data_copy, cp40_dict)
	return cp40_summary

# ��ʽ2�����10�����30�����������ĳǡ�������һCP40
# �������򡢼�����ݺͼ������ÿ������һ�У�����I/II/III�����
def format2(var_data, config):
	var_data_copy = copy.deepcopy(var_data)
	gene10_dict = CP40_SplitGene(config)["gene10"]
	gene30_dict = CP40_SplitGene(config)["gene30"]
	gene116_other = PAN116_other_gene(config)
	gene10_summary = getSum_forCP40(var_data_copy, gene10_dict)
	gene30_summary = getSum_forCP40(var_data_copy, gene30_dict)
	gene_116_other_summary = getSum_forCP40(var_data_copy, gene116_other)
	return gene10_summary, gene30_summary, gene_116_other_summary

# ����2���������򡢼�����ݡ���Ⱥͱ�������ÿ������һ�У���ϸ��/��Դ�����������I/II/III�࣬��ϵ����345��
# �����Ϊ�ѻ�������չʾʱ��Ҫ���Ǻţ�����������
# ��Ҫ������ϵ��Ʒ����HRR��HRD C�ȣ�
def getSum_forGermlinePanel(var_data, mlpa_data, gene_list_all, gene_list_appr):
	var_data_copy = copy.deepcopy(var_data)
	mlpa_data_copy = copy.deepcopy(mlpa_data)
	result = []
	var_list = [var for var in var_data_copy if (var["var_origin"] == "germline" and var["clinic_num_g"] in [3,4,5]) or\
				(var["var_origin"] != "germline" and var["clinic_num_s"] in [3,4,5])]
	detect_gene = [var["gene_symbol"] for var in var_list]
	# ����MLPA
	for k, v in mlpa_data_copy.items():
		if v:
			var_list += v
			detect_gene += [var["gene_symbol"] for var in v]
	for var in var_list:
		# ���Ʊ�������Ϊsnvindel����Ȼ�ᱨ��
		if var["gene_symbol"] in gene_list_all and var["bio_category"] == "Snvindel":
			# �²�/�����²�������ҩ�Ĺ�ΪIII��
			var["clinic_num_s"] = S_level(var)
			result.append(var)
		else:
			if "type" in var.keys() and var["type"] and re.search("Loss|Gain", var["type"]) and var["gene_symbol"] in gene_list_all:
				result.append(var)
	for gene in set(gene_list_all) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	# ��ǻ�������
	for i in result:
		if i["gene_symbol"] in  gene_list_appr:
			i["appr"] = "T"

	return sorted(result, key = lambda i:i["gene_symbol"])

# ��ʽ3: ������HRD C��ɽ����³
# �������򡢼�����ݡ���Ⱥͱ�������ÿ������һ�У�����I/II/III�����
# summary�����з�Ϊ�ѻ�������չʾʱ��Ҫ���Ǻţ�����������
def format3(var_data, mlpa_data):
	gene_list_all = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", "FANCA", "FANCL", "HDAC2", "PALB2", "PPP2R2A", \
					"PTEN", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "TP53"]
	gene_list_appr = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PALB2", "RAD51B", "RAD51C", "RAD51D", "RAD54L"]
	return getSum_forGermlinePanel(var_data, mlpa_data, gene_list_all, gene_list_appr)

# ��ʽ4: ������HRR��ɽ����³ HRRȫѪ����֯����ԣ�������ϸ��I/II/III��������ϵ3/4/5�����
# �������򡢼�����ݡ���Ⱥͱ�������ÿ������һ��
# summary�����з�Ϊ�ѻ�������չʾʱ��Ҫ���Ǻţ�����������
def format4(var_data, mlpa_data):
	gene_list_all = ["BRCA1", "BRCA2", "AR", "ATM", "ATR", "BARD1", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", "ESR1", \
					 "FANCA", "FANCL", "HDAC2", "HOXB13", "MRE11A", "NBN", "PALB2", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", \
					 "RAD51D", "RAD54L", "STK11", "TP53", "BRAF", "ERBB2", "KRAS", "NRAS", "PIK3CA"]
	gene_list_appr = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PALB2", "RAD51B", "RAD51C", "RAD51D", "RAD54L"]
	return getSum_forGermlinePanel(var_data, mlpa_data, gene_list_all, gene_list_appr)

# ����3��������116ϵ�в�Ʒ��CRC25��֯/��֯���/ѪҺ��ԡ�TC21��֯/��֯���/ѪҺ��ԡ�GA18��֯/��֯���/ѪҺ��ԣ�����������
# չʾ��ϸ��/��Դ���� I/II/III/����������չ��ر��죬��ϵ45�����
def getSum_forSomaticPanel(var_data, gene_list):
	var_data_copy = copy.deepcopy(var_data)
	result = []
	var_list = [var for var in var_data_copy if (var["var_origin"] == "germline" and var["clinic_num_g"] in [4,5]) or \
				(var["var_origin"] != "germline" and var["clinic_num_s"] in [3,4,5])]
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_list:
				# �²�/�����²�������ҩ�Ĺ�ΪIII��
				var["clinic_num_s"] = S_level(var)
				result.append(var)
	for gene in set(gene_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	return sorted(result, key = lambda i:i["gene_symbol"])

#��ʽ5��������116ϵ�в�Ʒ��CRC25��֯/��֯���/ѪҺ��ԡ�TC21��֯/��֯���/ѪҺ��ԡ�GA18��֯/��֯���/ѪҺ��ԣ�����������
def format5(var_data):
	crc25_list = ["AKT1", "ARAF", "BRAF", "EGFR", "EPCAM", "ERBB2", "FBXW7", "GNAS", "HRAS", "KDR", "KRAS", \
				  "MAP2K1", "MET", "MLH1", "MSH2", "MSH6", "NRAS", "NTRK1", "PIK3CA", "PMS2", "POLE", "PTEN", "SMAD4"]
	tc21_list = ["AKT1", "ALK", "BRAF", "CTNNB1", "EIF1AX", "GNAS", "HRAS", "KRAS", "NRAS", "NTRK1", "NTRK3", \
				 "PAX8", "PDGFRA", "PIK3CA", "PTEN", "RASAL1", "RET", "TERT", "TP53", "TSC2", "TSHR"]
	ga18_list = ["AKT1", "BRAF", "EGFR", "ERBB2", "FGF19", "FGFR1", "FGFR2", "FGFR3", "HRAS", "KIT", "KRAS", "MET", \
				 "MTHFR", "NRAS", "PDGFRA", "PIK3CA", "PTEN", "TP53"]
	lc10_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "RET", "ROS1", "NRAS", "PIK3CA"]
	CRC25_result = getSum_forSomaticPanel(var_data, crc25_list)
	TC21_result = getSum_forSomaticPanel(var_data, tc21_list)
	GA18_result = getSum_forSomaticPanel(var_data, ga18_list)
	LC10_result = getSum_forSomaticPanel(var_data, lc10_list)
	return CRC25_result, TC21_result, GA18_result, LC10_result

#��ʽ6����������������LC10
#ROS1/RET/ALK��I/II���ںϣ�MET��Ծͻ����MET I����죬����ƥ��hgvs_p
def format6(var_data):
	sv_gene = ["ALK", "ROS1", "RET"]
	spec_gene = ["MET"]
	snvindel_gene = {
		"BRAF" : ["V600E"],
		"EGFR" : ["E746_A750del", "L747_T751del", "E746_S752delinsV", "L747_P753delinsS", "L747_A750delinsP", \
				  "D770_N771insG", "L858R", "L861Q","G719A","T790M"],
		"ERBB2" : ["Y772_A775dup"],
		"KRAS" : ["G12S", "G12D", "G12A", "G12V", "G12C"],
		"NRAS" : ["G12D", "Q61R"],
		"PIK3CA" : ["H1047R"]
	}
	result = {}
	# �����ںϱ���, I/II�࣬������������չ��ΪIII�࣬�˴���չʾ����ע����������ܻ�Ϊ��gene1,gene2��
	for gene in sv_gene:
		if gene not in result.keys():
			result.setdefault(gene+"_sv", [])
		result[gene+"_sv"] = [var["var_hgvs"] for var in var_data if var["clinic_num_s"] in [4,5] and var["var_origin"] != "germline" and \
							  "evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and \
							  set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and \
							  var["bio_category"] == "Sv" and gene in re.split(",", var["gene_symbol"])]	
	# ����MET 14��Ծͻ�䣬I��
	result.setdefault("MET_14", [])
	result["MET_14"] = [(var["hgvs_c"],var["hgvs_p"]) for var in var_data if var["clinic_num_s"] == 5 and \
						var["var_origin"] != "germline" and "evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and \
						set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and \
						var["bio_category"] == "Snvindel" and var["gene_symbol"] == "MET"]
	# ����������ͻ��/����ȱʧͻ��
	for gene, info in snvindel_gene.items():
		for hgvs_p in info:
			if gene+"_"+hgvs_p not in result.keys():
				result.setdefault(gene+"_"+hgvs_p, [])
			result[gene+"_"+hgvs_p] = [var["hgvs_p"] for var in var_data if var["gene_symbol"] == gene and \
									   var["bio_category"] == "Snvindel" and \
									   hgvs_p == var["hgvs_p_ZJZL"].replace("p.", "")]

	return result

#��ʽ7����������������HRR
# ��ϵ4/5�࣬BRCA MLPA del��ģ�����ж��ͺã����ֻ����snvindel
# ��ϸ��I/II��
def format7(var_data):
	result = {}
	gene_list = ["AR","ATM","BRCA1","BRCA2","CDK12","CHEK1","CHEK2","ERBB2","FANCA","FANCL","KRAS","MRE11A",\
				 "NRAS","PALB2","PIK3CA","PTEN","RAD51B","RAD51C","RAD51D","RAD54L","STK11","TP53"]
	for gene in gene_list:
		if gene not in result.keys():
			result.setdefault(gene, [])
		result[gene] = [var["hgvs_c"] for var in var_data if ((var["clinic_num_s"] in [4,5] and \
						var["var_origin"] != "germline" and "evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and \
						set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys())) or \
						(var["clinic_num_g"] in [4,5] and var["var_origin"] == "germline")) and \
						gene in re.split(",", var["gene_symbol"]) and var["bio_category"] == "Snvindel"]
	return result

#��ʽ8����������������HRD C
# ��չʾI/II�����
# BRCA1/BRCA2��Ҫ���־��������Ӻţ�����������
# BRCA Exon2-3����չʾexon2 + intron2 + exon3�ı���
# BRCA Exon4����չʾexon4�ı���
def format8(var_data):
	gene_region = ["BRCA1", "BRCA2"]
	other_gene = ["ATM","BARD1","BRIP1","CDH1","CDK12","CHEK1","CHEK2","FANCA","FANCL","HDAC2","PALB2","PPP2R2A",\
				  "PTEN","RAD51B","RAD51C","RAD51D","RAD54L","TP53"]
	var_data_copy = copy.deepcopy(var_data)
	var_list = [var for var in var_data_copy if var["clinic_num_s"] in [4,5] and \
				var["var_origin"] != "germline" and \
				"evi_sum" in var.keys() and \
				var["evi_sum"]["evi_split"] and \
				set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and \
				var["bio_category"] == "Snvindel"]
	result = {}
	# ����BRCA����
	for var in var_list:
		if var["gene_symbol"] in gene_region:
			for region in re.split("_", var["gene_region"]):
				region = region.replace("5'UTR", "5UTR").replace("3'UTR", "3UTR")
				if var["gene_symbol"]+"_"+region not in result.keys():
					result.setdefault(var["gene_symbol"]+"_"+region, [])
				result[var["gene_symbol"]+"_"+region].append((var["hgvs_c"], var["hgvs_p"]))
	# ���������������
	for gene in other_gene:
		if gene not in result.keys():
			result.setdefault(gene, [])
		result[gene] = [(var["hgvs_c"],var["hgvs_p"]) for var in var_list if var["gene_symbol"] == gene]
	return result

#��ʽ9�������ڸ�������LC10����������-2022.11.03
#��һ�Σ�չʾ�������
#�ڶ��Σ����ñ�֤�ݽ���
#�����Σ�������
##ע�⣺ǰ���������췢������ͬ�����ϣ����һ����дһ�飬�����κϲ�д��
#���ĶΣ�ҩ���ܽᣬ���б��������Ϣ�����ܵ����һ��
def format9(var_data):
	var_level = s_var_rule(var_data)
	### ����ǰ��������
	# ���ݻ��򽫱�����з��ദ��
	var_gene_group = {}
	for var in var_level["level_I"] + var_level["level_II"] + var_level["level_onco_nodrug"] + var_level["level_III"]:
		for gene in re.split(",", var["gene_symbol"]):
			if gene not in var_gene_group.keys():
				var_gene_group.setdefault(gene, [])
			var_gene_group[gene].append(var)
	### ���������
	def stran_var_inter(i):
		inter = ""
		if i["bio_category"] == "Snvindel":
			if i["hgvs_p"] != "p.?":
				inter = i["gene_symbol"]+"�����"+i["hgvs_c"]+" "+i["hgvs_p"]+i["variant_desc_cn"].replace("�ñ���","")+str(i["variant_interpret_cn"])
			else:
				inter = i["gene_symbol"]+"�����"+i["hgvs_c"]+i["variant_desc_cn"].replace("�ñ���","")+str(i["variant_interpret_cn"])
		elif i["bio_category"] == "Cnv":
			inter = "����ʵ���⵽"+i["gene_symbol"]+"������"+str(i["variant_interpret_cn"])
		elif i["bio_category"] == "Sv":
			if i["variant_interpret_cn"]:
				inter = "����ʵ���⵽"+i["five_prime_gene"]+"-"+i["three_prime_gene"]+"�ںϣ��ں�����Ϊ"+i["five_prime_gene"]+":"+i["five_prime_cds"]+\
						"-"+i["three_prime_gene"]+":"+i["three_prime_cds"]+"��"+i["variant_interpret_cn"]
			else:
				inter = "����ʵ���⵽"+i["five_prime_gene"]+"-"+i["three_prime_gene"]+"�ںϣ��ں�����Ϊ"+i["five_prime_gene"]+":"+i["five_prime_cds"]+\
						"-"+i["three_prime_gene"]+":"+i["three_prime_cds"]+"��"
		return inter
	
	### ����������
	def stran_gene_inter(i):
		inter = {}
		if not re.search(",", i["gene_symbol"]):
			inter[i["gene_symbol"]] = i["gene_function"]
		else:
			if i["bio_category"] == "Sv" or i["bio_category"] == "PSeqRnaSv":
				inter[i["five_prime_gene"]] = i["five_prime_gene_function"]
				inter[i["three_prime_gene"]] = i["three_prime_gene_function"]
		
		return inter

	# ǰ������Ϣ����
	result = {}
	result["var"] = []
	for gene, var in var_gene_group.items():
		var_info = []
		gene_info = {}
		rpt_info = {}
		for i in var:
			var_info.append(stran_var_inter(i))
			gene_info.update(stran_gene_inter(i))
			if i["var_info_forFJZL"]:
				rpt_info.update(i["var_info_forFJZL"])
		result["var"].append(
			{
				"gene_info" : gene_info.get(gene),
				"rpt_info" : rpt_info.get(gene),
				"var_info" : "".join(var_info)
			}
		)

	## �������һ������
	# ����ҩ���ܽ��еı�������
	def stran_var_for_drug(i):
		inter = ""
		if i["bio_category"] == "Snvindel":
			type_cn = i["type_cn"] if i["type_cn"] != "--" else i["type"]
			type_cn = type_cn.replace("Intronic", "�ں���ͻ��").replace("3'UTR","3'UTRͻ��").replace("5'UTR", "5'UTRͻ��")
			if i["hgvs_p"] != "p.?":
				inter = "���"+i["gene_symbol"]+"�����"+i["hgvs_c"]+" "+i["hgvs_p"]+type_cn
			else:
				inter = "���"+i["gene_symbol"]+"�����"+i["hgvs_c"]+type_cn
		elif i["bio_category"] == "Cnv":
			inter = "���"+i["gene_symbol"]+"����"
		elif i["bio_category"] == "Sv":
			inter = "���"+i["five_prime_gene"]+"-"+i["three_prime_gene"]+"�ں�"

		return inter

	# ����ҩ���ܽ�
	result["drug"] = ""
	regimen_sum = []
	for var in var_level["level_I"] + var_level["level_II"]:
		var_regimen_S = []
		var_regimen_R = []
		var_regimen = []
		for i in ["A","B","C","D"]:
			if var["var_regimen_forFJZL"]["regimen_"+str(i)+"_S"]:
				var_regimen_S.append("��".join(var["var_regimen_forFJZL"]["regimen_"+str(i)+"_S"])+"��֤�ݵȼ�Ϊ"+str(i)+"��")
			if var["var_regimen_forFJZL"]["regimen_"+str(i)+"_R"]:
				var_regimen_R.append("��".join(var["var_regimen_forFJZL"]["regimen_"+str(i)+"_R"])+"��֤�ݵȼ�Ϊ"+str(i)+"��")

		if var_regimen_S:
			var_regimen.append("���ܶ�"+"��".join(var_regimen_S)+"����")
		if var_regimen_R:
			var_regimen.append("���ܶ�"+"��".join(var_regimen_R)+"��ҩ")
		regimen_sum.append(stran_var_for_drug(var)+"����ͻ����ʾ"+"��".join(var_regimen)+"��")
	result["drug"] = "".join(regimen_sum)
	return result

#��ʽ10����������������CP40/CRC12/TC21/GA18
# CP40չʾI/II/III������������չ��ر���
# CRC12չʾI/II����죬����Ҫ��ϸ�����ͣ���C12
# TC21/GA18չʾI/II����죨����������������չ��أ�
def format10(var_data):
	var_data_copy = copy.deepcopy(var_data)
	CP40_gene = ["AKT1","ALK","BRAF","CDK4","CTNNB1","DDR2","EGFR","ERBB2","FGFR1","FGFR2","FGFR3","KIT","KRAS","MAP2K1",\
				 "MET","MYC","NRAS","NTRK1","NTRK2","NTRK3","PDGFRA","PIK3CA","PTEN","POLE","RB1","RET","ROS1","TP53"]
	TC21_gene = ["BRAF","KRAS","NRAS","HRAS","GNAS","TERT","RET","ALK","NTRK1","PAX8","PIK3CA","AKT1","EIF1AX","TSC2",\
				 "CTNNB1","PDGFRA","RASAL1","TP53","PTEN","TSHR","NTRK3"]
	GA18_gene = ["AKT1","BRAF","EGFR","ERBB2","FGF19","FGFR1","FGFR2","FGFR3","HRAS","KIT","KRAS","MET","MTHFR",\
				 "NRAS","PDGFRA","PIK3CA","PTEN","TP53"]
	CP40_var_list = [var for var in var_data_copy if var["var_origin"] != "germline" and var["clinic_num_g"] in [3,4,5]]
	Other_var_list = [var for var in var_data_copy if (var["var_origin"] == "germline" and var["clinic_num_g"] in [4,5]) or\
					 (var["var_origin"] != "germline" and var["clinic_num_s"] in [4,5] and \
					  "evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and \
					  set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()))]
	# ����CP40��TC21��GA18
	def summary(var_list, gene_list):
		result = {}
		for gene in gene_list:
			if gene not in result.keys():
				result.setdefault(gene, [])
			result[gene] = [var for var in var_list if gene in re.split(",", var["gene_symbol"])]
		return result

	CP40_result = summary(CP40_var_list, CP40_gene)
	TC21_result = summary(Other_var_list, TC21_gene)
	GA18_result = summary(Other_var_list, GA18_gene)

	# CRC12-��ȫƥ�䰱����仯
	CRC12_snvindel_complete = {
		"BRAF" : ["V600E"],
		"POLE" : ["P286R", "V411L", "S297F"]
	}
	# CRC12-ƥ�䰱����仯��ͷ
	CRC12_snvindel_simple = {
		"KRAS" : ["G12", "G13", "Q61", "A146"],
		"NRAS" : ["G12", "G13", "Q61", "A146"],
		"PIK3CA" : ["E542", "E545", "H1047"]
		}
	# CRC12-�����������б���ʽ�������������
	CRC12_cnv = ["ERBB2"]
	# CRC12-�ں�
	CRC12_sv = ["NTRK1", "NTRK3", "ALK", "RET", "FGFR3"]

	def summary_CRC12(var_list):
		result = {}
		# �����ںϱ���, I/II�࣬ע����������ܻ�Ϊ��gene1,gene2��
		for gene in CRC12_sv:
			if gene not in result.keys():
				result.setdefault(gene+"_sv", [])
			result[gene+"_sv"] = [var["var_hgvs"] for var in var_list if var["bio_category"] == "Sv" and gene in re.split(",", var["gene_symbol"])]
		# ��������
		for gene in CRC12_cnv:
			if gene not in result.keys():
				result.setdefault(gene+"_cnv", [])
			result[gene+"_cnv"] = ["����" for var in var_list if var["bio_category"] == "Cnv" and gene == var["gene_symbol"]]
		# ����snvindel�� ��ȫƥ��hgvs_p
		for gene, hgvs_p_list in CRC12_snvindel_complete.items():
			for hgvs_p in hgvs_p_list:
				if gene+"_"+hgvs_p not in result.keys():
					result.setdefault(gene+"_"+hgvs_p, [])
				result[gene+"_"+hgvs_p] = [var["hgvs_p"] for var in var_list if var["bio_category"] == "Snvindel" and \
										   var["gene_symbol"] == gene and hgvs_p == (var["hgvs_p_ZJZL"].replace("p.", ""))]
		# ����snvindel������ȫƥ��hgvs_p
		for gene, hgvs_p_list in CRC12_snvindel_simple.items():
			for hgvs_p in hgvs_p_list:
				if gene+"_"+hgvs_p not in result.keys():
					result.setdefault(gene+"_"+hgvs_p, [])
				result[gene+"_"+hgvs_p] = [var["hgvs_p"] for var in var_list if var["bio_category"] == "Snvindel" and \
										   var["gene_symbol"] == gene and re.search(hgvs_p, var["hgvs_p"])]
		return result
	CRC12_result = summary_CRC12(Other_var_list)

	return CP40_result, TC21_result, GA18_result, CRC12_result

# ��ʽ11�������ڸ�����ɽCP40�����������ֵ�������ˣ��³����⡣-2022.11.28
# 1. �������б�
# 2. �ܽ᣺��⵽XX��XX����ͻ�䡣
def format11(var_data):
	var_data_copy = copy.deepcopy(var_data)
	cp40_list = ["AKT1","ALK","BRAF","CDK4","CTNNB1","DDR2","EGFR","ERBB2","ESR1","FGFR1",\
				"FGFR2","FGFR3","FGFR4","HRAS","IDH1","IDH2","KEAP1","KIT","KRAS","MAP2K1",\
				"MET","MYC","NFE2L2","NKX2-1","NRAS","NRG1","NTRK1","NTRK2","NTRK3","PDGFRA",\
				"PIK3CA","POLE","PTEN","RB1","RET","ROS1","STK11","TP53"]
	result = []
	detect_gene = []
	var_list = [var for var in var_data_copy if var["var_origin"] != "germline" and var["clinic_num_s"] in [3,4,5]]
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in cp40_list:
				detect_gene.append(gene)
				var["clinic_num_s"] = S_level(var)
				var["gene_symbol"] = gene
				result.append(var)
	for gene in set(cp40_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})	
	return sorted(result, key = lambda i:i["gene_symbol"]), "��".join(sorted(set(detect_gene)))

# ��ʽ12�������ڹ㶫����gHRR��չʾ��ϵ5/4/3�����
def format12(var_data, mlpa_data):
	var_data_copy = copy.deepcopy(var_data)
	mlpa_data_copy = copy.deepcopy(mlpa_data)
	gene_list = ["BRCA1","BRCA2","AR","ATM","ATR","BARD1","BRIP1","CDH1","CDK12","CHEK1",\
				 "CHEK2","ESR1","FANCA","FANCL","HDAC2","HOXB13","MRE11","NBN","PALB2","PPP2R2A",\
				 "PTEN","RAD51B","RAD51C","RAD51D","RAD54L","STK11","TP53","BRAF","ERBB2","KRAS",\
				 "NRAS","PIK3CA"]
	var_list = [var for var in var_data_copy if var["var_origin"] == "germline" \
												and var["clinic_num_g"] in [3,4,5] \
												and var["gene_symbol"] in gene_list \
												and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for k, v in mlpa_data_copy.items():
		if v:
			var_list += v
			detect_gene+= [var["gene_symbol"] for var in v]

	for gene in set(gene_list) - set(detect_gene):
		var_list.append({
			"gene_symbol" : gene
		})
	
	return sorted(var_list, key=lambda i:gene_list.index(i["gene_symbol"]))

#��ʽ13����������������BRCA
# ��չʾI/II/III�����ϵ3/4/5�����
# BRCA1/BRCA2��Ҫ���־��������Ӻ�
# BRCA Exon2-3����չʾexon2 + intron2 + exon3�ı���
# BRCA Exon4����չʾexon4�ı���
def format13(var_data, mlpa_data):
	var_data_copy = copy.deepcopy(var_data)
	var_list = [var for var in var_data_copy if ((var["clinic_num_s"] in [3,4,5] and var["var_origin"] != "germline") or \
												(var["clinic_num_g"] in [3,4,5] and var["var_origin"] == "germline")) and \
												var["bio_category"] == "Snvindel" and var["gene_symbol"] in ["BRCA1", "BRCA2"]]
	result = {}
	for var in var_list:
		for region in re.split("_", var["gene_region"]):
			region = region.replace("5'UTR", "5UTR").replace("3'UTR", "3UTR")
			if var["gene_symbol"]+"_"+region not in result.keys():
				result.setdefault(var["gene_symbol"]+"_"+region, [])
			result[var["gene_symbol"]+"_"+region].append((var["hgvs_c"], var["hgvs_p"]))
	
	# MLPA���
	# value�еĸ�ʽΪ��exon1a, exon1b, exon1-exon3, exon5-exon16
	mlpa_data_copy = copy.deepcopy(mlpa_data["B1_Loss"]+mlpa_data["B2_Loss"]+mlpa_data["B1_Gain"]+mlpa_data["B2_Gain"])
	for var in mlpa_data_copy:
		for i in re.split(",", var["value"].replace(" ","")):
			# ����Ϊ[1, 1, 1-3, 5-16]
			i = i.replace("exon","").replace("a","").replace("b", "").replace("up", "")
			exon_num = re.split("-", i)
			# �漰������ֻ��1��ʱ
			if len(exon_num) == 1:
				if var["gene_symbol"]+"_exon"+str(exon_num[0]) not in result.keys():
					result.setdefault(var["gene_symbol"]+"_exon"+str(exon_num[0]), [])
				result[var["gene_symbol"]+"_exon"+str(exon_num[0])].append("mlpa")
			# �漰���������ʱ
			elif len(exon_num) == 2:
				for num in range(int(exon_num[0]), int(exon_num[1])+1):
					if var["gene_symbol"]+"_exon"+str(num) not in result.keys():
						result.setdefault(var["gene_symbol"]+"_exon"+str(num), [])
					result[var["gene_symbol"]+"_exon"+str(num)].append("mlpa")
			else:
				pass

	return result

#��ʽ14����������������BPTM��������POLE��TP53��BRCA�����ø�ʽ13��
# ��չʾI/II/III��POLE��TP53
# ��Ҫ���־��������Ӻ�
# BRCA Exon2-3����չʾexon2 + intron2 + exon3�ı���
# BRCA Exon4����չʾexon4�ı���
def format14(var_data):
	var_data_copy = copy.deepcopy(var_data)
	var_list = [var for var in var_data_copy if ((var["clinic_num_s"] in [3,4,5] and var["var_origin"] != "germline") or \
												(var["clinic_num_g"] in [3,4,5] and var["var_origin"] == "germline")) and \
												var["bio_category"] == "Snvindel" and var["gene_symbol"] in ["POLE", "TP53"]]
	result = {}
	for var in var_list:
		for region in re.split("_", var["gene_region"]):
			region = region.replace("5'UTR", "5UTR").replace("3'UTR", "3UTR")
			if var["gene_symbol"]+"_"+region not in result.keys():
				result.setdefault(var["gene_symbol"]+"_"+region, [])
			result[var["gene_symbol"]+"_"+region].append((var["hgvs_c"], var["hgvs_p"]))
	
	return result