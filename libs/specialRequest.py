#-*- coding:gbk -*-
import re
from libs.AAstrans import strans
from libs.rule import s_var_rule
#from libs.getConfig import get_FDZS_database
#from libs.getConfig import get_FJZL_database
from libs.getConfig import getconfigxlsx, get_shrj_geneinfo
from libs.rule import S_level
import copy
from collections import defaultdict
'''
	����ģ���и�������Ҫ�� 
'''

# �ú����ɸ���
def varInter_FJZL(hgvs_p, var_type, gene_region, hgvs_c):
	'''
	����ʡ����ҽԺ��ֻ�е�����֯�͵���ȫѪ��
	1. �����������û�������Ѫ������DNA��Ʒ����BRCA2���� < 14�������ӷ���c.7409_7410insT����ͻ�䵼�»�����뵰�׵�2471λ���������հ���ͻ��Ϊ�鰱�Ტ��2474λ������ǰ��ֹ��
				�����γɹ������˻�ʧ��ĵ��ס�>
	< >�ڵ������ɸú���ʵ��
	'''
	#ת��gene_region����exonת��Ϊ12��������
	if re.search("exon", gene_region):
		gene_region_cn = gene_region.replace("exon", "")+"��������"
	elif re.search("intron", gene_region):
		gene_region_cn = gene_region.replace("intron", "")+"���ں���"
	else:
		gene_region_cn = gene_region

	inter1 = gene_region_cn+"����"+hgvs_c
	inter2 = strans(hgvs_p, var_type)
	inter = inter1+"��"+inter2 if inter2 else inter1
	return inter

# �ú������ڱ���ģ����ʵ��
def var_summary_FJZL(var_brca):
	'''
	����ʡ����ҽԺ��ֻ�е�����֯�͵���ȫѪ��
	1. �û�������Ѫ������DNA��Ʒ����< BRCA1�����²��Ա��� >
	'''	
	# ͳ���²��ԡ����²��ԡ�BRCA1��BRCA
	gene_rule = ["BRCA1", "BRCA2"]
	gene_list = [i["gene_symbol"] for i in var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"]+var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]+var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"]]
	level_rule = ["�����²��Ա���", "�²��Ա���", "���岻��ȷ"]
	level_dict = {4 : "�����²��Ա���", 5 : "�²��Ա���", 3 : "���岻��ȷ"}
	level_list = [level_dict.get(i["clinic_num_g"]) for i in var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"]+var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]]
	if var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
		level_list.append(level_dict.get(4))
	gene_str = "��".join(sorted(list(set(gene_list)), key=gene_rule.index))
	level_str = "��".join(sorted(list(set(level_list)), key=level_rule.index))

	return gene_str+"����"+level_str

# �ú������ڱ���ģ����ʵ��
def var_summary_GDRM(var_brca):
	'''
	�㶫����ֻ�е�����֯�͵���ȫѪ��
	ȫѪ��
		���ԣ����ڱ����ͼ����������μ����⵽X���²��Ա����X�������²��Ա���
		���ԣ����ڱ����ͼ�������δ��⵽�²��Ի������²��Ա���
	��֯��
		���ԣ����ڱ����ͼ����������μ����⵽X�������ٴ�����ı���
		���ԣ����ڱ����ͼ�������δ��⵽�����ٴ�����ı���
	'''
	level45_num_list_G = []
	level5_num = len(var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"])
	level4_num = len(var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]+var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"])
	if level5_num > 0:
		level45_num_list_G.append(str(level5_num)+"���²��Ա���")
	if level4_num > 0:
		level45_num_list_G.append(str(level4_num)+"�������²��Ա���")
	
	summary_G = "���ڱ����ͼ����������μ����⵽"+"��".join(level45_num_list_G) if level45_num_list_G else "���ڱ����ͼ�������δ��⵽�²��Ի������²��Ա���"
	summary_S = "���ڱ����ͼ����������μ����⵽"+str(level5_num+level4_num)+"�������ٴ�����ı���" if level5_num+level4_num > 0 else "���ڱ����ͼ�������δ��⵽�����ٴ�����ı���"
	
	return summary_G, summary_S
	
# ���ڱ���ģ����ʵ��
def var_GZZL(var_brca):
	'''
	����������ֻ�е�����֯�͵���ȫѪ��
	ȫѪ��
		����ʵ���⵽һЩ�����岻��ȷ���졢�������Ա��졢���Ա��죩
	��֯��
		����ʵ���⵽һЩ���ٴ����岻�����������Ա��졢���Ա��죩
	'''
	level1_3_G = []
	level1_3_S = []
	if var_brca["snv_s"]["B1_L3"] or var_brca["snv_s"]["B2_L3"]:
		level1_3_G.append("���岻��ȷ����")
		level1_3_S.append("�ٴ����岻��")
	if var_brca["snv_s"]["B1_L2"] or var_brca["snv_s"]["B2_L2"]:
		level1_3_G.append("�������Ա���")
		level1_3_S.append("�������Ա���")
	if var_brca["snv_s"]["B1_L1"] or var_brca["snv_s"]["B2_L1"]:
		level1_3_G.append("���Ա���")
		level1_3_S.append("���Ա���")
	level13_str_G = "��".join(level1_3_G)	
	level13_str_S = "��".join(level1_3_S)

	return level13_str_G, level13_str_S

def var_summary_SHRJ(var_brca):
	'''
	�Ϻ��ʼã�������֯������ȫѪ��
	�������Ҫƴ��
	'''
	varInter_list = []
	for var in var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"]+var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]:
		if var["hgvs_p"] != "p.?":
			inter = "��⵽"+var["gene_symbol"]+"�����"+var["hgvs_c"]+"��"+var["hgvs_p"]+"��ͻ�䣬"+var["variant_desc_cn"]
		else:
			inter = "��⵽"+var["gene_symbol"]+"�����"+var["hgvs_c"]+"ͻ�䣬"+var["variant_desc_cn"]
		if re.search("Splicing|Intronic", var["type"]):
			inter = inter[0:-1]+"����������쳣���ӡ�"
		if re.search("FrameShift|Nonsense", var["type"]):
			inter = inter[0:-1]+"�������γɹ������˻�ʧ��ĵ��ס�"
		varInter_list.append(inter)

	for var in var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"]:
		varInter_list.append("��⵽"+var["gene_symbol"]+"����Ĵ�Ƭ��ȱʧ���������Ϊ�����²��Ա��죬����Ӱ��"+var["gene_symbol"]+"�ĵ��׹��ܡ�")
	
	return "".join(varInter_list)

def varInfo_XAJDY(gene_region, type_cn):
	'''
	��������һ
	���������ӣ���XX��������/�ں���XXͻ��
	'''
	region_dict = {
		"exon" : "������",
		"intron" : "�ں���",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "�Ǳ�����",
		"5'FLANKING" : "�Ǳ�����"
		}

	region_list_en = re.split("_", gene_region)
	region_list_cn = []
	for i in region_list_en:
		if re.search("exon", i):
			region_list_cn.append("��"+i.replace("exon", "")+"��������")
		elif re.search("intron", i):
			region_list_cn.append("��"+i.replace("intron", "")+"���ں���")
		else:
			region_cn = region_dict[i] if i in region_dict.keys() else i
			region_list_cn.append(region_cn)
	# �ȼ���Ϊ�յ����
	type_cn = type_cn if type_cn else ""
	return "��".join(region_list_cn)+type_cn

# ����������ģ����ʵ��
def var_SD(var_brca):
	'''
	ɽ��ʡ��/ɽ����³��������֯������ȫѪ����ԣ�
	ȫѪ��
		����ʵ���⵽һЩ���������Ա��졢���Ա��죩
	��֯��
		����ʵ���⵽һЩ������/�������Ա��죩
	��ԣ�
		����ʵ���⵽һЩ���������Ա��졢���Ա��졢����/�������Ա��죩
	'''

	level1_2_G = []
	level1_2_S = []
	g_stran = {2 : "�������Ա���", 1 : "���Ա���"}
	g_rule = ["�������Ա���", "���Ա���"]
	s_stran = {2 : "����/�������Ա���", 1 : "����/�������Ա���"}
	s_rule = ["�������Ա���", "����/�������Ա���", "���Ա���"]

	# ��ϵ �����ڵ���ȫѪ������е�ȫѪ����
	for var in var_brca["snv_m"]["B1_G_L2"]+var_brca["snv_m"]["B2_G_L2"]+var_brca["snv_m"]["B1_G_L1"]+var_brca["snv_m"]["B2_G_L1"]:
		level1_2_G.append(g_stran.get(var["clinic_num_g"]))
	level1_2_G_str = "��".join(sorted(list(set(level1_2_G)), key=g_rule.index))

	# ��ϸ�� �����ڵ�����֯������е���֯����
	for var in var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"]:
		if var["var_origin"] == "germline":
			level1_2_S.append(g_stran.get(var["clinic_num_g"]))
		else:
			level1_2_S.append(s_stran.get(var["clinic_num_s"]))
	level1_2_S_str = "��".join(sorted(list(set(level1_2_S)), key=s_rule.index))

	return level1_2_G_str, level1_2_S_str

def var_summary_WHXH(var_brca):
	'''
	�人Э�ͣ�������֯������ȫѪ����ԣ���������ϵ�ĵȼ�
	1. �ܽ᣺
		���ԣ����ڱ����ͼ�������δ����к�ͻ��������к�ͻ�䡣
		���ԣ�����ʵ���⵽X��XXͻ�䡣
	2. level_3�ܽ᣺�������ݿ���Ϣ��Ϊ < ���岻��ȷͻ�䡢�����޺�ͻ�䡢�޺�ͻ�䡢��̬�Ըı� >
	'''
	# 4/5���ܽ�
	level4_5 = []
	level_4_count = len(var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]+var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"])
	level_5_count = len(var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"])
	if level_5_count:
		level4_5.append(str(level_5_count)+"���к�ͻ��")
	if level_4_count:
		level4_5.append(str(level_4_count)+"�������к�ͻ��")
	level4_5_summary = "��".join(level4_5)

	#1/2/3���ܽ�-��һ������ϵ�����1/2/3�ܽᣬ����ɽ��ģ��
	level_stran = {3 : "���岻��ȷͻ��", 2 : "�����޺�ͻ��", 1 : "�޺�ͻ��"}
	sort_rule = ["���岻��ȷͻ��", "�����޺�ͻ��", "�޺�ͻ��", "��̬�Ըı�"]
	def level_sum(var1_3_list):
		level1_3 = []
		for var in var_brca["snv_s"]["B1_L3"]+var_brca["snv_s"]["B2_L3"]+var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"]:
			if var["clinic_num_g"] in [2,3]:
				level1_3.append(level_stran.get(var["clinic_num_g"]))
			else:
				if "tag" in var.keys() and re.search("Polymorphism", var["tag"]):
					level1_3.append("��̬�Ըı�")
				else:
					level1_3.append(level_stran.get(var["clinic_num_g"]))
		return "��".join(sorted(list(set(level1_3)), key=sort_rule.index))
	all_origin_var = var_brca["snv_s"]["B1_L3"]+var_brca["snv_s"]["B2_L3"]+var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"]
	germline_var = var_brca["snv_m"]["B1_G_L3"]+var_brca["snv_m"]["B2_G_L3"]+var_brca["snv_m"]["B1_G_L2"]+var_brca["snv_m"]["B2_G_L2"]+var_brca["snv_m"]["B1_G_L1"]+var_brca["snv_m"]["B2_G_L1"]
	level1_3_summary = level_sum(all_origin_var)
	level1_3_summary_G = level_sum(germline_var)

	return level4_5_summary, level1_3_summary, level1_3_summary_G

# ����������ģ����ʵ��	
def var_level123_SDZL(var_brca):
	'''
	ɽ��������������֯������ȫѪ����ԣ�
	level1_3�ܽ�
	����ȫѪ����Ϊ < ���岻��ȷ���졢�������Ա��졢���Ա��� >
	������֯�����岻��ȷ���졢����/�������Ա���
	��ԣ����岻��ȷ���졢�������Ա��졢���Ա��졢����/�������Ա���

	'''
	g_stran = {3 : "���岻��ȷ����", 2 : "�������Ա���", 1 : "���Ա���"}
	s_stran = {3 : "���岻��ȷ����", 2 : "����/�������Ա���", 1 : "����/�������Ա���"}
	sort_rule = ["���岻��ȷ����", "�������Ա���", "����/�������Ա���", "���Ա���"]
	# ������BRCAȫѪ
	level1_3_G = [g_stran.get(var["clinic_num_g"]) for var in var_brca["snv_m"]["B1_G_L3"]+var_brca["snv_m"]["B2_G_L3"]+var_brca["snv_m"]["B1_G_L2"]+var_brca["snv_m"]["B2_G_L2"]+var_brca["snv_m"]["B1_G_L1"]+var_brca["snv_m"]["B2_G_L1"]]
	level1_3_G_str = "��".join(sorted(list(set(level1_3_G)), key=sort_rule.index))
	# ������BRCA��֯�����
	level1_3_S = [s_stran.get(var["clinic_num_s"]) for var in var_brca["snv_s"]["B1_L3"]+var_brca["snv_s"]["B2_L3"]+var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"] if var["var_origin"] != "germline"]
	level1_3_S += [g_stran.get(var["clinic_num_g"]) for var in var_brca["snv_s"]["B1_L3"]+var_brca["snv_s"]["B2_L3"]+var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"] if var["var_origin"] == "germline"]
	level1_3_S_str = "��".join(sorted(list(set(level1_3_S)), key=sort_rule.index))

	return level1_3_G_str, level1_3_S_str

def brca_judge_inApprovalTumor(var_brca):
	'''
	�ж�BRCA��Ŀ�����Ƿ��ڻ����������Ҫ������ϸ�����죬�ж���ߵȼ���I�໹��II�ࣩ
	�жϹ���
	�����Ʒ�������ߵȼ�ΪA��B���ж�Ϊ�ڻ��������C��D�ж�Ϊ��������
	'''
	var = var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L4"]+var_brca["snv_s"]["B2_L5"]
	level_list = []
	for i in var:
		# ��������Ʒ���Ϊ�գ��ᱻ�ж���Ǳ���ٴ����壡��
		if "evi_sum" in i.keys() and i["evi_sum"] and "evi_split" in i["evi_sum"].keys() and i["evi_sum"]["evi_split"] and "Predictive" in i["evi_sum"]["evi_split"].keys():
			level_list += [j["evi_conclusion_simple"] for j in i["evi_sum"]["evi_split"]["Predictive"]]
	if level_list:
		if set(["A", "B"]) & set(level_list):
			return "intumor"
		else:
			return "outtumor"

def judge_GA_tumor_KNB(var_data, tumor_list):
	'''
	������ɽCP40������Ϊ�������������ʱ��������ж�KNB
	'''
	# ����KNB�ж�����KRAS/NRAS��Ϊ�°�/�����°����Ҳ���exon2/3/4, δ��⵽BRAF V600E-2022.11.21
	GA_tumor_list = ["ʳ�ܰ�","θ��","θʳ�ܽ��簩","���Ű�","С����","������","���ܰ�","���Ұ�","�ΰ�","θ����������","���ٰ�","��β�ٰ�",\
					 "�������ٰ�","ʮ��ָ���ٰ�","ʳ����״ϸ����","ʳ���ٰ�","θ�ٰ�","���Ű�"]
	level1_2_var = [var for var in var_data if var["clinic_num_s"] in [5, 4] and var["bio_category"] == "Snvindel"]
	for var in level1_2_var:
		if any([var["gene_symbol"] in ["KRAS", "NRAS"] and var["gene_region"] in ["exon2", "exon3", "exon4"],
				var["gene_symbol"] == "BRAF" and var["hgvs_p"] == "p.(V600E)"]):
			break
	else:
		if set(tumor_list) & set(GA_tumor_list):
			return "yes"
		else:
			return 0

def sv_shfk(var_data):
	'''
	�Ϻ��ο�CP40����⵽�ں�ʱ������˵����Ҫд�����μ�⵽gene1-gene2�ںϡ�������ں���Ϊgene1:exon-gene2:exon��
	2022.07.27 ������MET 14 skippingҲҪչʾ��DNA��RNA����⵽����չʾDNA������+������ʵ����RNAˮƽҲ��⵽MET exon14 skipping������ֻ��RNA��⵽ʱ���RNA�������DNA��⵽��ʱ�޷��ж������Ž����
	'''
	sv_var = [var for var in var_data if var["bio_category"] == "Sv" and not (var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET")]
	RNA_met = [var for var in var_data if var["bio_category"] == "Sv" and var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET"]
	DNA_RNA_met = [var for var in var_data if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "MET" and "judge_mergeMET" in var.keys() and var["judge_mergeMET"]]
	sv_str_list = []
	for var in sv_var:
		if "var_desc_merge" in var.keys():
			sv_str_list.append("��⵽"+var["five_prime_gene"]+"-"+var["three_prime_gene"]+"�ںϣ�������ں���Ϊ"+var["var_desc_merge"])
	rna_str_list = []
	for var in RNA_met:
		rna_str_list.append(var["variant_desc_cn"]+var["variant_interpret_cn"])
	dna_rna_str_list = []
	for var in DNA_RNA_met:
		if var["hgvs_p"] == "p.?":
			dna_rna_str_list.append("����ʵ���⵽"+var["gene_symbol"]+" "+var["hgvs_c"]+"���죬"+var["variant_desc_cn"]+var["variant_interpret_cn"]+"����ʵ����RNAˮƽҲ��⵽MET exon14 skipping��")
		else:
			dna_rna_str_list.append("����ʵ���⵽"+var["gene_symbol"]+" "+var["hgvs_p"]+"���죬"+var["variant_desc_cn"]+var["variant_interpret_cn"]+"����ʵ����RNAˮƽҲ��⵽MET exon14 skipping��")
	
	sv_str = "����"+"��".join(sv_str_list)+"��" if sv_str_list else ""
	rna_str = ";".join(rna_str_list)
	dna_rna_str = ";".join(dna_rna_str_list)
	return sv_str + rna_str + dna_rna_str

def CP40_split_gene(var_data):
	'''
	�ĳǡ�������һCP40����������Ҫչʾ��10gene+30gene��ֻ����I/II����죩
	special["var_cp40_split"] = {
		gene10_level_I_II : [],
		gene30_level_I_II : []
	}
	'''
	gene10_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	var_levelI_II = var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]
	var_levelIII = var_data["var_somatic"]["level_onco_nodrug"] + var_data["var_somatic"]["level_III"]
	# ����III��-���ø�����ɽ����ҽԺ
	gene10_level_I_II = [var for var in var_levelI_II if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_I_II = [var for var in var_levelI_II if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene10_level_III = [var for var in var_levelIII if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_III = [var for var in var_levelIII if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	# III����Ϊ����������չ��غ�III�����-�����ڽ���ԺCP40-2023.06.09
	gene10_level_nj_onco_nodrug = [var for var in var_data["var_somatic"]["level_onco_nodrug"] if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_nj_onco_nodrug = [var for var in var_data["var_somatic"]["level_onco_nodrug"] if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene10_level_nj_III = [var for var in var_data["var_somatic"]["level_III"] if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_nj_III = [var for var in var_data["var_somatic"]["level_III"] if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	return gene10_level_I_II, gene30_level_I_II, gene10_level_III, gene30_level_III, gene10_level_nj_onco_nodrug, gene30_level_nj_onco_nodrug, gene10_level_nj_III, gene30_level_nj_III

def CP40_FJFY_summary(level_I, level_II, level_III, level_onco_nodrug):
	'''
	����ҽ�Ƹ�һCP40����������ҪչʾI/II���������
	I�ࣺXX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��ͻ����ΪX����XX��XX��������/������ҩ��A����
	II�ࣺXX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��ͻ����ΪX

	snvindel��XX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��ͻ����ΪX
	cnv��XX�����⵽������������ΪX
	sv��XX�����⵽XX-XX�ںϣ�������ΪX copies
	'''
	def var_info_stran(var):
		var_info = ""
		region_dict = {
		"exon" : "������",
		"intron" : "�ں���",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "�Ǳ�����",
		"5'FLANKING" : "�Ǳ�����"
		}
		if var["bio_category"] == "Cnv":
			var_info = var["gene_symbol"]+"�����⵽������������Ϊ"+str(var["cn_mean"])
		elif var["bio_category"] == "Sv":
			var_info = var["gene_symbol"]+"�����⵽"+var["five_prime_gene"]+"-"+var["three_prime_gene"]+"�ںϣ�������Ϊ"+str(var["copies"])+" copies"
		elif var["bio_category"] == "Snvindel":
			# ��һ��MET exon 14��Ծͻ��-2023.07.13
			if var["gene_symbol"] == "MET" and "exon14 skipping" in var["variant_interpret_cn"]:
				if var["hgvs_p_FJFY"] != "p.?":
					var_info  = "MET����14����������Ծͻ��"+var["hgvs_c"]+"("+var["hgvs_p_FJFY"]+")��ͻ����Ϊ"+str(var["freq_str"])
				else:
					var_info = "MET����14����������Ծͻ��"+var["hgvs_c"]+"��ͻ����Ϊ"+str(var["freq_str"])
			# MET exon14��Ծ��ӽ���-2023.07.13
			else:
				region_list_en = re.split("_", var["gene_region"])
				region_list_cn = []
				for i in region_list_en:
					if re.search("exon", i):
						region_list_cn.append(i.replace("exon", "")+"��������")
					elif re.search("intron", i):
						region_list_cn.append(i.replace("intron", "")+"���ں���")
					else:
						region_cn = region_dict[i] if i in region_dict.keys() else i
						region_list_cn.append(region_cn)
				if var["hgvs_p_FJFY"] != "p.?":
					var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽" + var["type_cn"] + var["hgvs_p_FJFY"]+"��ͻ����Ϊ"+str(var["freq_str"])
				else:
					if var["type_cn"] != "--":
						var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽"+ var["type_cn"] + var["hgvs_c"]+"ͻ�䣬ͻ����Ϊ"+str(var["freq_str"])
					else:
						var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽"+ var["hgvs_c"]+"ͻ�䣬ͻ����Ϊ"+str(var["freq_str"])
		return var_info

	def regimen_info(var):
		regimen_info = ""
		sense_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if \
							re.search("Sensitive", i["clinical_significance"]) and \
							i["evi_conclusion_simple"]=="A"] if \
								"Predictive" in var["evi_sum"]["evi_split"].keys() else \
								[]
		resis_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if \
							re.search("Resistant", i["clinical_significance"]) and \
							i["evi_conclusion_simple"]=="A"] if \
								"Predictive" in var["evi_sum"]["evi_split"].keys() else \
								[]

		sense_regimen_str = "��" + "��".join(sense_regimen) + "�������У�A����" if sense_regimen else ""
		resis_regimen_str = "��" + "��".join(resis_regimen) + "������ҩ��A����" if resis_regimen else ""

		return sense_regimen_str, resis_regimen_str
	# I�����
	var_levelI_sum = []

	# ������ʼ-2022.11.08
	lc10_gene = ["BRAF", "EGFR", "ERBB2", "KRAS", "MET", "ALK", "ROS1", "RET", "NRAS", "PIK3CA"]
	var_levelI_sum_lc10_list = []
	var_levelI_sum_other_list = []
	# ��������-2022.11.08

	for var in level_I:
		tmp_list = []
		var_info = var_info_stran(var)
		sense_regimen_str, resis_regimen_str = regimen_info(var)
		tmp_list.append(var_info)
		if sense_regimen_str:
			tmp_list.append(sense_regimen_str)
		if resis_regimen_str:
			tmp_list.append(resis_regimen_str)
		var_levelI_sum.append("��".join(tmp_list))

		# ������ʼ-2022.11.08
		for gene in re.split(",",var["gene_symbol"]):
			if gene in lc10_gene:
				if "��".join(tmp_list) and "��".join(tmp_list) not in var_levelI_sum_lc10_list:
					var_levelI_sum_lc10_list.append("��".join(tmp_list))
			else:
				if "��".join(tmp_list) and "��".join(tmp_list) not in var_levelI_sum_other_list:
					var_levelI_sum_other_list.append("��".join(tmp_list))
		# ��������-2022.11.08

	level_I_sum = "��".join(var_levelI_sum) if var_levelI_sum else ""
	var_levelI_sum_lc10 = "��".join(var_levelI_sum_lc10_list) if var_levelI_sum_lc10_list else ""
	var_levelI_sum_other = "��".join(var_levelI_sum_other_list) if var_levelI_sum_other_list else ""

	# II�����
	level_II_sum = "��".join([var_info_stran(var) for var in level_II]) if level_II else ""
	# ������ʼ-2022.11.08
	level_II_sum_lc10_list = []
	level_II_sum_other_list = []
	for var in level_II:
		for gene in re.split(",", var["gene_symbol"]):
			if gene in lc10_gene:
				if var_info_stran(var) and var_info_stran(var) not in level_II_sum_lc10_list:
					level_II_sum_lc10_list.append(var_info_stran(var))
			else:
				if var_info_stran(var) and var_info_stran(var) not in level_II_sum_other_list:
					level_II_sum_other_list.append(var_info_stran(var))

	level_II_sum_lc10 = "��".join(level_II_sum_lc10_list) if level_II_sum_lc10_list else ""
	level_II_sum_other = "��".join(level_II_sum_other_list) if level_II_sum_other_list else level_II_sum_other_list
	# ��������������Ҫ�������������������չ���
	level_III_sum_other = []
	for var in level_III+level_onco_nodrug:
		for gene in re.split(",", var["gene_symbol"]):
			if gene not in lc10_gene:
				level_III_sum_other.append(var)
	# ��������-2022.11.08

	result = {
		"level_I_sum" : level_I_sum,
		"level_II_sum" : level_II_sum,
		"level_I_lc10" : var_levelI_sum_lc10,
		"level_II_lc10" : level_II_sum_lc10,
		"level_I_other" : var_levelI_sum_other,
		"level_II_other" : level_II_sum_other,
		"level_III_other" : level_III_sum_other
	}
	
	return result
		
def Master_sum(var_data):
	# ͳ����ϸ�����������DNA/RNA�ֿ�ͳ�ƣ�ȡ��������
	result = s_var_rule(var_data)

	#level_I = [var for var in var_data if var["clinic_num_s"] == 5 and var["var_origin"] != "germline" and "evi_sum" in var.keys() and var["evi_sum"]["evi_split"]]
	#level_II = [var for var in var_data if var["clinic_num_s"] == 4 and var["var_origin"] != "germline" and "evi_sum" in var.keys() and var["evi_sum"]["evi_split"]]
	#level_onco_nodrug = [var for var in var_data if var["clinic_num_s"] in [5, 4] and var["var_origin"] != "germline" and "evi_sum" in var.keys() and not var["evi_sum"]["evi_split"]]
	#level_III = [var for var in var_data if var["clinic_num_s"] == 3 and var["var_origin"] != "germline"]

	def var_count(var_list):
		num = 0
		for var in var_list:
			num += 1
			# DNA/RNA�ں�������һ��������Ϊ1-2023.02.08
	#		if "rna_detect" in var.keys():
	#			num += 1 
		return num
	return var_count(result["level_I"]), var_count(result["level_II"]), var_count(result["level_onco_nodrug"]), var_count(result["level_III"])

def varInfo_SYX(gene_region, type_cn, type, hgvs_c, hgvs_p):
	'''
	������
	���������ӣ�XX��������/�ں���XX hgvs_p/hgvs_c XXͻ��
	'''
	region_dict = {
		"exon" : "������",
		"intron" : "�ں���",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "�Ǳ�����",
		"5'FLANKING" : "�Ǳ�����"
		}

	type_stran = {
		"Intronic" : "�ں���ͻ��",
		"3'UTR" : "3'UTRͻ��",
		"5'UTR" : "5'UTRͻ��",
		"3'FLANKING" : "�Ǳ�����ͻ��",
		"5'FLANKING" : "�Ǳ�����ͻ��"
	}
	region_list_en = re.split("_", gene_region)
	region_list_cn = []
	for i in region_list_en:
		if re.search("exon", i):
			region_list_cn.append(i.replace("exon", "")+"��������")
		elif re.search("intron", i):
			region_list_cn.append(i.replace("intron", "")+"���ں���")
		else:
			region_cn = region_dict[i] if i in region_dict.keys() else i
			region_list_cn.append(region_cn)
	# �ȼ���Ϊ�յ����
	type_cn = type_cn if type_cn and type_cn != "--" else type_stran.get(type, "ͻ��") 
	var_info = hgvs_p if hgvs_p != "p.?" else hgvs_c
	return "��".join(region_list_cn) + var_info + type_cn

def varInfo_FDZS(gene_symbol, gene_region, variant_desc_cn, config):
	'''
	������ɽ
	1����������ʹ�����ñ�
	2���Ŵ�����ʹ�����ñ�
	3���Ŵ�����������Ϣʹ�����ñ�
	4������ʹ�����ñ�
	5��gene_regionת��Ϊ���ģ��硰exon1��ת��Ϊ1��������
	6��������������ȡ������Ϣ����չʾ�����»�����뵰�׵�1056λ�������ɹȰ�����ͻ��Ϊ��ֹ�����ӡ�
	'''
	# 1. ��ȡ���ñ���Ϣ
	#FDZS_gene_dict = get_FDZS_database(config)
	FDZS_gene_dict, FJZL_dict, io_dict = getconfigxlsx(config)
	result = FDZS_gene_dict[gene_symbol] if gene_symbol in FDZS_gene_dict.keys() and FDZS_gene_dict[gene_symbol] else {}
	# 2. ��ȡgene_region
	result["gene_region_cn"] = ""
	if gene_region:
		region_dict = {
			"exon" : "������",
			"intron" : "�ں���",
			"3'UTR" : "3'UTR",
			"5'UTR" : "5'UTR",
			"3'FLANKING" : "�Ǳ�����",
			"5'FLANKING" : "�Ǳ�����"
			}

		region_list_en = re.split("_", gene_region)
		region_list_cn = []
		for i in region_list_en:
			if re.search("exon", i):
				region_list_cn.append(i.replace("exon", "")+"��������")
			elif re.search("intron", i):
				region_list_cn.append(i.replace("intron", "")+"���ں���")
			else:
				region_cn = region_dict[i] if i in region_dict.keys() else i
				region_list_cn.append(region_cn)
		result["gene_region_cn"] = "��".join(region_list_cn)
	# 3. �����������
	var_desc = re.split("��",variant_desc_cn) if variant_desc_cn else []
	result["var_desc"] = var_desc[1].replace("��","") if var_desc and len(var_desc) >= 2 else ""

	return result

def var_summary_CQFY(var_brca):
	'''
	���츽һ��ԺBRCA��������ȫѪ��
	***������MLPA***
	1. ������������м���Ҳ�ͬ����ȼ����ֿ�д
	       ��������⵽BRCA1����Ϊ2��--�������ԣ�likely benign���ı��죻BRCA2����Ϊ1��--���ԣ�benign���ı��졣
	2. ������������м����Ϊͬһ����ȼ���дһ��
	       ��������⵽BRCA1��BRCA2����Ϊ2��--�������ԣ�likely benign���ı��졣
	3. ��ͬһ��������ڶ��3�����ϱ���λ�����������г�����4,5����չʾ
	       ��������⵽BRCA1������5��--�²���pathogenic����3��--���岻����uncertain significance���ı��죻BRCA2����Ϊ1��--���ԣ�benign���ı��졣
	'''	
	# �ȼ�ת��
	level_dict = {
		"5" : "5��--�²���pathogenic��",
		"4" : "4��--�����²���likely pathogenic��",
		"3" : "3��--����δ��(uncertain significance)",
		"2" : "2��--�������ԣ�likely benign��",
		"1" : "1��--���ԣ�benign��"
		}

	result = {}

	# ��ȡBRCA1��������еȼ���������MLPA��
	brca1_level_list = [str(int(i["clinic_num_g"])) for i in var_brca["snv_s"]["B1_L5"] + \
															 var_brca["snv_s"]["B1_L4"] + \
															 var_brca["snv_s"]["B1_L3"] + \
															 var_brca["snv_s"]["B1_L2"] + \
															 var_brca["snv_s"]["B1_L1"]]
	# ��ȡBRCA2��������еȼ���������MLPA��
	brca2_level_list = [str(int(i["clinic_num_g"])) for i in var_brca["snv_s"]["B2_L5"] + \
															 var_brca["snv_s"]["B2_L4"] + \
															 var_brca["snv_s"]["B2_L3"] + \
															 var_brca["snv_s"]["B2_L2"] + \
															 var_brca["snv_s"]["B2_L1"]]
	# L21����ȡ12�е���ߵȼ����ַ���
	result["B1_L21"] = level_dict["2"] if "2" in brca1_level_list else level_dict["1"]
	result["B2_L21"] = level_dict["2"] if "2" in brca2_level_list else level_dict["1"]
	# L3 ����ȡ3���ַ���
	result["B1_L3"] = level_dict["3"] if "3" in brca1_level_list else ""
	result["B2_L3"] = level_dict["3"] if "3" in brca2_level_list else ""
	# L54 ��ȡ�����45���ַ�����ȥ�غ�ƴ�ӣ�
	result["B1_L54"] = "��".join(sorted(set([level_dict[i] for i in brca1_level_list if i in ["4","5"]]), reverse=True))
	result["B2_L54"] = "��".join(sorted(set([level_dict[i] for i in brca2_level_list if i in ["4","5"]]), reverse=True))
	# L543��ȡ345���б�ȥ�أ��������ж��������á��С����ǡ�Ϊ������1���á�Ϊ��������á��С�
	result["B1_L543"] = list(set([level_dict[i] for i in brca1_level_list if i in ["3","4","5"]]))
	result["B2_L543"] = list(set([level_dict[i] for i in brca2_level_list if i in ["3","4","5"]]))

	return result

def varInfo_FJZL(var, tumor_list,config):
	'''
	�����������������ڶ��Σ���Ҫչʾ���췢��Ƶ��
	���ĶΣ���Ҫ�ϲ����Ʒ���
	'''
	# ��ȡ���췢��Ƶ����Ϣ
	#fjzl_database = get_FJZL_database(config)
	fdzs_dict, fjzl_database, io_dict = getconfigxlsx(config)
	var_freq_info = {}
	for gene in re.split(",", var["gene_symbol"]):
		for tumor in tumor_list:
			# ��ƥ�䵽���� ���֣���ʹ�ø���Ϣ
			if (gene, tumor) in fjzl_database.keys():
				if gene not in var_freq_info.keys():
					var_freq_info.setdefault(gene, "")
				var_freq_info[gene] = fjzl_database[(gene, tumor)]
				
		# ��δƥ�䵽���� ���֣���ƥ�����
		if gene not in var_freq_info.keys() and (gene, "") in fjzl_database.keys():
			var_freq_info.setdefault(gene, "")
			var_freq_info[gene] = fjzl_database[(gene, "")]
	
	# ���Ʒ���С�ᴦ��
	# ��ҩ�ﰴ����/��ҩ���ȼ����в�֣����ں�������
	regimen_info = {}
	def getinfo_S(regimen_list, conclusion_list):
		return [a["regimen_name"] for a in regimen_list if a["evi_conclusion_simple"] in conclusion_list]

	for i in ["A", "B", "C", "D"]:
		regimen_info["regimen_"+str(i)+"_S"] = getinfo_S(var["evi_sum"]["regimen_S"], [i])
		regimen_info["regimen_"+str(i)+"_R"] = getinfo_S(var["evi_sum"]["regimen_R"], [i])
	
	return var_freq_info, regimen_info

def BRCA_FJFY_summary(var_data):
	'''
	����ҽ�Ƹ�һBRCA����������ҪչʾI/II���������-2022.11.11
	*BRCAȫѪ
		5�ࣺBRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p(P.XX)/hgvs_c��������Ϊ����/�Ӻϣ����ΪXX������XX��XX��������/������ҩ��A����
		4�ࣺBRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p(P.XX)/hgvs_c��������Ϊ����/�Ӻϣ����ΪXX������XX��XX��������/������ҩ��A����
		3�ࣺ����
	*BRCA��֯
		I�ࣺBRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p(P.XX)/hgvs_c��ͻ����ΪXX����XX��XX��������/������ҩ��A����
		II�ࣺBRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p(P.XX)/hgvs_c��ͻ����ΪXX
		III�ࣺBRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p(P.XX)/hgvs_c��ͻ����ΪXX

	snvindel��XX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��ͻ����ΪX
	mlpa��XX����value��⵽��Ƭ��ȱʧ/��Ƭ���ظ�
	'''
	def var_info_stran(var):
		var_info = ""
		region_dict = {
		"exon" : "������",
		"intron" : "�ں���",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "�Ǳ�����",
		"5'FLANKING" : "�Ǳ�����"
		}
		if var["type"] == "Loss":
			var_info = var["gene_symbol"]+"����"+var["value"]+"��⵽��Ƭ��ȱʧ"
		elif var["type"] == "Gain":
			var_info = var["gene_symbol"]+"����"+var["value"]+"��⵽��Ƭ���ظ�"
		elif var["bio_category"] == "Snvindel":
			region_list_en = re.split("_", var["gene_region"])
			region_list_cn = []
			for i in region_list_en:
				if re.search("exon", i):
					region_list_cn.append(i.replace("exon", "")+"��������")
				elif re.search("intron", i):
					region_list_cn.append(i.replace("intron", "")+"���ں���")
				else:
					region_cn = region_dict[i] if i in region_dict.keys() else i
					region_list_cn.append(region_cn)
			if var["hgvs_p"] != "p.?":
				if var["var_origin"] == "germline":
					gene_type = "����" if float(var["freq"]) >= 0.85 else "�Ӻ�"
					var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽" + var["type_cn"] + var["hgvs_p"]+"��������Ϊ"+gene_type+"�����Ϊ"+var["freq_str"]+"��"
				else:
					var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽" + var["type_cn"] + var["hgvs_p"]+"��ͻ����Ϊ"+str(var["freq_str"])
			else:
				if var["var_origin"] == "germline":
					gene_type = "����" if float(var["freq"]) >= 0.85 else "�Ӻ�"
					var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽" + var["type_cn"] + var["hgvs_c"]+"��������Ϊ"+gene_type+"�����Ϊ"+var["freq_str"]+"��"
				else:
					var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽"+ var["type_cn"] + var["hgvs_c"]+"��ͻ����Ϊ"+str(var["freq_str"])
		return var_info

	def regimen_info(var):
		sense_regimen = []
		resis_regimen = []
		if var["evi_sum"]:
			regimen_info = ""
			#print (var["evi_sum"])
			sense_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] \
							if re.search("Sensitive", i["clinical_significance"]) and i["evi_conclusion_simple"]=="A"] \
							if var["evi_sum"]["evi_split"] and "Predictive" in var["evi_sum"]["evi_split"].keys() else []
			resis_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] \
							if re.search("Resistant", i["clinical_significance"]) and i["evi_conclusion_simple"]=="A"] \
							if var["evi_sum"]["evi_split"] and "Predictive" in var["evi_sum"]["evi_split"].keys() else []

		sense_regimen_str = "��" + "��".join(sense_regimen) + "�������У�A����" if sense_regimen else ""
		resis_regimen_str = "��" + "��".join(resis_regimen) + "������ҩ��A����" if resis_regimen else ""

		return sense_regimen_str, resis_regimen_str

	# ����λ���ܽ�
	def var_inter(var):
		var_info = var_info_stran(var)
		sense_regimen_str, resis_regimen_str = regimen_info(var)
		tmp_list = []
		tmp_list.append(var_info)
		if sense_regimen_str:
			tmp_list.append(sense_regimen_str)
		if resis_regimen_str:
			tmp_list.append(resis_regimen_str)
		return "��".join(tmp_list)

	result = {}
	# BRCAȫѪ
	result["gBRCA_5"] = "��".join([var_inter(var) for var in var_data["snv_s"]["B1_L5"] + var_data["snv_s"]["B2_L5"]])
	result["gBRCA_4"] = "��".join([var_inter(var) for var in var_data["snv_s"]["B1_L4"] + var_data["snv_s"]["B2_L4"] + var_data["mlpa"]["B1_Loss"]+var_data["mlpa"]["B2_Loss"]])
	result["gBRCA_3"] = "��".join([var_inter(var) for var in var_data["snv_s"]["B1_L3"] + var_data["snv_s"]["B2_L3"] + var_data["mlpa"]["B1_Gain"]+var_data["mlpa"]["B2_Gain"]])
	# BRCA��֯
	result["tBRCA_I"] = "��".join([var_inter(var) for var in var_data["snv_s"]["B1_L5"] + var_data["snv_s"]["B2_L5"]+ var_data["snv_s"]["B1_L4"] + var_data["snv_s"]["B2_L4"] if var["clinic_num_s"] == 5])
	result["tBRCA_II"] = "��".join([var_inter(var) for var in var_data["snv_s"]["B1_L5"] + var_data["snv_s"]["B2_L5"]+ var_data["snv_s"]["B1_L4"] + var_data["snv_s"]["B2_L4"] if var["clinic_num_s"] == 4])
	result["tBRCA_III"] = "��".join([var_inter(var) for var in var_data["snv_s"]["B1_L3"] + var_data["snv_s"]["B2_L3"]])
		
	return result

def HRR_FJFY_summary(var_data):
	'''
	����ҽ�Ƹ�һHRR����������ҪչʾI/II���������-2022.11.17
	*HRRȫѪ
		BRCA1/BRCA2 5�ࣺXXX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��������Ϊ����/�Ӻϣ���XX��XX��������/������ҩ��A����
		BRCA1/BRCA2 4�ࣺXXX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��������Ϊ����/�Ӻϣ���XX��XX��������/������ҩ��A����
		BRCA1/BRCA2 3�ࣺXXX����X��������/�ں��Ӽ���XXͻ��hgvs_p/hgvs_c��������Ϊ����/�Ӻ�
		�������� 5/4/3�ࣺ��չʾ�Ƿ��м�⵽��������ϸ��Ϣ
	*HRR��֯
		BRCA1/BRCA2 I�ࣺBRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��ͻ����ΪXX����XX��XX��������/������ҩ��A����
		BRCA1/BRCA2 II�ࣺBRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p(P.XX)/hgvs_c��ͻ����ΪXX
		BRCA1/BRCA2 III��+����������չ��أ�BRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p(P.XX)/hgvs_c��ͻ����ΪXX
		�������� I/II/III��+����������չ��أ���չʾ�Ƿ��м�⵽��������ϸ��Ϣ
	*HRR���
	*��ֻ�迪�����������ȫѪ����֯�Ŀ�����ԵĽ����
		BRCA1/BRCA2 5�ࣺXXX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��������Ϊ����/�Ӻϣ���XX��XX��������/������ҩ��A����
		BRCA1/BRCA2 4�ࣺXXX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��������Ϊ����/�Ӻϣ���XX��XX��������/������ҩ��A����
		BRCA1/BRCA2 I�ࣺBRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��ͻ����ΪXX����XX��XX��������/������ҩ��A����
		BRCA1/BRCA2 II�ࣺBRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p(P.XX)/hgvs_c��ͻ����ΪXX
		BRCA1/BRCA2 3�ࣺXXX����X��������/�ں��Ӽ���XXͻ��hgvs_p/hgvs_c��������Ϊ����/�Ӻ�
		BRCA1/BRCA2 III��+����������չ��أ�BRCAX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p(P.XX)/hgvs_c��ͻ����ΪXX
		�������� 5/4/3���I/II/III+����������չ��أ���չʾ�Ƿ��м�⵽��������ϸ��Ϣ

	snvindel��XX����X��������/�ں��Ӽ�⵽XXͻ��hgvs_p/hgvs_c��ͻ����ΪX
	mlpa��XX����value��⵽��Ƭ��ȱʧ/��Ƭ���ظ�
	'''
	def var_info_stran(var):
		var_info = ""
		region_dict = {
		"exon" : "������",
		"intron" : "�ں���",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "�Ǳ�����",
		"5'FLANKING" : "�Ǳ�����"
		}
		if var["type"] == "Loss":
			var_info = var["gene_symbol"]+"����"+var["value"]+"��⵽��Ƭ��ȱʧ"
		elif var["type"] == "Gain":
			var_info = var["gene_symbol"]+"����"+var["value"]+"��⵽��Ƭ���ظ�"
		elif var["bio_category"] == "Snvindel":
			region_list_en = re.split("_", var["gene_region"])
			region_list_cn = []
			for i in region_list_en:
				if re.search("exon", i):
					region_list_cn.append(i.replace("exon", "")+"��������")
				elif re.search("intron", i):
					region_list_cn.append(i.replace("intron", "")+"���ں���")
				else:
					region_cn = region_dict[i] if i in region_dict.keys() else i
					region_list_cn.append(region_cn)
			if var["hgvs_p"] != "p.?":
				if var["var_origin"] == "germline":
					gene_type = "����" if float(var["freq"]) >= 0.85 else "�Ӻ�"
					var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽" + var["type_cn"] + var["hgvs_p"]+"��������Ϊ"+gene_type
				else:
					var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽" + var["type_cn"] + var["hgvs_p"]+"��ͻ����Ϊ"+str(var["freq_str"])
			else:
				if var["var_origin"] == "germline":
					gene_type = "����" if float(var["freq"]) >= 0.85 else "�Ӻ�"
					var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽" + var["type_cn"] + var["hgvs_c"]+"��������Ϊ"+gene_type
				else:
					var_info = var["gene_symbol"]+"����"+"��".join(region_list_cn)+"��⵽"+ var["type_cn"] + var["hgvs_c"]+"��ͻ����Ϊ"+str(var["freq_str"])
		return var_info

	def regimen_info(var):
		sense_regimen = []
		resis_regimen = []
		if var["evi_sum"]:
			regimen_info = ""
			#print (var["evi_sum"])
			sense_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if re.search("Sensitive", i["clinical_significance"]) and i["evi_conclusion_simple"]=="A"] if var["evi_sum"]["evi_split"] and "Predictive" in var["evi_sum"]["evi_split"].keys() else []
			resis_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if re.search("Resistant", i["clinical_significance"]) and i["evi_conclusion_simple"]=="A"] if var["evi_sum"]["evi_split"] and "Predictive" in var["evi_sum"]["evi_split"].keys() else []

		sense_regimen_str = "��" + "��".join(sense_regimen) + "�������У�A����" if sense_regimen else ""
		resis_regimen_str = "��" + "��".join(resis_regimen) + "������ҩ��A����" if resis_regimen else ""

		return sense_regimen_str, resis_regimen_str

	# ����λ���ܽ�
	def var_inter(var):
		var_info = var_info_stran(var)
		sense_regimen_str, resis_regimen_str = regimen_info(var)
		tmp_list = []
		tmp_list.append(var_info)
		if sense_regimen_str:
			tmp_list.append(sense_regimen_str)
		if resis_regimen_str:
			tmp_list.append(resis_regimen_str)
		return "��".join(tmp_list)

	result = {}
	# ��ϵ���
	result["G_B_5"] = "��".join([var_inter(var) for var in var_data["var_germline"]["level_5"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel"])

	result["G_B_4"] = "��".join([var_inter(var) for var in var_data["var_germline"]["level_4"] + var_data["mlpa"]["B1_Loss"]+var_data["mlpa"]["B2_Loss"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and (var["bio_category"] == "Snvindel" or ("type" in var.keys() and var["type"] == "Loss"))])

	result["G_B_3"] = "��".join([var_inter(var) for var in var_data["var_germline"]["level_3"] + var_data["mlpa"]["B1_Gain"]+var_data["mlpa"]["B2_Gain"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and (var["bio_category"] == "Snvindel" or ("type" in var.keys() and var["type"] == "Gain"))])

	result["G_other"] = "��".join([var_inter(var) for var in var_data["var_germline"]["level_5"]+var_data["var_germline"]["level_4"]+var_data["var_germline"]["level_3"] \
															if var["gene_symbol"] not in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel"])
	# ��ϸ�����
	result["S_B_I"] = "��".join([var_inter(var) for var in var_data["var_somatic"]["level_I"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel" ])

	result["S_B_II"] = "��".join([var_inter(var) for var in var_data["var_somatic"]["level_II"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel" ])

	result["S_B_III"] = "��".join([var_inter(var) for var in var_data["var_somatic"]["level_III"] + var_data["var_somatic"]["level_onco_nodrug"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel" ])

	result["S_other"] = "��".join([var_inter(var) for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]+\
															 var_data["var_somatic"]["level_III"] + var_data["var_somatic"]["level_onco_nodrug"] \
															if var["gene_symbol"] not in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel" ])														

		
	return result

def varS3_splitfor301(var_data):
	'''
	��3�������Ϊһ����������չʾ
	'''
	# ��ΪDNA/RNA�����ںϣ���RNA���ж��ʱ����Ҫ�ֿ�չʾ-2023.04.13
	var_list = []
	for var in var_data:
		var_list.append(var)
		if "dup_rna_detect" in var.keys() and len(var["dup_rna_detect"]) > 1 :
			for rna_var in var["dup_rna_detect"]:
				var_list.append(rna_var)
		#else:
		#	var_list.append(var)

	result = []
	if var_list:
		for i in range(0, len(var_list)-2, 2):
			tmp_dict = {}
			for j in range(1, 3):
				tmp_dict["var"+str(j)] = var_list[i+j-1]
			result.append(tmp_dict)

		rest_num= len(var_list) % 2
		rest_tmp_dict = {}
		for j in range(1, 3):
			rest_tmp_dict["var"+str(j)] = ""

		num = 1
		last_row_num = len(var_list)-rest_num if rest_num != 0 else len(var_list)-rest_num-2
		for j in range(last_row_num, len(var_list)):
			rest_tmp_dict["var"+str(num)] = var_list[j]
			num += 1
		result.append(rest_tmp_dict)
	return result

def sort_var_for301(var_data):
	'''
	301����MP�������С��ͽ�����ֱ�������
	���� �����򡢱���ȼ���I��II�������͡�Ƶ�ʣ�Ƶ���ڸ����ִ����ʱ���Ѿ������ˣ�
	10������ǰ��
	'''
	LC10_gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	for var in var_data:
		gene = re.split(",", var["gene_symbol"])
		var["sort_301"] = 0 if set(gene) & set(LC10_gene_list) else 1
	return sorted(var_data, key=lambda i : (i["sort_301"], i["gene_symbol"]))

def sort_var_for301_2(var_data):
	'''
	301����MP������������
	1. �����򣨰�I/II/III������ţ�����ȼ�Խ�ߵĻ�������ǰ�棬�ȼ���ͬ�İ�����ĸ�ţ�
	2. ��ͻ�䡢�ںϡ�CNV
	'''
	# 1. �������򣬸���I/II/III�ࣻ��ͻ�䡢�ںϡ�CNV
	var_data_copy = copy.deepcopy(var_data)
	type_rule = {"Snvindel" : 0, "Sv" : 1, "PSeqRnaSv" : 2, "Cnv" : 3}
	for var in var_data_copy:
		# s_level, I���Ӧ5��II���Ӧ4������������չ��ض�Ӧ3�����ڱ������
		var["s_level"] = S_level(var)
		# ֤����ߵȼ�ת��Ϊ���֣���������I����Ϊ0��II����Ϊ1�� ����������չ��ر��Ϊ2
		var["top_level_forsort"] = 0 if var["top_level"] in ["A", "B"] else 1 if var["top_level"] in ["C", "D"] else 2
		# sv��gene_symbol�������������ڼ�ⷶΧ��򷵻�five_prime_gene,three_prime_gene (Ŀǰ���ص���three_prime_gene,five_prime_gene)
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and "," in var["gene_symbol"]:
			var["gene_symbol"] = var["five_prime_gene"]+","+var["three_prime_gene"]
	var_data_copy_sort = sorted(var_data_copy, key=lambda i : (i["top_level_forsort"], type_rule.get(i["bio_category"])))

	# 2. ��������в�֣����ں�����װ
	var_split = {}
	for var in var_data_copy_sort:
		if var["gene_symbol"] not in var_split.keys():
			var_split.setdefault(var["gene_symbol"], [])
		var_split[var["gene_symbol"]].append(var)

	# 3. ��ȡ�����Ӧ�������ߵȼ�����˳���ȡ�����б�
	gene_level_sort = sorted([{"gene_symbol" : var["gene_symbol"], "top_level_forsort" : var["top_level_forsort"]} for var in var_data_copy],\
							 key=lambda i:(i["top_level_forsort"], i["gene_symbol"]))
	gene_list = []
	for i in gene_level_sort:
		if i["gene_symbol"] not in gene_list:
			gene_list.append(i["gene_symbol"])
	
	# 4. ����������Ļ����б������������װ
	result = []
	for gene in gene_list:
		result.extend(var_split[gene])

	return result

# ������ǿ��ڱ���ģ����ʵ��
def get_summary_GDRM_gHRR(var_data, somatic_var_data):
	'''
	���ڹ㶫����gHRR������ܽ�
	���ڱ����ͼ���������⵽XXX����gene_region hgvs_c hgvs_p�²�/�����²��Ա��졣
	2023.04.26������ϸ��/��Դ�������죺���ڱ����ͼ���������⵽XXX����gene_region hgvs_c hgvs_p I��/II�����
	'''
	clinic_num_g_stran = {5 : "�²��Ա���", 4 : "�����²��Ա���"}
	var_list = []
	for var in var_data:
		if var["type"] == "Loss":
			var_list.append("{0}����{1} del {2}".format(var["gene_symbol"], var["value"], "�����²��Ա���")) 
		else:
			if var["hgvs_p"] != "p.?":
				var_list.append("{0}����{1} {2} {3}{4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"],clinic_num_g_stran.get(var["clinic_num_g"])))
			else:
				var_list.append("{0}����{1} {2}{3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"],clinic_num_g_stran.get(var["clinic_num_g"])))
		
	# 2023.04.26������ϸ��/��Դ�����������ܽ�
	clinic_num_s_stran = {5 : "I�����", 4 : "II�����"}
	somatic_var_list = []
	for var in somatic_var_data:
		if var["hgvs_p"] != "p.?":
			somatic_var_list.append("{0}����{1} {2} {3} {4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"],clinic_num_s_stran.get(var["clinic_num_s"])))
		else:
			somatic_var_list.append("{0}����{1} {2} {3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"],clinic_num_s_stran.get(var["clinic_num_s"])))

	return ";".join(var_list), ";".join(somatic_var_list), ";".join(somatic_var_list+var_list)

def PAN116_split_gene(var_data):
	'''
	������һ116����������Ҫչʾ��10gene+��������ֻ����I/II����죩�Լ���ϵ4/5�����
	special["var_pan116_split"] = {
		gene10_level_I_II : [],
		other_level_I_II : [],
		gene10_germline_45 : [],
		other_germline_45 : []
	}
	'''
	gene10_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
#	gene10_level_I_II = [var for var in var_data if set(re.split(",",var["gene_symbol"])) & set(gene10_list) and var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and var["clinic_num_s"]==5] + [var for var in var_data if set(re.split(",",var["gene_symbol"])) & set(gene10_list) and var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and var["clinic_num_s"]==4]
#	gene30_level_I_II = [var for var in var_data if not set(re.split(",",var["gene_symbol"])) & set(gene10_list) and var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and var["clinic_num_s"]==5] + [var for var in var_data if not set(re.split(",",var["gene_symbol"])) & set(gene10_list) and var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and var["clinic_num_s"]==4]
	# ����III��-���ø�����ɽ����ҽԺ
	gene10_level_I_II = [var for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"] if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_I_II = [var for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"] if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene10_level_III = [var for var in var_data["var_somatic"]["level_onco_nodrug"] + var_data["var_somatic"]["level_III"] if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_III = [var for var in var_data["var_somatic"]["level_onco_nodrug"] + var_data["var_somatic"]["level_III"] if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]

	gene10_germline_45 = [var for var in var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"] if set(re.split(",", var["gene_symbol"])) & set(gene10_list)]
	other_germline_45 = [var for var in var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"] if not set(re.split(",", var["gene_symbol"])) & set(gene10_list)]


	return gene10_level_I_II, gene30_level_I_II, gene10_level_III, gene30_level_III, gene10_germline_45, other_germline_45

# ͳ����ϸ��/��Դδ�����Ƕ�̬�Ա�����-����֣��һHRR-2023.07.03
def varnum_ZDY_HRR(var_list):
	num = 0
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			var["tag"] = var["tag"] if var["tag"] else ""
			if not re.search("Polymorphism", var["tag"]):
				num += 1
	return num

# �Ϻ��ʼ�CP40����Ҫ������������ܣ�Ҫ���̶��Ļ����б�ͻ���˳��չʾ-2023.09.18
# summary: �����ؽ�����ܣ�����I/II/III�����
# var_inter������I/II�� + ��������I/II����
# gene_inter������I/II�� + ��������I/II�࣬��Ҫȥ��
def shrj_cp40_var(var_list, config):
	var_dict = {
		"KRAS" : "��ͻ��/����/ȱʧ",
		"NRAS" : "��ͻ��/����/ȱʧ",
		"BRAF" : "��ͻ��/����/ȱʧ",
		"POLE" : "��ͻ��/����/ȱʧ",
		"ERBB2" : "��ͻ��/����/ȱʧ/����������",
		"NTRK1" : "��ͻ��/����/ȱʧ/�ں�",
		"NTRK2" : "��ͻ��/����/ȱʧ/�ں�",
		"NTRK3" : "��ͻ��/����/ȱʧ/�ں�",
	}
	gene_list = ["KRAS", "NRAS", "BRAF", "POLE", "ERBB2", "NTRK1", "NTRK2", "NTRK3"]
	gene_info = get_shrj_geneinfo(config)
	detect_gene = []

	summary_result = []
	var_inter = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			var["level"] = str(int(S_level(var)))
			#var["gene_inter_shrj"] = gene_info.get(var["gene_symbol"])
			var["gene_num"] = gene_list.index(gene) if gene in gene_list else 999
			# summary���
			if gene in var_dict.keys():
				detect_gene.append(gene)
				summary_result.append(var)
			# I/II��������
			if var["level"] in ["5", "4"]:
				if var not in var_inter:
					var_inter.append(var)

	# δ����������Ҫ����
	for gene in set(gene_list) - set(detect_gene):
		summary_result.append({
			"gene_symbol" : gene,
			"var_type" : var_dict[gene],
			"gene_num" : gene_list.index(gene)
		})

	# ����
	summary_result_sort = sorted(summary_result, key=lambda i:i["gene_num"])
	var_inter_sort = sorted(var_inter, key=lambda i:i["gene_num"])
	gene_inter_sort = []
	# ����һ����gene�Ļ������-2023.09.28
	gene_inter_sort2 = []
	# �������-2023.09.28
	for var in var_inter_sort:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene_info.get(gene) not in gene_inter_sort:
				gene_inter_sort.append(gene_info.get(gene))
		#if var["gene_inter_shrj"] not in gene_inter_sort:
		#	gene_inter_sort.append(var["gene_inter_shrj"])
			# ����һ����gene�Ļ������-2023.09.28
			if {"gene" : gene, "gene_info" : gene_info.get(gene)} not in gene_inter_sort2:
				gene_inter_sort2.append({
					"gene" : gene,
					"gene_info" : gene_info.get(gene)
				})
			# �������-2023.09.28

	return summary_result_sort, var_inter_sort, gene_inter_sort, gene_inter_sort2

# XW0417
def sort_for_xw0417(var_list):
	result = {
		'BICC1' : '',
		'CASP7' : '',
		'TACC3v1' : '',
		'TACC3v3' : '',
		'BAI' : '',
		'R248' : '',
		'G370' : '',
		'S249' : '',
		'Y373' : ''
	}
	for var in var_list:
		if var['bio_category'] in ['Sv', 'PSeqRnaSv']:
			if var['gene_symbol'] == 'FGFR2' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'BICC1' and var['three_prime_cds'] == 'exon3':
				result["BICC1"] = 'BICC1'
			elif var['gene_symbol'] == 'FGFR2' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'CASP7' and var['three_prime_cds'] == 'exon2':
				result['CASP7'] = 'CASP7'
			elif var['gene_symbol'] == 'FGFR3' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'TACC3' and var['three_prime_cds'] == 'exon11':
				result['TACC3v1'] = 'v1'
			elif var['gene_symbol'] == 'FGFR3' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'TACC3' and var['three_prime_cds'] == 'exon10':
				result['TACC3v3'] = 'v3'
			elif var['gene_symbol'] == 'FGFR3' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'BAIAP2L1' and var['three_prime_cds'] == 'exon2':
				result['BAI'] = 'BAI'
		
#		if var['gene_symbol'] == 'FGFR3' and var['bio_category'] == 'Snvindel':
#			if var['hgvs_c'] == 'c.742C>T':
#				result['R248'] = 'T'
#			elif var['hgvs_c'] == 'c.1108G>T':
#				result['G370'] = 'T'
#			elif var['hgvs_c'] == 'c.746C>G':
#				result['S249'] = 'T'
#			elif var['hgvs_c'] == 'c.1118A>G':
#				result['Y373'] = 'T'
		# V4����
		if var['gene_symbol'] == 'FGFR3' and var['bio_category'] == 'Snvindel':
			if var['hgvs_p'] == 'p.(R248C)':
				result["R248"] = var['hgvs_c']
			elif var['hgvs_p'] == 'p.(G370C)':
				result["G370"] = var['hgvs_c']
			elif var['hgvs_p'] == 'p.(S249C)':
				result["S249"] = var['hgvs_c']
			elif var['hgvs_p'] == 'p.(Y373C)':
				result["Y373"] = var['hgvs_c']
	
	tmp_dict = {key:value for key, value in result.items() if value}
	result['final'] = 'T' if tmp_dict else 'F'

	return result

def getSum_for_tXW6002(var_list):
	t_XW6002_HRR_genelist = ["ABRAXAS1", "AKT1", "AKT2", "AR", "ATM", "ATR", "ATRX", "AURKA", "BARD1", "BRAF", "BRCA1", "BRCA2", "BRIP1", "CD274", "CHEK1", "CHEK2", "EGFR", "ERBB2",
    "ERCC1", "FANCA", "FANCC", "FANCD2", "FANCL", "FGFR1", "GEN1", "GSTP1", "IFNGR1", "KRAS", "MRE11", "MSH6", "MYB", "NBN", "PALB2", "POLD1", "POLE", "PRKACA", "RAD50", "RAD51B",
	"RAD51C", "RAD51D", "RAD54L", "RBM10", "TP53", "TSC1", "TSC2", "VEGFA"]
	hrr_gene = []
	other_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in t_XW6002_HRR_genelist:
				hrr_gene.append(gene)
			else:
				other_gene.append(gene)
	hrr_gene = sorted(set(hrr_gene))
	other_gene = sorted(set(other_gene))
	result = {}
	if hrr_gene and "AR" in hrr_gene:
		result['AR'] = 'AR'
		hrr_gene.remove('AR')
		result['HRR'] = '��'.join(hrr_gene) if hrr_gene else ''
		result['other'] = '��'.join(other_gene) if other_gene else ''
	else:
		result['AR'] = ''
		result['HRR'] = '��'.join(hrr_gene) if hrr_gene else ''
		result['other'] = '��'.join(other_gene) if other_gene else ''
	return result

def getSum_for_gXW6002(var_list):
	g_XW6002_HRR_genelist = ["ABRAXAS1", "AKT1", "AKT2", "AR", "ATM", "ATR", "AURKA", "BARD1", "BRAF", "BRCA1", "BRCA2", "BRIP1", "CD274", "CHEK1", "CHEK2", "EGFR", "ERBB2", "ERCC1",
    "ERCC3", "ERCC4", "FANCA", "FANCD2", "FANCL", "FANCM", "FGFR1", "GEN1", "KRAS", "MAPK1", "MRE11", "MSH6", "NBN", "NPM1", "PALB2", "POLD1", "POLE", "RAD50","RAD51", "RAD51B",
	"RAD51C", "RAD51D", "RAD52", "RAD54L", "SETD2", "SLX4", "TP53", "TSC1", "TSC2", "XRCC1", "XRCC2"]
	hrr_gene = []
	other_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in g_XW6002_HRR_genelist:
				hrr_gene.append(gene)
			else:
				other_gene.append(gene)
	hrr_gene = sorted(set(hrr_gene))
	other_gene = sorted(set(other_gene))
	result = {}
	if hrr_gene and "AR" in hrr_gene:
		result['AR'] = 'AR'
		hrr_gene.remove('AR')
		result['HRR'] = '��'.join(hrr_gene) if hrr_gene else ''
		result['other'] = '��'.join(other_gene) if other_gene else ''
	else:
		result['AR'] = ''
		result['HRR'] = '��'.join(hrr_gene) if hrr_gene else ''
		result['other'] = '��'.join(other_gene) if other_gene else ''
	return result

def getSum_for_XW5902(var_list, hd):
	result = ''
	if not var_list and not hd:
		result = 'δ���'
	elif var_list and hd:
		result = '���SNVIndel/HD'
	elif hd:
		result = '���HD'
	elif var_list:
		result = '���SNVIndel'
	
	return result

def var_ad3101_summary(var):
	somatic = var['var_somatic']
	germline = var['var_germline']
	variants = [
		*somatic['level_I'],
		*germline['level_5'],
		*somatic['level_II'],
		*germline['level_4'],
		*somatic['level_onco_nodrug'],
		*somatic['level_III'],
		*germline['level_3']
	]

	counter = defaultdict(int)
	for item in variants:
		g = item.get('clinic_ad3101_g', 0)
		s = item.get('clinic_ad3101_s', 0)

		if g >= 3 and not s:
			key = f'g{g}'
		elif s >= 3:
			key = f's{s}'
		else:
			continue
		counter[key] += 1

	return {
		'all' : len(variants),
		'num_5' : counter.get('g5', 0) + counter.get('s5', 0),
		'num_4' : counter.get('g4', 0) + counter.get('s4', 0),
		'num_3' : counter.get('g3', 0) + counter.get('s3', 0)
	}

def XW6701_summary(var_list, main_gene, hrr_genelist):
	'''
	����XW6701������ܽ�
	'''
	result = {'BRCA': [], 'hrr': [], 'other': []}
	allowed_bio = ['Snvindel', 'Cnv', 'Sv', 'PSeqRnaSv']
	seen = {'BRCA': set(), 'hrr': set(), 'other': set()}  # �����������Ŀ�������ظ�
	clin_map = {'Pathogenic': '�²��Ա���', 'Likely pathogenic': '�����²��Ա���'}
	func_map = {'Oncogenic': '�°��Ա���', 'Likely oncogenic': '�����°��Ա���'}
	sort_order = {
        '�²��Ա���': 0,
        '�°��Ա���': 1,
        '�����²��Ա���': 2,
        '�����°��Ա���': 3
    }

	def get_var_info(item):
		bio_category = item['bio_category']
		if bio_category == 'Snvindel':
			return f"{item['hgvs_c']} {item['hgvs_p']}"
		elif bio_category == 'Cnv':
			return item['cnv_type']
		elif bio_category in ['Sv', 'PSeqRnaSv']:
			return f"{item['five_prime_gene']}-{item['three_prime_gene']}"
		else:
			return ''

	for item in var_list:
		gene_symbol = item['gene_symbol']
		bio_category = item['bio_category']
		clin_sig = item['clinical_significance']
		func_class = item['function_classification']

		if bio_category not in allowed_bio:
			continue

		processed = False
		for category, genes in [('BRCA', main_gene), ('hrr', hrr_genelist)]:
			if gene_symbol in genes and clin_sig in clin_map:
				inter = clin_map[clin_sig]
				var_info = get_var_info(item)
				entry_key = (gene_symbol, inter, var_info)

				if entry_key not in seen[category]:
					seen[category].add(entry_key)
					result[category].append({
						'gene_symbol': gene_symbol,
						'inter': inter,
						'var_info': var_info
					})
				processed = True
				break
		if processed:
			continue
		
		inter = None
		if clin_sig != '-' and func_class == '-' and clin_sig in clin_map:
			inter = clin_map[clin_sig]
		elif func_class in func_map:
			inter = func_map[func_class]
		if inter:
			var_info = get_var_info(item)
			entry_key = (gene_symbol, inter, var_info)
			if entry_key not in seen['other']:
				seen['other'].add(entry_key)
				result['other'].append({
					'gene_symbol': gene_symbol,
					'inter': inter,
					'var_info': var_info
				})
	
	for category in result:
		result[category] = sorted(result[category], key=lambda x: (sort_order[x['inter']], x['gene_symbol'].upper()))
		
	return result

def getSum_for_tXW6701(var_list):
	var = var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug']
	main_gene = ['BRCA1', 'BRCA2', 'CDK12', 'PALB2']
	hrr_genelist = [
    "ARID1A", "ARID1B", "ARID2", "ATM", "ATR", "BAP1", "BARD1","BRIP1", "CCND1", "CCNE1", "CDK4",
	"CHEK1", "CHEK2", "ERCC1", "ERCC2", "FANCA", "FANCC", "FANCD2","FANCL", "GEN1", "HDAC2", "MDM2", "MDM4", "MLH1", "MLH3",
    "MSH2", "MSH3", "MSH6", "MTOR", "MUTYH", "NBN","PBRM1", "PMS2", "POLD1", "POLE", "RAD50", "RAD51B", "RAD51C",
    "RAD51D", "RAD54L", "SMARCA2", "SMARCA4", "SMARCB1", "TERT"]
	result = XW6701_summary(var, main_gene, hrr_genelist)

	return result

def getSum_for_gXW6701(var_list):
	var = var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug']
	main_gene = ['BRCA1', 'BRCA2', 'CDK12', 'PALB2']
	hrr_genelist = [
    "ARID1A", "ATM", "ATR", "BAP1", "BARD1", "BRIP1", "CCND1", "CCND2", "CCNE1", "CDK4", "CDKN2A",
    "CHEK1", "CHEK2", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "FANCA", "FANCD2", "FANCL", "FANCM", "GEN1", "HDAC2", "MDM2", "MDM4",
    "MLH1", "MSH2", "MSH6", "MTOR", "MUTYH", "NBN", "NPM1", "PMS2", "POLD1", "POLE", "RAD50", "RAD51", "RAD51B",
    "RAD51C", "RAD51D", "RAD52", "RAD54L", "SLX4", "SMARCA4", "TERT", "TP53", "XRCC1", "XRCC2"]
	result = XW6701_summary(var, main_gene, hrr_genelist)

	return result

def getSum_for_XW6003(var_list, hd):
	result = {'MTAP' : set(), 'CDKN2A' : set(), 'CDKN2B' : set(), 'other' : set()}
	target_genes = {'MTAP', 'CDKN2A', 'CDKN2B'}
	var_list = var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug'] + var_list['var_somatic']['level_III']

	for var in var_list:
		gene = var['gene_symbol']
		if gene in target_genes:
			result[gene].add(var['bio_category'])
		else:
			result['other'].add(gene)
	if hd:
		for item in hd:
			gene = item['gene_symbol']
			if gene in target_genes:
				result[gene].add('HD������ȱʧ��')
	
	return {
		'MTAP': '/'.join(sorted(result['MTAP'])),
		'CDKN2A': '/'.join(sorted(result['CDKN2A'])),
		'CDKN2B': '/'.join(sorted(result['CDKN2B'])),
		'other': '��'.join(sorted(result['other']))
	}