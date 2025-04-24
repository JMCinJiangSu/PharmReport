#-*- coding:gbk -*-

'''
Discription
	
	�ϲ�MET 14��Ծͻ�䡣
	RNA MET�ںϺ�DNA MET 14��Ծ 

'''
# �ж��Ƿ����������֤��
def judge_Predictive(var):
	if var["evi_sum"]["evi_split"] and "Predictive" in var["evi_sum"]["evi_split"]:
		return 1
	else:
		return 0

def mergeVar(var_data):
	judge_MET14 = []
	judge_METSV = []
	# ����һ������ֵ�������ж������Ƿ����DNA RNA����MET 14-2022.09.20
	judge_mergeMET = ""
	# 1. �ж��Ƿ����RNA�ں�MET-MET�������ռ���ͻ������Ʒ����͵ȼ�
	for var in var_data:
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET" and judge_Predictive(var):
			judge_METSV = var
			METSV_regimen = [(i["regimen_name"], i["evi_conclusion"]) for i in var["evi_sum"]["evi_split"]["Predictive"]]
	# 2. ������RNA�ں�MET-MET�����ж�Snvindel���Ƿ����MET 14
	if judge_METSV:
		for var in var_data:
			regimen = [(i["regimen_name"], i["evi_conclusion"]) for i in var["evi_sum"]["evi_split"]["Predictive"]] if judge_Predictive(var) else []
			# Snvindel��������ΪMET�������Ʒ����͵ȼ���MET-METһ�£����ж�ΪMET 14
			if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "MET" and regimen == METSV_regimen:
				judge_MET14 = var
				# ԭ�����У�ɾ��Snvindel��SV�е�MET 14�������������ٲ�����
				var_data.remove(var)
				var_data.remove(judge_METSV)
				# ����MET 14
				judge_MET14["hgvs_p_2"] = "MET-MET�ں�"
				judge_MET14["freq_2"] = str(judge_METSV["copies"]) + " copies" if "copies" in judge_METSV.keys() else \
										str(int(judge_METSV["reads"])) + " copies" if "reads" in judge_METSV.keys() else \
										""
				judge_MET14["judge_mergeMET"] = "yes"
				judge_mergeMET = "yes"
				# ����gene:exon-gene:exon-2022.09.27
				five_prime_region = judge_METSV["five_prime_region"] if "five_prime_region" in judge_METSV.keys() and judge_METSV["five_prime_region"] else \
									judge_METSV["five_prime_cds"]
				three_prime_region = judge_METSV["three_prime_region"] if "three_prime_region" in judge_METSV.keys() and judge_METSV["three_prime_region"] else \
									 judge_METSV["three_prime_cds"]
				judge_MET14["var_hgvs_2"] = "{0}:{1}-{2}:{3}".format(judge_METSV["five_prime_gene"], five_prime_region, judge_METSV["three_prime_gene"], three_prime_region)
				# �������������ͽ��-20220929
				judge_MET14["variant_desc_cn_2"] = judge_METSV["variant_desc_cn"]
				judge_MET14["variant_interpret_cn_2"] = judge_METSV["variant_interpret_cn"]
				judge_MET14["evi_sum_2"] = judge_METSV["evi_sum"]
	
	# ��MET 14�����������б��һλ
	if judge_MET14:
		var_data.insert(0, judge_MET14)

	return var_data, judge_mergeMET