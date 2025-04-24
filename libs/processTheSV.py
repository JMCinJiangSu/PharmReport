#-*- coding:gbk -*-
import re
from libs.getConfig import clinicalNumStran, functionNumStran
from libs.getEvi import varRegimen
import copy
from itertools import groupby
from libs.specialRequest import varInfo_FJZL
from libs.rule import S_function
from libs.rule import decimal_float, decimal_percen

'''
Discription	
	����sv��ʽ��
	����sv ��rna_sv
	�����ʽ�������������������ݣ�
	1. �ں��������ж�
	   ͨ����ⷶΧ��5�˻����ڷ�Χ�ڣ���3�˲��ڣ���������Ϊ5�˻������������Ĭ����3��Ϊ������(ȡ����ֱ��ʹ��json���ص�gene_symbol��2022.07.26)
	2. ����ϲ���������CP40 sv
	3. RNA SV ��DNA SVƥ�� 
'''

# ����ƴ���ںϣ�gene1:exon1-gene2:exon2
def var_info(var):
	return var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]

def process_sv(jsonDict, config):
	sv = copy.deepcopy(jsonDict["sv"])
	rna_sv = copy.deepcopy(jsonDict["rna_sv"])
	# sv����
	for var in sv:
		# ���������£��°��Ժ��²���ֻ����1��-2023.05.22
		#var["clinic_num_g"] = clinicalNumStran().get(var["clinical_significance"], 3) 
		#var["clinic_num_s"] = functionNumStran().get(var["function_classification"], 3) 
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
		#AD3101 ����һ���ֶΣ�ֻ���²��Խ���ķ����²��Եȼ�������������°��Խ���ȼ����,2025��3��6��
		var['clinic_ad3101_g'] = var['clinic_num_g'] if var["clinical_significance"] != "-" and var['function_classification'] == '-' else 3
		var['clinic_ad3101_s'] = var['clinic_num_s'] if var['function_classification'] != '-' else 3
		
		# 2023.05.22 �������
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["clinic_num_s"], var["top_level"] = S_function(var)
		# three_prime_cds �� five_prime_cds������ȡ������json��ɾ���ˡ�-ins��������
		var["five_prime_cds"] = "-".join(re.split("-", (re.split(":", var["var_hgvs"])[2]))[:-1]) \
								if not re.search("--", var["var_hgvs"]) else re.split("_", (re.split("--", var["var_hgvs"])[0]))[-1]
		var["three_prime_cds"] = re.split(":", var["var_hgvs"])[-1] if not re.search("--", var["var_hgvs"]) else re.split("_", (re.split("--", var["var_hgvs"])[1]))[-1]

		# ��һ������-2023.10.19
		# var_hgvs�¸�ʽ��gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, �ɵ�Ϊgene1:NM_xxx_exon1--gene2:NM_xxx_exon2
		# cds����xxx:exon1��xxx:exon2
		var["five_prime_cds"] = re.split(":", var["five_prime_cds"])[-1] if re.search(":", var["five_prime_cds"]) else var["five_prime_cds"]
		var["three_prime_cds"] = re.split(":", var["three_prime_cds"])[-1] if re.search(":", var["three_prime_cds"]) else var["three_prime_cds"]
		# �������-2023.10.19

		# ������ת¼��
		tmp_dict = {var["five_prime_gene"] : var["five_prime_transcript"], var["three_prime_gene"] : var["three_prime_transcript"]}
		var["transcript_primary"] = tmp_dict.get(var["gene_symbol"], "")
		# ����Ƶ�� CP40��copies��������Ŀ��freq
#		var["freq_str"] = "{:.2%}".format(float(var["freq"])) if "freq" in var.keys() and var["freq"] else ""

		# freq ���ݰٷ�����ʽ����ת��ΪС����ʽ-����v4-��쿷ң�2023.11.09
		freq_key = ["freq"]
		for key in freq_key:
			if key in var.keys():
				var[key] = float(var[key].replace("%", ""))/100 if re.search("%", str(var[key])) else var[key]
		# ���ݽ���-2023.11.09

		var["freq_str"] = decimal_percen(var["freq"]) if "freq" in var.keys() and var["freq"] else ""
		var["copies"] = var["copies"] if "copies" in var.keys() and var["copies"] else ""
		# �������������ر���Ƶ�������Ϣ����Դ���ñ������Ʒ�������
		var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
		# ���������ԭ�ԡ�,���ָ��������ԡ�\n���ָ�����������Ǵ�����������-2023.03.31
		var["gene_symbol"] = var["gene_symbol"].replace("\n", ",")
	# sv�ϲ���� ��CP40/CRC12Ҫ�ϲ���������Ŀֱ��ʹ��sv������
	sv_combination = combinationSV(sv) if re.search("Classic|CRC12", jsonDict["sample_info"]["prod_names"]) else sv

	# rna_sv����
	for var in rna_sv:
		# ���������£��°��Ժ��²���ֻ����1��-2023.05.22
		#var["clinic_num_g"] = clinicalNumStran().get(var["clinical_significance"], 3) 
		#var["clinic_num_s"] = functionNumStran().get(var["function_classification"], 3) 
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
		#AD3101 ����һ���ֶΣ�ֻ���²��Խ���ķ����²��Եȼ�������������°��Խ���ȼ����,2025��3��6��
		var['clinic_ad3101_g'] = var['clinic_num_g'] if var["clinical_significance"] != "-" and var['function_classification'] == '-' else 3
		var['clinic_ad3101_s'] = var['clinic_num_s'] if var['function_classification'] != '-' else 3
		# 2023.05.22 �������
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["clinic_num_s"], var["top_level"] = S_function(var)
		# ����three_prime_cds �� five_prime_cds
		var["five_prime_cds"] = re.split(":", var["five_prime_cds"])[0]
		var["three_prime_cds"] = re.split(":", var["three_prime_cds"])[0]
		# Ƶ�ʣ�����ģ����Ҫչʾ��
		var["freq"] = int(float(var["supp_splt_reads"]) + float(var["supp_span_reads"]))
		var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
		# ���������ԭ�ԡ�,���ָ��������ԡ�\n���ָ�����������Ǵ�����������-2023.03.31
		var["gene_symbol"] = var["gene_symbol"].replace("\n", ",")
		# ����v4��v3��reads�ֶ�v4��������-2023.11.06
		var["reads"] = int(float(var["supp_splt_reads"]) + float(var["supp_span_reads"])) if "reads" not in var.keys() else var["reads"]
		# �������-2023.11.06

	# rna sv ��dna svƥ��
	sv_combination_match, rna_sv_only = matchRNA_DNA_sv(sv_combination, rna_sv)

	# ��Ƶ��������
	# Myeloid����CP 2023��12��27�� jmc
	# xw1402 OncoPro(��֯)����ͬhandle jmc 2024��8��16��
	if not re.search("Classic|CRC12|Myeloid|OncoPro����֯��", jsonDict["sample_info"]["prod_names"]):
#		print (sv_combination_match)
		sv_combination_match = sorted(sv_combination_match, key=lambda i:float(str(i["freq"]).replace("%","")), reverse=True)
	else:
		sv_combination_match = sorted(sv_combination_match, key=lambda i:float(i["copies"]), reverse=True)
	
	rna_sv = sorted(rna_sv, key=lambda i:float(i["freq"]), reverse=True)
	rna_sv_only = sorted(rna_sv_only, key=lambda i:float(i["freq"]), reverse=True)

	return sv_combination_match, rna_sv, rna_sv_only

def combinationSV(sv_data):
	# �Ӹ�����ʹ���˻�����ͬ�ı����ŵ�һ��groupby����ABBA����ı��촦���е�����-2023.05.25
	sv_data = sorted(sv_data, key=lambda i:(i["five_prime_gene"], i["three_prime_gene"]))
	# 2023.05.25-�������
	# ��������CP40��CRC12
	for var in sv_data:
		var["hgvs_p"] = var["five_prime_gene"] + "-" + var["three_prime_gene"]
	sv_sum = groupby(sv_data, key=lambda x : x["hgvs_p"])
	sv_combination = []
	for i, j in sv_sum:
		sv_group = list(j)
		major_sv = sv_group[0]
		# ��Ӹ����ں��ͼ���Ӧ�Ŀ���������������������CP40-2022.08.23
		major_sv["sv_YNZL"] = []
		for k in sv_group:
			major_sv["sv_YNZL"].append({
				"five_prime_gene" : k["five_prime_gene"],
				"five_prime_cds" : k["five_prime_cds"],
				"three_prime_gene" : k["three_prime_gene"],
				"three_prime_cds" : k["three_prime_cds"],
				"copies" : k["copies"] if "copies" in k.keys() and k["copies"] else ""
				})
		# ������-2022.08.23
		major_sv['copies'] = sum(int(k["copies"]) for k in sv_group if k["copies"])
		# �ں�ģʽƴ��
		major_sv["var_desc_merge"] = "��".join([var_info(k) for k in sv_group])
		# �����ں������б����ʽ�棬���㱨��չʾ
		major_sv["merge_sv_list"] = [var_info(k) for k in sv_group]
	
		sv_combination.append(major_sv)

	return sv_combination

def matchRNA_DNA_sv(sv, rna_sv):
	rna_sv_dict = {}
	# 301 MP����DNA1��+RNA������exonͬ���ϵ㲻ͬ�����������������ݣ��������ֵ�-2023.04.12
	dup_rna_sv_dict = {}
	for var in rna_sv:
		rna_sv_dict.setdefault(var_info(var), {})
		rna_sv_dict[var_info(var)] = var
	# �����ֵ�-2023.04.12
		dup_rna_sv_dict.setdefault(var_info(var), [])
		dup_rna_sv_dict[var_info(var)].append(var)
	
	# ������RNA exon��ͬ���ϵ㲻ͬ�����ʱ����Ӹ���ǩ�����ж��Ƿ�չʾ����ϵ�-2023.04.12
	for var in rna_sv:
		key = var_info(var)
		if key in dup_rna_sv_dict.keys():
			if len(dup_rna_sv_dict[key]) > 1:
				var["dup_rna"] = "yes"
	
	rna_sv_only = copy.deepcopy(rna_sv)
	rna_sv_pop_key = []
	for var in sv:
		if var_info(var) in rna_sv_dict.keys():
			var["rna_detect"] = rna_sv_dict[var_info(var)]
			# �����ֵ�-2023.04.12
			var["dup_rna_detect"] = dup_rna_sv_dict[var_info(var)]
			rna_sv_pop_key.append(var_info(var))

	# ƥ����ģ�rna_sv_only��DNA/RNA�������ɾ��
	# ������£�ԭremove�޷����ظ���ȫ��ɾ������Ϊ����Ĵ���-2023.02.28
	rna_sv_only_reduce = []
	for var in rna_sv_only:
		if var_info(var) not in rna_sv_pop_key:
			rna_sv_only_reduce.append(var)

#	for var in rna_sv_only:
#		if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"] in rna_sv_pop_key:
#			rna_sv_only.remove(var)

	return sv, rna_sv_only_reduce
	