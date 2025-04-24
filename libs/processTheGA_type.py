#-*- coding:gbk -*-
import re
from libs.getEvi import varRegimen
from libs import listResultToDict
from libs.rule import judge_var, get_varsimpleinfo
import copy

'''
Discription
	
	����θ�����ӷ��͸�ʽ��
	��Ҫebv_type��MSI�Ͳ��ֻ������� 

'''

def process_ga_type(jsonDict, var_data):
	ga_result = {}
	ebv_gene_list = ["PIK3CA", "ARID1A", "BCOR"]
	gs_gene_list = ["CDH1", "RHOA", "CLDN18"]
	cin_gene_list = ["TP53", "ERBB2", "KRAS", "EGFR", "CDK6", "CCNE1", "APC", "CTNNB1", "SMAD2", "SMAD4", "PTEN"]
	ebv_type_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["ebv_type"])) if "ebv_type" in jsonDict.keys() and jsonDict["ebv_type"] else {}
	# θ�����ӷ����У�EBV��MSI��GS��CIN���֣�θ�����ӷ�����ر�־�ﲿ�ּ�⵽�ķ��Ͷ�չʾ
	# ���С�Ჿ��
	# ��⵽EBV��MSI��չʾEBV��MSI��ͬʱ������չʾ��
	# δ��⵽EBV��MSIʱ��������GS��CIN���ͣ���չʾ��ͬʱ������չʾ��
	### ����δ�޶���������
	# չʾ��ϸ��I/II/����������չ��ر������ϵ�²�/�����²�����
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# ��������б����ؼ����
	def get_var(gene_list):
		result_dict = {}
		for var in level_12_var:
			for gene in re.split(",", var["gene_symbol"]):
				if gene in gene_list:
					if gene not in result_dict.keys():
						result_dict.setdefault(gene, [])
					var_info = get_varsimpleinfo(var)
					result_dict[gene].append(var_info)
		return result_dict, ", ".join(["{0} {1}".format(k, v) for k, v in result_dict.items()])

	ebv_gene_result, ebv_gene_str = get_var(ebv_gene_list)
	gs_gene_list, gs_gene_str = get_var(gs_gene_list)
	cin_gene_list, cin_gene_str = get_var(cin_gene_list)
	ga_result = {
		"ebv_type" : ebv_type_dict,
		"ebv_gene" : ebv_gene_result,
		"ebv_sum" : ebv_gene_str,
		"gs_gene" : gs_gene_list,
		"gs_sum" : gs_gene_str,
		"cin_gene" : cin_gene_list,
		"cin_sum" : cin_gene_str
	}
	
	return ga_result