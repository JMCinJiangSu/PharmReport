#-*- coding:gbk -*-
import re
from libs.getConfig import clinicalNumStran, functionNumStran
from libs.getEvi import varRegimen
from libs.specialRequest import varInfo_FJZL
from libs.rule import S_function
import copy
from libs.rule import decimal_float, decimal_percen
from libs.getConfig import get_gene_class

'''
Discription
	
	����cnv��ʽ�� 

'''

def process_cnv(jsonDict, config):
	cnv = copy.deepcopy(jsonDict["cnv"])
	for var in cnv:
		#var["cn_mean"] = format(round(float(var["cn_mean"]), 2)+0.00, ".2f") if "cn_mean" in var.keys() and var["cn_mean"] else 0
		var["cn_mean"] = decimal_float(var["cn_mean"]) if "cn_mean" in var.keys() and var["cn_mean"] else 0
		
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
		
		#AD3101 ����һ���ֶΣ�ֻ���²��Խ���ķ����²��Եȼ�������������°��Խ���ȼ����,2025��3��6��
		var['clinic_ad3101_g'] = var['clinic_num_g'] if var["clinical_significance"] != "-" and var['function_classification'] == '-' else 3
		var['clinic_ad3101_s'] = var['clinic_num_s'] if var['function_classification'] != '-' else 3

		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["clinic_num_s"], var["top_level"] = S_function(var)
		# �������������ر���Ƶ�������Ϣ����Դ���ñ������Ʒ�������
		var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
	# ��cn_mean������
	cnv = sorted(cnv, key=lambda i:float(i["cn_mean"]), reverse=True)

	return cnv