#-*- coding:gbk -*-
import copy
from libs import listResultToDict
from libs.rule import decimal_float, is_number

'''
Discription
	
	�ýű�������ȡTME�����
	��ȡ���ݣ�
	1. ���
	2. ͼƬ 

'''

def getTME(jsonDict):
	tme_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["tme_type"])) if "tme_type" in jsonDict.keys() else {}
	tme_score = copy.deepcopy(jsonDict["tme_score"]) if "tme_score" in jsonDict.keys() else {}
	if tme_score:
		for k,v in tme_score.items():
			#tme_score[k] = "{:.2f}".format(float(v))
			#tme_score[k] = decimal_float(v)
			# tme_score�����������˲�����ֵ�͵��ֶΣ�sample_id�ȣ����Ӹ�����-2023.09.20
			tme_score[k] = decimal_float(v) if is_number(v) else v
			# �������-2023.09.20

	return tme_dict, tme_score