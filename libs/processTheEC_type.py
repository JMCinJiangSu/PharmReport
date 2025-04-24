#-*- coding:gbk -*-
import re
from libs.getEvi import varRegimen
from libs import listResultToDict
import copy

'''
Discription
	
	�����ӹ���Ĥ�����ӷ��͸�ʽ�� 

'''

def process_ec_type(jsonDict, config):
	ec_type_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["ec_type"]))
	ec_type_dict["evi_sum"] = varRegimen(jsonDict, ec_type_dict["evi_sum"], config, ec_type_dict)
	
	return ec_type_dict