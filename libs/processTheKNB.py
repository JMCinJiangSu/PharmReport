#-*- coding:gbk -*-
from libs.getEvi import varRegimen
import copy

'''
Discription
	
	����knb��ʽ�� 

'''

def process_knb(jsonDict, config):
	knb = copy.deepcopy(jsonDict["knb"])[0] if jsonDict["knb"] else {}
	if knb:
		knb["evi_sum"] = varRegimen(jsonDict, knb["evi_sum"], config, knb)

	return knb