#-*- coding:gbk -*-
import copy
from libs import listResultToDict
from libs.getEvi import varRegimen

'''
Discription
	
	�ýű�������ȡMSI�����
	��ȡ���ݣ�
	1. ���
	2. ͼƬ
	3. ���Ʒ���
	����MSI���⴦������ 

'''

def getMSI(jsonDict, config):
	msi_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["msi"]))
	# MSI_num json�޷��أ�ͨ���ű�����
	if msi_dict:
		msi_dict["msi_num_cp40"] = int(float(msi_dict["msi_score"]) * 55 + 0.5) if msi_dict["msi_score"] or msi_dict["msi_score"]==0 else ""
	if "evi_sum" in msi_dict.keys() and msi_dict["evi_sum"]:
		msi_dict["evi_sum"] = varRegimen(jsonDict, msi_dict["evi_sum"], config, msi_dict)
	return msi_dict