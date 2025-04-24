#-*- coding:gbk -*-
import copy
from libs import listResultToDict
'''
Discription
	
	�ýű�������ȡPD-L1�����
	��ȡ���ݣ�
	1. ���������/���ԣ��ΰ�TPS/��������CPS��
	2. ͼƬ
	3. ���Ʒ�����δȷ���費��Ҫ��ȡ��
	����PD-L1���⴦������ 

'''

def getPDL1(jsonDict):
	pdl1_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["pdl1"]))
	if pdl1_dict:
		pdl1_dict["value"] = str(int(float(pdl1_dict["value"])*100))+"%" if pdl1_dict["result"] == "����" and pdl1_dict["type"] == "TPS" else int(pdl1_dict["value"]) if pdl1_dict["result"] == "����" and pdl1_dict["type"] == "CPS" else pdl1_dict["value"]

	return pdl1_dict