#-*- coding:gbk -*-
import datetime
import copy
import re

'''
Discription
	
	��ȡjson�ļ��е�sample info��ת��Ϊ����ģ�巽�����ĸ�ʽ�� 

'''

def getSample(jsonDict):
	data = copy.deepcopy(jsonDict["sample_info"])
	data["report_date"] = str(datetime.date.today())
	for k, v in data.items():
		if not v:
			data[k] = ""
	
	#�������ڸ�ʽ��XX��XX��XX��
	receive_date = re.split("-", data["receive_data"])
	data["receive_date_special_1"] = receive_date[0] + "��" + receive_date[1] + "��" + receive_date[2] + "��" if receive_date and len(receive_date) >= 3 else ""
	report_date = re.split("-", data["report_date"])
	data["report_date_special_1"] = report_date[0] + "��" + report_date[1] + "��" + report_date[2] + "��" if report_date and len(report_date) >= 3 else ""

	#�������ŷ���ʱ��
	json_name_list = re.split("_", data["json_batch_name"])
	json_date = json_name_list[0] if json_name_list else ""
	data["json_date"] = json_date[:4] + "��" + json_date[4:6] + "��" + json_date[6:] + "��" if json_date and len(json_date) == 8 else ""

	#�ж��ټ����������Ž��ջ����Ϻ����յģ�������ĸ�Ƿ��S�����ж�
	data["locate"] = "SH" if re.match("S", data["sample_id"]) else "XM"

	# �㽭���������ͱ�ע����
	data["ZJFB_mark_dict"] = {}
	if data["mark"]:
		for i in re.split("��|;", data["mark"]):
			tmp = re.split("��|:", i)
			if len(tmp) == 2:
				data["ZJFB_mark_dict"][tmp[0]] = ""
				data["ZJFB_mark_dict"][tmp[0]] = tmp[1]
	
	# ����ϸ������������ֵ�ֶΣ����ڱȶ�
	data["tumor_content_num"] = data["tumor_content"].replace("%", "") if "tumor_content" in data.keys() and data["tumor_content"] else ""

	# ���⴦���ټ���Ŀ������Ϊ0�ķſգ��Ա�Ϊδ֪�ķſ�-20220923
	if jsonDict["sample_info"]["report_module_type"] != "hospital":
		data["age"] = "" if data["age"] and data["age"] == "0" else data["age"]
		data["gender"] = "" if data["gender"] and data["gender"] == "δ֪" else data["gender"]

	# ��۴�ѧ����ҽԺ�����ͱ�ע����-2022.10.10
	data["XGDX_mark_dict"] = {}
	if data["mark"]:
		for i in re.split("��", data["mark"]):
			tmp = re.split("��", i)
			if tmp:
				data["XGDX_mark_dict"][tmp[0]] = tmp[1] if len(tmp) >= 2 else ""
	
	# �ɼ����ڣ���֯�ɼ����ڻ���ԭ�����ֶΣ�gather_data
	# ����һ��ѪҺ�ɼ����ڣ����������-2022.11.29
#	data["blood_collection_date"] = data["blood_collection_date"] if "blood_collection_date" in data.keys() and data["blood_collection_date"] else ""
	data["blood_collection_date"] = re.split(" ", data["blood_collection_date"])[0] if "blood_collection_date" in data.keys() and data["blood_collection_date"] else ""

	# �������ڷ��ظ�ʽ�еĻ������ʱ�䣬��߼Ӹ����ݴ���-��쿷�-20231025
	date_key = ["receive_data", "tissue_collection_date", "blood_collection_date", "gather_data", "submission_date", "application_date", "sampling_time", "analysis_date", "section_date", "blood_date_received", "tissue_date_received"]
	for key in date_key:
		if key in data.keys():
			data[key] = re.split(" ", data[key])[0] if data[key] else ""
	# ������-20231025

	# ��������-��֯�ɼ����ڣ�Ŀǰ�󲿷�ģ�嶼��ʹ��	gather_data�ֶΣ�������һ��tissue_collection_date-��쿷ң�2023.11.14
	# ������gather_data����gather_dataʹ�ø��ֶ����ݣ���������gather_data����ʹ��tissue_collection_date
	data["gather_data"] = data["gather_data"] if "gather_data" in data.keys() and data["gather_data"] else \
						  data["tissue_collection_date"] if "tissue_collection_date" in data.keys() and data["tissue_collection_date"] else \
						  ""
	# �������-2023.11.14
	
	return data