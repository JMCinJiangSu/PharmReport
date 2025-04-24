#-*- coding:gbk -*-
import re
import copy
from functools import reduce
from libs.rule import decimal_float, decimal_percen

'''
Discription
	
	��ȡjson�ļ��е�qc��ת��Ϊ����ģ�巽�����ĸ�ʽ�� 

'''

# ʶ���Ƿ�Ϊ��ֵ
def is_number(i):
	try:
		float(i)
		return True
	except:
		pass 
	if i.isnumeric():
		return True
	return False

# �����жϲ����ǣ�Miseq/NextSeq������-2023.07.03
# ����1����λ����ΪMiseq��2��S��ͷ��Ϊ��������ǣ�3��A��ͷ��ΪCN500
# ������£��ⲿ����ı�ſ��ܻ���ɡ�SJMXXX��������CN500��������ж�Ϊ������߼Ӹ��ж�-2023.11.01
# ��ͷΪS+һ�����ֵ��ж�Ϊ������������"huada/cn500"�������ɫ
def judge_Sequencer(flowcell_lane):
	flowcell = re.split("_", flowcell_lane)[0] if flowcell_lane else ""
	if len(flowcell) == 5:
		return "miseq"
	elif re.match("S", flowcell):
	#	return "huada"
		if re.match("S[0-9]", flowcell):
			return "huada"
		else:
			return ""
	elif re.match("A", flowcell):
		return "cn500"
	else:
		return ""

def QCStran_dict(qcDict):
	QC_result = {}
	for k, v in qcDict.items():
		if k != "qc_type" and v:
			QC_result[k+"_num"] = float(v[:-1]) if re.search("cleandata_size|dv200|tumor_content|tumor_content_macrodissection_performed", k) else float(v) if is_number(v) else v
			#QC_result[k] = "{:.2%}".format(float(v)) if re.search("q30|q20|ratio|uni20", k) and not re.search("concentration", k) else v if re.search("cleandata_size", k) else "{:.2f}".format(float(v)) if is_number(v) else v
			QC_result[k] = decimal_percen(v) if re.search("q30|q20|ratio|uni20", k) and not re.search("concentration|registration_certificate", k) else \
						   v if re.search("cleandata_size", k) else \
						   decimal_float(v) if is_number(v) else \
						   v

		else:
			QC_result[k+"_num"] = 0
			QC_result[k] = ""
		# �����������ֶ�-2023.07.30
		if k == "flowcell_lane":
			QC_result["Sequencer"] = judge_Sequencer(v)
		# ��������-2023.07.03
	return QC_result

def QCStran_list(qclist):
	QC_result = {}
	for i in qclist:
		QC_result.setdefault(i["qc_type"], {})
		for k, v in i.items():
			if k != "qc_type" and v:
				QC_result[i["qc_type"]][k+"_num"] = float(v[:-1]) if re.search("cleandata_size|dv200|tumor_content|tumor_content_macrodissection_performed", k) else float(v) if is_number(v) else v
				#QC_result[i["qc_type"]][k] = "{:.2%}".format(float(v)) if re.search("q30|q20|ratio|uni20", k) and not re.search("concentration", k) else v if re.search("cleandata_size", k) else "{:.2f}".format(float(v)) if is_number(v) else v
				QC_result[i["qc_type"]][k] = decimal_percen(v) if re.search("q30|q20|ratio|uni20", k) and not re.search("concentration", k) else \
											 v if re.search("cleandata_size", k) else \
											 decimal_float(v) if is_number(v) else \
											 v
			else:
				QC_result[i["qc_type"]][k+"_num"] = 0
				QC_result[i["qc_type"]][k] = ""	
			# �����������ֶ�-2023.07.30
			if k == "flowcell_lane":
				QC_result["Sequencer"] = judge_Sequencer(v)
			# ��������-2023.07.03
	return QC_result
	
# ������
def getJsonQC(jsonDict):
	qc = copy.deepcopy(jsonDict["qc"])
	# ����XW0248 �ж�CNV�ʿ��Ƿ�ϸ�
	#qc['dna_data_qc'][0]['cnv_noise'] = jsonDict['cnv'][0]['cnv_noise'] if jsonDict['cnv'] else ''
	
	qc_items = [k for k in qc]
	data = {}
	for item in qc_items:
		# qc_gradient��������
		if item == "qc_gradient" and qc["qc_gradient"]:
			data["qc_gradient"] = {}
			for i in qc["qc_gradient"]:
				data["qc_gradient"][i["qc_source"]+"_"+i["gradient_num"]] = "{:.2%}".format(float(i["gradient_ratio"])) 
		else:
			# ���ʿ�����ֻ��һ����Ϊ�ֵ��ʽʱ��ʹ�ø�ת����ʽ
			if type(qc[item]).__name__=="dict":
				data[item] = QCStran_dict(qc[item])
			# ���ʿ�����Ϊ�б�
			else:
				# �б���Ϊ1ʱ��ת��Ϊ�ֶΣ�ʹ���ֵ��ת����ʽ
				if len(qc[item])==1:
					data[item] = QCStran_dict(qc[item][0])
				# �б��к��ж��Ԫ��ʱ��ʹ���б��ת����ʽ
				elif len(qc[item]) >= 2:
					# �Ƚ���ȥ�أ���ֹ�����ظ����ʿ����ݵ��±���
					qc_data = reduce(lambda x, y:x if y in x else x + [y], [[],]+qc[item])
					# ȥ�غ󻹱����ԭ��������������ʱû�а��޹�����F������
					if len(qc_data) ==1:
						data[item] = QCStran_dict(qc_data[0])
					else:
						data[item] = QCStran_list(qc[item])
	
	# ����ʪʵ���ʿ�-2022.08.10
	qc_lib = copy.deepcopy(jsonDict["lib_quality_control"]) if "lib_quality_control" in jsonDict.keys() and jsonDict["lib_quality_control"] else {}
	qc_lib_items = [k for k in qc_lib]
	lib_data = {}
	for item in qc_lib_items:
		if type(qc_lib[item]).__name__=="dict":
			lib_data[item] = QCStran_dict(qc_lib[item][0])
		else:
			if len(qc_lib[item]) == 1:
				lib_data[item] = QCStran_dict(qc_lib[item][0])
			elif len(qc_lib[item]) >= 2:
				qc_data = reduce(lambda x, y:x if y in x else x + [y], [[],]+qc_lib[item])
				if len(qc_data) == 1:
					lib_data[item] = QCStran_dict(qc_data[0])
				else:
					lib_data[item] = QCStran_list(qc_lib[item])

	return data, lib_data