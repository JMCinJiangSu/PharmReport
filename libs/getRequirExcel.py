#-*- coding:gbk -*-
import os
import re
import xlrd

'''
Discription
	
	��ȡconfig�ļ����еı���ģ�������ļ���ת��Ϊ��������������õĸ�ʽ�� 
	{
		�ټ� : {
			���� : {
				*���� (ҽԺ����Ʒ1����Ŀ) ��ģ����1,
				(ҽԺ, ��Ʒ1) : ģ����1
				},
			ͨ��-��� : {
				��Ʒ1 : ģ����1
				},
			ͨ��-���� : {
				��Ʒ1 : ģ����1
				}
			},
		��Ժ : {
			���� : {
				(ҽԺ, ����, ��Ʒ1) : ģ����1,
				(ҽԺ, ��Ʒ1) : ģ����1
 				},
			ͨ�� : {
				��Ʒ1 : ģ����1
				}
		}
	
'''

#BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#requirenment_path = os.path.join(BASE_DIR, "config/report_requirenment.xlsx")

key_stran = {
	"��Ʒ����" : "prod_name",
	"���ҵ������" : "business_type",
	"ģ������" : "report_type",
	"��λ����" : "company",
	"����" : "hosp_depart",
	"ģ����" : "report_name",
	"״̬" : "status",
	"ģ�忪����" : "developer",
	"�����" : "auditors",
	"���ʱ��" : "auditors_time",
	"��ע" : "note",
	"���¼�¼" : "update",
	"�ټ�" : "rummage",
	"��Ժ" : "hospital",
	"����" : "CustomEdition",
	"ͨ��" : "Universal",
	"ͨ��-���" : "Universal_simple",
	"ͨ��-����" : "Universal_complete",
	"����ģ��ƴ��" : "merge_template",
	"����ģ��-��̬" : "temp_dynamic",
	"����ģ��-�̶�" : "temp_fixed",
	"ҩ����Ŀ����" : "product_name"
	}

def getRequirenment(config):
	#BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	requirenment_path = os.path.join(config, "report_requirenment.xlsx")
	xls = xlrd.open_workbook(requirenment_path)
	requirenment_sheet = xls.sheet_by_name("report_requirenment")
	key = requirenment_sheet.row_values(0)
	Data = []
	for num in range(1, requirenment_sheet.nrows):
		rows = requirenment_sheet.row_values(num)
		tmpdict = {}
		for i in range(len(key)):
			tmpdict[key_stran.get(key[i])] = key_stran[rows[i]] if rows[i] in key_stran.keys() else rows[i]
		Data.append(tmpdict)

	#������Ϣ�ֵ�ṹ	
	requir_dict = {
		"rummage" : {
			"CustomEdition" : {},
			"Universal_simple" : {},
			"Universal_complete" : {}
			},
		"hospital" : {
			"CustomEdition" : {},
			"Universal" : {}
			}
		}

	# ����Ŀǰ���ڽ�Ժ��������Ҫ�����֣��ټ��ҩ��Ͳ�������
	for i in Data:
		if str(int(i["status"])) == "0":
			if i["company"]:
				if i["hosp_depart"]:
					requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["hosp_depart"], i["prod_name"]), "")
					requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["hosp_depart"], i["prod_name"])] = i["report_name"]
				else:
					# ����ҩ������-2023.11.17
					if i["product_name"]:
						requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"], i["product_name"]), "")
						requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"], i["product_name"])] = i["report_name"]
					else:
					# ��������-2023.11.17
						#print (i)
						requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"]), "")
						requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"])] = i["report_name"]
			else:	
				requir_dict[i["business_type"]][i["report_type"]].setdefault(i["prod_name"], "")
				requir_dict[i["business_type"]][i["report_type"]][i["prod_name"]] = i["report_name"]

	# �Ƿ���Ҫƴ�ӻ���ģ����ձ�
	report_merge_dict = {}
	for i in Data:
		if str(int(i["status"])) == "0":
			if i["report_name"] not in report_merge_dict.keys():
				report_merge_dict.setdefault(i["report_name"],{})
			report_merge_dict[i["report_name"]] = {
				"merge_template" : str(int(i["merge_template"])),
				"temp_dynamic" : re.split(";", i["temp_dynamic"].replace("\n","").replace(" ","")) if i["temp_dynamic"] else "",
				"temp_fixed" : re.split(";", i["temp_fixed"].replace("\n","").replace(" ","")) if i["temp_fixed"] else ""
			}

	return requir_dict, report_merge_dict
