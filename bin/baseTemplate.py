#-*- coding:gbk -*-
import os
import re
import jinja2
from docxtpl import DocxTemplate
import datetime
import io
#import threading
import customize_filters

jinja_env = jinja2.Environment()

'''
Discription
	
	�������ģ�壬��Ϊ�̶�ģ��Ͷ�̬���ģ��

'''

def BaseReport(data, tpl, merge_template, report_template, json_name, outfile):
	dynamic_dir = os.path.join(report_template, "template_dynamic")
	fixed_dir = os.path.join(report_template, "template_fixed")
	data["base_report"] = {}

	result = {}

	## ��̬���ģ��
	# ��Ҫƴ�ӣ�merge_template == 0
	if merge_template["merge_template"] == "0":
		temp_dynamic_list = []
		# ��̬ģ�壬ָ���б�Ϊnoʱ�򲻴���
		if "no" not in merge_template["temp_dynamic"]:
			# ���ܶ�̬ģ��Ͷ�Ӧ����ʽ
			if merge_template["temp_dynamic"]:
				# ��ʽ��module_name.type����ֻ��module_name�Ļ���type��ΪĬ��ֵdefault
				temp_dynamic_list = [{"module_name" : re.split("\.", i)[0], "type" : re.split("\.", i)[1] if len(re.split("\.", i)) == 2 else "default"} for i in merge_template["temp_dynamic"]]
			# ��̬ģ�壬Ĭ�ϣ���Ҫ�����ļ��У�ģ����ʽҲΪĬ�ϵ�
			else:
				temp_dynamic_list = [{"module_name" : re.split("\.", i)[0], "type" : ""} for i in os.listdir(dynamic_dir)]
		
		for temp in temp_dynamic_list:
			# ����̬ģ��Ͷ�Ӧ��ʽID���ص�data���ݽṹ�У�����ָ����ʽѡ��
			data["base_report"][temp["module_name"]+"_type"] = temp["type"]
#			print (data["base_report"])
			# ��̬ģ�����
			tpl_name = DocxTemplate(os.path.join(dynamic_dir, temp["module_name"]+".docx"))
			tpl_name.render(data, jinja_env)
			# ʹ��IO BytesIOд���ڴ��У����������������
			file_stream = io.BytesIO()
			tpl_name.save(file_stream)
			result[temp["module_name"]] = tpl.new_subdoc(file_stream)
#			tpl_name.save("./tmp/"+temp["module_name"]+"subdoc.docx")
			#print (temp, datetime.datetime.now())
		

		'''
		# ���߳�-���Խ�����ٶ�û������
		def newsub(temp):
			# ����̬ģ��Ͷ�Ӧ��ʽID���ص�data���ݽṹ�У�����ָ����ʽѡ��
			data["base_report"][temp["module_name"]+"_type"] = temp["type"]
			# ��̬ģ�����
			tpl_name = DocxTemplate(os.path.join(dynamic_dir, temp["module_name"]+".docx"))
			tpl_name.render(data, jinja_env)
			# ʹ��IO BytesIOд���ڴ��У��������ʵ���ļ�
			file_stream = io.BytesIO()
			tpl_name.save(file_stream)
			result[temp["module_name"]] = tpl.new_subdoc(file_stream)
			print (temp, datetime.datetime.now())
		
		for temp in temp_dynamic_list:
			t = threading.Thread(target = newsub, args = (temp, ))
			t.start()
			print (threading.active_count())
		# ���Զ��߳�
		'''

		# �̶��б���Ҫ�������Ƿ�Ҫ��ָ���б���Ĭ�ϵĴ��������ٶȲ����Σ�
		# �̶�ģ�壬ָ���б�Ϊnoʱ������
		if merge_template["temp_fixed"]:
			temp_fixed_list = merge_template["temp_fixed"] if "no" not in merge_template["temp_fixed"] else []
		else:
			temp_fixed_list = [re.split("\.", i)[0] for i in os.listdir(fixed_dir)]
		
		for temp in temp_fixed_list:
			result[temp] = tpl.new_subdoc(fixed_dir+"/"+temp+".docx")

	return result