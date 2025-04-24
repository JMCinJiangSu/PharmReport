#-*- coding:gbk -*-
import re
from libs.getRequirExcel import getRequirenment
from libs.getProdName_alias import alias_name

'''
Discription
	
	ƥ�䱨��ģ�����ñ���ȡģ�����ƣ����ڱ�����䡣 

'''

def MatchReport(jsonDict, config):
	requir_dict, report_merge_dict = getRequirenment(config)
	i = jsonDict.get("sample_info")
	# ��ȡ��Ʒ������ת��Ϊ��Ҫ��Ʒ��������ƥ��ģ��
	alias_name_dict = alias_name(config)
	i["prod_names"] = alias_name_dict[i["prod_names"]] if i["prod_names"] in alias_name_dict.keys() else i["prod_names"]
	#ƥ�䱨��ģ��
	# 1. �ټ�
	 # ����ҩ������-ƥ���ǰ��-2023.11.17
	 # 1.0 ��ƥ��"�ټ�:����:ҽԺ-��Ʒ-��Ŀ"
	 # 1.1 ��ƥ��"�ټ�:����:ҽԺ-��Ʒ"
	 # 1.2 �޶��ƵĻ�����ƥ��"�ټ�:ͨ��-���-��Ʒ"�� ����ƥ��"�ټ�:ͨ��-����-��Ʒ"
	# 2. ��Ժ
	 # 2.1 ��ƥ�䶨��"��Ժ:����:ҽԺ-����-��Ʒ"���޵Ļ�ƥ��"��Ժ:����:ҽԺ-��Ʒ"
	 # 2.2 �޶��ƵĻ�ƥ��"��Ժ:ͨ��:��Ʒ"
	report_name = ""
	if i["report_module_type"] == "rummage":
		# 20220627-����������ݣ��ټ�������Ҫƥ�����׺�ͼ쵥λ���ƣ����ֶ��ϴ��Ķ���ȱʧ�������ݣ�����Ҫ���ֶ��ϴ��������ټ챨�棬�ᱨ��
		# ��߼��ݳ�ȱ�ٴ���׺�ͼ쵥λ����ʱ��ֱ�ӳ��ټ������汨��
		i["origin_company"] = i["origin_company"] if i["origin_company"] else i["company"]
		# �������-20220627
		# ����ҩ������-2023.11.17����λ��������company�ɣ����ڻ����ֶ��ϴ��Ķ�����������Խ�LIMS�����ʲô��
		if (i["company"], i["prod_names"], i["product_name"]) in requir_dict["rummage"]["CustomEdition"].keys():
			report_name = requir_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["product_name"]))
		else:
		# ����ҩ�����ô�����½���-2023.11.17
			if (i["origin_company"], i["prod_names"]) in requir_dict["rummage"]["CustomEdition"].keys():
				report_name = requir_dict["rummage"]["CustomEdition"].get((i["origin_company"], i["prod_names"]))
			else:
				report_name = requir_dict["rummage"]["Universal_simple"].get(i["prod_names"], "") if re.search("����", i["origin_company"]) else requir_dict["rummage"]["Universal_complete"].get(i["prod_names"], "")
	elif i["report_module_type"] == "hospital":
		if (i["company"], i["hosp_depart"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
			report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["hosp_depart"], i["prod_names"]))
		elif (i["company"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
			report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["prod_names"]))
		else:
			report_name = requir_dict["hospital"]["Universal"].get(i["prod_names"], "")
	
	merge_template = report_merge_dict.get(report_name)

	return report_name, merge_template