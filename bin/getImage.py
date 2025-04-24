#-*- coding:gbk -*-
from docxtpl import InlineImage
from docx.shared import Mm
import os
import re
from libs.getImageSize import image_file

'''
Discription
	
	�ýű�������������ͼƬ����Ҫ����ʹ�ò��뷽����InlineImage����ͼƬ

'''
def render_image(tpl, data, jsonDict, report_name, image, config):
	result = {}
	# ����ͼƬ
	# 1. BRCA IGVͼ
	result["brca_all"] = [InlineImage(tpl, i, width=Mm(150)) for i in data["var_brca"]["igvplot"]["all"] if i != None and os.path.exists(i)]
#	print (result["brca_all"])
	result["brca_g"] = [InlineImage(tpl, i, width=Mm(150)) for i in data["var_brca"]["igvplot"]["germline"] if i != None and os.path.exists(i)]
	# 2. MSIͼ
	if data["msi"] and "img_path" in data["msi"].keys() and data["msi"]["img_path"] and os.path.exists(data["msi"]["img_path"]):
		result["msi"] = InlineImage(tpl, data["msi"]["img_path"], width=Mm(90))
	# 3. TMEͼ
	if data["qc"] and "rna_data_qc" in data["qc"].keys() and data["qc"]["rna_data_qc"] and "tmeplot" in data["qc"]["rna_data_qc"].keys() and \
		data["qc"]["rna_data_qc"]["tmeplot"] and os.path.exists(data["qc"]["rna_data_qc"]["tmeplot"]):
		result["tme"] = InlineImage(tpl, data["qc"]["rna_data_qc"]["tmeplot"], width=Mm(110))
	# 4. GEPͼ ���ж�gep.img_path���еĻ�ֱ��չʾ�������ж�qc.rna_data_qc.gepplot
	if data["gep"] and "img_path" in data["gep"].keys() and data["gep"]["img_path"] and os.path.exists(data["gep"]["img_path"]):
		result["gep"] = InlineImage(tpl, data["gep"]["img_path"], width=Mm(90))
	elif data["qc"] and "rna_data_qc" in data["qc"].keys() and data["qc"]["rna_data_qc"] and "gepplot" in data["qc"]["rna_data_qc"].keys() and \
		data["qc"]["rna_data_qc"]["gepplot"] and os.path.exists(data["qc"]["rna_data_qc"]["gepplot"]):
		result["gep"] = InlineImage(tpl, data["qc"]["rna_data_qc"]["gepplot"], width=Mm(90))
	# 5. ����I/II�����Snvindel IGVͼ-2022.10.20
	result["igv_I_II"] = [InlineImage(tpl, i, width=Mm(180)) for i in data["var"]["igv_I_II"] if os.path.exists(i)]
	# �������ø�����һCP40��IGVͼ����11cm����5cm-2023.10.20
	result["igv_I_II_FJFY"] = [InlineImage(tpl, i, width=Mm(110), height=Mm(50)) for i in data["var"]["igv_I_II"] if os.path.exists(i)]
	# 6. ����4/5�����Snvindel IGVͼ-2022.11.18
	result["igv_4_5"] = [InlineImage(tpl, i, width=Mm(180)) for i in data["var"]["igv_4_5"] if os.path.exists(i)]
	# 7. TMBͼ
	if data["tmb"] and "img_path" in data["tmb"].keys() and os.path.exists(data["tmb"]["img_path"]):
		result["tmb"] = InlineImage(tpl, data["tmb"]["img_path"], width=Mm(100))
	# 8. PD-L1ͼ
	if data["pdl1"] and "file_pdl1" in data["pdl1"].keys() and os.path.exists(data["pdl1"]["file_pdl1"]):
		result["pdl1"] = InlineImage(tpl, data["pdl1"]["file_pdl1"], width=Mm(80))
	# 9. CNVͼ
	if "cnv_file_path" in jsonDict.keys() and jsonDict["cnv_file_path"] and "abs_path" in jsonDict["cnv_file_path"].keys() and \
		jsonDict["cnv_file_path"]["abs_path"] and os.path.exists(jsonDict["cnv_file_path"]["abs_path"]):
		result["cnv"] = InlineImage(tpl, jsonDict["cnv_file_path"]["abs_path"], width=Mm(160))
	# 10. MLPAͼ����del��
	result["mlpa_image_del"] = [InlineImage(tpl, i, width=Mm(150)) for i in data["var_brca"]["mlpa_image_del"] if os.path.exists(i)]
	# 11. MRDͼ
	if data['mrd_last'] and 'mrd_img_path' in data['mrd_last'].keys() and os.path.exists(data['mrd_last']['mrd_img_path']):
		result['mrd'] = InlineImage(tpl, data['mrd_last']['mrd_img_path'], width=Mm(180))

	# ��������������չ���Snindel IGVͼ-2023.06.27
	result["igv_onconodrug"] = [InlineImage(tpl, i, width=Mm(180)) for i in data["var"]["igv_onconodrug"] if os.path.exists(i)]
	# ��������-2023.06.27

	# ��ȡ���ñ���ͼƬָ���ߴ磬Ĭ��Ϊ100
	image_sise_dict = image_file(config)
	# ģ���й̶���ͼƬ��ֱ��ƴ�ӻ���ʾ����������߰ѹ̶�ͼƬ���������ʽ�������
	image_list = [i for i in os.listdir(image)]
	for i in image_list:
		image_name = re.split("\.", i)[0]
		width = image_sise_dict.get(image_name, 100) 
		result["fixed_"+image_name] = InlineImage(tpl, image+"/"+i, width=Mm(width))


	# ������滻ͼƬ
	# 7. TMBͼ
#	if data["tmb"] and "img_path" in data["tmb"].keys() and os.path.exists(data["tmb"]["img_path"]):
#		# ����򿪱���TMBͼ���뷽ʽ�޸�-���ټ�ͨ��Master������2022.09.27����������ɽ��
#		if re.search("Master|90|87", data["sample"]["prod_names"]) and re.search("rummage|FDZS|SDSL", report_name):
#			tpl.replace_pic("test_TMB.jpg", data["tmb"]["img_path"])
#			print ("replace")
#		else:
#			result["tmb"] = InlineImage(tpl, data["tmb"]["img_path"], width=Mm(100))
#			print (inline)

	# �滻ͼƬ
	# 8. PD-L1ͼ ͼƬ̫����ܻᱨ��PDL1���滻
#	if data["pdl1"] and "file_pdl1" in data["pdl1"].keys() and os.path.exists(data["pdl1"]["file_pdl1"]):
#		tpl.replace_pic("test_PDL1.jpg", data["pdl1"]["file_pdl1"])
	# 9. CNVͼ
	# ������ͼ-���ԭ������ң�����
	# ʹ�ò���ͼƬ�ķ��������ز����������������д�word�ᱨ�����ܵ�ԭ����1��word��ʽ���⣬2��ģ��汾���죻�Ų�̫��ʱ�ˣ���ʱ����
	# ����replace����������ͼƬ�������ز����ϣ������Ҳ������Ų��ʱ���Ȳ����ˡ���ʱʹ��֮ǰ�����滻ͼƬ��ͼƬģ�壬�ɹ����У������Űɣ��п������Ų�
#	if re.search("Pan116����֯��|LC76����֯��|CRC25����֯��|GA18����֯��|TC21����֯��", data["sample"]["prod_names"]) and re.search("rummage", report_name):
#		if "cnv_file_path" in jsonDict.keys() and jsonDict["cnv_file_path"] and "abs_path" in jsonDict["cnv_file_path"].keys() and jsonDict["cnv_file_path"]["abs_path"] and os.path.exists(jsonDict["cnv_file_path"]["abs_path"]):
#			tpl.replace_pic("test_MSI.png", jsonDict["cnv_file_path"]["abs_path"])
	
	return result
