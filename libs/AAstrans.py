#-*- coding:gbk -*-
import re

'''
	����hgvs_p��ͻ������
	���ر���������
		����ͻ�䵼�»�����뵰�׵�2471λ���������հ���ͻ��Ϊ�鰱�Ტ��2474λ������ǰ��ֹ�� 

'''

def strans(hgvs_p,var_type):
	AA_dict = {"A" : "������",
			   "C" : "���װ���",
			   "D" : "�춬����",
			   "E" : "�Ȱ���",
			   "F" : "��������",
			   "G" : "�ʰ���",
			   "H" : "�鰱��",
			   "I" : "��������",
			   "K" : "������",
			   "L" : "������",
			   "M" : "������",
			   "N" : "�춬����",
			   "P" : "������",
			   "Q" : "�Ȱ�����",
			   "R" : "������",
			   "S" : "˿����",
			   "T" : "�հ���",
			   "V" : "�Ӱ���",
			   "W" : "ɫ����",
			   "Y" : "�Ұ���",
			   "*" : "��ֹ����"}
	outstr = ""
	mat = re.compile(r"\d+")
	#����ͻ�������ͻ��
	if (var_type == "nonSynonymous_Substitution" or var_type == "Nonsense_Mutation") and not re.search("del|ins|dup", hgvs_p):
		pos_list = [i for i in mat.findall(hgvs_p)]
		hgvs_p_stran = [AA_dict[i] for i in list(hgvs_p) if i in AA_dict.keys()]
		outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ��������"+hgvs_p_stran[0]+"ͻ��Ϊ"+hgvs_p_stran[-1]
	#����ͻ��
	if re.search("fs*", hgvs_p):
		for i in mat.findall(re.split("fs\*", hgvs_p)[0]):
			pos = i
		endpos = int(pos) + int(re.split("fs\*", hgvs_p)[-1].replace(")", "")) -1
		hgvs_p_stran = [AA_dict[i] for i in re.split("fs*", hgvs_p)[0] if i in AA_dict.keys()]
		outstr = "��ͻ�䵼�»�����뵰�׵�"+pos+"λ��������"+hgvs_p_stran[0]+"ͻ��Ϊ"+hgvs_p_stran[-1]+"����"+str(endpos)+"λ������ǰ��ֹ"
	#����ͻ��
	if var_type == "Extension":
		pos_list = [i for i in mat.findall(hgvs_p)]
		hgvs_p_stran = [AA_dict[i] for i in re.split("ext", hgvs_p)[0].replace("*", "") if i in AA_dict.keys()]
		outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ��ֹ����ͻ��Ϊ"+hgvs_p_stran[0]+"֮���γ��µ���ֹ����"
	#������ͻ��
	#ȱʧ
	if re.search("del", hgvs_p) and not re.search("delins", hgvs_p):
		pos_list = [i for i in mat.findall(hgvs_p)]
		hgvs_p_stran = [AA_dict[i] for i in hgvs_p if i in AA_dict.keys()]
		if len(pos_list) == 1:
			outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ������ȱʧ"
		else:
			outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ��"+pos_list[-1]+"λ������ȱʧ"
	#����
	if re.search("ins", hgvs_p) and not re.search("delins", hgvs_p):
		pos_list = [i for i in mat.findall(re.split("ins", hgvs_p)[0])]
		insAA = re.split("ins", hgvs_p)[-1].replace(")","")
		if re.search("^[0-9]*$", insAA):
			outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ��"+pos_list[1]+"λ������֮�����"+insAA+"��������"
		else:
			hgvs_p_stran = [AA_dict[i] for i in insAA if i in AA_dict.keys()]
			outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ��"+pos_list[1]+"λ������֮�����"+"��".join(hgvs_p_stran)
	#�ظ�
	if re.search("dup", hgvs_p):
		pos_list = [i for i in mat.findall(hgvs_p)]
		hgvs_p_stran = [AA_dict[i] for i in hgvs_p if i in AA_dict.keys()]
		if len(pos_list) == 1:
			outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ�������ظ�"
		else:
			outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ��"+pos_list[1]+"λ�������ظ�"
	#ȱʧ������
	if re.search("delins", hgvs_p):
		pos_list = [i for i in mat.findall(re.split("delins", hgvs_p)[0])]
		delinsAA = re.split("delins", hgvs_p)[-1].replace(")","")
		if len(pos_list) > 1:
			if re.search("^[0-9]*$", delinsAA):
				outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ��"+pos_list[1]+"λ������ȱʧ������"+delinsAA+"��������"
			else:
				hgvs_p_stran = [AA_dict[i] for i in delinsAA if i in AA_dict.keys()]
				outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ��"+pos_list[1]+"λ������ȱʧ������"+"��".join(hgvs_p_stran)
		else:
			if re.search("^[0-9]*$", delinsAA):
				outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ������ȱʧ������"+delinsAA+"��������"
			else:
				hgvs_p_stran = [AA_dict[i] for i in delinsAA if i in AA_dict.keys()]
				outstr = "��ͻ�䵼�»�����뵰�׵�"+pos_list[0]+"λ������ȱʧ������"+"��".join(hgvs_p_stran)
	#�ں��Ӻͼ���ͻ��
#	if var_type == "Intronic" or var_type == "Splicing":
#		outstr = "��ͻ���������쳣����"

	return outstr
