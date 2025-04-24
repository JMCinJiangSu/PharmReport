#-*- coding:gbk -*-
import copy
from libs import listResultToDict

'''
Discription
	
	�ýű�������ȡRNA�������������ģ��������з��д��� 

'''

def getRNA_exp(jsonDict):
	rna_exp_data = copy.deepcopy(jsonDict["rna_exp"])

	# ����һ������-����v4-2023.11.06
	# tpmΪ��ʱ��v3���ء�0.0����v4���ر���˿�ֵ��������˵��Ҫ�ı���ű�
	for i in rna_exp_data:
		if not i["tpm"]:
			i["tpm"] = "0.0"
	# �������-2023.11.06
	
	rna_exp = {}
	if rna_exp_data:
		# �ټ�ͨ�ð�չʾ��������
		rna_exp["column_4"] = process_data(rna_exp_data, 4) 
		# ����չʾ����
		rna_exp["column_5"] = process_data(rna_exp_data, 5)
		# ������ҪչʾGEP��ػ���ı��
		rna_exp["gep"] = gep_tpm(rna_exp_data)
		# xw5101 չʾ25�������tpm
		rna_exp['xw5101'] = process_data(exp_xw5101(rna_exp_data), 2)

	return rna_exp

# numָ����ģ����Ҫ������չʾ������+TPMΪһ��
def process_data(rna_exp_data, num):
	rna_exp = []
	for i in range(0, len(rna_exp_data)-num, num):
		tmp_dict = {}
		for j in range(1, num + 1):
			tmp_dict["gene"+str(j)] = rna_exp_data[i+j-1]["gene_symbol"]
			tmp_dict["tpm"+str(j)] = rna_exp_data[i+j-1]["tpm"]
		rna_exp.append(tmp_dict)
	
	rest_gene = len(rna_exp_data) % num
	rest_tmp_dict = {}
	for j in range(1, num+1):
		rest_tmp_dict["gene"+str(j)] = ""
		rest_tmp_dict["tpm"+str(j)] = ""

	rest_num = 1
	last_row_num = len(rna_exp_data)-rest_gene if rest_gene != 0 else len(rna_exp_data)-rest_gene-num
	for j in range(last_row_num, len(rna_exp_data)):
		rest_tmp_dict["gene"+str(rest_num)] = rna_exp_data[j]["gene_symbol"]
		rest_tmp_dict["tpm"+str(rest_num)] = rna_exp_data[j]["tpm"]
		rest_num += 1
	rna_exp.append(rest_tmp_dict)

	return rna_exp

# GEP ��ػ���ı����-��������Master
def gep_tpm(rna_exp_data):
	gene_list = ["CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6","HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT"]
	gep_exp = {}
	for i in rna_exp_data:
		if i["gene_symbol"] in gene_list:
			gep_exp.setdefault(i["gene_symbol"].replace("-",""), "")
			gep_exp[i["gene_symbol"].replace("-","")] = i["tpm"]

	return gep_exp

# XW5101 v3 RNAseqֻ��25������
def exp_xw5101(rna_exp_data):
	gene_list = ['BCL2', 'BCL2L1', 'CD14', 'CD34', 'CD38', 'CDK2', 'CDK4', 'CDK6', 'CDK7', 'CDK9', 'FLT3', 'HLA-A', 'HLA-B', 'HLA-C', 'HOXA9', 'ITGAM', 'MCL1', 'MEIS1', 'MEN1', 'MNDA', 'MYC', 'PBX3', 'PRMT5', 'PTPRC', 'TP53']
	exp_xw5101 = []
	for i in rna_exp_data:
		if i['gene_symbol'] in gene_list:
			exp_xw5101.append(i)
	
	return exp_xw5101