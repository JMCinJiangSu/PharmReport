#-*- coding:gbk -*-
import copy
from functools import reduce

'''
Discription
	
	�Ľű�������ȡ������Ϣ��������Ϊ���ڱ������ĸ�ʽ��Ŀǰ�еĸ�ʽ���£�
	1. ����������ֱ��ʹ��json���ص���Ϣ�����У�ȫ��չʾ��������صļӴ���ʾ
	2. ��λ��Ϊkey����������CP40��������ʱ��չʾ
	3. ��ҩ��Ϊkey�����������ְ��֣���δ��ӣ�
	4. ������������Ҫ���ְ��֣�ʹ��json������Ϣ��is_same_tumor=1�����ݣ��������ȫ��

'''

def getchemo(jsonDict):
	chemo_data = copy.deepcopy(jsonDict["PGx"])
	chemo_result = {}
	# 1. complete������MP�ټ�ͨ��
	chemo_result["complete"] = sorted(chemo_data, key=lambda i:(i["gene_symbol"], i["dbsnp"]))
	same_tumor_data = [i for i in chemo_result["complete"] if i["is_same_tumor"] == 1]
	# 2. dbsnp_simple_nosplitdrug������CP40-δ���ְ���
	chemo_result["dbsnp_simple_nosplitdrug"] = process_dbsnp_nosplitdrug(chemo_data)
	# 2.2 �����ڱ����������HRR-���ְ��֣�����Ӧ�����޻�����Ϣ�����ȫ��-2022.09.19����
	chemo_result["dbsnp_simple_splittumor"] = process_dbsnp_nosplitdrug(same_tumor_data) if same_tumor_data else process_dbsnp_nosplitdrug(chemo_result["complete"])
	# 3. HRR,����������������Ҫ���ְ��֣�����Ӧ�����޻�����Ϣ�����ȫ��chemo_result["complete"] ��is_same_tumor==1 Ϊͬ���֣�==0Ϊ��ͬ���֡�
	chemo_result["complete_split_tumor"] = reduce_chemo(same_tumor_data) if same_tumor_data else reduce_chemo(chemo_result["complete"])
	# 4. 116���ƣ�ȥ�ع���
	chemo_result["reduce_116"] = process_116(chemo_data)
	# 5. ��ҩ��Ϊ�����ҷְ��֣�����HRR������������-�ɵ�ģ�壩
	chemo_result["drug_splittumor"] = process_drug_splittumor(chemo_data)
	return chemo_result

def process_dbsnp_nosplitdrug(chemo_data):
	# ������CP40
	effict_drug_dict = {}
	for i in chemo_data:
		key = (i["gene_symbol"], i["dbsnp"], i["genotype"], i["evi_level"])
		effict_drug_dict.setdefault(key, [])
		effict_drug_dict[key].append((i["drug_name_cn"], i["impact_type_cn"]))
	dbsnp_list = set([(i["gene_symbol"], i["dbsnp"], i["genotype"], i["evi_level"]) for i in chemo_data])

	result = []
	for info in dbsnp_list:
		tmpdict = {}
		drug_inter = []
		if info in effict_drug_dict.keys():
			for i in effict_drug_dict[info]:
				tmpdict.setdefault(i[1], [])
				tmpdict[i[1]].append(i[0])
		for k, v in tmpdict.items():
			drug_inter.append("��".join(v)+str(k))
		result.append({
			"gene_symbol" : info[0],
			"dbsnp" : info[1],
			"genotype" : info[2],
			"evi_level" : info[3],
			"inter" : drug_inter
		})

	return sorted(result, key = lambda i:(i["gene_symbol"], i["dbsnp"]))

def process_116(chemo_data):
	# 116���Ʋ���Ҫȥ�أ��ø�֪ʶ�ⲻ�ģ���Ҫ����ȥ�أ�����ֻ����д��������������
	# �ֵ仹Ҫ�ع���ԭʼ�����滹��ҩ�����ƣ���ȥ�صĻ����ܻ�����������治׼ȷ�ķ��ա�
	chemo_116 = [{"gene_symbol" : i["gene_symbol"], "dbsnp" : i["dbsnp"], "genotype" : i["genotype"], "clin_anno_cn" : i["clin_anno_cn"], "evi_level" : i["evi_level"]} for i in chemo_data]
	chemo_116 = reduce(lambda x, y:x if y in x else x + [y], [[],]+chemo_116)
	return sorted(chemo_116, key=lambda i:(i["gene_symbol"], i["dbsnp"]))

def process_drug_splittumor(chemo_data):
	# ��ҩ��Ϊ���������ݴ�������HRR�����������ƣ���ͬ���ֵ��޻��ƽ������ŷ�����
	split_tumor_data = [i for i in chemo_data if i["is_same_tumor"] == 1]
	split_tumor_data = split_tumor_data if split_tumor_data else chemo_data
	result = []
	result_tmp = {}
	for i in split_tumor_data:
		result_tmp.setdefault(i["drug_name_cn"], [])
		result_tmp[i["drug_name_cn"]].append(i)
	
	for k,v in result_tmp.items():
		v_sort = sorted(v, key=lambda i:(i["evi_level"], i["gene_symbol"]))
		result.append({
			"drug_name_cn" : k,
			"info" : v_sort
		})
	result = sorted(result, key=lambda i:i["drug_name_cn"])
	return result

# ȥ�أ�����HRR
def reduce_chemo(chemo_data):
	result = []
	chemo_list = [{"gene_symbol" : i["gene_symbol"], "dbsnp" : i["dbsnp"], "genotype" : i["genotype"], "clin_anno_cn" : i["clin_anno_cn"], "evi_level" : i["evi_level"]} for i in chemo_data]
	chemo_list = reduce(lambda x, y:x if y in x else x + [y], [[],]+chemo_list)
	return sorted(chemo_list, key=lambda i:(i["gene_symbol"], i["dbsnp"]))