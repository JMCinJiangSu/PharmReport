#-*- coding:gbk -*-

'''
Discription
	
	��MSI��PD-L1��TMB��Ӧ��Ϊ�ֵ�ȴ���Ϊ�б����Ϣת���ֵ䣬���ж������Ĭ��ѡ���б��һ��Ԫ�� 
	
'''

def ListToDict(json_result):
	return json_result if type(json_result).__name__=="dict" else json_result[0] if json_result else {}