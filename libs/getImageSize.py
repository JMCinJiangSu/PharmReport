#-*- coding:gbk -*-
import os
#import xlrd
import pandas as pd 

'''
Discription
	
	��ȡͼƬ���ض�չʾ�ߴ� 
	
'''

def image_file(config):
	df = pd.read_excel(os.path.join(config, "report_requirenment.xlsx"), sheet_name="image_size", header=0)
	return df.set_index(df.columns[0]).to_dict()[df.columns[1]]