#-*- coding:gbk -*-
#from pypinyin import lazy_pinyin
from itertools import chain
from functools import lru_cache
from pypinyin import pinyin, Style

'''
Discription
	
	������ת��Ϊƴ������������ 
	
'''

@lru_cache(maxsize=128)
def topinyin(instr):
	# ʹ��pinyin��������꣬�硰���ġ����ء�zhong1wen2��������lazy_pinyin�����ء�zhongwen��
	# 2023.06.07-���ǸĻش�����ģ������Ļ����ܻ��С��������ᡱ���ڡ��������ᡱǰ�������������������ᡱ���ڰ�������ǰ���е����
	return "".join(chain.from_iterable(pinyin(instr, Style.TONE3)))
	#return "".join(chain.from_iterable(lazy_pinyin(instr)))
