import subprocess
import pandas as pd
import numpy as np
import os
import gzip
import shutil
import wget
import tkinter as tk
from tkinter import ttk
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
from pydeseq2.default_inference import DefaultInference
from deseq2_functions import process_folders
from deseq2_functions import get_counts
from deseq2_functions import extract_transcript_gene_mapping
from deseq2_functions import deseq2
from deseq2_functions import graph_gene
from display_dataset import display_dataset
from tkinter import *
import customtkinter
from deseq2_functions import umapp
import random

def run_deseq2(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, CRISPR, window):
	global output_genes, input_genes, saveas_

	platform_counts_exp = complete_exp['Platform'].value_counts()
	platform_counts_end = complete_end['Platform'].value_counts()
	rank_exp = platform_counts_exp.rank(method='max')
	rank_end = platform_counts_end.rank(method='max')
	ranks_df = pd.DataFrame({'exp_rank': rank_exp, 'end_rank': rank_end})
	ranks_df['average_rank'] = ranks_df.mean(axis=1)
	ranks_df = ranks_df.dropna()
	if ranks_df.empty:
	  warning = customtkinter.CTkLabel(master=window, text="*COMPARING ACROSS PLATFORMS*", fg_color=("red"))
	  warning.grid(row=0, column=4, padx=10)
	  end_set=complete_end.copy()
	  end_labs = list(set(end_set.index.tolist()))
	  end_tag= [small_end+"_"+big_end]*len(end_labs)
	  exp_set=complete_exp.copy()
	  exp_labs = list(set(exp_set.index.tolist()))
	  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)	  
	else:
	  most_common_platform = ranks_df['average_rank'].idxmax()
#	  most_common_platform='HiSeq2000'
	  end_set=complete_end[(complete_end['Platform']==most_common_platform)].copy()
	  end_labs = list(set(end_set.index.tolist()))
	  end_tag= [small_end+"_"+big_end]*len(end_labs)
	  exp_set=complete_exp[(complete_exp['Platform']==most_common_platform)].copy()
	  exp_labs = list(set(exp_set.index.tolist()))
	  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)

	if CRISPR == "Select CRISPR Method" or not CRISPR:
	  lbl_space2 = customtkinter.CTkLabel(master=window, text="*Please Select CRISPR Type*", fg_color=("red"))
	  lbl_space2.grid(row=2, column=0, columnspan=2)
	  return
	min_len=100
	if man_input_list != "NA":
		end_len=len(end_labs)
		exp_len=len(exp_labs)
		if end_len < min_len or exp_len < min_len:
		  min_len=min(end_len, exp_len)
		end_labs = random.sample(end_labs, min_len)
		exp_labs = random.sample(exp_labs, min_len)
		outputs= end_labs+exp_labs
		output_tag=[end_tag[0]]*min_len+[exp_tag[0]]*min_len
		input_genes=man_input_list 
	elif small_mid !='NA':
		platform_counts_exp = complete_exp['Platform'].value_counts()
		platform_counts_end = complete_end['Platform'].value_counts()
		platform_counts_mid = complete_mid['Platform'].value_counts()
		rank_exp = platform_counts_exp.rank(method='max')
		rank_end = platform_counts_end.rank(method='max')
		rank_mid = platform_counts_mid.rank(method='max')
		ranks_df = pd.DataFrame({'exp_rank': rank_exp, 'end_rank': rank_end, 'mid_rank': rank_mid})
		ranks_df['average_rank'] = ranks_df.mean(axis=1)
		ranks_df = ranks_df.dropna()
		if ranks_df.empty:
		  print('Data is empty')
		  warning = customtkinter.CTkLabel(master=window, text="*COMPARING ACROSS PLATFORMS*", fg_color=("red"))
		  warning.grid(row=0, column=3, padx=10)
		  end_set=complete_end.copy()
		  end_labs = list(set(end_set.index.tolist()))
		  end_tag= [small_end+"_"+big_end]*len(end_labs)
		  exp_set=complete_exp.copy()
		  exp_labs = list(set(exp_set.index.tolist()))
		  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
		  mid_set=complete_mid.copy()
		  mid_labs = list(set(mid_set.index.tolist()))
		  mid_tag= [small_mid+"_"+big_mid]*len(mid_labs)	  
		else:
		  most_common_platform = ranks_df['average_rank'].idxmax()
		  end_set=complete_end[(complete_end['Platform']==most_common_platform)].copy()
		  end_labs = list(set(end_set.index.tolist()))
		  end_tag= [small_end+"_"+big_end]*len(end_labs)
		  exp_set=complete_exp[(complete_exp['Platform']==most_common_platform)].copy()
		  exp_labs = list(set(exp_set.index.tolist()))
		  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
		  mid_set=complete_mid[(complete_mid['Platform']==most_common_platform)].copy()
		  mid_labs = list(set(mid_set.index.tolist()))
		  mid_tag= [small_mid+"_"+big_mid]*len(mid_labs)
		  end_len=len(end_labs)
		  exp_len=len(exp_labs)
		  mid_len=len(mid_labs)
		  if end_len < min_len or exp_len < min_len or mid_len < min_len:
		    min_len=min(end_len, exp_len, mid_len)
		  end_labs = random.sample(end_labs, min_len)
		  end_tag= [small_end+"_"+big_end]*len(end_labs)
		  exp_labs = random.sample(exp_labs, min_len)
		  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
		  mid_labs = random.sample(mid_labs, min_len)
		  mid_tag= [small_mid+"_"+big_mid]*len(mid_labs)
		outputs= end_labs+exp_labs 
		output_tag=end_tag+exp_tag
		inputs=mid_labs+exp_labs
		input_tag=mid_tag+exp_tag
		input_int, input_tag=process_folders(inputs, input_tag)
		input_int2=get_counts(input_int, input_tag, transcript_to_gene)
		if CRISPR != 'CRISPRa/i':
		    input_genes=deseq2(input_int2, "input", input_tag, CRISPR)
		else:
		    input_genes, CRISPRa_part, CRISPRi_part=deseq2(input_int2, "input", input_tag, CRISPR)
	else:
		end_len=len(end_labs)
		exp_len=len(exp_labs)
		if end_len < min_len or exp_len < min_len:
		  min_len=min(end_len, exp_len)
		end_labs = random.sample(end_labs, min_len)
		exp_labs = random.sample(exp_labs, min_len)
		inputs=end_labs+exp_labs
		input_tag = [end_tag[0]] * min_len + [exp_tag[0]] * min_len
		outputs= end_labs+exp_labs
		output_tag = [end_tag[0]] * min_len + [exp_tag[0]] * min_len
		input_int, input_tag=process_folders(inputs, input_tag)
		input_int2=get_counts(input_int, input_tag, transcript_to_gene)
		if CRISPR != 'CRISPRa/i':
		  input_genes=deseq2(input_int2, "input", input_tag, CRISPR)
		else:
		  input_genes, CRISPRa_part, CRISPRi_part=deseq2(input_int2, "input", input_tag, CRISPR)

	
	output_int, output_tag=process_folders(outputs, output_tag)
	print(output_int)
	print(output_tag)
	output_int2=get_counts(output_int, output_tag, transcript_to_gene)
	output_genes=deseq2(output_int2, "output", output_tag, CRISPR)
	style = ttk.Style()
	sty = ttk.Style(window)
	sty.configure("TRadiobutton", background="navy",
              foreground="deep sky blue", anchor="w",
              justify='left')
	current_directory = os.getcwd()
	input_csv_path = current_directory + '/test_input.csv'
	if os.path.exists(input_csv_path):
	   os.remove(input_csv_path)
	output_csv_path = current_directory + '/test_output.csv'
	if os.path.exists(output_csv_path):
	   os.remove(output_csv_path)
	input_genes.to_csv(input_csv_path)
	output_genes.to_csv(output_csv_path)

	def saver(input_genes, output_genes, current_directory, end_labs, end_tag, exp_labs, exp_tag):
              saveas_ = saveas.get()
              os.makedirs(current_directory+"/"+saveas_, exist_ok=True)
              save_path=current_directory+"/"+saveas_+"/"+saveas_
              input_genes.to_csv(save_path+ '_input_genes.csv')
              output_genes.to_csv(save_path+ '_output_genes.csv')
              output_genes2 = pd.DataFrame(output_genes.index, columns=['outputs'])
              output_genes2.to_csv(save_path + '_output_list.csv')
              input_genes2 = pd.DataFrame(input_genes.index, columns=['inputs'])
              input_genes2.to_csv(save_path + '_input_list.csv')
              data = {
                end_tag[0]: end_labs,
                exp_tag[0]: exp_labs
              }
              df = pd.DataFrame(data)
              df.to_csv(save_path+'_SRRs_used.csv')
              saveas_ = ', '.join(saveas_)
              savior=current_directory+"/saveas_.txt"
              with open(savior, "w") as file:
                   file.write(saveas_)
              return saveas_
	def saver2(input_genes, output_genes, CRISPRa_part, CRISPRi_part, current_directory):
              saveas_ = saveas.get()
              os.makedirs(current_directory+"/"+saveas_, exist_ok=True)
              save_path=current_directory+"/"+saveas_+"/"+saveas_
              input_genes.to_csv(save_path+ '_input_genes.csv')
              output_genes.to_csv(save_path+ '_output_genes.csv')
              output_genes2 = pd.DataFrame(output_genes.index, columns=['outputs'])
              output_genes2.to_csv(save_path + '_output_list.csv')
              CRISPRa_genes=pd.DataFrame(CRISPRa_part.index, columns=['CRISPRa'])
              CRISPRa_genes.to_csv(save_path+ '_a_list.csv')
              CRISPRi_genes=pd.DataFrame(CRISPRi_part.index, columns=['CRISPRi'])
              CRISPRi_genes.to_csv(save_path+ '_i_list.csv')
              input_genes2 = pd.DataFrame(input_genes.index, columns=['inputs'])
              input_genes2.to_csv(save_path + '_input_list.csv')
              saveas_ = ', '.join(saveas_)
              savior=current_directory+"/saveas_.txt"
              with open(savior, "w") as file:
                   file.write(saveas_)
              return saveas_
	new_widg=[]
	lbl_save = customtkinter.CTkLabel(master=window, text="File Prefix: ")
	lbl_save.grid(row=4, column=4,pady=5)
	new_widg.append(lbl_save)
	saveas = customtkinter.CTkEntry(master=window, width=150)
	saveas.grid(row=5, column=4, pady=5, sticky="e", padx=5)
	new_widg.append(saveas)
	if CRISPR !='CRISPRa/i':
	  btn_save = customtkinter.CTkButton(window, text="Save", command=lambda: saver(input_genes, output_genes, current_directory, end_labs, end_tag, exp_labs, exp_tag))
	  btn_save.grid(row=6, column=4, padx=5, pady=5)
	  new_widg.append(btn_save)
	else:
	  btn_save = customtkinter.CTkButton(window, text="Save", command=lambda: saver2(input_genes, output_genes, CRISPRa_part, CRISPRi_part, current_directory))
	  btn_save.grid(row=11, column=5, padx=5, pady=5)
	  new_widg.append(btn_save)
	btn_display_dataset1 = customtkinter.CTkButton(window, text="Display Input Genes", command=lambda: display_dataset('inputs', window))
	btn_display_dataset1.grid(row=2, column=4, padx=5, pady=5)
	new_widg.append(btn_display_dataset1)
	btn_display_dataset2 = customtkinter.CTkButton(window, text="Display Output Genes", command=lambda: display_dataset('outputs', window))
	btn_display_dataset2.grid(row=3, column=4, padx=5, pady=5)
	new_widg.append(btn_display_dataset2)
	return output_genes, input_genes, new_widg


def run_graph(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, genename, fitted, df1, processed_data, window):
	global output_genes, input_genes, df2, processed_data1
	platform_counts_exp = complete_exp['Platform'].value_counts()
	platform_counts_end = complete_end['Platform'].value_counts()
	rank_exp = platform_counts_exp.rank(method='max')
	rank_end = platform_counts_end.rank(method='max')
	ranks_df = pd.DataFrame({'exp_rank': rank_exp, 'end_rank': rank_end})
	ranks_df['average_rank'] = ranks_df.mean(axis=1)
	ranks_df = ranks_df.dropna()
	if ranks_df.empty:
	  print('yoooooooooooo its empty')
	  warning = customtkinter.CTkLabel(master=window, text="*COMPARING ACROSS PLATFORMS*", fg_color=("red"))
	  warning.grid(row=0, column=3, padx=10)
	  end_set=complete_end.copy()
	  end_labs = list(set(end_set.index.tolist()))
	  end_tag= [small_end+"_"+big_end]*len(end_labs)
	  exp_set=complete_exp.copy()
	  exp_labs = list(set(exp_set.index.tolist()))
	  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)	  
	else:
	  most_common_platform = ranks_df['average_rank'].idxmax()
	  end_set=complete_end[(complete_end['Platform']==most_common_platform)].copy()
	  end_labs = list(set(end_set.index.tolist()))
	  end_tag= [small_end+"_"+big_end]*len(end_labs)
	  exp_set=complete_exp[(complete_exp['Platform']==most_common_platform)].copy()
	  exp_labs = list(set(exp_set.index.tolist()))
	  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)

	genename=genename.get()
	new_widg1=[]
	min_len=100
	if man_input_list != "NA":
		end_len=len(end_labs)
		exp_len=len(exp_labs)
		if end_len < min_len or exp_len < min_len:
		  min_len=min(end_len, exp_len)
		end_labs = random.sample(end_labs, min_len)
		exp_labs = random.sample(exp_labs, min_len)
		inputs=end_labs+exp_labs
		input_tag=[end_tag[0]]*min_len+[exp_tag[0]]*min_len
		input_int=df1
		input_int2=processed_data
		if input_int is None or input_int2 is None:
		  input_int, input_tag=process_folders(inputs, input_tag)
		  input_int2=get_counts(input_int, input_tag, transcript_to_gene)
		a=graph_gene(input_int2, input_tag, genename, fitted)
	elif small_mid !='NA':
		platform_counts_exp = complete_exp['Platform'].value_counts()
		platform_counts_end = complete_end['Platform'].value_counts()
		platform_counts_mid = complete_mid['Platform'].value_counts()
		rank_exp = platform_counts_exp.rank(method='max')
		rank_end = platform_counts_end.rank(method='max')
		rank_mid = platform_counts_mid.rank(method='max')
		ranks_df = pd.DataFrame({'exp_rank': rank_exp, 'end_rank': rank_end, 'mid_rank': rank_mid})
		ranks_df['average_rank'] = ranks_df.mean(axis=1)
		ranks_df = ranks_df.dropna()
		if ranks_df.empty:
		  warning = customtkinter.CTkLabel(master=window, text="*COMPARING ACROSS PLATFORMS*", fg_color=("red"))
		  warning.grid(row=0, column=3, padx=10)
		  end_set=complete_end.copy()
		  end_labs = list(set(end_set.index.tolist()))
		  end_tag= [small_end+"_"+big_end]*len(end_labs)
		  exp_set=complete_exp.copy()
		  exp_labs = list(set(exp_set.index.tolist()))
		  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
		  mid_set=complete_mid.copy()
		  mid_labs = list(set(mid_set.index.tolist()))
		  mid_tag= [small_mid+"_"+big_mid]*len(mid_labs)	  
		else:
		  most_common_platform = ranks_df['average_rank'].idxmax()
		  end_set=complete_end[(complete_end['Platform']==most_common_platform)].copy()
		  end_labs = list(set(end_set.index.tolist()))
		  end_tag= [small_end+"_"+big_end]*len(end_labs)
		  exp_set=complete_exp[(complete_exp['Platform']==most_common_platform)].copy()
		  exp_labs = list(set(exp_set.index.tolist()))
		  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
		  mid_set=complete_mid[(complete_mid['Platform']==most_common_platform)].copy()
		  mid_labs = list(set(mid_set.index.tolist()))
		  mid_tag= [small_mid+"_"+big_mid]*len(mid_labs)
		  end_len=len(end_labs)
		  exp_len=len(exp_labs)
		  mid_len=len(mid_labs)
		  if end_len < min_len or exp_len < min_len or mid_len < min_len:
		    min_len=min(end_len, exp_len, mid_len)
		  end_labs = random.sample(end_labs, min_len)
		  end_tag= [small_end+"_"+big_end]*len(end_labs)
		  exp_labs = random.sample(exp_labs, min_len)
		  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
		  mid_labs = random.sample(mid_labs, min_len)
		  mid_tag= [small_mid+"_"+big_mid]*len(mid_labs)
		inputs=exp_labs+mid_labs+end_labs
		input_tag=exp_tag+mid_tag+end_tag
		input_int=df1
		input_int2=processed_data
		if input_int is None or input_int2 is None:
		  input_int, input_tag=process_folders(inputs, input_tag)
		  input_int2=get_counts(input_int, input_tag, transcript_to_gene)
		a=graph_gene(input_int2, input_tag, genename, fitted)
	else:
		end_len=len(end_labs)
		exp_len=len(exp_labs)
		if end_len < min_len or exp_len < min_len:
		  min_len=min(end_len, exp_len)
		end_labs = random.sample(end_labs, min_len)
		exp_labs = random.sample(exp_labs, min_len)
		inputs=end_labs+exp_labs
		input_tag = [end_tag[0]] * min_len + [exp_tag[0]] * min_len
		print(input_tag)
		input_int=df1
		input_int2=processed_data
		if input_int is None or input_int2 is None:
		  input_int, input_tag=process_folders(inputs, input_tag)
		  input_int2=get_counts(input_int, input_tag, transcript_to_gene)
		a=graph_gene(input_int2, input_tag, genename, fitted)
	fitted1=a
	df2=input_int
	processed_data1=input_int2
	return fitted1, df2, processed_data1, new_widg1


#def run_umap(value, exp_labs, exp_tag, end_labs, end_tag, mid_labs, mid_tag, man_input_list, transcript_to_gene, fitted, df1, processed_data, window):
def run_umap(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, fitted, df1, processed_data, window):
	global output_genes, input_genes, df2, processed_data1
	platform_counts_exp = complete_exp['Platform'].value_counts()
	platform_counts_end = complete_end['Platform'].value_counts()
	rank_exp = platform_counts_exp.rank(method='max')
	rank_end = platform_counts_end.rank(method='max')
	ranks_df = pd.DataFrame({'exp_rank': rank_exp, 'end_rank': rank_end})
	ranks_df['average_rank'] = ranks_df.mean(axis=1)
	ranks_df = ranks_df.dropna()
	if ranks_df.empty:
	  print('Data empty')
	  end_set=complete_end.copy()
	  end_labs = list(set(end_set.index.tolist()))
	  end_tag= [small_end+"_"+big_end]*len(end_labs)
	  exp_set=complete_exp.copy()
	  exp_labs = list(set(exp_set.index.tolist()))
	  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)	  
	else:
	  most_common_platform = ranks_df['average_rank'].idxmax()  
	  end_set=complete_end[(complete_end['Platform']==most_common_platform)].copy()
	  end_labs = list(set(end_set.index.tolist()))
	  end_tag= [small_end+"_"+big_end]*len(end_labs)
	  exp_set=complete_exp[(complete_exp['Platform']==most_common_platform)].copy()
	  exp_labs = list(set(exp_set.index.tolist()))
	  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
	min_len=100
	if man_input_list != "NA":
		end_len=len(end_labs)
		exp_len=len(exp_labs)
		if end_len < min_len or exp_len < min_len:
		  min_len=min(end_len, exp_len)
		end_labs = random.sample(end_labs, min_len)
		exp_labs = random.sample(exp_labs, min_len)
		inputs=end_labs+exp_labs
		input_tag=[end_tag[0]]*min_len+[exp_tag[0]]*min_len
		input_int=df1
		input_int2=processed_data
		if input_int is None or input_int2 is None:
		  input_int, input_tag=process_folders(inputs, input_tag)
		  input_int2=get_counts(input_int, input_tag, transcript_to_gene)
		a=umapp(input_int2, input_tag, fitted)
	elif small_mid !='NA':
		platform_counts_exp = complete_exp['Platform'].value_counts()
		platform_counts_end = complete_end['Platform'].value_counts()
		platform_counts_mid = complete_mid['Platform'].value_counts()
		rank_exp = platform_counts_exp.rank(method='max')
		rank_end = platform_counts_end.rank(method='max')
		rank_mid = platform_counts_mid.rank(method='max')
		ranks_df = pd.DataFrame({'exp_rank': rank_exp, 'end_rank': rank_end, 'mid_rank': rank_mid})
		ranks_df['average_rank'] = ranks_df.mean(axis=1)
		ranks_df = ranks_df.dropna()
		if ranks_df.empty:
		  print('yoooooooooooo its empty')
		  warning = customtkinter.CTkLabel(master=window, text="*COMPARING ACROSS PLATFORMS*", fg_color=("red"))
		  warning.grid(row=0, column=3, padx=10)
		  end_set=complete_end.copy()
		  end_labs = list(set(end_set.index.tolist()))
		  end_tag= [small_end+"_"+big_end]*len(end_labs)
		  exp_set=complete_exp.copy()
		  exp_labs = list(set(exp_set.index.tolist()))
		  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
		  mid_set=complete_mid.copy()
		  mid_labs = list(set(mid_set.index.tolist()))
		  mid_tag= [small_mid+"_"+big_mid]*len(mid_labs)	  
		else:
		  most_common_platform = ranks_df['average_rank'].idxmax()
		  end_set=complete_end[(complete_end['Platform']==most_common_platform)].copy()
		  end_labs = list(set(end_set.index.tolist()))
		  end_tag= [small_end+"_"+big_end]*len(end_labs)
		  exp_set=complete_exp[(complete_exp['Platform']==most_common_platform)].copy()
		  exp_labs = list(set(exp_set.index.tolist()))
		  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
		  mid_set=complete_mid[(complete_mid['Platform']==most_common_platform)].copy()
		  mid_labs = list(set(mid_set.index.tolist()))
		  mid_tag= [small_mid+"_"+big_mid]*len(mid_labs)
		  end_len=len(end_labs)
		  exp_len=len(exp_labs)
		  mid_len=len(mid_labs)
		  if end_len < min_len or exp_len < min_len or mid_len < min_len:
		    min_len=min(end_len, exp_len, mid_len)
		  end_labs = random.sample(end_labs, min_len)
		  end_tag= [small_end+"_"+big_end]*len(end_labs)
		  exp_labs = random.sample(exp_labs, min_len)
		  exp_tag= [small_exp+"_"+big_exp]*len(exp_labs)
		  mid_labs = random.sample(mid_labs, min_len)
		  mid_tag= [small_mid+"_"+big_mid]*len(mid_labs)
		inputs=exp_labs+mid_labs+end_labs
		input_tag=exp_tag+mid_tag+end_tag
		input_int=df1
		input_int2=processed_data
		if input_int is None or input_int2 is None:
		  input_int, input_tag=process_folders(inputs, input_tag)
		  input_int2=get_counts(input_int, input_tag, transcript_to_gene)
		a=umapp(input_int2, input_tag, fitted)
	else:
		end_len=len(end_labs)
		exp_len=len(exp_labs)
		if end_len < min_len or exp_len < min_len:
		  min_len=min(end_len, exp_len)
		end_labs = random.sample(end_labs, min_len)
		print(['exp_labs']+exp_labs)
		print(end_labs)
		exp_labs = random.sample(exp_labs, min_len)
		inputs=end_labs+exp_labs
		print(['end_tag']+end_tag)
		print(['exp_tag']+exp_tag)
		input_tag = [end_tag[0]] * min_len + [exp_tag[0]] * min_len
		print(input_tag)
		input_int=df1
		input_int2=processed_data
		if input_int is None or input_int2 is None:
		  input_int, input_tag=process_folders(inputs, input_tag)
		  print(input_tag)
		  input_int2=get_counts(input_int, input_tag, transcript_to_gene)
		a=umapp(input_int2, input_tag, fitted)
	fitted1=a
	df2=input_int
	processed_data1=input_int2
	return fitted1, df2, processed_data1

