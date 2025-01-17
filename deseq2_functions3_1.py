import subprocess
import pandas as pd
import numpy as np
import os
import gzip
import seaborn as sns
import shutil
import wget
from tkinter import *
import customtkinter
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
from pydeseq2.default_inference import DefaultInference
from sklearn.preprocessing import StandardScaler
import umap
import umap.umap_ as umap


def process_folders(srr_list, tag, GNAR):
	current_directory = os.getcwd()
#	df1=pd.read_csv(current_directory+"/big.csv", index_col=0)
#	df1=pd.read_csv("E:/CRISPR GEM/Retired_Datasets/data_est_4_29.csv", index_col=0)
	print(srr_list)
	srr_list3 = [srr for srr in srr_list if srr in GNAR.columns]
	print(srr_list3)
	tag_dict = {srr_list[i]: tag[i] for i in range(len(srr_list))}
	new_tag = [tag_dict[item] for item in srr_list3]
	print(new_tag)
	df2=GNAR[srr_list3]
	return df2, new_tag

def get_counts(df, tag, transcript_to_gene):
    df.index=df.index.map(transcript_to_gene)
    df.columns = [col.replace("output_", "") for col in df.columns]
    df2=df.groupby(df.index, axis=0).max()
    df2.fillna(0, inplace=True)
    df2=df2.astype(int)
    counts_per_row = df2.sum(axis=1)
    threshold = len(df2.columns)*2
    df2.columns= tag
    processed_data = df2[counts_per_row > threshold]
    return processed_data

def deseq_get_counts(df):
    df2=df.T
    Condition=df2.index
    metadata = pd.DataFrame({
        "Condition": Condition,
        "Group": Condition})
    metadata.index=df2.index
    counts_df=df2.copy()
    samples_to_keep = ~metadata.Condition.isna()
    counts_df = counts_df.loc[samples_to_keep]
    metadata = metadata.loc[samples_to_keep]
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="Condition",
        refit_cooks=True)
    dds.vst()
    vst_count=dds.layers["vst_counts"]
    fitted=pd.DataFrame(vst_count,index=counts_df.index, columns=counts_df.columns)
    return fitted

import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def umapp(df, tag, fitted):
    print(df.sum())
    print(tag)
    if fitted is None:
        df2 = (df.T)
        df3 = df2[df2.index.notna()]
        df2=df3.fillna(0).astype(int)
        Condition = df2.index
        metadata = pd.DataFrame({
            "Condition": Condition})
        metadata.index = df2.index
        counts_df = df2.copy()
        samples_to_keep = ~metadata.Condition.isna()
        counts_df = counts_df.loc[samples_to_keep]
        metadata = metadata.loc[samples_to_keep]
        genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
        counts_df = counts_df[genes_to_keep]
        print(counts_df)
        unique_columns = []
        for col in tag:
            if col not in unique_columns:
                unique_columns.append(col)
        dds = DeseqDataSet(
            ref_level=['Condition', unique_columns[1]],
            counts=counts_df,
            metadata=metadata,
            design_factors="Condition"
            )
        dds.vst()
        peters_a_bitch = dds.layers["vst_counts"]
        fitted = pd.DataFrame(peters_a_bitch, index=counts_df.index, columns=counts_df.columns)
    
    sns.set(style='white', context='notebook', rc={'figure.figsize': (14, 10)})
    df2 = df.T
    plot_window = customtkinter.CTkToplevel()
    fig, ax = plt.subplots(figsize=(8, 6))
    scaled_fitted = StandardScaler().fit_transform(df2)
    reducer = umap.UMAP(random_state=23, n_neighbors=5)
    embedding = reducer.fit_transform(scaled_fitted)
    if len(set(df2.index)) == 3:
        label_colors = {list(set(df2.index))[0]: 'red', list(set(df2.index))[1]: 'green', list(set(df2.index))[2]: 'blue'}
    elif len(set(df2.index)) == 2:
        label_colors = {list(set(df2.index))[0]: 'red', list(set(df2.index))[1]: 'blue'}
    colors = df2.index.map(label_colors)
    scatter = ax.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=colors,
        alpha=0.5,
        s=100
    )
    plt.gca().set_aspect('equal', 'datalim')
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label)
                   for label, color in label_colors.items()]
    ax.legend(handles=legend_elements, fontsize='large')
    canvas = FigureCanvasTkAgg(fig, master=plot_window)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    return fitted




def graph_gene(df, tag, genename, fitted):
    if fitted is None:
      df2 = (df.T)
      df3 = df2[df2.index.notna()]
      df2=df3.fillna(0).astype(int)
      Condition = df2.index
      metadata = pd.DataFrame({
            "Condition": Condition})
      metadata.index=df2.index
      print('df2 in graph gene')
      print(df2)
      print('Condition')
      print(Condition)
      counts_df=df2.copy()
      samples_to_keep = ~metadata.Condition.isna()
      counts_df = counts_df.loc[samples_to_keep]
      metadata = metadata.loc[samples_to_keep]
      genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
      counts_df = counts_df[genes_to_keep]
      unique_columns = []
      for col in tag:
        if col not in unique_columns:
            unique_columns.append(col)
      dds = DeseqDataSet(
#        ref_level=['Condition',unique_columns[1]],
        counts=counts_df,
        metadata=metadata,
        design_factors="Condition",
        refit_cooks=True)
      dds.vst(use_design=True)
      peters_a_bitch = dds.layers["vst_counts"]
      fitted = pd.DataFrame(peters_a_bitch, index=counts_df.index, columns=counts_df.columns)
    set1_palette = sns.color_palette("Set1")
    start_index = 1
    rearranged_palette = set1_palette[start_index:] + set1_palette[:start_index]
    sns.set_palette(rearranged_palette)
    Data = fitted.copy()
    Data['Labs']=fitted.index
    print('Data')
    print(Data['IL6'])
    print(Data)
    means = Data.groupby(Data['Labs'], sort=False)[genename].mean().reset_index()
    print('means')
    print(means)
    std = Data.groupby(Data['Labs'], sort=False)[genename].std().reset_index()
    error_kw = {'capsize': 5, 'capthick': 1, 'ecolor': 'black'}
    plot_window = customtkinter.CTkToplevel()
  
    fig = plt.Figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    print(means)
    print(std)
    sns.barplot(x=means['Labs'], y=means.iloc[:, 1], data=means, capsize=0.5,
              edgecolor='0.2', lw=2.5, errwidth=2.5, errcolor='0.2', ax=ax)
    kwargs = {'edgecolor': '0.2', 'linewidth': 2.5, 'fc': 'none'}
    sns.swarmplot(data=Data, x=Data['Labs'], y=Data[genename], marker='s', **kwargs, ax=ax)
    ax.set_ylabel(genename + " Expression", fontsize=14, fontweight='bold')
    canvas = FigureCanvasTkAgg(fig, master=plot_window)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    return fitted



def deseq2(df, type, tag, CRISPR):
    df2 = (df.T)
    df3 = df2[df2.index.notna()]
    df2=df3.fillna(0).astype(int)
    Condition = df2.index
    metadata = pd.DataFrame({
            "Condition": Condition})
    metadata.index=df2.index
    counts_df=df2.copy()
    print('counts_df')
    print(counts_df)
    samples_to_keep = ~metadata.Condition.isna()
    counts_df = counts_df.loc[samples_to_keep]
    print('counts_df')
    print(counts_df)
    metadata = metadata.loc[samples_to_keep]
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    print('counts_df')
    print(counts_df)
    counts_df = counts_df[genes_to_keep]
    print('counts_df')
    print(counts_df)
    unique_columns = []
    for col in tag:
        if col not in unique_columns:
            unique_columns.append(col)
    print('tag')
    print(tag)
    print('unique_columns')
    print(unique_columns)
    dds = DeseqDataSet(
        ref_level=['Condition',unique_columns[1]],
        counts=counts_df,
        metadata=metadata,
        design_factors="Condition",
        refit_cooks=True)
    print('matadata')
    print(metadata)
    print('counts_df.index')
    print(counts_df.index)
    print('counts_df.shape')
    print(counts_df.shape)
    dds.deseq2()
    b=dds.varm['LFC']
    c=b.copy()
    c2=c.copy()
    c_sort=c.sort_values(by='Condition_'+unique_columns[0].replace("_","-")+'_vs_'+unique_columns[1].replace("_","-"), ascending=False)
    c_sort.insert(loc=c_sort.columns.get_loc('Condition_'+unique_columns[0].replace("_","-")+'_vs_'+unique_columns[1].replace("_","-")) + 1, column="Rank", value=range(1, len(c) + 1))
    inference = DefaultInference()
    stat_res = DeseqStats(dds, inference=inference)
    stat_res.summary()
    data=stat_res.results_df
    data.insert(loc=data.columns.get_loc("padj") + 1, column="FoldChange", value=c_sort['Condition_'+unique_columns[0].replace("_","-")+'_vs_'+unique_columns[1].replace("_","-")])
    if type == 'output':
        keep3 = data[(data['padj'] < 0.01)]
        if keep3.shape[0] >= 5000:
          keep2=keep3.sort_values(by='padj', ascending=True)
          keep=keep2.head(5000)
        else:
          keep2=keep3.sort_values(by='padj', ascending=True)
          keep1=keep2.head(5000)
          num_r=keep1.shape[0]
          keep4=data[(data['padj'] < 0.05)]
          keep2=keep4.sort_values(by='padj', ascending=True)
          keep5=keep2.head(5000-num_r) 
          keep=pd.concat([keep1, keep5])
        CRISPR_candidates = data[data.index.isin(keep.index)]
        return CRISPR_candidates       
    elif type =='input' and (CRISPR=='CRISPRa' or CRISPR=='CRISPR-ki'):   
        keep = data[((data['log2FoldChange']>0.7)&(data['padj']<0.01))]
        CRISPR_candidates = data[data.index.isin(keep.index)]
        return CRISPR_candidates       
    elif type =='input' and (CRISPR=='CRISPRi' or CRISPR=='CRISPR-ko'):   
        keep = data[((data['log2FoldChange']<1)& (data['padj'] < 0.01))]
        CRISPR_candidates = data[data.index.isin(keep.index)]
        return CRISPR_candidates


def extract_transcript_gene_mapping(gtf_file):
    transcript_to_gene = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'transcript':
                attributes = dict(item.strip().split(' ') for item in fields[8].split(';') if item.strip())
                transcript_id = attributes['transcript_id'].strip('"')
                gene_symbol = attributes['gene_name'].strip('"')
                transcript_to_gene[transcript_id] = gene_symbol
    return transcript_to_gene

def download_annotation():
	url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz"
	gzip_filename = "gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz"
	output_filename = "gencode.v45.chr_patch_hapl_scaff.annotation.gtf"
	if not os.path.exists(output_filename):
	    wget.download(url, gzip_filename)
	    with gzip.open(gzip_filename, 'rb') as f_in:
 	       with open(output_filename, 'wb') as f_out:
 	           shutil.copyfileobj(f_in, f_out)            
	    print("File downloaded and decompressed successfully as", output_filename)
	else:
	    print("File already exists:", output_filename)
	transcript_to_gene = extract_transcript_gene_mapping(output_filename)
	return transcript_to_gene