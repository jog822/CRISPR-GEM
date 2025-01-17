import subprocess
import importlib.util
def install_if_needed(package_name):
    if importlib.util.find_spec(package_name) is None:
        subprocess.run(["pip", "install", package_name, "--quiet"], check=True)
    else:
        print(f"{package_name} is already installed.")

packages = [
    "pydeseq2",
    "wget",
    "pandas",
    "seaborn",
    "customtkinter",
    "umap",
    "umap-learn"
]

for package in packages:
    install_if_needed(package)
import tkinter as tk
import pandas as pd
import numpy as np
import os
import gzip
import shutil
from tkinter import ttk
import wget
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
from pydeseq2.default_inference import DefaultInference
from deseq2_functions import extract_transcript_gene_mapping
from deseq2_functions import download_annotation
from run_deseq3_1 import run_deseq2
from run_deseq3_1 import run_graph
from run_deseq3_1 import run_umap
from deseq2_functions3_1 import graph_gene
from deseq2_functions3_1 import get_counts
from deseq2_functions3_1 import umapp
from display_dataset import display_dataset
import tkinter.font as font
from PIL import Image, ImageTk

transcript_to_gene=download_annotation()
added_widgets =[]
CRISPR=[]

def on_select1(value):
    global CRISPR
    CRISPR=value
    print("Selected:", value)

#def update_srr_input(srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid):
#    # Update labels based on entered data
#    srr_input = manual_widg.get()  # Read from the Entry widget
#    srr_input_list = srr_input.split(",")
#    srr_input_list = [item.strip() for item in srr_input_list]
#    
#    entry_lab_end.config(text=f"Entered SRR IDs: {', '.join(srr_input_list)}")
#    
#    # Update the other fields accordingly
#    # Similarly update the other fields like labels or widgets when needed

def input_type(value1, window):
    global value
    value = value1
    current_directory = os.getcwd()
    labs = pd.read_csv(current_directory + '/Labels.csv', index_col=0)
    labs = labs.astype(str)
    labs_list = list(set(labs['Label'].tolist()))
      
    def end_select(value3, labb, labs, x, y, value):
        labs2 = labs[labs['Label'] == value3]
        lab_list2 = list(set(labs2['Description']))
        lab_list3 = lab_list2.copy()
        lab_list3.insert(0, 'Select ID')
        def end_select2(value2, labb, labs, value, value3):
            end_labs_whole = labs[(labs['Label'] == value3) & (labs['Description'] == value2)].copy()
            globals() ['complete_'+labb] = end_labs_whole
            globals() ['big_'+labb] = value3
            globals() ['small_'+labb] = value2
            return 
        selected_option_dos = customtkinter.StringVar(value=lab_list3[0])
        dropdown_end2 = customtkinter.CTkOptionMenu(master=window, values=lab_list3,
                                                    command=lambda value2: end_select2(value2, labb, labs, value, value3),
                                                    variable=selected_option_dos)
        new_x = int(x + 1)
        dropdown_end2.grid(row=y, column=new_x, padx=5, pady=10)
        added_widgets.append(dropdown_end2)
        return

    global man_input_list, mid_labs, mid_tag, big_exp   
    if value == "End-to-End":
        labs_list_exp = labs_list.copy()
        labs_list_exp.insert(0, 'Select Experimental Group')
        x = int(2)
        y = int(2)
        entry_lab_exp = customtkinter.CTkLabel(master=window, text="Experimental Group: ")
        entry_lab_exp.grid(row=y - 1, column=x, padx=10)
        added_widgets.append(entry_lab_exp)
        selected_option_exp = customtkinter.StringVar(value=labs_list_exp[0])
        dropdown_exp = customtkinter.CTkOptionMenu(master=window, values=labs_list_exp,
                                                   command=lambda value3: end_select(value3, 'exp', labs, x, y, value),
                                                   variable=selected_option_exp)
        dropdown_exp.grid(row=y, column=x, padx=5, pady=10)
        added_widgets.append(dropdown_exp)
        labs_list_end = labs_list.copy()
        labs_list_end.insert(0, 'Select Experimental Group')
        x_end=int(2)
        y_end=int(4)
        entry_lab_end = customtkinter.CTkLabel(master=window, text="Target Group: ")
        entry_lab_end.grid(row=y_end - 1, column=x_end, padx=10)
        selected_option_end = customtkinter.StringVar(value=labs_list_end[0])
        dropdown_end = customtkinter.CTkOptionMenu(master=window, values=labs_list_end,
                                                   command=lambda value4: end_select(value4, 'end', labs, x_end, y_end, value),
                                                   variable=selected_option_end)
        dropdown_end.grid(row=y_end, column=x_end, padx=5, pady=10)
        added_widgets.append(dropdown_end)
        added_widgets.append(entry_lab_end)

        globals() ['man_input_list'] = 'NA' 
        globals() ['big_mid'] = 'NA' 
        globals() ['small_mid'] = 'NA' 
        globals() ['complete_mid'] = 'NA' 
        globals() ['srr_lab_exp'] = 'NA'
        globals() ['srr_input_list'] = 'NA' 
        globals() ['srr_lab_end'] = 'NA' 
        globals() ['srr_input_list_exp'] = 'NA' 
        globals() ['srr_lab_mid'] = 'NA' 
        globals() ['srr_list_mid'] = 'NA'
        
    elif value == "Intermediate Cell-Type":
        labs_list_exp = labs_list.copy()
        labs_list_exp.insert(0, 'Select Experimental Group')
        x = int(2)
        y = int(2)
        entry_lab_exp = customtkinter.CTkLabel(master=window, text="Experimental Group: ")
        entry_lab_exp.grid(row=y - 1, column=x, padx=10)
        added_widgets.append(entry_lab_exp)
        selected_option_exp = customtkinter.StringVar(value=labs_list_exp[0])
        dropdown_exp = customtkinter.CTkOptionMenu(master=window, values=labs_list_exp,
                                                   command=lambda value3: end_select(value3, 'exp', labs, x, y, value),
                                                   variable=selected_option_exp)
        dropdown_exp.grid(row=y, column=x, padx=5, pady=10)
        added_widgets.append(dropdown_exp)
        labs_list_end = labs_list.copy()
        labs_list_end.insert(0, 'Select Target Group')
        x_end=int(2)
        y_end=int(6)
        entry_lab_end = customtkinter.CTkLabel(master=window, text="Target Group: ")
        entry_lab_end.grid(row=y_end - 1, column=x_end, padx=10)
        added_widgets.append(entry_lab_end)
        selected_option_end = customtkinter.StringVar(value=labs_list_end[0])
        dropdown_end = customtkinter.CTkOptionMenu(master=window, values=labs_list_end,
                                                   command=lambda value4: end_select(value4, 'end', labs, x_end, y_end, value),
                                                   variable=selected_option_end)
        dropdown_end.grid(row=y_end, column=x_end, padx=5, pady=10)
        added_widgets.append(dropdown_end)
        labs_list_mid = labs_list.copy()
        labs_list_mid.insert(0, 'Select Intermediate Group')
        x_mid=int(2)
        y_mid=int(4)
        entry_lab_mid = customtkinter.CTkLabel(master=window, text="Intermediate Group: ")
        entry_lab_mid.grid(row=y_mid - 1, column=x_mid, padx=10)
        added_widgets.append(entry_lab_mid)
        selected_option_mid = customtkinter.StringVar(value=labs_list_mid[0])
        dropdown_mid = customtkinter.CTkOptionMenu(master=window, values=labs_list_mid,
                                                   command=lambda value4: end_select(value4, 'mid', labs, x_mid, y_mid, value),
                                                   variable=selected_option_mid)
        dropdown_mid.grid(row=y_mid, column=x_mid, padx=5, pady=10)
        added_widgets.append(dropdown_mid)
        globals() ['man_input_list'] = 'NA' 
        globals() ['srr_lab_exp'] = 'NA'
        globals() ['srr_input_list'] = 'NA' 
        globals() ['srr_lab_end'] = 'NA' 
        globals() ['srr_input_list_exp'] = 'NA' 
        globals() ['srr_lab_mid'] = 'NA' 
        globals() ['srr_list_mid'] = 'NA'
    elif value == "Manual":
        entry_lab_end=customtkinter.CTkLabel(master=window, text="Please enter ENSEMBLE gene names for inputs")
        entry_lab_end.grid(row=1, column=2, padx=10)
        added_widgets.append(entry_lab_end)
        manual_widg = customtkinter.CTkTextbox(window, wrap="word", height=100, width=500)
        manual_widg.grid(row=2, column=2, columnspan=2, padx=5, pady=10)
        added_widgets.append(manual_widg)
        man_input = manual_widg.get("1.0", tk.END)  # Read from a Text widget
        man_input_list= man_input.split(",")
        man_input_list=[item.strip() for item in man_input_list]
        labs_list_exp = labs_list.copy()
        labs_list_exp.insert(0, 'Select Experimental Group')
        x = int(2)
        y = int(4)
        entry_lab_exp = customtkinter.CTkLabel(master=window, text="Experimental Group: ")
        entry_lab_exp.grid(row=y - 1, column=x, padx=10)
        added_widgets.append(entry_lab_exp)
        selected_option_exp = customtkinter.StringVar(value=labs_list_exp[0])

        dropdown_exp = customtkinter.CTkOptionMenu(master=window, values=labs_list_exp,
                                                   command=lambda value3: end_select(value3, 'exp', labs, x, y, value),
                                                   variable=selected_option_exp)
        dropdown_exp.grid(row=y, column=x, padx=5, pady=10)
        added_widgets.append(dropdown_exp)
        labs_list_end = labs_list.copy()
        labs_list_end.insert(0, 'Select Experimental Group')
        x_end=int(2)
        y_end=int(6)
        entry_lab_end = customtkinter.CTkLabel(master=window, text="Target Group: ")
        entry_lab_end.grid(row=y_end - 1, column=x_end, padx=10)
        added_widgets.append(entry_lab_end)
        selected_option_end = customtkinter.StringVar(value=labs_list_end[0])

        dropdown_end = customtkinter.CTkOptionMenu(master=window, values=labs_list_end,
                                                   command=lambda value4: end_select(value4, 'end', labs, x_end, y_end, value),
                                                   variable=selected_option_end)
        dropdown_end.grid(row=y_end, column=x_end, padx=5, pady=10)
        added_widgets.append(dropdown_end)
        globals() ['mid_tag'] = 'NA'
        globals() ['mid_labs'] = 'NA'
        globals() ['srr_lab_exp'] = 'NA'
        globals() ['srr_input_list'] = 'NA' 
        globals() ['srr_lab_end'] = 'NA' 
        globals() ['srr_input_list_exp'] = 'NA' 
        globals() ['srr_lab_mid'] = 'NA' 
        globals() ['srr_list_mid'] = 'NA'
        globals() ['big_mid'] = 'NA'
        globals() ['small_mid'] = 'NA'
        globals() ['complete_mid'] = 'NA'

    elif value == "Enter SRRs":
        global srr_input, srr_lab_end, srr_input_mid, srr_lab_mid, srr_input_exp, srr_lab_exp, srr_input_list, srr_input_list_exp, srr_list_mid
        entry_values={}
        def shred(event, widget_name, widget):
            entry_values[widget_name] = widget.get()
            print(f"{widget_name} updated: {entry_values[widget_name]}")

#        def shred(event, entry_name):
#         entry_values[entry_name] = event.widget.get()
#         print(f"{entry_name} updated: {entry_values[entry_name]}")
        entry_lab_end=customtkinter.CTkLabel(master=window, text="Target Group SRR ID numbers")
        entry_lab_end.grid(row=5, column=3, padx=10, columnspan=2)
        added_widgets.append(entry_lab_end)
#        entry.bind("<KeyRelease>", lambda event, name=entry_name: on_entry_change(event, name))
        manual_widg = customtkinter.CTkEntry(window, width=300)
        manual_widg.grid(row=6, column=3, columnspan=2, padx=5, pady=10)
        added_widgets.append(manual_widg)
#        srr_input = manual_widg.get("1.0", tk.END)  # Read from a Text widget
#        srr_input = manual_widg.get()  # Read from a Text widget


        entry_lab_end1=customtkinter.CTkLabel(master=window, text="Target Label")
        entry_lab_end1.grid(row=5, column=6, padx=10)
        added_widgets.append(entry_lab_end1)
        manual_widg1 = customtkinter.CTkEntry(window)
        manual_widg1.grid(row=6, column=6, padx=5, pady=10)
        added_widgets.append(manual_widg1)
#        srr_lab_end = manual_widg1.get()  # Read from a Text widget

        entry_lab_mid=customtkinter.CTkLabel(master=window, text="Intermediate Group SRR ID numbers")
        entry_lab_mid.grid(row=3, column=3, padx=10, columnspan=2)
        added_widgets.append(entry_lab_mid)
        manual_widg3 = customtkinter.CTkEntry(window, width=300,  placeholder_text="NA if not applicable")
        manual_widg3.grid(row=4, column=3, columnspan=2, padx=5, pady=10)
        added_widgets.append(manual_widg3)
#        srr_input_mid = manual_widg3.get()  # Read from a Text widget


        entry_lab_mid1=customtkinter.CTkLabel(master=window, text="Intermediate Label")
        entry_lab_mid1.grid(row=3, column=6, padx=10)
        added_widgets.append(entry_lab_mid1)
        manual_widg4 = customtkinter.CTkEntry(window)
        manual_widg4.grid(row=4, column=6, columnspan=1, padx=5, pady=10)
        added_widgets.append(manual_widg4)
#        srr_lab_mid = manual_widg4.get()  # Read from a Text widget


        entry_lab_exp=customtkinter.CTkLabel(master=window, text="Experimental Group SRR ID numbers")
        entry_lab_exp.grid(row=1, column=3, padx=10, columnspan=2)
        added_widgets.append(entry_lab_exp)
        srr_widg_exp = customtkinter.CTkEntry(window, width=300)
        srr_widg_exp.grid(row=2, column=3, columnspan=2, padx=5, pady=10)
        added_widgets.append(srr_widg_exp)
#        srr_input_exp = srr_widg_exp.get()  # Read from a Text widget


        entry_lab_exp1=customtkinter.CTkLabel(master=window, text="Experimental Label")
        entry_lab_exp1.grid(row=1, column=6, padx=10)
        added_widgets.append(entry_lab_exp1)
        manual_widg2 = customtkinter.CTkEntry(window)
        manual_widg2.grid(row=2, column=6, padx=5, pady=10)
        added_widgets.append(manual_widg2)
#        srr_lab_exp = manual_widg1.get()  # Read from a Text widget
#        entry_names = ["manual_widg", "manual_widg1", "manual_widg3", "manual_widg4", "srr_widg_exp", "manual_widg2"]
        entry_names = [
            ("manual_widg", manual_widg),
            ("manual_widg1", manual_widg1),
            ("manual_widg3", manual_widg3),
            ("manual_widg4", manual_widg4),
            ("srr_widg_exp", srr_widg_exp),
            ("manual_widg2", manual_widg2),
        ]

        for widget_name, widget in entry_names:
            widget.bind("<KeyRelease>", lambda event, wn=widget_name, w=widget: shred(event, wn, w))

        # Updating global variables from entry_values
        def update_globals():
            global srr_input, srr_lab_end, srr_input_mid, srr_lab_mid, srr_input_exp, srr_lab_exp, mid_tag, mid_labs, man_input_list, big_mid, small_mid, complete_mid, big_exp, big_end, small_end, small_exp, complete_end, complete_exp, srr_input_list_exp, srr_input_list, srr_list_mid
            srr_input = entry_values.get("manual_widg", "")
            srr_lab_end1 = entry_values.get("manual_widg1", "")
            srr_input_mid = entry_values.get("manual_widg3", "")
            srr_lab_mid1 = entry_values.get("manual_widg4", "")
            srr_input_exp = entry_values.get("srr_widg_exp", "")
            srr_lab_exp1 = entry_values.get("manual_widg2", "")
            mid_tag='NA'
            mid_labs='NA'
            man_input_list='NA'
            big_mid='NA'
            small_mid='NA'
            complete_mid='NA'
            big_exp='NA'
            big_end='NA'
            small_end='NA'
            small_exp='NA'
            complete_end='NA'
            complete_exp='NA'
            srr_lab_end = "".join(srr_lab_end1)
            print('SRR_LAB_END')
            print(srr_lab_end)
            srr_lab_mid = "".join(srr_lab_mid1)
            srr_lab_exp = "".join(srr_lab_exp1)
            srr_input_list_exp= srr_input_exp.split(",")
            srr_input_list_exp=[item.strip() for item in srr_input_list_exp]
            if srr_input_mid != "NA":
             srr_list_mid= srr_input_mid.split(",")
             srr_list_mid=[item.strip() for item in srr_list_mid]
            else:
             srr_list_mid = "NA"
            srr_input_list= srr_input.split(",")
            srr_input_list=[item.strip() for item in srr_input_list]

        # Add a button to save and print values
        save_button = customtkinter.CTkButton(window, text="Save Values", command=update_globals)
        save_button.grid(row=7, column=6, pady=10)
        added_widgets.append(save_button)
#        for entry in entry_names:
#         entry_name = entry.winfo_name()  # Use the widget's name as entry name
#         entry.bind("<KeyRelease>", partial(shred, entry_name=entry_name))
#        print(entry_values)        
#        srr_input = []        
#        srr_input = entry_values["manual_widg"]
#        srr_lab_end = []
#        srr_lab_end = entry_values["manual_widg1"]
#        srr_input_mid = []
#        srr_input_mid = entry_values["manual_widg3"]
#        srr_lab_mid = []
#        srr_lab_mid = entry_values["manual_widg4"]
#        srr_input_exp = []
#        srr_input_exp = entry_values["srr_widg_exp"]
#        srr_lab_exp = []
#        srr_lab_exp = entry_values["manual_widg2"]
            
    return

#added_widgets = ['lbl_pop','dropdown', 'entry_lab_exp', 'entry_exp', 'entry_lab_end', 'entry_end', 'tag_lab_exp', 'tag_exp', 'tag_lab_end', 'tag_end', #'btn_update', 'entry_mid', 'entry_lab_end3', 'tag_lab_mid', 'tag_mid', 'manual_widg']
def reset_selections(added_widgets):
    global exp_input_list, end_input_list, man_input_list, mid_input_list, exp_tag_list, end_tag_list, mid_tag_list, CRISPR, fitted
    exp_input_list = []
    end_input_list = []
    man_input_list = []
    mid_input_list = []
    exp_tag_list = []
    end_tag_list = []
    mid_tag_list = []
    CRISPR = ""
    fitted=None
    selected_option_CRISPR_method.set(options_CRISPR_method[0])
#    selected_option_data_cat.set(options_data_cat[0])
    selected_option_input.set(options_input[0])
    # Remove dynamically added widgets
    for widget in added_widgets:
        widget.grid_forget()
    added_widgets.clear()

from tkinter import *
import customtkinter
#from customtkinter.windows.widgets import ctk_image
window = customtkinter.CTk()
customtkinter.set_appearance_mode("dark")
#customtkinter.set_default_color_theme("orange")
current_directory = os.getcwd()
customtkinter.set_default_color_theme(current_directory+"/orange2.json")
window.iconbitmap(current_directory+"/Easton_icon1.ico")
window.title("CRISPR-GEM")
window.resizable(width=True, height=True)
window.configure()

sidebar_frame = customtkinter.CTkFrame(window, width=140, corner_radius=0)
sidebar_frame.grid(row=0, column=0, rowspan=15, columnspan=2, sticky="nsew")
sidebar_frame.grid_rowconfigure(4, weight=1)


title_lab=customtkinter.CTkLabel(master=window, text="CRISPR-GEM", font=('Roboto', 36))
title_lab.grid(row=0, column=0, columnspan=2, padx=5, pady=10)

options_CRISPR_method = ["Select CRISPR Method", "CRISPR-ko", "CRISPR-ki", "CRISPRa", "CRISPRi", "CRISPRa/i"]
selected_option_CRISPR_method = customtkinter.StringVar(value=options_CRISPR_method[0])
dropdown_CRISPR_method = customtkinter.CTkOptionMenu(master=window, values=options_CRISPR_method, command=on_select1, variable=selected_option_CRISPR_method)
dropdown_CRISPR_method.grid(row=1, column=0, columnspan=2, padx=5, pady=10)

try:
    GNAR
    print("GNAR is already loaded.")
except NameError:
    GNAR = pd.read_csv(current_directory+"/big.csv", index_col=0)
    print("Gnarly")

lbl_space1 = customtkinter.CTkLabel(master=window, text="")
lbl_space1.grid(row=3, column=0)
new_widg=[]
def runner(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, CRISPR, window, GNAR, srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid):
  global new_widg
  output_genes, input_genes, new_widg = run_deseq2(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, CRISPR, window, GNAR, srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid)
  return output_genes, input_genes, new_widg

btn_update = customtkinter.CTkButton(window, text="Get Inputs", command=lambda: run_deseq2(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, CRISPR, window, GNAR, srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid))
btn_update.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

get_gene = customtkinter.CTkEntry(master=window, width=50)
get_gene.grid(row=5, column=0, sticky="s", padx=5, pady=5)

fitted=None
df1=None
new_widg1=[]
processed_data=None
def show_graph(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, get_gene, fitted1, df1, processed_data, window, GNAR, srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid):
    global result_df
    if 'fitted' not in globals():
      global fitted
    fitted, df1, processed_data = run_graph(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, get_gene, fitted1, df1, processed_data, window, GNAR, srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid)
    return fitted, df1, processed_data

btn_graph = customtkinter.CTkButton(window, text="Graph Gene", command=lambda: show_graph(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, get_gene, fitted, df1, processed_data, window, GNAR, srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid))
btn_graph.grid(row=5, column=1, padx=5, pady=5)

def show_umap(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, fitted1, df1, processed_data, window, GNAR, srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid):
    global result_df
    if 'fitted' not in globals():
      global fitted
    fitted, df1, processed_data = run_umap(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, fitted1, df1, processed_data, window, GNAR, srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid)
    return fitted, df1, processed_data

btn_umap = customtkinter.CTkButton(window, text="Umap Plot", command=lambda: show_umap(big_exp, small_exp, big_end, small_end, big_mid, small_mid, complete_exp, complete_end, complete_mid, man_input_list, transcript_to_gene, fitted, df1, processed_data, window, GNAR, srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid))
btn_umap.grid(row=6, column=0, columnspan=2, padx=5, pady=5)

if not os.path.exists(current_directory+"/saveas_.txt"): 
	open('saveas_.txt', 'w').close()
with open(current_directory+"/saveas_.txt", "r") as file:
    saveas_ = file.read()
print(saveas_)
save_path=os.path.join(current_directory, saveas_)

options_input = ["Select Data", "End-to-End", "Intermediate Cell-Type", "Manual", "Enter SRRs"]
selected_option_input = customtkinter.StringVar(value="Select Data")
def hello(value, window):
#   global exp_labs, exp_tag, end_labs, end_tag, mid_labs, mid_tag, man_input_list
#   value, exp_labs, exp_tag, end_labs, end_tag, mid_labs, mid_tag, man_input_list=input_type(value, window)
   global srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid
   srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid = input_type(value, window)
#   return srr_lab_exp, srr_input_list, srr_lab_end, srr_input_list_exp, srr_lab_mid, srr_list_mid
   return
#   return exp_labs, exp_tag, end_labs, end_tag, mid_labs, mid_tag, man_input_list
dropdown_input = customtkinter.CTkOptionMenu(master=window, values=options_input, command=lambda value: input_type(value, window), variable=selected_option_input)
dropdown_input.grid(row=3, column=0, columnspan=2, padx=5, pady=10)

btn_reset = customtkinter.CTkButton(window, text="Reset", command=lambda: reset_selections(added_widgets+new_widg))
btn_reset.grid(row=7, column=0, columnspan=2, padx=5, pady=5)

window.mainloop()