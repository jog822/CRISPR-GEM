import tkinter as tk
from tkinter import ttk
import pandas as pd
import os

def display_dataset(inn, window):
    import pandas as pd
    current_directory = os.getcwd()
    if inn == 'inputs':
          dataset=input_genes=pd.read_csv(current_directory+'/test_input.csv', index_col=0)
    elif inn=='outputs':
          dataset=pd.read_csv(current_directory+'/test_output.csv', index_col=0)
    dataset_with_index = dataset.reset_index()

    window2 = tk.Toplevel(window)
    window2.title("Interactive Dataset")
    window2.configure(bg='navy')

    tree = ttk.Treeview(window2)
    tree["columns"] = list(dataset_with_index.columns)
    tree["show"] = "headings"

    tree.heading("index", text="Index")
    for column in dataset_with_index.columns:
        tree.heading(column, text=column)

    for _, row in dataset_with_index.iterrows():
        tree.insert("", "end", values=list(row))

    tree.grid(row=1, column=0, sticky="nsew")
    window2.grid_rowconfigure(1, weight=1)
    window2.grid_columnconfigure(0, weight=1)

    def sort_dataset():
        selected_column = comobox_column.get()
        ascending = int(radiobutton_asc.get())
        sorted_dataset = dataset.reset_index().sort_values(by=selected_column, ascending=ascending, ignore_index=True)
        # Clear existing rows
        for child in tree.get_children():
            tree.delete(child)
        # Insert sorted data
        for _, row in sorted_dataset.iterrows():
            tree.insert("", "end", values=list(row))

    ttk.Label(window2, text = "Selected Genes",  
          background = 'navy', foreground ="deep sky blue", anchor="w",  
          font = ("Times New Roman", 15)).grid(row = 0, column = 0) 
  
    ttk.Label(window2, text = "Select Column to filter:", foreground="deep sky blue", background="navy", anchor="w",
          font = ("Times New Roman", 10)).grid(column = 0, 
          row = 2, pady = 5) 

    n = tk.StringVar() 
    comobox_column = ttk.Combobox(window2, width = 15, textvariable = n, justify="left") 
  
    comobox_column['values'] = dataset.columns.tolist()
    comobox_column.grid(column = 0, row = 3, pady=5) 
    comobox_column.current() 

    ttk.Label(window2, text = "Search for genes:",anchor="w", foreground="deep sky blue", background="navy", font = ("Times New Roman", 10)).grid(column = 0, 
          row = 6, pady = 10) 

    radiobutton_asc = tk.StringVar(value="1")  # Default to ascending order

    style = ttk.Style()
    sty = ttk.Style(window)
    sty.configure("TRadiobutton", background="navy",
              foreground="deep sky blue", anchor="w",
              justify='left')
    ttk.Radiobutton(window2, text="Ascending", style="white.TRadiobutton", variable=radiobutton_asc, value="1", command=sort_dataset).grid(row=4, column=0)

    ttk.Radiobutton(window2, text="Descending", style="white.TRadiobutton", variable=radiobutton_asc, value="0", command=sort_dataset).grid(row=5, column=0)

    def filter_rows(event=None):
        search_query = entry_search.get()
        filtered_dataset = dataset[dataset.index.str.contains(search_query, case=False)]
        for child in tree.get_children():
            tree.delete(child)
        for idx, row in filtered_dataset.iterrows():
            tree.insert("", "end", values=[idx] + list(row))    
    entry_search = tk.Entry(window2, width=10)
    entry_search.grid(row=7, column=0, pady=10)
    entry_search.bind("<KeyRelease>", filter_rows)


def input_type(value):
    global exp_input_list, end_input_list, man_input_list, mid_input_list, exp_tag_list, end_tag_list, mid_tag_list
    def get_inputs():
        global exp_input_list, end_input_list, man_input_list, mid_input_list, exp_tag_list, end_tag_list, mid_tag_list
        if value == "Manual":
            man_input = manual_widg.get("1.0", tk.END)  # Read from a Text widget
            man_input_list= man_input.split(",")
            man_input_list=[item.strip() for item in man_input_list]
            mid_input_list="NA"
            mid_tag_list="NA"
            exp_input = entry_exp.get()
            end_input = entry_end.get()
            exp_input_list = exp_input.split(",")
            exp_input_list=[item.strip() for item in exp_input_list]
            end_input_list = end_input.split(",")
            end_input_list=[item.strip() for item in end_input_list]
            print(exp_input_list)
            print(end_input_list)
            print(man_input_list)
            exp_tag = tag_exp.get()
            end_tag = tag_end.get()
            exp_tag_list = exp_tag.split(",")
            exp_tag_list=[item.strip() for item in exp_tag_list]
            end_tag_list = end_tag.split(",")
            end_tagt_list=[item.strip() for item in end_tag_list]
            return man_input_list, exp_input_list, end_input_list, exp_tag_list, end_tag_list
        elif value == "End-to-End":
            exp_input = entry_exp.get()
            end_input = entry_end.get()
            exp_input_list = exp_input.split(",")
            exp_input_list=[item.strip() for item in exp_input_list]
            end_input_list = end_input.split(",")
            end_input_list=[item.strip() for item in end_input_list]
            mid_input_list="NA"
            mid_tag_list="NA"
            man_input_list="NA"
            print(exp_input_list)
            print(end_input_list)
            exp_tag = tag_exp.get()
            end_tag = tag_end.get()
            exp_tag_list = exp_tag.split(",")
            exp_tag_list=[item.strip() for item in exp_tag_list]
            end_tag_list = end_tag.split(",")
            end_tag_list=[item.strip() for item in end_tag_list]
            return exp_input_list, end_input_list, exp_tag_list, end_tag_list
        elif value == "Intermediate Cell-Type":
            exp_input = entry_exp.get()
            end_input = entry_end.get()
            mid_input=entry_mid.get()
            exp_input_list = exp_input.split(",")
            exp_input_list=[item.strip() for item in exp_input_list]
            end_input_list = end_input.split(",")
            end_input_list=[item.strip() for item in end_input_list]
            mid_input_list=mid_input.split(",")
            mid_input_list=[item.strip() for item in mid_input_list]
            man_input_list="NA"
            print(exp_input_list)
            print(end_input_list)
            print(mid_input_list)
            exp_tag = tag_exp.get()
            end_tag = tag_end.get()
            mid_tag=tag_mid.get()
            exp_tag_list = exp_tag.split(",")
            exp_tag_list=[item.strip() for item in exp_tag_list]
            end_tag_list = end_tag.split(",")
            end_tag_list=[item.strip() for item in end_tag_list]
            mid_tag_list=mid_tag.split(",")
            mid_tag_list=[item.strip() for item in mid_tag_list]
            print(exp_tag_list)
            print(end_tag_list)
            print(mid_tag_list)
            return exp_input_list, end_input_list, mid_input_list, exp_tag_list, end_tag_list, mid_tag_list
            
    if value == "End-to-End":
        entry_lab_exp=tk.Label(master=window, text="Please Enter Experimental Group SRR #'s (at least 3)'", bg='orange')
        entry_lab_exp.grid(row=6, column=2, padx=5)
        entry_exp = tk.Entry(master=window, width=30)
        entry_exp.grid(row=7, column=2, sticky="e", padx=5, pady=5)
        
        entry_lab_end
