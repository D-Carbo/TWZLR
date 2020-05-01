# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 13:00:08 2020

@author: C
"""

from tkinter import *
from tkinter import filedialog
from Bio import Entrez

Entrez.email = 'A.N.Other@example.com'

def get_gene_name(gene_ID)
    handle = Entrez.efetch(db = "gene", id = gene_ID, retmode ='text') #get information of specific NCBI ID
    raw_file = handle.read() #read fasta file from handle
    parsed_data = raw_file.splitlines()[0:]
    gene_name = parsed_data[2]
    return gene_name

def get_fasta(database_ID):
    
    Entrez.email = 'A.N.Other@example.com'
    handle = Entrez.efetch(db = "gene", id = database_ID, rettype = "fasta", retmode = "text") #get information of specific NCBI ID
    raw_file = handle.read() #read fasta file from handle
    sequence = raw_file.splitlines()[0:]
    return sequence


class TWIZLR_int:
    def __init__(self,main):
        
        self.insert_prompt = Label(main, text = "Enter the basepair number to insert the sequence after:")
        self.insert_prompt.place(x=50,y=150)
        self.insert_var = Entry(main, bd=2)
        self.insert_var.place(x=50+300,y=150)
            
        self.chrom_prompt = Label (main, text = "Which V. nat Chromosome would you like to edit:")
        self.chrom_prompt.place(x=50,y=50)
        chromosome = IntVar()
        self.chrom_select_1 = Radiobutton(main,text='1',variable=chromosome,value=1).place(x=50, y=75)
        self.chrom_select_2 = Radiobutton(main,text='2',variable=chromosome,value=2).place(x=50, y=100)
        
        self.homology_prompt = Label (main, text = "Enter the length of the homology arms for the primers:")
        self.homology_prompt.place(x=50, y=325)
        self.homology_var = Entry(main,bd=2)
        self.homology_var.place(x=50+300, y=325)        
        
        polymerase = StringVar()
        self.poly_prompt = self.chrom_prompt = Label (main, text = "Which Polymerase would you like to use:")
        self.poly_prompt.place(x=50, y=200)
        self.poly_select_1 = Radiobutton(main,text='Polymerase 1', variable=polymerase, value='Polymerase 1').place(x=50,y=225)
        self.poly_select_2 = Radiobutton(main,text='Polymerase 2', variable=polymerase, value='Polymerase 2').place(x=50,y=250)
        self.poly_select_3 = Radiobutton(main,text='Polymerase 3', variable=polymerase, value='Polymerase 3').place(x=50,y=275)
        
        self.output_prompt = Button(main,text = "Select output path", command = self.select_output_path, bd=2)
        self.output_prompt.place(x=175,y=433)
        
        self.input_prompt = Button(main,text = "Select input file", command = self.select_input_file, bd=2)
        self.input_prompt.place(x=175,y=400)
        
        
    
    
    def select_input_file(self):
        spec_int = Tk()
        TWIZLR_int.selection_int(spec_int)
        spec_int.geometry('768x512')
        spec_int.mainloop()               
        
    def select_output_path(self):
        output_path = filedialog.askdirectory()
        print (output_path)            
           

    
    
class selection_int:
    
    
    def __init__(self,aux_1):
                  
        self.prompt_id = Label(aux_1, text='Insert ID for gene to be inserted here').place (x=350,y=50)
        self.insert_id = Entry(aux_1,bd=2).place(x=350, y=75)
        
        self.prompt_range = Label(aux_1, text='Insert Range for gene to be inserted here (beginning,end) ').place (x=350,y=125)
        self.insert_range = Entry(aux_1,bd=2).place(x=350, y=150)
        
        self.add_gene = Button(aux_1, text='Add Gene', command=self.add_gene, bd=2).place(x=375, y=200)
        self.prompt_index = Label(aux_1, text='at index (0-indexing)').place (x=450,y=200)
        self.insert_index = Entry(aux_1,bd=2).place(x=575, y=200)
        
        self.gene_list_disp = StringVar()
        self.gene_list_disp.set('Here is some text') 
        self.gene_display = Label(aux_1, textvariable=self.gene_list_disp)
        self.gene_display.pack()
        
    def add_gene(self):
        gene_list.append(0)
        self.gene_list_disp.set(str(gene_list) )
                        
    
main_int = Tk()
TWIZLR_int(main_int)
main_int.geometry('512x768')

global gene_list
gene_list=[0] 

main_int.mainloop()