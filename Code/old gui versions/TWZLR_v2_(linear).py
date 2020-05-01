# -*- coding: utf-8 -*-
"""
TODO: add homology arms, combinatorial, integrate with sequencer (tomorrow)
Add: -index labels on gene list, pull homology arm, load data
Errors: -THROWN if gene ranges not included and gene added -negative numbers inserted anywhere -insertion or gene ranges inverted -ask Chris about insertion circularity (not important) 

Author: C
"""
import os
from tkinter import *
from tkinter import filedialog
from Bio import Entrez

Entrez.email = 'A.N.Other@example.com'


def get_gene_name(gene_ID):
    handle = Entrez.efetch(db = "gene", id = gene_ID, retmode ='text') #get information of specific NCBI ID
    raw_file = handle.read() #read fasta file from handle
    parsed_data = raw_file.splitlines()[0:]
    gene_name = parsed_data[1][3:]
    return gene_name



def get_gene_sequence(gene_ID,gene_range): #this function is definietly some black magic 
    gene_handle = Entrez.efetch(db = "gene", id = gene_ID, retmode ='text') #get information of specific NCBI ID
    gene_raw = gene_handle.read() #read fasta file from handle
    gene_data = gene_raw.splitlines()[0:]
    
    name = gene_data[2]
    
    # ignore this demonic conjuring, it works
    annotation_indices = [i for i, s in enumerate(gene_data) if 'Annotation' in s]
    annotation = gene_data[annotation_indices[0]]
    annotation = re.split('\s',annotation)
    annotation_id_indices = [i for i, s in enumerate(annotation) if 'NC' in s]
    annotation_id = annotation[annotation_id_indices[0]]
    
    genome_handle = Entrez.efetch(db = 'nucleotide', id = annotation_id, retmode='text', rettype='fasta')
    genome_raw = genome_handle.read()
    genome_data = genome_raw.splitlines()[1:]
    
    genome_sequence = ''.join(genome_data)
    gene_sequence = genome_sequence[gene_range[0]-1 : gene_range[1]]
    return gene_sequence



def check_params(param_config):
    
    if not param_config.get('Homology Length'):
        param_config['Homology Length'] = 25 
    if not param_config.get('Polymerase'):
        param_config['Polymerase'] = 'Polymerase 1'
    if param_config.get('Chromosome') == 0: #POTENTIAL BUG IF CHANGED MANUALLY IN CONSOLE BUT WILL THROW ERROR COME SEQUENCING STEP 
        return 0 
    if not param_config.get('Cut Range'):
        return 0
    print (param_config)


def get_vnat_genome(chromosome):
    chromo_sequences = []
    chromo_IDs = ['NZ_CP016347.1','NZ_CP016348.1']
    chromo_handle = Entrez.efetch(db = "nucleotide", id = chromo_IDs[chromosome-1], retmode ='fasta', rettype='fasta') #get information of specific NCBI ID
    raw_file = chromo_handle.read() #read fasta file from handle
    raw_sequence = raw_file.splitlines()[1:]
    chromo_sequence=(''.join(raw_sequence))
    return chromo_sequence


def get_homology_arms(chromosome, chromo_sequence, cut_range, homology_length):
    homology_arms = []
    if type(cut_range) == int:
        homology_arms.append(chromo_sequence[cut_range-homology_length:cut_range])
        homology_arms.append(chromo_sequence[cut_range:cut_range+homology_length])
    else:    
        homology_arms.append(chromo_sequence[cut_range[0]-homology_length:cut_range[0]])
        homology_arms.append(chromo_sequence[cut_range[-1]:cut_range[-1]+homology_length])
    return homology_arms


def get_gene_sequences(gene_IDs,gene_ranges):
    fluorescent_sequence = 'ATATATATATATATATATTATATATAT' #INSERT Fluorescent SEQUENCES HERE
    for  i in range (0,len(gene_IDs)):
        if gene_IDs[i] == '0000000' :
            gene_sequences.append(fluorescent_sequence)
        else:
            gene_sequences.append(get_gene_sequence(gene_IDs[i],gene_ranges[i]))
    
    


class TWIZLR_int:



    def __init__(self,main_int):
        
        self.prompt_id = Label(main_int, text='Insert gene ID for gene')
        self.prompt_id.place (x=250+350,y=50)
        self.insert_id = Entry(main_int,bd=2)
        self.insert_id.place(x=250+350, y=75)
        
        self.prompt_range = Label(main_int, text='Insert base range for gene in parent genome (beginning,end) ')
        self.prompt_range.place (x=250+350,y=125)
        self.insert_range = Entry(main_int,bd=2)
        self.insert_range.place(x=250+350, y=150)
        
        self.show_gene = Button(main_int, text='Add Gene', command = lambda: self.add_gene(self.insert_id.get(), self.insert_index.get(), self.insert_range.get() ) )
        self.show_gene.place(x=250+350, y=205) #need to lambda this for the actual addition of parameters (add id, range, and where in list to be inserted )
        self.prompt_index = Label(main_int, text='at index (0-indexing) \n in current gene list')
        self.prompt_index.place (x=250+450,y=200)
        self.insert_index = Entry(main_int,bd=2)
        self.insert_index.place(x=250+575, y=210)
        self.insert_index.insert(0,0)
        
        self.remove_entry = Button(main_int, text='Remove Entry', command = lambda: self.remove_gene(self.insert_removal_index.get() ) )
        self.remove_entry.place(x=250+350, y=265) #need to lambda this for the actual addition of parameters (add id, range, and where in list to be inserted )
        self.prompt_removal_index = Label(main_int, text='at index (0-indexing) \n in current gene list')
        self.prompt_removal_index.place (x=250+450,y=260)
        self.insert_removal_index = Entry(main_int,bd=2)
        self.insert_removal_index.place(x=250+575, y=270)

        self.gene_list_title = Label(main_int, text = 'Genes to be added in sequential order:')
        self.gene_list_title.place(x=25,y=25)

        self.gene_list_disp = StringVar()
        self.gene_display = Label(main_int, textvariable=self.gene_list_disp)
        self.gene_display.place(x=50,y=50)
            
        self.set_config = Button(main_int, text='Open Settings', command=self.get_config,bd=2)
        self.set_config.place(x=250+450,y=375)

        self.run_tool = Button(main_int, text ='Run tool', command = self.run_master, bd=2)
        self.run_tool.place(x=250+450,y=420)

    
        
    def add_gene(self,gene_ID,index,base_range):
        index = int(index)
        gene_name = get_gene_name(gene_ID)
        gene_IDs.insert(index,gene_ID)
        gene_ranges.insert(index,eval(base_range))
        gene_list.insert(index,gene_name)
        self.gene_list_disp.set("\n".join(gene_list) )
        self.insert_index.delete(0,END)
        self.insert_index.insert(0,index+1)


    def remove_gene(self,index):
        index = int(index)
        del gene_IDs[index]
        del gene_list[index]
        del gene_ranges[index]
        self.gene_list_disp.set("\n".join(gene_list) )


    def get_config(self):
        global config_int
        config_int=Tk()
        config_int.geometry('512x512')
        init_config(config_int)
        config_int.mainloop()
        
    def run_master(self):
        gene_sequences = get_gene_sequences(gene_IDs, gene_ranges)
        



class init_config:
        
        
    def __init__(self,config_int):
      
        self.chrom_prompt = Label (config_int, text = "V. nat Chromosome to edit:")
        self.chrom_prompt.place(x=25,y=25)
        chromosome = IntVar()
        self.chrom_select_1 = Radiobutton(config_int,text='1',variable=chromosome,value=1).place(x=25, y=50)
        self.chrom_select_2 = Radiobutton(config_int,text='2',variable=chromosome,value=2).place(x=25, y=75)
                    
        self.cut_range_prompt = Label(config_int, text = "Enter range of basepairs to cut (beginning,end) or a single insertion point")
        self.cut_range_prompt.place(x=25,y=125)
        self.cut_range_var = Entry(config_int, bd=2)
        self.cut_range_var.place(x=25+400,y=125)
                   
        self.homology_prompt = Label (config_int, text = "Enter primer homology arm length:")
        self.homology_prompt.place(x=25, y=175)
        self.homology_var = Entry(config_int,bd=2)
        self.homology_var.place(x=25+300, y=175)        
        
        polymerase = StringVar()
        self.poly_prompt = self.chrom_prompt = Label (config_int, text = "Select Polymerase:")
        self.poly_prompt.place(x=25, y=225)
        self.poly_select_1 = Radiobutton(config_int,text='Phusion', variable=polymerase, value='Phusion').place(x=25,y=250)
        self.poly_select_2 = Radiobutton(config_int,text='Q5', variable=polymerase, value='Q5').place(x=25,y=275)
        #self.poly_select_3 = Radiobutton(config_int,text='Polymerase 3', variable=polymerase, value='Polymerase 3').place(x=25,y=300)
                
        fluorescence = IntVar()
        self.prompt_fluroscence = Checkbutton(config_int,text = 'Check to use fluorescence cartridges', variable = fluorescence)
        self.prompt_fluroscence.place(x=25,y=310)

        self.output_prompt = Button(config_int,text = "Select results output path", command = self.select_output_path, bd=2)
        self.output_prompt.place(x=175,y=350)

        self.save_prompt = Button(config_int,text='Save Settings', command = lambda: self.save_config(chromosome.get(), self.homology_var.get(), self.cut_range_var.get(), polymerase.get() ) ).place(x=200, y=400)

    

    def save_config(self,chromosome,homology_length,cut_range_var,polymerase):
        global param_config
        param_config = {}
        param_config['Chromosome'] = chromosome
        param_config['Homology Length'] = homology_length
        param_config['Cut Range'] = cut_range_var
        param_config['Polymerase'] = polymerase
        param_config['Path'] = output_path
        param_config['Use Fluorescnece'] = fluorescence
        config_int.destroy()
   

    
    def select_output_path(self):
        output_path = filedialog.askdirectory()
             

            

global gene_list
gene_list = []

global gene_IDs
gene_IDs = []

global gene_ranges
gene_ranges = []

global gene_sequences
gene_sequences = []

global output_path
output_path = os.getcwd()

global config_int
config_int=Tk()
config_int.geometry('640x512')
init_config(config_int)
config_int.mainloop()


run_status = check_params(param_config)

if run_status != 0:
    
    param_config['Cut Range'] = eval(param_config['Cut Range'])
    
    main_int = Tk()
    TWIZLR_int(main_int)
    main_int.geometry('1024x512')
    main_int.mainloop()    
    
    chromo_sequence = get_vnat_genome(param_config['Chromosome'])
    homology_arms = get_homology_arms(param_config['Chromosome'],chromo_sequence,param_config['Cut Range'],param_config['Homology Length'])
    
    
    
else:
    print('Error. The chromosomal ,or cut range, parameters were missing, aborting run.')