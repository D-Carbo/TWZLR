# -*- coding: utf-8 -*-
"""

Author: C.
"""
import os
from tkinter import *
from tkinter import filedialog
from Bio import Entrez

import sequencer as TwzSeq
import graph

import json


Entrez.email = 'A.N.Other@example.com' #set emailto pull NCBI data


def get_gene_name(gene_ID):
    handle = Entrez.efetch(db = "gene", id = gene_ID, retmode ='text') #get gene information of specific NCBI ID
    raw_file = handle.read()                                           #read data  from handle
    parsed_data = raw_file.splitlines()[0:]                            #split text into a list
    gene_name = parsed_data[1][3:]                                     #get name of gene
    return gene_name 



def get_gene_sequence(gene_ID,gene_range): #this function is definietly some black magic 
    
    gene_handle = Entrez.efetch(db = "gene", id = gene_ID, retmode ='text') #get gene information of specific NCBI ID
    gene_raw = gene_handle.read()                                           #read data from handle
    gene_data = gene_raw.splitlines()[0:]
    
    name = gene_data[2]
    
    # list comprehension to extract specific genome annotation ID from gene data
    annotation_indices = [i for i, s in enumerate(gene_data) if 'Annotation' in s]
    annotation = gene_data[annotation_indices[0]]
    annotation = re.split('\s',annotation)
    annotation_id_indices = [i for i, s in enumerate(annotation) if 'NC' in s]
    annotation_id = annotation[annotation_id_indices[0]]
    
    genome_handle = Entrez.efetch(db = 'nucleotide', id = annotation_id, retmode='text', rettype='fasta') #get genome information of specific NCBI ID
    genome_raw = genome_handle.read()                                                                     #read data  from handle
    genome_data = genome_raw.splitlines()[1:]                                                             #remove title line from genome data and split text to list
    
    genome_sequence = ''.join(genome_data)                                                                #join genome into single string
    gene_range = eval(gene_range)                                                                         #force user specified range into type: tuple
    gene_sequence = genome_sequence[gene_range[0]-1 : gene_range[1]]                                      #extract gene range from genome
    return gene_sequence



def check_params(param_config): #set default parameters if not specified, return error if vital parameter is missing
    
    if not param_config.get('Homology Length'):
        param_config['Homology Length'] = 500 
    if not param_config.get('Polymerase'):
        param_config['Polymerase'] = 'Polymerase 1'
    if param_config.get('Chromosome') == 0: #POTENTIAL BUG IF CHANGED MANUALLY IN CONSOLE BUT WILL THROW ERROR COME SEQUENCING STEP 
        return 0 
    if not param_config.get('Cut Range'):
        return 0
    print (param_config)


def get_vnat_genome(chromosome): #get vibrio natriegen genome for specified chromosome (1 or 2)
    chromo_sequences = []
    chromo_IDs = ['NZ_CP016347.1','NZ_CP016348.1']
    chromo_handle = Entrez.efetch(db = "nucleotide", id = chromo_IDs[chromosome-1], retmode ='fasta', rettype='fasta') 
    raw_file = chromo_handle.read()
    raw_sequence = raw_file.splitlines()[1:]
    chromo_sequence=(''.join(raw_sequence))
    return chromo_sequence


def get_homology_arms(chromosome, chromo_sequence, cut_range, homology_length): #get left and right homology arms from v.nat chromosome according to homology length and user specified location
    homology_arms = []
    #pull specific range from chromosome
    if type(cut_range) == int:
        homology_arms.append(chromo_sequence[cut_range-homology_length:cut_range])
        homology_arms.append(chromo_sequence[cut_range:cut_range+homology_length])
    else:    
        homology_arms.append(chromo_sequence[cut_range[0]-homology_length:cut_range[0]])
        homology_arms.append(chromo_sequence[cut_range[-1]:cut_range[-1]+homology_length])
    return homology_arms


def get_gene_sequences(gene_IDs,gene_ranges): #get sequence for all genes inserted
    for i in range (0,len(gene_IDs)):
        gene_sequences.append([])
        for j in range(0,len(gene_IDs[i])):
            gene_sequences[i].insert(j,get_gene_sequence(gene_IDs[i][j],gene_ranges[i][j]))
    del gene_sequences[-1]

    

#class for main interface
class TWZLR_int:



    def __init__(self,main_int):
        
        #insert gene ID
        self.prompt_id = Label(main_int, text='Insert gene ID')
        self.prompt_id.place (x=1+350,y=50)
        self.insert_id = Entry(main_int,bd=2)
        self.insert_id.place(x=1+350, y=75)
        
        #insert base range for gene of specific ID
        self.prompt_range = Label(main_int, text='Insert base range for gene in parent genome (beginning,end) ')
        self.prompt_range.place (x=1+350,y=125)
        self.insert_range = Entry(main_int,bd=2)
        self.insert_range.place(x=1+350, y=150)
        
        #add gene to list and display
        self.show_gene = Button(main_int, text='Add Gene', command = lambda: self.add_gene(self.insert_id.get(), self.insert_layer.get(), self.insert_index.get(), self.insert_range.get() ), bd=2) #add gene at layer and subindex displayed 
        self.show_gene.place(x=1+350, y=205) 
        self.prompt_index = Label(main_int, text='at index (0-indexing) \n in layer shown below')
        self.prompt_index.place (x=1+450,y=200)
        self.insert_index = Entry(main_int,bd=2)
        self.insert_index.place(x=1+575, y=210)
        self.insert_index.insert(0,0) #initialize entry variable to 0
        
        #add layer to list and display
        self.show_layer = Button(main_int, text='Add Layer', command = lambda: self.add_layer(self.insert_layer.get() ), bd=2) #add layer at index displayed
        self.show_layer.place(x=1+350, y=260)
        self.prompt_layer = Label(main_int,text = 'at index (0-indexing) \n in current list ')
        self.prompt_layer.place(x=1+450,y=255)
        self.insert_layer = Entry(main_int, bd=2)
        self.insert_layer.place(x=1+575,y=260)
        self.insert_layer.insert(0,0) #initialize entry variable to 0

        #remove entry from list and display
        self.remove_entry = Button(main_int, text='Remove Entry', command = lambda: self.remove_gene(self.insert_layer.get(), self.insert_index_removal.get() ), bd=2) #remove entry at layer and subindex displayed
        self.remove_entry.place(x=1+350, y=305) 
        self.prompt_index_removal = Label(main_int, text='at index (0-indexing) \n in layer shown above')
        self.prompt_index_removal.place (x=1+450,y=300)
        self.insert_index_removal = Entry(main_int,bd=2)
        self.insert_index_removal.place(x=1+575, y=305)

        #delte layer
        self.remove_layer = Button(main_int, text = 'Delete Layer', command = lambda: self.delete_layer(self.insert_layer_removal.get() ), bd=2) #delete layer at index displayed
        self.remove_layer.place(x=1+350, y=352)
        self.prompt_layer_removal = Label(main_int, text = 'at index (0-index) \n in current list')
        self.prompt_layer_removal.place(x=1+450, y=345)
        self.insert_layer_removal = Entry (main_int, bd=2)
        self.insert_layer_removal.place(x=1+575, y=352)

        #list title 
        self.gene_list_title = Label(main_int, text = 'Genes to be added in sequential order:')
        self.gene_list_title.place(x=25,y=25)

        #display gene list
        self.gene_list_disp = StringVar()
        self.gene_display = Label(main_int, textvariable=self.gene_list_disp)
        self.gene_display.place(x=50,y=50)

        #open setting configuration menu
        self.set_config = Button(main_int, text='Open Settings', command=self.get_config, bd=2)
        self.set_config.place(x=1+450,y=475)

        self.run_tool = Button(main_int, text ='Run tool', command = self.run_master)
        self.run_tool.place(x=1+462,y=525)

    
    #function to add gene and relevant information to data lists     
    def add_gene(self,gene_ID,index, subindex, base_range):
        index = int(index)                                                                  #index of layer
        subindex = int(subindex)                                                            #index of gene
        gene_name = get_gene_name(gene_ID)                                                  
        gene_IDs[index].insert(subindex,gene_ID)                                            #add gene ID to list
        gene_ranges[index].insert(subindex,base_range)                                      #add gene base range to list
        gene_list[index].insert(subindex,gene_name)                                         #add gene name to list
        self.gene_list_disp.set('\n'.join(', '.join(map(str,sl) ) for sl in gene_list) )    #update list of gene names
        self.insert_index.delete(0,END)                                                     #clear entry for gene insertion index 
        self.insert_index.insert(0,subindex+1)                                              #make entry previous+1 to allow for immediate insertion


    def add_layer(self, index):
        index = int(index)
        #add empty sublists for inseriton to prevent indexing error
        gene_list.insert(index,[])
        gene_IDs.insert(index,[])
        gene_ranges.insert(index,[])
        self.gene_list_disp.set('\n'.join(', '.join(map(str,sl) ) for sl in gene_list) ) 
        self.insert_index.delete(0,END)
        self.insert_index.insert(0,0)


    def remove_gene(self,index,subindex):
        index = int(index)
        subindex = int(subindex)
        #delete relevant data from lists
        del gene_IDs[index][subindex]
        del gene_ranges[index][subindex]
        del gene_list[index][subindex]
        self.gene_list_disp.set('\n'.join(', '.join(map(str,sl) ) for sl in gene_list) )
    

    def delete_layer(self,index):
        index = int(index)
        #delete layer from lists
        del gene_IDs[index]
        del gene_ranges[index]
        del gene_list[index]
        self.gene_list_disp.set('\n'.join(', '.join(map(str,sl) ) for sl in gene_list) )

    def get_config(self): #initialize and run setting configuration interface
        global config_int
        config_int=Tk()
        config_int.geometry('640x512')
        init_config(config_int)
        config_int.mainloop()
        
    def run_master(self): #master run tool 
        main_int.destroy()
        gene_sequence = get_gene_sequences(gene_IDs, gene_ranges) #get base sequences for final gene list
        


#class for settings interface
class init_config:
        
        
    def __init__(self,config_int):
        
        #chromosome selection
        self.chrom_prompt = Label (config_int, text = "V. nat Chromosome to edit:")
        self.chrom_prompt.place(x=25,y=25)
        chromosome = IntVar()
        self.chrom_select_1 = Radiobutton(config_int,text='1',variable=chromosome,value=1).place(x=25, y=50)
        self.chrom_select_2 = Radiobutton(config_int,text='2',variable=chromosome,value=2).place(x=25, y=75)
        
        #construct insertion point 
        self.cut_range_prompt = Label(config_int, text = "Enter range of basepairs to cut (beginning,end) or a single insertion point")
        self.cut_range_prompt.place(x=25,y=125)
        self.cut_range_var = Entry(config_int, bd=2)
        self.cut_range_var.place(x=25+400,y=125)

        #homology arm length           
        self.homology_prompt = Label (config_int, text = "Enter primer homology arm length:")
        self.homology_prompt.place(x=25, y=175)
        self.homology_var = Entry(config_int,bd=2)
        self.homology_var.place(x=25+300, y=175)        
        
        #polymerase selection
        polymerase = StringVar()
        self.poly_prompt = self.chrom_prompt = Label (config_int, text = "Select Polymerase:")
        self.poly_prompt.place(x=25, y=225)
        self.poly_select_1 = Radiobutton(config_int,text='Phusion', variable=polymerase, value='Phusion').place(x=25,y=250)
        self.poly_select_2 = Radiobutton(config_int,text='Q5', variable=polymerase, value='Q5').place(x=25,y=275)
        #self.poly_select_3 = Radiobutton(config_int,text='Polymerase 3', variable=polymerase, value='Polymerase 3').place(x=25,y=300)
        polymerase.set('Phusion')

        #toggle adaptors
        adaptor = IntVar()
        self.prompt_adaptor = Checkbutton(config_int,text = 'Check to use adaptors', variable = adaptor)
        self.prompt_adaptor.place(x=25,y=310)

        #set results output path
        self.output_prompt = Button(config_int,text = "Select results output path", command = self.select_output_path, bd=2)
        self.output_prompt.place(x=175,y=350)

        #save setting configuration
        self.save_prompt = Button(config_int,text='Save Settings', command = lambda: self.save_config(chromosome.get(), self.homology_var.get(), self.cut_range_var.get(), polymerase.get(), adaptor.get() ) ).place(x=200, y=400)

    
    #get paramaeters and wrtie to dictionary 
    def save_config(self,chromosome,homology_length,cut_range_var,polymerase, adaptor):
        global param_config
        param_config = {}
        param_config['Chromosome'] = chromosome
        param_config['Homology Length'] = homology_length
        param_config['Cut Range'] = cut_range_var
        param_config['Polymerase'] = polymerase
        param_config['Path'] = output_path
        param_config['Use adaptors'] = adaptor
        config_int.destroy()
        
   
    #get user selected output path 
    def select_output_path(self):
        output_path = filedialog.askdirectory()
             

            

global gene_list
gene_list = [[]]

global gene_IDs
gene_IDs = [[]]

global gene_ranges
gene_ranges = [[]]

global gene_sequences
gene_sequences = [[]]

global output_path
output_path = os.getcwd()

#initial setting configuration run
global config_int
config_int=Tk()
config_int.geometry('640x512')
init_config(config_int)
config_int.mainloop()

#check if parameters meet conditions for run
run_status = check_params(param_config)

if run_status != 0:
    
    #set parameter types as numerical for indexing 
    param_config['Cut Range'] = eval(param_config['Cut Range'])
    if type(param_config['Homology Length'])== str:
        param_config['Homology Length'] = eval(param_config['Homology Length'])
    
    #run main interface
    main_int = Tk()
    TWZLR_int(main_int)
    main_int.geometry('768x768')
    main_int.mainloop()    
    
    #get chromosomal and homology sequences
    chromo_sequence = get_vnat_genome(param_config['Chromosome'])
    homology_arms = get_homology_arms(param_config['Chromosome'],chromo_sequence,param_config['Cut Range'],param_config['Homology Length'])
    homology_length = param_config['Homology Length']

    #append homology sequneces to parts lists
    gene_list.insert(0,['Homology left'])
    gene_list.append(['Homology Right'])
   
    gene_sequences.insert(0, [homology_arms[0]])
    gene_sequences.append([homology_arms[1]])
   
    #get primers for each possible cirucit combination
    
    if param_config['Use adaptors'] == 1: #check if user wants adaptors in circuit (work in progress)
        '69696969696969'#need some function here 
    
    else:
        
        primers_dict = {}
        for i in range(1,len(gene_sequences)-1):
                    primers_dict['Layer '+str(i)] = {}
                    
                    for j in range (0,len(gene_sequences[i]) ):
        
                        for k in range (0,len(gene_sequences[i-1]) ):
                            
                            for l in range (0,len(gene_sequences[i+1]) ):
                                
                                #get primers per combination and add to dictionary
                                primers =  TwzSeq.get_primers(gene_sequences[i-1][k][-homology_length:], gene_sequences[i+1][l][0:homology_length], gene_sequences[i][j] )
                                primers_dict['Layer '+str(i)][gene_list[i-1][k]+'-'+gene_list[i][j]+'-'+gene_list[i+1][l]] = primers     

        #output primers as data file1
        output_file = open(output_path+'\TWZLR_primers.json','w')   
        json_data = json.dumps(primers_dict, indent=5)
        output_file.write(json_data)
        output_file.close()

    #graph combinatorial possibilities of designed circuit  
    graph.crudegraph(gene_list, [], param_config['Use adaptors'])        

else:
    print('Error. The chromosomal ,or cut range, parameters were missing, aborting run.')
