# -*- coding: utf-8 -*-
"""
TODO: Errors if gene parameters are not there for calculation, side index list, break button (talk to chris), should add complement functionality , add 0 as default insertion element, negative numbers, ranges being inverted, add functionality to load file

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
    gene_name = parsed_data[2]
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
    for i in range (0,len(gene_IDs)):
        for j in range(0,len(gene_IDs[i])):
            if gene_IDs[i][j] == str('0000000') :
                gene_sequences[i].append(fluorescent_sequence)
            else:
                gene_sequences[i].append(get_gene_sequence(gene_IDs[i],gene_ranges[i]))
        
    


class TWIZLR_int:



    def __init__(self,main_int):
        
        self.prompt_id = Label(main_int, text='Insert gene ID for gene')
        self.prompt_id.place (x=725+350,y=50)
        self.insert_id = Entry(main_int,bd=2)
        self.insert_id.place(x=725+350, y=75)
        
        self.prompt_range = Label(main_int, text='Insert base range for gene in parent genome (beginning,end) ')
        self.prompt_range.place (x=725+350,y=125)
        self.insert_range = Entry(main_int,bd=2)
        self.insert_range.place(x=725+350, y=150)
        
        self.show_gene = Button(main_int, text='Add Gene', command = lambda: self.add_gene(self.insert_id.get(), self.insert_layer.get(), self.insert_index.get(), self.insert_range.get() ), bd=2)
        self.show_gene.place(x=725+350, y=205) 
        self.prompt_index = Label(main_int, text='at index (0-indexing) \n in layer shown below')
        self.prompt_index.place (x=725+450,y=200)
        self.insert_index = Entry(main_int,bd=2)
        self.insert_index.place(x=725+575, y=210)
        self.insert_index.insert(0,0)
        
        self.show_layer = Button(main_int, text='Add Layer', command = lambda: self.add_layer(self.insert_layer.get() ), bd=2)
        self.show_layer.place(x=725+350, y=260)
        self.prompt_layer = Label(main_int,text = 'at index (0-indexing) \n in current list ')
        self.prompt_layer.place(x=725+450,y=255)
        self.insert_layer = Entry(main_int, bd=2)
        self.insert_layer.place(x=725+575,y=260)
        self.insert_layer.insert(0,0)        

        self.remove_entry = Button(main_int, text='Remove Entry', command = lambda: self.remove_gene(self.insert_layer.get(), self.insert_removal_index.get() ), bd=2)
        self.remove_entry.place(x=725+350, y=305) 
        self.prompt_index_removal = Label(main_int, text='at index (0-indexing) \n in layer shown above')
        self.prompt_index_removal.place (x=725+450,y=300)
        self.insert_index_removal = Entry(main_int,bd=2)
        self.insert_index_removal.place(x=725+575, y=305)

        self.remove_layer = Button(main_int, text = 'Delete Layer', command = lambda: self.delete_layer(self.insert_layer_removal.get() ), bd=2)
        self.remove_layer.place(x=725+350, y=352)
        self.prompt_layer_removal = Label(main_int, text = 'at index (0-index) \n in current list')
        self.prompt_layer_removal.place(x=725+450, y=345)
        self.insert_layer_removal = Entry (main_int, bd=2)
        self.insert_layer_removal.place(x=725+575, y=352)

        self.gene_list_title = Label(main_int, text = 'Genes to be added in sequential order:')
        self.gene_list_title.place(x=25,y=25)

        self.gene_list_disp = StringVar()
        self.gene_display = Label(main_int, textvariable=self.gene_list_disp)
        self.gene_display.place(x=50,y=50)

        self.floro_break = Button(main_int,text = 'Insert fluorescence cartridge',command = lambda: self.insert_fluoro(self.insert_layer.get(), self.insert_index.get() ), bd=2)
        self.floro_break.place(x=725+350, y=420)

        self.set_config = Button(main_int, text='Open Settings', command=self.get_config, bd=2)
        self.set_config.place(x=725+450,y=475)

        self.run_tool = Button(main_int, text ='Run tool', command = self.run_master)
        self.run_tool.place(x=725+462,y=525)

    
        
    def add_gene(self,gene_ID,index, subindex, base_range):
        index = int(index)
        subindex = int(subindex)
        gene_name = get_gene_name(gene_ID)
        gene_IDs[index].insert(subindex,gene_ID)
        gene_ranges[index].insert(subindex,base_range)
        gene_list[index].insert(subindex,gene_name)
        self.gene_list_disp.set('\n'.join(', '.join(map(str,sl) ) for sl in gene_list) ) 
        self.insert_index.delete(0,END)
        self.insert_index.insert(0,subindex+1)


    def add_layer(self, index):
        index = int(index)
        gene_list.insert(index,[])
        gene_IDs.insert(index,[])
        gene_ranges.insert(index,[])
        gene_sequences.insert(index,[])
        self.gene_list_disp.set('\n'.join(', '.join(map(str,sl) ) for sl in gene_list) ) 
        self.insert_index.delete(0,END)
        self.insert_index.insert(0,0)


    def insert_fluoro(self,index,subindex):
        index=int(index)
        subindex=int(subindex)
        gene_IDs[index].insert(subindex,str('0000000'))
        gene_ranges[index].insert(subindex,(0,0))
        gene_list[index].insert(subindex,'Fluorescent Reporter Cartridge')
        self.gene_list_disp.set('\n'.join(', '.join(map(str,sl) ) for sl in gene_list) )
        self.insert_index.delete(0,END)
        self.insert_index.insert(0,subindex+1)


    def remove_gene(self,index,subindex):
        index = int(index)
        subindex = int(subindex)
        del gene_IDs[index][subindex]
        del gene_ranges[index][subindex]
        del gene_list[index][subindex]
        self.gene_list_disp.set('\n'.join(', '.join(map(str,sl) ) for sl in gene_list) )
    

    def delete_layer(self,index):
        index = int(index)
        del gene_IDs[index]
        del gene_ranges[index]
        del gene_list[index]
        self.gene_list_disp.set('\n'.join(', '.join(map(str,sl) ) for sl in gene_list) )

    def get_config(self):
        global config_int
        config_int=Tk()
        config_int.geometry('512x512')
        init_config(config_int)
        config_int.mainloop()
        
    def run_master(self):
        gene_sequences = get_gene_sequences(gene_IDs, gene_ranges)


global gene_list
gene_list = [[]]

global gene_IDs
gene_IDs = [[]]

global gene_ranges
gene_ranges = [[]]

global gene_sequences
gene_sequences = [[]]


param_config = {'Chromosome': 1, 'Homology Length': 25, 'Cut Range': '(34,60)', 'Polymerase': 'Phusion', 'Path': 'C:\\Users\\C\\Cindustries\\Engineering Work & History\\BU Classes\\BE 552\\Project'}


main_int = Tk()
TWIZLR_int(main_int)
main_int.geometry('1440x768')
main_int.mainloop()    

chromo_sequence = get_vnat_genome(param_config['Chromosome'])
homology_arms = get_homology_arms(param_config['Chromosome'],chromo_sequence,eval(param_config['Cut Range']),param_config['Homology Length'])        