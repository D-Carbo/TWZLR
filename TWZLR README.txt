TWZLR README

TWZLR dependencies (all of which can be installed with pip):
- Biopython
- networkx
- matplotlib
- primer3

All other dependencies may be imported from github repo directory

To use TWZLR:
0. pip install above dependencies if not already on machine
1. clone the repo into a directory on your local machine
2. change the working python directory to the cloned directory
3. run TWZLR.py

Settings Configuration:

Chromosome:  user selects Vibrio natriegen chromosome (1 or 2) to be modified (this parameter is critical and the tool will not run without a selection) DEFAULT: N/A
Insertion range: (this parameter is critical and tool will not run without proper input)                                                                 DEFAULT: N/A
    user may enter a single integer as an insertion point where the entire construct of combinations will be inserted after
    user may enter a tuple (format: (beginning,end) ) as the range to be cut from the chromosome and the construct inserted
Homology Length: user enters single integer to use as length for homology arms from chromosome for primer design                                         DEFAULT:500
Polymerase: user selects polymerase from list to use in assembly method                                                                                  DEFAULT: Phusion 
Linkers: user opts to use linkers in the construct                                                                                                       DEFAULT:off
Output path: file dialog box will open to prompt user to select the output for the data from the tool                                                    DEFAULT: current working directory
Save settings: confirm setting configuration and proceed to main window



TWZLR Main Interface:

Insert Gene ID entry: GENE ID from https://www.ncbi.nlm.nih.gov/gene/ for gene to be inserted
Insert base range entry: base range for the gene from parent genome for gene to be inserted (also found at https://www.ncbi.nlm.nih.gov/gene/) (if types tuple or list are not entered for each gene an error will be thrown)
Add gene: insert the gene with entered ID and base range at the subindex shown immediately to the right in current layer
Add layer: insert layer at index shown immediately to the right of button, (index shown is also the working layer for gene insertion)
Remove entry: remove gene at subindex specified immediately to the right in current layer
Delete layer: delete entire layer at index specified immediately to the right
Open settings: allows for user to reopen settings interface to modify parameters (if parameter window is reaccessed, all settings must be reconfigurated to desired value before pressing 'save settings')
Run tool: run tool for selected gene array, generate list of primers and combinatorial graph


Output of tool:
Upon running, user will see combinatorial graph of construct displayed

To acces results: 
Navigate to selected output directory to find TWZLR_primers.json file and open it with preferred json reader

Results format: 
The results are outputted in a json file of format
{
    "Layer n": {
            "component to the left- component - component to the right": [
                "CAGGCCTGCCAGGCAA",
                "AGGCGCGGGGGAGATAAGAGTGCCGGATCTTTTTCCTGA",
                "GCACTCTTATCTCCCCCGCGCC",
                "CAGTGAAAAGTTCTTCTCCTTTACTCATTTGGTTCGACCCCATTTGCC",
                "GGCAAATGGGGTCGAACCAAATGAGTAAAGGAGAAGAACTTTTCACTG",
                "TCAACAAGAATTGGGACAACTCC"
                ]
            }
}

for every possible combination of components between adjacent layers and the DNA base sequences shown are the primers