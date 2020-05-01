def revcomplemen(seq):
    if isinstance(seq,str):
        seq=list(seq)
    revseq=seq[::-1]
    revcom=[0]*len(seq)
    for s in range(len(seq)):
        if revseq[s]=="C":
            revcom[s]="G"
        elif revseq[s]=="G":
            revcom[s]="C"
        elif revseq[s]=="A":
            revcom[s]="T"
        elif revseq[s]=="T":
                revcom[s]="A"
    return "".join(revcom)

def complemen(seq):
    com=[0]*len(seq)
    for s in range(len(seq)):
        if seq[s]=="C":
            com[s]="G"
        elif seq[s]=="G":
            com[s]="C"
        elif seq[s]=="A":
            com[s]="T"
        elif seq[s]=="T":
                com[s]="A"
    return "".join(com)