import revcomplement as rev
import primer3

base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE, PTR_RESET = 0, 1, 2, 3, 4

def local(seq1, seq2, S, gap_penalty):
    """
    Return the score of the optimal semiglobal alignment for seq1
    and seq2.
    Note: gap_penalty should be positive (it is subtracted)
    """
    F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

    # initialize dynamic programming table for local alignment
    for i in range(1, len(seq1)+1):
        F[i][0] = 0 #local: start and end gap penalties not penalized on either sequence
        TB[i][0] = PTR_GAP2  # indicates a gap in seq2
    for j in range(1, len(seq2)+1):
        F[0][j] = 0 #local: start and end gap penalties not penalized on either sequence
        TB[0][j] = PTR_GAP1  # indicates a gap in seq1

    
    #Make F
    for i in range(1,len(seq1)+1):
        # all rows except for 0
        for j in range(1,len(seq2)+1):
            #all columns except 0
            #note seq1[0] is the first base of seq 1 
            #so we use seq1[i-1] for 0 to its full length
            F[i][j]=max((F[i-1][j-1]+S[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]),
             (F[i-1][j]-gap_penalty), (F[i][j-1]-gap_penalty))
            if F[i][j]<0:
                F[i][j]=0
            
    #Make TB
    #0 row/column of TB has already been set up
    #First, assign TB values for when sequences match
    for i in range(1,len(seq1)+1):
        # all rows except for 0
        for j in range(1,len(seq2)+1):
            #Match or mismatch w/o gap
            if F[i][j] == F[i-1][j-1]+S[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]] and F[i][j]>0:
                TB[i][j]=PTR_BASE
            #Gap in seq 2
            elif F[i][j] == F[i-1][j]-gap_penalty and F[i][j]>0:
                TB[i][j]=PTR_GAP2
            #Gap in seq 1
            elif F[i][j] == F[i][j-1]-gap_penalty and F[i][j]>0:
                TB[i][j]=PTR_GAP1
            elif F[i][j]==0:
                TB[i][j]=PTR_RESET
    return F[len(seq1)][len(seq2)], F, TB



def localtraceback(seq1, seq2, TB, loc1, loc2):
    s1 = ""
    s2 = ""

    i = loc1
    j = loc2

    while TB[i][j] != PTR_NONE:
        #not = 0
        if TB[i][j] == PTR_BASE:
            #=3 (alignment)
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i = i - 1
            j = j - 1
        elif TB[i][j] == PTR_GAP1:
            #=1 (one type of gap)
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j = j - 1
        elif TB[i][j] == PTR_GAP2:
            #=2 (other type of gap)
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i = i - 1
        elif TB[i][j]==PTR_RESET:
            break
        else:
            assert False

    return s1, s2

def localtm(sa,sb,Slocal,gap_penalty_local):
    score, F, TB = local(sa, sb, Slocal, gap_penalty_local)
    #get traceback of max score and calculate Tm
    maxlist=[0]*len(F)
    for l in range(len(F)):
        maxlist[l]=max(F[l])
    maxlocation=maxlist.index(max(maxlist)) #row of F with the best alignment
    maxloc2=F[maxlocation].index(max(F[maxlocation]))
    #get traceback of the best alignment sequence!
    s1, s2 = localtraceback(sa, sb, TB, maxlocation, maxloc2)
    s1list=list(s1)
    s2list=list(s2)
    while ('-') in s1list:
        s1list.remove('-')
    while ('-') in s2list:
        s2list.remove('-')
    s1_nogap="".join(s1list[0:60])#primer3 only goes to 60
    s2_nogap="".join(s2list[0:60])
    
    parts_tm=(primer3.bindings.calcHeterodimer(s1_nogap, rev.revcomplemen(s2_nogap)).tm)
    return parts_tm