base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE, PTR_RESET = 0, 1, 2, 3, 4

def semiglobal(seq1, seq2, S, gap_penalty):
    """
    Return the score of the optimal semiglobal alignment for seq1
    and seq2.
    Note: gap_penalty should be positive (it is subtracted)
    """
    F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

    # initialize dynamic programming table for Needleman-Wunsch alignment
    for i in range(1, len(seq1)+1):
        F[i][0] = 0 #semiglobal: start and end gap penalties not penalized
        TB[i][0] = PTR_GAP2  # indicates a gap in seq2
    for j in range(1, len(seq2)+1):
        F[0][j] = 0 - j*gap_penalty
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
            
    #Make TB
    #0 row/column of TB has already been set up
    #First, assign TB values for when sequences match
    for i in range(1,len(seq1)+1):
        # all rows except for 0
        for j in range(1,len(seq2)+1):
            #Match or mismatch w/o gap
            if F[i][j] == F[i-1][j-1]+S[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]:
                TB[i][j]=PTR_BASE
            #Gap in seq 2
            elif F[i][j] == F[i-1][j]-gap_penalty:
                TB[i][j]=PTR_GAP2
            #Gap in seq 1
            elif F[i][j] == F[i][j-1]-gap_penalty:
                TB[i][j]=PTR_GAP1
    return F[len(seq1)][len(seq2)], F, TB