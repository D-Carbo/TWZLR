import local #contains a function for local alignment and traceback
import semiglobal #contains a function for semi global alignment
import primer3
import revcomplement as rev


def get_primers(homology1, homology2, insert):
    
    Slocal = [
        # A   G   C   T
        [ 3, -3, -6, -6],  # A
        [-3,  3, -6, -6],  # G
        [-6, -6,  3, -3],  # C
        [-6, -6, -3,  3]   # T
    ]
    gap_penalty_local = 4


    homology1 = list(homology1)
    homology2 = list(homology2)
    insert = list (insert)
    tm=[]
    tm.append(local.localtm(homology1, homology2, Slocal, gap_penalty_local))
    tm.append(local.localtm(homology1, insert, Slocal, gap_penalty_local))
    tm.append(local.localtm(insert, homology2, Slocal, gap_penalty_local))
    seqstm=max(tm)

    Ssemi = [
        # A   G   C   T
        [ 3, -1, -2, -2],  # A
        [-1,  3, -2, -2],  # G
        [-2, -2,  3, -1],  # C
        [-2, -2, -1,  3]   # T
    ]
    gap_penalty_semi = 30


    #design insert primers, part1

    primermin=10
    primermax=40
    forwards=[0]*(primermax-primermin)
    reverses=[0]*(primermax-primermin)
    for n in range(primermin,primermax):
        forwards[n-primermin]=insert[:n];
    #List of potential forward and reverse primers in increasing size
    #make reverse complement of sequence and start doing alignments
    revcom=rev.revcomplemen(insert)
    for n in range(primermin,primermax):
        reverses[n-primermin]=revcom[:n];

    #semi global alignmment to forward and reverse strands
    cutoff=20
    for p in range(0,len(forwards)):
        primerseq=list(forwards[p])
        score, F, TB = semiglobal.semiglobal(insert, primerseq, Ssemi, gap_penalty_semi)
        #max score should be at F[len(seq2)][len[seq2]]
        scores=[0]*(len(F)-len(primerseq))
        for n in range(len(primerseq),len(F)):
            scores[n-len(primerseq)]=F[n][len(primerseq)]
        fullscore=scores[0]#get top score
        secondscore=max(scores[1::])
        #print(fullscore-secondscore)
        if fullscore-secondscore>cutoff:
            break
    potential_fwd="".join(forwards[p])

    for q in range(0,len(reverses)):
        primerseq=list(reverses[q])
        score, F, TB = semiglobal.semiglobal(revcom, primerseq, Ssemi, gap_penalty_semi)
        #max score should be at F[len(seq2)][len[seq2]]
        scores=[0]*(len(F)-len(primerseq))
        for n in range(len(primerseq),len(F)):
            scores[n-len(primerseq)]=F[n][len(primerseq)]
        fullscore=scores[0]#get top score
        secondscore=max(scores[1::])
        if fullscore-secondscore>cutoff:
            break
    potential_rev="".join(reverses[q])


    #calculate Tm of the two primers. if it is not close, add another base to cold one
    #use pip install primer3-py
    tf=primer3.bindings.calcTm(potential_fwd)
    tr=primer3.bindings.calcTm(potential_rev)
    tr=50
    while tr > (tf+4) or tf> (tr+4) or tr<53 or tf<53:
        if tr>tf:
            p=p+1
            potential_fwd=''.join(forwards[p])
            tf=primer3.bindings.calcTm(potential_fwd)
        else:
            q=q+1
            potential_rev=''.join(reverses[q])
            tr=primer3.bindings.calcTm(potential_rev)
            
    #design homology primers
            
    #default 
    h1f1=0
    h1f2=5
    h1r1=0
    h1r2=5      

    h2f1=0
    h2f2=5
    h2r1=0
    h2r2=5           

    #design the external primers h1f1 and h2r2 first
    #they need to be able to amplify the whole assembly
    assembledpart="".join(homology1+insert+homology2)

    primermin=10
    primermax=30
    h1fwds = [[0 for x in range(h1f2-h1f1)] for y in range(primermax-primermin)] 
    h2revs= [[0 for x in range(h1f2-h1f1)] for y in range(primermax-primermin)] 
    for m in range(h1f1,h1f2):
        for n in range(primermin,primermax):
            h1fwds[n-primermin][m]=homology1[m:m+n];
    #List of potential forward and reverse primers in increasing size
    #make reverse complement of sequence and start doing alignments
    revcom=rev.revcomplemen(assembledpart)
    for m in range(h2r1,h2r2):
        for n in range(primermin,primermax):
            h2revs[n-primermin][m]=revcom[m:m+n];

    #semi global alignmment to forward and reverse strands
    cutoff=20
    for p in range(0,len(h1fwds)):
        for p2 in range(0,len(h1fwds[1])):
            primerseq=list(h1fwds[p][p2])
            score, F, TB = semiglobal.semiglobal(assembledpart, primerseq, Ssemi, gap_penalty_semi)
            #max score should be at F[len(seq2)][len[seq2]]
            scores=[0]*(len(F)-len(primerseq))
            for n in range(len(primerseq),len(F)):
                scores[n-len(primerseq)]=F[n][len(primerseq)]
            fullscore=scores[0]#get top score
            secondscore=max(scores[1::])
            if fullscore-secondscore>cutoff:
                break
        #print(fullscore-secondscore)
        if fullscore-secondscore>cutoff:
            break
    potential_h1fwd="".join(h1fwds[p][p2])

    for q in range(0,len(h2revs)):
        for q2 in range(0,len(h2revs[1])):
            primerseq=list(h2revs[q][q2])
            score, F, TB = semiglobal.semiglobal(revcom, primerseq, Ssemi, gap_penalty_semi)
            #max score should be at F[len(seq2)][len[seq2]]
            scores=[0]*(len(F)-len(primerseq))
            for n in range(len(primerseq),len(F)):
                scores[n-len(primerseq)]=F[n][len(primerseq)]
            fullscore=scores[0]#get top score
            secondscore=max(scores[1::])
            if fullscore-secondscore>cutoff:
                break
        if fullscore-secondscore>cutoff:
            break
    potential_h2rev="".join(h2revs[q][q2])

    #calculate Tm of the two primers. if it is not close, add another base to cold one
    #use pip install primer3-py
    h1tf=primer3.bindings.calcTm(potential_h1fwd)
    h2tr=primer3.bindings.calcTm(potential_h2rev)
    while h2tr > (h1tf+4) or h1tf> (h2tr+4) or h2tr<53 or h1tf<53:
        if h2tr>h1tf:
            p=p+1
            potential_h1fwd=''.join(h1fwds[p][p2])
            h1tf=primer3.bindings.calcTm(potential_h1fwd)
        else:
            q=q+1
            potential_h2rev=''.join(h2revs[q][q2])
            h2tr=primer3.bindings.calcTm(potential_h2rev)
        
    #now h1tf and h2tr are fixed, so get values for h1tr and then h2tf
    #we can start with h1tr:
    primermin=10
    primermax=30
    h1revs = [[0 for x in range(h1f2-h1f1)] for y in range(primermax-primermin)] 
    revcom=rev.revcomplemen(homology1)
    for m in range(h1r1,h1r2):
        for n in range(primermin,primermax):
            h1revs[n-primermin][m]=revcom[m:m+n];
    #semi global alignmment to forward and reverse strands
    cutoff=20
    for q in range(0,len(h1revs)):
        for q2 in range(0,len(h1revs[1])):
            primerseq=list(h1revs[q][q2])
            score, F, TB = semiglobal.semiglobal(revcom, primerseq, Ssemi, gap_penalty_semi)
            #max score should be at F[len(seq2)][len[seq2]]
            scores=[0]*(len(F)-len(primerseq))
            for n in range(len(primerseq),len(F)):
                scores[n-len(primerseq)]=F[n][len(primerseq)]
            fullscore=scores[0]#get top score
            secondscore=max(scores[1::])
            if fullscore-secondscore>cutoff:
                break
        if fullscore-secondscore>cutoff:
            break
    potential_h1rev="".join(h1revs[q][q2])
    h1tr=primer3.bindings.calcTm(potential_h1rev)
    while h1tr<(h1tf):
        q=q+1
        potential_h1rev="".join(h1revs[q][q2])
        h1tr=primer3.bindings.calcTm(potential_h1rev)
        
        
    #now for h2tf:
    primermin=10
    primermax=30
    h2fwds = [[0 for x in range(h1f2-h1f1)] for y in range(primermax-primermin)] 
    for m in range(h2f1,h2f2):
        for n in range(primermin,primermax):
            h2fwds[n-primermin][m]=homology2[m:m+n];
    cutoff=20
    for p in range(0,len(h2fwds)):
        for p2 in range(0,len(h2fwds[1])):
            primerseq=list(h2fwds[p][p2])
            score, F, TB = semiglobal.semiglobal(homology2, primerseq, Ssemi, gap_penalty_semi)
            #max score should be at F[len(seq2)][len[seq2]]
            scores=[0]*(len(F)-len(primerseq))
            for n in range(len(primerseq),len(F)):
                scores[n-len(primerseq)]=F[n][len(primerseq)]
            fullscore=scores[0]#get top score
            secondscore=max(scores[1::])
            if fullscore-secondscore>cutoff:
                break
        #print(fullscore-secondscore)
        if fullscore-secondscore>cutoff:
            break
    potential_h2fwd="".join(h2fwds[p][p2])

    h2tf=primer3.bindings.calcTm(potential_h2fwd)
    while h2tf<(h2tr):
        p=p+1
        potential_h2fwd="".join(h2fwds[p][p2])
        h2tf=primer3.bindings.calcTm(potential_h2fwd)



    #finally calculate potential ovarlaps
    primer_sites=[potential_h1fwd, potential_h1rev, potential_fwd, potential_rev, potential_h2fwd, potential_h2rev]
    #rev primers are typed as their reverse complements
    #stitch together some primers to build initial overlaps, trim overlaps
    overlap1=rev.revcomplemen(potential_h1rev)+(potential_fwd)
    overlap2=rev.revcomplemen(potential_rev)+(potential_h2fwd)
    #recall seqstm
    temps=[h1tf,h1tr,tf,tr,h2tf,h2tr,seqstm]
    #trim overlaps to that overlap tm exceeds othet tms
    ol1tm=primer3.bindings.calcTm(overlap1)
    ol2tm=primer3.bindings.calcTm(overlap2)
    while ol1tm>(h1tf+10) and ol1tm>(h2tr+10) and ol1tm>(seqstm+10):
        overlap1=overlap1[1:]#trim away overlap to homology1
        ol1tm=primer3.bindings.calcTm(overlap1)
    while ol2tm>(h1tf+10) and ol2tm>(h2tr+10) and ol2tm>(seqstm+10):
        overlap1=overlap1[:-1]#trim away overlap to homology2
        ol2tm=primer3.bindings.calcTm(overlap1)

    #get final primers
    final_h1f=potential_h1fwd
    final_h2r=potential_h2rev
    final_fwd=overlap1
    final_rev=rev.revcomplemen(overlap2)
    final_h2f=rev.revcomplemen(potential_rev)+potential_h2fwd#full overlap2
    final_h1r=rev.revcomplemen(potential_fwd)+potential_h1rev #reverse complement of full overlap1
    primer_list=[final_h1f,final_h1r,final_fwd,final_rev,final_h2f,final_h2r]
    return primer_list
