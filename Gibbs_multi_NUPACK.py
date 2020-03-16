# -*- coding: utf-8 -*-
import numpy as np
import nupack

def Gibbs(probes, targets, send_end):
    G_list = np.zeros((np.shape(probes)[0], np.shape(targets)[0]))
    for p in range(np.shape(probes)[0]):
        for t in range(np.shape(targets)[0]):
            #print(p,t)
            s1 = probes[p,:]
            s2 = targets[t,:]

            Seq1, Seq2 = translate(s1,s2)
            rna_seqs = [Seq1[0], Seq2[0]]

            G = float(nupack.mfe(rna_seqs, ordering=None, material='rna', dangles='some', T=37, multi=0, pseudo=False, sodium=1.0, magnesium=0.0, degenerate=False)[0][1])
            #print(nupack.mfe(rna_seqs, ordering=None, material='rna', dangles='some', T=37, multi=0, pseudo=False, sodium=1.0, magnesium=0.0, degenerate=False))

            G_list[p,t] = G
    #print(G_list)
    send_end.send(G_list)

    return G_list


def translate(seq1, seq2):
    Seq1 = []
    Seq2 = []
    for i in range(np.shape(seq1)[0]):
        if seq1[i] == 1:
            Seq1.append("A")
        elif seq1[i] == 2:
            Seq1.append("U")
        elif seq1[i] == 3:
            Seq1.append("G")
        elif seq1[i] == 4:
            Seq1.append("C")

        Seq1 = ["".join(Seq1)]

        if seq2[i] == 1:
            Seq2.append("A")
        elif seq2[i] == 2:
            Seq2.append("U")
        elif seq2[i] == 3:
            Seq2.append("G")
        elif seq2[i] == 4:
            Seq2.append("C")

        Seq2 = ["".join(Seq2)]

    #print Seq1, Seq2

    return Seq1, Seq2
