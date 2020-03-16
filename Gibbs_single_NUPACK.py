# -*- coding: utf-8 -*-
import numpy as np
import nupack

def Gibbs(probes, send_end):
    G_list = np.zeros((np.shape(probes)[0],1))
    for p in range(np.shape(probes)[0]):
        s1 = probes[p,:]

        Seq = translate(s1)
        rna_seqs = [Seq[0]]

        G = float(nupack.mfe(rna_seqs, ordering=None, material='rna', dangles='some', T=37, multi=0, pseudo=False, sodium=1.0, magnesium=0.0, degenerate=False)[0][1])
        #print(nupack.mfe(rna_seqs, ordering=None, material='rna', dangles='some', T=37, multi=0, pseudo=False, sodium=1.0, magnesium=0.0, degenerate=False))

        G_list[p,0] = G
    #print(G_list)
    send_end.send(G_list)
    #print(G_list)
    return G_list


def translate(seq):
    Seq = []
    for i in range(np.shape(seq)[0]):
        if seq[i] == 1:
            Seq.append("A")
        elif seq[i] == 2:
            Seq.append("U")
        elif seq[i] == 3:
            Seq.append("G")
        elif seq[i] == 4:
            Seq.append("C")

        Seq = ["".join(Seq)]

    return Seq
