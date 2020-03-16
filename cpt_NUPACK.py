# -*- coding: utf-8 -*-

import numpy as np
from scipy import optimize

def Conc(Keqs, Khp_target, Khp_probe, TargetsConc, ProbeConc):

    Cp = np.zeros((np.shape(Khp_probe)[0],np.shape(TargetsConc)[1]))
    Ct = np.zeros((np.shape(TargetsConc)[1], np.shape(Khp_probe)[0], np.shape(Khp_target)[1]))

    Cp_hp = np.zeros((np.shape(Khp_probe)[0],np.shape(TargetsConc)[1]))
    Ct_hp = np.zeros((np.shape(TargetsConc)[1], np.shape(Khp_probe)[0], np.shape(Khp_target)[1]))

    for cell in range(np.shape(TargetsConc)[1]):
        for i in range(np.shape(Keqs)[0]):
            def f(x):
                y = 0
                for j in range(np.shape(Keqs)[1]):
                    #print(i,j, Keqs[i,j], TargetsConc[j,cell],Khp_target[0,j])
                    sigma = (Keqs[i, j] * TargetsConc[j, cell]) / (1 + (Keqs[i, j] * x) + Khp_target[0, j])
                    y += sigma
                return ((x * (1 + y + Khp_probe[i,0])) - ProbeConc)

            #뉴튼법 오류 극복 해야함. 초기값에 따라 결과 경항도 달라짐.
            Cp[i,cell] = optimize.newton(f, x0= ProbeConc, tol= 0.0000000000000000001)

            for k in range(np.shape(Keqs)[1]):
                Ct[cell, i, k] = TargetsConc[k, cell] / (1 + (Keqs[i, k] * Cp[i, cell]) + Khp_target[0,k])
                Ct_hp[cell, i, k] = Ct[cell, i, k] * Khp_target[0,k]
                Cp_hp[i,cell] = Cp[i,cell] * Khp_probe[i,0]

    return Cp, Ct, Ct_hp, Cp_hp