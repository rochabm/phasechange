# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 18:16:06 2017
@author: gessica
"""

import solucao_huang as sol
import matplotlib.pyplot as plt
import numpy as np
import os, sys, shutil

def Position(x,T,Tm):
    n = len(x)
    for i in range (1,n-1):
        if(T[i+1] != T[i]):
            qsi = 2.0*((Tm - T[i])/(T[i+1] - T[i])) -1.0
            print('%d' %i)
            print(qsi)
            if (qsi >= -1.0 and qsi <= 1.0):
                print('Entrou2 \n')
                a =  (T[i] - T[i-1])/(x[i]-x[i-1])
                b = T[i-1] - a*x[i-1]
                print('a = %f' %a)
                print('b = %f' %b)
                #y = a*x + b
                x_Inter = lambda y: (y-b)/a
                return x_Inter(Tm)
    return None

Tm = -0.15
dt = 0.01

tt = np.arange(0.0,4.0+dt,dt)
xIn = np.zeros(len(tt))

# solucao analitica da frente de solidificacao
st,TT,front = sol.solucao(4.0)

# solucao analitica da distribuicao de temperatura (no espaco)
xIn[0] = 0.0
for t in range(1,len(tt)):
    x,T = sol.solucaoEspaco(t*dt,24,4.0)
    xIn[t] = Position(x,T,Tm)

# plota graficos
plt.plot(tt,xIn,'.',label="Exato (dist. temperatura)")
plt.plot(st,front,'-',label="Exato")
plt.xlabel(u"Tempo(s)",fontsize=11)
plt.ylabel("Posicao da Interface",fontsize=11)
plt.legend(loc="best")
plt.grid()
