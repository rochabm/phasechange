# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 14:33:03 2016
@author: gessica
Problema de solidificacao usando metodo da CapEff.
"""

import matplotlib.pyplot as plt
import numpy as np
import solucao_huang as sol
plt.style.use('classic')  # comando pra usar a configuração antiga do matplotlib

P1 = np.loadtxt('saida/dados_temp_x.txt')
T_g1  = P1[:]
xx = np.linspace(0.0,4.0,25) # malha com 24 elementos

I1 = np.loadtxt('saida/Interface.txt')
xIn = I1[:,1]
xt = I1[:,0]

Q1 = np.loadtxt('saida/dados_pontoX1.txt')
T_1  = Q1[:,1]
x_1  = Q1[:,0]

st,T,front = sol.solucao(4.0)
#plt.plot(x,T,'-',label=u"Sol. Analítica")
x,TT = sol.solucaoEspaco(4.0,24,4.0) # malha com 24 elementos

plt.figure(1) # Solucao no espaco
plt.plot(x,TT,'-',label="Exato")
plt.plot(xx,T_g1,'.-',label=r"Numerico")
plt.xlabel(u'Espaco(m)',fontsize=11)
plt.ylabel('Temperatura(C)',fontsize=11)
plt.legend(loc="best")
plt.grid()

plt.figure(2)
plt.plot(st,front,'-',label="Exato")
plt.plot(xt,xIn,'.-',label="Numerico")
plt.xlabel(u"Tempo(s)",fontsize=11)
plt.ylabel("Posicao da Interface",fontsize=11)
plt.legend(loc="best")
plt.grid()

plt.figure(3)
plt.plot(st,T,'-',label=u"Exato")
plt.plot(x_1,T_1,'.-',label=r"Numerico")
plt.xlabel(u"Tempo(s)",fontsize=11)
plt.ylabel('Temperatura(C)',fontsize=11)
plt.legend(loc="best")
plt.grid()

plt.show()
