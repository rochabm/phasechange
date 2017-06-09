# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 14:33:03 2016
@author: gessica
Problema de solidificacao usando metodo da CapEff.
"""

import matplotlib.pyplot as plt
import numpy as np
import solucao_vila_real as sol
from math import sqrt
plt.style.use('classic')  # comando pra usar a configuração antiga do matplotlib

def Erro_Relativo(T_n,T_a,h):
    dif = np.array(T_n - T_a)
    aux = 0.0
    aux2 = 0.0
    for i in range (len(T_n)):
        aux = aux + dif[i]*dif[i]
        aux2 = aux2 + T_a[i]*T_a[i]
    Erro = (sqrt(h*aux)/sqrt(aux2))
    return Erro

tfim = 2.0
tam = int(tfim)*10 + 1
st,T,front = sol.solucao(tfim)
# Esquema para selecionar 41 ou 21 pontos da Temperatura , tempo e frontPosition da solução analítica
dt = tfim/(len(st)-1)
tr = 0.1
s = tr/dt
sT = np.zeros(tam)
sx = np.zeros(tam)
sf = np.zeros(tam)
i=0
for t in range (len(st)):
    if (t % s == 0):
        sT[i] = T[t]
        sx[i] = st[t]
        sf[i] = front[t]
        i=i+1
#------------------------------------------------------------------------------

P1 = np.loadtxt('saida/dados_temp_x.txt')
T_g1  = P1[:,1]
x_g1 = P1[:,0]
h = x_g1[1] - x_g1[0]


I1 = np.loadtxt('saida/Interface.txt')
xIn = I1[:,1]
xt = I1[:,0]


Q1 = np.loadtxt('saida/dados_pontoX1.txt')
T_1  = Q1[:,1]
x_1  = Q1[:,0]


#plt.plot(x,T,'-',label=u"Sol. Analítica")
x,TT = sol.solucaoEspaco(tfim,len(x_g1)-1,4.0)
#----------------------------------------------------------------------
erro1 = Erro_Relativo(T_g1,TT,h)
erro2 = Erro_Relativo(xIn,sf,h)
erro3 = Erro_Relativo(T_1,sT,h)
print('Erro relativo da Figura 1 - Temperatura x espaco: %e\n ' %erro1)
print('Erro relativo da Figura 2 - Pos_Interface x tempo: %e\n ' %erro2)
print('Erro relativo da Figura 3 - Temperatura x tempo: %e\n ' %erro3)



#----------------------------------------------------------------------

plt.figure(1) # Solucao no espaco
plt.plot(x,TT,'-',label="Exato")
plt.plot(x_g1,T_g1,'.-',label=r"Numerico")
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
