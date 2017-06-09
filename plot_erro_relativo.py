# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 14:33:03 2016
@author: gessica
Problema de solidificacao usando metodo da Entalpia. VR
"""

import matplotlib.pyplot as plt
import numpy as np
import solucao_huang as sol
import solucao_vila_real as sol
from math import sqrt
plt.style.use('classic')  # comando pra usar a configuração antiga do matplotlib

plot_x = [0.001, 0.0025, 0.004, 0.005, 0.01, 0.025, 0.05]

#--------------------------------------------------------------------------------

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
#--------------------------------------------------------------------------------

def Erro_Relativo(T_n,T_a,h):
    dif = np.array(T_n - T_a)
    aux = 0.0
    aux2 = 0.0
    for i in range (len(T_n)):
        aux = aux + dif[i]*dif[i]
        aux2 = aux2 + T_a[i]*T_a[i]
    Erro = (sqrt(h*aux)/sqrt(aux2))
    return Erro


#================================================================================
# MALHA 24 ELEMENTOS
# Dados para o Erro Relativo da Temperatura x espaco
P1_m24_e = np.loadtxt('VR_ME_dt_0001_m24/dados_temp_x.txt')
P2_m24_e = np.loadtxt('VR_ME_dt_00025_m24/dados_temp_x.txt')
P3_m24_e = np.loadtxt('VR_ME_dt_0004_m24/dados_temp_x.txt')
P4_m24_e = np.loadtxt('VR_ME_dt_0005_m24/dados_temp_x.txt')
P5_m24_e = np.loadtxt('VR_ME_dt_0025_m24/dados_temp_x.txt')
P6_m24_e = np.loadtxt('VR_ME_dt_001_m24/dados_temp_x.txt')
P7_m24_e = np.loadtxt('VR_ME_dt_005_m24/dados_temp_x.txt')

x_g1 = P1_m24_e[:,0]
h = x_g1[1] - x_g1[0]
x,TT = sol.solucaoEspaco(tfim,len(x_g1)-1,4.0)

T_1_m24_e  = P1_m24_e[:,1]
T_2_m24_e  = P2_m24_e[:,1]
T_3_m24_e  = P3_m24_e[:,1]
T_4_m24_e  = P4_m24_e[:,1]
T_5_m24_e  = P5_m24_e[:,1]
T_6_m24_e  = P6_m24_e[:,1]
T_7_m24_e  = P7_m24_e[:,1]

T_m24_e = [T_1_m24_e,T_2_m24_e,T_3_m24_e,T_4_m24_e,T_6_m24_e,T_5_m24_e,T_7_m24_e]

plot_y_m24_e = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m24_e[i]= Erro_Relativo(T_m24_e[i],TT,h)


# Dados para o Erro Relativo da Temperatura x tempo
P1_m24 = np.loadtxt('VR_ME_dt_0001_m24/dados_pontoX1.txt')
P2_m24 = np.loadtxt('VR_ME_dt_00025_m24/dados_pontoX1.txt')
P3_m24 = np.loadtxt('VR_ME_dt_0004_m24/dados_pontoX1.txt')
P4_m24 = np.loadtxt('VR_ME_dt_0005_m24/dados_pontoX1.txt')
P5_m24 = np.loadtxt('VR_ME_dt_0025_m24/dados_pontoX1.txt')
P6_m24 = np.loadtxt('VR_ME_dt_001_m24/dados_pontoX1.txt')
P7_m24 = np.loadtxt('VR_ME_dt_005_m24/dados_pontoX1.txt')

T_1_m24  = P1_m24[:,1]
T_2_m24  = P2_m24[:,1]
T_3_m24  = P3_m24[:,1]
T_4_m24  = P4_m24[:,1]
T_5_m24  = P5_m24[:,1]
T_6_m24  = P6_m24[:,1]
T_7_m24  = P7_m24[:,1]

T_m24 = [T_1_m24,T_2_m24,T_3_m24,T_4_m24,T_6_m24,T_5_m24,T_7_m24]

plot_y_m24 = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m24[i]= Erro_Relativo(T_m24[i],sT,h)

# Dados para o Erro Relativo da Interface x tempo
P1_m24_i = np.loadtxt('VR_ME_dt_0001_m24/Interface.txt')
P2_m24_i = np.loadtxt('VR_ME_dt_00025_m24/Interface.txt')
P3_m24_i = np.loadtxt('VR_ME_dt_0004_m24/Interface.txt')
P4_m24_i = np.loadtxt('VR_ME_dt_0005_m24/Interface.txt')
P5_m24_i = np.loadtxt('VR_ME_dt_0025_m24/Interface.txt')
P6_m24_i = np.loadtxt('VR_ME_dt_001_m24/Interface.txt')
P7_m24_i = np.loadtxt('VR_ME_dt_005_m24/Interface.txt')

xt = P1_m24_i[:,0]

T_1_m24_i  = P1_m24_i[:,1]
T_2_m24_i  = P2_m24_i[:,1]
T_3_m24_i  = P3_m24_i[:,1]
T_4_m24_i  = P4_m24_i[:,1]
T_5_m24_i  = P5_m24_i[:,1]
T_6_m24_i  = P6_m24_i[:,1]
T_7_m24_i  = P7_m24_i[:,1]

T_m24_i = [T_1_m24_i,T_2_m24_i,T_3_m24_i,T_4_m24_i,T_6_m24_i,T_5_m24_i,T_7_m24_i]

plot_y_m24_i = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m24_i[i]= Erro_Relativo(T_m24_i[i],sf,h)

#================================================================================
# MALHA 48 ELEMENTOS
# Dados para o Erro Relativo da Temperatura x espaco
P1_m48_e = np.loadtxt('VR_ME_dt_0001_m48/dados_temp_x.txt')
P2_m48_e = np.loadtxt('VR_ME_dt_00025_m48/dados_temp_x.txt')
P3_m48_e = np.loadtxt('VR_ME_dt_0004_m48/dados_temp_x.txt')
P4_m48_e = np.loadtxt('VR_ME_dt_0005_m48/dados_temp_x.txt')
P5_m48_e = np.loadtxt('VR_ME_dt_0025_m48/dados_temp_x.txt')
P6_m48_e = np.loadtxt('VR_ME_dt_001_m48/dados_temp_x.txt')
P7_m48_e = np.loadtxt('VR_ME_dt_005_m48/dados_temp_x.txt')

x_g1 = P1_m48_e[:,0]
h = x_g1[1] - x_g1[0]
x,TT = sol.solucaoEspaco(tfim,len(x_g1)-1,4.0)

T_1_m48_e  = P1_m48_e[:,1]
T_2_m48_e  = P2_m48_e[:,1]
T_3_m48_e  = P3_m48_e[:,1]
T_4_m48_e  = P4_m48_e[:,1]
T_5_m48_e  = P5_m48_e[:,1]
T_6_m48_e  = P6_m48_e[:,1]
T_7_m48_e  = P7_m48_e[:,1]

T_m48_e = [T_1_m48_e,T_2_m48_e,T_3_m48_e,T_4_m48_e,T_6_m48_e,T_5_m48_e,T_7_m48_e]

plot_y_m48_e = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m48_e[i]= Erro_Relativo(T_m48_e[i],TT,h)


# Dados para o Erro Relativo da Temperatura x tempo
P1_m48 = np.loadtxt('VR_ME_dt_0001_m48/dados_pontoX1.txt')
P2_m48 = np.loadtxt('VR_ME_dt_00025_m48/dados_pontoX1.txt')
P3_m48 = np.loadtxt('VR_ME_dt_0004_m48/dados_pontoX1.txt')
P4_m48 = np.loadtxt('VR_ME_dt_0005_m48/dados_pontoX1.txt')
P5_m48 = np.loadtxt('VR_ME_dt_0025_m48/dados_pontoX1.txt')
P6_m48 = np.loadtxt('VR_ME_dt_001_m48/dados_pontoX1.txt')
P7_m48 = np.loadtxt('VR_ME_dt_005_m48/dados_pontoX1.txt')

T_1_m48  = P1_m48[:,1]
T_2_m48  = P2_m48[:,1]
T_3_m48  = P3_m48[:,1]
T_4_m48  = P4_m48[:,1]
T_5_m48  = P5_m48[:,1]
T_6_m48  = P6_m48[:,1]
T_7_m48  = P7_m48[:,1]

T_m48 = [T_1_m48,T_2_m48,T_3_m48,T_4_m48,T_6_m48,T_5_m48,T_7_m48]

plot_y_m48 = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m48[i]= Erro_Relativo(T_m48[i],sT,h)

# Dados para o Erro Relativo da Interface x tempo
P1_m48_i = np.loadtxt('VR_ME_dt_0001_m48/Interface.txt')
P2_m48_i = np.loadtxt('VR_ME_dt_00025_m48/Interface.txt')
P3_m48_i = np.loadtxt('VR_ME_dt_0004_m48/Interface.txt')
P4_m48_i = np.loadtxt('VR_ME_dt_0005_m48/Interface.txt')
P5_m48_i = np.loadtxt('VR_ME_dt_0025_m48/Interface.txt')
P6_m48_i = np.loadtxt('VR_ME_dt_001_m48/Interface.txt')
P7_m48_i = np.loadtxt('VR_ME_dt_005_m48/Interface.txt')

xt = P1_m48_i[:,0]

T_1_m48_i  = P1_m48_i[:,1]
T_2_m48_i  = P2_m48_i[:,1]
T_3_m48_i  = P3_m48_i[:,1]
T_4_m48_i  = P4_m48_i[:,1]
T_5_m48_i  = P5_m48_i[:,1]
T_6_m48_i  = P6_m48_i[:,1]
T_7_m48_i  = P7_m48_i[:,1]

T_m48_i = [T_1_m48_i,T_2_m48_i,T_3_m48_i,T_4_m48_i,T_6_m48_i,T_5_m48_i,T_7_m48_i]

plot_y_m48_i = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m48_i[i]= Erro_Relativo(T_m48_i[i],sf,h)

#================================================================================
# MALHA 96 ELEMENTOS
# Dados para o Erro Relativo da Temperatura x espaco
P1_m96_e = np.loadtxt('VR_ME_dt_0001_m96/dados_temp_x.txt')
P2_m96_e = np.loadtxt('VR_ME_dt_00025_m96/dados_temp_x.txt')
P3_m96_e = np.loadtxt('VR_ME_dt_0004_m96/dados_temp_x.txt')
P4_m96_e = np.loadtxt('VR_ME_dt_0005_m96/dados_temp_x.txt')
P5_m96_e = np.loadtxt('VR_ME_dt_0025_m96/dados_temp_x.txt')
P6_m96_e = np.loadtxt('VR_ME_dt_001_m96/dados_temp_x.txt')
P7_m96_e = np.loadtxt('VR_ME_dt_005_m96/dados_temp_x.txt')

x_g1 = P1_m96_e[:,0]
h = x_g1[1] - x_g1[0]
x,TT = sol.solucaoEspaco(tfim,len(x_g1)-1,4.0)

T_1_m96_e  = P1_m96_e[:,1]
T_2_m96_e  = P2_m96_e[:,1]
T_3_m96_e  = P3_m96_e[:,1]
T_4_m96_e  = P4_m96_e[:,1]
T_5_m96_e  = P5_m96_e[:,1]
T_6_m96_e  = P6_m96_e[:,1]
T_7_m96_e  = P7_m96_e[:,1]

T_m96_e = [T_1_m96_e,T_2_m96_e,T_3_m96_e,T_4_m96_e,T_6_m96_e,T_5_m96_e,T_7_m96_e]

plot_y_m96_e = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m96_e[i]= Erro_Relativo(T_m96_e[i],TT,h)


# Dados para o Erro Relativo da Temperatura x tempo
P1_m96 = np.loadtxt('VR_ME_dt_0001_m96/dados_pontoX1.txt')
P2_m96 = np.loadtxt('VR_ME_dt_00025_m96/dados_pontoX1.txt')
P3_m96 = np.loadtxt('VR_ME_dt_0004_m96/dados_pontoX1.txt')
P4_m96 = np.loadtxt('VR_ME_dt_0005_m96/dados_pontoX1.txt')
P5_m96 = np.loadtxt('VR_ME_dt_0025_m96/dados_pontoX1.txt')
P6_m96 = np.loadtxt('VR_ME_dt_001_m96/dados_pontoX1.txt')
P7_m96 = np.loadtxt('VR_ME_dt_005_m96/dados_pontoX1.txt')

T_1_m96  = P1_m96[:,1]
T_2_m96  = P2_m96[:,1]
T_3_m96  = P3_m96[:,1]
T_4_m96  = P4_m96[:,1]
T_5_m96  = P5_m96[:,1]
T_6_m96  = P6_m96[:,1]
T_7_m96  = P7_m96[:,1]

T_m96 = [T_1_m96,T_2_m96,T_3_m96,T_4_m96,T_6_m96,T_5_m96,T_7_m96]

plot_y_m96 = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m96[i]= Erro_Relativo(T_m96[i],sT,h)

# Dados para o Erro Relativo da Interface x tempo
P1_m96_i = np.loadtxt('VR_ME_dt_0001_m96/Interface.txt')
P2_m96_i = np.loadtxt('VR_ME_dt_00025_m96/Interface.txt')
P3_m96_i = np.loadtxt('VR_ME_dt_0004_m96/Interface.txt')
P4_m96_i = np.loadtxt('VR_ME_dt_0005_m96/Interface.txt')
P5_m96_i = np.loadtxt('VR_ME_dt_0025_m96/Interface.txt')
P6_m96_i = np.loadtxt('VR_ME_dt_001_m96/Interface.txt')
P7_m96_i = np.loadtxt('VR_ME_dt_005_m96/Interface.txt')

xt = P1_m96_i[:,0]

T_1_m96_i  = P1_m96_i[:,1]
T_2_m96_i  = P2_m96_i[:,1]
T_3_m96_i  = P3_m96_i[:,1]
T_4_m96_i  = P4_m96_i[:,1]
T_5_m96_i  = P5_m96_i[:,1]
T_6_m96_i  = P6_m96_i[:,1]
T_7_m96_i  = P7_m96_i[:,1]

T_m96_i = [T_1_m96_i,T_2_m96_i,T_3_m96_i,T_4_m96_i,T_6_m96_i,T_5_m96_i,T_7_m96_i]

plot_y_m96_i = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m96_i[i]= Erro_Relativo(T_m96_i[i],sf,h)

#================================================================================
# MALHA 192 ELEMENTOS
# Dados para o Erro Relativo da Temperatura x espaco
P1_m192_e = np.loadtxt('VR_ME_dt_0001_m192/dados_temp_x.txt')
P2_m192_e = np.loadtxt('VR_ME_dt_00025_m192/dados_temp_x.txt')
P3_m192_e = np.loadtxt('VR_ME_dt_0004_m192/dados_temp_x.txt')
P4_m192_e = np.loadtxt('VR_ME_dt_0005_m192/dados_temp_x.txt')
P5_m192_e = np.loadtxt('VR_ME_dt_0025_m192/dados_temp_x.txt')
P6_m192_e = np.loadtxt('VR_ME_dt_001_m192/dados_temp_x.txt')
P7_m192_e = np.loadtxt('VR_ME_dt_005_m192/dados_temp_x.txt')

x_g1 = P1_m192_e[:,0]
h = x_g1[1] - x_g1[0]
x,TT = sol.solucaoEspaco(tfim,len(x_g1)-1,4.0)

T_1_m192_e  = P1_m192_e[:,1]
T_2_m192_e  = P2_m192_e[:,1]
T_3_m192_e  = P3_m192_e[:,1]
T_4_m192_e  = P4_m192_e[:,1]
T_5_m192_e  = P5_m192_e[:,1]
T_6_m192_e  = P6_m192_e[:,1]
T_7_m192_e  = P7_m192_e[:,1]

T_m192_e = [T_1_m192_e,T_2_m192_e,T_3_m192_e,T_4_m192_e,T_6_m192_e,T_5_m192_e,T_7_m192_e]

plot_y_m192_e = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m192_e[i]= Erro_Relativo(T_m192_e[i],TT,h)


# Dados para o Erro Relativo da Temperatura x tempo
P1_m192 = np.loadtxt('VR_ME_dt_0001_m192/dados_pontoX1.txt')
P2_m192 = np.loadtxt('VR_ME_dt_00025_m192/dados_pontoX1.txt')
P3_m192 = np.loadtxt('VR_ME_dt_0004_m192/dados_pontoX1.txt')
P4_m192 = np.loadtxt('VR_ME_dt_0005_m192/dados_pontoX1.txt')
P5_m192 = np.loadtxt('VR_ME_dt_0025_m192/dados_pontoX1.txt')
P6_m192 = np.loadtxt('VR_ME_dt_001_m192/dados_pontoX1.txt')
P7_m192 = np.loadtxt('VR_ME_dt_005_m192/dados_pontoX1.txt')

T_1_m192  = P1_m192[:,1]
T_2_m192  = P2_m192[:,1]
T_3_m192  = P3_m192[:,1]
T_4_m192  = P4_m192[:,1]
T_5_m192  = P5_m192[:,1]
T_6_m192  = P6_m192[:,1]
T_7_m192  = P7_m192[:,1]

T_m192 = [T_1_m192,T_2_m192,T_3_m192,T_4_m192,T_6_m192,T_5_m192,T_7_m192]

plot_y_m192 = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m192[i]= Erro_Relativo(T_m192[i],sT,h)

# Dados para o Erro Relativo da Interface x tempo
P1_m192_i = np.loadtxt('VR_ME_dt_0001_m192/Interface.txt')
P2_m192_i = np.loadtxt('VR_ME_dt_00025_m192/Interface.txt')
P3_m192_i = np.loadtxt('VR_ME_dt_0004_m192/Interface.txt')
P4_m192_i = np.loadtxt('VR_ME_dt_0005_m192/Interface.txt')
P5_m192_i = np.loadtxt('VR_ME_dt_0025_m192/Interface.txt')
P6_m192_i = np.loadtxt('VR_ME_dt_001_m192/Interface.txt')
P7_m192_i = np.loadtxt('VR_ME_dt_005_m192/Interface.txt')

xt = P1_m192_i[:,0]

T_1_m192_i  = P1_m192_i[:,1]
T_2_m192_i  = P2_m192_i[:,1]
T_3_m192_i  = P3_m192_i[:,1]
T_4_m192_i  = P4_m192_i[:,1]
T_5_m192_i  = P5_m192_i[:,1]
T_6_m192_i  = P6_m192_i[:,1]
T_7_m192_i  = P7_m192_i[:,1]

T_m192_i = [T_1_m192_i,T_2_m192_i,T_3_m192_i,T_4_m192_i,T_6_m192_i,T_5_m192_i,T_7_m192_i]

plot_y_m192_i = np.zeros(len(plot_x))
for i in range(len(plot_x)):
    plot_y_m192_i[i]= Erro_Relativo(T_m192_i[i],sf,h)

#================================================================================

plt.figure(1)
plt.plot(plot_x,plot_y_m24,'-*',label=u"MCE h = 0.167m")
plt.plot(plot_x,plot_y_m48,'-*',label=u"MCE h = 0.083m")
plt.plot(plot_x,plot_y_m96,'-*',label=u"MCE h = 0,042m")
plt.plot(plot_x,plot_y_m192,'-*',label=u"MCE h = 0,021m")
plt.xlabel(u"Tamanho do passo de tempo",fontsize=11)
plt.ylabel(u'Erro Relativo',fontsize=11)
plt.legend(loc='upper right',prop={'size':11})
plt.grid()

plt.figure(2)
plt.plot(plot_x,plot_y_m24_e,'-*',label=u"MCE h = 0.167m")
plt.plot(plot_x,plot_y_m48_e,'-*',label=u"MCE h = 0.083m")
plt.plot(plot_x,plot_y_m96_e,'-*',label=u"MCE h = 0,042m")
plt.plot(plot_x,plot_y_m192_e,'-*',label=u"MCE h = 0,021m")
plt.xlabel(u"Tamanho do passo de tempo",fontsize=11)
plt.ylabel(u'Erro Relativo',fontsize=11)
plt.legend(loc='upper right',prop={'size':11})
plt.grid()

plt.figure(3)
plt.plot(plot_x,plot_y_m24_i,'-*',label=u"MCE h = 0.167m")
plt.plot(plot_x,plot_y_m48_i,'-*',label=u"MCE h = 0.083m")
plt.plot(plot_x,plot_y_m96_i,'-*',label=u"MCE h = 0,042m")
plt.plot(plot_x,plot_y_m192_i,'-*',label=u"MCE h = 0,021m")
plt.xlabel(u"Tamanho do passo de tempo",fontsize=11)
plt.ylabel(u'Erro Relativo',fontsize=11)
plt.legend(loc='upper right',prop={'size':11})
plt.grid()

plt.show()
