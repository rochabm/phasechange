# -*- coding: utf-8 -*-

"""
Created on Thu Jun  1 15:42:46 2017
@author: gessica
"""

import matplotlib.pyplot as plt
import numpy as np
import solucao_vila_real as sol
from math import sqrt
plt.style.use('classic')

#
# Estilo dos graficos
#
import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['legend.fontsize'] = 18
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times'
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode']=True

#-------------------------------------------------------------------------------

def Erro_Relativo(T_n,T_a,h): 
    dif = np.array(T_n - T_a)
    aux = 0.0
    aux2 = 0.0
    for i in range (len(T_n)):
        aux = aux + dif[i]*dif[i]
        aux2 = aux2 + T_a[i]*T_a[i]
    Erro = (sqrt(h*aux)/sqrt(aux2))
    return Erro    

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_x = [0.001, 0.0025, 0.004, 0.005, 0.01, 0.025, 0.05]  
    plot_x1 = ['0001', '00025', '0004', '0005', '001', '0025', '005']  

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

    # MALHA 24 ELEMENTOS
    T_m24_e, T_m24, T_m24_i = [], [], []
    for t in plot_x1:
        P1_m24_e = np.loadtxt('VR_MCE_dt_%s_m24/dados_temp_x.txt' %t )[:,1]
        P1_m24 = np.loadtxt('VR_MCE_dt_%s_m24/dados_pontoX1.txt' %t )[:,1] 
        P1_m24_i = np.loadtxt('VR_MCE_dt_%s_m24/Interface.txt' %t )[:,1]
    
        T_m24_e.append(P1_m24_e)
        T_m24.append(P1_m24)
        T_m24_i.append(P1_m24_i)

    P1_m24_e = np.loadtxt('VR_MCE_dt_0001_m24/dados_temp_x.txt')
    x_g1 = P1_m24_e[:,0]
    h = x_g1[1] - x_g1[0]
    x,TT = sol.solucaoEspaco(tfim,len(x_g1)-1,4.0) 

    plot_y_m24_e = np.zeros(len(plot_x))
    plot_y_m24 = np.zeros(len(plot_x))
    plot_y_m24_i = np.zeros(len(plot_x))
    
    for i in range(len(plot_x)):  
        plot_y_m24_e[i]= Erro_Relativo(T_m24_e[i],TT,h)
        plot_y_m24[i]= Erro_Relativo(T_m24[i],sT,h)
        plot_y_m24_i[i]= Erro_Relativo(T_m24_i[i],sf,h)
        

    # MALHA 48 ELEMENTOS
    T_m48_e, T_m48, T_m48_i = [], [], []
    for t in plot_x1:
        P1_m48_e = np.loadtxt('VR_MCE_dt_%s_m48/dados_temp_x.txt' %t )[:,1]
        P1_m48 = np.loadtxt('VR_MCE_dt_%s_m48/dados_pontoX1.txt' %t )[:,1] 
        P1_m48_i = np.loadtxt('VR_MCE_dt_%s_m48/Interface.txt' %t )[:,1]
    
        T_m48_e.append(P1_m48_e)
        T_m48.append(P1_m48)
        T_m48_i.append(P1_m48_i)

    P1_m48_e = np.loadtxt('VR_MCE_dt_0001_m48/dados_temp_x.txt')
    x_g1 = P1_m48_e[:,0]
    h = x_g1[1] - x_g1[0]
    x,TT = sol.solucaoEspaco(tfim,len(x_g1)-1,4.0) 

    plot_y_m48_e = np.zeros(len(plot_x))
    plot_y_m48 = np.zeros(len(plot_x))
    plot_y_m48_i = np.zeros(len(plot_x))
    for i in range(len(plot_x)):  
        plot_y_m48_e[i]= Erro_Relativo(T_m48_e[i],TT,h)
        plot_y_m48[i]= Erro_Relativo(T_m48[i],sT,h)
        plot_y_m48_i[i]= Erro_Relativo(T_m48_i[i],sf,h)

    # MALHA 96 ELEMENTOS
    T_m96_e, T_m96, T_m96_i = [], [], []
    for t in plot_x1:
        P1_m96_e = np.loadtxt('VR_MCE_dt_%s_m96/dados_temp_x.txt' %t )[:,1]
        P1_m96 = np.loadtxt('VR_MCE_dt_%s_m96/dados_pontoX1.txt' %t )[:,1] 
        P1_m96_i = np.loadtxt('VR_MCE_dt_%s_m96/Interface.txt' %t )[:,1]
        
        T_m96_e.append(P1_m96_e)
        T_m96.append(P1_m96)
        T_m96_i.append(P1_m96_i)

    P1_m96_e = np.loadtxt('VR_MCE_dt_0001_m96/dados_temp_x.txt')
    x_g1 = P1_m96_e[:,0]
    h = x_g1[1] - x_g1[0]
    x,TT = sol.solucaoEspaco(tfim,len(x_g1)-1,4.0) 

    plot_y_m96_e = np.zeros(len(plot_x))
    plot_y_m96 = np.zeros(len(plot_x))
    plot_y_m96_i = np.zeros(len(plot_x))
    for i in range(len(plot_x)):  
        plot_y_m96_e[i]= Erro_Relativo(T_m96_e[i],TT,h)
        plot_y_m96[i]= Erro_Relativo(T_m96[i],sT,h)
        plot_y_m96_i[i]= Erro_Relativo(T_m96_i[i],sf,h)   

    # MALHA 192 ELEMENTOS
    T_m192_e, T_m192, T_m192_i = [], [], []
    for t in plot_x1:
        P1_m192_e = np.loadtxt('VR_MCE_dt_%s_m192/dados_temp_x.txt' %t )[:,1]
        P1_m192 = np.loadtxt('VR_MCE_dt_%s_m192/dados_pontoX1.txt' %t )[:,1] 
        P1_m192_i = np.loadtxt('VR_MCE_dt_%s_m192/Interface.txt' %t )[:,1]
        T_m192_e.append(P1_m192_e)
        T_m192.append(P1_m192)
        T_m192_i.append(P1_m192_i)

    P1_m192_e = np.loadtxt('VR_MCE_dt_0001_m192/dados_temp_x.txt')
    x_g1 = P1_m192_e[:,0]
    h = x_g1[1] - x_g1[0]
    x,TT = sol.solucaoEspaco(tfim,len(x_g1)-1,4.0) 

    plot_y_m192_e = np.zeros(len(plot_x))
    plot_y_m192 = np.zeros(len(plot_x))
    plot_y_m192_i = np.zeros(len(plot_x))
    for i in range(len(plot_x)):  
        plot_y_m192_e[i]= Erro_Relativo(T_m192_e[i],TT,h)
        plot_y_m192[i]= Erro_Relativo(T_m192[i],sT,h)
        plot_y_m192_i[i]= Erro_Relativo(T_m192_i[i],sf,h)


    plt.figure(1) #Erro da temperatura_x_tempo
    plt.plot(plot_x,plot_y_m24,'-*',label=u"MCE h = 0.167 m")
    plt.plot(plot_x,plot_y_m48,'-*',label=u"MCE h = 0.083 m")
    plt.plot(plot_x,plot_y_m96,'-*',label=u"MCE h = 0.042 m")
    plt.plot(plot_x,plot_y_m192,'-*',label=u"MCE h = 0.021 m")
    plt.xlabel(u"Tamanho do passo de tempo")
    plt.ylabel(u'Erro Relativo')
    plt.legend(loc='best')
    plt.grid()
    
    plt.figure(2) #Erro da temperatura_x_espaco
    plt.plot(plot_x,plot_y_m24_e,'-*',label=u"MCE h = 0.167 m")
    plt.plot(plot_x,plot_y_m48_e,'-*',label=u"MCE h = 0.083 m")
    plt.plot(plot_x,plot_y_m96_e,'-*',label=u"MCE h = 0.042 m")
    plt.plot(plot_x,plot_y_m192_e,'-*',label=u"MCE h = 0.021 m")
    plt.xlabel(u"Tamanho do passo de tempo")
    plt.ylabel(u'Erro Relativo')
    plt.legend(loc='best')
    plt.grid()
    
    plt.figure(3) #Erro da interface_x_tempo
    plt.plot(plot_x,plot_y_m24_i,'-*',label=u"MCE h = 0.167 m")
    plt.plot(plot_x,plot_y_m48_i,'-*',label=u"MCE h = 0.083 m")
    plt.plot(plot_x,plot_y_m96_i,'-*',label=u"MCE h = 0.042 m")
    plt.plot(plot_x,plot_y_m192_i,'-*',label=u"MCE h = 0.021 m")
    plt.xlabel(u"Tamanho do passo de tempo")
    plt.ylabel(u'Erro Relativo')
    plt.legend(loc='best')
    plt.grid()
    
    plt.show()
