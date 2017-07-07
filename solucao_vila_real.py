# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import optimize
from math import erf, exp, pi, sqrt

# Posicao da frente (numerica)
def FrontPosition(x,T,Tf):
    n = len(x)
    for i in range (n-1):
        qsi = 2.0*((Tf - T[i])/(T[i+1] - T[i])) -1.0
        if (qsi >= -1.0 and qsi <= 1.0):
            xInter = (x[i])*((1 - qsi)/2.0) + (x[i+1])*((1 + qsi)/2.0)
            return xInter

def PositionInterface(x,T,Tm):
    n = len(x)
    for i in range (1,n-1):
        if(T[i+1] != T[i]):
            qsi = 2.0*((Tm - T[i])/(T[i+1] - T[i])) -1.0
            if (qsi >= -1.0 and qsi <= 1.0):
                a = (T[i] - T[i-1])/(x[i]-x[i-1])
                b = T[i-1] - a*x[i-1]
                #y = a*x + b
                x_Inter = lambda y: (y-b)/a
                return x_Inter(Tm)
    return None

def X(lamb,Ks,dt,Nt,t):
    posicao = 2*lamb*sqrt(Ks*t*dt)
    return posicao

def temperaturaZonaSolida(Tf,Ks,lamb,Nt,dt,erf,x,t,Tw):
    p = erf(lamb)
    p1 = erf(x/(2*(Ks*t*dt)**(1./2.)))
    result = Tw + (Tf-Tw)*(p1/p)
    return result

def temperaturaZonaLiquida(Tinf,Tf,Ks,Kl,lamb,Nt,dt,x,t):
    p = 1 - erf(lamb*(Ks/Kl)**(1./2.))
    p1 = 1 - erf(x/(2*(Kl*t*dt)**(1./2.)))
    result = Tinf - (((Tinf - Tf)/p)*p1)
    return result

def X_espaco(lamb,Ks,t):
    posicao = 2*lamb*sqrt(Ks*t)
    return posicao

def temperaturaZonaSolida_espaco(Tf,Ks,lamb,x,t,Tw):
    p = erf(lamb)
    p1 = erf(x/(2*(Ks*t)**(1./2.)))
    result = Tw + (Tf-Tw)*(p1/p)
    return result

def temperaturaZonaLiquida_espaco(Tinf,Tf,Ks,Kl,lamb,x,t):
    p = 1 - erf(lamb*(Ks/Kl)**(1./2.))
    p1 = 1 - erf(x/(2*(Kl*t)**(1./2.)))
    result = Tinf - (((Tinf - Tf)/p)*p1)
    return result

def funcao(lamb,Tinf,Tf,Ks,Kl,L,cs,Tw):
    p  = erf(lamb)
    p1 = 1 - erf(lamb*(Ks/Kl)**(1./2.))
    m = (Kl/Ks)**(1./2.) * (Tinf - Tf)/(Tf-Tw) * (exp(-(Ks/Kl)*(lamb**2))/p1)
    result = exp(-lamb**2)/p - m - (lamb*L*(pi**(1./2.)))/(cs*(Tf-Tw))
    return result

def calculaLambda(funcao,lam,Tinf,Tf,Ks,Kl,L,cs,Tw):
    l = sp.optimize.newton(funcao,lam,
                           args=(Tinf,Tf,Ks,Kl,L,cs,Tw),
                           maxiter=100,
                           tol=1e-4)
    return l

def solucao(Tfim):

    parametro = open('dados/arqHowCheng.txt','r')
    lista = []
    for linha in parametro:
        valores = linha.split()
        lista.append( valores[1] )

    Ks   = float(lista[0]) # Conductibilidade na fase solida
    Kl   = float(lista[1]) # Conductibilidade na fase liquida
    Tf   = float(lista[2]) # temperatura de mudanca de fase
    Tinf = float(lista[3]) # Temperatura inicial na fase liquida T0
    Tw   = float(lista[4]) # Temperatura na fronteira no instante inicial
    L    = float(lista[5]) # calor latente
    cs   = float(lista[6]) # calor especifico na zona solida
    cl   = float(lista[7]) # calor especifico na zona liquida
    cf   = float(lista[8]) # calor especifico na zona de mudanca de fase
    p    = float(lista[9]) # pho = massa especifica

    parametro.close()
    Tf = -1.0
    Ks = 1.08
    Kl = 1.08

    #----------------------------------------------------------------------
    x    = 1.0       # x = 1 ponto onde supostamente ocorre  mudança de fase
    lam  = 0.5064    # chute inicial

    Nt = 2000
    T  = np.zeros(Nt+1)   #vetor de temperatura
    xt = np.linspace(0.0,Tfim,Nt+1)
    dt = xt[1]-xt[0]

    frenteS = np.zeros(Nt+1)

    for t in range(Nt+1):
        lamb = calculaLambda(funcao,lam,Tinf,Tf,Ks,Kl,L,cs,Tw)
        frenteS[t] = X(lamb,Ks,dt,Nt,t)
        if (t==0):
            T[t] = Tinf
        else:
            if (frenteS[t] > x):
                T[t] = temperaturaZonaSolida(Tf,Ks,lamb,Nt,dt,x,t,Tw)
            elif (frenteS[t] < x):
                T[t] = temperaturaZonaLiquida(Tinf,Tf,Ks, Kl,lamb,Nt,dt,x,t)
            else:
                T[t] = Tf
            lam = lamb

    return xt,T,frenteS

def solucaoEspaco(arq,t,n,xf):
    """
    Retorna vetores com a distribuicao de tempetura em um instante de tempo t
    """
    
    parametros = open(arq,'r')
    lista = []
    for linha in parametros:
        valores = linha.split()
        lista.append( valores[1] )
    Ks   = float(lista[0]) # Conductibilidade na fase solida
    Kl   = float(lista[1]) # Conductibilidade na fase liquida
    Tf   = float(lista[2]) # temperatura de mudanca de fase
    Tinf = float(lista[3]) # Temperatura inicial na fase liquida T0
    Tw   = float(lista[4]) # Temperatura na fronteira no instante inicial
    L    = float(lista[5]) # calor latente
    cs   = float(lista[6]) # calor especifico na zona solida
    cl   = float(lista[7]) # calor especifico na zona liquida
    cf   = float(lista[8]) # calor especifico na zona de mudanca de fase
    p    = float(lista[9]) # pho = massa especifica
    parametros.close()    

    lam = 0.5064                 # chute inicial    
    T  = np.zeros(n+1)           # temperaturas
    x  = np.linspace(0.0,xf,n+1) # espaco

    T[0] = Tw
    for i in range(1,n+1):
        xx = x[i]
        lamb = calculaLambda(funcao,lam,Tinf,Tf,Ks,Kl,L,cs,Tw)
        frente = X_espaco(lamb,Ks,t)
        if (frente > xx):
            T[i] = temperaturaZonaSolida_espaco(Tf,Ks,lamb,xx,t,Tw)
        elif (frente < xx):
            T[i] = temperaturaZonaLiquida_espaco(Tinf,Tf,Ks,Kl,lamb,xx,t)
        else:
            T[i] = Tf
        lam = lamb
    return x,T

def solucaoEspacoX(arq,t,x):
    """
    Retorna vetores com a distribuicao de tempetura em um instante de tempo t
    """
    
    parametros = open(arq,'r')
    lista = []
    for linha in parametros:
        valores = linha.split()
        lista.append( valores[1] )
    Ks   = float(lista[0]) # Conductibilidade na fase solida
    Kl   = float(lista[1]) # Conductibilidade na fase liquida
    Tf   = float(lista[2]) # temperatura de mudanca de fase
    Tinf = float(lista[3]) # Temperatura inicial na fase liquida T0
    Tw   = float(lista[4]) # Temperatura na fronteira no instante inicial
    L    = float(lista[5]) # calor latente
    cs   = float(lista[6]) # calor especifico na zona solida
    cl   = float(lista[7]) # calor especifico na zona liquida
    cf   = float(lista[8]) # calor especifico na zona de mudanca de fase
    p    = float(lista[9]) # pho = massa especifica
    parametros.close()
    
    lam = 0.5064      # chute inicial    
    n = len(x)        # tamanho da malha
    T  = np.zeros(n)  # temperaturas

    T[0] = Tw
    for i in range(1,n):
        xx = x[i]
        lamb = calculaLambda(funcao,lam,Tinf,Tf,Ks,Kl,L,cs,Tw)
        frente = X_espaco(lamb,Ks,t)
        if (frente > xx):
            T[i] = temperaturaZonaSolida_espaco(Tf,Ks,lamb,xx,t,Tw)
        elif (frente < xx):
            T[i] = temperaturaZonaLiquida_espaco(Tinf,Tf,Ks,Kl,lamb,xx,t)
        else:
            T[i] = Tf
        lam = lamb
    return T
    
if __name__ == "__main__":

    arq = "dados/arqVilaReal.txt"
    
    #
    # testa a solucao analitica para a distribuicao espacial da temperatura
    #
    x,u1 = solucaoEspaco(arq,2.0,40,4.0)
    plt.plot(x,u1)
    plt.show()

    x = np.linspace(0,4.0,41)
    u2 = solucaoEspacoX(arq,2.0,x)
    plt.plot(x,u2)
    plt.show()
    