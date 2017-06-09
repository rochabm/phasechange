# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Solucao Analitica problema de Mudanca de Fase
# Bernardo M. Rocha
# Gessica L. Siqueira
# 2016
# -----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import optimize
from math import erf, exp, pi, sqrt

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

# Solucao analitica para solidificacao
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

def temperaturaZonaSolida_espaco(Tf,Ks,lamb,erf,x,t,Tw):
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
                           args=(Tinf,Tf,Ks,Kl,L,cs,Tw),maxiter=100,tol=1e-4)
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
            T[t] = Tinf #condição inicial
            #print'frenteSolidificacao',frenteS
        else:
            #print'frenteSolidificacao',frenteS
            if (frenteS[t] > x): #or T[t] < Tf):
                T[t] = temperaturaZonaSolida(Tf, Ks, lamb, Nt, dt,erf,x,t,Tw)
                #print 'T zona solida',T
            elif (frenteS[t] < x): #or T[t] > Tf):
                T[t] = temperaturaZonaLiquida(Tinf,Tf,Ks, Kl,lamb,Nt,dt,x,t)
                #print 'T zona liquida',T
            else:
                T[t] = Tf
            lam = lamb
    return xt,T,frenteS

def solucaoEspaco(t,nquad,xf):

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

    #----------------------------------------------------------------------

    lam  = 0.5064    # chute inicial
    T  = np.zeros(nquad+1)   #vetor de temperatura
    x = np.linspace(0.0,xf,nquad+1)
    dx = x[1]-x[0]

    for i in range(nquad+1):
        xx = x[i]
        lamb = calculaLambda(funcao,lam,Tinf,Tf,Ks,Kl,L,cs,Tw)
        frente = X_espaco(lamb,Ks,t)
        if (i==0):
            T[i] = Tw
            #condição inicial
            #print'frenteSolidificacao',frenteS
        else:
            #print'frenteSolidificacao',frenteS
            if (frente > xx): #or T[t] < Tf):
                T[i] = temperaturaZonaSolida_espaco(Tf, Ks, lamb, erf, xx, t, Tw)
                #print 'T zona solida',T
            elif (frente < xx): #or T[t] > Tf):
                T[i] = temperaturaZonaLiquida_espaco(Tinf,Tf,Ks, Kl,lamb, xx, t)
                #print 'T zona liquida',T
            else:
                T[i] = Tf
            lam = lamb

    return x,T

if __name__ == "__main__":
    pass
