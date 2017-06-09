# Funcoes para calculo da Entalpia
# Gessica, 2016

import numpy as np

def Entalpia(x,T,L,Tf,Ts,Tl,cf,cs,cl,rho):
    n = len(x)
    H = np.zeros(n)
    for i in range(n):
        # fase liquida
        if (T[i] >= Tf):
            H[i] = rho*L + rho*cl*(T[i] - Tf)
        # fase solida
        else:
            H[i] = cs*(T[i] - Tf)
    return H

def EntalpiaVilaReal(x,T,L,Tf,Ts,Tl,cf,cs,cl,rho):
    n = len(x)
    H = np.zeros(n)
    for i in range (n):
        # fase solida
        if (T[i] < Ts):
            H[i] = cs*(T[i] - Tf)
        # zona mushy
        elif(Ts <= T[i] <= Tl):
            H[i] = cs*(T[i]-Tf) + (L/(Tl-Ts))*(T[i]-Ts) + cf*(T[i]-Ts)
        # fase liquida
        else:
            H[i] = cs*(Ts-Tf) + L + cf*(Tl-Ts) + cl*(T[i]-Tl)
    return H

def EntalpiaL(x,T,L,Tf,Ts,Tl,cf,cs,cl,rho):
    n = len(x)
    H = np.zeros(n)
    for i in range (n):
        # fase solida
        if (T[i] < Ts):
            H[i] = cs*T[i]
        # zona mushy
        elif(Ts <= T[i] <= Tl):
            H[i] = cs*Ts + L*((T[i]-Ts)/(Tl-Ts))
        # fase liquida
        else:
            H[i] = cs*Ts + L + cl*(T[i]-Tl)
    return H
    
def EntalpiaL2D(npts,T,L,Ts,Tl,cs,cl):
    #print("EntalpiaL2D")
    H = np.zeros((npts))
    for i in range (npts):
        # fase solida
        if (T[i] < Ts):
            #if(i==5): print("%d solida %f %f" % (i,Ts,Tl))
            H[i] = cs*(T[i])
        # zona mushy
        elif(Ts <= T[i] <= Tl):
            #if(i==5): print("%d mushy %f %f" % (i,Ts,Tl))
            H[i] = cs*(Ts) + L*((T[i]-Ts)/(Tl-Ts)) + cl*(T[i]-Ts)
        # fase liquida
        else:
            #if(i==5): print("%d liq %f %f" % (i,Ts,Tl))
            H[i] = cs*(Ts) + L + cl*(T[i]-Tl) + cl*(Tl-Ts)
    return H
    
def EntalpiaL2D_Tref(npts,T,L,Ts,Tl,cs,cl,Tref):
    #print("EntalpiaL2D")
    H = np.zeros((npts))
    for i in range (npts):
        # fase solida
        if (T[i] < Ts):
            #if(i==5): print("%d solida %f %f" % (i,Ts,Tl))
            H[i] = cs*(T[i] - Tref)
        # zona mushy
        elif(Ts <= T[i] and T[i]<= Tl):
            #if(i==5): print("%d mushy %f %f" % (i,Ts,Tl))
            H[i] = cs*(Ts - Tref) + L*((T[i]-Ts)/(Tl-Ts)) + cl*(T[i]-Ts)
        # fase liquida
        else:
            #if(i==5): print("%d liq %f %f" % (i,Ts,Tl))
            H[i] = cs*(Ts - Tref) + L + cl*(T[i]-Tl) + cl*(Tl-Ts)
    return H

def Entalpia2DIso(npts,T,L,Ts,Tl,cs,cl):#_antigo1
    H = np.zeros((npts))
    for i in range (npts):
        # fase solida
        if (T[i] < Ts):
            #if(i==5): print("%d solida %f %f" % (i,Ts,Tl))
            H[i] = cs*(T[i]-Ts)            
        # fase liquida
        else:
            #if(i==5): print("%d liq %f %f" % (i,Ts,Tl))
            H[i] = cl*(T[i]-Ts) + L
    return H
   
def CapEff(npts,T,L,Ts,Tl,cs,cl):
    DeltaT = Tl - Ts
    Cp = np.zeros((npts))
    for i in range (npts):
        if(T[i] < Ts):
            Cp[i] = cs
        elif(T[i] >= Ts and T[i] <= Tl):
            Cpeff = (L/DeltaT) + cl
            Cp[i] = Cpeff
        else:
            Cp[i] = cl
    return Cp

if __name__ == "__main__":

    from pylab import *

    L = 70.26
    Tm = -0.15
    cs,cl = 1.0, 1.0

    n = 2000
    T = np.linspace(-45, 40, n)

    DT = 2.0
    Ts = Tm - DT/2
    Tl = Tm + DT/2
    H = EntalpiaL2D(n, T, L, Ts, Tl, cs, cl)
    plot(T,H,label="DT=2.0")
    
    DT = 1.0
    Ts = Tm - DT/2
    Tl = Tm + DT/2
    H2 = EntalpiaL2D(n, T, L, Ts, Tl, cs, cl)
    plot(T,H2,label="DT=1.0")

    DT = 0.5
    Ts = Tm - DT/2
    Tl = Tm + DT/2
    H3 = EntalpiaL2D(n, T, L, Ts, Tl, cs, cl)
    plot(T,H3,'*',label="DT=0.5")

    DT = 0.25
    Ts = Tm - DT/2
    Tl = Tm + DT/2
    H4 = EntalpiaL2D(n, T, L, Ts, Tl, cs, cl)
    plot(T,H4,'--',label="DT=0.25")

    Ts = Tm 
    Tl = Tm 
    Hi = Entalpia2DIso(n, T, L, Ts, Tl, cs, cl)
    plot(T,Hi,"-^",label="isotermica")
    
    DT = 1.0
    Ts = Tm - DT/2.0
    Tl = Tm + DT/2.0
    Ce = CapEff(n, T, L, Ts, Tl, cs, cl)
    plot(T,Ce,'-',label="DT = 1.0")

    DT = 0.5
    Ts = Tm - DT/2.0
    Tl = Tm + DT/2.0
    Ce2 = CapEff(n, T, L, Ts, Tl, cs, cl)
    plot(T,Ce2,'-',label="DT = 0.5")    

    legend(loc="best")
    xlabel("Temperatura")
    ylabel("Entalpia")
    show()
