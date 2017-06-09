# Modulos com as funcoes do Metodo dos Elementos Finitos (MEF)
# Bernardo, Gessica
# 2016

import sys
import numpy as np
from math import sqrt
from numpy.polynomial.legendre import leggauss

# -----------------------------------------------------------------------------

def ShapeP1(r,s):
    S = np.array([1.0-r-s, r, s])
    dSdr = np.array([-1.0, 1.0, 0.0])
    dSds = np.array([-1.0, 0.0, 1.0])
    return S, dSdr, dSds

# -----------------------------------------------------------------------------

def ShapeQ1(r,s):
    S = np.array([0.25*(1-r)*(1-s),
                  0.25*(1+r)*(1-s),
                  0.25*(1+r)*(1+s),
                  0.25*(1-r)*(1+s)])
    dSdr = np.array([-0.25*(1-s), 0.25*(1-s), 0.25*(1+s), -0.25*(1+s)])
    dSds = np.array([-0.25*(1-r),-0.25*(1+r), 0.25*(1+r),  0.25*(1-r)])
    return S, dSdr, dSds

# -----------------------------------------------------------------------------

def ShapeQ2(r,s):
    S = np.array([0.25*(r*r-r)*(s*s-s),
                  0.25*(r*r+r)*(s*s-s),
                  0.25*(r*r+r)*(s*s+s),
                  0.25*(r*r-r)*(s*s+s),
                  0.5*(1- r*r)*(s*s-s),
                  0.5*(r*r+r)*(1- s*s),
                  0.5*(1- r*r)*(s*s+s),
                  0.5*(r*r-r)*(1- s*s),
		      (1 - r*r)*(1 - s*s)])
    dSdr = np.array([0.25*(s*s-s)*(2*r -1),
                     0.25*(s*s-s)*(2*r +1),
                     0.25*(s*s+s)*(2*r +1),
                     0.25*(s*s+s)*(2*r -1),
                     0.5*(s*s-s)*(-2*r),
                     0.5*(1 -s*s)*(2*r+1),
                     0.5*(s*s+s)*(-2*r),
                     0.5*(1-s*s)*(2*r-1),
                     (1-s*s)*(-2*r)])
    dSds = np.array([0.25*(r*r-r)*(2*s -1),
                     0.25*(r*r+r)*(2*s -1),
                     0.25*(r*r+r)*(2*s +1),
                     0.25*(r*r-r)*(2*s+1),
                     0.5*(1-r*r)*(2*s-1),
                     0.5*(r*r+r)*(-2*s),
                     0.5*(1-r*r)*(2*s+1),
                     0.5*(r*r-r)*(-2*s),
                     (1-r*r)*(-2*s)])
    return S, dSdr, dSds
    
# -----------------------------------------------------------------------------

def Isoparametric(x,y,r,s,shape):
    S, dSdr, dSds = shape(r,s)
    j11 = np.dot(dSdr,x)
    j12 = np.dot(dSdr,y)
    j21 = np.dot(dSds,x)
    j22 = np.dot(dSds,y)
    det = j11*j22 - j12*j21

    # ordering fix
    if(det<0): det*=-1

    dSdx = ( j22*dSdr - j12*dSds)/det
    dSdy = (-j21*dSdr + j11*dSds)/det
    return S, dSdx, dSdy, det

# -----------------------------------------------------------------------------

def GaussTri(precision):
    if(precision==1):
        w = np.array([1])
        q = np.array([[1./3., 1./3.]])
    elif(precision==2):
        w = np.array([1./3, 1./3., 1./3.])
        q = np.array([[1./6., 1./6.],
                      [2./3., 1./6.],
                      [1./6., 2./3.]])
    elif(precision==3):
        w = np.array([-27./28., 25./48., 25./48., 25./48.])
        q = np.array([[1./3., 1./3.],
                      [1./5., 1./5.],
                      [1./5., 3./5.],
                      [3./5., 1./5.]])
    else:
        print("erro GaussTri")
    return w, q
    
# -----------------------------------------------------------------------------

def GaussQuad(n):
    """
    n: numero de pontos e pesos
    precisao da regra em 1d: 2n-1 ou seja, integra polinomio de 
    grau ate 2n-1 exatamente     
    """
    q, w = leggauss(n)
    q2d = np.zeros((n*n, 2))
    w2d = np.zeros((n*n, 1))
    k = 0
    for i in range(n):
        for j in range(n):
            w2d[k] = w[i]*w[j]
            q2d[k,0] = q[i]
            q2d[k,1] = q[j]
            k = k + 1
    return w2d, q2d

# -----------------------------------------------------------------------------

def IntegracaoGaussDescontinua(shape,nen,ordem,cs,cl,ceff,yy,vT,Tm,DT,x,Tl,Ts):
    qsi   = 2.0*((Tm - vT[1])/(vT[0] - vT[1])) - 1.0
    qsi_l = 2.0*((Tl - vT[1])/(vT[0] - vT[1])) - 1.0
    qsi_s = 2.0*((Ts - vT[1])/(vT[0] - vT[1])) - 1.0
    pInterface = (qsi + 1.0)*0.5*(x[0] - x[1]) + x[1]
    pInterface_l = (qsi_l + 1.0)*0.5*(x[0] - x[1]) + x[1]
    pInterface_s = (qsi_s + 1.0)*0.5*(x[0] - x[1]) + x[1]
    dxi = (DT/2.0)*((x[0] - x[1])/(vT[0] - vT[1]))
    vm = np.zeros(nen)
    vs = np.zeros(nen)
    vl = np.zeros(nen)
    vTm = np.zeros(nen)
    vTs = np.zeros(nen)
    vTl = np.zeros(nen)

    if (nen == 4):
        # CASO 1---------------------------------------------------------------
        if ((pInterface_s > x[1]) and (pInterface_l < x[0]) ):
            caso = 1
            # Coordenada em x do Retangulo da parte mushy
            vm[0] = pInterface + dxi
            vm[1] = pInterface - dxi
            vm[2] = pInterface - dxi
            vm[3] = pInterface + dxi
            # Vetor de temperatura do Retangulo da parte mushy
            vTm[0] = Tm + DT/2.0
            vTm[1] = Tm - DT/2.0
            vTm[2] = Tm - DT/2.0
            vTm[3] = Tm + DT/2.0
            Im = IntegracaoParteMushy(shape,nen,ordem,cs,cl,ceff,yy,vm)
            
            # Coordenada em x do Retangulo da parte solida 
            vs[0] = vm[1] 
            vs[1] = x[1]
            vs[2] = x[2]
            vs[3] = vm[2]
            # Vetor de temperatura do Retangulo da parte solida
            vTs[0] = vTm[1] 
            vTs[1] = vT[1]
            vTs[2] = vT[2]
            vTs[3] = vTm[2]
            Is = IntegracaoParteSolida(shape,nen,ordem,cs,cl,ceff,yy,vs)
            
            # Coordenada em x do Retangulo da parte liquida
            vl[0] = x[0]
            vl[1] = vm[0]
            vl[2] = vm[3] 
            vl[3] = x[3]
            # Vetor de temperatura do Retangulo da parte liquida
            vTl[0] = vT[0]
            vTl[1] = vTm[0]
            vTl[2] = vTm[3]
            vTl[3] = vT[3]
            Il = IntegracaoParteLiquida(shape,nen,ordem,cs,cl,ceff,yy,vl)

            return Im + Is + Il,caso#,caso,(Im + Is + Il)

        # CASO 2 --------------------------------------------------------------
        elif ((pInterface_l < x[0] and pInterface_l > x[1]) and (pInterface_s < x[1]) ):
            caso = 2
            # Coordenada em x do Retangulo da parte mushy
            vm[0] = pInterface + dxi
            vm[1] = x[1]
            vm[2] = x[2]
            vm[3] = pInterface + dxi
            # Vetor de temperatura do Retangulo da parte mushy
            vTm[0] = Tm + DT/2.0
            vTm[1] = vT[1]
            vTm[2] = vT[2]
            vTm[3] = Tm + DT/2.0
            Im = IntegracaoParteMushy(shape,nen,ordem,cs,cl,ceff,yy,vm)
            
            # Coordenada em x do Retangulo da parte liquida
            vl[0] = x[0]
            vl[1] = vm[0]
            vl[2] = vm[3] 
            vl[3] = x[3]
            # Vetor de temperatura do Retangulo da parte liquida
            vTl[0] = vT[0]
            vTl[1] = vTm[0]
            vTl[2] = vTm[3]
            vTl[3] = vT[3]
            Il = IntegracaoParteLiquida(shape,nen,ordem,cs,cl,ceff,yy,vl)

            return Im + Il,caso#,caso,(Im + Il)

        # CASO 3 --------------------------------------------------------------
        elif ((pInterface_s > x[1] and pInterface_s < x[0]) and (pInterface_l > x[0]) ):
            caso = 3
            # Coordenada em x do Retangulo da parte mushy
            vm[0] = x[0]
            vm[1] = pInterface - dxi
            vm[2] = pInterface - dxi
            vm[3] = x[0]
            # Vetor de temperatura do Retangulo da parte mushy
            vTm[0] = vT[0]
            vTm[1] = Tm - DT/2.0
            vTm[2] = Tm - DT/2.0
            vTm[3] = vT[3]
            Im = IntegracaoParteMushy(shape,nen,ordem,cs,cl,ceff,yy,vm)
            
            # Coordenada em x do Retangulo da parte solida 
            vs[0] = vm[1] 
            vs[1] = x[1]
            vs[2] = x[2]
            vs[3] = vm[2]
            # Vetor de temperatura do Retangulo da parte solida
            vTs[0] = vTm[1] 
            vTs[1] = vT[1]
            vTs[2] = vT[2]
            vTs[3] = vTm[2]
            Is = IntegracaoParteSolida(shape,nen,ordem,cs,cl,ceff,yy,vs)
            return Im + Is, caso
        
        elif(vT[0] <= Ts and vT[1] <= Ts and vT[2] <= Ts and vT[3] <= Ts):
            # solido
            caso = 0.1
            vs[:] = x[:]
            I_s = IntegracaoParteSolida(shape,nen,ordem,cs,cl,ceff,yy,vs)
            return I_s, caso

        elif(vT[0] >= Tl and vT[1] >= Tl and vT[2] >= Tl and vT[3] >= Tl):
            # liquido
            caso = 0.2
            vl[:] = x[:]
            I_l = IntegracaoParteLiquida(shape,nen,ordem,cs,cl,ceff,yy,vl)
            return I_l,caso

        elif((vT[0] <= Tl ) and (vT[1] >= Ts ) and (vT[2] >= Ts) and (vT[3] <= Tl)):
            # mushy
            caso = 0.3
            vm[:] = x[:]
            I_m = IntegracaoParteMushy(shape,nen,ordem,cs,cl,ceff,yy,vm)

            return I_m,caso
        else:
            print("erro: tratamento da IntegracaoGaussDescontinua")
            sys.exit(1)
    else:
        print("erro: IntegracaoGaussDescontinua nen=9 nao implementado")

# fim IntegracaoGaussDescontinua

# -----------------------------------------------------------------------------
    
def IntegracaoParteSolida(shape,nen,ordem,cs,cl,ceff,yy,vs):
    w, q = GaussQuad(ordem)
    qr, qs = q[:,0], q[:,1]
    # Integracao na parte solida 
    Igs = np.zeros((nen,nen))
    for iq in range (len(w)):
        r, s = qr[iq], qs[iq]
        if (nen == 4):
            xs = [vs[0],vs[1],vs[2],vs[3]]
        elif(nen == 9):
            xs = [vs[0],vs[1],vs[2],vs[3],vs[4],vs[5],vs[6],vs[7],vs[8]]
        Ss, dSdx_s, dSdy_s, detj_s = Isoparametric(xs,yy,r,s,shape)    
        detjxw_s = detj_s * w[iq]
        Igs = Igs + cs*(np.outer(Ss,Ss))*detjxw_s
    return Igs

# -----------------------------------------------------------------------------

def IntegracaoParteMushy(shape,nen,ordem,cs,cl,ceff,yy,vm):
    w, q = GaussQuad(ordem)
    qr, qs = q[:,0], q[:,1]
    # Integracao na parte mushy
    Igm = np.zeros((nen,nen))
    for iq in range (len(w)):
        r, s = qr[iq], qs[iq]
        if (nen == 4):
            xm = [vm[0],vm[1],vm[2],vm[3]]
        elif(nen == 9):
            xm = [vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8]]            
        Sm, dSdx, dSdy, detj_m = Isoparametric(xm,yy,r,s,shape)
        detjxw_m = detj_m * w[iq]
        Igm = Igm + ceff*(np.outer(Sm,Sm))*detjxw_m
    return Igm

# -----------------------------------------------------------------------------

def IntegracaoParteLiquida(shape,nen,ordem,cs,cl,ceff,yy,vl):
    w, q = GaussQuad(ordem)
    qr, qs = q[:,0], q[:,1]
    # Integracao na parte liquida 
    Igl = np.zeros((nen,nen))
    for iq in range (len(w)):
        r, s = qr[iq], qs[iq]
        if (nen == 4):
            xl = [vl[0],vl[2],vl[3],vl[1]]
        elif(nen == 9):
            xl = [vl[0],vl[1],vl[2],vl[3],vl[4],vl[5],vl[6],vl[7],vl[8]]
        Sl, dSdx, dSdy, detj_l = Isoparametric(xl,yy,r,s,shape)
        detjxw_l = detj_l * w[iq]
        Igl = Igl + cl*(np.outer(Sl,Sl))*detjxw_l
    return Igl

# -----------------------------------------------------------------------------

def IntegracaoEntalpiaGauss(shape,nen,x,y,ordem,Te,He,cs,cl,Ts,Tl,ceff,DT,Tm):
    w, q = GaussQuad(ordem)
    qr, qs = q[:,0], q[:,1]
    Me = np.zeros((nen,nen))
    for iq in range(len(w)):
        r, s = qr[iq], qs[iq]
        S, dSdx, dSdy, detj = Isoparametric(x,y,r,s,shape)
        detjxw = detj * w[iq] 
        
	  # interpolate temp and enthalpy
        Tq = np.dot(S,Te)
        Hq = np.dot(S,He)
        
        # interpolar deriv temp x e y  - dSdx e dSdy
        dTqx = np.dot(dSdx,Te)    
        dTqy = np.dot(dSdy,Te)
    
        # interpolar deriv entalpia x - dSdx
        dHqx = np.dot(dSdx,He)
        dHqy = np.dot(dSdy,He)
        
        dTdxdTdy = (dTqx*dTqx + dTqy*dTqy)
        
        if(abs(dTdxdTdy) > 1.0e-10):
            # Morgan
            dHdT = sqrt((dHqx*dHqx + dHqy*dHqy)/(dTdxdTdy))
            # DelGuidice
            #dHdT = sqrt((dHqx*dTqx + dHqy*dTqy)/(dTdxdTdy))
            Cp = dHdT
            Me = IntegracaoGaussDescontinua(shape,nen,ordem,cs,cl,Cp,y,Te,Tm,DT,x,Tl,Ts)
                
        else:
            Me = IntegracaoGaussEntalpia(nen,cs,cl,ceff,Ts,Tl,Tq,detjxw,S)
    return Me

# -----------------------------------------------------------------------------

def IntegracaoGaussEntalpia(nen,cs,cl,ceff,Ts,Tl,Tq,detjxw,S):
    Ig = np.zeros((nen,nen))
    if(Tq < Ts):
        # elemento totalmente solido
        cq = cs
    elif(Ts <= Tq <= Tl):
        # elemento da zona mushy
        cq = ceff
    else:
        # elemento totalmente liquido
        cq = cl      
    Ig = Ig + cq*(np.outer(S,S))*detjxw
    return Ig                         

# -----------------------------------------------------------------------------

def IntegracaoGauss(shape,nen,xn,yn,Te,ordem,cs,cl,ceff,Ts,Tl):
    w, q = GaussQuad(ordem)
    qr, qs = q[:,0], q[:,1]
    Ig = np.zeros((nen,nen))
    for iq in range (len(w)):
        r, s = qr[iq], qs[iq]        
        S, dSdx, dSdy, detj = Isoparametric(xn,yn,r,s,shape)    
        detjxw = detj * w[iq]
        # interpola temperatura para o ponto de Gauss    
        Tq = np.dot(S,Te)
        # avalia Ceff
        if(Tq < Ts):
            # Elemento totalmente solido
            cq = cs
        elif(Ts <= Tq <= Tl):
            # Elemento da zona mushy
            cq = ceff
        else:
            # Elemento totalmente liquido
            cq = cl      
        Ig = Ig + cq*(np.outer(S,S))*detjxw
    return Ig
    
# -----------------------------------------------------------------------------

def BoundaryCondition(A, b, ua):
    n = len(b)
    # condicao Dirichlet em x=a
    A[0,0] = 1.0
    b[0] = ua
    for i in range(1,n):
        b[i] = b[i] - A[i,0]*ua
        A[i,0] = 0
        A[0,i] = 0

# -----------------------------------------------------------------------------

def BoundaryCondition(A, b, unodes, uvalues):
    n = len(b)
    for j in range(len(unodes)):
        k = unodes[j]
        # condicao Dirichlet em x=xk
        val = uvalues[j]
        A[k, k] = 1.0
        b[k] = val
        for i in range(n):
            if i != k:
                b[i] = b[i] - A[i, k] * val
                A[i, k] = 0
                A[k, i] = 0

# -----------------------------------------------------------------------------
        
def CondicaoInicial(T,Tinf):
    T[:] = Tinf

# -----------------------------------------------------------------------------    

def MatrizDeRigidez(x,k):
    n = len(x)
    A = np.zeros((n,n))
    h = x[2]-x[1]
    wg = [1.0, 1.0]
    xg = [-sqrt(3.0)/3.0, sqrt(3.0)/3.0]
    dphi1 = -0.5
    dphi2 = 0.5
    for i in range(n-1):
        ke = np.zeros((2,2))
        for ig in range(len(xg)):
            ke[0,0] += k * dphi1 * dphi1 * wg[ig] * (2.0/h)
            ke[0,1] += k * dphi1 * dphi2 * wg[ig] * (2.0/h)
            ke[1,0] += k * dphi2 * dphi1 * wg[ig] * (2.0/h)
            ke[1,1] += k * dphi2 * dphi2 * wg[ig] * (2.0/h)
        # montagem
        A[i,i]     += ke[0,0]
        A[i,i+1]   += ke[0,1]
        A[i+1,i]   += ke[1,0]
        A[i+1,i+1] += ke[1,1]
    return A

# -----------------------------------------------------------------------------

def MatrizDeRigidez2(x):
    n = len(x)
    A = np.zeros((n,n))
    h = x[2]-x[0]
    xg = np.array([-sqrt(3.0/5.0), 0.0,  sqrt(3.0/5.0)])
    wg = np.array([5.0/9.0, 8.0/9.0, 5.0/9.0])

    dphi1 = lambda x: x - 0.5
    dphi2 = lambda x: -2.0*x
    dphi3 = lambda x: x + 0.5

    for i in range(0,n-1,2):
        ke = np.zeros((3,3))
        for ig in range(len(xg)):
            xi = xg[ig]
            ke[0,0] += dphi1(xi) * dphi1(xi) * wg[ig] * (2.0/h)
            ke[0,1] += dphi1(xi) * dphi2(xi) * wg[ig] * (2.0/h)
            ke[0,2] += dphi1(xi) * dphi3(xi) * wg[ig] * (2.0/h)
            ke[1,0] += dphi2(xi) * dphi1(xi) * wg[ig] * (2.0/h)
            ke[1,1] += dphi2(xi) * dphi2(xi) * wg[ig] * (2.0/h)
            ke[1,2] += dphi2(xi) * dphi3(xi) * wg[ig] * (2.0/h)
            ke[2,0] += dphi3(xi) * dphi1(xi) * wg[ig] * (2.0/h)
            ke[2,1] += dphi3(xi) * dphi2(xi) * wg[ig] * (2.0/h)
            ke[2,2] += dphi3(xi) * dphi3(xi) * wg[ig] * (2.0/h)
        # montagem da matriz de rigidez global
        A[i,i]     += ke[0,0]
        A[i,i+1]   += ke[0,1]
        A[i,i+2]   += ke[0,2]
        A[i+1,i]   += ke[1,0]
        A[i+1,i+1] += ke[1,1]
        A[i+1,i+2] += ke[1,2]
        A[i+2,i]   += ke[2,0]
        A[i+2,i+1] += ke[2,1]
        A[i+2,i+2] += ke[2,2]
    return A

# -----------------------------------------------------------------------------
# End of MEF module
# -----------------------------------------------------------------------------