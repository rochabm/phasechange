#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Problema do Calor 2D
# MEF com retangulos lineares
# Bernardo M. Rocha
# Gessica L. Siqueira
# 2016
# -----------------------------------------------------------------------------

import os, sys, shutil
import numpy as np
from math import sqrt, sin, pi, exp, atan
from mef import *
from entalpia import *
import solucao_vila_real as sol
import matplotlib.pyplot as plt

# Some dirty globals
# Phase change

Tini = 0.0
Tm = -1.0#-0.15
Tw = -45.0
Lh = 70.26
Cps    = 1.0
Cpl    = 1.0
DeltaT = None
Ts = None     #Tm - (DeltaT/2.0)
Tl = None     #Tm + (DeltaT/2.0)
den = 1.0
ordem = 2

# *** Use or not LUMPING of the mass matrix ***
UseLump = False#True

# -----------------------------------------------------------------------------
# Matriz da Entalpia
def MassMatrixEntalpia(nen,p,t,Tn,Ts,Tl,lump=False):
    npts = np.size(p,1)
    ntri = np.size(t,1)
    M = np.zeros((npts,npts))

    Hn = EntalpiaL2D(npts,Tn,Lh,Ts,Tl,Cps,Cpl)

    if(nen==3):   shape = ShapeP1
    elif(nen==4): shape = ShapeQ1
    elif(nen==9): shape = ShapeQ2

    # get gaussian quadrature data
    w, q = GaussQuad(5)

    for k in range(ntri):
        loc2glb = t[:,k]
        x = p[0, loc2glb]
        y = p[1, loc2glb]
        Mk = np.zeros((nen,nen))

        qr, qs = q[:,0], q[:,1]

        # valores de T e H no elemento
        Te = Tn[loc2glb]
        He = Hn[loc2glb]

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
                # Morgan e Lemmon
                dHdT = sqrt((dHqx*dHqx + dHqy*dHqy)/(dTdxdTdy))
                # DelGuidice
                #dHdT = sqrt((dHqx*dTqx + dHqy*dTqy)/(dTdxdTdy))
                Cp = dHdT
            else:
                # cap effective
                if(Tq < Ts):
                    Cp = Cps
                elif(Tq >= Ts and Tq <= Tl):
                    Cpeff = (Lh/DeltaT) + Cpl
                    Cp = Cpeff
                else:
                    Cp = Cpl

            Mk = Mk + Cp * (np.outer(S,S)) * detjxw

        for il in range(nen):
            i = loc2glb[il]
            for jl in range(nen):
                j = loc2glb[jl]
                M[i,j] = M[i,j] + Mk[il,jl]

    # Mass Lumping da Matriz de Massa
    if(lump):
        for i in range(len(M)):
            aux = 0.0
            for j in range(len(M)):
                aux = aux + M[i,j]
            M[i,:] = 0.0
            M[i,i] = aux
    return M

#Matriz da Capacidade Efetiva ID
def MassMatrixCapEffID(nen,p,t,Tn,Ts,Tl,lump=False):
    npts = np.size(p,1)
    ntri = np.size(t,1)
    M = np.zeros((npts,npts))
    Cpeff =  1.0 + Lh/(DeltaT)

    if(nen==3):   shape = ShapeP1
    elif(nen==4): shape = ShapeQ1
    elif(nen==9): shape = ShapeQ2

    # get gaussian quadrature data
    w, q = GaussQuad(ordem)

    for k in range(ntri):
        loc2glb = t[:,k]
        x = p[0, loc2glb]
        y = p[1, loc2glb]
        Te = Tn[loc2glb]
        Mk = np.zeros((nen,nen))
        if ((Te[0] - Te[1]) > 1e-6):
            Mk,caso = IntegracaoGaussDescontinua(shape,nen,ordem,Cps,Cpl,Cpeff,y,Te,Tm,DeltaT,x,Tl,Ts)
        else:
            Mk = IntegracaoGauss(shape,nen,x,y,Te,ordem,Cps,Cpl,Cpeff,Ts,Tl)

        for il in range(nen):
            i = loc2glb[il]
            for jl in range(nen):
                j = loc2glb[jl]
                M[i,j] = M[i,j] + Mk[il,jl]

    # Mass Lumping
    if(lump):
        for i in range(len(M)):
            aux = 0.0
            for j in range(len(M)):
                aux = aux + M[i,j]
            M[i,:] = 0.0
            M[i,i] = aux

    return M

# -----------------------------------------------------------------------------

def MassMatrixCapEff(nen,p,t,Tn,Ts,Tl,lump=False):
    npts = np.size(p,1)
    ntri = np.size(t,1)
    M = np.zeros((npts,npts))

    if(nen==3):   shape = ShapeP1
    elif(nen==4): shape = ShapeQ1
    elif(nen==9): shape = ShapeQ2

    # get gaussian quadrature data
    w, q = GaussQuad(2)

    for k in range(ntri):
        loc2glb = t[:,k]
        x = p[0, loc2glb]
        y = p[1, loc2glb]

        Mk = np.zeros((nen,nen))
        qr, qs = q[:,0], q[:,1]
        Te = Tn[loc2glb]
        for iq in range(len(w)):
            r, s = qr[iq], qs[iq]
            S, dSdx, dSdy, detj = Isoparametric(x,y,r,s,shape)
            detjxw = detj * w[iq]
            # interpolate temp
            Tq = np.dot(S,Te)
            # cap effective
            if(Tq < Ts):
                Cp = Cps
            elif(Tq >= Ts and Tq <= Tl):
                Cpeff =  1.0 + Lh/(DeltaT)
                Cp = Cpeff
            else:
                Cp = Cpl
            Mk = Mk + Cp * (np.outer(S,S)) * detjxw
        for il in range(nen):
            i = loc2glb[il]
            for jl in range(nen):
                j = loc2glb[jl]
                M[i,j] = M[i,j] + Mk[il,jl]

    # Mass Lumping
    if(lump):
        for i in range(len(M)):
            aux = 0.0
            for j in range(len(M)):
                aux = aux + M[i,j]
            M[i,:] = 0.0
            M[i,i] = aux

    return M

# -----------------------------------------------------------------------------

def StiffnessMatrix(ordem,p,t,a):
    npts = np.size(p,1)
    ntri = np.size(t,1)
    A = np.zeros((npts,npts))
    nen = np.shape(t)[0]

    if(nen==3):   shape = ShapeP1
    elif(nen==4): shape = ShapeQ1
    elif(nen==9): shape = ShapeQ2

    w, q = GaussQuad(ordem)
    for k in range(ntri):
        loc2glb = t[:,k]
        x = p[0, loc2glb]
        y = p[1, loc2glb]
        # integracao numerica
        Ak = np.zeros((nen,nen))
        qr, qs = q[:,0], q[:,1]
        for iq in range(len(w)):
            r, s = qr[iq], qs[iq]
            S, dSdx, dSdy, detj = Isoparametric(x,y,r,s,shape)
            detjxw = detj * w[iq] * a(0,0)
            Ak = Ak + (np.outer(dSdx,dSdx) + np.outer(dSdy,dSdy)) * detjxw
        for il in range(nen):
            i = loc2glb[il]
            for jl in range(nen):
                j = loc2glb[jl]
                A[i,j] = A[i,j] + Ak[il,jl]
    return A

# -----------------------------------------------------------------------------

def FixedTemp(A, temperature, unodes, uvalues):
    n = len(b)
    for j in range(len(unodes)):
        k = unodes[j]
        # condicao Dirichlet em x=xk
        val = uvalues[j]
        A[k,k] = 1.0
        temperature[k] = val
        for i in range(n):
            if i != k:
                A[i,k] = 0
                A[k,i] = 0

# ----------------------------------------------------------------------------

def LoadMesh(basename):
    """
    Le malha no formato do Triangle - tri or quad
    nelem 3 0 -> triangulo
    nelem 4 0 -> quadrilatero (elemento bilinear 4 nos)
    #nelem 9 0 -> quadrilatero (elemento biquadratico 9 nos)
    """
    pfile = basename + ".node"
    tfile = basename + ".ele"

    pts = np.loadtxt(pfile,skiprows=1,comments='#')
    mpt = pts[:,1:3]

    # elementos
    #descobre tipo de elemento
    fele = open(tfile,"r")
    l = fele.readline()
    ls = l.split()
    tele = int(ls[1])
    fele.close()

    # le triangulos
    if(tele==3):
        ele = np.loadtxt(tfile,dtype=np.int32,skiprows=1,comments='#')
        ele[:] = ele[:] - 1
        ele = ele[:,1:4]
    elif(tele==4):
        ele = np.loadtxt(tfile,dtype=np.int32,skiprows=1,comments='#')
        ele[:] = ele[:] - 1
        ele = ele[:,1:5]
    elif(tele==9):
        ele = np.loadtxt(tfile,dtype=np.int32,skiprows=1,comments='#')
        ele[:] = ele[:] - 1
        ele = ele[:,1:10]

    mpt = mpt.transpose()
    ele = ele.transpose()

    return mpt, ele

# ----------------------------------------------------------------------------

def SaveVTK(fname, p, t, u):
    """
    Salva malha em formato VTK.
    Elementos Tri, Quad Bilinear e Quad Biquadratico
    """
    npts = np.size(p, 1)
    nele = np.size(t, 1)
    nen  = np.shape(t)[0]

    f = open(fname, 'w')
    f.write("# vtk DataFile Version 2.0\n")
    f.write('2D Unstructured Grid of Linear Triangles\n')
    f.write('ASCII\n')
    f.write('\n')

    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS %d float\n' % (npts))
    for i in range(npts):
        f.write('%f %f %f\n' % (p[0, i], p[1, i], 0.0))
    f.write('\n')

    if(nen==3): # triangulo
        cell_type = 5
        f.write('CELLS %d %d\n' % (nele, 4 * nele))
        for i in range(nele):
            f.write('3 %d %d %d\n' % (t[0, i], t[1, i], t[2, i]))
        f.write('\n')
    elif(nen==4): # quad bilinear
        cell_type = 9
        f.write('CELLS %d %d\n' % (nele, 5 * nele))
        for i in range(nele):
            f.write('4 %d %d %d %d\n' % (t[0,i], t[1,i], t[2,i], t[3,i]))
        f.write('\n')
    elif(nen==9): # quad biquadratico
        cell_type = 23
        f.write('CELLS %d %d\n' % (nele, 9 * nele))
        for i in range(nele):
            f.write('8 %d %d %d %d %d %d %d %d\n'
                    % (t[0,i], t[1,i], t[2,i], t[3,i], t[4,i], t[5,i], t[6,i], t[7,i]))  # ,t[8,i]))
        f.write('\n')
    else:
        print("erro: tipo de elemento nao reconhecido")
        sys.exit(1)

    f.write('CELL_TYPES %d\n' % nele)
    for i in range(nele):
        f.write('%d\n' % cell_type)
    f.write('\n')

    f.write('POINT_DATA %d\n' % npts)
    f.write('SCALARS scalar_field float\n')
    f.write('LOOKUP_TABLE default\n')
    for i in range(npts):
        f.write('%f\n' % u[i])
    f.write('\n')
    f.close()

# ----------------------------------------------------------------------------

def EscreveInfo(Tm,Tl,Ts,Tini,Tw,Lh,Cpl,Cps,den,
                UseLump,ft,dt,DeltaT,malha,tr,theta,matr):
    fsaida = open("saida/info.txt","w")
    if(matr == 1):
        fsaida.write('METODO ENTALPIA \n\n\n')
    elif(matr == 2):
        fsaida.write('METODO CAP. EFETIVA \n\n\n')
    else:
        fsaida.write('METODO CAP. EFETIVA ID \n\n\n')

    if (theta == 1.0):
        fsaida.write('Metodo de Euler  theta = %f\n\n' % theta)
    elif (theta == 0.5):
        fsaida.write('Metodo de Crank Nicolson  theta = %f\n\n' % theta)
    else:
        fsaida.write('Metodo de Galerkin theta = %f\n\n' % theta)
    fsaida.write('Temperatura de mudanca de fase Tm = %f\n\n' % Tm)
    fsaida.write('Temperatura fase liquida Tl = %f\n\n' % Tl)
    fsaida.write('Temperatura fase solida Ts = %f\n\n' % Ts)
    fsaida.write('Temperatura inicial em t=0 Tini = %f\n\n' % Tini)
    fsaida.write('Temperatura prescrita na face x=0 Tw = %f\n\n' % Tw)
    fsaida.write('Calor latente Lh = %f\n\n' % Lh)
    fsaida.write('Calor sensivel fase liquida cl = %f\n\n' % Cpl)
    fsaida.write('Calor sensivel fase solida cs = %f\n\n' % Cps)
    fsaida.write('Densidade pho = %f\n\n' % den)
    fsaida.write('Lumping na matriz de massa, UseLump = %s\n\n' % UseLump)
    fsaida.write('Tempo final tF = %f\n\n' % ft)
    fsaida.write('Tamanho do passo de tempo dt = %f\n\n' % dt)
    fsaida.write('Valor do delta de temperatura, DeltaT=%f\n\n' % DeltaT)
    fsaida.write('Malha utilizada = %s\n\n' % malha)
    fsaida.write("Salvando vtk a cada %f segundos =\n\n" % tr)
    fsaida.close()

# ------------------------------------------------------------------------------

if __name__ == "__main__":

    if(len(sys.argv) != 6):
        print("\n Usage: MetodoCapEfetiva2D [malha] [tempoFinal] [passoTempo] [deltaTemp] [1-ME 2-MCE 3-MCE ID]\n")
        sys.exit(1)

    # pega argumentos da linha de comando
    malha = str(sys.argv[1])
    ft = float(sys.argv[2])
    dt = float(sys.argv[3])
    DeltaT = float(sys.argv[4])
    matriz_metodo = float(sys.argv[5])

    metodo = 3

    # lambda functions - problem dependent
    a = lambda x,y: 1.08 # condutividade termica em x
    f = lambda x,y: 0.0 # condutividade termica em y

    # read mesh (nodes, elements)
    p,t  = LoadMesh(malha)
    npts = len(p[0])
    nen = np.shape(t)[0]
    nquad = np.size(t, 1)  # quantidade de elementos em x

    # discretizacao temporal e espacial
    tt = np.arange(0.0, ft + dt, dt)
    dist_x = abs(p[0,0] - p[0,1])
    dx = dist_x/nquad
    xx = np.linspace(0.0,dist_x,nquad+1)

    # temperatura solidus / liquidus
    Ts = Tm - (DeltaT / 2.0)
    Tl = Tm + (DeltaT / 2.0)

    # Cria diretorio saida
    outdir = "saida/"
    if(os.path.isdir(outdir)):
        shutil.rmtree(outdir)
    os.mkdir(outdir)

    # dados do ponto para analise
    pontoX1 = [1.0, 0.05]
    indxPontoX1 = None
    for i in range(npts):
        x,y = p[0,i], p[1,i]
        if(abs(x - pontoX1[0]) < 1e-8 and abs(y - pontoX1[1]) < 1e-8):
            indxPontoX1 = i
            break
    print("Indice do ponto para salvar dados %d" % indxPontoX1)

    # condicao inicial
    u0 = np.zeros((npts))
    for i in range(npts):
        x,y = p[0,i], p[1,i]
        u0[i] = Tini
        # regiao do contorno
        if(x < 1.0e-5):
            u0[i] = Tw

    if (matriz_metodo == 1):
        MassMatrix = MassMatrixEntalpia
        matr=1
    elif(matriz_metodo == 2):
        MassMatrix = MassMatrixCapEff
        matr=2
    else:
        MassMatrix = MassMatrixCapEffID
        matr=3


    # matrizes e vetores
    A = StiffnessMatrix(ordem,p,t,a)
    M = MassMatrix(nen,p,t,u0,Ts,Tl,UseLump)
    b = np.zeros(npts)
    R = np.zeros(npts)
    u1 = np.zeros((npts))
    ut = np.zeros((npts))
    du = np.zeros((npts))

    # matrizes iniciais
    K = (M + dt*A)
    rhs = M*u0 + dt*b

    # condicoes de contorno
    fixedIds = []
    for i in range(np.shape(p)[1]):
        x,y = p[0,i], p[1,i]
        if(x < 1e-5):
            fixedIds.append(i)
            u0[i] = Tw
    e = np.array(fixedIds)
    uval = np.ones(len(e)) * Tw
    print("Nos com condicao de Dirichlet:")
    print(e)

    #
    # TODO: conferir aqui...acho que tenho que mudar para FixTemp
    #
    BoundaryCondition(K, rhs, e, uval)

    # save initial conditions
    #SaveVTK('saida/solution_0.vtk',p,t,u0)

    fpontoX1 = open("saida/dados_pontoX1.txt","w")
    fpontoX1.write("%e %e\n" % (tt[0],u0[indxPontoX1]))

    xIn = np.zeros(len(tt))
    xIn[0] = 0.0
    Interface = open("saida/Interface.txt","w")
    Interface.write("%e %e\n" % (tt[0],xIn[0]))

    tr = 0.1
    s = tr/dt

    if(metodo == 3):
        # Esta fixo para Euler implicito
        # Theta nao muda o metodo
        # TODO: ajeitar isso
        theta = 2.0/3.0

        # loop in time
        for i in range(1,len(tt)):
            #
            # Newton method
            #
            u1 = u0.copy()

            u_ntheta = theta*u0 + (1-theta)*u0

            M = MassMatrix(nen,p,t,u0,Ts,Tl,UseLump)
            J = M + theta*dt*A
            G = M - (1.0 - theta)*dt*A
            R = np.dot(J,u1) - np.dot(G,u0)

            # modify system
            FixedTemp(J, u1, e, uval)
            R[e] = 0.0
            # solve

            #cnorm = np.linalg.norm(du)/np.linalg.norm(u0)
            #cnorm = np.linalg.norm(R)/np.linalg.norm(np.dot(K,u1))

            conv = False
            maxit = 100
            nit = 1
            ntol = 1e-3
            while(not conv):

                # solve
                du = np.linalg.solve(J, R)

                # update
                u1[:] = u1[:] - du[:]

                # continue
                u_ntheta = theta*u1 + (1-theta)*u0

                #M = MassMatrix(nen,p,t,u1,Ts,Tl,UseLump)
                M = MassMatrix(nen,p,t,u_ntheta,Ts,Tl,UseLump)

                J = M + theta*dt*A
                G = M - (1.0 - theta)*dt*A
                R = np.dot(J,u1) - np.dot(G,u0)

                # modify system
                FixedTemp(J, u1, e, uval)
                R[e] = 0.0

                # verify
                cnorm = np.linalg.norm(R) #/np.linalg.norm(np.dot(K,u1))
                if(cnorm < ntol or nit>=maxit):
                    conv = True
                    break
                    if(nit>=maxit):
                        print("maxit do newton %f" % cnorm)
                        sys.exit(0)
                if(nit>20):
                    ntol = 1.0e-2

                #du = np.linalg.solve(J, R)
                #cnorm = np.linalg.norm(du)/np.linalg.norm(u0)
                #cnorm = np.linalg.norm(R)/np.linalg.norm(np.dot(K,u1))

                nit = nit + 1
            print(nit,cnorm)
            #
            # end of Newton method
            #

            # salva dados
            if(i % s == 0): #influenciado pela escolha de dt
                print("Time %f" % (tt[i]) )
                name = "saida/solution_%d.vtk" %(i/s)
                #SaveVTK(name,p,t,u1)
                fpontoX1.write("%e %e\n" % (tt[i],u1[indxPontoX1]))

                # Salvando o vetor temperatura apenas da linha de baixo da malha
                down_u = np.zeros((nquad + 1))
                down_u[0] = u1[0]
                for j in range (4,nquad+3):
                    down_u[j-3] = u1[j]
                down_u[nquad] = u1[1]

                # criando vetor xIn - frente de solidificacao
                xIn[i] = sol.PositionInterface(xx,down_u,Tm)
                Interface.write("%e %e\n" % (tt[i],xIn[i]))

            # prepara para proximo passo
            u0 = u1.copy()

    print("fim")

    # Salva o vetor temperatura apenas da linha de baixo da
    # malha do ultimo passo de tempo pra malha de 24 elementos
    Temp_x = open("saida/dados_temp_x.txt", "w")
    Temp_x.write("%e %e\n" % (xx[0],u1[0]))
    for j in range(4, nquad+3):
        Temp_x.write("%e %e\n" % (xx[j-3],u1[j]))
    Temp_x.write("%e %e\n" % (xx[-1],u1[1]))
    Temp_x.close()

    # fecha arquivos
    fpontoX1.close()
    Interface.close()

    # escreve dados
    EscreveInfo(Tm,Tl,Ts,Tini,Tw,Lh,Cpl,Cps,den,UseLump,
                ft,dt,DeltaT,malha,tr,theta,matr)

    # plota os dados
#    ext,exTemp,exFrente = sol.solucao(4.0)
#    pdados = np.loadtxt("saida/dados_pontoX1.txt")
#    plt.plot(pdados[:,0], pdados[:,1])
#    plt.plot(ext,exTemp)
#    plt.grid()
#    plt.show()

# end of main
