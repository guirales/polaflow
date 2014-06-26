# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 09:50:39 2014

@author: memo
"""
from numpy import*
from scipy import*
from scipy.linalg import toeplitz
from numpy.linalg import inv

###Polariton una rama##    
def DF2L(Psi0,V0,g,a0,X,T,puntos=1):
    def S(Psi):
        return dot(A1-1.0j*g*diag(abs(Psi)**2),Psi)    
    
    Dx=X[1]-X[0]
    Nx=len(X)
    Dt=T[1]-T[0]
    Nt=len(T)
    
    F0j=Psi0
    PsiF=[F0j]
    d=Nx if puntos<2 else int(puntos)   
    B=array([[k**(2*r) for k in range(1,d)] for r in range(1,d)])
    BI=inv(B)
    C0=-2*sum(BI[::,0])
#    print C0," ",BI[::,0]
    Vj=-1.0j*diag(V0)
    gm0=1.0j*C0*a0*eye(Nx)/(2*Dx**2)
    Vh=append(BI[::,0],zeros(Nx-d))
    alfj=insert(1.0j*a0*Vh/(2*Dx**2),0,0,0)
#    print alfj
    A1=toeplitz(alfj,alfj)+gm0+Vj##


    for n in range(Nt):
        K1=Dt*S(F0j)
        K2=Dt*S(F0j+K1/2)
        K3=Dt*S(F0j+K2/2)
        F0j=F0j+K1/6.0+2.0*K2/3.0+K3/6.0
#        F0j=F0j+K1
        PsiF=append(PsiF,[F0j],axis=0)
    return PsiF

###Polariton dos ramas###
def DF2L2R(Psi0c,Psi0x,V0c,V0x,a0c,a0x,gR,g,X,T0,puntos=2,precision=10):
    def Sc(Psic,Psix):
        return dot(A1c,Psic)-1.0j*gR*Psix    
    def Sx(Psic,Psix):
        return dot(A1x,Psix)-1.0j*g*abs(Psix)**2*Psix-1.0j*gR*Psic
#        return dot(A1x-1.0j*g*diag(abs(Psix)**2),Psix)-1.0j*gR*Psic    
    
    Dx=X[1]-X[0]
    Nx=len(X)
    Nt=(len(T0)-1)*precision
    T=linspace(T0[0],T0[-1],Nt+1)
    Dt=T[1]-T[0]
    
#    d=Nx if puntos<2 else int(puntos)   
    d=puntos
    B=array([[k**(2*r) for k in range(1,d)] for r in range(1,d)])
    BI=inv(B)
    C0=-2*sum(BI[::,0]) 
    Vh=append(BI[::,0],zeros(Nx-d))

    
    #### Parte cavidad
    F0jc=Psi0c
    PsiFc=[F0jc]
    Vjc=-1.0j*diag(V0c)
    gm0c=1.0j*C0*a0c*eye(Nx)/(2*Dx**2)
    alfjc=insert(1.0j*a0c*Vh/(2*Dx**2),0,0,0)
    A1c=toeplitz(alfjc,alfjc)+gm0c+Vjc##    
    
    #####Parte exciton
    F0jx=Psi0x
    PsiFx=[F0jx]
    Vjx=-1.0j*diag(V0x)
    gm0x=1.0j*C0*a0x*eye(Nx)/(2*Dx**2)
    alfjx=insert(1.0j*a0x*Vh/(2*Dx**2),0,0,0)
    A1x=toeplitz(alfjx,alfjx)+gm0x+Vjx##


    for n in range(1,Nt+1):
        #primer orden   
        K1c=Dt*Sc(F0jc,F0jx)
        K1x=Dt*Sx(F0jc,F0jx)
        ##segundo orden
        K2c=Dt*Sc(F0jc+K1c/2,F0jx+K1x/2)
        K2x=Dt*Sx(F0jc+K1c/2,F0jx+K1x/2)
        ##Tercer orden
        K3c=Dt*Sc(F0jc+K2c/2,F0jx+K2x/2)
        K3x=Dt*Sx(F0jc+K2c/2,F0jx+K2x/2)
        ###Funciones
        F0jc=F0jc+K1c/6.0+2.0*K2c/3.0+K3c/6.0
        F0jx=F0jx+K1x/6.0+2.0*K2x/3.0+K3x/6.0
#

        PsiFc=append(PsiFc,[F0jc],axis=0) if n%precision==0 else PsiFc      
        PsiFx=append(PsiFx,[F0jx],axis=0) if n%precision==0 else PsiFx
    return PsiFc,PsiFx
    
########################################
########################################
###Polariton dos ramas en 2D###
def DF2L2R2D(Psi0c,Psi0x,V0c,V0x,a0c,a0x,gR,g,X,T0,puntos=2,precision=10):
    def Sc(Psic,Psix):
        return dot(A1c,Psic)+dot(Psic,alphaC)-1.0j*gR*Psix    
    def Sx(Psic,Psix):
        return dot(A1x,Psix)+dot(Psix,alphaX)-1.0j*g*abs(Psix)**2*Psix-1.0j*gR*Psic    
    
    Dx=X[1]-X[0]
    Nx=len(X)
#    Dy=Y[1]-Y[0]
#    Ny=len(Y)
    Nt=(len(T0)-1)*precision
    T=linspace(T0[0],T0[-1],Nt+1)
    Dt=T[1]-T[0]
    
#    d=Nx if puntos<2 else int(puntos)   
    d=puntos
    B=array([[k**(2*r) for k in range(1,d)] for r in range(1,d)])
    BI=inv(B)
    C0=-2*sum(BI[::,0]) 
    Vh=append(BI[::,0],zeros(Nx-d))

    
    #### Parte cavidad
    F0jc=Psi0c
    PsiFc=[F0jc]
    Vjc=-1.0j*diag(V0c)
    gm0c=1.0j*C0*a0c*eye(Nx)/(Dx**2)
    alfjc=insert(1.0j*a0c*Vh/(2*Dx**2),0,0,0)
    alphaC=toeplitz(alfjc,alfjc)
    A1c=alpC+gm0c+Vjc##    
    
    #####Parte exciton
    F0jx=Psi0x
    PsiFx=[F0jx]
    Vjx=-1.0j*diag(V0x)
    gm0x=1.0j*C0*a0x*eye(Nx)/(Dx**2)
    alfjx=insert(1.0j*a0x*Vh/(2*Dx**2),0,0,0)
    alphaX=toeplitz(alfjx,alfjx)
    A1x=alpX+gm0x+Vjx##


    for n in range(1,Nt+1):
        #primer orden   
        K1c=Dt*Sc(F0jc,F0jx)
        K1x=Dt*Sx(F0jc,F0jx)
        ##segundo orden
        K2c=Dt*Sc(F0jc+K1c/2,F0jx+K1x/2)
        K2x=Dt*Sx(F0jc+K1c/2,F0jx+K1x/2)
        ##Tercer orden
        K3c=Dt*Sc(F0jc+K2c/2,F0jx+K2x/2)
        K3x=Dt*Sx(F0jc+K2c/2,F0jx+K2x/2)
        ###Funciones
        F0jc=F0jc+K1c/6.0+2.0*K2c/3.0+K3c/6.0
        F0jx=F0jx+K1x/6.0+2.0*K2x/3.0+K3x/6.0
#

        PsiFc=append(PsiFc,[F0jc],axis=0) if n%precision==0 else PsiFc      
        PsiFx=append(PsiFx,[F0jx],axis=0) if n%precision==0 else PsiFx
    return PsiFc,PsiFx
    
