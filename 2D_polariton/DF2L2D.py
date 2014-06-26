# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 09:50:39 2014

@author: memo
"""
from numpy import*
from scipy import*
from scipy.linalg import toeplitz
from numpy.linalg import inv
import sys,time

# -*- coding: utf-8 -*-
"""
Created on Tue May 27 04:28:40 2014

@author: memo
"""
#######################################
########################################
###Polariton dos ramas en 2D###
def DF2L2R2D(Psi0c,Psi0x,V0c,V0x,f,a0c,a0x,gR,g,gammaC,gammaX,X,T0,puntos=2,precision=10):
    def CS(num):
        if abs(num)<10:
            if num<0:
                txt="-0%d"%int(abs(num))
            else:
                txt="0%d"%int(num)
        else:
            txt="%d"%int(num)
        return txt

    def hms(tf):
        hor=(int(tf/3600))  
        minu=int((tf-(hor*3600))/60)  
        seg=int(tf-((hor*3600)+(minu*60)))  
        tfx=CS(hor)+"h"+CS(minu)+"m"+CS(seg)+"s|"
        return tfx    
    
    
    def Sc(Psic,Psix,F):
        return A1c*Psic+dot(alphaC,Psic)+dot(Psic,alphaC)-1.0j*gR*Psix-1.0j*F    
    def Sx(Psic,Psix):
        return A1x*Psix+dot(alphaX,Psix)+dot(Psix,alphaX)-1.0j*g*abs(Psix)**2*Psix-1.0j*gR*Psic    
    
    Dx=X[1]-X[0]
    Nx=len(X)
#    X2,Y2=meshgrid(X,X)
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
    alfjc=insert(1.0j*a0c*Vh/(2*Dx**2),0,0,0)
    alphaC=toeplitz(alfjc,alfjc)
    A1c=1.0j*C0*a0c/(Dx**2)-1.0j*V0c##    
    
    #####Parte exciton
    F0jx=Psi0x
    PsiFx=[F0jx]
    alfjx=insert(1.0j*a0x*Vh/(2*Dx**2),0,0,0)
    alphaX=toeplitz(alfjx,alfjx)
    A1x=1.0j*C0*a0x/(Dx**2)-1.0j*V0x##
    
################# Barra de estado
    nl=int((Nt+1)/10)
    print "[__h__m__s|"+(10*"#")+"]"
    sys.stdout.write("[")
    start = time.time()
#################
    for n in range(1,Nt+1):
        if n==3:
            end=time.time()
            tfi=(end-start)*(Nt+1)
            sys.stdout.write(hms(tfi))        
        F=f(n*Dt)
        #primer orden   
        K1c=Dt*Sc(F0jc,F0jx,F)
        K1x=Dt*Sx(F0jc,F0jx)
        ##segundo orden
        K2c=Dt*Sc(F0jc+K1c/2,F0jx+K1x/2,F)
        K2x=Dt*Sx(F0jc+K1c/2,F0jx+K1x/2)
        ##Tercer orden
        K3c=Dt*Sc(F0jc+K2c/2,F0jx+K2x/2,F)
        K3x=Dt*Sx(F0jc+K2c/2,F0jx+K2x/2)
        ###Funciones
        F0jc=F0jc+K1c/6.0+2.0*K2c/3.0+K3c/6.0
        F0jc=(1-gammaC*Dt)*F0jc
        F0jx=F0jx+K1x/6.0+2.0*K2x/3.0+K3x/6.0
        F0jx=(1-gammaX*Dt)*F0jx
#

        PsiFc=append(PsiFc,[F0jc],axis=0) if n%precision==0 else PsiFc      
        PsiFx=append(PsiFx,[F0jx],axis=0) if n%precision==0 else PsiFx
        
        if (n%nl==0)&(n!=0): 
            sys.stdout.write("#")
    
    sys.stdout.write("]")
    
    return PsiFc,PsiFx
    
