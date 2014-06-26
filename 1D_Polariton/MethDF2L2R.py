# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 11:10:28 2014

@author: memo
"""
from numpy import*
from scipy import*
from DF2L import DF2L2R
import sys,os


print "inicia calculo"

def soliton(y,t):
    x=0.3*y
    asol=0.4
    vsol=0.0
    x0sol=0.0
    return asol*exp(1.0j*(vsol*x+(asol**2-vsol**2)*t/2))/cosh(asol*(x-vsol*t-x0sol))


def gauss(x,t,sigma=15.0):
    return exp(-0.5*x**2/sigma**2)/sqrt(sqrt(pi*sigma**2))

###****potencial
def potencial(x):
    omega=0.0
    return omega**2*x**2/2.0

####bombeo delta
def f2(x,t):
    return 5.0*exp(1j*0.5*t)*Ddelta(x,0.0,3*Dx)
#f2=vectorize(f2)

####Funcion cero
def fz(x,t):
    return 0.0

###########################

#x=arange(-7,7,0.01)
x=linspace(-50,50,2**7)
Dx=x[1]-x[0]
Nx=len(x)

#t=arange(0,7,0.05)
t=linspace(0,15,2**7)
Dt=t[1]-t[0]
Nt=len(t)

print "Dx=%.3f y Dt=%.3f"%(Dx,Dt)
gv=arange(0,5,1)
contador=1
for g in gv:
    ######  0.hbar 1.EmC 2.EmX 3.g   4.gR
    cons=[0.65731,0.560337,500,g,4.04943]
    hbar=cons[0] #meV.ps
    EmC=cons[1] #meV.[ps^2/um^2] foton efective mass
    EmX=cons[2]
    g=cons[3]/hbar
    gR=cons[4] ##ps^-1
    #500
    a0c=hbar/EmC
    a0x=hbar/EmX
    
    F0jc=soliton(x,0)
    #F0jc=gauss(x,0,sigma=20)
    #gauss(x,0,sigma=20)
    F0jx=zeros(Nx)
    V0c=zeros(Nx)/hbar
    V0x=zeros(Nx)/hbar
    
    #F0j=array([sol(xi,0) for xi in x])
    PsiNc,PsiNx=DF2L2R(F0jc,F0jx,V0c,V0x,a0c,a0x,gR,g,x,t,puntos=10,precision=200)
    
    
    prefix0="PsiNc_DFMoL"
    prefix1="PsiNx_DFMoL"
    
    sufix="soliton_g%.3f.npy"%g
    
    try:
        os.stat(sufix)
    except:
        os.mkdir(sufix)  
    
    
    #"DF2L_DX%.1E_Dt%.1E_g%.1E.npy"%(Dx,Dt,g)
    
    save(sufix+"/x"+sufix,x)
    save(sufix+"/t"+sufix,t)
    save(sufix+"/"+prefix0+sufix,PsiNc)
    save(sufix+"/"+prefix1+sufix,PsiNx)
    pru=str(cons).replace(" ","")
    os.system("python espectro.py "+prefix0+"  "+sufix+"  "+prefix1+"  "+sufix+"  "+pru)
    print "finalizado: ",contador,"/ ",len(gv)
    contador=contador+1





#os.system("python grafica.py "+prefix0+"  "+sufix+"  "+prefix1+"  "+sufix+"  100")

#Valores para el experimento:
#
#DK=2*pi/(Dx*J)
#DE=hbar*2*pi/(Dt*N)
#kmin = -pi/Dx
#kmax=-kmin
#Emin=-hbar*pi/Dt
#Emax=-Emin
#
#imshow(abs(FTC)**0.5,cmap=cm.spectral,extent=[kmin,kmax,Emin,Emax], aspect='auto')
#def egv(eg,k):
#    return (((EmC+EmX)*hbar**2*k**2+eg*sqrt(16.0*EmC**2*EmX**2*EOmegaR**2*hbar**2+ 
#    (EmC-EmX)**2*hbar**4*k**4))/(4.0*EmC*EmX))
