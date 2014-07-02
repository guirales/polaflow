# -*- coding: utf-8 -*-
"""
Editado por mi!

Created on Thu Apr 24 11:10:28 2014

@author: Polaflow
"""
from numpy import*
from scipy import*
from DF2L import DF2L2R
import sys,os,time

print "inicia calculo"
start2 = time.time()

################################Funciones de Formateo
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
    
############################3
formato="""##resulttado del codigo MethDF2L2DRec.py 
g={0} ##ps^-1 No linealidad
Nst={1} # presicion
hbar={2}    #meV.ps
EmC={3} #meV.[ps^2/um^2] foton efective mass
EmX={4} #meV.[ps^2/um^2] exciton efective mass
gR={5}  #rad/ps^-1 rabi freq
Sig0={6} #ancho cavidad
a0c={7} #cavidad hbar/2mc
a0x={8} #exciton hbar/2mx
W0c={9}  #cavity energy
W0x={10}  #exciton energy
gammaC={11} # Decaimiento cavidad
gammaX={12} #Decaimiento exciton
C0={13} # intensidad bombeo
wl={14}  # freq laser bomb
Sigl={15} #ancho espacial bomb
DT1={16} #duracion bomb
t1={17}  #maximo bombeo
prefix0="{18}"
prefix1="{19}"
dot={20}
"""

#############################3
############################################################
#############################Funciones matematicas

def soliton(x,t):
#    x=0.3*y
    asol=0.15
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
try:
    os.stat("Datos/")
except:
    os.mkdir("Datos")  
#x=arange(-7,7,0.01)
x=linspace(-500,500,2**12)
Dx=x[1]-x[0]
Nx=len(x)

#t=arange(0,7,0.05)
t=linspace(0,100,2**10)
Dt=t[1]-t[0]
Nt=len(t)
######################Constantes
gv=[0]##no linealidades
Nst=70##precision
hbar=0.65731 #meV.ps
EmC=0.025
#0.560337 #meV.[ps^2/um^2] foton efective mass
#EmX=EmC=1
EmX=2
#500
gR=2
#4.001 ##ps^-1
Sig0=0.5
a0c=hbar/EmC
a0x=hbar/EmX
W0c=0.0
W0x=0.0
##################DECAE
gammaC=0
#0.3 #ps^-1
gammaX=0
#0.001#ps^-1
#################BOMBEO
C0=0.0
wl=0.0
Sigl=10
DT1=1.0
#0.18
t1=1.0
#######################
prefix0="PsiNc_DF"
prefix1="PsiNx_DF" 
dot=10
#F0jc=soliton(x,0)
F0jc=gauss(x,0,sigma=Sig0)

F0jx=zeros(Nx)
V0c=zeros(Nx)/hbar+W0c/hbar
V0x=zeros(Nx)/hbar+W0x/hbar

def f(t):
    return 0
#    return (C0/(sqrt(2*pi)*DT1))*gauss(x,0,sigma=Sigl)*exp(-0.5*((t-t1)/DT1)**2+1j*wl*t)
##############################
print "Dx=%.3f y Dt=%.3f"%(Dx,Dt)
for gi in gv:
        
    g=gi/hbar

    PsiNc,PsiNx=DF2L2R(F0jc,F0jx,V0c,V0x,f,a0c,a0x,gR,g,gammaC,gammaX,x,t,puntos=dot,precision=Nst)  
    sufix="Gauss1D_David_"+time.strftime("%Y_%m_%d_%H_%M_%S")    
#    sufix="Soliton1D_"+time.strftime("%Y_%m_%d_%H_%M_%S")
    path="Datos/"+sufix
 ###########################   
    try:                    #
        os.stat(path)       #
    except:                 #
        os.mkdir(path)      #
######Bloque de respaldo de constantes
    const=open(path+"/constantes.py","w")
    const.write(formato.format(g,Nst,hbar,EmC,EmX,gR,Sig0,a0c,a0x,W0c,W0x,gammaC,gammaX,C0,wl,Sigl,DT1,t1,prefix0,prefix1,dot))
    const.close()    
###################    
    save(path+"/x"+sufix+".npy",x)
    save(path+"/t"+sufix+".npy",t)
    save(path+"/"+prefix0+sufix+".npy",PsiNc)
    save(path+"/"+prefix1+sufix+".npy",PsiNx)
print "\n Terminado en "+hms(time.time()-start2)
