# -*- coding: utf-8 -*-
"""
Created on Thu May 29 04:56:07 2014

@author: memo
"""

from numpy import*
from scipy import*
from DF2L2D import DF2L2R2D
from mpl_toolkits.mplot3d import Axes3D
import sys,os,time


print "inicia calculo"
start2 = time.time()

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

def soliton(x,y,t):    
    def soliton1D(x,t,asol,vsol,x0sol):
        #asol=0.4
        #vsol=0.0
        #x0sol=0.0
        return asol*exp(1.0j*(vsol*x+(asol**2-vsol**2)*t/2))/cosh(asol*(x-vsol*t-x0sol))

    return soliton1D(X,t,0.4,0,0)*soliton1D(Y,t,0.4,0,0)



def gauss(x,y,t,sigma=15.0):
    return exp(-0.5*(x**2+y**2)/sigma**2)/sqrt(pi*sigma**2)

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
w={14}  # freq laser bomb
Sig={15} #ancho espacial bomb
DT1={16} #duracion bomb
t1={17}  #maximo bombeo
prefix0="{18}"
prefix1="{19}"
"""

#############################3
#x=arange(-7,7,0.01)
x=linspace(-60,60,100)
Dx=x[1]-x[0]
Nx=len(x)

X,Y=meshgrid(x,x)

#t=arange(0,7,0.05)
t=linspace(0,20,80)
Dt=t[1]-t[0]
Nt=len(t)
############################33
gv=array([0])##no linealidad
Nst=80##precision
hbar=0.65731 #meV.ps
EmC=0.560337 #meV.[ps^2/um^2] foton efective mass
EmX=500
gR=4.001 ##ps^-1
Sig0=20
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
C0=5.0
wl=100.0
Sigl=20
DT1=0.07
#0.18
t1=1.48
#######################
prefix0="PsiNc_DFMoL"
prefix1="PsiNx_DFMoL" 

F0jc=gauss(X,Y,0,sigma=Sig0)
#F0jc=soliton(X,Y,0)
F0jx=zeros((Nx,Nx))
V0c=zeros((Nx,Nx))/hbar+W0c/hbar
V0x=zeros((Nx,Nx))/hbar+W0x/hbar
def f(t):
    return (C0/(sqrt(2*pi)*DT1))*gauss(X,Y,0,sigma=Sigl)*exp(-0.5*((t-t1)/DT1)**2+1j*wl*t)

print "Dx=%.3f y Dt=%.3f"%(Dx,Dt)
#gv=arange(10,,1)
contador=1
for gi in gv:

    g=gi/hbar

    PsiNc,PsiNx=DF2L2R2D(F0jx,F0jx,V0c,V0x,f,a0c,a0x,gR,g,gammaC,gammaX,x,t,puntos=10,precision=Nst)   
#    sufix="Guirales_m02_pulseA831_Nar_g%.3f.npy"%g
    sufix="Gauss2D_gR%.3g_g%.3g_gamm(%.3g,%.3g)_DT1%.3g.npy"%(gR,g,gammaC,gammaX,DT1)
    path="Datos/"+sufix
    
    try:
        os.stat(path)
    except:
        os.mkdir(path)  
######Bloque de respaldo de constantes
    const=open(path+"/constantes.py","w")
    const.write(formato.format(g,Nst,hbar,EmC,EmX,gR,Sig0,a0c,a0x,W0c,W0x,gammaC,gammaX,C0,wl,Sigl,DT1,t1,prefix0,prefix1))
    const.close()    
###################
    
    #"DF2L_DX%.1E_Dt%.1E_g%.1E.npy"%(Dx,Dt,g)
    
    save(path+"/x"+sufix,x)
    save(path+"/t"+sufix,t)
    save(path+"/"+prefix0+sufix,PsiNc)
    save(path+"/"+prefix1+sufix,PsiNx)
print "\n Terminado en "+hms(time.time()-start2)
#    pru=str(cons).replace(" ","")
#    os.system("python espectro.py "+prefix0+"  "+sufix+"  "+prefix1+"  "+sufix+"  "+pru)
#    print "finalizado: ",contador,"/ ",len(gv)
#    contador=contador+1





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
