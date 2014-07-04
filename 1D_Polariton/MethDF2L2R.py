"""
Created on Thu Apr 24 11:10:28 2014

@author: Polaflow
"""
from numpy import*
from scipy import*
from DF2L import DF2L2R
import os
import time
import inspect



print "Calculating...."
start2 = time.time()

################################Funciones de Formateo
def LineNo():
    return inspect.currentframe().f_back.f_lineno

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

####################################
#load("Datos/PsiCD.npy") 
#load("Datos/PsiXD.npy")
########################
try:
    os.stat("Data/")
except:
    os.mkdir("Data") 
Li=LineNo()
###########################################################################
########################ZERO FUNCTION######################################
###########################################################################
##For some special pattern you must modify the followin lines##############
###Cavity pattern:                                                      ###
def PsiC0(x):
    sigma=20                                                           ###
    return  exp(-0.5*x**2/sigma**2)/sqrt(sqrt(pi*sigma**2))                                    ###
##Exciton patern:                                                       ###
def PsiX0(x):                                                           ###
    return 0.0*x                                      ###
#Pumping (sourcing) pattern:                                            ###
def f(t):                                                               ###
    return 0.0*t                                                        ###
##Cavity potential pattern:                                             ###
def V0C(x):                                                             ###
    return 0.0*x+W0c/hbar                                               ###
##Exciton potential pattern:                                            ###
def V0X(x):                                                      ###
    return 0.0*x+W0x/hbar                                               ###
###########################################################################
################CORE OF PARAMETERS#########################################
###########################################################################
##SPACE GRID                                                            ###
x=linspace(-70,70,2**9)                                              ###
Dx=x[1]-x[0]                                                            ###
Nx=len(x)                                                               ###
##TIME GRID                                                             ###
t=linspace(0,29.3,300.3)                                                ###
Dt=t[1]-t[0]                                                            ###
Nt=len(t)                                                               ###
####################################################CONSTANS###############
gv=[0]##NO LIENALITY TERM                                               ###
Nst=200##TEMPORAL PRECISION                                             ###
dot=10##ESPATIAL PRECISION, DOTS TAKEN FOR DERIVATIVE                   ###
hbar=0.65731 #meV.ps                                                    ###
EmC= 0.560337# 0.025 #meV.[ps^2/um^2] foton efective mass                ###
EmX=500# 2 meV.[ps^2/um^2] exciton efective mass                         ###
gR=4.001# 2##ps^-1  rabii frequency                                     ###
a0c=hbar/EmC                                                            ###
a0x=hbar/EmX                                                            ###
####################################################DEPHASING##############
W0c=0.0 #depahasing cavity term                                         ###
W0x=0.0 #depahasing exciton term                                        ###
####################################################DECAY##################
gammaC=0 #0.3 #ps^-1 cavity decay                                       ###
gammaX=0 #0.001#ps^-1 exciton decay                                     ###
###########################################################################
#######################NAMES AND PLACE ADJUST
prefix0="PsiNc_DF"
prefix1="PsiNx_DF"
Lf=LineNo()
sufix="Gauss1D"+time.strftime("%Y_%m_%d_%H_%M_%S")    
#sufix="Soliton1D_"+time.strftime("%Y_%m_%d_%H_%M_%S")
path="Data/"+sufix 

#F0jc=soliton(x,0)
#F0jc=gauss(x,0,sigma=Sig0)
#F0jx=zeros(Nx)
###############################FUNCTIONS AND POTENTIALS INITIALIZATION
F0jc=PsiC0(x)
#[len(x)/2-400:len(x)/2+400]
F0jx=PsiX0(x)
#[len(x)/2-400:len(x)/2+400]
V0c=V0C(x)
#[len(x)/2-400:len(x)/2+400]
V0x=V0X(x)
#[len(x)/2-400:len(x)/2+400]
##############################
print "Dx=%.3f y Dt=%.3f"%(Dx,Dt)
for gi in gv: 
    g=gi/hbar
    PsiNc,PsiNx=DF2L2R(F0jc,F0jx,V0c,V0x,f,a0c,a0x,gR,g,gammaC,gammaX,x,t,puntos=dot,precision=Nst)  
 ###########################   
    try:                    #
        os.stat(path)       #
    except:                 #
        os.mkdir(path)      #
######Bloque de respaldo de constantes
    fpy=open(inspect.getfile(inspect.currentframe()),"r")
    #fpy=open("MethDF2L2RNV.py","r")
    content=fpy.readlines()
    fpy.close()
    const=open(path+"/constants.py","w")
    const.write("from numpy import *\n")
    const.write("from scipy import *\n")
    const.write("g="+str(g)+"\n")
    for li in content[Li:Lf-1]:
        const.write(li)    
    const.close()    
###################    
    save(path+"/x"+sufix+".npy",x)
    save(path+"/t"+sufix+".npy",t)
    save(path+"/"+prefix0+sufix+".npy",PsiNc)
    save(path+"/"+prefix1+sufix+".npy",PsiNx)
print "\n Ended in  "+hms(time.time()-start2)