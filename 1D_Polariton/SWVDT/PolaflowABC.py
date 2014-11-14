# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 18:46:46 2014
##prub
@author: Polaflow
"""
from __future__ import division

class Particle():
    from numpy import zeros,exp,log,arange    
    hbar=0.65731#meV.ps 
    rabi=2.62989731    ##ps^-1  rabii frequency
    x=zeros(50)
    t=zeros(50)
    name="Psi"
    mass=2.0    #meV.[ps^2/um^2] efective mass
    w0=0
    gamma=0.0
    g=0
    R=10
    abc=0.3    
    Psi0=zeros(len(x))
    Potential=zeros(len(x))
    Psi=zeros(len(x))
    Absor=exp(log(abc)*arange(2*R)/R)    
    
    
    def pump(self,x,t):
        return 0*x*t
    
    formato="""##code polaflow.py 
g={0} ## gi/hbar NOn LINEAR TERM                                               
hbar={1} #meV.ps                                                    
mass={2} #meV.[ps^2/um^2] efective mass                                       
rabi={3}# 4.001##ps^-1  rabii frequency 
################################################DECAY##################
gamma={4} #0.3 #ps^-1 decay                                        
#######################################################################
###################NAMES  ADJUST
name="{5}"
xi={6}
xf={7}
ti={8}
tf={9}
###################OTHERS ADJUSTABLE PARAMETERS
#all them must be spatial functions, except f that must be a (explicit) temporal function
#Psi0C
#Psi0X
#V0c
#V0x
#f
"""
    def Summary(self):
        print self.formato.format(self.g,self.hbar,self.mass,self.rabi,self.gamma,self.name,self.x[0],
                             self.x[-1],self.t[0],self.t[-1])
    def Gsummary(self):
        from numpy import abs
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2,2,figsize=(6,6))
        full_txt=r"""Name: %s
mass: %s
hbar: %s
rabi: %s
gamma: %s
w0: %s
g: %s"""%(self.name,self.mass,self.hbar,self.rabi,self.gamma,self.w0,self.g)
#        
        axes[0,0].set_title(self.name+"_Psi0", fontsize=14, weight='bold')
        axes[0,0].plot(self.x,abs(self.Psi0)**2)
        
        axes[0,1].text( 0.1, 0.1,full_txt,fontsize=14,ha='left')
        axes[0,1].set_xticklabels("", visible=False)
        axes[0,1].set_yticklabels("", visible=False)
#        
        axes[1,0].set_title("Potential", fontsize=14, weight='bold')
        axes[1,0].plot(self.x,self.Potential)
#        
        axes[1,1].set_title("Pump", fontsize=14, weight='bold')
        axes[1,1].plot(self.x,self.pump(self.x,0))
        
        return fig,axes
    
    def indices(self,lista,xa,xb):
        if xa>=lista[0]and xb<=lista[-1]:
            Nt=len(lista)
            L=lista[-1]-lista[0]
            na=(xa-lista[0])/L
            nb=(xb-lista[0])/L
            return int(Nt*na),int((Nt-1)*nb)
        else:
            print "index graph out of range"
            return 0,-1
            
    def PlotXT(self,contrast=2,xa=0,xb=-1,ta=0,tb=-1,**kargs):
        
        import matplotlib.pyplot as plt
        
        if xa==0:
            xa=self.x[0]    
        if xb==-1:
            xb=self.x[-1]        
        if ta==0:
            ta=self.t[0]
        if tb==-1:
            tb=self.t[-1]
        
        fig, axes = plt.subplots(1,1,**kargs)
        
        xia,xib=self.indices(self.x,xa,xb)
        tia,tib=self.indices(self.t,ta,tb)
        axes.set_title(self.name)        
        axes.imshow(abs(self.Psi[tia:tib,xia:xib])**contrast,extent=[self.x[xia],self.x[xib],self.t[tia],self.t[tib]],aspect='auto',origin='lower')
        return fig,axes
###########################################

########Polaflow1D                             

class Polaflow1D(object):
    def __init__(self,pho=None,exc=None):
        from numpy import zeros
        self.abc=0.3
        self.pho=pho
        self.exc=exc
        self.x=zeros(50)                                              ###
        self.Dx=self.x[1]-self.x[0]                                                            ###
        self.J=len(self.x)                                                               ###
    ##TIME GRID                                                             ###
        self.t=zeros(50)                                                ###
        self.Dt=self.t[1]-self.t[0]                                                            ###
        self.Nt=len(self.t)                                                               ###
    ####################################################CONSTANS###############
        #gi=0##NOn LINEAR TERM                                               ###
        self.M=10##TEMPORAL PRECISION (Pade)                                             ###
        self.R=10##ESPATIAL PRECISION, DOTS TAKEN FOR DERIVATIVE                   ###
        self.hbar=0.65731 #meV.ps                                                    ###
        self.rabi=2.62989731#4.001##ps^-1  rabii frequency                
        self.name='Polaflow'
        self.pathwd="DataPF/"
        self.path="DataPF/"
        self.error=(0,0)
        self.tolerance=1E-2
        self.sweeps=2
        self.apod=1.0E-3
###############################
#<CS############################
    def CS(self,num):
        if abs(num)<10:
            if num<0:
                txt="-0%d"%int(abs(num))
            else:
                txt="0%d"%int(num)
        else:
            txt="%d"%int(num)
        return txt
#CS>##########################
#<hms#########################
    def hms(self,tf):
        hor=(int(tf/3600))  
        minu=int((tf-(hor*3600))/60)  
        seg=int(tf-((hor*3600)+(minu*60)))  
        tfx=self.CS(hor)+"h"+self.CS(minu)+"m"+self.CS(seg)+"s|"
        return tfx   
#hms>#####################
#<Refresh###################
    def Refresh(self):
        self.Dx=self.x[1]-self.x[0]                                                            ###
        self.J=len(self.x)
        self.Dt=self.t[1]-self.t[0]                                                            ###
        self.Nt=len(self.t)
        self.pho.R=self.exc.R=self.R
        self.pho.abc=self.exc.abc=self.abc
        self.pho.x=self.exc.x=self.x
        self.pho.t=self.exc.t=self.t
        self.pho.rabi=self.exc.rabi=self.rabi
        self.pho.hbar=self.exc.hbar=self.hbar
#Refresh>######################
#<Save##########################            
    def Save(self,FileNameAuto=True):
        import time
        import os
        from numpy import save#,meshgrid
        if FileNameAuto:
            suffix=str(self.name)+time.strftime("_%Y_%m_%d_%H_%M_%S")
        else:
            suffix=str(self.name)
            
        self.path=self.pathwd+suffix
        formato="""
namex="{0}"
mx={1}
gammax={2}
gx={3}
w0x={4}
namec="{5}"
mc={6}
gammac={7}
gc={8}
w0c={9}
namep="{10}"
hbar={11}
rabi={12}
M={13}
R={14}
suffix="{15}"
sweeps={16}
abc={17}
apod={18}
        """  
        try:                    #
            os.stat(self.pathwd)       #
        except:                 #
            os.mkdir(self.pathwd)      #
            
        os.mkdir(self.path)
        const=open(self.path+"/constants.py","w")
        const.write(formato.format(self.exc.name,self.exc.mass,self.exc.gamma,self.exc.g,self.exc.w0,
                                   self.pho.name,self.pho.mass,self.pho.gamma,self.pho.g,self.pho.w0,
                                   self.name, self.hbar,self.rabi,self.M,self.R,suffix,self.sweeps,
                                   self.abc,self.apod))
        const.close()
        save(self.path+"/x"+suffix+".npy",self.x)
        save(self.path+"/t"+suffix+".npy",self.t)
        save(self.path+"/"+self.pho.name+suffix+".npy",self.pho.Psi)
        save(self.path+"/"+self.exc.name+suffix+".npy",self.exc.Psi)
        save(self.path+"/V0C"+suffix+".npy",self.pho.Potential)
        save(self.path+"/V0X"+suffix+".npy",self.exc.Potential)
        save(self.path+"/pumpx"+suffix+".npy",self.pho.pump(self.x,0))
        print "Saved in ",self.path
#Save>######################
#############################        
    def Calculate(self):
        from ShaoWangV3ABC import ShaoWang
        from numpy import shape
        import time
        start2 = time.time()
        V0C=self.pho.Potential/self.hbar+self.pho.w0
        V0X=self.exc.Potential/self.hbar+self.exc.w0
        EmC=self.pho.mass
        EmX=self.exc.mass
        
        if (shape(self.pho.Psi0)==shape(self.exc.Psi0) and shape(self.pho.Psi0)==shape(self.x) 
            and shape(self.pho.Psi0)==shape(self.pho.Potential)and shape(self.pho.Psi0)==shape(self.exc.Potential)):
                
            self.pho.Psi,self.exc.Psi,self.error=ShaoWang(self.pho.Psi0,self.exc.Psi0,V0C,V0X,self.pho.pump,EmC,
                                             EmX,self.rabi,self.exc.g,self.pho.gamma,self.exc.gamma,
                                             self.x,self.t,R=self.R,M=self.M,hbar=self.hbar,
                                             sweeps=self.sweeps,tolerance=self.tolerance,abc=self.abc,apod=self.apod)

        else:
            print "Warning! The shape of all vectors of functions, space, and potentials must be the same. Correct this and try again "

        print "\n Finished in  "+self.hms(time.time()-start2)
        
        
###################################      
####################################
    
    def Load(self):
        from numpy import load
        import os, imp
        wdir=os.getcwd()
        os.chdir(self.path)
        import constants as ct
        imp.reload(ct)
        self.exc.name=ct.namex
        self.exc.mass=ct.mx
        self.exc.gamma=ct.gammax
        self.exc.g=ct.gx
        self.exc.w0=ct.w0x
        self.pho.name=ct.namec
        self.pho.mass=ct.mc
        self.pho.gamma=ct.gammac
        self.pho.g=ct.gc
        self.pho.w0=ct.w0c
        self.name=ct.namep
        self.hbar=ct.hbar
        self.rabi=ct.rabi
        self.M=ct.M
        self.R=ct.R
        suffix=ct.suffix
        self.sweeps=ct.sweeps
        self.abc=ct.abc
        self.apod=ct.apod
        
        self.x=load("x"+suffix+".npy")
        self.t=load("t"+suffix+".npy")
        self.pho.Psi=load(self.pho.name+suffix+".npy")
        self.pho.Psi0=self.pho.Psi[0]
        self.pho.Potential=load("V0C"+suffix+".npy")
        self.pumpx=load("pumpx"+suffix+".npy")        
        self.exc.Psi=load(self.exc.name+suffix+".npy")
        self.exc.Psi0=self.exc.Psi[0]
        self.exc.Potential=load("V0X"+suffix+".npy")
        os.chdir(wdir)
        
        self.pho.Psi0=self.pho.Psi[0]
        self.exc.Psi0=self.exc.Psi[0]
        
        self.Dx=self.x[1]-self.x[0]                                                            ###
        self.J=len(self.x)                                                               ###                                              ###
        self.Dt=self.t[1]-self.t[0]                                                            ###
        self.Nt=len(self.t)        
        os.chdir(wdir)
        self.Refresh()
        print "Data %s loaded..."%self.name
####################################        
#############################
    def PlotKE(self,contrast=2):
        from numpy import pi,arange, sqrt,flipud
        from scipy.fftpack import fft2,fftshift
        import matplotlib.pyplot as plt
        
        DK=2*pi/(self.Dx*self.J)
        kf = pi/self.Dx
        kl=arange(-kf,kf,DK)
        DE=self.hbar*2*pi/(self.Dt*self.Nt)
        Ef=self.hbar*pi/self.Dt
        En=arange(-Ef,Ef,DE)
        #
        AFT1=flipud(abs(fftshift(fft2(self.pho.Psi))))
        AFT2=flipud(abs(fftshift(fft2(self.exc.Psi))))          
        
        def egv(eg,k):
            EmC=self.pho.mass
            EmX=self.exc.mass
            return (((EmC+EmX)*self.hbar**2*k**2+eg*sqrt(16.0*EmC**2*EmX**2*self.rabi**2*self.hbar**2+ 
                    (EmC-EmX)**2*self.hbar**4*k**4))/(4.0*EmC*EmX))
        
        fig, axes = plt.subplots(1,2)
        axes[0].set_title(self.pho.name)
        axes[0].imshow(abs(AFT1)**contrast,extent=[kl[0],kl[-1],En[0]+DE,En[-1]+DE],aspect='auto',origin='lower')
        axes[0].plot(kl,egv(1,kl),"black")
        axes[0].plot(kl,egv(-1,kl),"black")
        axes[0].set_ylim(En[0],En[-1])
        axes[1].set_title(self.exc.name)
        axes[1].imshow(abs(AFT2)**contrast,extent=[kl[0],kl[-1],En[0]+DE,En[-1]+DE],aspect='auto',origin='lower')
        axes[1].plot(kl,egv(1,kl),"black")
        axes[1].plot(kl,egv(-1,kl),"black");
        axes[1].set_ylim(En[0],En[-1])
        return fig, axes
    
    def AnimateXT(self):
        from numpy import array
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation
        import sys#,os,time

        
        pl1=abs(self.pho.Psi)**2
        pl2=abs(self.exc.Psi)**2
        maxi=(array([pl1[0].max(),pl2[0].max()])).max()
        
        fig2, axes = plt.subplots(1,2)
        axes[0].plot(self.x,pl1[0])
        axes[1].plot(self.x,pl2[0])
        
        nl=int(self.Nt/10)
        
        print "[__h__m__s|"+(10*"##")+"]"
        sys.stdout.write("[__h__m__s|")
        #start = time.time()
        
        def data_gen2(fr):
            
            lt="%d/%d"%(fr,self.Nt)
            axes[0].clear()
            axes[1].clear()
            
            axes[0].set_title(self.pho.name+" "+lt)
            axes[0].set_ylim(0,maxi)
            axes[0].plot(self.x,pl1[fr])
            
            axes[1].set_title(self.exc.name+" "+lt) 
            axes[1].set_ylim(0,maxi)
            axes[1].plot(self.x,pl2[fr])
        #    axes[1].plot(x,abs(soliton(x,Dt*fr))**2)
            if (fr%nl==0)&(fr!=0):
                sys.stdout.write("##")
        
        #    axes[1].set_ylim(0,maxi)
        #picosC=array([],dtype=int)
        ani = animation.FuncAnimation(fig2, data_gen2, interval=100,frames=self.Nt, blit=False)
        
        ani.save(self.name+'ABCtemporal.mp4',writer="avconv", codec="libx264")
        sys.stdout.write("]")
        plt.close(fig2)
        ##############################
        from IPython.display import HTML
        
        video = open(self.name+'ABCtemporal.mp4', "rb").read()
        video_encoded = video.encode("base64")
        video_tag = '<video controls alt="test" src="data:video/x-m4v;base64,{0}">'.format(video_encoded)
        ## width="600" height="400" 
        return HTML(video_tag)

   
   
   
##########################################
##################Next Upgrades, don't uncomment:
#    
#    class ListTable(list):
#        """ Overridden list class which takes a 2-dimensional list of 
#            the form [[1,2,3],[4,5,6]], and renders an HTML Table in 
#            IPython Notebook. """
#            
#        def _repr_html_(self):
#            html = ["<table>"]
#            for row in self:
#                html.append("<tr>")
#                
#                for col in row:
#                    html.append("<td>{0}</td>".format(col))
#                
#                html.append("</tr>")
#            html.append("</table>")
#            return ''.join(html)
#    
#    
#    def Ls(self):
#                
#        os.chdir("Data")
#        directorios=filter(os.path.isdir, os.listdir(os.getcwd()))
#        directorios.sort()
#        #directorios=os.listdir(os.getcwd()+"/Datos")
#        os.chdir(wdir)
#        #list(enumerate(directorios))
#        tabla=ListTable()
#        tabla.append(['i','$g/\hbar$','Nst','$\hbar$','$EmC$','$EmX$','$gR$','$t_i$','$t_f$','$\Delta t$','$x_i$','$x_f$','$\Delta x$','$a_{0c}$','$a_{0x}$','$\omega_{0c}$','$\omega_{0x}$','$\gamma_C$','$\gamma_X$','nombre'])
#        for i in range(len(directorios)):
#            os.chdir("Data/"+directorios[i])
#            (g,Nst,hbar,EmC,EmX,gR,ti,tf,Dt,xi,xf,Dx,a0c,a0x,W0c,W0x,gammaC,gammaX)=18*('X',)
#            import  constants.py as ct
#            xi=x[0]
#            xf=x[-1]
#            ti=t[0]
#            tf=t[-1]
#            filaT=map(lambda x:x if str(type(x))=="<type 'str'>" else x if str(type(x))=="<type 'int'>" else round(x,3) ,
#                      [i,g,Nst,hbar,EmC,EmX,gR,ti,tf,Dt,xi,xf,Dx,a0c,a0x,W0c,W0x,gammaC,gammaX,directorios[i]])
#            tabla.append(filaT)
#            #%cd -0
#            os.chdir(wdir)
#        
#        tabla
#    
