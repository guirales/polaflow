# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 14:47:41 2014
@author: Polaflow
"""
from __future__ import division
from numpy import exp,append,array,bmat,dot,conjugate,sqrt,concatenate,ones,arange,log,zeros
from scipy.misc import pade
from scipy.linalg import inv,solve
from scipy.interpolate import approximate_taylor_polynomial  
import sys

##begin Shao-Wang-toyama method
def ShaoWang(Psi0c,Psi0x,V0c,V0x,f,mc,mx,Om,gx,gammaC,gammaX,x,t,R=10,M=10,hbar=0.6582,sweeps=2,tolerance=1E-2,abc=0.3,apod=1.0E-3):
    print "using the improved version  ShaoWang of ShaoWangV3ABC final apodized" 
    ##Begin formatting purpose functions
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
    
    ##End formatting purpose functions
    J = len(x)
    T = len(t)
    dx = x[1]-x[0]
    dt = t[1]-t[0]
    print "Dx: %.3f, Dt: %.3f, abc: %.2f"%(dx,dt,abc)
    ##############################################################################
    #Functions 
    ##############################################################################
    #begin function to calculate the ak and ck coefficients
    def H(i,j):
            if j<=(R-1):
                r=i+1
                k=j+1
                MAB=2*k**(2*r-2)
                return MAB
            else:
                r=i+1
                k=j-R+1
                MAB=2*k**(2*r-2)
                return MAB*k**2/((2*r)*(2*r-1)) 
    #end function to calculate the ak and ck coefficients
    #begin process to calculate the zs coefficients
    def zs(m):
        e_exp=approximate_taylor_polynomial(exp,0,2*m,0)
        p, q = pade(e_exp.c[::-1], m)
        return (p.r)[::-1]    
    #########################################
    zs=zs(M)    
    #end  process to calculates zs coefficients 
    #begin process to calculates the ak coefficients 
    Ar=array([[H(i,j) for j in range(2*R)] for i in range(2*R)])
    rl=-1*inv(Ar)[:,0]
    a=append(1,rl[:R])
    c=append(-2*sum(rl[R:]),rl[R:])
    #end process to calculate the ak coefficients   
    #begin functions that  evaluates the matrix element ij of A(s) for the s-th iteration
    def hc(i,j,zs):
        if abs(j-i) <= R:
            return a[abs(j-i)]*(1.0-1j*dt*(V0c[j])/zs)-1j*hbar*dt*c[abs(j-i)]/(2.0*mc*zs*dx**2)
        else:
            return 0.0           
    def hx(i,j,zs,psinl):
        if abs(j-i) <= R:
            return a[abs(j-i)]*(1.0-1j*dt*(V0x[j]+gx*abs(psinl[j])**2.0)/zs)-1j*hbar*dt*c[abs(j-i)]/(2.0*mx*zs*dx**2)
        else:
            return 0.0           
    def hD(i,j,zs):
        if abs(j-i) <= R:
            return a[abs(j-i)]*(-1j*dt*Om)/zs
        else:
            return 0.0                    
    # This function build up the matrix A(s) for the s-th iteration    
    def AC(zs):
        return array([[hc(i,j,zs) for j in range(J)] for i in range(J)])
        
    def AX(zs,psinl):
        return array([[hx(i,j,zs,psinl) for j in range(J)] for i in range(J)])
        
    def A2(zs):
        return array([[hD(i,j,zs) for j in range(J)] for i in range(J)])
    #end functions that  evaluates the matrix element ij of A(s) for the s-th iteration         
    psiin=concatenate((Psi0x,Psi0c))
    phix1=Psi0x[1:2*R+1]
    phix2=Psi0x[-2*R-1:-1]
    phic1=Psi0c[1:2*R+1]
    phic2=Psi0c[-2*R-1:-1]
    Intx0=sum(abs(Psi0x)**2)*dx
    Intc0=sum(abs(Psi0c)**2)*dx
    ##This control the division by zero
    if Intx0==0.0:
        Intx0=1
    if Intc0==0.0:
        Intc0=1
    
    Absor=exp(log(abc)*arange(2*R)/R)
    dissipa=-concatenate((dt*gammaX*ones(J),dt*gammaC*ones(J)))
    psinl=psiin
    MT=[]
    MT.append(psiin)  
    tn=1
    nl=int((T+1)/10)
    print "["+(10*"#")+"]"
    sys.stdout.write("[")
    ##begin s-matrices product  
    for tn in range(1,T):
        crit=1
        crit3=1.0        
        while (crit<sweeps) and (crit3>=tolerance):
            psiin1=psiin
            for i in range(M):
                A=AC(zs[i])
                D=AX(zs[i],psiin[:J])
                D2=AX(zs[i],psinl[:J])
                C=A2(zs[i])
                M1=(bmat([[D,C],[C,A]])).getA()
                M2=(bmat([[D2,C],[C,A]])).getA()
                b=dot(conjugate(M1),psiin1) 
                psiin1=solve(M2,b,overwrite_a=True, overwrite_b=True)
                
                
                
            crit3=sqrt(sum(abs(psiin1-psinl)**2.0))
            psinl=psiin1
            crit=crit+1
########begin ABCD boundaries process 
        #getting the boundaries
        phix1B=psiin1[1:2*R+1]
        phix2B=psiin1[J-2*R-1:J-1]
        phic1B=psiin1[J+1:J+2*R+1]
        phic2B=psiin1[-2*R-1:-1] 
        #putting the absorbing        
        psiin1[:2*R]=Absor[::-1]*phix1
        psiin1[J-2*R:J]=Absor*phix2
        psiin1[J:J+2*R]=Absor[::-1]*phic1
        psiin1[-2*R:]=Absor*phic2 
        #saving the absorbing boundaries
        phix1=phix1B
        phix2=phix2B
        phic1=phic1B
        phic2=phic2B
        #End ABCD process    
######Apodization
        Intxn=sum(abs(psiin1[:J])**2)*dx/Intx0 #exc
        Intcn=sum(abs(psiin1[J:])**2)*dx/Intc0 #cav
        #print "apod ",apod,"intxn ",Intxn," intx0 ", Intx0,"intcn ",Intcn," intc0 ", Intc0
        if Intxn<=apod:
            psiin1[:J]=zeros(J)
        if Intcn<=apod:
            psiin1[J:]=zeros(J)
        
###########New wavefunction
        psiin=psiin1
        psinl=psiin1
        psiin+=dissipa*psiin
        MT.append(psiin)        
        if (tn%nl==0)&(tn!=0): 
            sys.stdout.write("#")
    sys.stdout.write("]")
    MT=array(MT)   
    return MT[:,J:],MT[:,:J],crit3   ##Photon, Exciton 
