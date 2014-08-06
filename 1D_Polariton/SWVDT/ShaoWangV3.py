# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 14:47:41 2014

@author: Polaflow
"""
from __future__ import division
from numpy import exp,append,array,bmat,dot,conjugate,sqrt,concatenate
from scipy.misc import factorial,pade
#from numpy.linalg.linalg import inv, lu
from scipy.linalg import inv,solve

from scipy.interpolate import approximate_taylor_polynomial  
import sys
 
def ShaoWang(Psi0c,Psi0x,V0c,V0x,f,mc,mx,Om,gx,gammaX,gammaC,x,t,R=10,M=10,hbar=0.6582):
    print "using the improved version  ShaoWang of ShaoWangV3" 
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
    
    
    J = len(x)
    T = len(t)
    dx = x[1]-x[0]
    dt = t[1]-t[0]  
    ##############################################################################
    #Functions 
    ##############################################################################
    # to calculate the ak and ck coefficients
    # to calculate the ak and ck coefficients
    def H(i,j):        
            if j<=(R-1):
                r=i+1
                k=j+1
                A_rk=2.0*k**(2*r-2)/factorial(2.0*r-2)
                return A_rk 
            else:
                r=i+1
                k=j-R+1
                A_rk=2.0*k**(2*r-2)/factorial(2.0*r-2)
                return k**2*A_rk/(2.0*r*(2.0*r-1.0))
    
    # to calculate the zs coefficients
    def zs(m):
        e_exp=approximate_taylor_polynomial(exp,0,2*m,0)
        p, q = pade(e_exp.c[::-1], m)
        return (p.r)[::-1]    
    #########################################
    zs=zs(M)    
    # Calculates the ck coefficients 
    # Calculates the ak coefficients 
    Ar=array([[H(i,j) for j in range(2*R)] for i in range(2*R)])
    rl=-inv(Ar)[:,0]
    a=append(1,rl[:R])
    c=append(-2*sum(rl[R:]),rl[R:])
    
    
    # This function evaluates the matrix element ij of A(s) for the s-th iteration
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
        
    #This function uses the iterative approach to obtain the propagation matrix KT    
#    def KT(zs,psinl,psiendant):
#        for i in range(M):
#            A=AC(zs[i])
#            D=AX(zs[i],psiendant)
#            D2=AX(zs[i],psinl)
#            C=A2(zs[i])
#            M1=(bmat([[D,C],[C,A]])).getA()
#            M2=(bmat([[D2,C],[C,A]])).getA()
##            P,L,U=lu(M2)            
##            invM2=dot(inv(U),dot(inv(L),P))
#            AA=dot(inv(M2),conjugate(M1))
#            if i<1:
#                Kt=AA;
#            else:
#                Kt=dot(AA,Kt)
#        return Kt
         
    psiin=concatenate((Psi0x,Psi0c))
    psinl=psiin
    MT=[]
    MT.append(psiin)  
    tn=1
    nl=int((T+1)/10)
    print "["+(10*"#")+"]"
    sys.stdout.write("[")
#    start = time.time()    
#    while tn<=T:
    out=open("STOUT.dat","w")
    out.write("Data out ")
    out.close()
    for tn in range(1,T):
        crit=1
        while crit<=1:
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
            out=open("STOUT.dat","a")
            psinl=psiin1
            out.write("tn %s, crit %s, i  %s, crit3  %s \n"%(tn,crit,i,crit3))
            out.close()
#            print(tn," , ",crit3)
            crit+=1
        
        psiin=psiin1
        MT.append(psiin)
#        tn+=1
        if (tn%nl==0)&(tn!=0): 
            sys.stdout.write("#")
    
    
    sys.stdout.write("]")
    MT=array(MT)   
    
#    MT1=array([[MT[i,j] for j in range(J+1)] for i in range(T)])
#    MT2=array([[MT[i,j] for j in range(J+1,2*J+2)] for i in range(T)])
    return MT[:,J:],MT[:,:J]   ##Cavidad, Exciton 
            
            