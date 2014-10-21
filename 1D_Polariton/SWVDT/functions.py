# -*- coding: utf-8 -*-
"""
Created on Sat Jul  5 16:26:02 2014

@author: fabrice
"""
from __future__ import division

def egv(eg,k,EmX=0.6581,EmC=1.0E5,EOmegaR=2.66947,hbar=6.582E-1):
    from numpy import  sqrt
    return (((EmC+EmX)*hbar**2*k**2+eg*sqrt(16.0*EmC**2*EmX**2*EOmegaR**2+ 
        (EmC-EmX)**2*hbar**4*k**4))/(4.0*EmC*EmX))

def polasolitonD(x,Dx,Rabi=4,hbar=0.65731,mC=0.5,mX=500,k0=0,Delta=0,Sigma=0.5):
    from numpy import sqrt,fft,arange,exp,pi,array

    def polC(k):
        return -(((-2*mC*mX*Delta + k**2*(mC - mX)*hbar + sqrt(4*mC**2*mX**2*(Delta**2 + 4*Rabi**2) + 
                4*k**2*mC*mX*(-mC + mX)*Delta*hbar + k**4*(mC - mX)**2*hbar**2))*
                sqrt((16*Sigma + Sigma*abs((2*mC*mX*Delta + k**2*(-mC + mX)*hbar - 
                sqrt(4*mC**2*mX**2*(Delta**2 + 4*Rabi**2) + 4*k**2*mC*mX*(-mC + mX)*Delta*
                hbar + k**4*(mC - mX)**2*hbar**2))/(mC*mX*Rabi))**2)**(-1)))/
                (exp((k - k0)**2/(2*Sigma**2))*mC*mX*pi**(1/4)*Rabi))
    def polX(k):
        return ((4*sqrt((16*Sigma + Sigma*abs((2*mC*mX*Delta + k**2*(-mC + mX)*hbar - 
               sqrt(4*mC**2*mX**2*(Delta**2 + 4*Rabi**2) + 4*k**2*mC*mX*(-mC + mX)*Delta*
               hbar + k**4*(mC - mX)**2*hbar**2))/(mC*mX*Rabi))**2)**(-1)))/
               (exp((k - k0)**2/(2*Sigma**2))*pi**(1/4)))
    Nxp = sqrt(2*pi*len(x))/Dx
    Dkp = 2*pi/(Dx*Nxp)
    Kf = pi/Dkp
    Ki=arange(-Kf,Kf,Dkp)
    
    Psi1=fft.fftshift(fft.ifft(array([polC(k) for k in Ki]))*sqrt(len(x)))
    Psi2=fft.fftshift(fft.ifft(array([polX(k) for k in Ki]))*sqrt(len(x)))
    return Psi1[:len(x)],Psi2[:len(x)]
    
def gauss(x,sigma=20):
    from numpy import exp, sqrt, pi                                                     
    return  exp(-0.5*x**2/sigma**2)/sqrt(sqrt(pi*sigma**2))
    
def gauss2D(X,Y,sigma=15.0):
    from numpy import exp, sqrt, pi
    return exp(-0.5*(X**2+Y**2)/sigma**2)/sqrt(pi*sigma**2)

def gaussXT(x,t,C0=1.0,DT1=1.0,Sigl=10,t1=0.5,wl=10):
    from numpy import sqrt,pi,exp
    return (C0/(sqrt(2*pi)*DT1))*gauss(x,sigma=Sigl)*exp(-0.5*((t-t1)/DT1)**2+1j*wl*t)
    
    
def soliton1D(x,t,asol=0.4,vsol=0.0,x0sol=0.0):
    from numpy import exp,cosh
    return asol*exp(1.0j*(vsol*x+(asol**2-vsol**2)*t/2))/cosh(asol*(x-vsol*t-x0sol)) 
        
def soliton2D(X,Y,t):  
    return soliton1D(X,t)*soliton1D(Y,t)

def DCI(x):
    from numpy import linspace, load
    from scipy import interpolate
    yC=load("PsiCD.npy")
    yX=load("PsiXD.npy")    
    x0=linspace(-500,500,2**12)
    fC=interpolate.interp1d(x0,yC)
    fX=interpolate.interp1d(x0,yX)
    return fC(x),fX(x) 
    
def pot(x,A,W):
    from numpy import array
    Dx=x[1]-x[0]
    L=x[-1]
    
    def inte(xi):
        if xi<=(-L+A*Dx):
            v=W*(xi-(-L+A*Dx))**2
        if xi>=(L-A*Dx):
            v=W*(xi-(L-A*Dx))**2
        if xi>(-L+A*Dx) and xi<(L-A*Dx):
            v=0
        return v
    return array([inte(xi) for xi in x])
    
def Airyfun(x,t,hbar=0.658212,m=0.0658191,au=0.05,ki=0.0):
    from numpy import exp
    from scipy.special import airy,ai_zeros
    Xo=10.0
    Bo=(ai_zeros(2)[1][0]-ai_zeros(2)[1][1])/Xo
    B=Bo*hbar**(2.0/3.0)
    xxa=B/hbar**(2.0/3.0)*(x-(B**3.0/4.0*m**2)*t**2)
    xxe=(1j*B**3*t/(2*m*hbar))*(x-(B**3*t**2/(6*m**2)))   
    return airy(xxa)[0]*exp(xxe)*exp(1j*ki*x)*exp(au*x)
    
def CoWavePack(x,t,hbar=1,m=1,a=10,w=0.2,x0=-40):
    from scipy import pi,cos,sin,exp
    K1=m*w**2
    alpha=(m*K1/(hbar**2))**(1.0/4.0)
    zi=alpha*x
    zi0=alpha*a
    return (alpha**(1.0/2.0)/pi**(1.0/4.0))*exp(-(1.0/2.0)(zi-zi0*cos(w*t))**2-1j*((w*t/2)+zi*zi0*sin(w*t)-(1.0/4.0)*zi0**2*sin(2*w*t)))

