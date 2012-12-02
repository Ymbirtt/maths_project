from __future__ import division
from random import *
from pylab import *   
from scipy.cluster.vq import *


def poissonFun(lambMax, lamb, tmax):
    vs = [0]
    
    v = 0
    u = 0
    
    while v < tmax:
        v = v + expovariate(1/lambMax)
        u = uniform(0,lambMax)
        while u > lamb(v):
            v = v + expovariate(1/lambMax)
            u = uniform(0,lambMax)
            
        vs.append(v)
        
    return vs[:-1]
    
def poisson(lamb,t0 , tmax):
    vs = [0]   
    v = 0
    
    while v < tmax:
        v = v + expovariate(lamb)
        vs.append(v)
        
    return vs[:-1]
    
def poissonFun2(lambMax, lamb, tmax):
    vs = poisson(lambMax, tmax)
    xmax = len(vs)
    x=1
    
    while x<xmax:
        if uniform(0,1)>lamb(vs[x])/lambMax:
            vs.pop(x)
            xmax=xmax-1
        else:
            x=x+1
    return vs
        

def plotPoisson(xs, args=''):
    """Given a list of values generated by a poisson process, 
    adds them to a plot in the canonical fashion
    show() must be called afterwards to display the plot"""

    ys = [x for x in range(len(xs)) for _ in (0,1)][:-1]
    
    xs = [xs[0]]+[x for x in xs[1:] for _ in (0,1)]

    plot(xs,ys,args)

def argmin(args,f):
    m=args[0]
    
    for x in args[1:]:
        if f(x)<f(m):
            m=x
            
    return m

def fit(xs,clusters):
    means=[m*sqrt(var(xs)) for m in kmeans(whiten(xs),2)][0]
    ys = [0 for _ in range(len(xs))]
    
    for n in range(len(xs)):
        cluster = argmin(means,lambda x:abs(x-xs[n]))
        ys[n] = (xs[n],cluster)
        print ys[n]
    return (ys,means)

def lookup(x,xs,width):
    index = int(x/width)
    return xs[index][1]

def foo(x):
    return sin(x/5)+1
    
def baz(x):
    if x<100:
        return 1
    return 0.5
    
def main():
    ys = []
    
    samples = 1000
    bins = 100
    
    tmax = 200
    lmax = 2
    
    for _ in range(samples):    
        xs = poissonFun2(lmax,foo,tmax)
        #plotPoisson(xs)
        #show()
        ys=xs[1:]+ys
        
    #show()
    ys.sort()
    #plotPoisson(ys)
    #show()
 
    ys,_,_= hist(ys,bins)
    cla()
    
    
    ys = [(y*bins)/(samples*tmax) for y in ys]

    bar(arange (0,tmax,tmax/bins),ys,tmax/bins)
    plot([foo(x) for x in range(0,tmax)],'red')
    show()

    zs,ms = fit(ys,2)
       
    f=lambda x:lookup(x,zs,tmax/bins)
    
    plot([foo(x) for x in range(0,tmax)],'red')
    plot([f(x) for x in range(0,tmax)],'blue')
    show()
	
    
if __name__ == "__main__":
	main()

