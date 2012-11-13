from __future__ import division
from random import *
from pylab import *   

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
    
def poisson(lamb, tmax):
    vs = [0]   
    v = 0
    
    while v < tmax:
        v = v + expovariate(lamb)
        vs.append(v)
        
    return vs[:-1]
    
def poissonFun2(lambMax, lamb, tmax):
    vs = poisson(lambMax, tmax)
    
    for v in vs[1:]:
        if uniform(0,1) > lamb(v)/lambMax: vs.remove(v)
        
    return vs
        

def plotPoisson(xs, args=''):
    """Given a list of values generated by a poisson process, 
    adds them to a plot in the canonical fashion
    show() must be called afterwards to display the plot"""
	
    ys = [x for x in range(len(xs)) for _ in (0,1)][:-1]
    
    xs = [xs[0]]+[x for x in xs[1:] for _ in (0,1)]

    plot(xs,ys,args)

def foo(x):
    return 20*sin(x/5)+20
    
def baz(x):
    if x<20:
        return 5
    if x<60:
        return 40
    return 5
    
def main():
    ys = []
    
    samples = 10
    bins = 400
    
    tmax = 200
    lmax = 40
    
    for _ in range(samples):    
        xs = poissonFun2(lmax,baz,tmax)
        #plotPoisson(xs)
        #show()
        ys=xs[1:]+ys
        
    #show()
    ys.sort()
    #plotPoisson(ys)
    #show()
 
    ys,_,_= hist(ys,bins)
    show()
    
    
    ys = [(y*bins)/(samples*tmax) for y in ys]
	
    bar(arange (0,tmax,tmax/bins),ys,tmax/bins)
    plot([baz(x) for x in range(0,tmax)],'red')
    show()
	
	
    
if __name__ == "__main__":
	main()

