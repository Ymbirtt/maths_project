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
    
def poisson(lamb, tmax):
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
    
    l = len(xs)-1

    ys = [x/500 for x in range(len(xs)) for _ in (0,1)][:-1]
    
    xs = [xs[0]]+[x for x in xs[1:] for _ in (0,1)]

    plot(xs,ys,args, label = "Estimated integral")
    xlabel("Time")
    ylabel("Emissions")
    annotate(ys[-1], (xs[-1],ys[-1]), xytext = (xs[-1]+1,ys[-1]))

def argmin(args,f):
    m=args[0]
    
    for x in args[1:]:
        if f(x)<f(m):
            m=x
            
    return m

def fit(xs,clusters):
    means=[m*sqrt(var(xs)) for m in kmeans(whiten(xs),clusters)][0]
    ys = [0 for _ in range(len(xs))]
    
    for n in range(len(xs)):
        cluster = argmin(means,lambda x:abs(x-xs[n]))
        ys[n] = (xs[n],cluster)
    return (ys,means)

def lookup(x,xs,width):
    index = int(x/width)
    return xs[index][1]

def foo(x):
    return sin(x/5)+1
    
def baz(x):
    if 30<x<50:
        return 10
    else:
        return 5

def int_baz(x):
    if x<30:
        return 5*x
    if x<50:
        return 10*(x-30)+150
    return 5*(x-50)+350
    
def int_foo(x):
    return -5*cos(x/5)+x+5
    
def main():
    ys = []
    samples = 500
    bins = 100
    
    tmax = 100
    lmax = 2
    
    for _ in range(samples):
        xs = poissonFun2(lmax,foo,tmax)
        #plotPoisson(xs)
        #plot([int_foo(x) for x in range(0,tmax)],'red')
        #show()
        ys=xs[1:]+ys
        
    #show()
    ys.sort()
    plotPoisson(ys)
    print len(ys)
    plot([int_foo(x) for x in range(0,tmax)],'red', label = "True integral")
    legend()
    show()
 
    exit()
 
    ys,_,_= hist(ys,bins)
    cla()
    
    
    ys = [(y*bins)/(samples*tmax) for y in ys]

    bar(arange (0,tmax,tmax/bins),ys,tmax/bins)
    plot([baz(x) for x in range(0,tmax)],'red')
    show()

    zs,ms = fit(ys,2)

    f=lambda x:lookup(x,zs,tmax/bins)
    
    plot([baz(x) for x in range(0,tmax)],'red', label = "True rate function")
    plot([f(x) for x in range(0,tmax)],'blue', label = "Estimated rate function")
    legend()
    xlabel("Time")
    ylabel("Emission rate")
    show()
	
    
if __name__ == "__main__":
	main()

