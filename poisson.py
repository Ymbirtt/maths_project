from __future__ import division
import random
   

def poisson(lambMax, lamb, tmax):
    vs = []
    
    v = 0
    u=0
    
    while v < tmax:
        v = v + random.expovariate(1/lambMax)
        u = random.uniform(0,lambMax)
        while u > lamb(v):
            v = v + random.expovariate(1/lambMax)
            u = random.uniform(0,lambMax)
            
        vs.append(v)
        
    return vs[:-1]
    
if __name__ == "__main__":
    print "Hello World!"
