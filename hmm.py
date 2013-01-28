from __future__ import division
from random  import *
from pylab import *
from rpy2.robjects import *
from itertools import *
from datetime import *

def mindex(xs):
    x = xs[0]
    min = 0
    for i in range(1,len(xs)):
        if xs[i]<x: min = i

    return min

#Todo: make this work for 0-rate processes
def poisson(lamb, t0, tmax):
    vs = [t0]   
    v = t0

    while v < tmax:
        v = v + expovariate(lamb)
        vs.append(v)

    return vs[:-1]

def plotPoisson(xs, args=''):
    """Given a list of values generated by a poisson process, 
    adds them to a plot in the canonical fashion
    show() must be called afterwards to display the plot"""

    ys = [x for x in range(len(xs)) for _ in (0,1)][:-1]
    
    xs = [xs[0]]+[x for x in xs[1:] for _ in (0,1)]

    plot(xs,ys,args)

def fitmmpp(tau, states):
    
    tau   = FloatVector(tau)
    delta = FloatVector([1/states for _ in range(states)])
    lambd = FloatVector([1        for _ in range(states)])
    
    Q = [[random() for _ in range(states)] for _ in range(states)]
    
    for i in range(states):
        Q[i][i] = - (sum(Q[i][:i]+Q[i][i+1:]))
        
    Q = FloatVector([q for q_ in Q for q in q_])
    
    Q = r['matrix'](Q,states)
    
    r('''library(HiddenMarkov)''')

    model = r['mmpp'](tau,Q,delta,lambd)
    
    return r['BaumWelch'](model)
    
def getDates(f):
    
    f = open(f,'r')
    f = f.read().split('\n')
    contents = []
    for l in f:
        l=l.split(',')[1:]
        l = map(lambda x:datetime.strptime(x,"%Y-%m-%d %H:%M:%S"),l)
        l.sort()
        contents.append(l)
    
    
    return contents
    
    
    

############################################################

states = (0,0.5,2)
rates  = {
0  : {0.5: 1/60, 2:1/30},
0.5: {0: 1/10,   2:1/30},
2  : {0: 1/10, 0.5:1/60}
}


state = choice(states)

t_max = 10000
sims = 1

t=0
process = []
ts = [(0,state)]


for _ in range(sims):
    t = 0
    state = choice(states)
    ts = [(0,state)]
    while t<t_max:
        Ts,ss = zip(*[(expovariate(rates[state][s]),s) for s in states if s != state])
        i = mindex(Ts)
        if state != 0: process += poisson(state,t,t+Ts[i])
        t += Ts[i]
        state = ss[i]
        ts.append((t,state))
    print ts

process.sort()
process = [0] + [x for x in dropwhile(lambda x:x==0,process)]

print process
print len(process)
model = fitmmpp(process,3)
print model
