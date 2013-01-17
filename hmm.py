from __future__ import division
from random  import *
from pylab import *
from rpy2.robjects import *

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

############################################################

states = (0,0.5,2)
rates  = {
0  : {0.5: 1/60, 2:1/30},
0.5: {0: 1/10,   2:1/30},
2  : {0: 1/10, 0.5:1/60}
}


state = choice(states)

t_max = 10000
t=0
ts = [(0,state)]
process = []


while t<t_max:
    Ts,ss = zip(*[(expovariate(rates[state][s]),s) for s in states if s != state])
    i = mindex(Ts)
    if state != 0: process += poisson(state,t,t+Ts[i])
    t += Ts[i]
    state = ss[i]
    ts.append((t,state))

print ts
print process
#plotPoisson(process)
#show()

print len(process)

#ts = [0] + [process[i+1]-process[i] for i in range(len(process)-1)]

process = FloatVector(process)

r('''library(HiddenMarkov)''')

r('''mmpp(NULL,matrix(c(0,0,0,0,0,0,0,0,0),3),c(1/3,1/3,1/3),c(1,1,1))''')

r_mmpp = r['mmpp']
Q = r('''matrix(c(-1,0.5,0.5,0.6,-1.2,0.6,0.7,0.7,-1.4),3)''')


model = r_mmpp(process,Q,FloatVector([1/3,1/3,1/3]),FloatVector([1,1,1]))

print model

r_BaumWelch = r['BaumWelch']
model = r_BaumWelch(model)
r['print'](model)

#print r['summary'](model)



#for _ in range(100):
#    process = []
#    t=0
#    state = choice(states)
#    ts = [(0,state)]
#    while t<t_max:
#        Ts,ss = zip(*[(expovariate(rates[state][s]),s) for s in states if s != state])
#        i = mindex(Ts)
#        if state != 0: process += poisson(state,t,t+Ts[i])
#        t += Ts[i]
#        state = ss[i]
#        ts.append((t,state))

#    model[0] = FloatVector(process)
#    model = r_BaumWelch(model)  

#print model
