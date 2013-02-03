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

def poisson(lamb, t0, tmax):
    if lamb != 0:
        vs = [t0]   
        v = t0

        while v < tmax:
            v = v + expovariate(lamb)
            vs.append(v)

        return vs[:-1]
    else:
        return []

def shadeStates(xs):
    colours = {1:'r',2:'y'  ,3:'g'}
    ys = []

    for x in xs:
        ys.append(x)
        fst = [1]
        xs = dropwhile(lambda (x,y):y==fst,xs)
    
    for n in range(1,len(ys)):
        axvspan(ys[n-1][0],ys[n][0],facecolor=colours[ys[n-1][1]],alpha=0.25,lw=0)

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
    lambd = FloatVector(range(1,states+1))
    
    print delta
    
    Q = [[1 for _ in range(states)] for _ in range(states)]
    
    #for i in range(states):
    #    for j in range(i,states):
    #        Q[i][j] = random()
    
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

    return contents[1:]
    
def malleate(tau):
    f = lambda x: (x-tau[0]).total_seconds()/3600
    
    tau = map(f,tau)
    tau.sort()
    
    return tau
    
    
    

############################################################

#TODO: Play with ~/HiddenMarkov/R/forwardback.mmpp.R

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

print process[:100]
print len(process)
#model = fitmmpp(process,3)
#plotPoisson(process[:100])
#shadeStates([(x,y) for (x,y) in ts if x<process[100]])
#show()
#r['print'](model)
#V = r['Viterbi'](model)

tau = getDates("./twitterextract")
tau = tau[0][1:]
tau = malleate(tau)
print tau

print len(tau)

model = fitmmpp(tau,3)
r['print'](model)
r['print'](r['Viterbi'](model))
V = r['Viterbi'](model)
plotPoisson(tau)
shadeStates(zip(tau,V))
show()
