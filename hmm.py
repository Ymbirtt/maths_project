from __future__ import division
from random  import *
from pylab import *
from rpy2.robjects import *
from itertools import *
from datetime import *
from scipy.stats import *
from numpy import *
import pickle

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
    colours = {1:'r',2:'y'  ,3:'g', 4:'c', 5:'k', 0:'b', 0.5:'r', 0.01:'g'}
    ys = xs
    
    for n in range(1,len(ys)):
        print "shading ", n 
        axvspan(ys[n-1][0],ys[n][0],facecolor=colours[ys[n-1][1]],alpha=0.25,lw=0)

def plotPoisson(xs, args=''):
    """Given a list of values generated by a poisson process,
    adds them to a plot in the canonical fashion
    show() must be called afterwards to display the plot"""

    ys = [x for x in range(len(xs)) for _ in (0,1)][:-1]

    xs = [xs[0]]+[x for x in xs[1:] for _ in (0,1)]

    plot(xs,ys,args)
    xlabel("Time")
    ylabel("Emissions")

def plotPoisson2D(xs):
    """Given a list of values generated by a poisson process, adds them to a
    2d plot, with days up the y-axis and hours across the x-axis"""

    days = [x.day + x.month*31         for x in xs]
    xs   = [(x.second + x.minute*60 + x.hour*3600)/3600  for x in xs]

    scatter(xs,days,marker = 'x')
    xlabel("Hour")
    ylabel("Day")
    xlim(0,24)
    ylim(min(days),max(days))

def shadeStates2D(xs,states):
    """Given a list of times and states at each time, shades in the states at
    each time (badly)"""
    colours = {1:'r', 2:'y', 3:'g', 4:'c', 5:'k', 6:'m', 7:'b'}

    print len(states)
    print len(xs)

    days = [x.day + 31*x.month for x in xs]
    xs   = [(x.second + x.minute*60 + x.hour*3600)/3600  for x in xs]

    for n in range(1,len(days)-1):
        if days[n-1]==days[n]:
            #Shade in the area between these two states
            gca().add_patch(Rectangle((xs[n-1],days[n]-0.5),xs[n]-xs[n-1],1,alpha=0.5,color=colours[states[n-1]],linewidth=0))
        elif days[n]-days[n-1]==1:
            #shade in the area from this state to the end of the day, then from the start of the next to the next day
            gca().add_patch(Rectangle((xs[n-1],days[n-1]-0.5),24-xs[n-1],1,alpha=0.5,color=colours[states[n-1]],linewidth=0))
            gca().add_patch(Rectangle((0,days[n]-0.5),xs[n],1,alpha=0.5,color=colours[states[n-1]],linewidth=0))
        else:
            #shade in the area from this state to the end of the day, then from the start of the next to the next day
            gca().add_patch(Rectangle((xs[n-1],days[n-1]-0.5),24-xs[n-1],1,alpha=0.5,color=colours[states[n-1]],linewidth=0))
            gca().add_patch(Rectangle((0,days[n]-0.5),xs[n],1,alpha=0.5,color=colours[states[n-1]],linewidth=0))
            gca().add_patch(Rectangle((0,days[n-1]+0.5),24,days[n]-days[n-1]-2,alpha=0.5,color=colours[states[n-1]],linewidth=0))


def fitmmpp(tau, states = -1, maxiters = 500):

    if states == -1:
        xs = []
        mini = (-1,Inf,None)
        print "Estimating states..."

        for k in range(2,10):
            m = fitmmpp(tau,states=k, maxiters = maxiters)
            LL = m[-3][0]
            BIC = ((k*k+k-1)*log(len(tau))) - (2*LL)
            if BIC<mini[1]:
                mini = (k,BIC,m)
            print "k= " + str(k)
            print "Best fit for k=" + str(mini[0]) + ", BIC= " + str(mini[1])
            raw_input()

    else:

        r('''library(HiddenMarkov)''')

        control = r['bwcontrol'](maxiters)
        tau     = FloatVector(tau)
        delta   = FloatVector([1/states for _ in range(states)])
        lambd   = FloatVector(range(1,states+1))

        Q = [[1 for _ in range(states)] for _ in range(states)]

        for i in range(states):
            Q[i][i] = - (sum(Q[i][:i]+Q[i][i+1:]))

        Q = FloatVector([q for q_ in Q for q in q_])

        Q = r['matrix'](Q,states)

        model = r['mmpp'](tau,Q,delta,lambd)

        return r['BaumWelch'](model, control)

def fitdthmm(tau, states = -1, maxiters = 500):

    if states == -1:
        xs = []
        mini = (-1,Inf,None)
        print "Estimating states..."

        for k in range(2,10):
            m = fitdthmm(tau,states=k, maxiters = maxiters)
            LL = m[-3][0]
            print "LL=" + str(LL)
            BIC = ((k*k+k-1)*log(len(tau))) - (2*LL)
            if BIC<mini[1]:
                mini = (k,BIC,m)
            print "k= " + str(k)
            print "Best fit for k=" + str(mini[0]) + ", BIC= " + str(mini[1])
            raw_input()
    else:

        r('''library(HiddenMarkov)''')

        t = [(tau[n]-tau[n-1]) for n in range(1,len(tau))]
        delta = [1/states for _ in range(states)]
        distn = "exp"
        pm = [s+1 for s in range(states)]
        ones = [1 for _ in range(states)]
        Pi = [[1/states for _ in range(states)] for _ in range(states)]
        control = r['bwcontrol'](maxiters)

        Pi = FloatVector([q for q_ in Pi for q in q_])
        Pi = r['matrix'](Pi,states)
        print t

        t = FloatVector(t)
        delta = FloatVector(delta)
        pm = FloatVector(pm)
        ones = FloatVector(ones)
        model = r['dthmm'](t,Pi,delta,distn,r['list'](rate=pm))
        r['print'](model)
        raw_input()
        return r['BaumWelch'](model, control)


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

    #The first few results seem to be outliers - let's lop some of them off
    return tau[50:]

def testGoodness(ss,model):
    lambdas = list(model[4][0])
    N = len(lambdas)
    ts = model[0]
    ss = map(lambda x:int(x)-1,ss)

    print ss
    #raw_input()
    print N
    #raw_input()
    print lambdas
    #raw_input()
    print len(ts)
    #raw_input()
    print len(ss)
    #raw_input()

    ts_for_states = [[] for _ in lambdas]

    #print ss

    #for x in range(1,len(ts)):
    #      if ss[x] == ss[x-1]:
    #          ts_for_states[ss[x]].append(ts[x])

    for (t,s) in zip(ts,ss):
          ts_for_states[s].append(t)

    for x in range(N):
        if ts_for_states[x]:
            print "Testing fit for state " + str(x)
            print str(len(ts_for_states[x])) + " events"
            mu = mean(ts_for_states[x])
            print "mu = ", mu
            print "1/mu = ", 1/mu
            print "var = " + str(std(ts_for_states[x]))
            print "lambda = " + str(lambdas[x])
            print "1/lambda = " + str(1/lambdas[x])
            #print(ts_for_states[x])
            
            #xs = FloatVector([log(t) for t in ts_for_states[x] if t!=0])
            xs = FloatVector(ts_for_states[x])
            d = r['density'](xs)
            
            #r['svg']('./write-up/images/density_mmpp1_state' + str(x+1) + '.svg')
            #r['plot'](d, main="Estimated probability density of state " + str(x+1), col = "blue", xlab = "x")
            #r('curve(dexp(x, rate ='+ str(1/mu) + '), add=TRUE, col = "red")')
            #r['dev.off']()
            
            r['svg']('./write-up/images/density_mmpp2_state' + str(x+1) + '.svg')
            r['plot'](d, main="Estimated probability density of state " + str(x+1), col = "blue", xlab = "log x")
            r('curve(dnorm(x, mean ='+ str(mean(xs)) + ', sd =' + str(std(xs)) +'), add=TRUE, col = "red")')
            r['dev.off']()
            #r['print'](r['ks.test'](xs,'pexp',lambdas[x]))
            #r['print'](r['ks.test'](FloatVector(ts_for_states[x]),'pnorm',lambdas[x]))
        else:
            print "State " + str(x) + " is empty"
        raw_input()

############################################################
print "Doing nothing of any importance..."

states = (0.01,0.5,2)
rates  = {
0.01  : {0.5: 1/60, 2:1/30},
0.5: {0.01: 1/10,   2:1/30},
2  : {0.01: 1/10, 0.5:1/60}
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

#model = fitmmpp(process,3)
#plotPoisson(process)
#shadeStates(ts)
#xlim(0,10000)
#savefig("./write-up/images/trace_mmpp1.svg")
#show()
#V = r['Viterbi'](model)
#r['print'](model)
#raw_input()
#plotPoisson(process)
#print "shading states"
#print V
#xlim(0,10000)
#shadeStates(zip(process,list(V)))
#print "done!"
#savefig("./write-up/images/fit_mmpp1.svg")
#show()
print "Reading data..."

tau = getDates("./twitterextract")

print "Preprocessing data..."

model = None
tau = tau[0][1:]
times = tau
tau = malleate(tau)

print "Fitting mmpp..."

#model = fitdthmm(tau,5, maxiters = 1500)

#r['print'](r['summary'](model))

#raw_input()

print "Pickling model..."

model = pickle.load(open('./mmpp.stash','r'))
#pickle.dump(model,open('./mmpp.stash','w'))

print "Predicting state transitions..."

r('''library(HiddenMarkov)''')
V = list(r['Viterbi'](model))

testGoodness(V,model)

print "Plotting process..."

#plotPoisson(tau)
plotPoisson2D(times[50:])

#print "Shading states..."

#shadeStates(zip(tau,V))
shadeStates2D(times[50:],V)

print "OH MY GOD IT'S A GRAPH"

show()


