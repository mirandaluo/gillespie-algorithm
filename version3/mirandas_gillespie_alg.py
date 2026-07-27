#single gene molecule x sits there and random stuff happens to it either a new molecule is produced or an existing on degrades. time keeps passing
# -> toggle swithc of 2 genes that muatually repress eachother
#each gene's production is blocked by other's protein (hill function)
#result of bistable system where it picks one of two stable states
#general-purpose gillespie simulator
import numpy as np
import matplotlib.pyplot as plt


#parameters
k_prod = 10.0 #rate of producing new molecule (per unit time)
k_deg = 1.0 #rate of degrading molecule
alpha = 10.0 #max production rate
K = 2.0 #half repression threshold
n = 4.0 #hill coeffiecient and has to be >1
k_deg = 1.0 #degradation rate


#expect x to fluctuate around 10
#k_prod/k_deg = 10 steady state




#stochiometry matrix rows = species, columns = reactions
#each entry is change in species when reaction happens




#hill function adding new bio circuit
#2 species of a and b with 4 reactions: a degraded, b degraded, a produced, b produced
#a and b repress the others production using hill function
#if protein blocks a production, b will also be blocked
#switch of coefficient(larger is sharper)
# when b <k denominate is around 1 and production rate is alpha (full strengith)
#when b is large b>k denominator is huge production rate is around 0(shut off)


#row for a and b
stoichiometry = np.array([[+1,-1,0,0],[0,0,+1,-1]]) #+1 when production and -1 when degradation


def propensities(state):
   #new propensities
   A = state[0]
   B=state[1]


   a_A_prod = alpha/(1+(B/K)**n) #a is made but repressed by b
   a_A_deg = k_deg *A
   a_B_prod = alpha/(1+(A/K)**n) #b is made but repressed by a
   a_B_deg = k_deg * B
   return np.array([a_A_prod,a_A_deg,a_B_prod,a_B_deg])






#simulation settings
t_max = 1000.0
initial_state = np.array([0,0])#x starts at 0


#state
#record for trajectory
times = [0.0]
states = [initial_state.copy()] #list of state vectors
state = initial_state.copy() #numpy is mutable copy dso it doesnt effect both
#current state (change over time)
t= 0.0


#number generator
rng = np.random.default_rng(seed =42) #note: 42 as a placeholder, number chosen doesnt matter seed is indentical
#gillespie main loop
while t < t_max:
   #propensities
   a = propensities(state)
   a_total = a.sum()


   if a_total == 0: #edgecase
       break


   #if events occur at a total rate a_total, the waiting time until the next one is exponentially distrubuted with mean 1/a_total
   time_increment = rng.exponential(1.0/a_total) #poisson processes
   t = t + time_increment #random


   reaction_index = rng.choice(len(a), p=a/a_total)


   state = state + stoichiometry[:, reaction_index]
   #probability of production = a_prod/a_total
   #probability of degradation = a_deg/a_total


   #rng random gives a uniform number in [0,1), if it falls below a_prod /a_total, production wins


   times.append(t) #add new time
   states.append(state.copy())
states_array = np.array(states)
#plot
plt.figure(num = 'two gene gillespie stochastic simulator')
plt.step(times, states_array[:,0], where = 'post', label = 'A')
plt.step(times, states_array[:,1], where = 'post', label = 'B')
plt.xlabel('time')
plt.ylabel('x(molecule count)')
plt.title('two gene gillespie toggle switch simulator gardner et al (2000)')
# no longer used plt.axhline(k_prod/k_deg, color = 'red', linestyle='--', label = 'steady state') a and b will populate from the steps
plt.legend()
plt.show()


#re\sult of bistability















