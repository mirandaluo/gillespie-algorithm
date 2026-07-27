#single gene molecule x sits there and random stuff happens to it either a new molecule is produced or an existing on degrades. time keeps passing
# -> toggle swithc of 2 genes that muatually repress eachother
#each gene's production is blocked by other's protein (hill function)
#result of bistable system where it picks one of two stable states
#general-purpose gillespie simulator wrapped as a function
import numpy as np
import matplotlib.pyplot as plt


#parameters
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
#stoichiometry doesnt depend on parameters so it can live outside
stoichiometry = np.array([[+1,-1,0,0],[0,0,+1,-1]]) #+1 when production and -1 when degradation
def simulate(alpha, K, n, k_deg, t_max,seed, initial_state = None):
   if initial_state is None: # if user doesn't provide initial state start with 0 molecules of A and B
       initial_state = np.array([0,0])
   def propensities(state): #probability per unit time
       A, B = state
       prod_A = alpha/(1+(B/K)**n) #productoin rate of A and b represses using Hill function
       deg_A = k_deg * A #more a molecules more chance to degrade
       prod_B = alpha/(1+(A/K)**n)
       deg_B = k_deg * B


       #4 reaction rates
       return np.array([prod_A,deg_A,prod_B,deg_B])








   #state
   #record for trajectory
   t= 0.0
   times = [t]
   states = [initial_state.copy()] #list of state vectors
   state = initial_state.copy() #numpy is mutable copy dso it doesnt effect both
   #current state (change over time)




   #number generator
   rng = np.random.default_rng(seed) #note: 42 as a placeholder, number chosen doesnt matter seed is indentical
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
       if t >t_max:
           break


       reaction_index = rng.choice(len(a), p=a/a_total)


       state = state + stoichiometry[:, reaction_index]
       #probability of production = a_prod/a_total
       #probability of degradation = a_deg/a_total


       #rng random gives a uniform number in [0,1), if it falls below a_prod /a_total, production wins


       times.append(t) #add new time
       states.append(state.copy())
   return np.array(times), np.array(states)


#run simulation
times, states_array = simulate(alpha = alpha, K=K, n=n, k_deg = k_deg, t_max = 1000.0, seed = 42, initial_state = np.array ([0,0]))


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















