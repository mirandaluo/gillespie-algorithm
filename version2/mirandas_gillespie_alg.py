#molecule x sits there and random stuff happens to it either a new molecule is produced or an existing on degrades. time keeps passing
#general-purpose gillespie simulator
import numpy as np
import matplotlib.pyplot as plt


#parameters
k_prod = 10.0 #rate of producing new molecule (per unit time)
k_deg = 1.0 #rate of degrading molecule
#expect x to fluctuate around 10
#k_prod/k_deg = 10 steady state




#stochiometry matrix rows = species, columns = reactions
#each entry is change in species when reaction happens


stoichiometry = np.array([[+1,-1]]) #+1 when production and -1 when degradation


def propensities(state):
   x= state[0]#extract x from state vector index 0
   a_prod=k_prod
   a_deg = k_deg * x
   return np.array([a_prod,a_deg])






#simulation settings
t_max = 20.0
initial_state = np.array([0])#x starts at 0


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
plt.figure(num = 'single gene medium rates')
plt.step(times, states_array[:,0], where = 'post')
plt.xlabel('time')
plt.ylabel('x(molecule count)')
plt.title('single gene gillespie simulation')
plt.axhline(k_prod/k_deg, color = 'red', linestyle='--', label = 'steady state')
plt.legend()
plt.show()

















