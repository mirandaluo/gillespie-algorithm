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


#ensemble experiment on how does the population split


#how many independent runs to do each run is one virtual cell
n_runs = 50


#2 empty list to collect final molecuef counts from each run and use append inside the loop to add one number per run
final_A = []
final_B = []


#loop over seeds and range of n_runs generate the numbers 0 through n_runs-1
#each iteration runs one simulation with a diff random seed whcih means each run gets a diff stream of nunbers
#so the trajectory is diff even through parameters are the same
for seed in range(n_runs):
   #call simulate with the same parameters but w diff seed each time and the function returns 2 arrays of event times and the state of each event
   #Wmight need to change names because they shadow the ones from the single run but its ok cus we dont need old values anaymore
   times,states= simulate(alpha = 10.0, K=2.0,n=4.0, k_deg = 1.0, t_max = 1000.0, seed = seed)
   #states is a 2d array of shape [num_events,2]
   #last row and final count of a is [-1,0] and b is [-1,1] and append
   final_A.append(states[-1,0]) #final A
   final_B.append(states[-1,1]) #final B count from this run
#after the loop both lists have 50 entries and convert them to numpy arrays so we can do math easily
final_A = np.array(final_A)
final_B = np.array(final_B)


#make a new figure window and separate from the single run plot
plt.figure(num = 'ensemble final states')
plt.scatter(final_A,final_B, alpha = 0.6)
plt.xlabel('final count of A')
plt.ylabel('final count of B')
#
plt.title(f'final states from {n_runs} independent runs')
plt.axhline(0, color = 'gray', lw = 0.5)
plt.axvline(0,color = 'gray', lw=0.5)
plt.show()


#quick stats


#np.sum() shows how many true values tehre are and final_A > final_B is boolean array
#how many runs ended with a higher than b
n_A_wins = np.sum(final_A > final_B)
n_B_wins = np.sum(final_B > final_A)




#print counts and percentages
#new sytax is .0f formats with 0 decimal places
print(f"A wins: {n_A_wins} runs ({100*n_A_wins/n_runs:.0f}%)")
print(f"B wins: {n_B_wins} runs ({100*n_B_wins/n_runs:.0f}%)")




#explore how bistability works with alpha
alpha_values = [5.0,10.0,20.0,50.0]
n_runs = 50


#make figure with 2x2 subplots


fig, axes = plt.subplots(2,2,figsize = (10,10), num = 'exploring alpha!')
#axes is a 2d array of subplot objects and flatten it for easier looping
axes = axes.flatten()
for i, alpha_val in enumerate(alpha_values):
   final_A = []
   final_B = []


   for seed in range(n_runs):
       times,states = simulate(alpha=alpha_val, K = 2.0, n = 4.0, k_deg = 1.0,t_max = 1000.0, seed=seed)
       final_A.append(states[-1,0])
       final_B.append(states[-1,1])


   final_A = np.array(final_A)
   final_B = np.array(final_B)


   #plot into subplot i
   ax = axes[i]
   ax.scatter(final_A, final_B, alpha=0.6)
   ax.set_xlabel('final A')
   ax.set_ylabel('final B')
   ax.set_title(f'a = {alpha_val}')
   ax.axhline(0, color = 'gray', lw = 0.5)
   ax.axvline(0, color = 'gray', lw = 0.5)
plt.tight_layout()
plt.show()


# summary stats function which takes trajectory (time+states) and returns a fixed length feature vector
def compute_summary(times,states):
   """args: times: 1d numpy array, shape(T,). the time of each event
    states: 2d numpy array, shape (T,2).
     column 0 = count of A, column 1 = count of B
      returns:
       features:1d numpy array of lenght 9 """
   #hold long each state held np.diff(times) computes [times[1]-times[0], times[2]-times[1],..]
   #each entry is the duration the state stayed constant between events and if times has length T, durations has length T-1
   durations = np.diff(times)
   #total sim time (should be close to t_max)
   total_time = durations.sum()
   #pair each duration with the state held during it
   #the state at row i was held from times[i] to times[i+1] which is a duration of durations[i]
   #so use states [:-1] (everything except last row) to align with durations
   #0 = column 0 (species A) and 1 = column 1 (species B)
   A = states[:-1,0]
   B = states[:-1,1] #1d array of b values t-1


   #time weighted means
   # time weighted mean: sum(value * duration it was held)/total time
   #accounts for some states held longer than others
   #(A*durations) is element-wise multiplication for each event, multiply A value by how long it lasted and sum and divide
  
   mean_A = (A*durations).sum()/total_time
   mean_B = (B*durations).sum()/total_time


   #time weighted variances
   #vaariance = avg squared deviation from the mean
   #same time-weighting idea: each deviation is weighted by how long that state was held
   #(A-mean_A) is element wise subtraction(each entry of A minus the scalar mean_A)
   #**2 squares each element then multiply by durations, sum, divide


   var_A=((A-mean_A)**2*durations).sum()/total_time
   var_B=((B-mean_B)**2*durations).sum()/total_time


   #faction of time A>B
   #boolean array where if a exeeded b, its true but false other wise
   #.astype(float) converts true to 1.0 and false to 0.0
   #time weight how much time was spent with A>B
   A_winning = (A>B).astype(float)
   frac_A_high = (A_winning *durations).sum()/total_time
   #count switches
   #swithc happens when A_winning changes between conscutive events
   switch_events = np.diff(A_winning)    #np.diff(A_winning) give +1 when A started winning, -1 where A stops winning, 0 otherwise
   n_switches = np.abs(switch_events).sum()  #np.abs makes both +1 and -1 become 1 so summing counts both kinds of transitions


   #time weighted mean of total protein
   #just a+b time averaged the same way


   total = A + B
   mean_total = (total*durations).sum()/total_time


   #maximum values
   #no time weight only the largest value ever observed for species
   #states [:,0] is all rows, column 0 (full A trajectory)
   max_A = states[:,0].max()
   max_B = states[:,1].max()


   #pack it all in one fixed lenght numpy array


#feauture vector which is 9 numbers
   features = np.array([
       mean_A,
       mean_B,
       var_A,
       var_B,
       frac_A_high,
       n_switches,
       mean_total,
       max_A,
       max_B,
   ])
   return features
#test to compute summary stats for one trajectory and print


times, states = simulate(alpha = 10.0, K=2.0, n = 4.0, k_deg = 1.0, t_max = 1000.0, seed = 42)
#use the parameters that show spontaneous switching(3 images before)
features = compute_summary(times,states)
print("summary stats for one toggle switch run (seed = 42)")
print(f"  mean_A = {features[0]:.2f}")
print(f"  mean_B = {features[1]:.2f}")
print(f"  var_A = {features[2]:.2f}")
print(f"  var_B = {features[3]:.2f}")
print(f"  frac_A>B = {features[4]:.2f}")
print(f"  n_switches = {features[5]:.0f}")
print(f"  mean_total = {features[6]:.2f}")
print(f"  max_A = {features[7]:.0f}")
print(f"  max_B  = {features[8]:.0f}")















