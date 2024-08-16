######################################################################
#%%####### Problem 1.c)
import scipy.stats
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.stats import bernoulli
N = 1000
p = 0.3
lambda_ = 10

outcome = np.random.poisson(lambda_, N)
X = []
Y = []
for num_toss in outcome:
    numhead = np.sum(bernoulli.rvs(p = p, size = num_toss))
    X.append(numhead)
    Y.append(num_toss - numhead)
cor_coef, p_val = scipy.stats.pearsonr(X,Y)
plt.figure(dpi = 300, figsize = (6,4))
plt.scatter(X,Y, alpha = 0.3, c = 'k', label = f'cor_coef = {round(cor_coef,4)}\n p_value = {round(p_val,4)}')
plt.xlabel('# of heads')
plt.ylabel('# of tails')
plt.legend()



######################################################################
#%%####### Problem 2.b)
N = 1000
n = 20
c = 20000

from scipy.stats import bernoulli
money_fin = []
for i in range(N):
    outcome = 2*bernoulli.rvs(p = 0.5, size = n).astype('float')
    np.place(outcome, (outcome == 0), 0.5)
    money = c
    for j in outcome:
        money *= j
    money_fin.append(money) 
print(f'Expectation = {int(c*(5/4)**n)}, Simulation = {int(np.mean(money_fin))}')



######################################################################
#%%####### Problem 3.b)

#%% Data construction
n = 10000
c = 100
data = 2*bernoulli.rvs(p = 0.5,  size = (c, n))-1
cum_data = np.cumsum(data, axis = 1)
#%% Graph with 90% confidence region
time = np.arange(n)
plt.figure(dpi = 300, figsize = (10,8))
for sole_data in cum_data:
    plt.scatter(time, sole_data, marker = '.', s = 1)
z = norm.ppf(0.95)
plt.plot(time, z * np.sqrt(time), color='k')
plt.plot(time, -z * np.sqrt(time), color='k')
plt.fill_between(time, z * np.sqrt(time), -z * np.sqrt(time), color='gray', alpha=0.3)
plt.xlabel('n')
plt.ylabel('Stock price')
# %% Outlier distribution
counts = []
for i in range(n):
    counts.append(sum(np.abs(cum_data[:,i-1]) < z*np.sqrt(time[i-1])))
plt.figure(dpi = 300, figsize = (10,3))
plt.plot(np.array(counts)/c, c = 'k')
plt.ylim(0.75, 1.1)
plt.xlabel('n')
plt.ylabel('Probability within gray region')
#%% Mean & variance
mean = np.mean(cum_data, axis= 0)
variance = np.sum(cum_data**2 - mean, axis = 0)/c
#%% Mean
plt.figure(dpi = 300, figsize = (10,3))
plt.scatter(time, mean , c = 'k', marker = '.', s = 1)
plt.ylim(-20,20)
plt.plot([0,10000],[0,0],linestyle = '--', linewidth = 1, c= 'k')
plt.ylabel('Mean')
plt.xlabel('n')
# %% Variance
plt.figure(dpi = 300, figsize = (10,3))
plt.scatter(time, variance , c = 'r', marker = '.', s = 1)
plt.plot([0,10000],[0,10000],linestyle = '--', linewidth = 1, c= 'r')
plt.ylim(-500,10500)
plt.ylabel('Variance')
plt.xlabel('n')



######################################################################
#%%####### Problem 5
#%% Problem 5.a)
time_unit = 5
initial_case = 5
lambda_ = 10
p_infect = 0.2
duration = 100
timeslot = np.arange(5, duration+5, time_unit)
N = 100

all_result = []
for n_sim in range(N):
    current_case = initial_case
    num_infect =[]
    for i in range(len(timeslot)):
        n_contact = np.random.poisson(lambda_, current_case)
        new_case = np.sum(np.random.binomial(n_contact, p_infect))
        current_case = new_case
        num_infect.append(new_case)
    all_result.append(num_infect)
#%% Problem 5.b)
time_unit = 5
initial_case = 5
lambda_ = 10
p_infect = 0.2
p_infect_wash = 0.14
duration = 100
timeslot = np.arange(5, duration+5, time_unit)
N = 100

all_result = []
for n_sim in range(N):
    current_case = initial_case
    num_infect =[]
    for i in range(len(timeslot)//2):
        n_contact = np.random.poisson(lambda_, current_case)
        new_case = np.sum(np.random.binomial(n_contact, p_infect))
        current_case = new_case
        num_infect.append(new_case)
    for i in range(len(timeslot)//2):
        n_contact = np.random.poisson(lambda_, current_case)
        new_case = np.sum(np.random.binomial(n_contact, p_infect_wash))
        current_case = new_case
        num_infect.append(new_case)
    all_result.append(num_infect)
    
    
#%% Problem 5.c)
time_unit = 5
initial_case = 5
lambda_ = 10
lambda_distancing = 3
p_infect = 0.2
p_infect_wash = 0.14
duration = 100
timeslot = np.arange(5, duration+5, time_unit)
N = 100

all_result = []
for n_sim in range(N):
    current_case = initial_case
    num_infect =[]
    for i in range(sum(timeslot<=50)): # ~day 50
        n_contact = np.random.poisson(lambda_, current_case)
        new_case = np.sum(np.random.binomial(n_contact, p_infect))
        current_case = new_case
        num_infect.append(new_case)
    for i in range(sum((timeslot>50)&(timeslot <=60))): # day 50 ~ day 60
        n_contact = np.random.poisson(lambda_, current_case)
        new_case = np.sum(np.random.binomial(n_contact, p_infect_wash))
        current_case = new_case
        num_infect.append(new_case)
    for i in range(sum((timeslot >60))): # day 60~
        n_contact = np.random.poisson(lambda_distancing, current_case)
        new_case = np.sum(np.random.binomial(n_contact, p_infect_wash))
        current_case = new_case
        num_infect.append(new_case)
    all_result.append(num_infect)
# %% Plot
from matplotlib import ticker
plt.figure(dpi = 300, figsize = (8,4))
plt.xlim(0,21)
plt.ylim(1e0,1e8)
plt.xlabel('Days')
plt.ylabel('# of new infection')

plt.minorticks_on()
ax = plt.subplot()
ax.boxplot(np.array(all_result), showfliers = False, labels = timeslot)
ax.set_yscale('log')
y_major = ticker.LogLocator(base = 10.0, numticks = 5)
ax.yaxis.set_major_locator(y_major)
y_minor = ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
ax.yaxis.set_minor_locator(y_minor)

