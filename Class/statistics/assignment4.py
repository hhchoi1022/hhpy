
#%%
from scipy.stats import norm, cauchy
import matplotlib.pyplot as plt
import numpy as np
#%% Problem 4 #####################################

n = 100
alpha = 0.05
#observations = norm.rvs(size = n) # For standard Normal distribution
observations = cauchy.rvs(size = n) # For standard Cauchy distribution
def calc_epsilon(n, alpha):
    return np.sqrt( (1/(2*n))*np.log((2/alpha)) ) # From the appendix page 99
def CDF(x, observations):
    fraction = np.sum(observations<x)/len(observations)
    return fraction
def CDF_upper(x, observations, epsilon):
    fraction = np.min([CDF(x, observations) + epsilon, 1])
    return fraction
def CDF_lower(x, observations, epsilon):
    fraction = np.max([CDF(x, observations) - epsilon, 0])
    return fraction
# %% Calc confidence bands
epsilon = calc_epsilon(n,alpha) 
xrange_ = np.linspace(np.min(observations), np.max(observations), 1000)
CDF_data = []
for x in xrange_:
    CDF_data.append([CDF(x, observations), CDF_upper(x, observations, epsilon), CDF_lower(x, observations, epsilon)])
CDF_data = np.array(CDF_data)
# %% Single CDF
plt.figure(dpi = 600)
plt.xlabel('value')
plt.ylabel('fraction')
plt.plot(xrange_, CDF_data[:,0], linestyle = '--', c = 'k', label ='CDF[observations]')
#plt.plot(xrange_, norm.cdf(xrange_), linestyle = '--', c = 'r', label ='CDF')
plt.plot(xrange_, cauchy.cdf(xrange_), linestyle = '--', c = 'r', label ='CDF')
plt.fill_between(xrange_, CDF_data[:,1], CDF_data[:,2], alpha = 0.3, color = 'k', label = '95% confidence')
plt.legend()

# %% 1000 iteration CDF
nn = 1000
alpha = 0.05
containlist = []
fraclist = []
plt.figure(dpi =600)
plt.xlabel('value')
plt.ylabel('fraction')
for i in range(nn):
    n = 100
    #observations = norm.rvs(size = n) # For standard normal distribution
    observations = cauchy.rvs(size = n) # For standard Cauchy distribution
    observations.sort()
    epsilon = calc_epsilon(n,alpha) 
    CDF_data = []
    for x in observations:
        #CDF_data.append([CDF(x, observations), CDF_upper(x, observations, epsilon), CDF_lower(x, observations, epsilon), norm.cdf(x)]) # For standard Normal distribution
        CDF_data.append([CDF(x, observations), CDF_upper(x, observations, epsilon), CDF_lower(x, observations, epsilon), cauchy.cdf(x)]) # For standard Cauchy distribution
    CDF_data = np.array(CDF_data)
    contain_key = ((CDF_data[:,3]<=CDF_data[:,1]) & (CDF_data[:,3]>=CDF_data[:,2])).all()
    containlist.append(contain_key)
    fraclist.append(sum(containlist)/len(containlist))
    plt.plot(observations, CDF_data[:,3], c = 'k', alpha = 0.1)
    
    if contain_key == True:
        plt.plot(observations, CDF_data[:,1], c = 'k', alpha = 0.005, linestyle = '--')
        plt.plot(observations, CDF_data[:,2], c= 'k', alpha = 0.005)

    if contain_key == False:
        plt.plot(observations, CDF_data[:,1], c = 'r', alpha = 0.03, linestyle = '--')
        plt.plot(observations, CDF_data[:,2], c = 'r', alpha = 0.03)
    plt.xlim(-10,10)
plt.text(5,0.2,'fraction = %.3f'%fraclist[-1])

#%% Fraction variation by iteration
plt.figure(dpi =600)
plt.xlabel('iteration')
plt.ylabel('fraction in bound')
plt.plot(fraclist, c = 'k')
# %%
#%% Problem 5 #####################################
#%% Problem 5.a)
def theta(mu):
    return np.exp(mu)
n = 10
mu = 5
stdev = 1
tot_sample = norm.rvs(loc = mu, size = n, scale= stdev)
tot_sample.sort()
bins = np.arange(0, np.exp(np.mean(tot_sample[-1:])), 1)
#%%
n_Boot = 100000
samplethetalist = []
z = 1.96
for i in range(n_Boot):
    sample = np.random.choice(tot_sample, size = n, replace = True)
    samplethetalist.append(theta(np.mean(sample)))
std_Boot = np.std(samplethetalist)
mean_Boot = theta(np.mean(tot_sample))
lower_lim = mean_Boot - z*std_Boot
upper_lim = mean_Boot + z*std_Boot
#%%
plt.figure(dpi =600)
plt.xlim(50,500)
plt.hist(samplethetalist, bins = bins, density = True, histtype = 'step', color ='k')
plt.axvline(lower_lim, linestyle = '--', c= 'r', linewidth = 1)
plt.axvline(upper_lim, linestyle = '--', c= 'r', linewidth = 1)
plt.xlabel(r'$\theta$')
plt.ylabel('fraction')
#%% Problem 5.b)
def theta_cdf(x):
    return norm.cdf(np.log(x), loc=5, scale=stdev/np.sqrt(n))
# %%
theta_cdf_bins = theta_cdf(bins)
theta_cdf_bins_delta = np.zeros(len(bins))
theta_cdf_bins_delta[1:] = np.diff(theta_cdf_bins)
# %%
true_samples = np.exp(norm.rvs(loc = 5, scale =stdev/np.sqrt(n), size = 10000000))
plt.figure(dpi =600, figsize = (10,6))
plt.xlim(50,500)
plt.hist(samplethetalist, bins, density = True, color = 'k', histtype = 'step', label = 'Bootstrap')
plt.hist(true_samples, bins, density = True, color = 'r', histtype = 'step', label = 'True sampling distribution')
plt.axvline(np.mean(samplethetalist), linestyle = '--', c= 'k')
plt.axvline(np.mean(true_samples), linestyle = '--', c= 'r')
plt.xlabel(r'$\theta$')
plt.ylabel('fraction')
plt.legend()

#%% Problem 7 #####################################
import pandas as pd
import os
import seaborn as sns
os.chdir('../../Documents/(Class)Math & Statistics/')
df = pd.read_csv('coris.txt', sep=',', skiprows=[0, 1])
#df = df[df['tobacco']>0]
df['logTob'] = np.log(df['tobacco']+0.1)
#%% Problem 7.a)
plt.figure(dpi = 600)
plt.hist(df['tobacco'], bins = 30, density = True, color = 'k', histtype = 'step', label = 'Bootstrap')
plt.xlabel('tobacco')
plt.ylabel('Frequency')
plt.title('Histogram of Tobacco')
# %% Problem 7.b)
plt.figure(dpi = 600)
plt.hist(df['logTob'], bins = 30, density = True, color = 'k', histtype = 'step', label = 'Bootstrap')
plt.xlabel(r'$\log_e({tobacco})$')
plt.ylabel('Frequency')
plt.title('Histogram of Tobacco')

# %% Problem 7.c)
n = len(df)
n_Boot = 1000
sammedlist_0 =[]
sammedlist_1 =[]
df_0 = df[df['chd'] == 0]
df_1 = df[df['chd'] == 1]
for i in range(n_Boot):
    sample_0 = np.random.choice(df_0['logTob'], size=len(df_0), replace=True)
    sample_1 = np.random.choice(df_1['logTob'], size=len(df_1), replace=True)
    median_sample_0 = np.mean(sample_0)
    median_sample_1 = np.mean(sample_1)
    sammedlist_0.append(median_sample_0)
    sammedlist_1.append(median_sample_1)
boot_chd = np.hstack([np.zeros(len(sammedlist_0)),np.ones(len(sammedlist_1))])
boot_logTob = np.hstack([sammedlist_0,sammedlist_1])

#%%
plt.figure(dpi = 600)
sns.violinplot(data=df, x='chd', y='logTob', inner = 'box');
plt.xlabel('CHD')
plt.ylabel(r'$\log_e({tobacco})$')
plt.title('Violin plot of all samples')

#%%
def Get_CI(T, B_Sample):
    # Percentile Interval
    CI_Percent = list(np.quantile(B_Sample, q=[0.025, 0.975]))
    # Normal Interval
    se_b = np.std(B_Sample, ddof=1)
    CI_Normal = [T - 1.96*se_b, T + 1.96*se_b]
    # Pivot Interval
    CI_Pivot = [2*T - CI_Percent[1], 2*T - CI_Percent[0]]

    re = {'CI_Percent(2.5%, 97.5%)':CI_Percent, 'CI_Normal(2.5%, 97.5%)':CI_Normal, 'CI_Pivot(2.5%, 97.5%)':CI_Pivot}
    print(re)
    return CI_Percent,CI_Normal,CI_Pivot, se_b

# %%
CI_perc_0, CI_norm_0, CI_pivot_0, se_0 = Get_CI(np.mean(df_0['logTob']), sammedlist_0)
CI_perc_1, CI_norm_1, CI_pivot_1, se_1 = Get_CI(np.mean(df_1['logTob']), sammedlist_1)
#%%
plt.figure(dpi = 600)
sns.violinplot(x=boot_chd, y=boot_logTob, inner = 'box');
for i in range(2):
    plt.axhline(y = np.mean(df_1['logTob']), xmin = 0.5, xmax = 0.8, linewidth = 0.7, linestyle = '-', c = 'k')
    plt.text(x = 0.4, y = np.mean(df_1['logTob'])-0.015, s = '%.3f'%np.mean(df_1['logTob']), size = 5, c = 'k')
    plt.axhline(y = np.mean(df_0['logTob']), xmin = 0.2, xmax = 0.5, linewidth = 0.7, linestyle = '-', c = 'k')
    plt.text(x = 0.55, y = np.mean(df_0['logTob'])-0.015, s = '%.3f'%np.mean(df_0['logTob']), size = 5, c = 'k')
    #Confidence intervals
    plt.axhline(y = CI_perc_0[i], xmin = 0.2, xmax = 0.3, linewidth = 0.7, linestyle = '--', c = 'r', label = 'Percent(95%)')
    plt.axhline(y = CI_norm_0[i], xmin = 0.2, xmax = 0.3, linewidth = 0.7, linestyle = '--', c = 'b', label  ='Noraml(95%)')
    plt.axhline(y = CI_pivot_0[i], xmin = 0.2, xmax = 0.3, linewidth = 0.7, linestyle = '--', c = 'g', label = 'Pivot(95%)')
    plt.text(x = -0.2, y = CI_perc_0[i]-0.015, s = '%.3f'%CI_perc_0[i], size = 5, c = 'r')
    plt.text(x = -0.2, y = CI_norm_0[i]-0.015, s = '%.3f'%CI_norm_0[i], size = 5, c = 'b')
    plt.text(x = -0.2, y = CI_pivot_0[i]-0.015, s = '%.3f'%CI_pivot_0[i], size = 5, c = 'g')
    plt.axhline(y = CI_perc_1[i], xmin = 0.7, xmax = 0.8, linewidth = 0.7, linestyle = '--', c = 'r')
    plt.axhline(y = CI_norm_1[i], xmin = 0.7, xmax = 0.8, linewidth = 0.7, linestyle = '--', c = 'b')
    plt.axhline(y = CI_pivot_1[i], xmin = 0.7, xmax = 0.8, linewidth = 0.7, linestyle = '--', c = 'g')
    plt.text(x = 0.8, y = CI_perc_1[i]-0.015, s = '%.3f'%CI_perc_1[i], size = 5, c = 'r')
    plt.text(x = 0.8, y = CI_norm_1[i]-0.015, s = '%.3f'%CI_norm_1[i], size = 5, c = 'b')
    plt.text(x = 0.8, y = CI_pivot_1[i]-0.015, s = '%.3f'%CI_pivot_1[i], size = 5, c = 'g')
    
plt.xlabel('CHD')
plt.ylabel(r'$\log_e({tobacco})$')
plt.title('Violin plot of bootstrapped samples')
#%% Problem 7.e)
T_meandiff = np.mean(df_1['logTob'])-np.mean(df_0['logTob'])
n = len(df)
n_Boot = 1000
sammed_diff =[]
df_0 = df[df['chd'] == 0]
df_1 = df[df['chd'] == 1]
for i in range(n_Boot):
    sample_0 = np.random.choice(df_0['logTob'], size=len(df_0), replace=True)
    sample_1 = np.random.choice(df_1['logTob'], size=len(df_1), replace=True)
    median_sample_0 = np.mean(sample_0)
    median_sample_1 = np.mean(sample_1)
    sammed_diff.append(median_sample_1-median_sample_0)
CI_perc, CI_norm, CI_pivot, se = Get_CI(T_meandiff, sammed_diff)
plt.figure(dpi = 600)
sns.histplot(x=sammed_diff, bins=15)
plt.axvline(x = CI_perc[0], linewidth = 0.7, linestyle = '--', c = 'r', label = 'Percent(95%)')
plt.axvline(x = CI_perc[1], linewidth = 0.7, linestyle = '--', c = 'r')
plt.axvline(x = CI_norm[0], linewidth = 0.7, linestyle = '--', c = 'b', label  ='Noraml(95%)')
plt.axvline(x = CI_norm[1], linewidth = 0.7, linestyle = '--', c = 'b')
plt.axvline(x = CI_pivot[0], linewidth = 0.7, linestyle = '--', c = 'g', label = 'Pivot(95%)')
plt.axvline(x = CI_pivot[1], linewidth = 0.7, linestyle = '--', c = 'g')
plt.axvline(x = T_meandiff, linewidth = 1.5, linestyle = '--', c = 'k')
plt.legend()
plt.xlabel(r'$\ln({tob})_{CHD=1} - \ln({tob})_{CHD=0}$')
# %%
