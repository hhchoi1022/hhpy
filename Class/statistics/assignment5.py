#%% Problem 1.c)
import numpy as np
from scipy.stats import norm, poisson
observations = np.array([ 3.23,-2.50, 1.88,-0.68, 4.43, 0.17,
                          1.03,-0.07,-0.01, 0.76, 1.76, 3.18,
                          0.33,-0.31, 0.30,-0.61, 1.52, 5.43,
                          1.54, 2.28, 0.42, 2.33,-1.03, 4.00,
                          0.39])
z_05 = norm.ppf(0.95)
z_025 = norm.ppf(0.975)

def tau(mean, std):
    return z_05 * std + mean

mean_hat = np.mean(observations)
se_hat = np.std(observations)
tau_hat = tau(mean_hat, se_hat)

se_delta = se_hat * np.sqrt( (2+z_05**2)/(2*len(observations)))

n_Boot = 10000
taulist = []
for i in range(n_Boot):
    samples = norm.rvs(size=len(observations), loc=tau_hat, scale=se_hat)
    tau_samples = np.std(samples) * z_05 + np.mean(samples)
    taulist.append(tau_samples)
se_boot = np.std(taulist)
print('tau_hat = %.3f, se_delta = %.3f, se_bootstrap = %.3f' %(tau_hat, se_delta, se_boot))
#%% Problem 3.a
import numpy as np
from scipy.stats import norm
n = 100
std = 1
mean = 5
observations = norm.rvs(loc=mean, scale=std, size=n)
#%% Delta method
def theta(mu):
    return np.exp(mu)
z_025 = norm.ppf(0.975)
mu_hat = np.mean(observations)
theta_hat = theta(mu_hat)
se_hat =  theta(np.mean(observations))*np.std(observations)/np.sqrt(n)
CI_delta = (theta_hat - z_025*se_hat, theta_hat + z_025*se_hat)
theta_delta_list = norm.rvs(size = n_Boot, loc = theta_hat, scale = se_hat)
theta_true_list = norm.rvs(size = 1000000, loc = theta_hat, scale = theta(np.mean(observations))*std/np.sqrt(n))
#%% Parametric Bootstrap
n_Boot = 10000
theta_param_list = []
for i in range(n_Boot):
    samples = norm.rvs(size=len(observations), loc=mu_hat, scale=std)
    theta_samples = theta(np.mean(samples))
    theta_param_list.append(theta_samples)
se_param_bootstrap = np.std(theta_param_list)
CI_bootparam= (theta_hat - z_025*se_param_bootstrap, theta_hat + z_025*se_param_bootstrap)
# %% Non-parametric Bootstrap
theta_nonparam_list = []
for i in range(n_Boot):
    samples = np.random.choice(observations, len(observations), replace = True)
    theta_samples = theta(np.mean(samples))
    theta_nonparam_list.append(theta_samples)
se_nonparam_bootstrap = np.std(theta_nonparam_list)
CI_bootnonparam = (theta_hat - z_025*se_nonparam_bootstrap, theta_hat + z_025*se_nonparam_bootstrap)
#%%
print('CI(delta) = (%.3f, %.3f)\nCI(bootstrap_param) = (%.3f, %.3f)\nCI(bootstrap_nonparam) = (%.3f, %.3f) '%(CI_delta[0],CI_delta[1],CI_bootparam[0],CI_bootparam[1],CI_bootnonparam[0],CI_bootnonparam[1]))

#%% Problem 3.b
import matplotlib.pyplot as plt
plt.figure(dpi = 500)
plt.hist(theta_true_list, bins = 300, color= 'r', histtype = 'step', density= True, label = 'True')
plt.hist(theta_delta_list, bins = 50, color= 'k', histtype = 'step', density= True, label = 'delta')
plt.hist(theta_nonparam_list, bins = 50, color= 'y', histtype = 'step', density= True, label = 'Bootstrap(nonparam)')
plt.hist(theta_param_list, bins = 50, color= 'b', histtype = 'step', density= True, label = 'Bootstrap(param)')
plt.xlabel(r'$\mu$')
plt.ylabel('fraction')
plt.legend()
# %% CDF
def CDF(value, observations):
    return np.sum(observations < value)/len(observations)
xrange_ = np.linspace(100, 240, 100)
fraclist_true = []
fraclist_delta = []
fraclist_nonparam = []
fraclist_param = []
for x in xrange_:
    fraclist_delta.append(CDF(x,theta_delta_list))
    fraclist_nonparam.append(CDF(x,theta_nonparam_list))
    fraclist_param.append(CDF(x,theta_param_list))
    fraclist_true.append(CDF(x,theta_true_list))
plt.figure(dpi = 500)
plt.title('CDF')
plt.plot(xrange_, fraclist_true, label = 'true', c ='r')
plt.plot(xrange_, fraclist_delta, label = 'delta', c ='k')
plt.plot(xrange_, fraclist_nonparam, label = 'Non-param bootstrap', c='y')
plt.plot(xrange_, fraclist_param, label = 'param bootstrap', c='b')
plt.xlabel(r'$\mu$')
plt.ylabel('cumulative fraction')
plt.legend()
# %% CDF residual 
comp_delta = np.array(fraclist_true)-np.array(fraclist_delta)
comp_nonparam = np.array(fraclist_true)-np.array(fraclist_nonparam)
comp_param = np.array(fraclist_true)-np.array(fraclist_param)
plt.figure(dpi = 500)
plt.title('CDF residual')
plt.plot(xrange_, comp_delta, label = 'delta', c ='k')
plt.plot(xrange_, comp_nonparam, label = 'Non-param bootstrap', c='y')
plt.plot(xrange_, comp_param, label = 'param bootstrap', c='b')
plt.xlabel(r'$\mu$')
plt.ylabel('cumulative fraction residual')
plt.legend()
# %%


#%% Problem 5.b
from scipy.stats import poisson, norm
alpha = 0.05
lambda_0 = 1
size = 20
n_simulation = 100000
z_025 = norm.ppf(0.975)
z_05 = norm.ppf(0.95)
n_reject = 0
waldlist = []
for i in range(n_simulation):
    sample = poisson.rvs(mu = lambda_0, size = size)
    wald = (np.mean(sample) - lambda_0)/(np.sqrt(lambda_0/size))  
    waldlist.append(wald)  
n_reject = np.sum(np.abs(waldlist) > z_025)
ratio_error1 = n_reject/n_simulation
print('Type I error = %.3f'%(ratio_error1))
#%% Problem 5.d
z_025 = norm.ppf(0.975)
size = 20
lambda_ = 1.5
n_simulation = 1000
def powerfunc(lambda_, n):
    power = 1 - norm.cdf(np.sqrt(n)*(lambda_-1)+z_025) + norm.cdf(np.sqrt(n)*(lambda_-1)-z_025)
    return power
# %%
vallist = []
for i in range(n_simulation):
    sample = poisson.rvs(mu = lambda_0, size = size)
    mean_sample = np.mean(sample)
    val = powerfunc(mean_sample, size)
    vallist.append(val)
n_reject = np.sum(np.array(vallist)>0.5)
print("Rejected counts for null hypothesiss = %d" %(n_reject))
# %%
import pandas as pd
import numpy as np
from scipy.stats import chi2
from astropy.table import Table
import pprint
df = pd.read_csv('RateMyProfessor.csv', sep=',')
features = list(df.columns)[3:]
# %% Problem 6.a)
def Pvalue(x1, x2):
    n1, n2 = len(x1), len(x2)
    Mu1, Mu2 = np.mean(x1), np.mean(x2)
    std1, std2 = np.std(x1, ddof=1), np.std(x2, ddof=1)
    Theta = Mu2 - Mu1
    SE = np.sqrt(std1**2/n1 + std2**2/n2)
    W = Theta/SE
    Pvalue = chi2.sf(W**2, df=1)
    return Pvalue
pvallist = []
for feature in features:
    x1 = df[df[feature] == 0]['student_difficult']
    x2 = df[df[feature] == 1]['student_difficult']
    pvallist.append(Pvalue(x1,x2))
result = Table()
result['feature'] = features
result['P-value'] = pvallist
result.sort('P-value')
# %% Problem 6.b)
result['Bonferroni_P-value'] = result['P-value'] * len(pvallist)
result[result['Bonferroni_P-value'] < 0.05]
# %% Problem 6.c)
"""Reference : https://bioinformaticsandme.tistory.com/129"""

result['Rank'] = np.arange(len(result))+1
pval_FDRlist = []
pval_adjusted = result[-1]['P-value']
for i, feature in enumerate(reversed(result)):
    rank = feature['Rank']
    pvalue = feature['P-value']
    pval_FDR = np.min([pval_adjusted, pvalue*(len(result)/rank)])
    pval_FDRlist.append(pval_FDR)
result['Adjusted_P-value_for_B-H'] = list(reversed(pval_FDRlist))
result[result['Adjusted_P-value_for_B-H'] < 0.05]
# %% Problem 6.d)
n_permutation = 1000000
from tqdm import tqdm
pval_perlist = []
for feature in result['feature']:
    x = df['star_rating']
    x1 = df[df[feature] == 0]['star_rating']
    x2 = df[df[feature] == 1]['star_rating']
    n1, n2 = len(x1), len(x2)
    Mu1, Mu2 = np.mean(x1), np.mean(x2)
    Theta = Mu2 - Mu1
    theta_perlist = []
    for i in tqdm(range(n_permutation)):
        idx_1 = np.random.choice(n1+n2, size=n1, replace=False)
        mu1 = np.mean(x[np.in1d(range(len(x)),idx_1)])
        mu2 = np.mean(x[~np.in1d(range(len(x)),idx_1)])
        theta_per = mu2 - mu1
        theta_perlist.append(theta_per)
    pval = (np.sum(np.array(theta_perlist)**2 >= Theta**2) + 1)/(n_permutation+1)
    pval_perlist.append(pval)
result['Permutation_P-value'] = pval_perlist
result[result['Permutation_P-value']<0.05]
