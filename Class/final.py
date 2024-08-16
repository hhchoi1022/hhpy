
#%%
# Conda environment (python 3.10)
from astropy.io import ascii
import numpy as np
from IPython.display import display, Latex
#%% Problem 1
# Load data
data = ascii.read('./halospin.txt')
data.rename_column('col1','idx') # Index
data.rename_column('col2','spin') # Spin parameter
data.rename_column('col3','M_v') # Virial mass

# Problem 1.a)
# Find median value of the virial mass & divide into two subsamples
data.sort('M_v')
idx = len(data)//2
small_M_v = data[:idx]
large_M_v = data[idx:]
# small_M_v has two samples with spin 0 which makes loss function diverges, so exclude from the samples assuming these samples have no measurements of spin parameter
small_M_v = small_M_v[small_M_v['spin'] != 0]
# Problem 1.b)
# Sample mean for each subsampels
lambda_0_small = np.exp(np.sum(np.log(small_M_v['spin']))/len(small_M_v))
lambda_0_large = np.exp(np.sum(np.log(large_M_v['spin']))/len(large_M_v))
sigma_lambda_small = np.sqrt(np.sum((np.log(small_M_v['spin'])-np.log(lambda_0_small))**2)/len(small_M_v))
sigma_lambda_large = np.sqrt(np.sum((np.log(large_M_v['spin'])-np.log(lambda_0_large))**2)/len(large_M_v))
display(Latex('$\lambda_{0,small}$ = %.5f '%(lambda_0_small)))
display(Latex('$\lambda_{0,large}$ = %.5f'%(lambda_0_large)))
display(Latex('$\sigma_{\lambda,small}$ = %.5f '%(sigma_lambda_small)))
display(Latex('$\sigma_{\lambda,large}$ = %.5f '%(sigma_lambda_large)))
#%%
# Problem 1.c)
# Fisher information matrix & Error
error_lambda_0_small = np.sqrt(sigma_lambda_small**2 * lambda_0_small**2 / len(small_M_v))
error_lambda_0_large = np.sqrt(sigma_lambda_large**2 * lambda_0_large**2 / len(large_M_v))
error_sigma_0_small = np.sqrt(sigma_lambda_small**2 /(2*len(small_M_v)))
error_sigma_0_large = np.sqrt(sigma_lambda_large**2 /(2*len(large_M_v)))
display(Latex('$\lambda_{0,small}$ = %.5f $\pm$ %.5f'%(lambda_0_small, error_lambda_0_small)))
display(Latex('$\lambda_{0,large}$ = %.5f $\pm$ %.5f'%(lambda_0_large, error_lambda_0_large)))
display(Latex('$\sigma_{\lambda,small}$ = %.5f $\pm$ %.4f'%(sigma_lambda_small, error_sigma_0_small)))
display(Latex('$\sigma_{\lambda,large}$ = %.5f $\pm$ %.4f'%(sigma_lambda_large, error_sigma_0_large)))
#%%
# Problem 1.d)
import math
sigma_sub = np.sqrt(sigma_lambda_small**2 + sigma_lambda_large**2)
lambda_sub = lambda_0_small-lambda_0_large
n_sub = len(small_M_v) + len(large_M_v)
nu = (n_sub*lambda_sub)/2
Qfunc = lambda x: 0.5*(1-math.erf(x/np.sqrt(2)))
P_FA = float(Qfunc(nu/np.sqrt(n_sub*sigma_sub**2)))
P_D = float(Qfunc((nu-n_sub*lambda_sub)/np.sqrt(n_sub*sigma_sub**2)))
display(Latex('$\sigma_{\lambda,sub}$ = %.5f'%(sigma_sub)))
display(Latex('$\lambda_{0,sub}$ = %.5f'%(lambda_sub)))
display(Latex('$P_{FA}$ = %.5f%%'%(P_FA*100)))
display(Latex('$P_{D}$ = %.5f%%'%(P_D*100)))








#%% Problem 2
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
import numpy as np
from scipy.integrate import quad

sn_data = ascii.read('./Perlmutter1999.data')
abs_peak_magnitude = -18 # setting the universal absolute B-band
H0 = (70 * u.km / u.s /(1e6*u.pc)).to(u.cm/u.s/u.Mpc).value # 70km/s/Mpc to cgs unit
DM = sn_data['Magnitude_sub'] - abs_peak_magnitude 
sn_data['distance'] = 10**((DM + 5)/5)
speec_light = const.c.cgs.value

# Luminosity distance fucntion 
def luminosity_distance(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation):
    component1 = (speec_light * (1+redshift))/H0
    component2 = lambda z: 1/(np.sqrt((1+z)**4*fraction_radiation + (1+z)**3*fraction_mass + (1+z)**2*fraction_curvature + fraction_cosmological))
    int_value = quad(component2, 0, redshift)[0]
    return (component1 * int_value) # Luminosity distance in Mpc

def model_magnitude(distance_Mpc, abs_peak_magnitude):
    return 5*np.log10(distance_Mpc) + 25 + abs_peak_magnitude

def log_likelihood(params):
    abs_peak_magnitude, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation = params
    lum_distance = np.array([luminosity_distance(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation) for redshift in sn_data['z']])
    mag_model = model_magnitude(lum_distance, abs_peak_magnitude)
    log_likelihood_values = -0.5 * np.log(2 * np.pi) - np.log(sn_data['Uncertainty_1']) - 0.5 * ((sn_data['Magnitude'] - mag_model)**2) / (sn_data['Uncertainty_1']**2)
    sum_log_likelihood_values = np.sum(log_likelihood_values)
    print(sum_log_likelihood_values, abs_peak_magnitude, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation)
    return sum_log_likelihood_values
#%%
from scipy.optimize import minimize, LinearConstraint
initial_params = [-19.31, 0.7, 0.3, 0.0, 0]
# Set bounds for parameters
param_bounds = [(-20, -17), (0, 1), (0, 1), (0, 1), (0,1)] 
# Set the cosmological parameters sum 1
equality_constraint = LinearConstraint([[0, 1, 1, 1, 1]], [1], [1]) 
result_powell = minimize(lambda params: -log_likelihood(params), initial_params, bounds=param_bounds, constraints = [equality_constraint], method='SLSQP')
#%% Global minimum
from scipy.optimize import differential_evolution
result_diff = differential_evolution(lambda params: -log_likelihood(params), bounds=param_bounds, constraints = [equality_constraint])

#%% Visualization
import matplotlib.pyplot as plt
redshift_grid = np.arange(0, 1, 0.01)
lum_distance_powell = np.array([luminosity_distance(H0, redshift, result_powell.x[1],result_powell.x[2],result_powell.x[3],result_powell.x[4]) for redshift in redshift_grid])
lum_distance_diff = np.array([luminosity_distance(H0, redshift, result_diff.x[1],result_diff.x[2],result_diff.x[3],result_diff.x[4]) for redshift in redshift_grid])
lum_distance_real = np.array([luminosity_distance(H0, redshift, 0.72, 0.28, 0, 0) for redshift in redshift_grid])
mag_model_powell = model_magnitude(lum_distance_powell, result_powell.x[0])
mag_model_diff = model_magnitude(lum_distance_diff, result_diff.x[0])
mag_model_real = model_magnitude(lum_distance_real, result_diff.x[0])

plt.figure(dpi = 300, figsize = (8,6))
plt.scatter(sn_data['z'], sn_data['Magnitude'], facecolor = None, edgecolor = 'k')
plt.errorbar(sn_data['z'], sn_data['Magnitude'], sn_data['Uncertainty_1'], fmt = 'none', ecolor = 'k')
plt.plot(redshift_grid, mag_model_powell, c = 'r', label = rf'Powell method($\Omega_\Lambda$ = %.2f, $\Omega_m$ = %.2f, $\Omega_k$ = %.2f (fixed), $\Omega_r$ = %.2f (fixed))'%(result_powell.x[1],result_powell.x[2],result_powell.x[3],result_powell.x[4]))
plt.plot(redshift_grid, mag_model_real, c = 'g', label = rf'Permutter1998($\Omega_\Lambda$ = %.2f, $\Omega_m$ = %.2f, $\Omega_k$ = %.2f, $\Omega_r$ = %.2f)'%(0.72, 0.28, 0, 0), linestyle = '--')
plt.legend(loc = 4)
plt.xlabel('z')
plt.ylabel('Peak magnitude')
#%%
# Error
def partial_derivative_omega_lambda(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation):
    component1 = (speec_light * (1+redshift))/H0
    component2 = lambda z: -0.5* 1/(np.sqrt((1+z)**4*fraction_radiation + (1+z)**3*fraction_mass + (1+z)**2*fraction_curvature + fraction_cosmological))**3
    int_value = quad(component2, 0, redshift)[0]
    return (component1 * int_value)# Luminosity distance in Mpc
def partial_derivative_omega_k(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation):
    component1 = (speec_light * (1+redshift))/H0
    component2 = lambda z: -0.5* (1+z)**2/(np.sqrt((1+z)**4*fraction_radiation + (1+z)**3*fraction_mass + (1+z)**2*fraction_curvature + fraction_cosmological))**3
    int_value = quad(component2, 0, redshift)[0]
    return (component1 * int_value) # Luminosity distance in Mpc
def partial_derivative_omega_m(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation):
    component1 = (speec_light * (1+redshift))/H0
    component2 = lambda z: -0.5* (1+z)**3/(np.sqrt((1+z)**4*fraction_radiation + (1+z)**3*fraction_mass + (1+z)**2*fraction_curvature + fraction_cosmological))**3
    int_value = quad(component2, 0, redshift)[0]
    return (component1 * int_value) # Luminosity distance in Mpc
def partial_derivative_omega_r(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation):
    component1 = (speec_light * (1+redshift))/H0
    component2 = lambda z: -0.5* (1+z)**4/(np.sqrt((1+z)**4*fraction_radiation + (1+z)**3*fraction_mass + (1+z)**2*fraction_curvature + fraction_cosmological))**3
    int_value = quad(component2, 0, redshift)[0]
    return (component1 * int_value) # Luminosity distance in Mpc
def partial2_derivative_omega_lambda(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation):
    component1 = (speec_light * (1+redshift))/H0
    component2 = lambda z: (3/4)* 1/(np.sqrt((1+z)**4*fraction_radiation + (1+z)**3*fraction_mass + (1+z)**2*fraction_curvature + fraction_cosmological))**5
    int_value = quad(component2, 0, redshift)[0]
    return (component1 * int_value) # Luminosity distance in Mpc
def partial2_derivative_omega_k(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation):
    component1 = (speec_light * (1+redshift))/H0
    component2 = lambda z: (3/4)* (1+z)**4/(np.sqrt((1+z)**4*fraction_radiation + (1+z)**3*fraction_mass + (1+z)**2*fraction_curvature + fraction_cosmological))**5
    int_value = quad(component2, 0, redshift)[0]
    return (component1 * int_value) # Luminosity distance in Mpc
def partial2_derivative_omega_m(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation):
    component1 = (speec_light * (1+redshift))/H0
    component2 = lambda z: (3/4)* (1+z)**6/(np.sqrt((1+z)**4*fraction_radiation + (1+z)**3*fraction_mass + (1+z)**2*fraction_curvature + fraction_cosmological))**5
    int_value = quad(component2, 0, redshift)[0]
    return (component1 * int_value) # Luminosity distance in Mpc
def partial2_derivative_omega_r(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation):
    component1 = (speec_light * (1+redshift))/H0
    component2 = lambda z: (3/4)* (1+z)**8/(np.sqrt((1+z)**4*fraction_radiation + (1+z)**3*fraction_mass + (1+z)**2*fraction_curvature + fraction_cosmological))**5
    int_value = quad(component2, 0, redshift)[0]
    return (component1 * int_value) # Luminosity distance in Mpc
#%%
def fisher_diagonal_component(H0, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation, absolute_magnitude_peak, partial_derivative, partial_derivative2):
    lum_distance = np.array([luminosity_distance(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation) for redshift in sn_data['z']])
    component1 = 5/sn_data['Uncertainty_1']
    component2 = 1/(lum_distance**2)
    component3 = 1/(np.log(10))**2
    component4 = np.array([partial_derivative(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation)for redshift in sn_data['z']])**2
    model_mag = model_magnitude(lum_distance, absolute_magnitude_peak)
    component5 = (5 + np.log(10)* (sn_data['Magnitude']- model_mag))
    component6 = -5
    component7 = (sn_data['Magnitude']- model_mag)/sn_data['Uncertainty_1']
    component8 = 1/(lum_distance)
    component9 = 1/np.log(10)
    component10 = np.array([partial_derivative2(H0, redshift, fraction_cosmological, fraction_mass, fraction_curvature, fraction_radiation)for redshift in sn_data['z']])
    return sum(component1 * component2 * component3 * component4 * component5 + component6 * component7 * component8 * component9 * component10)
# %%
omega_lambda = result_powell.x[1]
omega_k = result_powell.x[2]
omega_m = result_powell.x[3]
omega_r = result_powell.x[4]
error_lambda = np.sqrt(1/fisher_diagonal_component(H0, result_powell.x[1], result_powell.x[2], result_powell.x[3], result_powell.x[4], result_powell.x[0], partial_derivative_omega_lambda, partial2_derivative_omega_lambda))
error_k = np.sqrt(1/fisher_diagonal_component(H0, result_powell.x[1], result_powell.x[2], result_powell.x[3], result_powell.x[4], result_powell.x[0], partial_derivative_omega_k, partial2_derivative_omega_k))
error_m = np.sqrt(1/fisher_diagonal_component(H0, result_powell.x[1], result_powell.x[2], result_powell.x[3], result_powell.x[4], result_powell.x[0], partial_derivative_omega_m, partial2_derivative_omega_m))
error_r = np.sqrt(1/fisher_diagonal_component(H0, result_powell.x[1], result_powell.x[2], result_powell.x[3], result_powell.x[4], result_powell.x[0], partial_derivative_omega_r, partial2_derivative_omega_r))
# %%
display(Latex('$\Omega_\Lambda$ = %.5f $\pm$ %.5f'%(omega_lambda, error_lambda)))
display(Latex('$\Omega_k$ = %.5f $\pm$ %.5f'%(omega_k, error_k)))
display(Latex('$\Omega_m$ = %.5f $\pm$ %.4f'%(omega_m, error_m)))
display(Latex('$\Omega_r$ = %.5f $\pm$ %.4f'%(omega_r, error_r)))
