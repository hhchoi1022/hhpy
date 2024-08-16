
#%%

import astropy.constants as const
import numpy as np
import matplotlib.pyplot as plt
h = const.h.cgs
k = const.k_B.cgs
e = 4.8032e-10 # cm3/2 g1/2 s-1
v_e = 1e8 # cm/s
r_e = 2.1e-8 # cm
c = 3e10 #cm/s
pi = np.pi

h = const.h.cgs
m_e = 9.1094e-28 # g
#%% Problem 1
#1.1
phase_difference = 0.0913 #s
nu_1 = 1610e6 #Hz
nu_2 = 1660e6 #Hz
DM = phase_difference / (4.15e15) / ((1/nu_1**2)-(1/nu_2**2))
#1.2
distance = 6e3
n_e = DM/6e3
# 1.3
lamb_1 = 3e8/nu_1
lamb_2 = 3e8/nu_2
RM_min = (1)/(lamb_1**2 - lamb_2**2)
RM_min_2 = (2*np.pi - 1)/(lamb_1**2 - lamb_2**2)
#1.4
B_ll = 1.23*RM_min/DM
#%% Problem 2
# 2.1
n = 2
h_bar = h/2/np.pi
a0 = h_bar**2 / m_e / e**2
alpha = e**2 / h_bar / c
radius = n**2 * a0
velocity = alpha * c / n
energy = (1/2*alpha**2*m_e*c**2)/n**2
# 2.2
B = -e/c * velocity / radius**2
# 2.3
magnetic_moment = 2 * e / 2 /m_e/c * 1/2 * h_bar
delE = 2 *magnetic_moment*B
delE/1.6e-12
# 2.4
l = 1
j = 1/2
s = -1/2
a = alpha**2 * energy / n / l / (l+1/2) / (l+1)
E_so_1 = a / 2 * (j*(j+1) - l*(l+1) - s*(s+1))

#%% Problem 3 
def gamma_21(omega, g2, temp = 1e4):
    comp1 = 8.63e-6 * omega
    comp2 = g2 * temp**0.5
    return comp1/comp2
# 3.1.1
omega_311 = 0.58
g2_311 = 1
gamma_21_1 = gamma_21(omega = omega_311, g2 = g2_311)
A_311 = 1.6
rho_crit_1 = A_311 / gamma_21_1
# 3.1.2
omega_312 = 0.29
g2_312 = 1
gamma_21_2 = gamma_21(omega= omega_312, g2 = g2_312)
A_312 = 1/4 * (6.1e-4 + 2.3e-1)
rho_crit_2 = A_312 / gamma_21_2
# 3.1.3
omega_313 = 2.29
g2_313 = 5
gamma_21_3 = gamma_21(omega= omega_313, g2 = g2_313)
A_313 = 5/4 * (2.0e-2 + 6.8e-3 + 1.7e-6)
rho_crit_3 = A_313 / gamma_21_3
# 3.2
def ratio(n_e, gam_12, gam_21, a21):
    comp1 = n_e * gam_12
    comp2 = n_e * gam_21 + a21
    return comp1/comp2
n_e_range = np.logspace(2, 8, 1000)
temp = 1e4
g1_1 = 9
g2_1 = 5
E_1 = h*c/5000e-8
gamma_12_1 = g2_1/g1_1*gamma_21_3*np.exp((-E_1/k/temp).value)
result_3_1 = ratio(n_e_range, gamma_12_1, gamma_21_3, A_313)

g1_2 = 5
g2_2 = 1
E_2 = h*c/4363e-8
gamma_12_2 = g2_2/g1_2*gamma_21_1*np.exp((-E_2/k/temp).value)
result_3_2 = ratio(n_e_range, gamma_12_2, gamma_21_1, A_311)

plt.figure(dpi = 300)
plt.plot(n_e_range, result_3_1, c= 'k', label = r'$n(^1 D)/n(^3 P)$')
plt.plot(n_e_range, result_3_2, c= 'r', label = r'$n(^1 S)/n(^1 D)$')
plt.axvline(rho_crit_3, c= 'k', linestyle = '--', label = r'$n_{crit}$')
plt.axvline(rho_crit_1, c= 'r', linestyle = '--', label = r'$n_{crit}$')
plt.legend()
plt.xscale('log')
plt.ylabel('Ratio')
plt.xlabel(r'Electron density($n_e$)')
# 3.3
n_crit_1 = (A_311 + A_312) /(gamma_21_1 + gamma_21_2)
n_crit_2 = (A_313) /(gamma_21_3)

# 3.4
A21 = 0.034
g1 = 5
g2 = 1
E21 = h*c/4363e-8
t_range = np.linspace(6000, 20000, 100)
omega = lambda temp: 0.523 * (temp/1e4)**(-0.390-0.056*np.log((temp/1e4)))
comp1 = g2/g1
comp2 = 8.63e-6 * omega(t_range)
comp3 = t_range**0.5
comp4 = np.exp((-E21/k/t_range).value)
gam_12 = comp1 * comp2 / comp3 * comp4
R = A21/gam_12 / 1e2
plt.figure()
plt.plot(t_range, R, c= 'k')
plt.yscale('log')
plt.ylabel(r'$(j_{4959}+j_{5007})/j_{4363}$')
plt.xlabel('Temperature(K)')
plt.ylim(1e5, 1e8)
#%%
1.6/5.05e-8
rho_crit_1
# %% Problem 4
# 4.2
lamb = 1215.67e-8 # cm
nu = c /lamb
a0 = 5.292e-9 # cm
d21_sq = 2**15 * (e * a0)**2 / 3**10
comp1 = 64 * pi**4 / 3 / lamb**3 / h
A21 = comp1 * d21_sq
# 4.3

comp1 = 8 * pi**2 * e**2
comp2 = m_e * c**3
comp3 = 1/3*nu**2
f12 = A21 / comp1 * comp2 / comp3
f12 = 8*pi**2*m_e*nu/3/h/e**2*d21_sq*3

