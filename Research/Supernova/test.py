#%%
from companion_interaction_K10 import CompanionInteractionK10
from DOM_interaction_L17 import DOMInteractionL17
import numpy as np
#%%
td = np.arange(0.1, 15, 0.1)
DOM = DOMInteractionL17()
result_DOM, _, _ = DOM.calc_magnitude(td = td, filterset=  'UBVR', visualize= True)
solar_r = 6.955e10
td = np.arange(0.1, 15, 0.1)
COMP = CompanionInteractionK10(rstar = 2e12/2/solar_r, m_wd = 1.4, commonangle= False, kappa = 0.2)
result_comp, _, _ = COMP.calc_magnitude(td, filterset = 'UBVR', visualize= True)
# %%
import matplotlib.pyplot as plt
plt.gca().invert_yaxis()
plt.plot(td, result_DOM['U'])
plt.plot(td, result_comp['U'], linestyle = '--')

# %%


#%%
from astropy.io import ascii
# %%
