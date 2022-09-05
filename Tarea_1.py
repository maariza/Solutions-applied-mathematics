#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


# In[2]:


def derivative_wbacteria(b, t, params):
    # b es el vector [y,x,z]
    # params contiene k_a, V, k_e, rho, alpha
    k_a, V, k_e, rho, alpha = params
    db_dt = np.zeros_like(b) # No es estrictamente necesario. Por claridad se hace
    db_dt[0] = -k_a * b[0] # b[0] = y
    db_dt[1] = (k_a / V) * b[0] - k_e * b[1] #b[1] = x
    db_dt[2] = rho * b[2] * (1 - b[2]) - alpha * b[1] * b[2] #b[2]=z
    
    return db_dt

def derivative_wobacteria(b, t, params):
    # b es el vector [y,x]
    # params contiene k_a, V, k_e
    k_a, V, k_e = params
    db_dt = np.zeros_like(b) # No es estrictamente necesario. Por claridad se hace
    db_dt[0] = -k_a * b[0] # b[0] = y
    db_dt[1] = (k_a / V) * b[0] - k_e * b[1] #b[1] = x
    
    return db_dt
    


# In[3]:


def multidose(derivative, init_cond, params, treatment_time, dose_interval):
    results = [] #truco para que funcione en primera iteración
    n_dose = int(treatment_time // dose_interval)
    time_steps = np.arange(0, dose_interval, 0.01)
    new_dose = init_cond.copy()
    new_dose[1:] = 0 # La døsis se aplica solo al estømago 
    for dose in range(n_dose):
        if dose == 0: # Primera dosis
            init_cond_pass = init_cond
        else: # Segunda dosis en adelante
            init_cond_pass = results[-1][-1] + new_dose
        b = odeint(derivative, init_cond_pass, time_steps, args = (params,))
        results.append(b)
    return np.arange(0, treatment_time, 0.01), np.array(results)   


# In[4]:


init_cond_wbacteria = np.array([1., 0., 1.])
 # params contiene k_a, V, k_e, rho, alpha
params_wbacteria = (0.6, 12, 0.2, 1., 8.)
treatment_time = 25.
dose_interval = 5.
time_steps, results_3_dose = multidose(derivative_wbacteria, 
                                       init_cond_wbacteria, 
                                       params_wbacteria, 
                                       treatment_time, 
                                       dose_interval)
results_3_dose = np.concatenate(results_3_dose, axis = 0)


# In[6]:


fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(5,10))
for i in range(3):
    ax[i].plot(time_steps, results_3_dose[:,i])

[a.set_xlabel("Tiempo tratamiento") for a in ax]
ax[0].set_ylabel("Concentración en estómago")
ax[1].set_ylabel("Concentración en sangre")
ax[2].set_ylabel("Densidad de bacterias en órgano")
fig.tight_layout()
fig.savefig("graph_3.png", dpi=300)


# In[ ]:





# In[6]:


# Punto 2
init_cond_wobacteria = np.array([1., 0.])
params_wobacteria = (0.05, 0.06, 12)
treatment_time = 30
dose_interval = 30.
time_steps, results = multidose(derivative_wobacteria, 
                                       init_cond_wobacteria, 
                                       params_wobacteria, 
                                       treatment_time, 
                                       dose_interval)
results = np.concatenate(results, axis = 0)


# In[7]:


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
for i in range(2):
    ax[i].plot(time_steps, results[:,i])

[a.set_xlabel("Treatment time") for a in ax]
ax[0].set_ylabel("Concentration in Stomach")
ax[1].set_ylabel("Concentration in Blood")
fig.tight_layout()


# In[8]:


# Safety band
treatment_time = 2000. # n_dose -> infinito
time_steps_sb, results_sb = multidose(derivative_wobacteria, 
                                       init_cond_wobacteria, 
                                       params_wobacteria, 
                                       treatment_time, 
                                       dose_interval)
results_sb = np.concatenate(results_sb, axis = 0)


# In[9]:


max_sb = results_sb[-10000:].max(axis=0)
min_sb = results_sb[-10000:].min(axis=0)
print(min_sb, max_sb)

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
for i in range(2):
    ax[i].plot(time_steps_sb, results_sb[:,i])
    ax[i].fill_between(time_steps_sb[[0,-1]], min_sb[i], max_sb[i], color = "r", alpha = 0.5, zorder = 0)

[a.set_xlabel("Treatment time") for a in ax]
ax[0].set_ylabel("Concentration in Stomach")
ax[1].set_ylabel("Concentration in Blood")
fig.tight_layout()


# In[10]:


def single_dose_analytical(t, params, init_cond):
    ka, ke, v = params
    d = init_cond[0]
    return (d * ka) / (v * (ka - ke)) * (np.exp(-ke*t) - np.exp(-ka*t))


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
for i in range(2):
    ax[i].plot(time_steps, results[:,i])
    ax[i].fill_between(time_steps[[0,-1]], min_sb[i], max_sb[i], color = "r", alpha = 0.5, zorder = 0)
[a.set_xlabel("Treatment time") for a in ax]
ax[0].set_ylabel("Concentration in Stomach")
ax[1].set_ylabel("Concentration in Blood")
fig.tight_layout()


# In[ ]:





# In[ ]:




