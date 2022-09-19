#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


# In[2]:


def drug_model(b, t, params):
    # b es el vector [y,x]
    # params contiene k_a, V, k_e
    k_a, V, k_e = params
    db_dt = np.zeros_like(b) # No es estrictamente necesario. Por claridad se hace
    db_dt[0] = -k_a * b[0] # b[0] = y
    db_dt[1] = (k_a / V) * b[0] - k_e * b[1] #b[1] = x
    
    return db_dt


# In[15]:


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


# In[16]:


init_cond= np.array([1., 0.])
# params contiene k_a, V, k_e
params = (3, 6, 6) #L/h, L, L/h
treatment_time = 2.
dose_interval = 2.
time_steps, results = multidose(drug_model,
                                init_cond,
                                params,
                                treatment_time, 
                                dose_interval)
results = np.concatenate(results, axis = 0)


# In[17]:


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
for i in range(2):
    ax[0].plot(time_steps, results[:,i])

[a.set_xlabel("Treatment time") for a in ax]
ax[0].set_ylabel("Concentration in Stomach")
ax[1].set_ylabel("Concentration in Blood")
fig.tight_layout()


# In[6]:


def simulation(b, t):
    
    """
    b es el vector v,x
    b[0]=x
    b[1]=v
    """
    db_dt = np.zeros_like(b)
    db_dt[0] = b[1]
    db_dt[1] = -k * b[0]/m
        
    return db_dt


# In[7]:


# Vector condiciones iniciales
# donde: x(t=0)=0.2 m y(t=0)=0. m/s

init_conditions= [0.2, 0.]
time_steps = np.arange(0,2, 0.01)
k = 0.5 #N/m
m = 10e-3 #kg


# In[8]:


results = odeint(simulation, init_conditions, time_steps)    


# In[9]:


fig,ax=plt.subplots(2, 1, figsize=[6,8])
#figura 1 
ax[0].plot(time_steps, results[:,0], lw=3, label = "Sim.")
ax[1].plot(time_steps, results[:,1], lw=3)
omega_0 = np.sqrt(k / m) # rad / s
f_0 = omega_0 / (2 * np.pi) # Hz
T_0 = 1 / f_0 # s
[a.axvline(T_0, c = "k", ls = ":", label = "Period") for a in ax]


x = init_conditions[0] * np.cos(omega_0 * time_steps)
v = -init_conditions[0] * omega_0 * np.sin(omega_0 * time_steps)
ax[0].plot(time_steps, x, ls="--", label = "Analytic")
ax[1].plot(time_steps, v, ls="--")
[a.grid() for a in ax]
ax[0].set_ylabel("Position [m]")
ax[1].set_ylabel("Velocity [m/s]")
ax[1].set_xlabel("Time [s]")
ax[0].legend(loc = "best")
ax[0].set_title(f"Period = {T_0:.2f} s")
fig.savefig("graph_t2.png", dpi=300)


# In[10]:


#b = 0.05e-3 #kg/s
b = 0.05
def simulation(y, t):
    
    """
    y es el vector v,x
    y[0]=x
    y[1]=v
    """
    dy_dt = np.zeros_like(y)
    dy_dt[0] = y[1]
    dy_dt[1] = -k * y[0]/m - b * y[1]/m
        
    return dy_dt


# In[11]:


init_conditions= [0.2, 0.]
time_steps = np.arange(0,2, 0.01)
k = 0.5 #N/m
m = 10e-3 #kg


# In[12]:


results = odeint(simulation, init_conditions, time_steps)    


# In[13]:


fig,ax=plt.subplots(2, 1, figsize=[6,8])
#figura 1 
ax[0].plot(time_steps, results[:,0], lw=3, label = "Sim.")
ax[1].plot(time_steps, results[:,1], lw=3)
omega_0 = np.sqrt(k / m) # rad / s
f_0 = omega_0 / (2 * np.pi) # Hz
T_0 = 1 / f_0 # s
[a.axvline(T_0, c = "k", ls = ":", label = "Period") for a in ax]

gamma = b / m
r_plus = 0.5 * (-gamma + np.sqrt(gamma**2 - 4 * omega_0**2 + 0j))
r_minus = 0.5 * (-gamma - np.sqrt(gamma**2 - 4 * omega_0**2 + 0j))

B = init_conditions[0] / (1-r_minus/r_plus)
A = init_conditions[0] - B

print(r_plus)
x = A * np.exp(r_plus * time_steps) + B * np.exp(r_minus * time_steps) # Faltan terminos
v= A * r_plus * np.exp(r_plus * time_steps) + B * r_minus * np.exp(r_minus * time_steps)
#v = init_conditions[0] * omega_0 * - np.cos(omega_0 * time_steps)
ax[0].plot(time_steps, x, ls="--", label = "Analytic")
ax[1].plot(time_steps, v, ls="--")
[a.grid() for a in ax]
ax[0].set_ylabel("Position [m]")
ax[1].set_ylabel("Velocity [m/s]")
ax[1].set_xlabel("Time [s]")
ax[0].legend(loc = "best")
ax[0].set_title(f"Period = {T_0:.2f} s")
fig.savefig("graph_t3_2.png", dpi=300)


# In[43]:


import sympy as sy
x, rp, rm, t, A, B, x0 = sy.symbols("x, r_+ r_- t A B x_0")
gamma, omega, m, b, k = sy.symbols("gamma omega m b k", real=True, positive=True)


# In[53]:


position = A * sy.exp(rp * t) + B * sy.exp(rm * t)
position = position.replace(A, x0-B)
position = position.replace(B, x0 / (1 - rm/rp)).simplify()
position = position.replace(rp, (-gamma + sy.sqrt(gamma**2 - 4 * omega**2)) / 2)
position = position.replace(rm, (-gamma - sy.sqrt(gamma**2 - 4 * omega**2)) / 2).simplify()
position = position.replace(gamma, b/(m))
position = position.replace(omega, k/m).simplify()

display(position)
vel = position.diff(t)
acc = position.diff(t,t)


# In[52]:


sy.solve(sy.exp(-b*t/2/m) - 1e-3, t) #tiempo para que la amplitud disminuya a un 1e-3 veces la original.


# In[ ]:




