import  numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# Definition of parameters

beta1  = 0.060    #  Effective contact rate between individuals
beta2  = 0.22     #  Per capita contact rate for humans and contaminated water
alpha  =  100     #  Shedding rate
k      = 100000   #  Half-saturation constant
gamma0 = 0.05     #  Minimum recovery rate
gamma1 = 0.09     #  Maximum recovery rate
b      = 2        #  Hospital-bed population ratio
delta  = 30       #  Bacterial net death rate
mu     = 0.06     #  Natural death rate of humans

#define reproduction number
# R_{0} < 1

R0 = beta1/(gamma1 + mu) + (alpha*beta2)/(k*delta*(gamma1 + mu))
print R0

def dYdt(Y,t=0):
    """ Return the human population and bacteria concentration in environment """
    return np.array([mu - beta1*Y[0]*Y[1] - beta2*(Y[0]*Y[3])/(Y[3]+k) - mu*Y[0],
                  beta1*Y[0]*Y[1] + beta2*(Y[0]*Y[3])/(Y[3]+k) - \
                   (gamma0 + ((gamma1 - gamma0) * b )/ (Y[1] + b)) * Y[1] - mu*Y[1],
                  (gamma0 + ((gamma1 - gamma0)* b) / (Y[1] + b))* Y[1] - mu*Y[2],
	           alpha*Y[1] - delta*Y[3] ])


t = np.linspace(0, 150, 1000)              # time
Y0 = np.array([0.80, 0.15, 0.04, 0.40])        # initials conditions
Y, infodict = integrate.odeint(dYdt, Y0, t, full_output=True)
infodict['message']  
s, i, r, x = Y.T


plt.figure(1)
plt.plot(t, s, 'g-',lw = 1, label='Susceptible$\quad$individuals')
plt.xlabel("Time$(t)$$\quad$in$\quad$days",fontsize = 16)
plt.ylabel("$s\quad$$($Population$)$",fontsize = 16)
plt.legend(loc = '$upper right$')
plt.grid(True)
plt.savefig('Susc_DFE.png')

plt.figure(2)
plt.plot(t, r, 'g-',lw = 1, label='Recovered$\quad$individuals')
plt.xlabel("Time$(t)$$\quad$in$\quad$days",fontsize = 16)
plt.ylabel("$r\quad$$($Population$)$",fontsize = 16)
plt.legend(loc = '$upper right$')
plt.grid(True)
plt.savefig('Recov_DFE.png')

plt.figure(3)
plt.plot(t, i, 'g-',lw = 1, label= 'Infected$\quad$individuals')
plt.xlabel("Time$(t)$$\quad$in$\quad$days",fontsize = 16)
plt.ylabel("$i\quad$$($Population$)$",fontsize = 16)
plt.legend( loc ='$upper right$' , fontsize=10)
plt.grid(True)
plt.savefig('Infect_DFE.png')

#plt.figure(4)
#plt.xlabel("Time$(t)$$\quad$in$\quad$days",fontsize = 16)
#plt.ylabel("$x\quad$$($Population$)$",fontsize = 16)
#plt.legend(loc = '$upper right$')
#plt.grid(True)
#plt.savefig('VCholerae_DFE.png')
plt.show()
