import  numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# Definition of parameters

beta1  = 0.08     # Effective contact rate between individuals
beta2  = 0.04     # Per capita contact rate for humans and contaminated water
alpha  = 10       # Shedding rate
k      = 1000000   # Half-saturation constant
gamma0 = 0.001      # Minimum recovery rate
gamma1 = 0.04      # Maximum recovery rate
delta  = 0.8         # Bacterial net death rate
mu     = 0.00524  # Natural death rate of humans

# different values for the hospital-bed population ratio
b1 = 0.302
b2 = 0.521
b3 = 1.990


#define reproduction number
# R_{0} > 1
R0 = beta1/(gamma1 + mu) + (alpha*beta2)/(k*delta*(gamma1 + mu))
print R0


def dYdt1(Y1,t):
	""" Return the human population and bacteria concentration in environment """
	return np.array([mu - (beta1*Y1[0]*Y1[1]) - ((beta2*(Y1[0]*Y1[3]))/(Y1[3]+k)) - (mu*Y1[0]),
				  beta1*Y1[0]*Y1[1] + beta2*(Y1[0]*Y1[3])/(Y1[3]+k) - 
				   (gamma0 + ((gamma1 - gamma0) * b1 )/ (Y1[1] + b1)) * Y1[1] - mu*Y1[1],
				  (gamma0 + ((gamma1 - gamma0)* b1) / (Y1[1] + b1))* Y1[1] - mu*Y1[2],
			   alpha*Y1[1] - delta*Y1[3] ])


t = np.linspace(0, 800, 1000)              # time
Y01 = np.array([0.8, 0.15, 0.05, 0.40])        # initials conditions
Y1, infodict = integrate.odeint(dYdt1, Y01, t, full_output=True)
infodict['message']  
s1, i1, r1, x1 = Y1.T
def dYdt2(Y2,t):
	""" Return the human population and bacteria concentration in environment """
	return np.array([mu - (beta1*Y2[0]*Y2[1]) - ((beta2*(Y2[0]*Y2[3]))/(Y2[3]+k)) - (mu*Y2[0]),
				  beta1*Y2[0]*Y2[1] + beta2*(Y2[0]*Y2[3])/(Y2[3]+k) - 
				   (gamma0 + ((gamma1 - gamma0) * b2 )/ (Y2[1] + b2)) * Y2[1] - mu*Y2[1],
				  (gamma0 + ((gamma1 - gamma0)* b2) / (Y2[1] + b2))* Y2[1] - mu*Y2[2],
			   alpha*Y2[1] - delta*Y2[3] ])



Y02 = np.array([0.8, 0.15, 0.05, 0.40])        # initials conditions
Y2, infodict = integrate.odeint(dYdt2, Y02, t, full_output=True)
infodict['message']  
s2, i2, r2, x2 = Y2.T

def dYdt3(Y3,t):
	""" Return the human population and bacteria concentration in environment """
	return np.array([mu - (beta1*Y3[0]*Y3[1]) - ((beta2*(Y3[0]*Y3[3]))/(Y3[3]+k)) - (mu*Y3[0]),
				  beta1*Y3[0]*Y3[1] + beta2*(Y3[0]*Y3[3])/(Y3[3]+k) - 
				   (gamma0 + ((gamma1 - gamma0) * b3 )/ (Y3[1] + b3)) * Y3[1] - mu*Y3[1],
				  (gamma0 + ((gamma1 - gamma0)* b3) / (Y3[1] + b3))* Y3[1] - mu*Y3[2],
			   alpha*Y3[1] - delta*Y3[3] ])



Y03 = np.array([0.8, 0.15, 0.05, 0.40])        # initials conditions
Y3, infodict = integrate.odeint(dYdt3, Y03, t, full_output=True)
infodict['message']  
s3, i3, r3, x3 = Y3.T

plt.figure(1)
plt.plot(t, r1, 'g--',lw = 2, label='Recovered$\quad$individuals$\quad$ b = 0.302')
plt.plot(t, r2, 'b--',lw = 2, label='Recovered$\quad$individuals$\quad$ b = 0.521')
plt.plot(t, r3, 'r-.',lw = 2, label='Recovered$\quad$individuals$\quad$ b = 1.99')
plt.xlabel("Time (days)",fontsize = 16)
plt.ylabel("Recovered population",fontsize = 16)
plt.legend(loc = 'upper right')
plt.grid(True)
plt.savefig('RecovdESS_Vary.png')

plt.figure(2)
plt.plot(t, i1, 'g--',lw = 2, label= 'Infected$\quad$individuals$\quad$ b = 0.302')
plt.plot(t, i2, 'b--',lw = 2, label= 'Infected$\quad$individuals$\quad$ b = 0.521')
plt.plot(t, i3, 'r-.',lw = 2, label= 'Infected$\quad$individuals$\quad$ b = 1.990')
plt.xlabel("Time (days)",fontsize = 16)
plt.ylabel("Infected population",fontsize = 16)
plt.legend( loc ='upper right' , fontsize=10)
plt.grid(True)
plt.savefig('InfectESS_Vary.png')

plt.show()
