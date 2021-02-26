from pandemic import *
from matplotlib import pyplot as plt

# Creation of instances with example parameters
corona = Pandemic() # name of the pandemic
COVID19 = Disease(3,0.1,365) # instance disease 
Germany = Population(83.9e6, 1e3, 0) # instance population
BNT = Vaccine(2e5,2,0.9) # instance vaccination

# Simulation
simulation = corona.SIRS(COVID19, Germany, BNT, 1800)

# Rename the output of the simulation
time = simulation[0]
incidence = simulation[1]
susceptibles = simulation[2]
infectious = simulation[3]
removed = simulation[4]

# Plot of the dynamics of the model
plt.grid(True)
plt.title('Dynamics in the SIRS-model with vaccination')
plt.xlabel('Time [days]')
plt.ylabel('People in Million')
plt.plot(time, susceptibles, 'blue', label='susceptible')
plt.plot(time, infectious, 'red', label='infectious')
plt.plot(time, removed, 'green', label='removed')
plt.legend()
plt.show()

# Plot of the incidence curve after the first wave
start = 300

plt.grid(True)
plt.title('Incidence curve in the SIRS-model with vaccination')
plt.xlabel('Time [days] after the first peak')
plt.ylabel('7-days incidence per 100k people')
plt.plot(time[start:], incidence[start:], label='incidence')
plt.legend()
plt.show()
