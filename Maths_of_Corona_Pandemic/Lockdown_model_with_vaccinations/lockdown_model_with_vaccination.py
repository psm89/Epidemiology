from containment_strategy import *
from matplotlib import pyplot as plt

# Definition of instances
pandemic = Containment_Strategy()
COVID19 = Disease(3, 0.1, 1e10) # note: last parameter = 1e10 to approximate eternal immunity
Germany = Population(83.9e6, 60.8, 129203, 2243892*2, 0.3, 0.7)
BNT = Vaccine(6e5, 2, 0.95, 28, 0, 1.7e6)
shutdown = Lockdown(215, 35, 1.6, 0.7, 'lockdown')

# Run the simulation
simulation = pandemic.dynamics_with_lockdown(COVID19,Germany,BNT,shutdown,360)

# Choice of suitable names for output of the simulation
time = simulation[0]
susceptibles = simulation[1]
incidence = simulation[2]
infectious = simulation[3]
removed = simulation[4]
vaccinated = simulation[5]
days_in_lockdown = simulation[6]

# Plot of the incidence curve. 
# Notice: Discontinuity occurs due to jump between effective reproduction numbers.
plt.grid(True)
plt.title('Incidence curve')
plt.xlabel('Time [Days]')
plt.ylabel('7-days Incidence per 100k')
plt.bar(time, incidence)
plt.show()

# Computation of further days in shutdown
print('Expected number of days in shutdown: ', sum(days_in_lockdown))
