"""
Definition of class Pandemic
"""
# The class "Pandemic" contains the following subclasses: 
#   - "Disease" with the following parameters: 
#       - basic reproduction number 
#       - average duration of infectivity 
#       - average duration of immunity
#   - "Population" with the following parameters: 
#       - total population size
#       - number of actively infectious people
#       - number of removed people 
#   - "Vaccine" with the following parameters: 
#       - vaccination frequency, i.e. how many people are getting vaccinated each day
#       - number of required doeses for immunity
#       - the efficacy of the vaccination
#
# The class "Pandemic" contains a function "SIRS" which implements the SIRS-model 
# including vaccination.
# The input of the SIRS model are the instances of "Disease", "Population" and "Vaccine" and
# the duration of the simulation in days.     

class Pandemic(object):
    
    def SIRS(self, Disease, Population, Vaccine, duration_of_simulation):
        
        # Parameters of the instance of Disease
        R0 = Disease.R0
        gamma = Disease.gamma
        delta = Disease.delta
        
        # Parameters of the instance of Population
        N = Population.N 
        I = Population.I
        R = Population.R
        S = N - I - R # number of susceptibles
        
        # Defining lists for plots
        time = []
        S_list = []
        I_list = []
        R_list = []
        incidence = []
        
        # Parameters of the instance of 
        vaccination_frequency = Vaccine.vaccination_frequency 
        shots = Vaccine.number_of_doses 
        efficacy = Vaccine.efficacy 
        
        # Duration of the simulation
        duration = duration_of_simulation
        
        # Implementation of the SIRS-model with vaccination
        for i in range(0,duration):
            
            # Full vaccination possible only for S > vaccination_frequency
            daily_vaccination = min(S, vaccination_frequency)
            
            # Differential equations
            dS = -gamma*R0*I*S/N + delta*R - daily_vaccination/shots*efficacy
            dI = gamma*(R0*S/N-1)*I 
            I_new = gamma*R0*S*I/N # number of new infections for incidence
            # S, I and R must be non-negative
            S = max(S + dS,0)
            I = max(I + dI,0)
            R = max(N-S-I,0)
            # filling the lists
            time.append(i)
            S_list.append(S/1e6)
            I_list.append(I/1e6)
            R_list.append(R/1e6)
            incidence.append(I_new*7*1e5/N)
        
        return time, incidence, S_list, I_list, R_list
    
class Disease(Pandemic):
    
    def __init__(self, reproduction_number, duration_infectivity, duration_immunity):
        self.R0 = reproduction_number # basic reproduction number 
        self.gamma = duration_infectivity # average duration of infectivity in days
        self.delta = 1/duration_immunity # average duration of immunity in days
     
class Population(Pandemic):
    
    def __init__(self, population, infectious, removed):
        self.N = population # total population size
        self.I = infectious # number of actively infectious people
        self.R = removed # number of removed people
        
class Vaccine(Pandemic):
    
    def __init__(self, vaccination_frequency, number_of_doses, efficacy):
        self.vaccination_frequency = vaccination_frequency # vaccination frequency, i.e. how many people are getting vaccinated each day
        self.number_of_doses = number_of_doses # number of required doeses for immunity
        self.efficacy = efficacy # the efficacy of the vaccination