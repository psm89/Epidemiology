#---------------------------------------------------
# Implementation of a containment strategy
#---------------------------------------------------

# Idea: 
# Containment strategy is a shutdown whenever the 7-days incidence per 100k exceeds some critical value.
# The shutdown terminates once the incidence has dropped below some target value.
# The class Containment_Strategy contains a function dynamics_with_lockdown which implements the strategy.
# The subclasses Disease, Population, Vaccine, Lockdown are used to define corresponding instances with desired parameters.

class Containment_Strategy(object):
    
    def dynamics_with_lockdown(self, Disease, Population, Vaccine, Lockdown, duration_of_simulation):
        
        # Parameters of the instance of Disease
        R0 = Disease.R0
        gamma = Disease.gamma
        delta = Disease.delta
        
        # Parameters of the instance of Population
        N = Population.N 
        incidence = Population.incidence
        I = Population.I
        R = Population.R
        S = N - I - R # number of susceptibles
        high_risk_group = Population.high_risk_group
        vaccination_readiness = Population.vaccination_readiness
        
        # Parameters of the instance of Vaccine
        vaccination_frequency = Vaccine.vaccination_frequency 
        shots = Vaccine.number_of_doses 
        efficacy = Vaccine.efficacy 
        delay = Vaccine.delay
        begin_vaccination = Vaccine.begin_vaccination
        V = Vaccine.initially_vaccinated
        
        V_max = S*vaccination_readiness
        vacc = vaccination_frequency/shots 
        S_hr_fail = S*high_risk_group*((1-vaccination_readiness) + vaccination_readiness*(1-efficacy))
        
        # Parameters of Lockdown
        critical_incidence = Lockdown.critical_incidence*(N/7e5)
        subcritical_incidence = Lockdown.subcritical_incidence*(N/7e5)
        R_normal = Lockdown.R_normal
        R_lockdown = Lockdown.R_lockdown
        mode = Lockdown.mode
        
        # effective critical incidence 
        critical_incidence_eff = critical_incidence*high_risk_group
        
        # Duration of the simulation
        duration = duration_of_simulation
        
        # definition of lists for the purpose of plotting
    
        V_list = [V]
        S_list = [S]
        I_list = [I]
        R_list = [R]
        New_list = [incidence]
        time = [0]
        
        # for counting days in lockdown 
        count_lockdowns = 0
        count_days_in_lockdown = 0
        list_lockdowns = []
        days_in_lockdown = []
        
        for i in range (1, duration):
        
            if mode == 'lockdown':
                R0 = R_lockdown # effective reproduction number during lockdown
                if count_lockdowns == 0:
                    count_lockdowns = 1 
            else:
                R0 = R_normal # (initial) basis reproduction number
            
            New = gamma*R0*S*I/N # new cases 
            New_hr = New*high_risk_group # new cases in the high risk group
            
            # criterion for lockdown
            if (New_hr > critical_incidence_eff) and (mode != 'lockdown'):
                mode = 'lockdown'
                count_lockdowns += 1
            
            # counting the total days of lockdown
            if mode == 'lockdown':
                count_days_in_lockdown += 1 
            
            # criterion for switching from lockdown to normal mode
            if (New < subcritical_incidence) and (mode == 'lockdown'):
                R0 = R_normal
                mode = 'normal'
                list_lockdowns.append(count_days_in_lockdown)
            
            # counting number of people receiving vaccination for each day 
            # the criterions accounts for temporal delay in the protection
            if (V < V_max) and (i > delay + begin_vaccination): 
                dV = min(vacc,S)
            else:
                dV = 0
            
            dS = - gamma*R0*I*S/N - efficacy*dV + delta*(R+V)# change in the susceptibles
            dI = gamma*(R0*S/N-1)*I # change in the number of infectious people 
            dR = gamma*I - delta*R# change in the number of recovered people accounting for deaths
            # updating those variables
            V = V + dV - delta*V
            S = max(S + dS, 0)
            I = max(I + dI, 0)
            R = max(R + dR, 0)
            # inserting the variables to the corresponding lists
            V_list.append(V)
            S_list.append(S)
            I_list.append(I)
            R_list.append(R)
            New_list.append(New*7*1e5/N)
            time.append(i)
                    
            S_hr = high_risk_group*S # calculating the size of susceptible population among high risk group at the beginning of day i
            S_hr = S_hr - dI*high_risk_group # subtracting hr people who got infected before vaccination
    
            # this criterion ensures that the high risk group gets vaccination first
            if S_hr > S_hr_fail:
                S_hr = S_hr - dV + delta*V + delta*R # size of the susceptible high risk group at the end of day i
            
            high_risk_group = max(S_hr,0)/S # updating the proportion of the high risk group
        
        # counting the days for each lockdown
        for i in range(len(list_lockdowns)):
            if i == 0:
                days_in_lockdown.append(list_lockdowns[0])
            else:
                days_in_lockdown.append(list_lockdowns[i]-list_lockdowns[i-1])
                
        return time, S_list, New_list, I_list, R_list, V_list, days_in_lockdown
    
class Disease(Containment_Strategy):
    
    def __init__(self, reproduction_number, duration_infectivity, duration_immunity):
        self.R0 = reproduction_number # basic reproduction number 
        self.gamma = duration_infectivity # average duration of infectivity in days
        self.delta = 1/duration_immunity # average duration of immunity in days
     
class Population(Containment_Strategy):
    
    def __init__(self, population, incidence, infectious, removed, high_risk_group, vaccination_readiness):
        self.N = population # total population size
        self.incidence = incidence # 7-days incidence per 100k
        self.I = infectious # number of actively infectious people
        self.R = removed # number of removed people
        self.high_risk_group = high_risk_group # proportion of the high risk group
        self.vaccination_readiness = vaccination_readiness # proportion of people willing to get vaccinated
        
class Vaccine(Containment_Strategy):
    
    def __init__(self, vaccination_frequency, number_of_doses, efficacy, delay, begin_vaccination, initially_vaccinated):
        self.vaccination_frequency = vaccination_frequency # vaccination frequency, i.e. how many people are getting vaccinated each day
        self.number_of_doses = number_of_doses # number of required doeses for immunity
        self.efficacy = efficacy # the efficacy of the vaccination
        self.delay = delay # time until full protection in days
        self.begin_vaccination = begin_vaccination # counter for starting day of vaccination (zero if vaccination has already started)
        self.initially_vaccinated = initially_vaccinated # number of people who are initially vaccinated
        
class Lockdown(Containment_Strategy):
    
    def __init__(self, critical_incidence, subcritical_incidence, R_normal, R_lockdown, mode):
        self.critical_incidence = critical_incidence # critical incidence at which a shutdown is imposed
        self.subcritical_incidence = subcritical_incidence # target incidence at which shutdown terminates 
        self.R_normal = R_normal # effective reproduction number outside the shutdown
        self.R_lockdown = R_lockdown # effective reproduction number during the shutdown
        self.mode = mode # mode = 'lockdown' if currently in shutdown, otherwise 'normal'