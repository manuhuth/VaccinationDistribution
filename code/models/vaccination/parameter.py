import numpy as np

def general_set_up():
    non_vaccination_state = "vac0"
    virus_states = ["virW", "virM"]
    areas = ["countryA", "countryB"]
    species_comp = ["susceptible", "infectious", "recovered", "dead"]
    distances = np.array([[0, 3], [3, 0]])
    
    out = {"non_vaccination_state" : non_vaccination_state,
           "virus_states" : virus_states,
           "areas" : areas,
           "species_compartments" : species_comp,
           "distances" : distances,}
    return out
    
    
def fixed_parameter():
    parameter = {
    "beta": 10**-1.1519, #infection rate; parameter from Simon
    "lambda1": 1/22.61993, #rate of transition from infectious to recovered/dead; parameter from Simon
    "p": 0.04, #probability of dying p*lambda is fraction that dies in vac0
    "eta_virW": 1, #factor with which beta is scaled for the wild type
    "eta_virM": 1.3, #factor with which beta is scaled for the mutant type
    }
    return parameter


def parameters_vaccine_one():
    parameter = {"number_vac1": 400, #inflow of vaccine 1
    "omega_vac1_virW": 0.99, #if 1, vaccine1 gives full death protection against wild type if infected
    "delta_vac1_virW": 0.9, #if 1, vaccine1 gives full infection protection against wild type
    "omega_vac1_virM": 0.99, #if 1, vaccine1 gives full death protection against mutant type if infected
    "delta_vac1_virM": 0.9, #if 1, vaccine1 gives full infection protection against mutant type
    }
    
    return parameter

def parameters_vaccine_two():
    parameter = {
    "number_vac2": 100, #inflow of vaccine 2
    "omega_vac2_virW": 0.5, #if 1, vaccine2 gives full death protection against wild type if infected
    "delta_vac2_virW": 0.6, #if 1, vaccine2 gives full infection protection against wild type
    "omega_vac2_virM": 0.5, #if 1, vaccine1 gives full death protection against mutant type if infected
    "delta_vac2_virM": 0.6, #if 1, vaccine2 gives full infection protection against mutant type
    }
    
    return parameter

def start_parameter_one():
    parameter = {"susceptible_countryA_vac0_t0": 80000,
    "susceptible_countryB_vac0_t0": 80000,
    "infectious_countryA_vac0_virW_t0": 10,
    "infectious_countryA_vac0_virM_t0": 10,
    "infectious_countryB_vac0_virW_t0": 10,
    "infectious_countryB_vac0_virM_t0": 10,
    "infectious_countryA_vac1_virW_t0": 0.1,
    "infectious_countryA_vac1_virM_t0": 0.1,
    "infectious_countryB_vac1_virW_t0": 0.1,
    "infectious_countryB_vac1_virM_t0": 0.1,}
    
    return parameter

def start_parameter_two():
    parameter = {"infectious_countryA_vac2_virW_t0": 0.1,
    "infectious_countryA_vac2_virM_t0": 0.1,
    "infectious_countryB_vac2_virW_t0": 0.1,
    "infectious_countryB_vac2_virM_t0": 0.1,}
    
    return parameter