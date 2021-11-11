import numpy as np


def general_set_up():
    non_vaccination_state = "vac0"
    virus_states = ["virW", "virM"]
    areas = ["countryA", "countryB"]
    species_comp = ["susceptible", "infectious", "recovered", "dead"]
    distances = np.array([[0, 4999], [4999, 0]])

    out = {
        "non_vaccination_state": non_vaccination_state,
        "virus_states": virus_states,
        "areas": areas,
        "species_compartments": species_comp,
        "distances": distances,
    }
    return out


def fixed_parameter():
    parameter = {
        "beta": 0.3,  # infection rate;
        "lambda1": 1
        / 10,  # rate of transition from infectious to recovered/dead; RKI https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Steckbrief.html;jsessionid=66D23F33F076AD1AA611E1C2159F45F7.internet052?nn=13490888#doc13776792bodyText10
        "prob_deceasing": 0.025,  # probability of dying p*lambda is fraction that dies in vac0; https://www.nature.com/articles/s41598-021-84055-6
        "eta_virW": 1,  # factor with which beta is scaled for the wild type
        "eta_virM": 1.2,  # factor with which beta is scaled for the mutant type
    }
    return parameter


def parameters_vaccine_one():
    parameter = {
        "omega_vac1_virW": 0.99,  # if 1, vaccine1 gives full death protection against wild type if infected
        "delta_vac1_virW": 0.94,  # if 1, vaccine1 gives full infection protection against wild type https://www.nejm.org/doi/full/10.1056/nejmoa2034577, https://www.mdpi.com/2076-393X/9/1/30/pdf
        "omega_vac1_virM": 0.99,  # if 1, vaccine1 gives full death protection against mutant type if infected
        "delta_vac1_virM": 0.9,  # if 1, vaccine1 gives full infection protection against mutant type
    }

    return parameter


def parameters_vaccine_two():
    parameter = {
        "omega_vac2_virW": 0.99,  # if 1, vaccine2 gives full death protection against wild type if infected
        "delta_vac2_virW": 0.7,  # if 1, vaccine2 gives full infection protection against wild type
        "omega_vac2_virM": 0.99,  # if 1, vaccine2 gives full death protection against mutant type if infected
        "delta_vac2_virM": 0.65,  # if 1, vaccine2 gives full infection protection against mutant type
    }

    return {**parameters_vaccine_one(), **parameter}


def parameters_vaccine_one_toy():
    parameter = {
        "omega_vac1_virW": 0.99,  # if 1, vaccine1 gives full death protection against wild type if infected
        "delta_vac1_virW": 0.9,  # if 1, vaccine1 gives full infection protection against wild type https://www.nejm.org/doi/full/10.1056/nejmoa2034577, https://www.mdpi.com/2076-393X/9/1/30/pdf
        "omega_vac1_virM": 0,  # if 1, vaccine1 gives full death protection against mutant type if infected
        "delta_vac1_virM": 0,  # if 1, vaccine1 gives full infection protection against mutant type
    }

    return parameter


def parameters_vaccine_two_toy():
    parameter = {
        "omega_vac2_virW": 0,  # if 1, vaccine2 gives full death protection against wild type if infected
        "delta_vac2_virW": 0,  # if 1, vaccine2 gives full infection protection against wild type
        "omega_vac2_virM": 0.99,  # if 1, vaccine2 gives full death protection against mutant type if infected
        "delta_vac2_virM": 0.9,  # if 1, vaccine2 gives full infection protection against mutant type
    }

    return {**parameters_vaccine_one_toy(), **parameter}


def start_parameter_one():
    parameter = {
        "susceptible_countryA_vac0_t0": 80000000,
        "susceptible_countryB_vac0_t0": 80000000,
        "infectious_countryA_vac0_virW_t0": 10,
        "infectious_countryA_vac0_virM_t0": 0,
        "infectious_countryB_vac0_virW_t0": 0,
        "infectious_countryB_vac0_virM_t0": 10,
        "infectious_countryA_vac1_virW_t0": 0,
        "infectious_countryA_vac1_virM_t0": 0,
        "infectious_countryB_vac1_virW_t0": 0,
        "infectious_countryB_vac1_virM_t0": 0,
    }

    return parameter


def start_parameter_two():
    parameter = {
        "infectious_countryA_vac2_virW_t0": 0,
        "infectious_countryA_vac2_virM_t0": 0,
        "infectious_countryB_vac2_virW_t0": 0,
        "infectious_countryB_vac2_virM_t0": 0,
    }

    return {**start_parameter_one(), **parameter}
