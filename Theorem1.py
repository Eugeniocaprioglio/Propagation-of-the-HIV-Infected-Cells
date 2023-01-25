import numpy as np
import math
import array as arr
import pickle

# Defining Variables into a Dictionary

max_dict = {
    "r": [0, 100],
    "n": [1, 300],
    "alpha": [2, 1.5],
    "gamma": [3, 0.001],
    "s_max": [4, 1500],
    "mus": [5, 0.1],
    "mui": [6, 0.5],
    "muv": [7, 10]
    }

# SAVING DICTIONARY TO A FILE  --------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

np.save('GlobalVarDictionary.npy', max_dict) # max_dict.npy
txt_filename = 'GlobalVariablesDictionary.txt' # max_dict.txt (Human eye)

with  open(txt_filename, 'w') as f :

    f.write(str(max_dict))

# BUILDING GENERATOR FUNCTION (used for simulation) ---------------------------------------------
# ------------------------------------------------------------------------------------------------------------


def gen(a,k): # k = max_dict[...][1]
    return abs(math.sin(a**2*k+a*k*3-k**4*a+36))

# ------------------------------------------ SYSTEM STATE FUNCTIONS
# ------------------------------------------------------------------------------------------------------------

def uninfected_steady_state(a):

    alpha = max_dict["alpha"][1]*gen(a,max_dict["alpha"][0])
    mus = max_dict["mus"][1]*gen(a,max_dict["mus"][0])
    r = max_dict["r"][1]*gen(a,max_dict["r"][0])
    s_max = max_dict["s_max"][1]*gen(a,max_dict["s_max"][0])

    return [(r - mus + math.sqrt((r - mus)**2 + (4*alpha*r/s_max)))*s_max/(2*r),
    0,
    0]
# ------------------------------------------------------------------------------------------------------------

def infected_steady_state(a):

    alpha = max_dict["alpha"][1]*gen(a,max_dict["alpha"][0])
    mus = max_dict["mus"][1]*gen(a,max_dict["mus"][0])
    r = max_dict["r"][1]*gen(a,max_dict["r"][0])
    s_max = max_dict["s_max"][1]*gen(a,max_dict["s_max"][0])
    gamma = max_dict["gamma"][1]*gen(a,max_dict["gamma"][0])
    mui = max_dict["mui"][1]*gen(a,max_dict["mui"][0])
    muv = max_dict["muv"][1]*gen(a,max_dict["muv"][0])
    n = max_dict["n"][1]*gen(a,max_dict["n"][0])

    return [muv/(gamma*n),
    alpha/mui - (mus*muv)/(gamma*mui*n) + ((muv*r)/(gamma*mui*n))*(1-(muv/(gamma*n*s_max))),
    alpha*n/muv - mus/gamma + ((r/gamma)*(1-(muv/(gamma*n*s_max)))) ]

# ------------------------------------------------------------------------------------------------------------

def reproduction_ratio(a):

    gamma = max_dict["gamma"][1]*gen(a,max_dict["gamma"][0])
    n = max_dict["n"][1]*gen(a,max_dict["n"][0])
    muv = max_dict["muv"][1]*gen(a,max_dict["muv"][0])
    s_0 = uninfected_steady_state(a)[0]

    return (gamma*n*s_0)/muv

# ---------------------------------------------------------------------------------------
