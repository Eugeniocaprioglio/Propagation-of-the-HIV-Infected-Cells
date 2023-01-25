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

# SAVING DICTIONARY TO A FILE  -------------------------------------------------
# ------------------------------------------------------------------------------

np.save('GlobalVarDictionary.npy', max_dict) # max_dict.npy
txt_filename = 'GlobalVariablesDictionary.txt' # max_dict.txt (Human eye)

with  open(txt_filename, 'w') as f :

    f.write(str(max_dict))

# BUILDING GENERATOR FUNCTION (used for simulation) ----------------------------


def gen(a,k): # k = max_dict[...][1]
    return abs(math.sin(a**2*k+a*k*3-k**4*a+36))

# ------------------------------------------ SYSTEM STATE FUNCTIONS
# ------------------------------------------------------------------------------

def uninfected_steady_state(a):

    alpha = max_dict["alpha"][1]*gen(a,max_dict["alpha"][0])
    mus = max_dict["mus"][1]*gen(a,max_dict["mus"][0])
    r = max_dict["r"][1]*gen(a,max_dict["r"][0])
    s_max = max_dict["s_max"][1]*gen(a,max_dict["s_max"][0])

    return [(r - mus + math.sqrt((r - mus)**2 + (4*alpha*r/s_max)))*s_max/(2*r),
    0,
    0]
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------

def reproduction_ratio(a):

    gamma = max_dict["gamma"][1]*gen(a,max_dict["gamma"][0])
    n = max_dict["n"][1]*gen(a,max_dict["n"][0])
    muv = max_dict["muv"][1]*gen(a,max_dict["muv"][0])
    s_0 = uninfected_steady_state(a)[0]

    return (gamma*n*s_0)/muv

# ------------------------------------------------------------------------------

def n_parameters(a):

    alpha = max_dict["alpha"][1]*gen(a,max_dict["alpha"][0])
    mus = max_dict["mus"][1]*gen(a,max_dict["mus"][0])
    s_max = max_dict["s_max"][1]*gen(a,max_dict["s_max"][0])
    gamma = max_dict["gamma"][1]*gen(a,max_dict["gamma"][0])
    mui = max_dict["mui"][1]*gen(a,max_dict["mui"][0])
    muv = max_dict["muv"][1]*gen(a,max_dict["muv"][0])

    return [(gamma**2)*mui*s_max - 4*alpha*(mui + muv),
            - 2*gamma*mui*muv*((mui**2 + 3*mui*muv + muv**2)*s_max - 4*alpha*(mui + muv)),
            muv*(mui**4 + 2*(mui**2)*muv*(3*mui - 2*mus) + 6*mui*(muv**3) + mui*(muv**2)*(11*mui - 4*mus) + muv**4)
           ]
# ------------------------------------------------------------------------------
def r_parameters(a):

    alpha = max_dict["alpha"][1]*gen(a,max_dict["alpha"][0])
    mus = max_dict["mus"][1]*gen(a,max_dict["mus"][0])
    s_max = max_dict["s_max"][1]*gen(a,max_dict["s_max"][0])
    gamma = max_dict["gamma"][1]*gen(a,max_dict["gamma"][0])
    mui = max_dict["mui"][1]*gen(a,max_dict["mui"][0])
    muv = max_dict["muv"][1]*gen(a,max_dict["muv"][0])
    n = max_dict["n"][1]*gen(a,max_dict["n"][0])

    return [(muv**4)*(mui + muv),
    gamma*(muv**2)*n*s_max*(-gamma*mui*muv*n*s_max + 2*alpha*gamma*mui*n + 2*alpha*gamma*muv*n + (mui**2)*muv + muv**3),
    (n**2)*(gamma**2)*(s_max**2)*(alpha*gamma*(muv**3)*n + alpha*gamma*mui*muv*n*(mui + muv) + (alpha**2)*(gamma**2)*(n**2)*(mui + muv) + mui*mus*muv**3),
    ]

# ------------------------------------------------------------------------------

    def state_is_in_P(a):

    n = max_dict["n"][1]*gen(a,max_dict["n"][0])
    r = max_dict["r"][1]*gen(a,max_dict["r"][0])
    [an, bn, cn] = n_parameters(a)
    [ar, br, cr] = r_parameters(a)
    root_1 = (-bn + (bn**2 -4*an*cn)**0.5)/(2*an)
    root_2 = (-bn - (bn**2 -4*an*cn)**0.5)/(2*an)

    max_root = max(root_1, root_2)

    a = (n > max_root) # boolean: true/false

    b = (ar*r**2+br*r+cr < 0)

    if (a == True and b==True):
        return True

    return False

def choose_in_p():

    arr = []
    for a in range(12):
        if state_is_in_P(a):
            arr.append(a)

    return arr

#print(choose_in_p()) # print "a" parameters for which (N,r) belongs to set P

state_is_in_p_array = choose_in_p()

# SYSTEM STATE MATRICES --------------------------------------------------------
# and data saving --------------------------------------------------------------

# BUILDING MATRIX FOR STEADY STATE S_0, I_0=0, V_0=0 ---------------------------

row = len(choosen_a_array)
col = 3

steady_matrix = np.zeros((row, col))

for i in range(0, row):
    steady_matrix[i][0] = uninfected_steady_state(choosen_a_array[i])[0]

# BUILDING MATRIX FOR STEADY STATE S_1, I_1, V_1 only for r > 1 (state is in P)

rowP = len(state_is_in_p_array)
colP = 3

steady_P_matrix = np.zeros((rowP, colP))

for i in range (0, len(state_is_in_p_array)):
    steady_P_matrix [i, :] = infected_steady_state(state_is_in_p_array[i])

# Saving values for which state is (not) in P into a file ------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

arr = np.asarray(choosen_a())
np.save("state_is_not_in_P.npy", arr) # binary data
np.savetxt('state_is_not_in_P.txt', arr) # human readable data

arr2 = np.asarray(choose_in_p())
np.save("state_is_in_P.npy", arr2) # binary data
np.savetxt('state_is_in_P.txt', arr2) # human readable data

a_used = np.concatenate((arr, arr2))

print (a_used)

# Trying to save matrix S_0/I_0/V_0 to file ------------------------------------
# ------------------------------------------------------------------------------

#Binary data
np.save('S_0_matrix.npy', steady_matrix)

#Human readable data
np.savetxt('S_0_matrix.txt', steady_matrix)

# Trying to save matrix S_1/I_1/V_1 to binary file -----------------------------
# ------------------------------------------------------------------------------

#Binary data
np.save('S_P_matrix.npy', steady_P_matrix)

#Human readable data
np.savetxt('S_P_matrix.txt', steady_P_matrix)

# Saving dictionary values for all used a (5 values used) ----------
# ------------------------------------------------------------------

alpha_array = np.zeros(len(a_used))
mus_array = np.zeros(len(a_used))
r_array = np.zeros(len(a_used))
s_max_array = np.zeros(len(a_used))
gamma_array = np.zeros(len(a_used))
mui_array = np.zeros(len(a_used))
muv_array = np.zeros(len(a_used))
n_array = np.zeros(len(a_used))

    
for i in range(len(a_used)): #col

        alpha_array[i] = max_dict["alpha"][1]*gen(a_used[i],max_dict["alpha"][0])
        mus_array[i] = max_dict["mus"][1]*gen(a_used[i],max_dict["mus"][0])
        r_array[i] = max_dict["r"][1]*gen(a_used[i],max_dict["r"][0])
        s_max_array[i] = max_dict["s_max"][1]*gen(a_used[i],max_dict["s_max"][0])
        gamma_array[i] = max_dict["gamma"][1]*gen(a_used[i],max_dict["gamma"][0])
        mui_array[i] = max_dict["mui"][1]*gen(a_used[i],max_dict["mui"][0])
        muv_array[i] = max_dict["muv"][1]*gen(a_used[i],max_dict["muv"][0])
        n_array[i] = max_dict["n"][1]*gen(a_used[i],max_dict["n"][0])

params_used_array = np.zeros(len(max_dict)*len(a_used))

row_used = len(max_dict)
col_used = len(a_used)
params_used_matrix = np.zeros ((row_used,col_used))

params_used_matrix[0,:] = alpha_array
params_used_matrix[1,:] = mus_array
params_used_matrix[2,:] = r_array
params_used_matrix[3,:] = s_max_array
params_used_matrix[4,:] = gamma_array
params_used_matrix[5,:] = mui_array
params_used_matrix[6,:] = muv_array
params_used_matrix[7,:] = n_array

# Binary data
np.save('params_used_matrix.npy', params_used_matrix)

# Human readable data
np.savetxt('params_used_matrix.txt', params_used_matrix)
