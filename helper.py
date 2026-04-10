import qutip
import numpy as np
import random
from collections import defaultdict,Counter
import itertools
from itertools import product
from collections import defaultdict, Counter
import cvxpy as cp
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from scipy.optimize import minimize
import warnings

def generate_binary_dict(n_min=1, n_max=20):
    result = {}
    for n in range(n_min, n_max + 1):
        all_keys = [''.join(p) for p in itertools.product('01', repeat=n)]
        result[f"all_keys_var_{n}"] = all_keys
    return result

def ket(*args):
    return qutip.tensor(*[qutip.basis(2, i) for i in args])

def bra(*args):
    return ket(*args).dag()

def kb(ket_str_bra_str):
    mid = len(ket_str_bra_str) // 2
    ket_str = ket_str_bra_str[:mid]
    bra_str = ket_str_bra_str[mid:]

    ket_state = ket(*[int(i) for i in ket_str])
    bra_state = bra(*[int(i) for i in bra_str])

    return ket_state * bra_state

def bk(bra_str_ket_str):
    mid = len(bra_str_ket_str) // 2
    ket_str = bra_str_ket_str[:mid]
    bra_str = bra_str_ket_str[mid:]
    ket_state = ket(*[int(i) for i in ket_str])
    bra_state = bra(*[int(i) for i in bra_str])

    return bra_state*ket_state

def project_onto_prob_simplex(p):
    n = len(p)
    def objective(x):
        return np.sum((x - p) ** 2)
    constraints = [
        {'type': 'eq', 'fun': lambda x: np.sum(x) - 1}
    ]
    bounds = [(0, None)] * n 
    x0 = np.maximum(p, 0)
    result = minimize(objective, x0, bounds=bounds, constraints=constraints)
    return result.x if result.success else None

def renormalize_counts(count_dict, target_total_count=10000):
    keys = list(count_dict.keys())
    raw_values = np.array([count_dict[k] for k in keys])
    total_raw = np.sum(raw_values)
    if total_raw == 0:
        raise ValueError("Total raw count is zero. Cannot normalize.")
    probs = raw_values / total_raw
    projected_probs = project_onto_prob_simplex(probs)
    if projected_probs is None:
        raise RuntimeError("Optimization failed during projection.")

    final_counts = {k: int(round(p * target_total_count)) for k, p in zip(keys, projected_probs)}
    return final_counts

def Counts_Readout_Dic(count_dic, num_qubit, readout_mat):
    if num_qubit < 2 or num_qubit > 20:
        raise ValueError("num_qubit must be between 2 and 20.")
    
    all_keys_var_dict=generate_binary_dict(2,6)

    if num_qubit == 2:
        all_keys = all_keys_var_dict[f"all_keys_var_{num_qubit}"]
    else:
        all_keys = all_keys_var_dict[f"all_keys_var_{num_qubit}"]  
    
    count_noisyarr = np.array([count_dic[key] for key in all_keys])
    count_idealarr = np.linalg.inv(readout_mat) @ count_noisyarr
    count_mit_dic = {}
    for idx, count in enumerate(count_idealarr):
        key = format(idx, f'0{num_qubit}b')
        count_mit_dic[key] = round(np.real(count))
    
    return count_mit_dic

def ideal_4q_GHZ():
    return qutip.Qobj(((kb('00000000')+kb('00001111')+kb('11110000')+kb('11111111'))/2).full().reshape(16,16))

def ideal_4q_0():
    return qutip.Qobj((kb('00000000')).full().reshape(16,16))

def ideal_4q_p():
    plus = (qutip.basis(2, 0) + qutip.basis(2, 1)).unit()
    plus4 = qutip.tensor(plus, plus, plus, plus)
    rho = plus4 * plus4.dag()
    return qutip.Qobj(rho.full().reshape(16,16))

def modify_last_4_bits(state, mask):
    # seperate first/last four bits
    front_bits = state[:4]
    back_bits = state[4:]
    
    # Flip last 4 bits
    new_back_bits = ''
    for b, m in zip(back_bits, mask):
        if m == 'X':  # Bit flip for X
            new_back_bits += '1' if b == '0' else '0'
        elif m == 'I':
            new_back_bits += b

    return front_bits + new_back_bits

