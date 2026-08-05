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

