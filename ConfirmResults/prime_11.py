from mmgroup import *
import timeit
import sys
import itertools as it
from functools import reduce
import numpy as np
import time
import random
import math
from itertools import islice
import multiprocessing as mp

def group_generated_by(L, n, order_only=False):     
    """Generates group from array of elements. 
    Written by Dietrich, Lee, and Popiel while writing the paper on the completion of the classification of the maximal subgroups."""
    start = time.time()
    orb = [L[0]]
    orbset = {tuple(L[0].as_tuples())}
    os = 0
    for el in L:
        eltup = tuple(el.as_tuples())
        if not eltup in orbset:
            orb.append(el)
            orbset.add(eltup)
            os = os+1;
             
    j = 0
    while j <= os-1:
        for g in L:
            el= orb[j]*g
            eltup = tuple(el.as_tuples())
            if not eltup in orbset:
                orb.append(el)
                orbset.add(eltup)
                os = os+1;
                    
        j = j+1
        end = time.time()       
        print("Limit", n, "; have", os, "in time ", round(end-start,4), end='\r')
        if len(orbset)>n:
            print("Group is larger than imposed limit -- abort in time ", round(end-start,4) )
            return []
    
    print("Limit", n, "; have", len(orb), "in time ", round(end-start,4))
    if order_only:
        return len(orb)
    return orb

gens_11 = [
    MM("M<y_79h*x_0cc9h*d_3a0h*p_191654485*l_2*p_2880*l_1*p_24240*l_1*p_624000>"),
    MM("M<y_107h*x_1a47h*d_0a7h*p_44352107*l_2*p_2830080*l_2*p_32128794*l_1*t_2*l_2*p_1900800*l_2*p_23202560*l_2*t_1*l_2*p_2386560*l_2*p_33015896*l_2*t_2*l_2*p_1985280*l_1*p_85329047*t_2*l_1*p_1499520*l_2*p_43158128*t_1*l_2*p_2597760*l_1*p_37670*l_2*t_2*l_1*p_1499520*l_2*p_42664548*t_1*l_2*p_1943040*l_2*p_64002560>")
]

print("The generators commute:", gens_11[0] * gens_11[1] == gens_11[1] * gens_11[0])

G = group_generated_by(gens_11, n=11**2)

print("The order of the generated group is correct:", len(G) == 11**2)



# Output:

# The generators commute: True
# Limit 121 ; have 121 in time  3.4509
# The order of the generated group is correct: True