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

x3 = MM("M<y_400h*x_19ddh*d_0c5eh*p_94056107*l_1*p_1415040*l_2*p_1296*l_1*p_497280>")

gens_normal_3 = [
    MM("M<y_4c6h*x_47dh*d_62ah*p_210482411*l_2*p_1457280*l_1*p_2415843*t_1*l_1*p_23040*l_1*p_12103056*t_1*l_1*p_2640000*l_1*p_11539*t_2*l_1*p_182328960*t_2*l_2*p_1900800*l_2*p_14425*t_2*l_2*p_2787840*l_2*p_1018898>"), 
    MM("M<y_4f1h*x_920h*d_16ah*p_185054066*l_2*p_2344320*l_2*p_11601816*l_1*t_2*l_2*p_2830080*l_2*p_23197940*l_1*t_2*l_1*p_960*l_2*p_465888*l_2*p_338880*t_1*l_2*p_2830080*l_2*p_21334177*t_2*l_2*p_2344320*l_2*p_22767713*l_1*t_2*l_2*p_2386560*l_2*p_42716530*t_2*l_1*p_1499520*l_1*p_42670273>"), 
    MM("M<y_0d3h*x_0c05h*d_0a34h*p_169910595*l_1*p_1499520*l_1*p_2791974*t_1*l_1*p_1858560*l_2*p_21416304*t_1*l_1*p_1499520*l_2*p_64043947*t_2*l_2*p_2344320*l_2*p_32551155*l_2*t_1*l_1*p_2027520*l_1*p_4802*t_1*l_1*p_6042240*l_1>"), 
    MM("M<y_28ch*x_16dch*d_3d7h*p_166890611*l_1*p_1858560*l_2*p_3120*l_1*p_109440*t_1*l_2*p_1393920*l_2*p_10709088*t_1*l_2*p_1985280*l_1*p_21374144*t_2*l_2*p_2597760*l_1*p_13060339*l_2*t_2*l_2*p_2386560*l_2*p_42675072*t_1*l_2*p_2956800*l_1*p_64016147>"), 
    MM("M<y_123h*x_0b01h*d_6f1h*p_109083059*l_1*p_360000*l_2*t_2*l_2*p_2597760*l_1*p_21424*l_1*t_2*l_1*p_466560*l_2*p_21923040*l_1*t_1*l_2*p_2386560*l_2*p_32955362*l_1*t_2*l_1*p_1457280*l_2*p_9601*t_1*l_2*p_1985280*l_1*p_42736709>"), 
    MM("M<y_7e3h*x_1fb0h*d_24fh*p_144712389*l_2*p_1394880*l_2*p_24240*l_1*p_1907520*t_1*l_2*p_1900800*l_2*p_32507011*l_2*t_2*l_2*p_2956800*l_1*p_32956343*l_2*t_2*l_2*p_2787840*l_2*p_501256*t_1*l_1*p_3840*l_1*p_22416*l_1*p_3888960*t_2*l_2*p_1985280*l_1*p_42664514*t_2*l_2*p_1985280*l_1*p_53350562>"),
    MM("M<y_93h*x_2a2h*d_0acch*p_66086744*l_2*p_1393920*l_2*p_951840*t_2*l_2*p_2830080*l_2*p_21889554*t_2*l_1*p_2027520*l_1*p_6738*t_1*l_1*p_1457280*l_2*p_43621*t_2*l_2*p_2880*l_2*p_21360*l_2*p_6127680*t_1*l_1*p_1499520*l_2*p_42706969*t_2*l_1*p_3400320*l_2>"),
    MM("M<y_332h*x_19bch*d_0b5h*p_196521696*l_2*p_255360*l_1*t_2*l_1*p_1457280*l_2*p_1487563*l_2*t_1*l_2*p_2386560*l_2*p_32532131*l_1*t_2*l_2*p_2597760*l_1*p_43600690*t_1*l_1*p_2027520*l_1*p_482029*l_1*t_2*l_1*p_1499520*l_2*p_42834827*t_2*l_1*p_1866240>"),
    MM("M<y_7d7h*x_0d09h*d_0d17h*p_133224044*l_1*p_2027520*l_1*p_32003699*l_1*t_2*l_2*p_1457280*l_1*p_2860305*l_2*t_2*l_2*p_2830080*l_2*p_465882*l_1*t_2*l_2*p_2956800*l_1*p_2355408*l_1*t_1*l_2*p_1900800*l_2*p_10581*t_1*l_2*p_1943040*l_2*p_42717409>"),
    MM("M<y_4h*x_59ch*d_285h*p_117194863*l_2*p_2597760*l_1*p_32068419*l_1*t_1*l_1*p_1499520*l_2*p_64004448*t_1*l_1*p_1920*l_2*p_528384*t_2*l_1*p_14657280*l_2*p_53708160*t_1*l_1*p_1394880*l_2*p_24240*l_1*p_2031360*t_1*l_2*p_2386560*l_2*p_21336098>"),
    MM("M<y_0c6h*x_0f59h*d_8c7h*p_214008324*l_2*p_1985280*l_1*p_32027747*l_2*t_2*l_2*p_2597760*l_1*p_27940*t_1*l_2*p_1457280*l_1*p_66387*t_1*l_1*p_1499520*l_1*p_10693537*t_2*l_1*p_2880*l_1*p_21408*l_1*p_360000*t_1*l_2*p_2386560*l_2*p_33434373>"),
    MM("M<y_360h*x_78dh*d_0a3ah*p_102516690*l_1*p_2640000*l_1*p_43592096*t_1*l_1*p_464640*l_2*p_1858800*l_2*t_2*l_1*p_1393920*l_2*p_22272*l_1*p_492480*t_1*l_1*p_31680*l_1*t_2*l_1*p_1499520*l_2*p_106703620*t_2*l_2*p_2956800*l_1*p_53821843*t_2*l_2*p_1943040*l_2*p_85816850*t_2*l_2*p_2830080*l_2*p_106703619>")
]

cyclic = []
cyclic_as_int = []

curr = x3
for _ in range(3):
    cyclic.append(curr)
    cyclic_as_int.append(curr.as_int())
    curr *= x3


def commutator(a,b):
    return a**-1 * b**-1 * a * b 


print("Checking the normal subgroup 3^{1+12}")
print()



Failed = False
for i in range(len(gens_normal_3)):         #Verify that all generators commute mod x5
    for j in range(i+1, len(gens_normal_3)):
        a = gens_normal_3[i]
        b = gens_normal_3[j]
        if commutator(a,b).as_int() not in cyclic_as_int:
            Failed = True

print("Verifying that all elements commute mod x3:", not Failed)
print()

gens_abc = gens_normal_3[:6]
gens_xyz = gens_normal_3[6:]


G_abc = group_generated_by(gens_abc, n=3**7)            #Generate partial groups
G_xyz = group_generated_by(gens_xyz, n=3**7)

G_abc_as_int = [g.as_int() for g in G_abc]              #int comparisons are much faster than mmgroup object comparisons
G_xyz_as_int = [g.as_int() for g in G_xyz]


Failed = False
for g_abc_as_int in G_abc_as_int:
    if g_abc_as_int in G_xyz_as_int and g_abc_as_int not in cyclic_as_int:
        Failed = True


print("The intersection in trivial modulo x3:", not Failed)


# Output:

# Checking the normal subgroup 3^{1+12}

# Verifying that all elements commute mod x3: True

# Limit 2187 ; have 2187 in time  251.2124
# Limit 2187 ; have 2187 in time  257.9827
# The intersection in trivial modulo x3: True