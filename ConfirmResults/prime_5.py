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


g2_G = MM("M<y_4f6h*x_1f98h*d_0b7h*p_67615847*l_1*p_2999040*l_1*p_86264262*l_2*p_11172480>")
g3_G = MM("M<y_4e1h*x_19cbh*d_9c8h*p_19643307*l_1*p_2999040*l_1*p_64003504*l_2*p_1478400>")
g5 = g2_G*g3_G

gens_normal_5 = [
    MM("M<y_368h*x_1352h*d_39ah*p_180576707*l_2*p_2830080*l_2*p_96040849>"),
    MM("M<y_77eh*x_104ah*d_386h*p_231348305*l_1*p_2640000*l_1*p_11323780*t_2*l_2*p_2830080*l_2*p_12579258*l_1*t_2*l_2*p_960*l_2*p_22368*l_2*p_617280*t_1*l_2*p_2787840*l_2*p_603280*t_2*l_1*p_1858560*l_1*p_22368*l_2*p_193920*t_1*l_2*p_1858560*l_2*p_85367712*t_2*l_2*p_1943040*l_2*p_86256660>"),
    MM("M<y_731h*x_15d2h*d_34dh*p_177049707*l_1*p_1499520*l_1*p_507205*t_1*l_2*p_1900800*l_2*p_23219987*t_1*l_1*p_1394880*l_1*p_42789696*t_2*l_2*p_2787840*l_2*p_7698*t_1*l_2*p_1394880*l_1*p_10665888*l_2*p_6105600*t_1*l_2*p_2597760*l_1*p_85328070*t_2*l_2*p_1943040*l_2*p_86277651>"),
    MM("M<y_9bh*x_809h*d_0b2ch*p_58396269*l_1*p_1499520*l_1*p_2261*l_2*t_2*l_1*p_2999040*l_1*p_32065323*t_1*l_2*p_1985280*l_1*p_63995944*t_1*l_2*p_960*l_2*p_466896*l_2*p_359040*t_1*l_1*p_1499520*l_2*p_32082597>"),
    MM("M<y_11bh*x_1edfh*d_974h*p_24562111*l_2*p_1457280*l_1*p_21340951*t_1*l_2*p_1985280*l_1*p_1463362*l_2*t_2*l_2*p_1393920*l_1*p_23280*l_2*p_10437120*t_1*l_1*p_1457280*l_2*p_946672*t_1*l_2*p_2830080*l_2*p_42838662*t_1*l_2*p_2956800*l_1*p_43239926*t_1*l_2*p_2956800*l_1*p_64123810>"),
    MM("M<y_50dh*x_15aah*d_1a9h*p_150943776*l_2*p_2344320*l_2*p_21859494*l_1*t_2*l_1*p_2999040*l_1*p_43182182*t_2*l_2*p_1457280*l_1*p_517586*t_2*l_2*p_2830080*l_2*p_21331876*t_2*l_2*p_1943040*l_2*p_85413728*t_1*l_2*p_2597760*l_1*p_11307408*t_2*l_2*p_1985280*l_1*p_53802562>")
]

cyclic = []
cyclic_as_int = []

curr = g5
for _ in range(5):
    cyclic.append(curr)
    cyclic_as_int.append(curr.as_int())
    curr *= g5


def commutator(a,b):
    return a**-1 * b**-1 * a * b 


print("Checking the normal subgroup 5^{1+6}")
print()



Failed = False
for i in range(len(gens_normal_5)):         #Verify that all generators commute mod x5
    for j in range(i+1, len(gens_normal_5)):
        a = gens_normal_5[i]
        b = gens_normal_5[j]
        if commutator(a,b).as_int() not in cyclic_as_int:
            Failed = True

print("Verifying that all generators commute mod x5:", not Failed)
print()

gens_abc = gens_normal_5[:3]
gens_xyz = gens_normal_5[3:]


G_abc = group_generated_by(gens_abc, n=5**4)            #Generate partial group s
G_xyz = group_generated_by(gens_xyz, n=5**4)

G_abc_as_int = [g.as_int() for g in G_abc]              #int comparisons are much faster than mmgroup object comparisons
G_xyz_as_int = [g.as_int() for g in G_xyz]


Failed = False
for g_abc_as_int in G_abc_as_int:
    if g_abc_as_int in G_xyz_as_int and g_abc_as_int not in cyclic_as_int:
        Failed = True


print("The intersection is trivial modulo x5:", not Failed)
print()


print("Checking the complement 2 * J_2:4")
print()

gens_comp_5 = [
    MM("M<y_387h*x_68ch*d_6c2h*p_48364112*l_1*p_1499520*l_2*p_23197937*t_2*l_2*p_464640*l_1*p_1858896*l_1*t_2*l_1*p_1858560*l_2*p_21408*l_1*p_340800*t_1*l_2*p_5641920*l_1*t_1*l_2*p_2386560*l_2*p_21334307*t_2*l_2*p_2597760*l_1*p_106706512*t_2*l_1*p_1457280*l_2*p_971769*t_2*l_2*p_2956800*l_1*p_43160098>"),
    MM("M<y_608h*x_1496h*d_0c11h*p_186954798*l_1*p_1499520*l_1*p_33019748*l_1*t_1*l_2*p_2830080*l_2*p_11665459*l_1*t_1*l_1*p_1858560*l_1*p_23328*l_2*p_3382080*t_1*l_1*p_1499520*l_2*p_10695456*t_1*l_1*p_1499520*l_1*p_42665472*t_1*l_2*p_1985280*l_1*p_53350595*t_2*l_2*p_2956800*l_1*p_85377137>")
]

Failed = False
for gen in gens_comp_5:             #Verify that no generators commute mod x5 with generators of 5^{1+6}
    if commutator(gen, gens_normal_5[0]).as_int() in cyclic_as_int:
        Failed = False

print("Checking generators are in the complement 2 * J_2:4:", not Failed)



print("The generators commute:", gens_comp_5[0] * gens_comp_5[1] == gens_comp_5[1] * gens_comp_5[0])

G_comp = group_generated_by(gens_comp_5, n=5**2)



# Output:

# Checking the normal subgroup 5^{1+6}

# Verifying that all generators commute mod x5: True

# Limit 625 ; have 625 in time  33.6535
# Limit 625 ; have 625 in time  35.4755
# The intersection is trivial modulo x5: True

# Checking the complement 2 * J_2:4

# Checking generators are in the complement 2 * J_2:4: True
# The generators commute: True
# Limit 25 ; have 25 in time  0.9813