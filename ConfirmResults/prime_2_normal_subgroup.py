from mmgroup import *
import time


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

z = MM("M<x_1000h>")            #z (that G_x0 is defined as the centralizer of) is the central element of G_x0.
                                # we use is as we did x3, x5 etc.

gens_normal_2 = [
    MM("M<d_1h>"),
    MM("M<d_2h>"),
    MM("M<d_4h>"),
    MM("M<d_8h>"),
    MM("M<d_10h>"),
    MM("M<d_20h>"),
    MM("M<d_40h>"),
    MM("M<d_80h>"),
    MM("M<d_100h>"),
    MM("M<d_200h>"),
    MM("M<d_400h>"),
    MM("M<d_800h>"),
    MM("M<x_1h*d_0eh>"),
    MM("M<x_2h*d_0dh>"),
    MM("M<x_4h*d_0bh>"),
    MM("M<x_8h*d_7h>"),
    MM("M<x_10h*d_0fh>"),
    MM("M<x_20h*d_0fh>"),
    MM("M<x_40h*d_0eh>"),
    MM("M<x_80h*d_0eh>"),
    MM("M<x_100h*d_401h>"),
    MM("M<x_200h*d_401h>"),
    MM("M<x_400h>"),
    MM("M<x_800h>"),
]

cyclic = []
cyclic_as_int = []

curr = z
for _ in range(2):
    cyclic.append(curr)
    cyclic_as_int.append(curr.as_int())
    curr *= z


def commutator(a,b):
    return a**-1 * b**-1 * a * b 


print("Checking the normal subgroup 2^{1+24}")
print()



Failed = False
for i in range(len(gens_normal_2)):         #Verify that all generators commute mod x5
    for j in range(i+1, len(gens_normal_2)):
        a = gens_normal_2[i]
        b = gens_normal_2[j]
        if commutator(a,b).as_int() not in cyclic_as_int:
            Failed = True

print("Verifying that all generators commute mod z:", not Failed)
print()

gens_abc = gens_normal_2[:12] + [z]         #z is added because it cannot be generated using the generators. 
gens_xyz = gens_normal_2[12:] + [z]         # in the other primes the central elements could be generated using the generators


G_abc = group_generated_by(gens_abc, n=2**13)            #Generate partial groups. The groups are huge however we are in G_x0.
G_xyz = group_generated_by(gens_xyz, n=2**13)            # Computations are much faster.

G_abc_as_int = [g.as_int() for g in G_abc]              #int comparisons are much faster than mmgroup object comparisons
G_xyz_as_int = [g.as_int() for g in G_xyz]


Failed = False
for g_abc_as_int in G_abc_as_int:
    if g_abc_as_int in G_xyz_as_int and g_abc_as_int not in cyclic_as_int:
        Failed = True


print("The intersection is trivial modulo z:", not Failed)
print()


# Output:

# Checking the normal subgroup 2^{1+24}

# Verifying that all generators commute mod z: True

# Limit 8192 ; have 8192 in time  1.2322
# Limit 8192 ; have 8192 in time  1.2418
# The intersection is trivial modulo z: True