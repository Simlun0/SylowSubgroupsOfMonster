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


g7 = MM("M<y_5d3h*x_0a6dh*d_8d4h*p_111142481*l_1*p_2999040*l_1*p_43234193>")


gens_normal_7 = [
MM("M<y_718h*x_13deh*d_249h*p_224166850*l_1*p_2027520*l_1*p_488181*t_1*l_2*p_2956800*l_1*p_32991828*l_2*t_2*l_2*p_2386560*l_2*p_1914710*l_2*t_1*l_2*p_1985280*l_1*p_22368501*t_1*l_1*p_1394880*l_1*p_3168*l_2*p_3867840*t_2*l_2*p_1900800*l_2*p_2880*t_1*l_2*p_1985280*l_1*p_85331925>"),
MM("M<y_457h*x_1be7h*d_0b13h*p_213988192*l_2*p_1943040*l_2*p_21800869*l_1*t_1*l_2*p_1943040*l_2*p_1907843*l_1*t_2*l_1*p_2640000*l_1*p_22757251*l_1*t_2*l_2*p_2344320*l_2*p_10562*t_1*l_2*p_1457280*l_1*p_486185*t_2*l_2*p_3840*l_2*p_10666800*l_1*p_6086400*t_1*l_2*p_2597760*l_1*p_85329029*t_2*l_2*p_1457280*l_1*p_17322>"),
MM("M<y_1cah*x_3h*d_1f9h*p_183743694*l_1*p_1457280*l_2*p_23219975*l_2*t_2*l_1*p_2027520*l_1*p_11684675*l_2*t_1*l_2*p_1985280*l_1*p_39536*l_1*t_2*l_2*p_2386560*l_2*p_53403281*t_2*l_1*p_2027520*l_1*p_12103029*l_2*t_2*l_1*p_1499520*l_1*p_42837731*t_2*l_2*p_2386560*l_2*p_85818768>"),
MM("M<y_91h*x_2c5h*d_0e8h*p_230757800*l_1*p_2999040*l_1*p_10709008*t_1*l_2*p_2787840*l_2*p_11595172*l_2*t_2*l_1*p_1920*l_1*p_4032*l_1*p_3360000*t_1*l_2*p_2830080*l_2*p_85415619*t_1*l_1*p_1499520*l_1*p_64015156*t_2*l_2*p_2956800*l_1*p_21445856*t_2*l_2*p_1985280*l_1*p_21445856>")
]

cyclic = []
cyclic_as_int = []

curr = g7
for _ in range(7):
    cyclic.append(curr)
    cyclic_as_int.append(curr.as_int())
    curr *= g7


def commutator(a,b):
    return a**-1 * b**-1 * a * b 


print("Checking the normal subgroup 7^{1+4}")
print()



Failed = False
for i in range(len(gens_normal_7)):         #Verify that all generators commute mod x7
    for j in range(i+1, len(gens_normal_7)):
        a = gens_normal_7[i]
        b = gens_normal_7[j]
        if commutator(a,b).as_int() not in cyclic_as_int:
            Failed = True

print("Verifying that all generators commute mod x7:", not Failed)
print()

gens_abc = gens_normal_7[:2]
gens_xyz = gens_normal_7[2:]


G_abc = group_generated_by(gens_abc, n=7**3)            #Generate partial group s
G_xyz = group_generated_by(gens_xyz, n=7**3)

G_abc_as_int = [g.as_int() for g in G_abc]              #int comparisons are much faster than mmgroup object comparisons
G_xyz_as_int = [g.as_int() for g in G_xyz]


Failed = False
for g_abc_as_int in G_abc_as_int:
    if g_abc_as_int in G_xyz_as_int and g_abc_as_int not in cyclic_as_int:
        Failed = True


print("The intersection is trivial modulo x7:", not Failed)
print()


print("Checking the complement 3 x 2S_7")
print()

gen_comp_7 = MM("M<y_184h*x_16b6h*d_50dh*p_203787511*l_2*p_1457280*l_1*p_32537859*l_2*t_1*l_1*p_1499520*l_1*p_32067200*l_1*t_1*l_2*p_1943040*l_2*p_32998679*l_1*t_2*l_2*p_2787840*l_2*p_10580*t_1*l_1*p_7941120*l_2*t_2*l_1*p_1499520*l_2*p_13059382*l_2*t_2*l_2*p_2344320*l_2*p_6725*t_2*l_2*p_1985280*l_1*p_64042011>")

Failed = False
if commutator(gen_comp_7, gens_normal_7[0]).as_int() in cyclic_as_int:    #Verify that the generator does not commute mod x7 with generators of 7^{1+4}
    Failed = False

print("Checking generator is in the complement 3 x 2S_7:", not Failed)

print("Order of generator:", gen_comp_7.order())


# Output:

# Checking the normal subgroup 7^{1+4}

# Verifying that all generators commute mod x7: True

# Limit 343 ; have 343 in time  15.0088
# Limit 343 ; have 343 in time  14.4152
# The intersection is trivial modulo x7: True

# Checking the complement 3 x 2S_7

# Checking generator is in the complement 3 x 2S_7: True
# Order of generator: 7