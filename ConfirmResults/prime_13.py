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

def commutator(a,b):
    return a**-1 * b**-1 * a * b 


x13 = MM("M<y_381h*x_9c3h*d_6b6h*p_161269430*l_2*p_1457280*l_1*p_96463266*l_2*p_528000>")



cyclic = []
cyclic_as_int = []

curr = x13
for _ in range(13):
    cyclic.append(curr)
    cyclic_as_int.append(curr.as_int())
    curr *= x13

gens_13 = [
    MM("M<y_40ch*x_11ech*d_7aah*p_56137889*l_1*p_2027520*l_1*p_32567525*l_1*t_1*l_1*p_2640000*l_1*p_1482712*l_1*t_2*l_1*p_1393920*l_2*p_1200*l_1*p_10275840*t_1*l_1*p_2640000*l_1*p_1017971*t_2*l_1*p_2999040*l_1*p_85835155*l_2*p_507840*t_2*l_1*p_50160000*l_1*p_1795200*l_1>"),
    MM("M<y_22fh*x_0eaeh*d_0bdh*p_163280938*l_1*p_23040*l_2*p_23193120*t_1*l_2*p_1457280*l_1*p_32533077*l_2*t_2*l_1*p_1920*l_2*p_465936*l_2*p_965760*t_1*l_2*p_1943040*l_2*p_21388998*t_2*l_1*p_1920*l_2*p_10665936*l_2*p_2398080*t_2*l_1*p_960*l_2*p_42726192*t_1*l_2*p_2386560*l_2*p_21429446>"),
]


Failed = False
for i in range(len(gens_13)):         #Verify that all generators commute mod x13
    for j in range(i+1, len(gens_13)):
        a = gens_13[i]
        b = gens_13[j]
        if commutator(a,b).as_int() not in cyclic_as_int:
            Failed = True

print("Verifying that all generators commute mod x13:", not Failed)


G = group_generated_by(gens_13, n=13**3)            #While we could use the same method as the code for prime_5 and prime_7
                                                    # it is easier to simply generate the whole group in this case.

print("The order of the generated group is correct:", len(G) == 13**3)


# Output:

# Verifying that all generators commute mod x13: True
# Limit 2197 ; have 2197 in time  92.9706
# The order of the generated group is correct: True