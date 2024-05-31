from mmgroup import *


a = MM("M<d_200h>")
b = MM("M<y_756h*x_1b0eh*d_29ah*p_142818544*l_2*p_1415040*l_1*p_12992112*l_1*t_1*l_2*p_2787840*l_2*p_31997178*l_1*t_1*l_2*p_960*l_2*p_10666896*l_1*p_4292160*t_1*l_1*p_1499520*l_1*p_21378434*t_1*l_2*p_1858560*l_2*p_464880*l_1*p_1927680*t_2*l_2*p_1393920*l_1*p_85409856*t_1*l_2*p_2956800*l_1*p_85837058>")

SLPs = [
    "ababababababbababababbabbabbabb",      #17
    "ababababababbabababababbabb",          #19
    "ababababababbabababbababbabbabb",      #23
    "ab",                                   #29
    "ababababbabababbababbabb",             #31
    "abababababbababbabbababbabb",          #41
    "ababababababababbababbabbabb",         #47
    "abababababababbabababbababb",          #59
    "abababababababbabbabababbababbabb"     #71
]

primes = [
    17,
    19,
    23,
    29,
    31,
    41,
    47,
    59,
    71
]

def word_to_mmgroup(word):
    ret = MM("M<1>")        #Unit element
    for letter in word:
        if letter == 'a':
            ret *= a
        if letter == 'b':
            ret *= b

    return ret

for index, word in enumerate(SLPs):
    g = word_to_mmgroup(word)
    print("Order of word:", g.order(), "Expected:", primes[index], "Equal:", g.order()==primes[index])



# Output:

# Order of word: 17 Expected: 17 Equal: True
# Order of word: 19 Expected: 19 Equal: True
# Order of word: 23 Expected: 23 Equal: True
# Order of word: 29 Expected: 29 Equal: True
# Order of word: 31 Expected: 31 Equal: True
# Order of word: 41 Expected: 41 Equal: True
# Order of word: 47 Expected: 47 Equal: True
# Order of word: 59 Expected: 59 Equal: True
# Order of word: 71 Expected: 71 Equal: True