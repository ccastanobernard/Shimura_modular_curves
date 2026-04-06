# Mon Apr  6 14:35:05 CST 2026

# The following Sage program computes the prime discriminant factorisation of a fundamental discriminant.

def decompose_fundamental_discriminant(D):
    if not D.is_fundamental_discriminant():
        raise ValueError(f"{D} is not a fundamental discriminant.")
    factors = abs(D).factor()
    prime_factors = [f[0] for f in factors]
    prime_discs = []
    if 2 in prime_factors:
        if D % 16 == 0:
            if (D // -4) % 8 == 0:
                pass
        if D % 8 == 1 or D % 8 == 5:
            pass
        elif D % 4 == 0:
            if D % 8 == 0:
                if (D // 8) % 4 == 1: prime_discs.append(8)
                else: prime_discs.append(-8)
            else:
                prime_discs.append(-4)
        prime_factors = [p for p in prime_factors if p != 2]
    for p in prime_factors:
        if p % 4 == 1:
            prime_discs.append(p)
        else:
            prime_discs.append(-p)
    return prime_discs

# The output is [N, group_order, i] [p_1^*, ... ,p_n^*], where:
# N is the level of the Shimura modular curve asociated with a real quadratic field of discriminant N with n distinct prime divisors, for N between minimum and maximum;
# group_order, i is the group ID of the Atkin-Lehner group of the Shimura modular curve;
# p_1^*, ... ,p_n^* is the prime discriminant decomposition of N.

n = 2
minimum = 5
maximum = 3997
fundamental = [D for D in range(minimum, maximum)  if Integer(D).is_fundamental_discriminant() and gp.omega(D) == n]
for N in fundamental:
    delta = [n for n in range(1,N) if kronecker(N, n) == 1]
    congruence_subgroup = GammaH(N, delta)
    C = congruence_subgroup.cusps()
    generators = []
    for p in prime_divisors(N):
        W = matrix(congruence_subgroup.atkin_lehner_matrix(p))
        P = [congruence_subgroup.reduce_cusp(alpha.apply(W.list())) for alpha in C]
        perm = Word(P).standard_permutation() / Word(C).standard_permutation()
        assert [C[i-1] for i in perm] == P
        generators.append(perm)
    atkin_lehner_group = PermutationGroup(generators)
    group_order = atkin_lehner_group.order()
    i = 1
    G = gap.Image(gap.IsomorphismPermGroup(gap.SmallGroup(group_order, i)))
    group = PermutationGroup(gap_group = G)
    while not group.is_isomorphic(atkin_lehner_group):
        i = i + 1
        G =   gap.Image(gap.IsomorphismPermGroup(gap.SmallGroup(group_order, i)))
        group = PermutationGroup(gap_group = G)
    print([N, group_order, i], decompose_fundamental_discriminant(Integer(N)))


