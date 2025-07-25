
from sage.all import matrix, QQ

import random

from src.planar_groups import (
    prepare_gap_env, to_L_basis,
    build_finite_group, from_indices_list, 
    SpaceGroup_Element, SpaceGroup_gap
    )


def test_build_finite():
    prepare_gap_env()

    for n in range(1, 17):
        G = to_L_basis(n, dim=2)
        P = G.PointGroup()
        gens = [matrix(QQ, el) for el in P.GeneratorsOfGroup()]

        triv = matrix(QQ, [[1, 0], [0, 1]])

        found = str(sorted(build_finite_group(gens, triv)))

        needed = str(sorted([str(matrix(QQ, el)) for el in P.AsList()]))

        assert found == needed 

        found = build_finite_group(gens, triv)
        for el, seq in found.items(): 

            res = from_indices_list(gens, triv, seq)
            assert str(el) == str(res), f"{el}, {res}, {seq}"


def test_space_group(): 
    G = SpaceGroup_gap.from_gap_cryst(12, dim=2)

    for el in G.G_gens: 
        assert el in G

    for x in G.snot: 
        assert x in G 

    for x in G.snot: 
        for y in G.snot: 
            assert (x * y) in G


    test1 = SpaceGroup_Element(
        matrix(QQ, [[-1, 0, 1/2], 
                    [0, 1, 0], 
                    [0, 0, 1]])
    )

    assert test1 not in G

    test2 = SpaceGroup_Element(
            matrix(QQ, [[1, 0, 1/2], 
                        [0, 1, 0], 
                        [0, 0, 1]])
        )
    assert test2 not in G


def test_word_space_group(): 
    G = SpaceGroup_gap.from_gap_cryst(12, dim=2)

    s_n = len(G.snot)

    for _ in range(100): 
        idx = random.randint(0, s_n-1)
        x = SpaceGroup_Element(G.snot[idx])

        k = random.randint(-100, 100)
        j = random.randint(-100, 100)

        x._body[0, 2] += k
        x._body[1, 2] += j

        word = ''.join(G.as_word(x, readable=True))
        if k:
            assert (_t := f'e_1^({k})') in word, f"expected: {_t}, got: {word}"
        if j: 
            assert (_t := f'e_2^({j})') in word, f"expected: {_t}, got: {word}"

    for _ in range(1000): 
        x = G.random_element()
        it_seq = G.as_word(x, readable=False)
        res = from_indices_list(G.G_sorted_gens, SpaceGroup_Element(G.G_triv), it_seq)

        assert str(res) == str(x), f"expected: \n{x}, got: \n{res}, {it_seq}"


if __name__ == "__main__": 
    test_build_finite()
    test_space_group()
    test_word_space_group()
    print('all good.✅')
