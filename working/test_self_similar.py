
from sage.all import matrix, QQ

from src.space_groups import SpaceGroup_gap, prepare_gap_env, is_simple
from src.srdegrees import SR_Degrees



def test_self_similar(): 

    G = SpaceGroup_gap.from_gap_cryst(13, dim=2)

    T = matrix(QQ, [
        [-1/3, 2/3, 0], 
        [1/3, 1/3, 0], 
        [0, 0, 1]
    ])

    act = G.self_similar(T, verbose=False)
    G._wreath_recursion(act)


def test_wreath_recursion(): 
    G = SpaceGroup_gap.from_gap_cryst(7, dim=2)

    T = matrix(QQ, [
        [1/3, 0, 0], 
        [0, 1/2, 0], 
        [0, 0, 1],
    ])
    act = G.self_similar(T, verbose=False)
    G._wreath_recursion(act)


def test_min_sr(): 
    for i in range(1, 18): 
        G = SpaceGroup_gap.from_gap_cryst(i, dim=2)
        G.minimal_self_replicating(verbose=False, build_wreath=False)


def test_naive_sr(): 
    for n in range(1, 18): 
        G = SpaceGroup_gap.from_gap_cryst(n, dim=2)
        G.naive_self_replicating(verbose=False, build_wreath=False)

    for n in [13, 17, 21, 100, 134, 200]: 
        G = SpaceGroup_gap.from_gap_cryst(n, dim=3)
        G.naive_self_replicating(verbose=False, build_wreath=False)


def test_srdegrees(): 
    x = SR_Degrees(7, verbose=False)
    x.algorithm()

    x = SR_Degrees(7, verbose=True)
    x.algorithm()


def test_simplicity(): 

    assert is_simple(matrix(QQ, [
        [1/2, 0], 
        [0, 1/2]
    ]))

    assert not is_simple(matrix(QQ, [
        [1, 2, 3], 
        [0, 1, 0], 
        [0, 0, 2],
    ]))

    assert is_simple(matrix(QQ, [
        [1/2, 0, 0, 0], 
        [0, 1/2, 0, 0], 
        [0, 0, 1/2, 0],
        [0, 0, 0, 1/2]
    ]))

    assert not is_simple(matrix(QQ, [
        [0, 0, 0], 
        [0, 0, 2], 
        [1, 2, 3]
    ]))

    assert not is_simple(matrix(QQ, [
        [1/2, 0], 
        [0, 1]
    ]))


if __name__ == '__main__': 
    prepare_gap_env()
    test_simplicity()
    test_self_similar() 
    test_wreath_recursion()
    test_min_sr()
    test_naive_sr()
    test_srdegrees()
    print('all good.✅')
