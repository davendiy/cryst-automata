
from sage.all import matrix, QQ

from src.space_groups import SpaceGroup_gap, prepare_gap_env
from src.srdegrees import SR_Degrees



def test_self_similar(): 

    G = SpaceGroup_gap.from_gap_cryst(13, dim=2)

    T = matrix(QQ, [
        [-1/3, 2/3, 0], 
        [1/3, 1/3, 0], 
        [0, 0, 1]
    ])

    act = G.self_similar(T, verbose=True)
    G._wreath_recursion(act)


def test_wreath_recursion(): 
    G = SpaceGroup_gap.from_gap_cryst(7, dim=2)

    T = matrix(QQ, [
        [1/3, 0, 0], 
        [0, 1/2, 0], 
        [0, 0, 1],
    ])
    act = G.self_similar(T, verbose=True)
    G._wreath_recursion(act)


def test_srdegrees(): 
    x = SR_Degrees(7)
    x.algorithm()
    

if __name__ == '__main__': 
    prepare_gap_env()
    test_self_similar() 
    # test_wreath_recursion()
    test_srdegrees()
    print('all good.✅')
