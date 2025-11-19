
from sage.all import MatrixGroup
from src.space_groups import SpaceGroup_gap, prepare_gap_env, check_div
# from src.srdegrees import _factorize

prepare_gap_env()

__f = open('normalizers_dim3_changed_basis', 'w')


def nprint(*args, **kwargs):
    print(*args, **kwargs)
    print(*args, **kwargs, file=__f)


found = set()
point_groups = set()

for i in range(2, 231):
    nprint(f'##################### {i} ############################')
    G = SpaceGroup_gap.from_gap_cryst(i, dim=3, change_basis=True)
    els = list(sorted(MatrixGroup(G.P_gens)))
    if str(els) in point_groups:
        continue
    point_groups.add(str(els))
    for A in G.point_group_normalizer():
        if str(A) in found:
            continue
        if check_div(A):
            continue
        found.add(str(A))
        nprint(A)
        nprint('\n')

__f.close()
