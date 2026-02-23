
from ..space_groups import SpaceGroup_gap, check_div
from ..srdegrees import tau

from .common import point_groups_representatives


def nprint(*args, **kwargs):
    print(*args, **kwargs)
    # print(*args, **kwargs, file=__f)


def bruteforce_normalizers():
    simple_found = set()

    for i in point_groups_representatives:
        nprint(f'##################### {i} ############################')
        G = SpaceGroup_gap.from_gap_cryst(i, dim=3, change_basis=True)

        norm_matrices = set()
        for A in G.point_group_normalizer():
            norm_matrices.add(str(A))
            if check_div(A):
                continue
            simple_found.add(str(A))
            nprint(A)
            nprint('tau(A) = ', tau(A).simplify_rational())
            nprint('\n')
        with open(f'norm_{i}.txt', 'w') as file:
            for el in norm_matrices:
                file.write(el)
                file.write('\n\n')

    return simple_found
