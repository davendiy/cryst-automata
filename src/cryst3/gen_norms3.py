
import os

import re

from base64 import b64encode, b64decode
from sage.all import dumps, loads

from ..space_groups import SpaceGroup_gap, check_div
from ..srdegrees import tau

from .common import point_groups_representatives, by_point_groups

root = os.path.realpath(__file__)
root = os.path.realpath(os.path.join(root, '../../../'))
datafolder = os.path.join(root, 'cached_normalizers_dim3')


cached_b64 = re.compile(r'######## INTERNAL_REPR #######\n(.*)\n#### END OF INTERNAL_REPR ####')


def nprint(*args, **kwargs):
    print(*args, **kwargs)
    # print(*args, **kwargs, file=__f)


def bruteforce_normalizers(folder=datafolder, enforce_integral=True):
    simple_found = set()

    # for i in [148]:
    for i in point_groups_representatives:
        nprint(f'##################### {i} ############################')
        G = SpaceGroup_gap.from_gap_cryst(i, dim=3, change_basis=True)

        norm_matrices = {}
        for A in G.point_group_normalizer(enforce_integral=enforce_integral):
            norm_matrices[str(A)] = A
            if check_div(A):
                continue
            simple_found.add(str(A))
            nprint(A)
            nprint('tau(A) = ', tau(A).simplify_rational())
            nprint('\n')
        with open(os.path.join(folder, f'norm_{i}.txt'), 'w') as file:
            for st_A, A in norm_matrices.items():
                file.write(st_A)
                file.write('\n\n')
                file.write('\n######## INTERNAL_REPR #######\n')
                file.write(b64encode(dumps(A)).decode('ascii'))
                file.write('\n#### END OF INTERNAL_REPR ####\n')

    return simple_found


__cached = {}
__cached_matrices = {}
all_matrices = []


for i in point_groups_representatives:
    filename = os.path.join(datafolder, f'norm_{i}.txt')
    if not os.path.exists(filename):
        print(f'couldnt find {filename}')
        continue
    with open(filename) as file:
        data = file.read()
        __cached[i] = data
        mts = [el.encode() for el in cached_b64.findall(data)]
        mts = [loads(b64decode(el)) for el in mts]
        __cached_matrices[i] = mts
        all_matrices.extend(mts)


all_unique_matrices = list({str(el): el for el in all_matrices}.values())


def load_cached_normalizers(group_num):
    point_group = [el[0] for el in by_point_groups if group_num in el]
    if len(point_group) != 1:
        raise ValueError(f'found {len(point_group)} point group representatives for {group_num}..')
    point_group = point_group[0]
    matrices = __cached_matrices[point_group]
    return matrices


def find_groups(A):
    res = []
    for i, val in __cached.items():
        if str(A) in val:
            res.append(i)
    return res
