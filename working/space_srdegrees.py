import os

from src.space_groups import prepare_gap_env
from src.srdegrees import solve_simple_mat, factorize
from src.cryst3.all_norm3 import all_unique_matrices


prepare_gap_env()

folder = 'cached_normalizers_dim3/'
cached = {}
for i in range(230):
    name = f'{folder}/norm_{i}.txt'
    if os.path.exists(name):
        with open(name) as file:
            matrices = file.read().split('\n\n')
            matrices = filter(bool, [m.strip() for m in matrices])
            cached.update({m: i for m in matrices})

count = 0
for A in all_unique_matrices:
    # print(A)
    # print()
    try:
        solve_simple_mat(A)
    except Exception:
        print(A)
        print()
        sub, v = factorize(A)
        print('sub determinant:', sub.det())
        print(f'group num: {cached[str(A)]}')
        count += 1

        print('---------------------------')

print(f'no solvable: {count}')
