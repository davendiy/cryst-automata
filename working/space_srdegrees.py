
from src.space_groups import prepare_gap_env
from src.srdegrees import solve_simple_mat, factorize
from src.cryst3.all_norm3 import all_unique_matrices


prepare_gap_env()

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
        count += 1

        print('---------------------------')

print(f'no solvable: {count}')
