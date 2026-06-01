from sage.all import latex

from src.space_groups import prepare_gap_env
from src.srdegrees import solve_simple_mat, factorize
from src.cryst3.gen_norms3 import all_unique_matrices, find_groups


prepare_gap_env()

count = 0
for A in all_unique_matrices:
    # print(A)
    # print()
    try:
        solve_simple_mat(A)
    except Exception:
        print("$$", latex(A), "$$")
        print()
        sub, v = factorize(A)
        print('sub determinant: $', latex(sub.det()), '$')
        print(f'group num: {find_groups(A)}')
        count += 1
        print()
        print()

print(f'no solvable: {count}')
