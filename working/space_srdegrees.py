import sys

from tqdm import tqdm
from sage.all import latex, table
from src.space_groups import prepare_gap_env

from src.srdegrees import SR_Degrees, solve_simple_mat, factorize, tau, standardize_sol_matrix
from src.cryst3.gen_norms3 import load_cached_normalizers


prepare_gap_env()

VER = False

filename = sys.argv[1]


def get_polys(M):
    if tau(M) == 0 and M.trace() == 0:
        return M.det().simplify_rational()

    t = factorize(M)
    if t is None:
        raise ValueError('something went wrong')

    subM, p = t
    q = subM.det().simplify_rational()
    g = subM.trace().simplify_rational()
    return p.simplify_rational(), q, g


def printv(*args, **kwargs):
    if VER:
        print(*args, **kwargs)


headers = ['group num', r'$p_\sigma$', r'$q_\sigma$', r'$g_\sigma$']
headers_simple = ['group num', r'$f_\sigma$']
res_table_simple = []
res_table = []

# for i in [24]:
# for i in [9]:
# for i in [148]:
for i in tqdm(range(2, 231)):
    printv('\n\n')
    printv(f'==================== group #{i} ======================')
    t = SR_Degrees(group_index=i, dim=3, verbose=(2 if VER else 0))
    # for A in t.G.point_group_normalizer(enforce_integral=False):
    for A in load_cached_normalizers(i):

        printv()
        printv('original A:')
        printv(A)

        printv('\nlattice compat A:')
        B = t.lattice_compat(A)
        assert not B.failed
        B = B.result
        printv(B)

        if t.G.is_symmorphic():
            printv(solve_simple_mat(B).result)
            B = standardize_sol_matrix(B, 'y')
            if tau(B) == 0 and B.trace() == 0:
                res_table_simple.append([i] + [get_polys(B)])
            else:
                res_table.append([i] + list(get_polys(B)))
            continue

        # t.construct_congruences_v2(A.inverse().simplify_rational(), A)

        m = t.cocycle_compat_v2(B)
        if m.failed:
            printv('no solutions')
            continue
        printv('cocycle compat A:')
        printv(m.result)
        printv()
        printv('det(A) =', m.result.det().simplify_rational())
        printv('tr(A) =', m.result.trace().simplify_rational())
        printv(solve_simple_mat(m.result))

        if tau(m.result) == 0 and m.result.trace() == 0:
            res_table_simple.append([i] + [get_polys(m.result)])
        else:
            res_table.append([i] + list(get_polys(m.result)))

# print(res_table)

M = 50

for i in range(len(res_table)):
    for j in range(1, 4):
        res_table[i][j] = '$' + latex(res_table[i][j]) + '$' if res_table[i][j] is not None else ''

with open(filename, 'w') as file:
    for i in range(0, len(res_table), M):
        print(r'\begin{table}[H]', file=file)
        print(r'\begin{center}', file=file)
        print(r'\begin{scriptsize}', file=file)
        chunk = res_table[i: i+M]
        print(latex(table([headers] + chunk, frame=True, header_row=True)), file=file)
        print(r'\end{scriptsize}', file=file)

        print(r'\caption{' + f'Case 1 for Groups {chunk[0][0]}..{chunk[-1][0]}' + r'}', file=file)
        print(r'\end{center}', file=file)
        print(r'\end{table}', file=file)
        print('', file=file)
        print(r'\newpage', file=file)
        print('', file=file)


for i in range(len(res_table_simple)):
    res_table_simple[i][1] = '$' + latex(res_table_simple[i][1]) + '$' if res_table_simple[i][1] is not None else ''

with open(filename, 'a') as file:
    for i in range(0, len(res_table_simple), M):
        print(r'\begin{table}[H]', file=file)
        print(r'\begin{center}', file=file)
        print(r'\begin{scriptsize}', file=file)
        chunk = res_table_simple[i: i+M]
        print(latex(table([headers_simple] + chunk, frame=True, header_row=True)), file=file)
        print(r'\end{scriptsize}', file=file)

        print(r'\caption{' + f'Case 2 for Groups {chunk[0][0]}..{chunk[-1][0]}' + r'}', file=file)
        print(r'\end{center}', file=file)
        print(r'\end{table}', file=file)
        print('', file=file)
        print(r'\newpage', file=file)
        print('', file=file)
