
from sage.all import table, latex

from src.space_groups import prepare_gap_env  # , check_div
# from src.cryst3.common import cached_simple_matrices
# from src.cryst3.all_norm3 import all_unique_matrices
from src.cryst3.gen_norms3 import find_groups, all_unique_lattice_matrices

# from src.srdegrees import tau

prepare_gap_env()


# for gr_num in cached_simple_matrices:
#     print(f'################# {gr_num} ######################')
#     for a in cached_simple_matrices[gr_num]:
#         # if tau(a) == 0:
#         #     continue
#         # print(a)
#         print("equation:", (a.det() - a.trace()).simplify_rational())
#         # print('tau(A):', tau(a))


# check either matrix is factorizable or is simple
# for a in all_unique_matrices:
#     if check_div(a):
#         continue
#     if tau(a) != 0 or a.trace() != 0:
#         print(a)
#         print(f'tau: {tau(a)}')
#         print(f'trace: {a.trace()}')
#         print()


# for el in all_unique_matrices:
#     print(find_groups(el))

N = 2
M = 17

header = []
for i in range(N):
    header.append('No.')
    header.append('Matrix')
    header.append('Groups No.')

matrices = [(find_groups(A, latt=True, all=False), A) for A in all_unique_lattice_matrices]
matrices = sorted(matrices)


def normalize_grouplist(groups):
    if len(groups) <= 4:
        return table([[','.join([str(el) for el in groups])]])
    elif len(groups) <= 16:
        return table(rows=[[','.join(str(el) for el in groups[i:i+4])] for i in range(0, len(groups), 4)])
    else:
        return f'{groups[0]},{groups[1]}...{groups[-2]},{groups[-3]}'


matrices = [(i, A, normalize_grouplist(groups)) for i, (groups, A) in enumerate(matrices)]


def unwrap(row):

    res = []
    for subrow in row:
        for el in subrow:
            res.append(el)
    return res


byNrow = []
for i in range(0, len(matrices), N):
    byNrow.append(unwrap(matrices[i:i+N]))

for i in range(0, len(byNrow), M):
    print(r'\begin{table}[H]')
    print(r'\begin{center}')
    print(r'\begin{scriptsize}')
    if i == 0:
        header[2] = r'Groups No.\footnotemark[1]{}'
    else:
        header[2] = r'Groups No.'
    print(latex(table([header] + byNrow[i: i+M], frame=True, header_row=True)))
    print(r'\end{scriptsize}')

    print(r'\caption{' + f'Unique matrices {i*N}..{(i+M)*N} satisfying PC condition for space groups' + r'}')
    print(r'\label{' + f'tab:norm_dim3_{i*N}' + r'}')
    print(r'\end{center}')
    print(r'\end{table}')
    print()
    if i == 0:
        print(r'\footnotetext[1]{Only symmorphic groups are provided as point group representatives.}')
    print(r'\newpage')
    print()
