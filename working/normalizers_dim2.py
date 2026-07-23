
from sage.all import table, latex, var, matrix

from src.space_groups import prepare_gap_env
from src.srdegrees import SR_Degrees

prepare_gap_env()


def unwrap(row):

    res = []
    for subrow in row:
        for el in subrow:
            res.append(el)
    return res


headers = (r'Group No.', r'$A \in \mathcal{Y}$', r'$\det{A}$')
table_end = r'\end{tabular}'
cur_t = []
first_i = 1
for i in range(1, 18):
    t = SR_Degrees(i, dim=2, verbose=0)
    G = t.G

    if i in [1, 2]:
        y0, y1, y2, y3 = var('y0 y1 y2 y3')
        C = matrix([[y0, y1], [y2, y3]])
        cur_t.append((str(i), C, C.det()))

    for A in G.point_group_normalizer():
        # TODO: add check whether A forms virtual endomorphism
        # compat = '+' if t.cocycle_compat(A) else ''

        l_compat = t.lattice_compat(A)
        if l_compat.failed:
            continue

        # print(l_compat.result)
        c_compat = t.cocycle_compat_v2(l_compat.result)
        if c_compat.failed:
            continue

        C = c_compat.result
        cur_t.append((str(i), C, C.det()))

    if i in [10, 15, 17]:

        double_cur = [headers + headers]
        if len(cur_t) % 2:
            cur_t.append(['' for _ in headers])
        half = len(cur_t) // 2
        for j in range(0,  half):
            double_cur.append(unwrap([cur_t[j], cur_t[half + j]]))
        cur_t = table(double_cur, header_row=True, frame=True)
        if first_i == i:
            caption = r'\caption{' + f'PC-CC solutions for plane group {i}' + r'}'
        else:
            caption = r'\caption{' + f'PC-CC solutions for plane groups {first_i}-{i}' + r'}'
        print(r'\begin{table}[H]')
        print(r'\begin{scriptsize}')
        print(r'\begin{center}')
        print(latex(cur_t))
        print(r'\end{center}')
        print(caption)
        print(r'\label{' + f'tab:norm_dim2_{i}' + r'}')
        print(r'\end{scriptsize}')
        print(r'\end{table}')
        print()
        cur_t = []
        first_i = i+1
