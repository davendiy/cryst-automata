
from sage.all import table, latex

from src.space_groups import prepare_gap_env
from src.srdegrees import SR_Degrees

prepare_gap_env()


headers = ('ITA num', r'$A \in \mathcal{X}$', r'$\det{A}$', r'$A \in \mathcal{Y}$')
table_end = r'\end{tabular}'
cur_t = [headers]
first_i = 1
for i in range(1, 18):
    t = SR_Degrees(i, dim=2, verbose=0)
    G = t.G

    for A in G.point_group_normalizer():
        # TODO: add check whether A forms virtual endomorphism
        compat = '+' if t.cocycle_compat(A) else ''
        cur_t.append((i, A, A.det().simplify_rational(), compat))

    if i in [10, 12, 15, 17]:
        cur_t = table(cur_t, header_row=True, frame=True)
        if first_i == i:
            caption = r'\caption{' + f'Point group compatibility solutions for group {i}' + r'}'
        else:
            caption = r'\caption{' + f'Point group compatibility solutions for groups {first_i}-{i}' + r'}'
        print(r'\begin{table}[H]')
        print(r'\begin{center}')
        print(latex(cur_t))
        print(r'\end{center}')
        print(caption)
        print(r'\end{table}')
        print()
        cur_t = [headers]
        first_i = i+1
