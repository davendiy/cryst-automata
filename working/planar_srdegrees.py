from sage.all import latex, table, matrix, ascii_art

from src.space_groups import SpaceGroup_gap, prepare_gap_env
from src.srdegrees import SR_Degrees


prepare_gap_env()

for n in range(2, 18):
    # for n in [12]:
    # G = SpaceGroup_gap.from_gap_cryst(n, dim=2)

    sr = SR_Degrees(n, method='latex', verbose=1)
    if n == 2:
        sr.texdoc_header()

    eig_table, res_table, cond_table = sr.sr_degrees()
    eig_table = table(eig_table, header_row=True)

    for el in res_table[1:]:
        if el[2]:
            el[2] = table(el[2])

    res_table = table(res_table, header_row=True, frame=True)
    cond_table = table(cond_table, header_row=True)
    print(latex(eig_table))
    print(r'\\')

    print(latex(res_table))
    print(r'\\')
    print(latex(cond_table))
    if n == 17:
        sr.texdoc_ending()
