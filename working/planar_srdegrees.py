from sage.all import latex, table

from src.space_groups import prepare_gap_env
from src.srdegrees import SR_Degrees

def uncover(el):
    if isinstance(el, list) or isinstance(el, tuple):
        if len(el) == 1:
            return el[0]

    return el


prepare_gap_env()

general_table = [["num", "srdegrees"]]

sr = SR_Degrees(1, method='latex', verbose=2)
sr.texdoc_header()

for n in range(2, 18):
    # for n in [12]:

    sr = SR_Degrees(n, method='latex', verbose=2)
    eig_table, res_table, cond_table, pred_col_table = sr.sr_degrees()
    eig_table = table(eig_table, header_row=True, frame=True)

    smaller_table = set()
    for el, pred_sols in zip(res_table[1:], pred_col_table[1:]):
        pred_sols = pred_sols[1]

        if (el_n := len(el[1].variables())) > 1:    # type: ignore
            degree = f"^{el_n}"
            lp = '('
            rp = ')'
        else:
            lp = rp = ''
            degree = ''

        field = r' \times '.join(r'\mathbb{Z}' if pred_sols.get(v,0) == 0 else r'OZ' for v in el[1].variables())  # type: ignore
        # (x_1, x_2) \in Z^2
        inz = lp + ",".join([latex(_x) for _x in el[1].variables()]) + rp + sr.in_sym + ' ' + field # type: ignore
        # {(x_1, x_2) \in Z^2}
        inz = r"\left\{" + latex(el[1]) + r'\, | \,' + inz + r'\right\}'

        if not el[2]:
            smaller_table.add("$" + inz + "$")
            continue

        # / {(-1, 1), (1, 1), (0, 0)}
        new_arr = r"\left\{"
        new_arr += ','.join(latex(uncover(tuple(sol))) for sol in el[2])
        new_arr += r'\right\}'
        el[2] = "$" + new_arr + "$"

        smaller_table.add(
                r"$" + inz.rstrip(r'\right\}') + r"\, / \," + new_arr + r'\right\}' + "$",
        )

    if smaller_table:
        general_table.append([n, table([[el] for el in smaller_table])])  # type: ignore

    res_table = table(res_table, header_row=True, frame=True)
    pred_col_table = table(pred_col_table, header_row=True, frame=True)

    print(latex(eig_table))
    print(r'\\')
    print(latex(res_table))
    print(r'\\')
    print()
    if sr.G.is_symmorphic():
        continue

    for row in cond_table[1:]:
        padded = []
        for i, el in enumerate(row[1]):
            if i % 4 == 0:
                padded.append(list())
            padded[-1].append("$" + latex(el) + r"\in \mathbb{Z}$")
        row[1] = table(padded)

    print(latex(pred_col_table))
    print(r'\\')
    cond_table = table(cond_table, header_row=True, frame=True)
    print()
    print(r'\begin{scriptsize}')
    print(latex(cond_table))
    print(r'\end{scriptsize}')

print(r'\newpage')
print(r'\scriptsize')
print(latex(table(general_table, header_row=True, frame=True)))

sr.texdoc_ending()
