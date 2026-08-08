from sage.all import latex, table, var, matrix

from src.space_groups import prepare_gap_env
from src.srdegrees import SR_Degrees


def stack_solutions(pairs):
    if len(pairs) <= 6:
        return r'\left\{' + ','.join(latex(par) for par in pairs) + r'\right\}'
    else:
        n = 3
        res = r'\left\{\begin{array}{' + 'r'*n + r'}' + '\n'

        for i in range(0, len(pairs), n):
            row = [latex(el) + ',' for el in pairs[i: i + n]]
            if len(row) != n:
                row = row + ['' for _ in range(n - len(row))]
            else:
                row[-1] = row[-1].strip(',')

            assert len(row) == n
            res += ' & '.join(row)
            res += r' \\' + '\n'

        return res + r'\end{array}\right\}'


def uncover(el):
    if isinstance(el, list) or isinstance(el, tuple):
        if len(el) == 1:
            return el[0]
    return el


prepare_gap_env()

min_scd = [2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 2, 9, 3, 4, 4, 3, 3]
min_srd = [2, 2, 4, 6, 4, 2, 6, 3, 3, 2, 2, 9, 3, 4, 4, 3, 3]

FULL = False
verbose = 2 if FULL else 0

y0, y1, y2, y3 = var('y0 y1 y2 y3')
gen_mat = matrix([[y0, y1], [y2, y3]])

general_table = [
    [r"No. \footnotemark[1]{}", r"$A \in \mathcal{Y}$", r"$f_\sigma$", "$MinSRD$", r"$SRD$ restriction"],
    ["1", gen_mat, r"$|y_0y_3 - y_1y_2|$", min_srd[0], r"$\dots$",],
    ["2", gen_mat, r"$|y_0y_3 - y_1y_2|$", min_srd[1], r"$\dots$",],
]

sr = SR_Degrees(1, method='latex', verbose=verbose)

if FULL:
    sr.texdoc_header()

for n in range(2, 18):

    sr = SR_Degrees(n, method='latex', verbose=verbose)
    eig_table, res_table = sr.sr_degrees()
    eig_table = table(eig_table, header_row=True, frame=True)

    smaller_table = dict()
    used_polys = set()
    for el in res_table[1:]:

        if (el_n := len(el[1].variables())) > 1:    # type: ignore
            degree = f"^{el_n}"
            lp = '('
            rp = ')'
            var = '(' + ','.join(latex(v) for v in el[1].variables()) + ')'
        else:
            lp = rp = ''
            degree = ''
            var = latex(el[1].variables()[0])

        if latex(el[1]) in used_polys:
            continue

        # field = r'\mathbb{Z}^2'
        # # (x_1, x_2) \in Z^2
        # inz = lp + ",".join([latex(_x) for _x in el[1].variables()]) + rp + sr.in_sym + ' ' + field  # type: ignore
        # # {(x_1, x_2) \in Z^2}
        # inz = r"\left\{" + latex(el[1]) + r'\, | \,' + inz + r'\right\}'

        inz = latex(el[1].abs())

        used_polys.add(latex(el[1]))
        used_polys.add(latex(-el[1]))

        if not el[2]:
            smaller_table[("$" + inz + "$", "")] = el[0]
            continue

        # / {(-1, 1), (1, 1), (0, 0)}
        new_arr = var + r'\notin' + r"\left\{"
        new_arr = var + r'\notin' + stack_solutions([uncover(tuple(sol)) for sol in el[2]])
        el[2] = "$" + new_arr + "$"

        smaller_table[
                ('$' + inz + '$', '$' + new_arr + '$')
                # r"$" + inz.rstrip(r'\right\}') + r"\, / \," + new_arr + r'\right\}' + "$",
        ] = el[0]

    if smaller_table:
        # general_table.append([n, table([[el] for el in smaller_table])])  # type: ignore
        general_table.append([
            n,
            table([[v] for v in smaller_table.values()]),
            table([[el[0]] for el in smaller_table]),
            min_srd[n-1],
            table([[el[1]] for el in smaller_table]),
        ])

    res_table = table(res_table, header_row=True, frame=True)

    if FULL:
        print(latex(eig_table))
        print(r'\\')
        print(latex(res_table))
        print(r'\\')
        print()

    if sr.G.is_symmorphic():
        continue


# print(r'\newpage')
print(r'\begin{table}[H]')
print(r'\begin{tiny}')
print(latex(table(general_table, header_row=True, frame=True)))
print(r'\end{tiny}')
print(r'\caption{Self-covering and self-replicating degrees for plane groups.}')
print(r'\label{tab:planar_srdegrees}')
print(r'\end{table}')

print(r'\footnotetext[1]{A positional number of the group in the International Tables for Crystallography, Volume A (ITA)}')

if FULL:
    sr.texdoc_ending()
