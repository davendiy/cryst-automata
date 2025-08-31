from sage.all import latex, table

from src.space_groups import prepare_gap_env
from src.srdegrees import SR_Degrees


prepare_gap_env()

giant_table = [["num", "srdegrees"]]

sr = SR_Degrees(1, method='latex', verbose=1)
sr.texdoc_header()

for n in range(2, 18):
    # for n in [12]:

    sr = SR_Degrees(n, method='latex', verbose=1)
    eig_table, res_table, cond_table = sr.sr_degrees()
    eig_table = table(eig_table, header_row=True, frame=True)

    smaller_table = []
    for el in res_table[1:]:
        inz = ",".join([latex(_x) for _x in el[1].variables()]) + sr.in_sym + sr.z  # type: ignore
        inz = r"\{" + latex(el[1]) + r'\, | ' + inz + r'\}'

        if not el[2]:
            smaller_table.append(["$" + inz + "$"])

        new_arr = r"\{"
        new_arr += ','.join(latex(tuple(sol)) for sol in el[2])
        new_arr += r'\}'
        el[2] = "$" + new_arr + "$"

        smaller_table.append(
            [
                r"$" + inz + " / " + new_arr + "$",
            ]
        )

    if smaller_table:
        giant_table.append([n, table(smaller_table)])  # type: ignore

    res_table = table(res_table, header_row=True, frame=True)

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

    cond_table = table(cond_table, header_row=True, frame=True)
    print()
    print(r'\begin{small}')
    print(latex(cond_table))
    print(r'\end{small}')

print(r'\newpage')
print(r'\scriptsize')
print(latex(table(giant_table, header_row=True, frame=True)))

sr.texdoc_ending()
