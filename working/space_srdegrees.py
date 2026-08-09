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


def collect_groups_case1(rowtable):
    unique_triples = {}
    result = []

    for row in rowtable:
        groupnum, p, q, g = row
        f_hash = str(p) + str(q) + str(g)
        f_altern = []
        for sp, sg in [[1, 1], [1, -1], [-1, 1], [-1, -1]]:
            # q can't alternate due to condition (q+-g) != -1
            f_altern.append(str(sp * p) + str(q) + str(sg * g))

        if f_hash not in unique_triples:
            base_struct = [set(), p, q, g]
            for h in f_altern:
                unique_triples[h] = base_struct
            result.append(base_struct)

        grouplist, *_ = unique_triples[f_hash]
        grouplist.add(groupnum)
    return result


def collect_groups_case2(rowtable):
    unique_polys = {}
    count = 1
    result = []
    for row in rowtable:
        groupnum, f = row
        f_hash = str(f)
        f_altern = str(-f)

        if f_hash not in unique_polys:
            base_struct = [count, f, set()]
            count += 1
            unique_polys[f_hash] = base_struct
            unique_polys[f_altern] = base_struct
            result.append(base_struct)

        _, _, grouplist = unique_polys[f_hash]
        grouplist.add(groupnum)
    return result


def collect_by_polys(rowtable):
    unique_triples = {}
    result = []

    for row in rowtable:
        groupnum, p, q, g = row
        for i, f in enumerate([p, q, g]):
            f_hash = str(f)
            f_alt = str(-f)

            if f_hash not in unique_triples:
                base_struct = [f, set(), set(), set()]
                unique_triples[f_hash] = base_struct
                unique_triples[f_alt] = base_struct
                result.append(base_struct)

            unique_triples[f_hash][i+1].add(groupnum)
    return sorted(result, key=lambda x: len(str(x[0])))


def normalize_grouplist(groups):
    groups = list(sorted(groups))
    if not groups:
        return ''
    elif len(groups) <= 4:
        return table([[','.join([str(el) for el in groups])]])
    elif len(groups) <= 40:
        return table(rows=[[','.join(str(el) for el in groups[i:i+4])] for i in range(0, len(groups), 4)])
    else:
        return f'{groups[0]},{groups[1]}...{groups[-2]},{groups[-3]}'


def printv(*args, **kwargs):
    if VER:
        print(*args, **kwargs)


headers = ['Groups No.', r'$p_\sigma$', r'$q_\sigma$', r'$g_\sigma$']
headers_simple = ['No.', r'$f_\sigma$', 'Groups No.']
res_table_simple = []
res_table = []

# for i in [24]:
# for i in [9]:
# for i in [148]:
for i in tqdm(range(2, 231)):
# for i in tqdm(range(2, 40)):
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

table_polys = collect_by_polys(res_table)
res_table = collect_groups_case1(res_table)
for i in range(len(res_table)):
    for j in range(1, 4):
        res_table[i][j] = '$' + latex(res_table[i][j]) + '$' if res_table[i][j] is not None else ''
    res_table[i][0] = normalize_grouplist(res_table[i][0])

res_table_simple = collect_groups_case2(res_table_simple)
for i in range(len(res_table_simple)):
    res_table_simple[i][1] = '$' + latex(res_table_simple[i][1]) + '$' if res_table_simple[i][1] is not None else ''
    res_table_simple[i][2] = normalize_grouplist(res_table_simple[i][2])

for i in range(len(table_polys)):
    f, grp, grq, grg = table_polys[i]
    f = '$' + latex(f) + '$' if f is not None else ''
    grp = normalize_grouplist(grp)
    grq = normalize_grouplist(grq)
    grg = normalize_grouplist(grg)

    table_polys[i] = [i,  f, grp, grq, grg]


with open(filename, 'w') as file:
    i = 0
    while i < len(res_table):
        chunk_size = 0
        j = i
        while chunk_size < M and j < len(res_table):
            grlist = res_table[j][0]
            gsize = 1 if isinstance(grlist, str) else len(latex(grlist).splitlines()) - 2
            if gsize > M+1:
                raise ValueError
            if chunk_size + gsize > M + 1:
                break
            chunk_size += gsize
            j += 1
        chunk = res_table[i: j]

        print(r'\begin{table}[H]', file=file)
        print(r'\begin{center}', file=file)
        print(r'\begin{scriptsize}', file=file)
        print(latex(table([headers] + chunk, frame=True, header_row=True)), file=file)
        print(r'\end{scriptsize}', file=file)

        print(r'\caption{\textit{SRD}' + f' unique triplets for the block-triangular case {i+1}..{i + len(chunk)}' + r'}', file=file)
        print(r'\label{tab:case1_srdegrees' + f'{i}' + '}', file=file)
        print(r'\end{center}', file=file)
        print(r'\end{table}', file=file)
        print('', file=file)
        print(r'\newpage', file=file)
        print('', file=file)

        i = j

    for i in range(0, len(res_table_simple), M):
        print(r'\begin{table}[H]', file=file)
        print(r'\begin{center}', file=file)
        print(r'\begin{scriptsize}', file=file)
        chunk = res_table_simple[i: i+M]
        print(latex(table([headers_simple] + chunk, frame=True, header_row=True)), file=file)
        print(r'\end{scriptsize}', file=file)

        print(r'\caption{\textit{SRD}' + f' unique polynomials for the reduced matrix case {i+1}..{i + len(chunk)}' + r'}', file=file)
        print(r'\label{tab:case2_srdegrees' + f'{i//M}' + '}', file=file)
        print(r'\end{center}', file=file)
        print(r'\end{table}', file=file)
        print('', file=file)
        print(r'\newpage', file=file)
        print('', file=file)

    # headers = ['num', 'polynomial', r'$p_\sigma$ groups', r'$q_\sigma$ groups', r'$g_\sigma$ groups']

    # while i < len(table_polys):
    #     chunk_size = 0
    #     j = i
    #     while chunk_size < M and j < len(table_polys):
    #         max_gsize = 0
    #         for p in [2, 3, 4]:
    #             grlist = table_polys[j][p]
    #             gsize = 1 if isinstance(grlist, str) else len(latex(grlist).splitlines()) - 2
    #             if gsize > M+1:
    #                 raise ValueError
    #             max_gsize = max(max_gsize, gsize)

    #         if chunk_size + max_gsize > M + 1:
    #             break
    #         chunk_size += gsize
    #         j += 1
    #     chunk = table_polys[i: j]

    #     print(r'\begin{table}[H]', file=file)
    #     print(r'\begin{center}', file=file)
    #     print(r'\begin{scriptsize}', file=file)
    #     print(latex(table([headers] + chunk, frame=True, header_row=True)), file=file)
    #     print(r'\end{scriptsize}', file=file)

    #     print(r'\caption{' + 'Unique polynomials ' + f'{i+1}..{i + len(chunk)}' + r'}', file=file)
    #     print(r'\end{center}', file=file)
    #     print(r'\end{table}', file=file)
    #     print('', file=file)
    #     print(r'\newpage', file=file)
    #     print('', file=file)

    #     i = j
