from sage.typeset.ascii_art import AsciiArt

from src.space_groups import prepare_gap_env
from src.srdegrees import SR_Degrees

prepare_gap_env()

AsciiArt._terminal_width = lambda self: 200

print('===================== planar case =========================')
for i in [4, 7, 8, 12]:
    print(f'\n\n\n ---------------------- Group {i} -------------------')
    t = SR_Degrees(i, verbose=2, method='markdown')
    ps = t.G.point_group_normalizer()

    for A in ps:
        print('---')
        print('matrix A:')
        print(A)
        A_inv = A.inverse().simplify_rational()

        eqs, base_vars = t.construct_congruences_v2(A_inv, A)
        print("eqs:")
        print(eqs, base_vars)
        print('solutions:')
        print(t.solve_congruences_v4(eqs, base_vars, list()))


print('==================== spatial case ==========================')

for i in [4, 12, 50, 55, 120, 220]:

    print(f'\n\n\n ---------------------- Group {i} -------------------')
    t = SR_Degrees(i, dim=3, verbose=2, method='markdown')
    if t.G.is_symmorphic():
        continue

    ps = t.G.point_group_normalizer()

    for A in ps:
        print('----')
        print('matrix A:')
        print(A)
        A_inv = A.inverse().simplify_rational()

        eqs, base_vars = t.construct_congruences_v2(A_inv, A)
        print('solutions:')
        print(t.solve_congruences_v4(eqs, base_vars, list()))
