from src.space_groups import prepare_gap_env
from src.srdegrees import SR_Degrees


prepare_gap_env()


for i in [4, 7, 8, 12]:
    print(f'\n\n\n         Group {i}             ')
    t = SR_Degrees(i, verbose=0)
    ps = t.G.point_group_normalizer()

    for A in ps:
        print('matrix A:')
        print(A)
        A_inv = A.inverse().simplify_rational()

        eqs, base_vars, variables = t.construct_congruences(A_inv, A)
        print('solutions:')
        t.solve_congruences_v3(eqs, base_vars, variables)
