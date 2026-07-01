
from src.space_groups import prepare_gap_env
from src.srdegrees import SR_Degrees

prepare_gap_env()

print('===================== planar case =========================')
for i in [4, 7, 8, 12]:
    print(f'\n\n\n ---------------------- Group {i} -------------------')
    t = SR_Degrees(i, dim=2, verbose=2, method='ascii')
    ps = t.G.point_group_normalizer()

    for A in ps:
        print('---')
        print('matrix A:')
        print(A)
        A_inv = A.inverse().simplify_rational()

        left, right = t.construct_congruences_v3(A_inv, A)
        print('--------------congruences-----------------')
        print("left part:")
        print(left)
        print("right:")
        print(right)

        print('----------------answer----------------')
        res = t.solve_congruences_v5(left, right)
        if not res.failed:
            ans = res.result
            A_new = A.subs({var: exp for var, exp in zip(list(ans.base_variables), list(ans.expressions))})
            print('result matrix:')
            print(A_new)
        else:
            print('no solutions')


# input()
print('==================== spatial case ==========================')
for i in [9, 11]:        # , 13, 21, 24, 64, 149]:

    print(f'\n\n\n ---------------------- Group {i} -------------------')
    t = SR_Degrees(i, dim=3, verbose=2, method='ascii')
    if t.G.is_symmorphic():
        continue

    ps = t.G.point_group_normalizer()

    for A in ps:
        print('---')
        print('matrix A:')
        print(A)
        A_inv = A.inverse().simplify_rational()

        left, right = t.construct_congruences_v3(A_inv, A)
        print('--------------congruences-----------------')
        print("left part:")
        print(left)
        print("right:")
        print(right)

        print('----------------answer----------------')
        res = t.solve_congruences_v5(left, right)
        if not res.failed:
            ans = res.result
            A_new = A.subs({var: exp for var, exp in zip(list(ans.base_variables), list(ans.expressions))})
            print('result matrix:')
            print(A_new)
        else:
            print('no solutions')
