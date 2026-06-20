
# from sage.all import latex
from src.space_groups import prepare_gap_env

from src.srdegrees import SR_Degrees
from src.cryst3.gen_norms3 import load_cached_normalizers


prepare_gap_env()


# for i in [24]:
# for i in [9]:
for i in range(2, 231):
    print('\n\n')
    print(f'==================== group #{i} ======================')
    t = SR_Degrees(group_index=i, dim=3, verbose=2)
    # for A in load_cached_normalizers(i):
    for A in t.G.point_group_normalizer(enforce_integral=True):
        print()
        print('original A:')
        print(A)

        if t.G.is_symmorphic():
            continue

        # t.construct_congruences_v2(A.inverse().simplify_rational(), A)

        m = t.cocycle_compat_v2(A)
        if m.status == m.status.Error:
            print('no solutions')
            continue
        print('cocycle compat A:')
        print(m.result)
        print()
        print('det(A) =', m.result.det().simplify_rational())
        print('tr(A) =', m.result.trace().simplify_rational())

