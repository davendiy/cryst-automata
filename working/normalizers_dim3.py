
from src.space_groups import prepare_gap_env
from src.cryst3.common import cached_simple_matrices


prepare_gap_env()


for gr_num in cached_simple_matrices:
    print(f'################# {gr_num} ######################')
    for a in cached_simple_matrices[gr_num]:
        # if tau(a) == 0:
        #     continue
        # print(a)
        print("equation:", (a.det() - a.trace()).simplify_rational())
        # print('tau(A):', tau(a))
