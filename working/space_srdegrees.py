
from src.space_groups import prepare_gap_env
from src.srdegrees import SR_Degrees, solve_simple_mat


prepare_gap_env()

sr = SR_Degrees(22, dim=3)

print(*sr.G.point_group_normalizer(), sep='\n\n')

# sc = sr.sc_degrees()

# for A, _, _, _ in sc:
    # print(solve_simple_mat(A))
