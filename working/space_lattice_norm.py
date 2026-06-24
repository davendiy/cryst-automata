
from base64 import b64encode

import os
import sys

from src.space_groups import prepare_gap_env

from src.srdegrees import SR_Degrees
from src.cryst3.gen_norms3 import load_cached_normalizers, dumps, point_groups_representatives

prepare_gap_env()


save_folder = sys.argv[1]
os.makedirs(save_folder, exist_ok=True)

for i in point_groups_representatives:
    print(f'==================== group #{i} ======================')

    t = SR_Degrees(group_index=i, dim=3, verbose=0)
    with open(os.path.join(save_folder, f'norm_{i}.txt'), 'w') as file:
        for A in load_cached_normalizers(i):
            B = t.lattice_compat(A).result
            file.write(str(B))
            file.write('\n\n')
            file.write('\n######## INTERNAL_REPR #######\n')
            file.write(b64encode(dumps(B)).decode('ascii'))
            file.write('\n#### END OF INTERNAL_REPR ####\n')
