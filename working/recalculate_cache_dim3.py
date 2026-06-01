import sys

from src.space_groups import prepare_gap_env, SpaceGroup_gap
from src.cryst3.gen_norms3 import bruteforce_normalizers

prepare_gap_env()

folder = sys.argv[1]

bruteforce_normalizers(folder)
