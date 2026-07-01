
from tqdm import tqdm
from src.space_groups import prepare_gap_env, SpaceGroup_gap

prepare_gap_env()

snots = set()

for i in tqdm(range(1, 231)):
    G = SpaceGroup_gap.from_gap_cryst(i, dim=3, change_basis=True)
    for el in G._alpha.values():
        snots.add(el.denominator())

print(*snots, sep='\n\n')
