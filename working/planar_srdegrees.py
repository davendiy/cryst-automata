

from src.space_groups import SpaceGroup_gap, prepare_gap_env
from src.srdegrees import SR_Degrees


prepare_gap_env()

for n in range(1, 18): 
#for n in [12]:
    #G = SpaceGroup_gap.from_gap_cryst(n, dim=2)

    sr = SR_Degrees(n, method='latex')
    if n == 1:
        sr.texdoc_header()

    sr.sr_degrees()
    
    if n == 17:
        sr.texdoc_ending()
    
