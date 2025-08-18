
from sage.all import latex

from src.space_groups import SpaceGroup_gap, prepare_gap_env
from src.srdegrees import SR_Degrees


prepare_gap_env()

for n in range(2, 18): 
#for n in [12]:
    #G = SpaceGroup_gap.from_gap_cryst(n, dim=2)

    sr = SR_Degrees(n, method='latex', verbose=1)
    if n == 2:
        sr.texdoc_header()

    sr_table = sr.sr_degrees()
    print(latex(sr_table))
    if n == 17:
        sr.texdoc_ending()
    
