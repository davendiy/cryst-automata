from src.srdegrees import SR_Degrees
from src.space_groups import prepare_gap_env

prepare_gap_env()

for i in range(2, 18):
    SR_Degrees(i, method='markdown').sr_degrees()
