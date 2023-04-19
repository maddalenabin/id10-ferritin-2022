import re

def get_scan_number(s):
    s = str(s)
    s = re.search('(?<=scan)\d{4}', s)
    if s is not None:
        return int(s.group(0))
    else:
        return -1

def get_datdir_number(s):
    s = str(s)
    s = re.search(r"(?<=conc120_gly_50_2_)\d{4}", s)
    if s is not None:
        return int(s.group(0))
    else:
        return    
    
def get_rep(x, reps_per_spot=1):
    scan = get_scan_number(x)
    return (scan % reps_per_spot) + 1