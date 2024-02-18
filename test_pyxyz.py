import os
from pyxyz import Confpool
import opt_one
def opt_and_compare(rpath):
    opt_one.find_TS(rpath, "to_opt.xyz",optimized_cap=0.0001, ratio=300, maxstep=1000)
    p = Confpool()
    p.include_from_file(os.path.join(rpath,"TS.xyz"))
    p.include_from_file(os.path.join(rpath,"xtbopt.xyz"))
    p.generate_connectivity(0, mult=1.3, ignore_elements=[])
    p.generate_isomorphisms()
    rmsd_value= p[0].rmsd(p[len(p)-1])[0]
    print(rmsd_value)
    return rmsd_value


rpaths=[ "da_test", "ep_test","sn2_test", "et_test","apw_test", "apw2_test","apw3_test"]
result={}
for rpath in rpaths:
    print("   ")
    print(rpath)
    result[rpath]=opt_and_compare(os.path.join(os.getcwd(),"tests",rpath))
print(result)