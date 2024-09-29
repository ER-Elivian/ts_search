import os
from pyxyz import Confpool
from TS_find import optTS
def opt_and_compare(rpath):
    got_TS=optTS(os.path.join(rpath,"to_opt.xyz"), threshold_rel=3, threshold_force=0.00004, print_output=True,mode="strict", maxstep=10**3, programm=dict(name="xtb", force_constant=6))
    p = Confpool()
    p.include_from_file(os.path.join(rpath,"TS.xyz"))
    p.include_from_file(os.path.join(rpath,"xtbopt.xyz"))
    p.generate_connectivity(0, mult=1.3, ignore_elements=[])
    p.generate_isomorphisms()
    rmsd_value= p[0].rmsd(p[len(p)-1])[0]
    print(rmsd_value)
    return rmsd_value, got_TS.settings["step"]

def print_result(result):
    color_by_value=lambda val :"\033"+ ("[92mgood"if val<0.001 else "[93mnot bad" if val<0.02 else "[91mbad") + "\033[00m"
    for key in result.keys():
        print(f'{key.ljust(20, " ")} {"{:9.7f}".format(result[key][0])} {color_by_value(result[key][0]).ljust(17, " ")} in {"{:6}".format(result[key][1])} steps')
        

rpaths=["da_test", "ep_test","sn2_test","bul_test","apw_test","apw2_test"]
result={}
for rpath in rpaths:
    print("   ")
    print(rpath)
    result[rpath]=opt_and_compare(os.path.join(os.getcwd(),"tests",rpath))
print_result(result)