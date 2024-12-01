
import os,shutil
import TS_find

BtS_base="0\nvacuum\n"

folder_work="tests_neb"
wd=os.getcwd()

wpath=os.path.join(wd,folder_work)

listworks=os.listdir(wpath)
for work in listworks:
    TS_find.optTS(xyz_path=os.path.join("tests_neb", work, "to_opt.xyz"), threshold_rel=8, threshold_force=0.00004, print_output=True,mode="strict", maxstep=10**3, programm=dict(name="xtb", force_constant= 6))
    