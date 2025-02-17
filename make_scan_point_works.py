import os
rpath=os.path.join(os.getcwd(),"da_scan")
structsname="scan_structs_points.xyz"
bts_content="0\nvacuum\nb 1 11 1\nb 4 12 1"
with open(os.path.join(rpath,structsname), "r") as xyzlog:
    lines_xyzs=xyzlog.readlines()
num=int(lines_xyzs[0])

wpath=os.path.join(os.getcwd(),"scan_opt_da")
n=0#string number in lines_xyzs
k=0#xyz number
while n+num+2<=len(lines_xyzs) and k<500:#k<500 just for not folder overflow
    xyz_lines=lines_xyzs[n:n+num+2]
    path_wf=os.path.join(wpath,f"work{k}")
    if not os.path.exists(path_wf):
        os.mkdir(path_wf)
    path_w=os.path.join(path_wf,"to_opt.xyz")
    with open(path_w,"w+") as work:
        work.writelines(xyz_lines)
    path_bts=os.path.join(path_wf,"bonds_to_search")
    with open(path_bts,"w+") as bts:
        bts.writelines(bts_content)
    k+=1
    n+=num+2
    
    