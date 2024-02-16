import os,subprocess,numpy as np

vec_len = lambda v:(np.matmul(v,v.T))**0.5
control_top = lambda chrg,control_strs: control_strs.extend([f"$chrg {chrg}\n","$constrain\n"]) #заполнить верхнюю часть control файла
control_end = lambda f_c,control_strs: control_strs.extend([f"    force constant = {f_c}\n","$end\n"])#заполнить нижнюю часть control файла

def extract_AB_dir(xyzs_strs, num_A,num_B):
    vec_A=xyzs_strs[num_A+1].split()[1:]
    for num,coord in enumerate(vec_A):
        vec_A[num]=float(coord)

    vec_B=xyzs_strs[num_B+1].split()[1:]
    for num,coord in enumerate(vec_B):
        vec_B[num]=float(coord)
    
    res=np.subtract(vec_B,vec_A)
    return res

def increment_bond(scans,counters):
    counters[len(counters)-1]+=1
    for i in range(len(scans)-1,-1,-1):
        if counters[i]>scans[i][4]:
            counters[i]=0
            if i>0:
                counters[i-1]+=1
            else:
                return "full"
            

initial_cwd = os.getcwd()
path=os.path.join(initial_cwd,"scan_DA")

print(path)

scans=[]#nA, nB, init_len, fin_len, steps

reaction="to_scan.xyz"
rpath=path
os.chdir(path)
scanname="scan_local"
with open(os.path.join(rpath,scanname+".txt"), "w+") as scanfile:
    a=0

with open(os.path.join(rpath, "scans"),"r") as scansfile:
    line=scansfile.readline()
    while(line!=""):
        linesplit=line.split()
        scans.append([int(linesplit[0]),int(linesplit[1]),float(linesplit[2]),float(linesplit[3]),int(linesplit[4])])
        line=scansfile.readline()

print(scans)
counters=[]
for i in range(len(scans)):
    counters.append(0)

line=0
while line!="full":
    control_strs=[]
    length_str=""

    control_top(0,control_strs)
    for i,scanval in enumerate(scans):
        length=scanval[2]+(scanval[3]-scanval[2])*counters[i]/scanval[4]
        control_strs.append(f"    distance: {scanval[0]}, {scanval[1]}, {length}\n")  
    control_end(6,control_strs)

    with open(os.path.join(rpath,"control"),"w+") as control:
        control.writelines(control_strs)
    with open(os.path.join(rpath,"xtbout"),"w+") as xtbout:
        subprocess.call(["xtb", reaction, "-I", "control","--opt","--vtight"],stdout=xtbout)
        #subprocess.call(["xtb", "xtbopt.xyz", "--chrg", "1","--alpb","water","--acc","0.5","--grad"],stdout=xtbout)
    
    xyzs_strs=[]
    with open(os.path.join(rpath,"xtbopt.xyz"), "r") as file:
        xyzs_strs=file.readlines()
    for i in range(len(scans)):    
        vec=extract_AB_dir(xyzs_strs,scans[i][0],scans[i][1])
        length=vec_len(vec)
        length_str+=f"{length} "
    
    energy=xyzs_strs[1].split()[1]
    length_str+=f"{energy}\n"
    
    with open(os.path.join(rpath,scanname+".txt"), "a") as scanfile:
        scanfile.write(length_str)
    print(counters)
    line=increment_bond(scans,counters)


