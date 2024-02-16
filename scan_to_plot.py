import os,numpy as np

initial_cwd = os.getcwd()
path=os.path.join(initial_cwd,"scan_DA")

scanname="scan_local"
curvename="way_log"
scan_vals=[]
curve_vals=[]
with open(os.path.join(path,scanname+".txt"),"r") as file:
    for line in file:
        linesplit=line.split()
        scan_vals.append([float(linesplit[0]), float(linesplit[1]), float(linesplit[2])])

with open(os.path.join(path,curvename+".txt"),"r") as file:
    for line in file:
        linesplit=line.split()
        curve_vals.append([float(linesplit[0]), float(linesplit[1]), float(linesplit[2])])

minvals=np.array(scan_vals[0])
maxvals=np.array(scan_vals[0])
for val in scan_vals:
    minvals=[min(minvals[0],val[0]), min(minvals[1],val[1]), min(minvals[2],val[2])]
    maxvals=[max(maxvals[0],val[0]), max(maxvals[1],val[1]), max(maxvals[2],val[2])]



abs_hs=np.subtract(maxvals,minvals)
medians=0.5*np.add(maxvals,minvals)
with open(os.path.join(path,scanname+"_plot.txt"),"w+") as file:
    for val in scan_vals:
        val=20*np.subtract(np.array(val),medians)/abs_hs
        file.write(f"{val[0]} {val[1]} {val[2]}\n")
#        print(val)
        
with open(os.path.join(path,curvename+"_"+scanname+"_plot.txt"),"w+") as file:
    for val in curve_vals:
        val=20*np.subtract(np.array(val),medians)/abs_hs
        file.write(f"{val[0]} {val[1]} {val[2]}\n")
        print(val)
