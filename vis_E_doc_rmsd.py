import numpy as np
import plotly.graph_objects as go
import os
import statistics

def makelist(name):
    with open(os.path.join(rpath, name),"r") as file:
        f_strs=file.readlines()
    aslist=[]
    for f_str in f_strs:
        aslist.append(float(f_str))
    return aslist
    
rpath=os.path.join(os.getcwd(),"tests","ep_test")
Energ=makelist("log_E")
rmsd=makelist("log_rmsd")
doc=makelist("log_doc")

with open(os.path.join(rpath, "TS.xyz"),"r") as TSfile:
    lines=TSfile.readlines()
E_TS=float(lines[1].split()[1])

fig = go.Figure()
fig.add_scatter3d(x=Energ+[-21.532561913404],y=rmsd,z=doc,  mode='markers+lines', marker=dict(size=1.5,color="red"), line=dict(color="black"), text=[str(i) for i in range(len(Energ))],name="way")
fig.add_scatter3d(x=[E_TS],y=[0],z=[0.00004],  mode='markers', marker=dict(size=2,color="blue"), name="TS")

fig.update_layout(width=900, height=800, showlegend=False,
                  scene=dict(
                            xaxis_title='Energy',
                            yaxis_title='rsmd to TS',
                            zaxis_title='div of forces',
                        ))
fig.show()