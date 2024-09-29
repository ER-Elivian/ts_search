import numpy as np
import plotly.graph_objects as go
import os
import statistics

def make3list(name):
    with open(os.path.join(rpath, name),"r") as file:
        f_strs=file.readlines()
    
    res_str=""
    for f_str in f_strs:
        res_str+=f_str
    f_strs=res_str.split("\n")

    x=[]
    f_str_spl=f_strs[0].split()
    for i in range(len(f_str_spl)):
        x.append(float(f_str_spl[i]))
    
    y=[]
    f_str_spl=f_strs[1].split()
    for i in range(len(f_str_spl)):
        y.append(float(f_str_spl[i]))
    
    z=[]
    f_str_spl=f_strs[2].split()
    for i in range(len(f_str_spl)):
        z.append(float(f_str_spl[i]))
    
    return x,y,z
    
rpath=os.path.join(os.getcwd(),"tests","ep_test")
x,y,z=make3list("bond_lens")

fig = go.Figure()
fig.add_scatter3d(x=x,y=y,z=z,  mode='markers+lines', marker=dict(size=1.5,color="red"), line=dict(color="black"),text=[str(i) for i in range(len(x))], name="way")
fig.add_scatter3d(x=[1.65117],y=[2.07427],z=[2.07378],  mode='markers', marker=dict(size=2,color="blue"), name="TS")

v_max=[max(x),max(y),max(z)]
v_min=[min(x),min(y),min(z)]

vec=[np.multiply(0.5,np.add(v_max,v_min))]
div=min(np.subtract(v_max,v_min))/10
vec.append(np.add(vec[0],[div,div,div]))

fig.add_scatter3d(x=[vec[0][0],vec[1][0]],y=[vec[0][1],vec[1][1]],z=[vec[0][2],vec[1][2]],  mode='lines', line=dict(width=1,color="green"), name="dir x+y+z=const")
fig.update_layout(width=900, height=800, showlegend=False,
                  scene=dict(
                            xaxis_title='8, 7',
                            yaxis_title='1, 8',
                            zaxis_title='4, 8',
                            aspectmode="data"
                        ))
fig.show()