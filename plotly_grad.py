import numpy as np
import plotly.graph_objects as go
import os, copy

b2a=0.529177208
def read_grad(grad_path):
    with open(os.path.join(grad_path,"gradient")) as file:
        gradlines=file.readlines()
    n_athoms=(len(gradlines)-3)//2
    print(n_athoms)
    athoms=[]#symbol, X, Y, Z
    for i in range(2, 2+n_athoms):
        linesplit=gradlines[i].split()
        athoms.append([linesplit[3], float(linesplit[0]), float(linesplit[1]), float(linesplit[2])])
    forces=[]
    for i in range(2+n_athoms, 2+n_athoms*2):
        linesplit=gradlines[i].split()
        forces.append([-float(linesplit[0]), -float(linesplit[1]), -float(linesplit[2])])
    return athoms, forces

def color_and_size_by_symb(symb):
    if symb=="O":
        return "red", 11
    if symb=="C":
        return "black", 11
    if symb=="H":
        return "gray", 6
    if symb=="N":
        return "blue", 10
    return "magneta", 16

def bond_val(v1, v2):
    if v1>v2:
        v1,v2=v2,v1
    if v1=="C":
        if v2=="C":
            return 1.54
        if v2=="H":
            return 1.08
        if v2=="N":
            return 1.47
        if v2=="O":
            return 1.44
    if v1=="H":
        if v2=="H":
            return 0.74
        if v2=="N":
            return 1.00
        if v2=="O":
            return 0.96
    if v1=="N":
        if v2=="N":
            return 1.45
        if v2=="O":
            return 1.41
    if v1=="O":
        if v2=="O":
            return 1.48

def find_connections(athoms):
    n_ath=len(athoms)
    dist=lambda a1,a2: ((a1[1]-a2[1])**2+(a1[2]-a2[2])**2+(a1[3]-a2[3])**2)**0.5
    connections=np.empty((n_ath,n_ath),dtype=int)
    for i in range (n_ath):
        for j in range(i+1,n_ath):
            if dist(athoms[i], athoms[j])*b2a < bond_val(athoms[i][0],athoms[j][0])*1.2:
                connections[i][j]=1
            else:
                connections[i][j]=0
    return connections

def ends_move(a1,a2):
    v=np.subtract(a1[1:],a2[1:])
    print(v)
    v_len=np.linalg.norm(v)
    print(v_len)
    v=np.multiply(1/v_len,v)
    print(np.linalg.norm(v))
    if a1[0]==0:
        mul_a1=0
    else:
        mul_a1=0.057*color_and_size_by_symb(a1[0])[1]

    if a2[0]==0:
        mul_a2=0
    else:
        mul_a2=0.057*color_and_size_by_symb(a2[0])[1]
    
    if a1[0]==0:
        p1=[a1[1]+v[0]*mul_a2, a1[2]+v[1]*mul_a2, a1[3]+v[2]*mul_a2]
        p2=[a2[1]+v[0]*mul_a2, a2[2]+v[1]*mul_a2, a2[3]+v[2]*mul_a2]
    elif a2[0]==0:
        p1=[a1[1]-v[0]*mul_a1, a1[2]-v[1]*mul_a1, a1[3]-v[2]*mul_a1]
        p2=[a2[1]-v[0]*mul_a1, a2[2]-v[1]*mul_a1, a2[3]-v[2]*mul_a1]
    else:        
        p1=[a1[1]-v[0]*mul_a1, a1[2]-v[1]*mul_a1, a1[3]-v[2]*mul_a1]
        p2=[a2[1]+v[0]*mul_a2, a2[2]+v[1]*mul_a2, a2[3]+v[2]*mul_a2]
    return p1,p2

def draw_on_figure(fig,athoms=None, forces=None, connections=None):
    if athoms!=None:
        for i,athom in enumerate(athoms):
            xx=[athom[1]]
            yy=[athom[2]]
            zz=[athom[3]]
            c, s=color_and_size_by_symb(athom[0])
            line_marker = dict(color=c, size=s)
        
            fig.add_scatter3d(x=xx, y=yy, z=zz, mode='markers',marker=line_marker,name="", text=f"{athom[0]}{i}")
            
    if forces!=None:
        max_force=0
        for force in forces:
            max_force=max(max_force, np.linalg.norm(force))
        for i in range(len(forces)):
            p1=copy.deepcopy(athoms[i])
            p2=[0,athoms[i][1]+forces[i][0]/max_force*2, athoms[i][2]+forces[i][1]/max_force*2, athoms[i][3]+forces[i][2]/max_force*2]
            
            p1,p2=ends_move(p1,p2)

            fig.add_scatter3d(x=[p1[0], p2[0]], y=[p1[1], p2[1]], z=[p1[2], p2[2]], mode='lines',
                              line = dict( color = "rgb(84,0,0)",width = 6),
                              hoverinfo="text+name",
                              text=f"{forces[i][0]} {forces[i][1]} {forces[i][2]}",
                              name=f'f{athoms[i][0]}{i}')
    if type(connections)!=None:
        n_ath=len(athoms)
        for i in range (n_ath):
            for j in range(i,n_ath):
                if connections[i][j]==1:
                    p1,p2=ends_move(athoms[i],athoms[j])
                    fig.add_scatter3d(x=[p1[0],p2[0]],
                                      y=[p1[1],p2[1]], 
                                      z=[p1[2],p2[2]], 
                                      mode='lines',
                                      line=dict( color = "rgb(50,50,50)",width = 3),
                                      hoverinfo='skip',
                                      name='')



initial_cwd = os.getcwd()
rpath=os.path.join("/media/user/D/MYDOCS/Projects/TS_search_old/visualization")
athoms, forces=read_grad(rpath)
connects=find_connections(athoms)
fig = go.Figure()
draw_on_figure(fig,athoms, forces, connects)
fig.update_layout(width=900, height=800, showlegend=False, scene=dict(aspectmode="data"))
fig.show()
#connections=connect(athoms)
