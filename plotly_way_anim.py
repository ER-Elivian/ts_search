import numpy as np
import plotly.graph_objects as go
import os, copy
import  io
from PIL import Image
b2a=0.529177208
def read_way(way_path):
    athom_frames=[]
    with open(os.path.join(way_path)) as file:
        xyz_log_lines=file.readlines()
    n_athoms=int(xyz_log_lines[0])
    struct_num=0
    while (n_athoms+2)*struct_num+2 < len(xyz_log_lines):
        athom_lines=xyz_log_lines[(n_athoms+2)*struct_num+2 : (n_athoms+2)*(struct_num+1)]
        athoms=[]#symbol, X, Y, Z
        for i in range(len(athom_lines)):
            linesplit=athom_lines[i].split()
            athoms.append([linesplit[0], float(linesplit[1]), float(linesplit[2]), float(linesplit[3])])
        athom_frames.append(athoms)
        struct_num+=1
    return athom_frames

def color_and_size_by_symb(symb):
    if symb=="O":
        return "red", 11*b2a
    if symb=="C":
        return "black", 11*b2a
    if symb=="H":
        return "gray", 6*b2a
    if symb=="N":
        return "blue", 10*b2a
    return "magenta", 16*b2a

def bond_val(v1, v2): #by chemcraft
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

def find_connections(athom_frames):
    connection_frames=[]
    n_ath=len(athom_frames[0])
    dist=lambda a1,a2: ((a1[1]-a2[1])**2+(a1[2]-a2[2])**2+(a1[3]-a2[3])**2)**0.5
    
    for frame_num in range(len(athom_frames)):
        connections=np.empty((n_ath,n_ath),dtype=int)
        for i in range(n_ath):
            for j in range(i+1,n_ath):
                if dist(athom_frames[frame_num][i], athom_frames[frame_num][j]) < bond_val(athom_frames[frame_num][i][0],athom_frames[frame_num][j][0])*1.2:
                    connections[i][j]=1
                else:
                    connections[i][j]=0
        connection_frames.append(connections)
    return connection_frames

def ends_move(a1,a2):
    v=np.subtract(a1[1:],a2[1:])
    v_len=np.linalg.norm(v)
    v=np.multiply(1/v_len,v)
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

def draw_on_figure(athom_frames=None, connection_frames=None):
    frames=[]
    frame_figures=[]
    base_fig = go.Figure()
    for frame_num in range(len(athom_frames)):
        frame_fig=go.Figure()
        athoms=athom_frames[frame_num]
        if athoms!=None:
            x_arr, y_arr, z_arr, m_arr = [],[],[],dict(color=[], size=[])
            for i,athom in enumerate(athoms):
                c, s=color_and_size_by_symb(athom[0])
                x_arr.append(athom[1])
                y_arr.append(athom[2])
                z_arr.append(athom[3])
                m_arr["color"].append(c)
                m_arr["size"].append(s*7)
                m_arr["opacity"]=1
            
            frame_fig.add_trace(go.Scatter3d(x=np.array(x_arr), 
                                             y=np.array(y_arr), 
                                             z=np.array(z_arr),name=''))
            if frame_num==0:
                base_fig.add_trace(go.Scatter3d(x=np.array(x_arr), y=np.array(y_arr), z=np.array(z_arr),mode="markers",marker=m_arr))
        
        connections=connection_frames[frame_num]
        if type(connections)!=None:
            n_ath=len(athoms)
            for i in range (n_ath):
                for j in range(i,n_ath):
                    if connections[i][j]==1:
                        p1,p2=ends_move(athoms[i],athoms[j])
                        if frame_num==0:
                            base_fig.add_trace(go.Scatter3d(x=[p1[0],p2[0]],
                                                   y=[p1[1],p2[1]],
                                                   z=[p1[2],p2[2]],
                                                   mode='lines',
                                                   line=dict( color = "rgb(50,50,50)",width = 3)
                                                   
                                                   ))
                        frame_fig.add_trace(go.Scatter3d(x=[p1[0],p2[0]],
                                      y=[p1[1],p2[1]], 
                                      z=[p1[2],p2[2]]))
                        
        frame=go.Frame(data=frame_fig["data"])
        frame_figures.append(frame_fig)
        frames.append(frame)
        print(f"frame {frame_num}")
    
    base_fig.update(frames=frames)
    return base_fig, frame_figures

def save_gif(name, fig):
    frame_images=[]
    print("converting to images")
    for slider_pos, frame in enumerate(fig.frames):
        if slider_pos%10==0:
            print(f"{slider_pos+1} of {len(fig.frames)}")
        fig.update(data=frame.data)
        frame_images.append(Image.open(io.BytesIO(fig.to_image(format="png"))))
    
    print("writing")
    frame_images[0].save(name,
               save_all=True,
               append_images=frame_images[1:],
               optimize=True,
               duration=650,
               loop=0)
    
initial_cwd = os.getcwd()
rpath="tests/da_test/TS_search_log.xyz"
print("reading")
athom_frames=read_way(rpath)
print("conneting")
connects=find_connections(athom_frames)

print("drawing")
fig, frame_figures=draw_on_figure(athom_frames, connects)

fig.update()
camera = dict(
    eye=dict(x=0, y=0.9, z=1.1)
)

fig.update_layout(#updatemenus=[dict(type="buttons",
                  #        buttons=[dict(label="Play",
                  #                      method="animate",
                  #                      args=[None, dict(frame=dict(redraw=True,fromcurrent=True, mode='immediate'))      ])])],
                  scene=dict(aspectmode="data",
                             xaxis = dict(visible=False),
                             yaxis = dict(visible=False),
                             zaxis = dict(visible=False)),
                  scene_camera=camera,
                  showlegend=False)


save_gif("out.gif", fig)

#fig.show()
