import numpy as np,os,subprocess,datetime,time, networkx as nx
projection = lambda va,vb: np.multiply( np.matmul(va,vb.T)/np.matmul(vb,vb.T), vb )#a to b projection
vec_len = lambda v:(np.matmul(v,v.T))**0.5

def save_control(rpath:str,chrg:int,f_c:float, control_strs:list):
    with open(os.path.join(rpath,"control"),"w+") as control:
        control.writelines([f"$chrg {chrg}\n","$constrain\n"])
        control.writelines(control_strs)
        control.writelines([f"    force constant = {f_c}\n","$end\n"])
            
def sign(x:float):
    if x>0:
        return 1
    if x<0:
        return -1
    return 0

def opt(rpath:str,xyz_name:str,solvent:str):
    with open(os.path.join(rpath,"xtbout"),"w+") as xtbout:
        if solvent=="vacuum":
            subprocess.call(["xtb", xyz_name, "-I", "control","--vtight","--opt"],stdout=xtbout)
        else:
            subprocess.call(["xtb", xyz_name, "-I", "control","--alpb",solvent,"--opt","--vtight"],stdout=xtbout)

def o_grad(rpath:str,chrg:int,solvent:str):
    with open(os.path.join(rpath,"xtbout"),"w+") as xtbout:
        if solvent=="vacuum":
            subprocess.call(["xtb", "xtbopt.xyz", "--chrg", str(chrg),"--grad"],stdout=xtbout)
        else:
            subprocess.call(["xtb", "xtbopt.xyz", "--chrg", str(chrg), "--alpb", solvent,"--grad"],stdout=xtbout)

def adaptive_opt_cap(grad_strs:list, xyzs_strs:list, search_bonds:list,settings:dict):
    search_atoms=set()
    for bond in search_bonds:
        search_atoms.add(bond[0])
        search_atoms.add(bond[1])
    sum_forces=0
    num_forces=0
    for i in range(1,settings["nAtoms"]+1):
        if i not in search_atoms:
            vec_force=extractGradient(i+settings["nAtoms"]+1,grad_strs)
            sum_forces+=vec_len(vec_force)
            num_forces+=1
    return sum_forces/num_forces

def change_fn(length:float, cap:float,settings:dict):
    len_sb=settings["nBonds"]
    if len_sb>2:
        return length/(len_sb*2)
    len_sign=sign(length)
    x=len_sign*length
    change=(x**1.05)/2#(x**(0.1+1.9/(x**0.1+1)))*5
    if change>cap:
        change=cap
    
    return change*len_sign

def log_xyz(path:str, xyzs_strs:list, mode=None):
    if mode=="new":
        open_prm_str="w+"
    else:
        open_prm_str="a"
    with open(os.path.join(path,"TS_search_log.xyz"),open_prm_str) as log:
        log.writelines(xyzs_strs)

def log_forces(path:str, forces:str, mode=None):
    if mode=="new":
        open_prm_str="w+"
    else:
        open_prm_str="a"
    with open(os.path.join(path,"force_log.xyz"),open_prm_str) as log:
        log.write(forces)

def extract_AB_dir(xyzs_strs:list, num_A:int,num_B:int):
    vec_A=xyzs_strs[num_A+1].split()[1:]
    for num,coord in enumerate(vec_A):
        vec_A[num]=float(coord)

    vec_B=xyzs_strs[num_B+1].split()[1:]
    for num,coord in enumerate(vec_B):
        vec_B[num]=float(coord)

    res=np.subtract(vec_B,vec_A)
    return res

def extractGradient(line_num:int,grad_strs:list):
    gradline=grad_strs[line_num]
    gradstr_arr=gradline.split()
    gradarr=[]
    for i in range(len(gradstr_arr)):
        gradarr.append(float(gradstr_arr[i]))
    return np.array(gradarr)

def move_bonds( xyzs_strs:list, search_bonds:list, lens:dict, init_bonds:dict, grad_strs:list, control_strs:list, settings:dict):
    MIN_BOND=0.8
    MAX_BOND=3.1
    sum_changes=0
    min_change=1000
    max_change=-1000
    changes={}
    settings["step"]+=1
    cur_c_m=settings["change_mode"]
    for bond in search_bonds:#найдём градиент (желание растянуться) вдоль каждой связи
        num_A=bond[0]
        num_B=bond[1]
        key=f"{bond[0]}, {bond[1]}"

        grad_A=extractGradient(num_A+settings["nAtoms"]+1,grad_strs)
        grad_B=extractGradient(num_B+settings["nAtoms"]+1,grad_strs)
        AB_dir=extract_AB_dir(xyzs_strs,num_A,num_B)
        summ_grad=np.subtract(grad_B,grad_A)

        if not key in lens.keys():
            if settings["print_output"]:
                print(f"key \"{key}\" not in lens.keys()") 
            lens[key] = vec_len(AB_dir)

        s_g_proj=projection(summ_grad, AB_dir)
        proj_len=vec_len(s_g_proj)

        s_g_p_sign=sign(s_g_proj[0]*AB_dir[0])

        changes[key]=s_g_p_sign*proj_len#удлиннение (если отрицательно - укорочение) связи
        max_change=max(changes[key], max_change)
        min_change=min(changes[key], min_change)
        sum_changes+=changes[key]
    mean_change=sum_changes/settings["nBonds"]#среднее "желание растянуться"
    div_of_changes=max_change-min_change#отклонение 

    if settings["DoC_cutoff"]==0:
        settings["DoC_cutoff"]=max(div_of_changes*1.001, 0.001)
    if settings["print_output"]:
        print(f'div of forces: \033[01mcur\033[00m {"{:.8f}".format(div_of_changes)}    \033[01mprev\033[00m {"{:.8f}".format(settings["prev_dc"])} (\033{"[92mless" if div_of_changes<settings["prev_dc"] else "[091mhigher"}\033[00m)') 
        print(f'change mode {settings["change_mode"]}') 

    change_if_less_DoC=settings["preferred_change_mode"]*(change_fn(mean_change,0.015,settings))
    if abs(change_if_less_DoC)<0.000001:#предел разрешения control, а вблизи ПС тут получается очень маленькое значение (10^-8)
        change_if_less_DoC=sign(change_if_less_DoC)*0.000001
    if settings["print_output"]:
        print(f'\033[01m\033[93mbonds:\033[00m')
    for bond in search_bonds:
        key=f"{bond[0]}, {bond[1]}"
        bond_len=lens[key]
        
        if div_of_changes<settings["DoC_cutoff"]: # нужно "отпущение" связей
            bond_change=change_if_less_DoC
            settings["change_mode"]=-settings["preferred_change_mode"]
            settings["pass_turns"]=2#1 снимется уже в конце этого хода, так что 2
        elif div_of_changes<settings["prev_dc"] or settings["pass_turns"]>0:#нужно выравнивание по силам
            cur_c_m=settings["change_mode"]
            bond_change=settings["change_mode"]*(change_fn(changes[key]-mean_change,0.01,settings))
        else:#если разница между силами не становится меньше - меняем знак в алгоритме
            cur_c_m=-settings["change_mode"]
            bond_change=cur_c_m*(change_fn(changes[key]-mean_change,0.01,settings))

        res_bond=bond_len+bond_change
        if settings["print_output"]:
            print(f'\033[93m{key}\033[00m\tchange {"{:14.10f}".format(bond_change)}, res {"{:14.10f}".format(res_bond)}')
        control_strs.append(f"    distance: {key}, {res_bond}\n")
        lens[key]=res_bond

        if res_bond>MAX_BOND or res_bond<MIN_BOND:
            settings=dict(step=0,prev_dc=100,change_mode=0, preferred_change_mode=0, pass_turns=0, DoC_cutoff=0,bond_reach_critical_len=True, nAtoms=settings["nAtoms"],nBonds=settings["nBonds"], print_output=settings["print_output"])
    
            #если в результате поиска ПС связь порвалась или замкнулась, то  результат сбрасывается, а в начальной геометрии соответствующая связь немного (на 0,02) растягивается или укорачивается, чтобы притяжение или отталкивание было меньше
            settings["bond_reach_critical_len"]=True

            if bond[0]<bond[1]:
                key=f"{bond[0]}, {bond[1]}"
            else:
                key=f"{bond[1]}, {bond[0]}"

            if res_bond>MAX_BOND:
                init_bonds[key]-=0.1
            else:
                init_bonds[key]+=0.1
            print(init_bonds)
    if settings["pass_turns"]>0:
        if abs(change_if_less_DoC)==0.000001:
            settings["DoC_cutoff"]*=0.99
        else:
            settings["DoC_cutoff"]*=0.997
        settings["pass_turns"]-=1

    settings["prev_dc"]=div_of_changes
    if settings["change_mode"]*cur_c_m<0:
        settings["change_mode"]=cur_c_m
    return max(abs(min_change),max_change)

def read_bonds(rpath:str):
    search_bonds=[]
    with open(os.path.join(rpath,"bonds_to_search"),"r") as bonds:
        #print(rpath)
        chrg=int(bonds.readline())
        solvent=bonds.readline().split()[0]
        line=bonds.readline()
        #print(line)
        while line != "":
            if line.startswith("b"):#это связь
                #print(line)
                line_split=line.split()
                if line_split[0]=='b':
                    search_bonds.append([int(line_split[1]), int(line_split[2]), int(line_split[3])])
            line=bonds.readline()
    return chrg,solvent,search_bonds

def find_reac_type_by_graph_and_phases__and__measure_init_bonds(settings:dict, search_bonds:list, xyzs_strs:list):#именно то, что написано на упаковке, более короткого, но осмысленного названия я придумать не смог
    init_bonds={}
    Reag_graph=nx.Graph()

    first_phase=0
    reac_type=2#как Дильс-Альдер
    for bond in search_bonds:
        if first_phase==0:
            first_phase=bond[2]
        if bond[0]<bond[1]:
            key= f"{bond[0]}, {bond[1]}" 
        else:
            key= f"{bond[1]}, {bond[0]}"  
        init_bonds[key]=vec_len(extract_AB_dir(xyzs_strs,bond[0],bond[1])) 

        if bond[0] not in Reag_graph.nodes():
            Reag_graph.add_node(bond[0])
        if bond[1] not in Reag_graph.nodes():
            Reag_graph.add_node(bond[1])
        if first_phase*bond[2]<0:#как sn2
            reac_type=1
        Reag_graph.add_edge(bond[0],bond[1],phase=bond[2])
    if settings["print_output"]:
        print(reac_type)

    if reac_type==2:    
        for node in Reag_graph.nodes():
            if len(list(Reag_graph.neighbors(node)))!=1:
                if settings["print_output"]:
                    print("strange reaction") 
                exit()
        settings["preferred_change_mode"]=1
    elif reac_type==1:    
        for node in Reag_graph.nodes():
            connected_nodes=list(Reag_graph.neighbors(node))
            if len(connected_nodes)>1:
                is_1=False
                is_rev=False
                for node_2 in connected_nodes:

                    if Reag_graph.get_edge_data(node, node_2)["phase"]==1:
                        is_1=True
                    elif Reag_graph.get_edge_data(node, node_2)["phase"]==-1:
                        is_rev=True
                if not (is_1 and is_rev):
                    print("strange reaction")
                    exit()
        settings["preferred_change_mode"]=-1
    settings["change_mode"]=settings["preferred_change_mode"]
    if settings["print_output"]:
        print(init_bonds)
    return init_bonds

def log(str:str,logname:str):
    with open(logname,"a") as file:
        file.write(str)

def check_tresholds_converged(proj_len:float, optimized_cap:float, ratio:float, grad_strs:list,xyzs_strs:list, search_bonds:list,settings:dict):
    trashold_template=lambda name,cur,target,conv:f'{name} trashold {"{:15.7f}".format(cur)} of {"{:15.7f}".format(target)}: \033{"[92m" if conv else "[91mnot "}converged\033[00m'
    
    converged=True
    if optimized_cap!=0:
        if proj_len>optimized_cap:
            if settings["print_output"]:
                print(trashold_template("force",proj_len,optimized_cap,False))
            converged=False
        elif settings["print_output"]:
            print(trashold_template("force",proj_len,optimized_cap,True))
            
    if ratio!=0:
        mean_not_opt=adaptive_opt_cap(grad_strs,xyzs_strs, search_bonds,settings)
        cur_ratio=proj_len/mean_not_opt
        if cur_ratio>ratio:
            if settings["print_output"]:
                print(trashold_template("ratio",cur_ratio,ratio,False))
            converged=False
        elif settings["print_output"]:
            print(trashold_template("ratio",cur_ratio,ratio,True))
            
    return converged

def read_file(rpath:str,file_name:str):
    file_strs=[]
    with open(os.path.join(rpath,file_name),"r") as file:
        line=file.readline()
        while line!="":
            file_strs.append(line)
            line=file.readline()
    return file_strs

    #------run------#
def find_TS(rpath:str, xyz_name:str,optimized_cap:float=0, ratio:float=0, maxstep:int=7000, print_output:bool=True):
    if optimized_cap==0 and ratio==0:
        print("please, enter optimized_cap or (and) ratio")
        return

    settings=dict(step=0,prev_dc=100,change_mode=0, preferred_change_mode=0, pass_turns=0, DoC_cutoff=0,bond_reach_critical_len=True, nAtoms=0,nBonds=0, print_output=print_output)
    
    with open(os.path.join(rpath,"way_log.txt"),"w+") as file:
        pass
    initial_cwd = os.getcwd()
    os.chdir(rpath)

    logname=os.path.join(initial_cwd,"grad_log")

    with open(logname,"w+") as file:
        file.write( (datetime.datetime.now()).strftime("%Y m%m d%d %H:%M:%S\n") )

    search_bonds=[] # список связей, для которых ищется ПС. [A,B,phase]
    if settings["print_output"]:
        print("reading inputs\n")

    chrg,solvent,search_bonds= read_bonds(rpath)
    settings["nBonds"]=len(search_bonds)
    log_xyz(rpath,[],"new") 
    #log_forces(rpath,"","new")

    #Все длины, над которыми производятся операции - в ангстремах
    init_bonds={}#["a, b"] (при этом a<b - номера атомов) - начальные длины связей. При разваливании ПС изменяются так, чтобы меньше разваливались (растягиваются при критичном сжатии, сжимаются при критичном растяжении..)
    lens={}#["a, b"] текущие длины связей (к которым применяется изменение длины по градиенту)

    xyzs_strs=[]#строки координат атомов (вида "A X Y Z\n", A - символ элемента, X,Y,Z - его координаты)
    xyzs_strs=read_file(rpath,xyz_name)
    settings["nAtoms"]=len(xyzs_strs)-2
    init_bonds=find_reac_type_by_graph_and_phases__and__measure_init_bonds(settings, search_bonds, xyzs_strs)
    
    not_completed=True
    while not_completed:
        if settings["bond_reach_critical_len"]==True:
            #if rpath in lens.keys():
            #print(lens)
            lens.clear()
            #print("lens is clear")
            control_strs=[]

            for bond_atoms, bond_value in zip(init_bonds.keys(), init_bonds.values()):
                control_strs.append(f"    distance: {bond_atoms}, {bond_value}\n")
            
            save_control(rpath,chrg,6,control_strs)
            opt(rpath,xyz_name,solvent)

            settings["bond_reach_critical_len"]=False

        else:
            xyzs_strs=read_file(rpath,"xtbopt.xyz")#readen
            log_xyz(rpath,xyzs_strs)

            #string_curve=f'{vec_len(extract_AB_dir(xyzs_strs,13,14))} {vec_len(extract_AB_dir(xyzs_strs,18,19))} {xyzs_strs[1].split()[1]}\r\n'
            #log(string_curve, "way_log.txt")

            grad_strs=read_file(rpath,"gradient")
            
            if optimized_cap=="auto":
                optimized_cap=adaptive_opt_cap(grad_strs,xyzs_strs, search_bonds,settings)
                if settings["print_output"]:
                    print(f"because optimized cap is \"auto\", calculated optimized_cap is {optimized_cap}")

            control_strs=[]

            proj_len=move_bonds( xyzs_strs, search_bonds, lens, init_bonds, grad_strs, control_strs,  settings)

            not_completed = not check_tresholds_converged(proj_len, optimized_cap, ratio, grad_strs,xyzs_strs, search_bonds,settings)
            
            if settings["step"]>=maxstep:
                not_completed=False
                 
            save_control(rpath,chrg,6,control_strs)
            if settings["print_output"]:
                print("opt geometry with new control") 
            opt(rpath,"xtbopt.xyz",solvent)

        if not_completed:
            if settings["print_output"]:
                print(f'\nstep {settings["step"]}') 
        else:
            log(f"completed at {(datetime.datetime.now()).strftime('%Y m%m d%d %H:%M:%S')}\n",logname)
        '''
        os.chdir(path)
        iteration_not_completed=True
        while iteration_not_completed:
            log(f"(stdout) "+os.path.join(initial_cwd,"squeue_stdout\n"),logname)
            with open(os.path.join(initial_cwd,"squeue_stdout"),"w") as file:
                subprocess.call(["squeue"],stdout=file)
            with open(os.path.join(initial_cwd,"squeue_stdout"),"r") as file:
                file.readline()
                if file.readline()=="":
                    iteration_not_completed=False
            time.sleep(1)'''

        if settings["print_output"]:
            print("gradient calculation")
        
        o_grad(rpath,chrg,solvent)
        #subprocess.call(["sxtb_grad", "xtbopt.xyz","control"])
    '''
        iteration_not_completed=True
        while iteration_not_completed:
            log(f"(stdout) "+os.path.join(initial_cwd,"squeue_stdout\n"),logname)
            with open(os.path.join(initial_cwd,"squeue_stdout"),"w") as file:
                subprocess.call(["squeue"],stdout=file)
            with open(os.path.join(initial_cwd,"squeue_stdout"),"r") as file:
                file.readline()
                if file.readline()=="":
                    iteration_not_completed=False
            time.sleep(1)#xtb runs about 5 seconds so there no reasonable to check more frequent 
    '''
    os.chdir(initial_cwd)
initial_cwd=os.getcwd()
path=os.path.join(initial_cwd,"tests","rad_test")
find_TS(path,"to_opt.xyz", ratio=4, optimized_cap=0.00002, print_output=True)