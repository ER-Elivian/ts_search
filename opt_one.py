import numpy as np,os,subprocess,datetime,time, networkx as nx
projection = lambda va,vb: np.multiply( np.matmul(va,vb.T)/np.matmul(vb,vb.T), vb )#a to b projection
extractdir = lambda xyzs,num_N,num_O:np.subtract( xyzs[num_N][1:],xyzs[num_O][1:])
vec_len = lambda v:(np.matmul(v,v.T))**0.5
control_top = lambda chrg,control_strs: control_strs.extend([f"$chrg {chrg}\n","$constrain\n"]) #заполнить верхнюю часть control файла
control_end = lambda f_c,control_strs: control_strs.extend([f"    force constant = {f_c}\n","$end\n"])#заполнить нижнюю часть control файла

def find_TS(path, xyz_name,optimized_cap=0, ratio=0, maxstep=7000):
    if optimized_cap==0 and ratio==0:
        print("please, enter optimized_cap or (and) ratio")
        return
    
    settings=dict(step=0,prev_dc=100,change_mode=0, preferred_change_mode=0, pass_turns=0, DoC_cutoff=0)
    MIN_BOND=0.8
    MAX_BOND=2.8

    def sign(x):
        if x>0:
            return 1
        if x<0:
            return -1
        return 0

    def opt(rpath,xyz_name,solvent):
        with open(os.path.join(rpath,"xtbout"),"w+") as xtbout:
            if solvent=="vacuum":
                subprocess.call(["xtb", xyz_name, "-I", "control","--vtight","--opt"],stdout=xtbout)
            else:
                subprocess.call(["xtb", xyz_name, "-I", "control","--alpb",solvent,"--opt","--vtight"],stdout=xtbout)

    def o_grad(rpath,chrg,solvent):
        with open(os.path.join(rpath,"xtbout"),"w+") as xtbout:
            if solvent=="vacuum":
                subprocess.call(["xtb", "xtbopt.xyz", "--chrg", str(chrg),"--grad"],stdout=xtbout)
            else:
                subprocess.call(["xtb", "xtbopt.xyz", "--chrg", str(chrg), "--alpb", solvent,"--grad"],stdout=xtbout)
    
    def adaptive_opt_cap(grad_strs, bonds_to_search):
        search_atoms=set()
        for bond in bonds_to_search:
            search_atoms.add(bond[0])
            search_atoms.add(bond[1])
        sum_forces=0
        num_forces=0
        for i in range(1,len(xyzs_strs)-1):
            if i not in search_atoms:
                vec_force=extractGradient(i+len(xyzs_strs)-1,grad_strs)
                sum_forces+=vec_len(vec_force)
                num_forces+=1
        return sum_forces/num_forces
        
    def change_fn(length, cap):
        return length/len(search_bonds)
        '''len_sign=sign(length)
        x=len_sign*length
        change=(x**1.4)*5#(x**(0.1+1.9/(x**0.1+1)))*5

        if change>cap:
            change=cap
        return change*len_sign'''

    def log_xyz(path, xyzs_strs, mode=None):
        if mode=="new":
            open_prm_str="w+"
        else:
            open_prm_str="a"
        with open(os.path.join(path,"TS_search_log.xyz"),open_prm_str) as log:
            log.writelines(xyzs_strs)

    def log_forces(path, forces, mode=None):
        if mode=="new":
            open_prm_str="w+"
        else:
            open_prm_str="a"
        with open(os.path.join(path,"force_log.xyz"),open_prm_str) as log:
            log.write(forces)

    def extract_AB_dir(xyzs_strs, num_A,num_B):
        vec_A=xyzs_strs[num_A+1].split()[1:]
        for num,coord in enumerate(vec_A):
            vec_A[num]=float(coord)

        vec_B=xyzs_strs[num_B+1].split()[1:]
        for num,coord in enumerate(vec_B):
            vec_B[num]=float(coord)

        res=np.subtract(vec_B,vec_A)
        return res

    def extractGradient(line_num,grad_strs):
        gradline=grad_strs[line_num]
        gradstr_arr=gradline.split()
        gradarr=[]
        for i in range(len(gradstr_arr)):
            gradarr.append(float(gradstr_arr[i]))
        return np.array(gradarr)
    
    def move_bonds(rpath, xyzs_strs, search_bonds, lens, not_optimized_in_iteration, init_bonds, grad_strs, control_strs, logname):
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

            grad_A=extractGradient(num_A+len(xyzs_strs)-1,grad_strs)
            grad_B=extractGradient(num_B+len(xyzs_strs)-1,grad_strs)
            AB_dir=extract_AB_dir(xyzs_strs,num_A,num_B)
            summ_grad=np.subtract(grad_B,grad_A)

            if not key in lens.keys():
                print(f"key \"{key}\" not in lens.keys()")
                lens[key] = vec_len(AB_dir)

            s_g_proj=projection(summ_grad, AB_dir)
            proj_len=vec_len(s_g_proj)

            s_g_p_sign=sign(s_g_proj[0]*AB_dir[0])

            changes[key]=s_g_p_sign*proj_len#удлиннение (если отрицательно - укорочение) связи
            max_change=max(changes[key], max_change)
            min_change=min(changes[key], min_change)
            sum_changes+=changes[key]
        mean_change=sum_changes/len(search_bonds)#среднее "желание растянуться"
        div_of_changes=max_change-min_change#отклонение 

        if settings["DoC_cutoff"]==0:
            settings["DoC_cutoff"]=max(div_of_changes*1.001, 0.001)
        
        #print(div_of_changes)
        #print(prev_dc)

        change_if_less_DoC=settings["preferred_change_mode"]*(change_fn(mean_change,0.01))
        if abs(change_if_less_DoC)<0.000001:#предел разрешения control, а вблизи ПС тут получается очень маленькое значение (10^-8)
            change_if_less_DoC=sign(change_if_less_DoC)*0.000001
        for bond in search_bonds:
            key=f"{bond[0]}, {bond[1]}"
            bond_len=lens[key]
            
            if div_of_changes<settings["DoC_cutoff"]: # нужно "отпущение" связей
                bond_change=change_if_less_DoC
                settings["change_mode"]=-settings["preferred_change_mode"]
                settings["pass_turns"]=2#1 снимется уже в конце этого хода, так что 2
            elif div_of_changes<settings["prev_dc"] or settings["pass_turns"]>0:#нужно выравнивание по силам
                #print(f"{div_of_changes<prev_dc} {settings["change_mode"]}")
                cur_c_m=settings["change_mode"]
                bond_change=settings["change_mode"]*(change_fn(changes[key]-mean_change,0.01))
            else:#если разница между силами не становится меньше - меняем знак в алгоритме
                cur_c_m=-settings["change_mode"]
                bond_change=cur_c_m*(change_fn(changes[key]-mean_change,0.01))

            res_bond=bond_len+bond_change
            #print(f"change {bond_change}, len {bond_len}, res {res_bond}")
            control_strs.append(f"    distance: {key}, {res_bond}\n")
            lens[key]=res_bond

            if res_bond>MAX_BOND or res_bond<MIN_BOND:
            #если в результате поиска ПС связь порвалась или замкнулась, то  результат сбрасывается, а в начальной геометрии соответствующая связь немного (на 0,02) растягивается или укорачивается, чтобы притяжение или отталкивание было меньше
                not_optimized_in_iteration[rpath]=True

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
                settings["DoC_cutoff"]*=0.995
            settings["pass_turns"]-=1

        settings["prev_dc"]=div_of_changes
        if settings["change_mode"]*cur_c_m<0:
            settings["change_mode"]=cur_c_m
        return max(abs(min_change),max_change)

    def read_bonds(rpath,search_bonds):
        with open(os.path.join(rpath,"bonds_to_search"),"r") as bonds:
            #print(rpath)
            chrg=int(bonds.readline())
            solvent=bonds.readline().split()[0]
            line=bonds.readline()
            #print(line)
            while line != "":
                if line.startswith("b"):#все остальные реакции обозначим просто связями без отступа
                    #print(line)
                    line_split=line.split()
                    if line_split[0]=='b':
                        search_bonds.append([int(line_split[1]), int(line_split[2]), int(line_split[3])])
                line=bonds.readline()
        return chrg,solvent

    def log(str,logname):
        with open(logname,"a") as file:
            file.write(str)

    #------run------#
    #method="--gfn1"
    #with open("way_log.txt","w+") as file:
    #    a=1
    initial_cwd = os.getcwd()
    os.chdir(path)
    #print(path) 

    logname=os.path.join(initial_cwd,"grad_log")

    with open(logname,"w+") as file:
        file.write( (datetime.datetime.now()).strftime("%Y m%m d%d %H:%M:%S\n") )

    not_small_g_p_dict={} #not small gradient prijection dict - если [rpath] true, то нужны доальнейшие шаги
    not_optimized_in_iteration={} # если [rpath] true, если в результате наших итераций связь растянулась или сжалась - в общем, ПС развалилось и надо заново делать начальный control с небольшими изменениями, и оптимизировать с самого начала

    search_bonds=[] # список связей, для которых ищется ПС. Связь может быть вида [a,b], а может [a,b,c] - тогда считается, что это SN2
    log("reading inputs\n",logname)

    rpath=path
    not_small_g_p_dict[rpath]=True
    not_optimized_in_iteration[rpath]=True

    chrg,solvent= read_bonds(rpath,search_bonds)

    log_xyz(rpath,[],"new") 
    log_forces(rpath,"","new")

    init_bonds={}#["a, b"] (при этом a<b - номера атомов) - начальные длины связей. При разваливании ПС изменяются так, чтобы меньше разваливались (растягиваются при критичном сжатии, сжимаются при критичном растяжении..)
    lens={}#["a, b"] текущие длины связей (к которым применяется изменение длины по градиенту)

    xyzs_strs=[]#строки координат атомов (вида "A X Y Z\n", A - символ элемента, X,Y,Z - его координаты)
    with open(os.path.join(rpath,xyz_name),"r") as xyz_file:
        line=xyz_file.readline()
        while line!="":
            xyzs_strs.append(line)
            line=xyz_file.readline()

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

        if bond[1] not in Reag_graph.nodes():
            Reag_graph.add_node(bond[1])
        if bond[2] not in Reag_graph.nodes():
            Reag_graph.add_node(bond[2])
        if first_phase*bond[2]<0:#как sn2
            reac_type=1

        Reag_graph.add_edge(bond[0],bond[1],phase=bond[2])
    if reac_type==2:    
        for node in Reag_graph.nodes():
            if len(list(Reag_graph.neighbors(node)))!=1:
                #print("strange reaction")
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
    print(init_bonds)

    not_completed=True
    while not_completed:
        if not not_small_g_p_dict[rpath]:
            continue

        if not_optimized_in_iteration[rpath]:
            #if rpath in lens.keys():
            #print(lens)
            lens.clear()
            step=0
            #print("lens is clear")
            control_strs=[]

            control_top(chrg,control_strs)
            #print(init_bonds)
            for bond_atoms, bond_value in zip(init_bonds.keys(), init_bonds.values()):
                control_strs.append(f"    distance: {bond_atoms}, {bond_value}\n")
            control_end(3,control_strs)

            with open(os.path.join(rpath,"control"),"w+") as control:
                control.writelines(control_strs)
            opt(rpath,xyz_name,solvent)

            not_optimized_in_iteration[rpath]=False
            not_optimized=True

        else:
            xyzs_strs=[]#readen
            with open(os.path.join(rpath,"xtbopt.xyz"),"r") as xyz_file:
                line=xyz_file.readline()
                while line!="":
                    xyzs_strs.append(line)
                    line=xyz_file.readline()
            log_xyz(rpath,xyzs_strs)

            #string_curve=f'{vec_len(extract_AB_dir(xyzs_strs,1,11))} {vec_len(extract_AB_dir(xyzs_strs,4,12))} {xyzs_strs[1].split()[1]}\r\n'
            #log(string_curve, "way_log.txt")

            grad_strs=[]
            with open(os.path.join(rpath,"gradient"),"r") as gradient:
                line=gradient.readline()
                while line!="":
                    grad_strs.append(line)
                    line=gradient.readline()
            if optimized_cap=="auto":
                optimized_cap=adaptive_opt_cap(grad_strs, search_bonds)
                print(f"because optimized cap is not defined, calculated optimized_cap is {optimized_cap}")
                exit()
            control_strs=[]
            control_top(chrg,control_strs)

            proj_len=1000
            not_optimized=False

            proj_len=move_bonds(rpath, xyzs_strs, search_bonds, lens, not_optimized_in_iteration, init_bonds, grad_strs, control_strs, logname)

            if optimized_cap!=0:
                if proj_len>optimized_cap:
                    not_optimized=True
                else:
                    log(f"{rpath} completed\n",logname)
            if ratio!=0:
                mean_not_opt=adaptive_opt_cap(grad_strs, search_bonds)
                print(proj_len/mean_not_opt)
                if proj_len/mean_not_opt>ratio:
                    not_optimized=True
            
            if step>=maxstep:
                not_optimized=False
                 
            not_small_g_p_dict[rpath]=not_optimized
            control_end(6,control_strs)
            with open(os.path.join(rpath,"control"),"w+") as control:
                control.writelines(control_strs)

            os.chdir(rpath)
            opt(rpath,"xtbopt.xyz",solvent)

        need_next_iteration=False
        for reaction_is_not_compl in not_small_g_p_dict.values():
            if reaction_is_not_compl==True:
                need_next_iteration=True
                break
        not_completed=need_next_iteration
        if not_completed:
            log("-------next iteration-------\n",logname)
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

        log("grad opt\n",logname)

        if not not_small_g_p_dict[rpath]:
            continue
        
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
#initial_cwd=os.getcwd()
#path=os.path.join(initial_cwd,"tests","sn2_test")
#find_TS(path,"to_opt.xyz", ratio=10)