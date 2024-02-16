import numpy as np,os,subprocess,datetime,time
projection = lambda va,vb: np.multiply( np.matmul(va,vb.T)/np.matmul(vb,vb.T), vb )#a to b projection
extractdir = lambda xyzs,num_N,num_O:np.subtract( xyzs[num_N][1:],xyzs[num_O][1:])
vec_len = lambda v:(np.matmul(v,v.T))**0.5
control_top = lambda chrg,control_strs: control_strs.extend([f"$chrg {chrg}\n","$constrain\n"]) #заполнить верхнюю часть control файла
control_end = lambda f_c,control_strs: control_strs.extend([f"    force constant = {f_c}\n","$end\n"])#заполнить нижнюю часть control файла

min_bond=1.75
max_bond=2.65
def sign(x):
    if x>0:
        return 1
    if x<0:
        return -1
    return 0

def change_fn(length, cap):
    change=(length**1.2)*10
    if change<cap:
        change=cap
    return change

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

def move_bond_3(rpath,xyzs_strs,num_A, num_center, num_B, lens, not_optimized_in_iteration,init_bonds, grad_strs, control_strs, logname): #SN2 реакция
    grad_center=extractGradient(num_center+len(xyzs_strs)-1,grad_strs) #градиент цетрального атома
    AB_dir=extract_AB_dir(xyzs_strs,num_A,num_B) #вектор A->B 
    print(AB_dir)
    print(grad_center)
    g_c_proj=projection(grad_center,AB_dir) #проекция градиента центра на  A->B
    g_c_p_len=vec_len(g_c_proj) #знак этой проекции
    g_c_len=vec_len(grad_center)
    Ac_dir=extract_AB_dir(xyzs_strs,num_A,num_center) #вектор A->center  
    Bc_dir=extract_AB_dir(xyzs_strs,num_B,num_center) #вектор B->center

    #Пояснение - констрейны в control не выдерживаются абсолютно точно. Т.е. если ставится констрейн на 2 ангстрема, то результирующее значение будет примерно 2+-0.00001. Поэому нужно запоминать длины связей, какие они должны быть, а не пересчитывать их из полученного xyz
    if not rpath in lens.keys():#если длины связи неизвестны - создадим такой элемент в словаре и занесём длины, которые прочитались из xyz
        lens[rpath]={}
    if not f"{num_A}_{num_center}_{num_B}" in lens[rpath].keys():
        lens[rpath][f"{num_A}_{num_center}_{num_B}"] = (vec_len(Ac_dir), vec_len(Bc_dir))
    Ac_len, Bc_len=lens[rpath][f"{num_A}_{num_center}_{num_B}"] #Считаем длинами связей длины, находящиеся в словаре
    
    log(f"{rpath} ({num_A}, {num_center}, {num_B})\n atom_A={xyzs_strs[num_A]}\natom_center={xyzs_strs[num_center]}\natom_B={xyzs_strs[num_B]}\n",logname)
    #Пояснение - вначале выставляется центральный атом так, чтобы на него не действовала сила вдоль отрезка, соединяющего крайние атомы (движение против силы)
    # затем выставляются крайние атомы (движение по силе)
    
    g_c_p_sign=sign(g_c_proj[0]*AB_dir[0])

    grad_A=extractGradient(num_A+len(xyzs_strs)-1,grad_strs) #градиант на атоме num_A
    grad_B=extractGradient(num_B+len(xyzs_strs)-1,grad_strs) #и на num_B
    g_A_proj=projection(grad_A,Ac_dir) #проекции на направления связей A с центральным атомом и B с ним же
    g_B_proj=projection(grad_B,Bc_dir)
    g_A_p_len=vec_len(g_A_proj) #их длины
    g_B_p_len=vec_len(g_B_proj)

    g_A_p_sign=sign(g_A_proj[0]*Ac_dir[0])
    g_B_p_sign=sign(g_B_proj[0]*Bc_dir[0])
    
    log_forces(rpath,f"{g_A_p_len*g_A_p_sign}, {g_c_p_len*g_c_p_sign}, {g_B_p_len*g_B_p_sign}\n")
    print(g_c_p_len/g_c_len)
    grad=max(g_A_p_len,g_B_p_len,g_c_p_len)# максимальный граиент (нужен ля определения - надо ли продолжать оптимизацию)
        
    if  g_c_p_len/g_c_len<0.1 or g_c_p_len<0.00001:
        
        
        
        # Пояснение - для нахождения смещения атома нужна функция, обладающая следующими свойствами: 
        # 1) монотонно возрастающая
        # 2) f(0)=0 
        # 3) имеющая уменьшающуюся при приближении к 0 производную - чтобы чем ближе сила к нулю, тем слабее было изменение и, соответственно для минимизации дрожания
        # 4) стремящаяся к константе при больших значениях - чтобы даже при большой силе атом сильно не улетал
        # Вариант "y=x^1.3, но если x^1.3>0.00015 y=0.00015" оказался работоспособным. Оба коэффициента могут изменяться в довольно широких пределах
        res_Ac=Ac_len+g_A_p_sign*change_fn(g_A_p_len, 0.07)
        res_Bc=Bc_len+g_B_p_sign*change_fn(g_B_p_len, 0.07)

        log(f"move A and B\nchange Ac {g_A_p_sign*g_A_p_len}\nchange Bc {g_B_p_sign*g_B_p_len}\n",logname)
        print(f"move A and B\nchange Ac {g_A_p_sign*g_A_p_len}\nchange Bc {g_B_p_sign*g_B_p_len}")
        #else:
        #    res_Ac=Ac_len
        #    res_Bc=Bc_len
    else:
        log(f"move center, Ac change {g_c_p_sign*g_c_p_len}\n",logname)
        print(f"move center, Ac change {g_c_p_sign*g_c_p_len}")
        change=change_fn(g_c_p_len,0.05)
        res_Ac=Ac_len+g_c_p_sign*change
        res_Bc=Bc_len-g_c_p_sign*change
        
    print(f"res_Ac {res_Ac}, res_Bc {res_Bc}")
    control_strs.append(f"    distance: {num_A}, {num_center}, {res_Ac}\n")
    control_strs.append(f"    distance: {num_B}, {num_center}, {res_Bc}\n")
    lens[rpath][f"{num_A}_{num_center}_{num_B}"]=(res_Ac,res_Bc)
    
    if res_Ac>max_bond or res_Ac<min_bond or res_Bc>max_bond or res_Bc<min_bond:
    #Пояснение- если в результате поиска ПС какая-то связь порвалась или замкнулась, то  результат сбрасывается, а в начальной геометрии соответствующая связь немного (на 0,02) растягивается или укорачивается, чтобы притяжение или отталкивание было меньше
        not_optimized_in_iteration[rpath]=True
        
        if num_A<num_center:
            key1=f"{num_A}, {num_center}"
        else:
            key1=f"{num_center}, {num_A}"
        if num_B<num_center:
            key2=f"{num_B}, {num_center}"
        else:
            key2=f"{num_center}, {num_B}"

        if res_Ac>max_bond:
            init_bonds[rpath][key1]-=0.1
        elif res_Ac<min_bond:
            init_bonds[rpath][key1]+=0.1
        if res_Bc>max_bond:
            init_bonds[rpath][key2]-=0.1
        elif res_Bc<min_bond:
            init_bonds[rpath][key2]+=0.1

    return grad
        
def move_bond_2(rpath,xyzs_strs, num_A, num_B, lens, not_optimized_in_iteration,init_bonds, grad_strs, control_strs, logname):
    grad_A=extractGradient(num_A+len(xyzs_strs)-1,grad_strs)
    grad_B=extractGradient(num_B+len(xyzs_strs)-1,grad_strs)
    AB_dir=extract_AB_dir(xyzs_strs,num_A,num_B)
    summ_grad=np.subtract(grad_B,grad_A)
    
    s_g_proj=projection(summ_grad, AB_dir)
    proj_len=vec_len(s_g_proj)

    s_g_p_sign=sign(s_g_proj[0]*AB_dir[0])
    
    bond_change=s_g_p_sign*change_fn(proj_len,0.07)
    if not f"{rpath}_{num_A}_{num_B}" in lens.keys():
        lens[f"{rpath}_{num_A}_{num_B}"] = vec_len(AB_dir)
    bond_len=lens[f"{rpath}_{num_A}_{num_B}"]
    res_bond=bond_len-bond_change
    print(f"change {bond_change}, len {bond_len}, res {res_bond}")
    log(f"{rpath} ({num_A}, {num_B})\n atom_A={xyzs_strs[num_A]}\natom_B={xyzs_strs[num_B]}\ngrad_sum={summ_grad}\nproj={s_g_proj}\nproj_len={proj_len}\n\n",logname)
    control_strs.append(f"    distance: {num_A}, {num_B}, {res_bond}\n")
    lens[f"{rpath}_{num_A}_{num_B}"]=res_bond
    return proj_len

def log(str,logname):
    with open(logname,"a") as file:
        file.write(str)



#------run------#

initial_cwd = os.getcwd()
path=os.path.join(initial_cwd,"reactionsamples")

print(path)

logname=os.path.join(initial_cwd,"grad_log")
with open(logname,"w+") as file:
    file.write( (datetime.datetime.now()).strftime("%Y m%m d%d %H:%M:%S\n") )

reactions=os.listdir(path)[:1]#список папок реакций. rpath - полный path папки реакции

not_small_g_p_dict={} #not small gradient prijection dict - если [rpath] true, то нужны доальнейшие шаги
not_optimized_in_iteration={} # если [rpath] true, если в результате наших итераций связь растянулась или сжалась - в общем, ПС развалилось и надо заново делать начальный control с небольшими изменениями, и оптимизировать с самого начала

search_bonds={} #[rpath] список связей, для которых ищется ПС. Связь может быть вида [a,b], а может [a,b,c] - тогда считается, что это SN2
log("reading inputs\n",logname)
for reaction in reactions:
    rpath=os.path.join(path,reaction)
    not_small_g_p_dict[rpath]=True
    not_optimized_in_iteration[rpath]=True
    
    search_bonds[rpath]=[]
    with open(os.path.join(rpath,"bonds_to_search"),"r") as bonds:
        print(rpath)
        line=bonds.readline()
        print(line)
        while line != "":
            line_split=line.split()
            len_l_s=len(line_split)
            if len_l_s == 2:
                search_bonds[rpath].append([int(line_split[0]), int(line_split[1])])
            elif len_l_s == 3:
                search_bonds[rpath].append([int(line_split[0]), int(line_split[1]), int(line_split[2])])
            else:
                print(line_split)
            line=bonds.readline()

    log_xyz(rpath,[],"new") 
    log_forces(rpath,"","new")


init_bonds={}#[rpath]["a, b"] (при этом a<b - номера атомов) - начальные длины связей. При разваливании ПС изменяются так, чтобы меньше разваливались (растягиваются при критичном сжатии, сжимаются при критичном растяжении..)
for reaction in reactions:
    rpath=os.path.join(path,reaction)
    xyzs_strs=[]#строки координат атомов (вида "A X Y Z\n", A - символ элемента, X,Y,Z - его координаты)
    with open(os.path.join(rpath,"xtbopt.xyz"),"r") as xyz_file:
        line=xyz_file.readline()
        while line!="":
            xyzs_strs.append(line)
            line=xyz_file.readline()
            
    for bond in search_bonds[rpath]:
        init_bonds[rpath]={}
        if len(bond)==2:
            if bond[0]<bond[1]:
                init_bonds[rpath][ f"{bond[0]}, {bond[1]}" ]=vec_len(extract_AB_dir(xyzs_strs,bond[0],bond[1]))     
            else:
                init_bonds[rpath][ f"{bond[1]}, {bond[0]}" ]=vec_len(extract_AB_dir(xyzs_strs,bond[0],bond[1])) 
        elif len(bond)==3:
            if bond[0]<bond[1]:
                init_bonds[rpath][ f"{bond[0]}, {bond[1]}" ]=vec_len(extract_AB_dir(xyzs_strs,bond[0],bond[1]))     
            else:
                init_bonds[rpath][ f"{bond[1]}, {bond[0]}" ]=vec_len(extract_AB_dir(xyzs_strs,bond[0],bond[1])) 
            if bond[1]<bond[2]:
                init_bonds[rpath][ f"{bond[1]}, {bond[2]}" ]=vec_len(extract_AB_dir(xyzs_strs,bond[1],bond[2])) 
            else:
                init_bonds[rpath][ f"{bond[2]}, {bond[1]}" ]=vec_len(extract_AB_dir(xyzs_strs,bond[1],bond[2])) 
print(init_bonds)

lens={}#[rpath]["{num_A}_{num_center}_{num_B}"] текущие длины связей (к которым применяется изменение длины по градиенту)
not_completed=True
while not_completed:
    for reaction in reactions:
        rpath=os.path.join(path,reaction)
        if not not_small_g_p_dict[rpath]:
            continue
            
        if not_optimized_in_iteration[rpath]:
            if rpath in lens.keys():
                del lens[rpath]
            control_strs=[]
            
            control_top(1,control_strs)
            for bond_atoms, bond_value in zip(init_bonds[rpath].keys(), init_bonds[rpath].values()):
                control_strs.append(f"    distance: {bond_atoms}, {bond_value}\n")
            control_end(6,control_strs)
            with open(os.path.join(rpath,"control"),"w+") as control:
                control.writelines(control_strs)
            
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
            
            grad_strs=[]
            with open(os.path.join(rpath,"gradient"),"r") as gradient:
                line=gradient.readline()
                while line!="":
                    grad_strs.append(line)
                    line=gradient.readline()
            
            control_strs=[]
            control_top(1,control_strs)
            
            proj_len=1000
            not_optimized=False
            print(search_bonds[rpath])
            for bond in search_bonds[rpath]:
                len_bond=len(bond)
                if len_bond == 2:
                    proj_len=move_bond_2(rpath,xyzs_strs,bond[0],bond[1], lens, grad_strs,control_strs,logname)
                elif len_bond == 3:
                    proj_len=move_bond_3(rpath,xyzs_strs,bond[0],bond[1],bond[2],lens, not_optimized_in_iteration, init_bonds, grad_strs, control_strs, logname)
                else:                
                    print("strange bond:")
                print(bond)
                if proj_len>0.00002:
                    not_optimized=True
                else:
                    log(f"{rpath} completed\n",logname)
            
            not_small_g_p_dict[rpath]=not_optimized
            control_end(6,control_strs)
            with open(os.path.join(rpath,"control"),"w+") as control:
                control.writelines(control_strs)
        
        os.chdir(rpath)
        with open(os.path.join(rpath,"xtbout"),"w+") as xtbout:
            subprocess.call(["xtb", "xtbopt.xyz", "-I", "control","--alpb","water","--acc","0.5","--opt","--tight"],stdout=xtbout)
            subprocess.call(["xtb", "xtbopt.xyz", "--chrg", "1","--alpb","water","--acc","0.5","--grad"],stdout=xtbout)
            #subprocess.call(["sxtb_opt", "xtbopt.xyz","control"])
            
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
        time.sleep(1)
        
    log("grad opt\n",logname)
    for reaction in reactions:
        rpath=os.path.join(path,reaction)
        if not not_small_g_p_dict[rpath]:
            continue
        os.chdir(rpath)
        subprocess.call(["sxtb_grad", "xtbopt.xyz","control"])

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
        time.sleep(1)#xtb runs about 5 seconds so there no reasonable to check more frequent 
        '''
