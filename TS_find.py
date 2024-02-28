import numpy as np,os,subprocess,datetime,time, networkx as nx
    
class optTS:
    def __init__(self,rpath:str, xyz_name:str,optimized_cap:float=0, ratio:float=0, maxstep:int=7000, print_output:bool=True, mode="strict"):
        
        if optimized_cap==0 and ratio==0:
            print("please, enter optimized_cap or (and) ratio")
            return
        self.const_settings=dict(rpath=rpath,xyz_name=xyz_name,print_output=print_output, optimized_cap=optimized_cap,ratio=ratio,maxstep=maxstep, mode=mode,preferred_change_mode=0)
        self.settings=dict(step=0,prev_dc=100,change_mode=0, pass_turns=0, DoC_cutoff=0,bond_reach_critical_len=True)
        with open(os.path.join(self.const_settings["rpath"],"way_log.txt"),"w+") as file:
            pass
        

        self.initial_cwd = os.getcwd()
        os.chdir(self.const_settings["rpath"])
    
        self.logname=os.path.join(self.initial_cwd,"grad_log")
    
        with open(self.logname,"w+") as file:
            file.write( (datetime.datetime.now()).strftime("%Y m%m d%d %H:%M:%S\n") )
    
        self.search_bonds=[] # список связей, для которых ищется ПС. [A,B,phase]
        if self.const_settings["print_output"]:
            print("reading inputs\n")
    
        self.read_bonds()
        self.const_settings["nBonds"]=len(self.search_bonds)
        if self.const_settings["nBonds"]<3 and self.const_settings["mode"]=="autostrict":
            self.const_settings["mode"]="strict"
        else:
            self.const_settings["mode"]=""
        self.xyzs_strs=[]
        self.log_xyz("new") 
        
        #Все длины, над которыми производятся операции - в ангстремах
        self.init_bonds={}#["a, b"] (при этом a<b - номера атомов) - начальные длины связей. При разваливании ПС изменяются так, чтобы меньше разваливались (растягиваются при критичном сжатии, сжимаются при критичном растяжении..)
        self.lens={}#["a, b"] текущие длины связей (к которым применяется изменение длины по градиенту)
    
        self.xyzs_strs=[]#строки координат атомов (вида "A X Y Z\n", A - символ элемента, X,Y,Z - его координаты)
        self.xyzs_strs=self.read_file(self.const_settings["xyz_name"])
        self.const_settings["nAtoms"]=len(self.xyzs_strs)-2
        self.init_bonds=self.find_reac_type_by_graph_and_phases__and__measure_init_bonds()
        
        if self.const_settings["optimized_cap"]=="auto":
            self.const_settings["optimized_cap"]=self.mean_force()
            if self.const_settings["print_output"]:
                print(f'because optimized cap is \"auto\", calculated optimized_cap is {self.const_settings["optimized_cap"]}')

        self.not_completed=True
         
        self.proceed()

    #mathematics
    @staticmethod  
    def projection(va,vb):
        return np.multiply( np.matmul(va,vb.T)/np.matmul(vb,vb.T), vb )#a to b projection
    @staticmethod  
    def vec_len(v):
        return (np.matmul(v,v.T))**0.5
    @staticmethod               
    def sign(x:float):
        if x>0:
            return 1
        if x<0:
            return -1
        return 0
    
    def change_fn(self,length:float, cap:float):
        len_sb=self.const_settings["nBonds"]
        if len_sb>2:
            return length/(len_sb)
        len_sign=self.sign(length)
        x=len_sign*length
        change=(x**1.05)/2#(x**(0.1+1.9/(x**0.1+1)))*5
        if change>cap:
            change=cap
        
        return change*len_sign
    #~mathematics

    #xtb
    def opt(self,xyz_name:str):
        with open(os.path.join(self.const_settings["rpath"],"xtbout"),"w+") as xtbout:
            if self.const_settings["solvent"]=="vacuum":
                subprocess.call(["xtb", xyz_name, "-I", "control","--vtight","--opt"],stdout=xtbout)
            else:
                subprocess.call(["xtb", xyz_name, "-I", "control","--alpb",self.const_settings["solvent"],"--opt","--vtight"],stdout=xtbout)
    
    def o_grad(self):
        with open(os.path.join(self.const_settings["rpath"],"xtbout"),"w+") as xtbout:
            if self.const_settings["solvent"]=="vacuum":
                subprocess.call(["xtb", "xtbopt.xyz", "--chrg", str(self.const_settings["chrg"]),"--grad"],stdout=xtbout)
            else:
                subprocess.call(["xtb", "xtbopt.xyz", "--chrg", str(self.const_settings["chrg"]), "--alpb", self.const_settings["solvent"],"--grad"],stdout=xtbout)
    #~xtb   
    #rw
    def save_control(self,f_c:float):
        with open(os.path.join(self.const_settings["rpath"],"control"),"w+") as control:
            control.writelines([f'$chrg {self.const_settings["chrg"]}\n',"$constrain\n"])
            control.writelines(self.control_strs)
            control.writelines([f"    force constant = {f_c}\n","$end\n"])

    def log_xyz(self, mode=None):
        if mode=="new":
            open_prm_str="w+"
        else:
            open_prm_str="a"
        with open(os.path.join(self.const_settings["rpath"],"TS_search_log.xyz"),open_prm_str) as log:
            log.writelines(self.xyzs_strs)
    
    @staticmethod
    def log(str:str,logname:str):
        with open(logname,"a") as file:
            file.write(str)
    
    def read_file(self,file_name:str):
        file_strs=[]
        with open(os.path.join(self.const_settings["rpath"],file_name),"r") as file:
            line=file.readline()
            while line!="":
                file_strs.append(line)
                line=file.readline()
        return file_strs
    #~rw
    #vector extract
    def extract_AB_dir(self, num_A:int,num_B:int):
        vec_A=self.xyzs_strs[num_A+1].split()[1:]
        for num,coord in enumerate(vec_A):
            vec_A[num]=float(coord)
    
        vec_B=self.xyzs_strs[num_B+1].split()[1:]
        for num,coord in enumerate(vec_B):
            vec_B[num]=float(coord)
    
        res=np.subtract(vec_B,vec_A)
        return res
    
    def extractGradient(self,line_num:int):
        gradline=self.grad_strs[line_num]
        gradstr_arr=gradline.split()
        gradarr=[]
        for i in range(len(gradstr_arr)):
            gradarr.append(float(gradstr_arr[i]))
        return np.array(gradarr)
    #~vector extract
    #init fns
    def read_bonds(self):
        self.search_bonds=[]
        with open(os.path.join(self.const_settings["rpath"],"bonds_to_search"),"r") as bonds:
            #print(rpath)
            self.const_settings["chrg"]=int(bonds.readline())
            self.const_settings["solvent"]=bonds.readline().split()[0]
            line=bonds.readline()
            #print(line)
            while line != "":
                if line.startswith("b"):#это связь
                    #print(line)
                    line_split=line.split()
                    if line_split[0]=='b':
                        self.search_bonds.append([int(line_split[1]), int(line_split[2]), int(line_split[3])])
                line=bonds.readline()
    
    def find_reac_type_by_graph_and_phases__and__measure_init_bonds(self):#именно то, что написано на упаковке, более короткого, но осмысленного названия я придумать не смог
        init_bonds={}
        Reag_graph=nx.Graph()
    
        first_phase=0
        reac_type=2#как Дильс-Альдер
        for bond in self.search_bonds:
            if first_phase==0:
                first_phase=bond[2]
            if bond[0]<bond[1]:
                key= f"{bond[0]}, {bond[1]}" 
            else:
                key= f"{bond[1]}, {bond[0]}"  
            init_bonds[key]=self.vec_len(self.extract_AB_dir(bond[0],bond[1])) 
    
            if bond[0] not in Reag_graph.nodes():
                Reag_graph.add_node(bond[0])
            if bond[1] not in Reag_graph.nodes():
                Reag_graph.add_node(bond[1])
            if first_phase*bond[2]<0:#как sn2
                reac_type=1
            Reag_graph.add_edge(bond[0],bond[1],phase=bond[2])
        if self.const_settings["print_output"]:
            print(reac_type)
    
        if reac_type==2:#как Дильс-Альдер   
            for node in Reag_graph.nodes():
                if len(list(Reag_graph.neighbors(node)))!=1:
                    if self.const_settings["print_output"]:
                        print("strange reaction") 
                    exit()
            self.const_settings["preferred_change_mode"]=1
        elif reac_type==1:#как sn2 
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
            self.const_settings["preferred_change_mode"]=-1
            
        self.settings["change_mode"]=self.const_settings["preferred_change_mode"]
        if self.const_settings["print_output"]:
            print(init_bonds)
        return init_bonds
    #~init fns

    #main loop fns
    def proceed(self):
    
        while self.not_completed:
            self.control_strs=[]
            if self.settings["bond_reach_critical_len"]==True:
                self.lens.clear()
                print("lens is clear")
                
    
                for bond_atoms, bond_value in zip(self.init_bonds.keys(), self.init_bonds.values()):
                    self.control_strs.append(f"    distance: {bond_atoms}, {bond_value}\n")
                
                self.save_control(6)
                self.opt(self.const_settings["xyz_name"])
    
                self.settings["bond_reach_critical_len"]=False
    
            else:
                self.xyzs_strs=self.read_file("xtbopt.xyz")#readen
                self.log_xyz(self.xyzs_strs)
    
                if self.const_settings["nBonds"]>1:
                    string_curve=f'{self.vec_len(self.extract_AB_dir(self.search_bonds[0][0],self.search_bonds[0][1]))} {self.vec_len(self.extract_AB_dir(self.search_bonds[1][0],self.search_bonds[1][1]))} {self.xyzs_strs[1].split()[1]}\r\n'
                    self.log(string_curve, "way_log.txt")
    
                self.grad_strs=self.read_file("gradient")
    
                proj_len=self.move_bonds()
    
                self.not_completed = not self.check_tresholds_converged(proj_len)
                
                if self.settings["step"]>=self.const_settings["maxstep"]:
                    if self.const_settings["print_output"]:
                        print("\033[91mnot optimized but reached maximum number of steps\033[00m") 
                    self.not_completed=False
                     
                self.save_control(6)
                if self.const_settings["print_output"]:
                    print("opt geometry with new control") 
                self.opt("xtbopt.xyz")
    
            if self.not_completed:
                if self.const_settings["print_output"]:
                    print(f'\nstep {self.settings["step"]}') 
            else:
                self.log(f"completed at {(datetime.datetime.now()).strftime('%Y m%m d%d %H:%M:%S')}\n",self.logname)
    
            if self.const_settings["print_output"]:
                print("gradient calculation")
            
            self.o_grad()
        os.chdir(self.initial_cwd)

    def move_bonds(self):
        MIN_BOND=0.8
        MAX_BOND=3.5
        sum_changes=0
        min_change=1000
        max_change=-1000
        changes={}
        self.settings["step"]+=1
        cur_c_m=self.settings["change_mode"]
        for bond in self.search_bonds:#найдём градиент (желание растянуться) вдоль каждой связи
            num_A=bond[0]
            num_B=bond[1]
            key=f"{bond[0]}, {bond[1]}"
    
            grad_A=self.extractGradient(num_A+self.const_settings["nAtoms"]+1)
            grad_B=self.extractGradient(num_B+self.const_settings["nAtoms"]+1)
            AB_dir=self.extract_AB_dir(num_A,num_B)
            summ_grad=np.subtract(grad_B,grad_A)
    
            if not key in self.lens.keys():
                if self.const_settings["print_output"]:
                    print(f"key \"{key}\" not in lens.keys()") 
                self.lens[key] = self.vec_len(AB_dir)
    
            s_g_proj=self.projection(summ_grad, AB_dir)
            proj_len=self.vec_len(s_g_proj)
    
            s_g_p_sign=self.sign(s_g_proj[0]*AB_dir[0])
    
            changes[key]=s_g_p_sign*proj_len#удлиннение (если отрицательно - укорочение) связи
            max_change=max(changes[key], max_change)
            min_change=min(changes[key], min_change)
            sum_changes+=changes[key]
        mean_change=sum_changes/self.const_settings["nBonds"]#среднее "желание растянуться"
        div_of_changes=max_change-min_change#отклонение 
    
        if self.settings["DoC_cutoff"]==0:
            self.settings["DoC_cutoff"]=min(div_of_changes*1.05, 0.01)
        if self.const_settings["print_output"]:
            print(f'div of forces: \033[01mcur\033[00m {"{:.8f}".format(div_of_changes)}    \033[01mprev\033[00m {"{:.8f}".format(self.settings["prev_dc"])} (\033{"[92mless" if div_of_changes<self.settings["prev_dc"] else "[091mhigher"}\033[00m)') 
            print(f'change mode {self.settings["change_mode"]}') 
    
        change_if_less_DoC=self.const_settings["preferred_change_mode"]*(self.change_fn(mean_change,0.015))
        if abs(change_if_less_DoC)<0.000001:#предел разрешения control, а вблизи ПС тут получается очень маленькое значение (10^-8)
            change_if_less_DoC=self.sign(change_if_less_DoC)*0.000001
        if self.const_settings["print_output"]:
            print(f'\033[01m\033[93mbonds:\033[00m')
        for bond in self.search_bonds:
            key=f"{bond[0]}, {bond[1]}"
            bond_len=self.lens[key]
            
            if div_of_changes<self.settings["DoC_cutoff"] or self.const_settings["nBonds"]==1: # нужно "отпущение" связей
                bond_change=change_if_less_DoC
                self.settings["change_mode"]=self.const_settings["preferred_change_mode"]
                self.settings["pass_turns"]=2#1 снимется уже в конце этого хода, так что 2
            elif div_of_changes<self.settings["prev_dc"] or self.settings["pass_turns"]>0 or self.const_settings["mode"]=="strict":#нужно выравнивание по силам
                cur_c_m=self.settings["change_mode"]
                bond_change=-self.settings["change_mode"]*(self.change_fn(changes[key]-mean_change,0.01))
            else:#если разница между силами не становится меньше - меняем знак в алгоритме
                cur_c_m=-self.settings["change_mode"]
                bond_change=-cur_c_m*(self.change_fn(changes[key]-mean_change,0.01))
    
            res_bond=bond_len+bond_change
            if self.const_settings["print_output"]:
                print(f'\033[93m{key}\033[00m\tchange {"{:14.10f}".format(bond_change)}, res {"{:14.10f}".format(res_bond)}')
            self.control_strs.append(f"    distance: {key}, {res_bond}\n")
            self.lens[key]=res_bond
    
            if res_bond>MAX_BOND or res_bond<MIN_BOND:
                self.settings["step"]=0
                self.settings["prev_dc"]=100
                self.settings["pass_turns"]=0
                self.settings["DoC_cutoff"]=0
                self.settings["bond_reach_critical_len"]=True  
                #если в результате поиска ПС связь порвалась или замкнулась, то  результат сбрасывается, а в начальной геометрии соответствующая связь немного (на 0,02) растягивается или укорачивается, чтобы притяжение или отталкивание было меньше
    
    
                if bond[0]<bond[1]:
                    key=f"{bond[0]}, {bond[1]}"
                else:
                    key=f"{bond[1]}, {bond[0]}"
    
                if res_bond>MAX_BOND:
                    self.init_bonds[key]-=0.1
                else:
                    self.init_bonds[key]+=0.1
                print(self.init_bonds)
        if self.settings["pass_turns"]>0:
            if self.const_settings["mode"]=="strict":
                self.settings["DoC_cutoff"]=min(div_of_changes,self.settings["DoC_cutoff"])*0.997
            elif abs(change_if_less_DoC)==0.000001:
                self.settings["DoC_cutoff"]*=0.99#Уменьшаем сильно, т.к. делаем не то, что надо скорее всего
            else:
                self.settings["DoC_cutoff"]*=0.997#Уменьшаем немного
            self.settings["pass_turns"]-=1
    
        self.settings["prev_dc"]=div_of_changes
        if self.settings["change_mode"]*cur_c_m<0:
            self.settings["change_mode"]=cur_c_m
        return max(abs(min_change),max_change)
    
    def mean_force(self):
        search_atoms=set()
        for bond in self.search_bonds:
            search_atoms.add(bond[0])
            search_atoms.add(bond[1])
        sum_forces=0
        num_forces=0
        for i in range(1,self.const_settings["nAtoms"]+1):
            if i not in search_atoms:
                vec_force=self.extractGradient(i+self.const_settings["nAtoms"]+1)
                sum_forces+=self.vec_len(vec_force)
                num_forces+=1
        return sum_forces/num_forces
    
    def check_tresholds_converged(self,proj_len:float):
        trashold_template=lambda name,cur,target,conv:f'{name} trashold {"{:15.7f}".format(cur)} of {"{:15.7f}".format(target)}: \033{"[92m" if conv else "[91mnot "}converged\033[00m'
        
        converged=True
        if self.const_settings["optimized_cap"]!=0:
            if proj_len>self.const_settings["optimized_cap"]:
                if self.const_settings["print_output"]:
                    print(trashold_template("force",proj_len,self.const_settings["optimized_cap"],False))
                converged=False
            elif self.const_settings["print_output"]:
                print(trashold_template("force",proj_len,self.const_settings["optimized_cap"],True))
                
        if self.const_settings["ratio"]!=0:
            mean_not_opt=self.mean_force()
            cur_ratio=proj_len/mean_not_opt
            if cur_ratio>self.const_settings["ratio"]:
                if self.const_settings["print_output"]:
                    print(trashold_template("ratio",cur_ratio,self.const_settings["ratio"],False))
                converged=False
            elif self.const_settings["print_output"]:
                print(trashold_template("ratio",cur_ratio,self.const_settings["ratio"],True))
                
        return converged
    #~main loop fns
#------run------#
'''
initial_cwd=os.getcwd()
path=os.path.join(initial_cwd,"tests","da_test")
optTS(path,"to_opt.xyz", ratio=8, optimized_cap=0.00004, print_output=True,mode="")
'''