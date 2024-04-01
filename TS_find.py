import numpy as np,os,subprocess,datetime,copy
    
class optTS:
    def __init__(self,rpath:str, xyz_name:str,optimized_cap:float=0, ratio:float=0, maxstep:int=7000, print_output:bool=True, mode="strict"):
        
        if optimized_cap==0 and ratio==0:
            print("please, enter optimized_cap or (and) ratio")
            return
        self.const_settings=dict(rpath=rpath,xyz_name=xyz_name,print_output=print_output, optimized_cap=optimized_cap,ratio=ratio,maxstep=maxstep, mode=mode)
        self.settings=dict(step=0,prev_dc=100, pass_turns=0, DoC_cutoff=0,bond_reach_critical_len=True)

        self.log("",os.path.join(self.const_settings["rpath"],"log_doc"))
        self.log("",os.path.join(self.const_settings["rpath"],"log_E"))
        self.log("",os.path.join(self.const_settings["rpath"],"way_log.txt"))

        self.initial_cwd = os.getcwd()
        os.chdir(self.const_settings["rpath"])
    
        self.logname=os.path.join(self.initial_cwd,"grad_log")
    
        with open(self.logname,"w+") as file:
            file.write( (datetime.datetime.now()).strftime("%Y m%m d%d %H:%M:%S\n") )
    
        self.search_bonds=[] # список связей, для которых ищется ПС. [A,B,phase]
        self.ifprint("reading inputs\n")
    
        self.read_bonds()
        self.const_settings["nBonds"]=len(self.search_bonds)
        if self.const_settings["nBonds"]<3 and self.const_settings["mode"]=="autostrict":
            self.const_settings["mode"]="strict"
        
        self.xyzs_strs=[]
        self.log_xyz("new") 
        
        #Все длины, над которыми производятся операции - в ангстремах
        self.init_bonds={}#["a, b"] (при этом a<b - номера атомов) - начальные длины связей. При разваливании ПС изменяются так, чтобы меньше разваливались (растягиваются при критичном сжатии, сжимаются при критичном растяжении..)
        self.lens={}#["a, b"] текущие длины связей (к которым применяется изменение длины по градиенту)
    
        self.xyzs_strs=[]#строки координат атомов (вида "A X Y Z\n", A - символ элемента, X,Y,Z - его координаты)
        self.xyzs_strs=self.read_file(self.const_settings["xyz_name"])
        self.const_settings["nAtoms"]=len(self.xyzs_strs)-2
        
        
        self.change_projections={}
        self.change_projections["vn_for_change"]=[[self.const_settings["nBonds"]**-0.5 for i in range(self.const_settings["nBonds"])]]#self.make_ort111(self.const_settings["nBonds"])
        for i in range(self.const_settings["nBonds"]-1):
            self.change_projections["vn_for_change"].append([0 for j in range(self.const_settings["nBonds"])])
        self.change_projections["reliability"]=[(2 if i>0 else 100) for i in range(self.const_settings["nBonds"])]
        self.change_projections["selected_vector"]=1
        
        self.init_bonds=self.find_reac_type_by_phases__and__measure_init_bonds()

        self.init_change_projections=copy.deepcopy(self.change_projections)
        

        if self.const_settings["optimized_cap"]=="auto":
            self.const_settings["optimized_cap"]=self.mean_force()
            self.ifprint(f'because optimized cap is \"auto\", calculated optimized_cap is {self.const_settings["optimized_cap"]}')

        self.not_completed=True
        self.series=dict(serie=[], last=[]) 
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
    @staticmethod
    def vsign(v1:list,v2:list):
        cos_v1v2=np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
        if cos_v1v2>0.9:
            return 1
        if cos_v1v2<-0.9:
            return -1
        return None
        
    @staticmethod               
    def cos2v(v1:list, v2:list):
        return abs(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
    
    def change_fn(self,length:float, cap:float):
        len_sb=self.const_settings["nBonds"]
        if len_sb>2:
            return length/((len_sb)**(1+len_sb*0.25))
        len_sign=self.sign(length)
        x=len_sign*length
        change=(x**1.05)/2#(x**(0.1+1.9/(x**0.1+1)))*5
        if change>cap:
            change=cap
        
        return change*len_sign
    
    def produce_new_vector(self,changes):#процесс Грамма-Шмидта
        init=copy.deepcopy(changes)
        self.ifprint(init)
        for vector in self.change_projections["vn_for_change"]:
            if np.linalg.norm(vector)>0.1:
                init=np.subtract(init,self.projection(np.array(init),np.array(vector)))
            else:
                break
        init=np.multiply(1/np.linalg.norm(init),init)
        self.ifprint(init)
        return init
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
    def ifprint(self,to_print):
        if self.const_settings["print_output"]:
            print(to_print)

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
        with open(logname,"a" if str!="" else "w+") as file:
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
            self.const_settings["chrg"]=int(bonds.readline())
            self.const_settings["solvent"]=bonds.readline().split()[0]
            line=bonds.readline()
            while line != "":
                if line.startswith("b"):#это связь
                    line_split=line.split()
                    if line_split[0]=='b':
                        self.search_bonds.append([int(line_split[1]), int(line_split[2]), int(line_split[3])])
                line=bonds.readline()
    
    def find_reac_type_by_phases__and__measure_init_bonds(self):#именно то, что написано на упаковке, более короткого, но осмысленного названия я придумать не смог
        init_bonds={}
        reac_type=2
        phases_vec=[]
        for bond in self.search_bonds:
            phases_vec.append(bond[2])
            if bond[0]<bond[1]:
                key= f"{bond[0]}, {bond[1]}" 
            else:
                key= f"{bond[1]}, {bond[0]}"  
            init_bonds[key]=self.vec_len(self.extract_AB_dir(bond[0],bond[1])) 
        for i in range(len(phases_vec)-1):
            if phases_vec[i]*phases_vec[i+1]<0:
                reac_type=1
                break
        self.ifprint(reac_type)

        self.change_projections["signs"]=[]
        if reac_type==2:#как Дильс-Альдер   
            self.change_projections["signs"].append(1)
            for i in range(1,len(phases_vec)):
                self.change_projections["signs"].append(-1)
        elif reac_type==1:#как sn2 
            self.change_projections["vn_for_change"][1]=self.produce_new_vector(phases_vec)
            self.change_projections["signs"].append(-1)
            self.change_projections["signs"].append(1)
            for i in range(2,len(phases_vec)):
                self.change_projections["signs"].append(-1)
            
        self.ifprint(init_bonds)
        return init_bonds
    #~init fns

    #reliability
    def increase_rel(self,vec_num:int):#увеличиваем уверенность в том, что это правильное изменение по знаку
        REL_MAX=60
        REL_INIT=2
        if self.change_projections["reliability"][vec_num]<=1:
            self.change_projections["reliability"][vec_num]=REL_INIT
        else:
            self.change_projections["reliability"][vec_num]+=1
            if self.change_projections["reliability"][vec_num]>REL_MAX:
                self.change_projections["reliability"][vec_num]=REL_MAX
    def decrease_rel(self,vec_num:int):#уменьшаем уверенность в том, что это правильное изменение по знаку
        if self.change_projections["reliability"][vec_num]<=1:
            if self.const_settings["mode"]!="strict":
                self.change_projections["signs"][vec_num]=-self.change_projections["signs"][vec_num]
                self.increase_rel(vec_num)
        else:
            if self.change_projections["reliability"][vec_num] >= 12:
                self.change_projections["reliability"][vec_num]-=2
            else:    
                self.change_projections["reliability"][vec_num]-=1
            if self.const_settings["nBonds"]==1:
                self.change_projections["reliability"][vec_num]-=1
    #~reliability

    #series 
    def add_series(self,vnum,length): #True, если началась новая серия
        if self.series["last"]==[]:#первый вызов
            self.series["last"]=[vnum,length]
            return False
        
        elif self.series["last"][0]==vnum:
            self.series["last"][1]+=length
            return False
        elif self.series["last"][1]!=0:
            v_last=copy.deepcopy(self.series["last"])
            self.series["last"]=[vnum,length]
            self.series["serie"].append(v_last)
            return True
        
    def search_for_cycle(self):
        if len(self.series["serie"])<4:
            return False, []
        for rel, vec in zip(self.change_projections["reliability"],self.change_projections["vn_for_change"]):
            if rel<4 and rel!=0 and np.linalg.norm(vec)>0.1:
                return False, []
        len_ser=len(self.series["serie"])
        
        #var1
        v1_a=[0 for c in range(self.const_settings["nBonds"])]
        v2_a=[0 for c in range(self.const_settings["nBonds"])]
        v1_a[self.series["serie"][len_ser-1][0]]+=self.series["serie"][len_ser-1][1]
        v1_a[self.series["serie"][len_ser-2][0]]+=self.series["serie"][len_ser-2][1]
        
        v2_a[self.series["serie"][len_ser-3][0]]+=self.series["serie"][len_ser-3][1]
        v2_a[self.series["serie"][len_ser-4][0]]+=self.series["serie"][len_ser-4][1]

        
        cos_v1v2=self.cos2v(v1_a,v2_a)
        #print(f"cos {cos_v1v2}")
        if cos_v1v2<0.96:
            return False,[]
        
        
        '''
        DEPTH=15 #рассматриваемая глубина истории
        num_first=max(1,len_ser-DEPTH)
        mid=(num_first-len_ser)//2
        v1= [0 for c in range(self.const_settings["nBonds"])] #сумма изменеий от num_first до mid не включительно
        v2= [0 for c in range(self.const_settings["nBonds"])] #сумма изменеий от mid до конца
        for i in range(num_first,mid):
            v1[self.series["serie"][i][0]]+=self.series["serie"][i][1]
        for i in range(mid,len_ser):
            v2[self.series["serie"][i][0]]+=self.series["serie"][i][1]
        flag_find=False
        print(1)

        for i in range(num_first,len_ser):
            NUM_IN_MID_TO_SEARCH=3
            mid_i=(i+len_ser)//2
            if mid_i!=mid:#отличие не может быть больше, чем на 1
                v1[self.series["serie"][mid][0]]+=self.series["serie"][mid][1]
                v2[self.series["serie"][mid][0]]-=self.series["serie"][mid][1]
                mid+=1
            
            v1_a=copy.deepcopy(v1)
            v2_a=copy.deepcopy(v2)
            for j in range(mid+1,min(mid+NUM_IN_MID_TO_SEARCH,len_ser-6)):
                print(f"cos {abs(np.dot(v1_a, v2_a)/(np.linalg.norm(v1_a)*np.linalg.norm(v2_a)))}")
                if abs(np.dot(v1_a, v2_a)/(np.linalg.norm(v1_a)*np.linalg.norm(v2_a)))>0.94:#угол меньше 20 градусов
                     flag_find=True
                     break
                v1_a[self.series["serie"][j][0]]+=self.series["serie"][j][1]
                v2_a[self.series["serie"][j][0]]-=self.series["serie"][j][1]
            if flag_find==True:
                break
            
            v1_a=copy.deepcopy(v1_a)
            v2_a=copy.deepcopy(v2_a)
            for j in range(mid-1,max(mid-NUM_IN_MID_TO_SEARCH,i),-1):
                v1_a[self.series["serie"][j][0]]-=self.series["serie"][j][1]
                v2_a[self.series["serie"][j][0]]+=self.series["serie"][j][1]
                if abs(np.dot(v1_a, v2_a)/(np.linalg.norm(v1_a)*np.linalg.norm(v2_a)))>0.94:#угол меньше 20 градусов
                     flag_find=True
                     break
            if flag_find==True:
                break
            
            v1[self.series["serie"][i][0]]-=self.series["serie"][i][1]
        if not flag_find:
            return False, []
        '''
        v1_a_len=np.linalg.norm(v1_a)
        v2_a_len=np.linalg.norm(v2_a)
        rel_decrease=v2_a_len/v1_a_len
        if rel_decrease<0.9:
            len_res=v2_a_len*rel_decrease*(1+rel_decrease)#аппроксимируем первыми 2 членами геометрической прогрессии
        else:
            len_res=(v2_a_len+v1_a_len)*2
        v_change=[0 for i in range(self.const_settings["nBonds"])]
        vsum=np.add(v1_a,v2_a)
        for i,coord in enumerate(vsum):
            v_change=np.add(v_change,np.multiply(coord,self.change_projections["vn_for_change"][i]))
        v_change=np.multiply(len_res/np.linalg.norm(v_change),v_change)
        #print(np.array(self.series["serie"]))
        self.series["serie"]=[]
        self.series["last"]=[]
        
        return True,v_change

    #~series
        
    #main loop fns
    def reset(self):
        self.settings["prev_dc"]=100
        self.settings["pass_turns"]=0
        self.settings["DoC_cutoff"]=0
        self.settings["bond_reach_critical_len"]=True
        self.change_projections=copy.deepcopy(self.init_change_projections)

        self.series["last"]=[]
        self.series["serie"]=[]

                
    def proceed(self):
        while self.not_completed:
            self.control_strs=[]
            if self.settings["bond_reach_critical_len"]==True:
                self.lens.clear()
                self.ifprint("lens is clear")
                
    
                for bond_atoms, bond_value in zip(self.init_bonds.keys(), self.init_bonds.values()):
                    self.control_strs.append(f"    distance: {bond_atoms}, {bond_value}\n")
                
                self.save_control(6)
                self.opt(self.const_settings["xyz_name"])
    
                self.settings["bond_reach_critical_len"]=False
    
            else:
                self.xyzs_strs=self.read_file("xtbopt.xyz")#readen
                self.log_xyz(self.xyzs_strs)
                self.log(f"{self.xyzs_strs[1].split()[1]}\n",os.path.join(self.const_settings["rpath"],"log_E"))
        
                if self.const_settings["nBonds"]>1:
                    string_curve=f'{self.vec_len(self.extract_AB_dir(self.search_bonds[0][0],self.search_bonds[0][1]))} {self.vec_len(self.extract_AB_dir(self.search_bonds[1][0],self.search_bonds[1][1]))} {self.xyzs_strs[1].split()[1]}\r\n'
                    self.log(string_curve, "way_log.txt")
    
                self.grad_strs=self.read_file("gradient")
    
                proj_len=self.move_bonds()
    
                self.not_completed = not self.check_tresholds_converged(proj_len)
                
                if self.settings["step"]>=self.const_settings["maxstep"]:
                    self.ifprint("\033[91mnot optimized but reached maximum number of steps\033[00m") 
                    self.not_completed=False
                     
                self.save_control(6)
                self.ifprint("opt geometry with new control") 
                self.opt("xtbopt.xyz")
    
            if self.not_completed:
                self.ifprint(f'\nstep {self.settings["step"]}') 
            else:
                self.log(f"completed at {(datetime.datetime.now()).strftime('%Y m%m d%d %H:%M:%S')}\n",self.logname)
    
            self.ifprint("gradient calculation")
            
            self.o_grad()
        os.chdir(self.initial_cwd)

    def move_bonds(self):
        MIN_BOND=0.8
        MAX_BOND=3.5
        sum_changes=0
        min_change=1000
        max_change=-1000
        changes=[]
        self.settings["step"]+=1
        for i,bond in enumerate(self.search_bonds):#найдём градиент (желание растянуться) вдоль каждой связи
            num_A=bond[0]
            num_B=bond[1]
            key=f"{bond[0]}, {bond[1]}"
    
            grad_A=self.extractGradient(num_A+self.const_settings["nAtoms"]+1)
            grad_B=self.extractGradient(num_B+self.const_settings["nAtoms"]+1)
            AB_dir=self.extract_AB_dir(num_A,num_B)
            summ_grad=np.subtract(grad_B,grad_A)
    
            if not key in self.lens.keys():
                self.ifprint(f"key \"{key}\" not in lens.keys()") 
                self.lens[key] = self.vec_len(AB_dir)
    
            s_g_proj=self.projection(summ_grad, AB_dir)
            proj_len=self.vec_len(s_g_proj)
    
            s_g_p_sign=self.sign(s_g_proj[0]*AB_dir[0])
    
            changes.append(s_g_p_sign*proj_len)#удлиннение (если отрицательно - укорочение) связи
            max_change=max(changes[i], max_change)
            min_change=min(changes[i], min_change)
            sum_changes+=changes[i]
        mean_change=sum_changes/self.const_settings["nBonds"]#среднее "желание растянуться"
        div_of_changes=max_change-min_change#отклонение 
    
        if self.settings["DoC_cutoff"]==0:
            self.settings["DoC_cutoff"]=min(div_of_changes*1.05, 0.05)
        self.ifprint(f'div of forces: \033[01mcur\033[00m {"{:.8f}".format(div_of_changes)}    \033[01mprev\033[00m {"{:.8f}".format(self.settings["prev_dc"])} (\033{"[92mless" if div_of_changes<self.settings["prev_dc"] else "[091mhigher"}\033[00m), cutoff {"{:.8f}".format(self.settings["DoC_cutoff"])}') 
            
    
        change_if_less_DoC=10*self.change_projections["signs"][0]*(self.change_fn(mean_change,0.015))
        if abs(change_if_less_DoC)<0.000001:#предел разрешения control, а вблизи ПС тут получается очень маленькое значение (10^-8)
            change_if_less_DoC=self.sign(change_if_less_DoC)*0.000001
            self.settings["DoC_cutoff"]=div_of_changes
        
        #блок того, что надо сделать 1 раз
        cycle_res=False
        if div_of_changes<self.settings["DoC_cutoff"] or self.const_settings["nBonds"]==1:
            cond_all=True#"условие, что двигаем все" чтобы ещё раз не проверять то же самое
            self.add_series(0,change_if_less_DoC*self.const_settings["nBonds"]**0.5)

            self.settings["DoC_cutoff"]=div_of_changes
        elif div_of_changes<self.settings["prev_dc"] or self.settings["pass_turns"]>0:#нужно выравнивание по силам
            cond_all=False
            self.increase_rel(self.change_projections["selected_vector"])
        else:#если разница между силами не становится меньше - берём следующий вектор
            cond_all=False
            self.decrease_rel(self.change_projections["selected_vector"])
            if self.const_settings["nBonds"]==2:
                self.decrease_rel(self.change_projections["selected_vector"])
            self.change_projections["selected_vector"] = self.change_projections["selected_vector"]+1 if self.change_projections["selected_vector"]+1<self.const_settings["nBonds"] else 1
            
        cycle_res,bond_jump=self.search_for_cycle()
                
        if not cond_all:
            if np.linalg.norm(self.change_projections["vn_for_change"][self.change_projections["selected_vector"]])<0.1:
                self.change_projections["vn_for_change"][self.change_projections["selected_vector"]]=self.produce_new_vector(changes)
            self.ifprint(f'\nselected vector {self.change_projections["selected_vector"]} {self.change_projections["vn_for_change"][self.change_projections["selected_vector"]]}') 
            
            changes_projected=self.projection(changes,self.change_projections["vn_for_change"][self.change_projections["selected_vector"]])
            changes_len=np.linalg.norm(changes_projected)#длина
            changes_real_len=self.change_fn(changes_len,0.01)*self.change_projections["signs"][self.change_projections["selected_vector"]]*(self.change_projections["reliability"][self.change_projections["selected_vector"]])**0.5#умножаем изменение на корень из ьуверенности в нём
            changes_real=np.multiply(changes_real_len/changes_len, changes_projected)#получаем вектор реальных изменений
            #print(changes_real)
            self.add_series(self.change_projections["selected_vector"], abs(changes_real_len) * self.vsign(changes_real, self.change_projections["vn_for_change"][self.change_projections["selected_vector"]] ) )#кладём его в историю изменений

        self.ifprint(f'reliabilites:{self.change_projections["reliability"]}')
        self.ifprint(f'signs:       {self.change_projections["signs"]}')
        
        
        #~блок того, что надо сделать 1 раз

        self.ifprint(f'\033[01m\033[93mbonds:\033[00m')
        
        for i,bond in enumerate(self.search_bonds):
            key=f"{bond[0]}, {bond[1]}"
            bond_len=self.lens[key]
            #блок того, что надо сделать для каждой связи
            if cycle_res:#если надо jump
                bond_change=bond_jump[i]
                self.ifprint("jump")
            elif cond_all: # нужно "отпущение" связей
                bond_change=change_if_less_DoC
                self.settings["pass_turns"]=2#1 снимется уже в конце этого хода, так что 2
            else:
                bond_change=changes_real[i]
            
            
            #~блок того, что надо сделать для каждой связи
            
            res_bond=bond_len+bond_change
            self.ifprint(f'\033[93m{key}\033[00m\tchange {"{:14.10f}".format(bond_change)}, res {"{:14.10f}".format(res_bond)}')
            self.control_strs.append(f"    distance: {key}, {res_bond}\n")
            self.lens[key]=res_bond
    
            if res_bond>MAX_BOND or res_bond<MIN_BOND:
                self.reset()
                #если в результате поиска ПС связь порвалась или замкнулась, то  результат сбрасывается, а в начальной геометрии соответствующая связь немного (на 0,02) растягивается или укорачивается, чтобы притяжение или отталкивание было меньше
    
                if bond[0]<bond[1]:
                    key=f"{bond[0]}, {bond[1]}"
                else:
                    key=f"{bond[1]}, {bond[0]}"
    
                if res_bond>MAX_BOND:
                    self.init_bonds[key]-=0.1
                else:
                    self.init_bonds[key]+=0.1
                self.ifprint(self.init_bonds)
        if self.settings["pass_turns"]>0:
            if self.const_settings["mode"]=="strict":
                self.settings["DoC_cutoff"]=min(div_of_changes,self.settings["DoC_cutoff"])*0.997
            elif abs(change_if_less_DoC)==0.000001:
                self.settings["DoC_cutoff"]*=0.99#Уменьшаем сильно, т.к. делаем не то, что надо скорее всего
            else:
                self.settings["DoC_cutoff"]*=0.997#Уменьшаем немного
            self.settings["pass_turns"]-=1
    
        self.settings["prev_dc"]=div_of_changes
            
        self.log(f"{div_of_changes}\n",os.path.join(self.const_settings["rpath"],"log_doc"))
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
            cond=proj_len<=self.const_settings["optimized_cap"]
            self.ifprint(trashold_template("force",proj_len,self.const_settings["optimized_cap"],cond))
            converged &= cond
                
        if self.const_settings["ratio"]!=0:
            mean_not_opt=self.mean_force()
            cur_ratio=proj_len/mean_not_opt
            cond=cur_ratio<=self.const_settings["ratio"]
            self.ifprint(trashold_template("ratio",cur_ratio,self.const_settings["ratio"],cond))
            converged &= cond
               
        return converged
    #~main loop fns
#------run------#
if __name__ == "__main__":
    initial_cwd=os.getcwd()
    path=os.path.join(initial_cwd,"tests","apw_test")
    optTS(path,"to_opt.xyz", ratio=8, optimized_cap=0.00004, print_output=True,mode="strict", maxstep=10e100)
