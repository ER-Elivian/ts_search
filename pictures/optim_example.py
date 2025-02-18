import numpy as np

import matplotlib.pyplot as plt
import random
e=2.718281828
func= lambda x,y :-x**2 + y**2
class optim():
    def __init__(self,entry_x,entry_y,vf_x,vf_y,func):
        self.grad=np.zeros((2))
        self.xy=np.array([float(entry_x),float(entry_y)])
        self.phase_vec=np.array([vf_x,vf_y])
        self.func=func
        self.coef_grad=0.4
        self.step=0
        self.prev_maxgrad=1e100

        self.vk=np.zeros((2))
        self.Gk=0
        

    def mirror(self):
        mul_res=np.sum(self.phase_vec*self.grad)
        sqr_res=np.sum(self.phase_vec*self.phase_vec)
        mirror_grad_cos=min(0.4,abs(mul_res/(sqr_res*np.sum(self.grad*self.grad))**0.5))
        #print( f"pre-mirr_g cos {abs(mul_res/(sqr_res*np.sum(self.grad*self.grad))**0.5)}")
        #print (f"mirrorgrad cos {mirror_grad_cos}")
        self.grad=np.subtract(self.grad,(1+mirror_grad_cos)*np.multiply(mul_res/sqr_res,self.phase_vec))
    
    def get_grad(self):
        dx=lambda x,y,f:(-f(x-0.001,y)+f(x+0.001,y))/0.002
        dy=lambda x,y,f:(-f(x,y-0.001)+f(x,y+0.001))/0.002
        grad=np.array([dx(self.xy[0],self.xy[1],self.func),dy(self.xy[0],self.xy[1],self.func)])
        return grad, max(abs(grad))
    
    def apply_grad(self):    
        self.xy-=self.grad*self.coef_grad
    
    def proceed(self):
        self.xs=[self.xy[0]]
        self.ys=[self.xy[1]]
        while(self.move_DoFs()>0.001):
            #print(f"step {self.step} ({self.xy[0]}, {self.xy[1]})")
            self.xs.append(self.xy[0])
            self.ys.append(self.xy[1])
            #input()
            if self.step>=100000:
                print("!!!!!")
                break
        print(f"step {self.step} ({self.xy[0]}, {self.xy[1]})")
        self.get_grad()
        print(f"get      {self.grad}")

        self.mirror()
        print(f"mirrored {self.grad}\nv = {self.phase_vec}")

        return self.xs,self.ys
  
    def move_DoFs(self):
        b1=0.33
        b2=0.99
        eps=1e-8
        TRUST_RAD=1000
        
        self.grad, maxgrad=self.get_grad()
        self.mirror()
        if(maxgrad*self.coef_grad>TRUST_RAD):
            self.coef_grad=TRUST_RAD/maxgrad

        if self.step>0 or 1:
            #ADAM
            self.vk = b1*self.vk + (1-b1)*self.grad#*self.coef_grad
            self.Gk = b2*self.Gk + (1-b2)*np.sum(self.grad*self.grad)#*self.coef_grad**2
            self.xy=self.xy - 4.e-1*(self.Gk+eps)**(-0.5) * self.vk
            
        else:#GD - потому что первые 20 происходит значительная смена параметров, и нечего давать её в инерционный алгоритм
            self.apply_grad()
        
        self.step+=1
        if maxgrad<self.prev_maxgrad:
            self.coef_grad*=1.01
        elif maxgrad>self.prev_maxgrad:
            self.coef_grad*=0.9
            if self.coef_grad<0.4:
                self.coef_grad=0.4
        
        #print(f"coef grad {self.coef_grad}")
        self.prev_maxgrad=maxgrad
        return maxgrad

def frange(x, y, jump):
  n=0
  while x+n*jump < y:
    yield x+n*jump
    n+=1

plt.axes().set_aspect(1)  
plt.grid(zorder=0)

xc, yc = np.meshgrid(np.linspace(-6, 6, 100),  
                   np.linspace(-6, 6, 100)) 
plt.contour(xc, yc, func(xc,yc),linewidths=2.4,zorder=2) 


J_MIN=-5
J_MAX=5
I_MIN=-5
I_MAX=5

for i in frange (I_MIN,I_MAX+1,3.333):
    for j in frange(J_MIN,J_MAX+1,3.333):
        trace=optim(i+(random.random()-0.5)/10,j+(random.random()-0.5)/10,1,0,func)
        xs,ys=trace.proceed()

        #color_tuple=(((i-I_MIN)/(I_MAX-I_MIN))**2/2+0.5,(1-((j-J_MIN)/(J_MAX-J_MIN)+(j-J_MIN)/(J_MAX-J_MIN))/2)**2,((j-J_MIN)/(J_MAX-J_MIN))**2)
        color_tuple=(1,0,0)
        plt.plot(xs, ys,color=color_tuple,linewidth=1,zorder=2*(3+(i-I_MIN)*(J_MAX-J_MIN)+(j-J_MIN)))
        plt.scatter([xs[0]], [ys[0]],color="black",linewidth=1,zorder=2*(3+(i-I_MIN)*(J_MAX-J_MIN)+(j-J_MIN))+1)
        plt.scatter([xs[-1]], [ys[-1]],color="b",marker="+",linewidth=1,zorder=2*(3+(i-I_MIN)*(J_MAX-J_MIN)+(j-J_MIN))+1)
        
        trace=0
 
'''
I_MIN=-3.1415+3.1415/32
I_MAX=3.1415

for i in frange (I_MIN,I_MAX+0.1,3.1415*2/16):
    trace=optim(np.sin(i)/1.2,np.cos(i)/1.2,np.sin(i),np.cos(i),func)
    xs,ys=trace.proceed()

    #color_tuple=(((i-I_MIN)/(I_MAX-I_MIN))**2/2+0.5,(1-((j-J_MIN)/(J_MAX-J_MIN)+(j-J_MIN)/(J_MAX-J_MIN))/2)**2,((j-J_MIN)/(J_MAX-J_MIN))**2)
    color_tuple=(random.random(),random.random()/10,random.random()/10)
    plt.plot(xs, ys,color=color_tuple,linewidth=1,zorder=2*(3+(i-I_MIN)))
    plt.scatter([xs[0]], [ys[0]],color="black",linewidth=1,zorder=2*(3+(i-I_MIN))+1)
    trace=0
'''
plt.xlim(-6.5, 6.5) 
plt.ylim(-6.5, 6.5) 

plt.savefig("pictures/fig_trace_multiple_points_saddle_roubst",dpi=300)