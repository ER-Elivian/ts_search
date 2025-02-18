def mirror_fn(grad,#list 3*N, gradient
              xyzs,#list 4*N, coord i.e."C, 0, 1, 1"
              search_DoFs#list of reaction dofs: type, athoms, value. i.e. "b 1 2 -1" for stratching (-1) bond (b) between 1 and 2 athoms
              ):

    import numpy as np
    from angle_3d_rev import find_vectors as f_v_a
    from dihedral_3d_rev import find_vectors as f_v_d
    m_grad=np.array(grad)
    nAthoms=len(xyzs)
    
    v0=[0.,0.,0.]
    mirror_vec=[]
    for i in range(nAthoms):
        mirror_vec.append(v0)
    mirror_vec=np.array(mirror_vec)
    
    for dof in search_DoFs:
        if dof[0]=='b':
            direction=np.subtract(xyzs[dof[1]-1], xyzs[dof[2]-1])
            direction=np.multiply(dof[3]/np.linalg.norm(direction),direction)
            mirror_vec[dof[1]-1]-=direction
            mirror_vec[dof[2]-1]+=direction
        
        if dof[0]=='a':
            vectors=f_v_a([xyzs[dof[1]-1],xyzs[dof[2]-1],xyzs[dof[3]-1]])
            mirror_vec[dof[1]-1]+=dof[4]*vectors[0]
            mirror_vec[dof[2]-1]+=dof[4]*vectors[1]
            mirror_vec[dof[3]-1]+=dof[4]*vectors[2]

        if dof[0]=='d':
            vectors=f_v_d([xyzs[dof[1]-1],xyzs[dof[2]-1],xyzs[dof[3]-1],xyzs[dof[4]-1]])
            mirror_vec[dof[1]-1]+=dof[4]*vectors[0]
            mirror_vec[dof[2]-1]+=dof[4]*vectors[1]
            mirror_vec[dof[3]-1]+=dof[4]*vectors[2]
            mirror_vec[dof[4]-1]+=dof[4]*vectors[3]



    for i in range(nAthoms):
        v_len=np.linalg.norm(mirror_vec[i])
        if v_len > 0.01:
            np.multiply(1/v_len, mirror_vec[i])


    mul_res=np.sum(mirror_vec*m_grad)
    sqr_res=np.sum(mirror_vec*mirror_vec)
    
    mirror_grad_cos=min(0.6,abs(mul_res/(sqr_res*np.sum(m_grad*m_grad))**0.5)**4)
    print(abs(mul_res/(sqr_res*np.sum(m_grad*m_grad))**0.5))
    print (f"mgcos {mirror_grad_cos}")
    m_grad=np.subtract(m_grad,(1+mirror_grad_cos)*np.multiply(mul_res/sqr_res,mirror_vec))
    m_grad_mean=np.array([0,0,0])
    for i in range(nAthoms):
        m_grad_mean=np.add(m_grad_mean,m_grad[i])
    m_grad_mean/=nAthoms
    
    for i in range(nAthoms):
        m_grad[i]-=m_grad_mean
    
    return m_grad