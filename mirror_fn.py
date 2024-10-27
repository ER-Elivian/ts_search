def mirror_fn(grad,#list 3*N, gradient
              xyzs,#list 4*N, coord i.e."C, 0, 1, 1"
              search_DoFs#list of reaction dofs: type, athoms, value. i.e. "b 1 2 -1" for stratching (-1) bond (b) between 1 and 2 athoms
              ):

    import copy, numpy as np
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

    for i in range(nAthoms):
        v_len=np.linalg.norm(mirror_vec[i])
        if v_len > 0.00001:
            np.multiply(1/v_len, mirror_vec[i])

    mul_res=np.sum(mirror_vec*m_grad)
    sqr_res=np.sum(mirror_vec*mirror_vec)
    
    mirror_grad_cos=min(0.5,abs(mul_res/(sqr_res*np.sum(m_grad*m_grad))**0.5))
    print (f"mgcos {mirror_grad_cos}")
    m_grad=np.subtract(m_grad,(1+mirror_grad_cos**2)*np.multiply(mul_res/sqr_res,mirror_vec))
    return m_grad