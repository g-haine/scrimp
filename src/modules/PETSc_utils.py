# -*- coding: utf-8 -*-
"""
Utility functions for petsc4py.

@author: Florian Monteghetti
"""
import numpy as np
import scipy.sparse as sp
from petsc4py import PETSc

def convert_scipy_to_PETSc(A,comm=None):
    """
    Convert from scipy sparse to PETSc sparse.

    Parameters
    ----------
    A : scipy.sparse matrix
        
    comm : MPI communicator PETSc.Comm, optional
        
    Returns
    -------
    A_petsc : PETSc.Mat
    
    """
    if A.format!="csr":
        A_csr = A.tocsr(copy=False)
    else:
        A_csr = A
        
    A_petsc = PETSc.Mat().createAIJ(size=A_csr.shape,
                                   csr=(A_csr.indptr, A_csr.indices,
                                       A_csr.data),comm=comm)
    return A_petsc 
    

def convert_PETSc_to_scipy(A):
    """
    Convert from PETSc.Mat to scipy sparse (csr).

    Parameters
    ----------
    A : PETSc.Mat


    Returns
    -------
    A_scipy : scipy.sparse.csr.csr_matrix
        

    """
    (indptr,indices,val) = A.getValuesCSR()
    A_scipy = sp.csr_matrix((val, indices, indptr), shape=A.size)
    return A_scipy

def set_diagonal_entries(A):
    """
    Set all diagonal entries of A without modifying existing entries.

    Parameters
    ----------
    A : PETSc.Mat

    """

    x=A.createVecRight(); x.setArray(0)
    A.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    A.setDiagonal(x,addv=True)
    A.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR,True)

def get_cleared_options_db():
    """ Clear the option database and return a pointer to it. """
    OptDB = PETSc.Options() # pointer to THE option database
    for key in OptDB.getAll(): OptDB.delValue(key) # Clear database
    return OptDB

######################### PETSc TS ###########################################
import time

class Fully_implicit_DAE(object):
    """ DAE under a fully implicit form F(t,y,yd) = 0 suitable for
    integration with PETSc TS fully implicit schemes. """ 
    
    def __init__(self,N,name,idx_alg=None):
        self.name = name        # str: name
        self.N = N              # int: number of unknowns
        self.idx_alg = idx_alg  # list(int): algebraic variables
        self.F = 0              # PETSc.Vec(): vector used during time-integration
        self.J = 0              # PETSc.Mat(): matrix used during time-integration
           
    def init(self):
        """ Initialize. Call before new time integration. """
            # History
        self.init_history()

    def IFunction(self, ts, t, y, yd, F):
        """ Evaluate residual vector F(t,y,yd)."""
        pass
        
    def IJacobian(self, ts, t, y, yd, c, Amat, Pmat):
        """ Evaluate jacobian matrix J = c*J_yd + J_y."""
        pass
    
    def init_history(self):
        """ Initialize history structure. """
        self.history = dict(t=list(),z=list())#,ally=list()) # z because y is for the output in the pHs framework
        self.t_start = 0
        
    def monitor(self, ts, i, t, y, dt_export=1, t0=0.0):
        """ Monitor to use during iterations. """
        if self.history['t']:
            lastt = self.history['t'][-1]            
        else:
            lastt = t-dt_export
            self.t_start = time.time()
        if (t >= lastt + 0.95*dt_export) or (i==-1) or (t==t0):
            print(f"i={i:8d} t={t:8g} * ({int(time.time()-self.t_start)}s)")
            self.history['t'].append(t)
            self.history['z'].append(y.copy()) # z because y is for the output in the pHs framework
        else: # (i<10): 
            print(f"i={i:8d} t={t:8g}   ({int(time.time()-self.t_start)}s)")
        
        # if self.history['ally']:
        #     lasty = self.history['ally'][-1]
        #     diffy = lasty.copy()
        #     diffy.aypx(-1.,y)
        #     diff = diffy.norm(norm_type=2)
        #     self.history['ally'].append(y.copy()) 
        #     print(f"Absolute difference: {diff}")
        # else:
        #     norm = y.norm()
        #     self.history['ally'].append(y.copy()) 
        #     print(f"Initial norm: {norm}")
           
    def get_vector(self):
        """ Get new residual-size vector. """
        x = PETSc.Vec().createWithArray(np.zeros(self.N,))
        return x

def TS_exclude_var_from_lte(dae,ts):
    """ Exclude algebraic variables from local truncation error. 
    
    Inputs
    -------
    dae: Fully_implicit_DAE
    
    ts: PETSc.TS
    
    """
    atol = ts.getTolerances()[1]
    atol_v = atol*np.ones((dae.N,))
    atol_v[dae.idx_alg] = np.inf # inf for algebraic variables
    atol_v_petsc = dae.get_vector()
    atol_v_petsc.setArray(atol_v)
    ts.setTolerances(atol=atol_v_petsc)
    
def TS_integration_dae(dae,x0,dt,tf,dt_export=None,
                        cinit=False,cinit_dt=1e-2,cinit_nstep=1,t0=0.0,interpolate_last=True):
    """ Time-integration of DAE under fully implicit form F(t,y,yd)=0.
    
    Inputs
    -------
    dae: Fully_implicit_DAE
    
    x0: PETSc.Vec
        Initial condition. This vector is used throughout the time integration.
        
    dt,tf: real
        Time step and final time.
    
    dt_export: real (Optional)
        Time step for saving solution.
        
    cinit, cinit_dt, cinit_nstep: bool, real, int
        Parameters related to consistent initialization.
    
    """
    tc = -cinit*1.e10
    # if dt_export is None:
    #     dt_export = dt
    dae.init()
    # monitor = lambda ts, i, t, x: dae.monitor(ts, i, t, x, dt_export=dt_export, t0=t0)
    OptDB = PETSc.Options()
    import time
        # -- Consistent initialization
    if cinit==True:
        # print(f"Consistent initialization: {cinit_nstep} step(s) of Backward Euler.")
        print(f"Consistent initialization: {cinit_nstep} step(s) max. of Pseudo TS.")
        if OptDB.hasName("ts_type"):
            tsType = OptDB.getString("ts_type")
        else:
            tsType = "ts_ssp"
        monitor_init = lambda ts, i, t, x: dae.monitor(ts, i, t, x, dt_export=cinit_dt*cinit_nstep, t0=t0)
        # OptDB["ts_type"] = "beuler"
        OptDB["ts_type"] = "pseudo"
        OptDB["ts_pseudo_increment"] = 1.1
        OptDB["ts_pseudo_fatol"] = 1e-6
        OptDB["ts_pseudo_frtol"] = 1e-9
        # if OptDB.hasName("ts_atol"):
        #     tsAtol = OptDB.getString("ts_atol")
        # else:
        #     tsAtol = "No_Atol"
        #     OptDB["ts_atol"] = 1e-6
        # if OptDB.hasName("ts_rtol"):
        #     tsRtol = OptDB.getString("ts_rtol")
        # else:
        #     tsRtol = "No_Rtol"
        #     OptDB["ts_Rtol"] = 1e-9
        ts = PETSc.TS().create()
        ts.setMonitor(monitor_init)
        ts.setIFunction(dae.IFunction, dae.F)
        ts.setIJacobian(dae.IJacobian, dae.J)
        ts.setMaxSNESFailures(-1)
        ts.setFromOptions()
        ts.setTime(tc)
        # ts.setMaxTime(cinit_nstep*cinit_dt)
        # ts.setTimeStep(cinit_dt)
        # ts.setMaxSteps(cinit_nstep)
        ts.setMaxTime(0.)
        ts.setTimeStep(cinit_dt)
        ts.setMaxSteps(cinit_nstep)
        ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)
        start = time.time()
        ts.solve(x0)
        tc = ts.getTime()
        print(f"Done (t={tc:8g}).\nElapsed time: {time.time()-start:1.4g}s")
        ts.destroy()
        del ts
        OptDB.delValue("ts_pseudo_increment")
        # OptDB.delValue("ts_atol")
        # OptDB.delValue("ts_rtol")
        if (tsType=="ts_ssp"):
            OptDB.delValue("ts_type")
        else:
            OptDB["ts_type"] = tsType
        # if (tsAtol=="No_Atol"):
        #     OptDB.delValue("ts_atol")
        # else:
        #     OptDB["ts_atol"] = tsAtol
        # if (tsRtol=="No_Rtol"):
        #     OptDB.delValue("ts_rtol")
        # else:
        #     OptDB["ts_rtol"] = tsRtol
        # -- Main integration
    if dt_export is None:
        dt_export = dt
    dae.init()
    monitor = lambda ts, i, t, x: dae.monitor(ts, i, t, x, dt_export=dt_export, t0=t0)
    ts = PETSc.TS().create()
    ts.setMonitor(monitor)
    ts.setIFunction(dae.IFunction, dae.F)
    ts.setIJacobian(dae.IJacobian, dae.J)
    ts.setTime(t0)
    ts.setMaxTime(tf)
    ts.setTimeStep(dt)
    ts.setMaxSNESFailures(-1)
    # ts.setExactFinalTime(PETSc.TS.ExactFinalTime.INTERPOLATE)
    ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)
    ts.setFromOptions()
    TS_exclude_var_from_lte(dae,ts)
    import time
    start = time.time()
    ts.solve(x0)
    print(f"Elapsed time: {time.time()-start:1.4g}s")
    print(f"Steps: {ts.getStepNumber()} ({ts.getStepRejections()} rejected, {ts.getSNESFailures()} Nonlinear solver failures)")
    print(f"Nonlinear iterations: {ts.getSNESIterations()}, Linear iterations: {ts.getKSPIterations()}")
    ts.reset()
    ts.destroy()
