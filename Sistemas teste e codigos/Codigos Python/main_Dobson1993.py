from numpy import zeros, ones, conj, exp, diag, asmatrix, asarray, concatenate, sum, linalg, round, array
from scipy.sparse import issparse, csr_matrix as sparse
from calcYbus import calcYbus
from copy import deepcopy

import array as arr
import numpy as np
import math


from importSystemData import importSystemData
from dSbus_dV import dSbus_dV
from d2Sbus_dV2 import d2Sbus_dV2


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#Leitura do arquivo 
system = 'ieee14cdf.txt'

[dbus,dbranch,nbus,MVAbase] = importSystemData(system)

#Calculo da matriz Ybus
ybus = calcYbus(nbus, dbus, dbranch)    
sYbus = sparse(ybus)

typeArr = arr.array('i',dbus['type'].values)
pl0 = dbus['lMW'].values / MVAbase
ql0  = dbus['lMVAR'].values / MVAbase
pg0 = dbus['gMW'].values / MVAbase
qg0  = dbus['gMVAR'].values / MVAbase
magV = dbus['magV'].values 
angV = dbus['angV'].values*math.pi/180

pv = [i for i in range(len(typeArr)) if typeArr[i] == 2]
pq = [i for i in range(len(typeArr)) if typeArr[i] == 0]
ref = [i for i in range(len(typeArr)) if typeArr[i] == 3]

npv = len(pv)
npq = len(pq)
nbus = npv + npq + 1

## ---- Initial Values ---- ##
#eigenvalues
wP = zeros((nbus,1)) 
wP[pv+pq] = 1
wQ = zeros((nbus,1))
wQ[pq] = 1
#load growth parameter
t = ones((1,1))
#maximum number of iterations
itMax = 20

pg = pg0; qg = qg0

espS = zeros(nbus, dtype = complex)    
espP = zeros(nbus, dtype = float)            
espQ = zeros(nbus, dtype = float)            

calcV = zeros(nbus, dtype = complex)   
calcI = zeros(nbus, dtype = complex)   

tol = 1e-7
idx = 0
lambP = zeros((nbus,1),)
lambQ = zeros((nbus,1),)
while True:
    if idx == 0:    
            n0 = zeros((nbus*2,1))
            n0[:nbus,0] = pl0/linalg.norm(pl0)   
            n0[nbus:,0] = ql0/linalg.norm(ql0)  
            n0v = deepcopy(n0)
    it = 0 
    while True:    
        it += 1    
        pl = pl0+t*n0[:nbus,0];   ql = ql0+t*n0[nbus:,0]
        
        espP = pg-pl; espQ = qg-ql
        espS.real = espP;  espS.imag = espQ

        calcV = magV*exp(angV*1j); calcI = sYbus@calcV; 
        calcS = calcV*conj(calcI)

        eS = calcS - espS
        eP = eS.real
        eQ = eS.imag

        
        f = concatenate([eP[pv+pq],eQ[pq]],axis=0)
        f = f.reshape(f.size,1)

        [dS_magV,dS_angV] = dSbus_dV(sYbus,calcV) #Complex Derivatives
        A11 = (dS_angV.A[pv+pq,:][:,pv+pq]).real  #dP_dAngV
        A12 = (dS_magV.A[pv+pq,:][:,pq]).real #dP_dMagV
        A21 = (dS_angV.A[pq,:][:,pv+pq]).imag #dQ_AngV
        A22 = (dS_magV.A[pq,:][:,pq]).imag #dQ_MagV

        df_dx = concatenate((concatenate((A11,A21),axis=0),concatenate((A12,A22),axis=0)),axis=1)
        
        dpl_dt = n0[:nbus,0]; dql_dt = n0[nbus:,0]
        df_dt = concatenate((dpl_dt[pv+pq],dql_dt[pq]),axis = 0)
        df_dt = df_dt.reshape(df_dt.size,1)

        df_dw = zeros((df_dx.shape[1],npv+2*npq))
    
        g = df_dx.T@concatenate((wP[pv+pq],wQ[pq]),axis=0) #g function 
        
        h = sum(concatenate((wP[pv+pq],wQ[pq]),axis=0)*concatenate((wP[pv+pq],wQ[pq]),axis=0))-1 #h function 
        h = h.reshape(1,1)
        
        F = concatenate((f,g,h),axis=0)
        
        ## ---- Function g(x,t,u) = Dxf.'w derivatives ---- ##  
        [GAA1, GAV1, GVA1, GVV1] = d2Sbus_dV2(sYbus, calcV, wP) #complex jacobian matrix Dxxf.*uP
        [GAA2, GAV2, GVA2, GVV2] = d2Sbus_dV2(sYbus, calcV, wQ) #complex jacobian matrix Dxxf.*uQ

        M1 = concatenate((concatenate((GAA1.A[pv+pq,:][:,pv+pq],GVA1.A[pq,:][:,pv+pq]),axis=0),concatenate((GAV1.A[pv+pq,:][:,pq],GVV1.A[pq,:][:,pq]),axis=0)),axis=1)
        M2 = concatenate((concatenate((GAA2.A[pv+pq,:][:,pv+pq],GVA2.A[pq,:][:,pv+pq]),axis=0),concatenate((GAV2.A[pv+pq,:][:,pq],GVV2.A[pq,:][:,pq]),axis=0)),axis=1)

        dg_dx = M1.real + M2.imag 
        dg_dt = zeros(df_dt.shape)
        dg_dw = df_dx.T

        ## ---- Function h(w) = w*w.'- 1 derivatives ---- ##
        dh_dx = zeros((1,df_dx.shape[1]))
        dh_dt = zeros((1,1))
        dh_dw = 2*concatenate((wP[pv+pq],wQ[pq]),axis=0).T

        ## ---- NR Loop ---- ##
        J = concatenate((concatenate((df_dx,dg_dx,dh_dx),axis=0), concatenate((df_dt,dg_dt,dh_dt),axis=0), concatenate((df_dw,dg_dw,dh_dw),axis=0)), axis = 1) #complete jacobian matrix
        angV_PQV = angV[pv+pq]; angV_PQV = angV_PQV.reshape(npv+npq,1) #%previous state vector
        magV_PQ = magV[pq]; magV_PQ = magV_PQ.reshape(npq,1)
        x = concatenate((angV_PQV,magV_PQ,t,wP[pv+pq],wQ[pq]),axis = 0)
        s = linalg.lstsq(J,F) #convergence error
        s = s[0].reshape(x.shape)
        s = s.reshape(x.shape)
        x = x-s # %next state vector
        x = x.reshape(x.size,)

        ## ---- State Variable Update ---- ##    
        angV[pv+pq] = x[0:(npv+npq),]
        magV[pq] = x[(npv+npq):npv+2*npq,]
        t = x[npv+2*npq]
        t = t.reshape((1,1))

        x = x.reshape(x.size,1)
        wP[pv+pq] = x[npv+2*npq+1:2*npv+3*npq+1,]
        wQ[pq] = x[2*npv+3*npq+1:x.size,]

        if it >= itMax or linalg.norm(s) <= tol:
            break   
    if idx == 0:
        print('-----------RESULTADO--------')       
    idx+=1    
    print('Iteracao: ', idx,'\n')
    n0[pv+pq] =wP[pv+pq]
    n0[array(pq)+nbus] =wQ[pq]
    n0v = concatenate((n0v, n0), axis=1)    
    if idx == 1: 
        lambP = round(t*n0[:nbus], 3)
        lambQ = round(t*n0[nbus:], 3)
    else: 
        lambP = concatenate((lambP, round(t*n0[:nbus], 3)), axis=1)
        lambQ = concatenate((lambQ, round(t*n0[nbus:], 3)), axis=1)
    idx_minP = 0
    value_minP = 0    
    for i in range (npv):
        if idx_minP == 0:
            idx_minP = pv[i]
            value_minP = lambP[pv[i],idx-1]
        else:
            if value_minP > lambP[pv[i],idx-1]:
                idx_minP = pv[i]
                value_minP = lambP[pv[i],idx-1] 
    for i in range (npq):
        if idx_minP == 0:
            idx_minP = pq[i]
            value_minP = lambP[pq[i],idx-1]
        else:
            if value_minP > lambP[pq[i],idx-1]:
                idx_minP = pq[i]
                value_minP = lambP[pq[i],idx-1] 
    print('Minimo P-P0 na barra ', dbus['name'][idx_minP], 'com valor igual a ', round(value_minP*MVAbase,3), 'MW', '\n')             
    idx_minQ = 0
    value_minQ = 0
    for i in range (npq):
        if idx_minQ == 0:
            idx_minQ = pq[i]
            value_minQ = lambQ[pq[i],idx-1]
        else:
            if value_minQ > lambQ[pq[i],idx-1]:
                idx_minQ = pq[i]
                value_minQ = lambQ[pq[i],idx-1] 
    print('Minimo Q-Q0 na barra ', dbus['name'][idx_minQ], 'com valor igual a ', round(value_minQ*MVAbase,3), 'Mvar', '\n')
   
    if idx >= 2 and linalg.norm(n0v[:,idx-1] - n0v[:,idx-2]) <= 1e-4:
        break    
   
teste = 1


