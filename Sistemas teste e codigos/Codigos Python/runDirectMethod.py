from numpy import zeros, ones, conj, exp, diag, asmatrix, asarray, concatenate, sum, linalg
from scipy.sparse import issparse, csr_matrix as sparse
import array as arr
import numpy as np
import math
from dSbus_dV import dSbus_dV
from d2Sbus_dV2 import d2Sbus_dV2


def runDirectMethod(wP, wQ, tol, pq, npq, pv, npv, sYbus, magV, angV, pg0, qg0, pl0, ql0, nbus, itMax, kp):

    t = zeros((1, 1))  # tamanho do crescimento

    pg = pg0
    qg = qg0

    espS = zeros(nbus, dtype=complex)
    espP = zeros(nbus, dtype=float)
    espQ = zeros(nbus, dtype=float)

    calcV = zeros(nbus, dtype=complex)
    calcI = zeros(nbus, dtype=complex)

    it = 0
    while True:

        it += 1

        pl = (1+t*kp)*pl0
        ql = ql0

        espP = pg-pl
        espQ = qg-ql
        espS.real = espP
        espS.imag = espQ

        calcV = magV*exp(angV*1j)
        calcI = sYbus@calcV
        calcS = calcV*conj(calcI)

        eS = calcS - espS
        eP = eS.real
        eQ = eS.imag

        f = concatenate([eP[pv+pq], eQ[pq]], axis=0)
        f = f.reshape(f.size, 1)

        [dS_magV, dS_angV] = dSbus_dV(sYbus, calcV)  # Complex Derivatives
        A11 = (dS_angV.A[pv+pq, :][:, pv+pq]).real  # dP_dAngV
        A12 = (dS_magV.A[pv+pq, :][:, pq]).real  # dP_dMagV
        A21 = (dS_angV.A[pq, :][:, pv+pq]).imag  # dQ_AngV
        A22 = (dS_magV.A[pq, :][:, pq]).imag  # dQ_MagV

        df_dx = concatenate((concatenate((A11, A21), axis=0),
                            concatenate((A12, A22), axis=0)), axis=1)

        dpl_dt = pl0
        dql_dt = zeros((nbus,))
        df_dt = concatenate((dpl_dt[pv+pq], dql_dt[pq]), axis=0)
        df_dt = df_dt.reshape(df_dt.size, 1)

        df_dw = zeros((df_dx.shape[1], npv+2*npq))

        g = df_dx.T@concatenate((wP[pv+pq], wQ[pq]), axis=0)  # g function

        h = sum(concatenate((wP[pv+pq], wQ[pq]), axis=0) *
                concatenate((wP[pv+pq], wQ[pq]), axis=0))-1  # h function
        h = h.reshape(1, 1)

        F = concatenate((f, g, h), axis=0)

        ## ---- Function g(x,t,u) = Dxf.'w derivatives ---- ##
        [GAA1, GAV1, GVA1, GVV1] = d2Sbus_dV2(
            sYbus, calcV, wP)  # complex jacobian matrix Dxxf.*uP
        [GAA2, GAV2, GVA2, GVV2] = d2Sbus_dV2(
            sYbus, calcV, wQ)  # complex jacobian matrix Dxxf.*uQ

        M1 = concatenate((concatenate((GAA1.A[pv+pq, :][:, pv+pq], GVA1.A[pq, :][:, pv+pq]), axis=0),
                         concatenate((GAV1.A[pv+pq, :][:, pq], GVV1.A[pq, :][:, pq]), axis=0)), axis=1)
        M2 = concatenate((concatenate((GAA2.A[pv+pq, :][:, pv+pq], GVA2.A[pq, :][:, pv+pq]), axis=0),
                         concatenate((GAV2.A[pv+pq, :][:, pq], GVV2.A[pq, :][:, pq]), axis=0)), axis=1)

        dg_dx = M1.real + M2.imag
        dg_dt = zeros(df_dt.shape)
        dg_dw = df_dx.T

        ## ---- Function h(w) = w*w.'- 1 derivatives ---- ##
        dh_dx = zeros((1, df_dx.shape[1]))
        dh_dt = zeros((1, 1))
        dh_dw = 2*concatenate((wP[pv+pq], wQ[pq]), axis=0).T

        ## ---- NR Loop ---- ##
        J = concatenate((concatenate((df_dx, dg_dx, dh_dx), axis=0), concatenate(
            # complete jacobian matrix
            (df_dt, dg_dt, dh_dt), axis=0), concatenate((df_dw, dg_dw, dh_dw), axis=0)), axis=1)
        angV_PQV = angV[pv+pq]
        angV_PQV = angV_PQV.reshape(npv+npq, 1)  # %previous state vector
        magV_PQ = magV[pq]
        magV_PQ = magV_PQ.reshape(npq, 1)
        x = concatenate((angV_PQV, magV_PQ, t, wP[pv+pq], wQ[pq]), axis=0)
        s = linalg.lstsq(J, F, rcond=None)  # convergence error
        s = s[0].reshape(x.shape)
        s = s.reshape(x.shape)
        x = x-s  # %next state vector
        x = x.reshape(x.size,)

        ## ---- State Variable Update ---- ##
        angV[pv+pq] = x[0:(npv+npq),]
        magV[pq] = x[(npv+npq):npv+2*npq,]
        t = x[npv+2*npq]
        t = t.reshape((1, 1))

        x = x.reshape(x.size, 1)
        wP[pv+pq] = x[npv+2*npq+1:2*npv+3*npq+1,]
        wQ[pq] = x[2*npv+3*npq+1:x.size,]

        if it >= itMax or linalg.norm(s) <= tol:
            dPtotal = sum(t*kp@pl0)
            break

    return (t, dPtotal, it)
