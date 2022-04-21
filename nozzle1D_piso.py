#!/usr/bin/env python3.3

import numpy as np
import scipy.sparse as sps
from scipy.sparse import linalg as spslinalg
from scipy import linalg

np.set_printoptions(linewidth=np.nan)

rho = 1
mass_flow = 1
L = 2
AA = 0.5
AE = 0.1
pA = 9.6
pE = 0
pNx = 5
uNx = pNx - 1

Ap = np.zeros(pNx)
for c in range(0, pNx):
    Ap[c] = AA - (AA - AE)/(pNx - 1) * c

Au = np.zeros(uNx)
for c in range(0, uNx):
    Au[c] = (Ap[c] + Ap[c +1])/2

u = mass_flow/(rho * Au)
p = np.zeros(pNx)
for c in range(0, pNx):
    p[c] = pA - (pA - pE)/(pNx - 1) * c

Flux = np.zeros(pNx)

aW = np.zeros(uNx)
aE = np.zeros(uNx)
aP = np.zeros(uNx)
Su = np.zeros(uNx)
uu = np.zeros(uNx)
uuu = np.zeros(uNx)
bla = np.zeros(uNx)

aWp = np.zeros(pNx)
aEp = np.zeros(pNx)
aPp = np.zeros(pNx)
b = np.zeros(pNx - 2)

aWpp = np.zeros(pNx)
aEpp = np.zeros(pNx)
aPpp = np.zeros(pNx)
bp = np.zeros(pNx - 2)

alpha = 0.5
alphap = 0.5
# alpha = 1.0
# alphap = 1.0


A = sps.coo_matrix((uNx, uNx), dtype=np.float64)
A.row = np.array(np.arange(uNx), dtype=np.intc)
A.col = np.array(np.arange(uNx), dtype=np.intc)
A.row = np.array(np.append(A.row, np.arange(1,uNx)), dtype=np.intc)
A.col = np.array(np.append(A.col, np.arange(0,uNx-1)), dtype=np.intc)

C = sps.coo_matrix((pNx - 2, pNx - 2), dtype=np.float64)
C.row = np.array(np.arange(pNx -2), dtype=np.intc)
C.col = np.array(np.arange(pNx - 2), dtype=np.intc)
C.row = np.array(np.append(C.row, np.arange(1,pNx -2)), dtype=np.intc)
C.col = np.array(np.append(C.col, np.arange(0,pNx-3)), dtype=np.intc)
C.row = np.array(np.append(C.row, np.arange(0,pNx -3)), dtype=np.intc)
C.col = np.array(np.append(C.col, np.arange(1,pNx-2)), dtype=np.intc)

Cp = sps.coo_matrix((pNx - 2, pNx - 2), dtype=np.float64)
Cp.row = np.array(np.arange(pNx -2), dtype=np.intc)
Cp.col = np.array(np.arange(pNx - 2), dtype=np.intc)
Cp.row = np.array(np.append(Cp.row, np.arange(1,pNx -2)), dtype=np.intc)
Cp.col = np.array(np.append(Cp.col, np.arange(0,pNx-3)), dtype=np.intc)
Cp.row = np.array(np.append(Cp.row, np.arange(0,pNx -3)), dtype=np.intc)
Cp.col = np.array(np.append(Cp.col, np.arange(1,pNx-2)), dtype=np.intc)


for it in range(0, 200):

    for c in range(1, pNx - 1):
        Flux[c] = rho * (u[c] + u[c - 1])/2 * Ap[c]
    Flux[0] = rho * u[0] * Au[0]
    Flux[-1] = rho * u[-1] * Au[-1]

    for c in range(1, uNx):
        aW[c] = Flux[c]

    for c in range(0, uNx):
        aP[c] = Flux[c +1]

    aP = aP/alpha
    d = Au/aP

    for c in range(0, uNx):
        Su[c] = (p[c] - p[c +1]) * Au[c] + (1 - alpha) * aP[c] * u[c]

    A.data = aP
    A.data = np.append(A.data, -aW[1:uNx])
    u = spslinalg.bicgstab(A, Su, tol=1e-6)[0]


    for c in range(1, pNx - 1):
        aWp[c] = rho * d[c - 1] * Au[c - 1]
        aEp[c] = rho * d[c] * Au[c]
        aPp[c] = aWp[c] + aEp[c]

    for c in range(0, pNx - 2):
        b[c] = rho * u[c] * Au[c] - rho * u[c + 1] * Au[c + 1]



    C.data = aPp[1:-1]
    C.data = np.append(C.data, -aWp[2:pNx -1])
    C.data = np.append(C.data, -aEp[1:pNx -2])

    pc = spslinalg.bicgstab(C, b, tol=1e-6)[0]


    for c in range(1, pNx - 1):
        p[c] = p[c] + alphap * pc[c - 1]

    for c in range(1, uNx - 1):
        uu[c] = u[c] + d[c] * (pc[c - 1] - pc[c])
    uu[0] = u[0] + d[0] * (0 - pc[0])
    uu[-1] = u[-1] + d[-1] * (pc[-1] - 0)

    bla = u.copy()

    du = uu -u

    aa = aW.copy()
    aa[0] = 0

    for c in range(1, uNx):
        aa[c] = aW[c] * du[c - 1]
    # aa[2] = aW[2] * du[1]
    # aa[3] = aW[3] * du[2]
    # aa[4] = aW[4] * du[3]
    # aa[5] = aW[5] * du[4]
    # aa[6] = aW[6] * du[5]
    # aa[7] = aW[7] * du[6]

#eigentlich sinnlos
    for c in range(1, pNx - 1):
        aWpp[c] = rho * d[c - 1] * Au[c - 1]
        aEpp[c] = rho * d[c] * Au[c]
        aPpp[c] = aWpp[c] + aEpp[c]





    for c in range(0, pNx - 2):
        bp[c] = rho * Au[c] / aP[c] * aa[c] - rho * Au[c + 1] / aP[c + 1] * aa[c + 1]



    Cp.data = aPpp[1:-1]
    Cp.data = np.append(Cp.data, -aWpp[2:pNx -1])
    Cp.data = np.append(Cp.data, -aEpp[1:pNx -2])

    ppc = spslinalg.bicgstab(Cp, bp, tol=1e-6)[0]

    for c in range(1, pNx - 1):
        p[c] = p[c] + alphap * ppc[c - 1]


    for c in range(1, uNx - 1):
        uuu[c] = uu[c] + d[c] * (ppc[c - 1] - ppc[c]) + aa[c]/aP[c]
    uuu[0] = uu[0] + d[0] * (0 - ppc[0])+ aa[0]/aP[0]
    uuu[-1] = uu[-1] + d[-1] * (ppc[-1] - 0)+ aa[-1]/aP[-1]



    u = uuu.copy()
    # print(it,":",uuu)
    # print()
    diff = rho * bla * Au - rho * u *Au
    if abs(np.max(diff)) < 10**-5:
        print("fertig")
        print(it)
        break
    bla = uuu.copy()
