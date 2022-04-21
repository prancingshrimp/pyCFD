#!/usr/bin/env python3.3

import numpy as np
import scipy.sparse as sps

from scipy import linalg
from scipy.sparse import linalg as spslinalg

np.set_printoptions(linewidth=np.nan)

def SetGridPoints(options):
    if options == 0:
        return 5, 5
        # return 4, 7
    if options == 1:
        return 10, 15
    if options == 2:
        return 60, 40
    if options == 3:
        return 200, 60
    if options == 4:
        return 300, 151
    if options == 5:
        return 1000, 600

def SetRelaxation(options):
    if options == 0:
        return 1.0, 1.0
    if options == 1:
        return 0.8, 0.9
    if options == 2:
        return 0.5, 0.5

def CreateFields(n, Nx, Ny):
    tmp = []
    for item in range(0, n):
        tmp.append(np.zeros([Nx,Ny]))
    return tmp

def CreateMatrix(Ny, Nx):
    num = 5*Nx*Ny-2*Ny-2*Nx
    M = sps.coo_matrix((Ny*Nx, Ny*Nx), dtype=np.float64)

    M.row = np.array(np.zeros(num), dtype=np.intc)
    M.col = np.array(np.zeros(num), dtype=np.intc)
    M.data = np.array(np.zeros(num), dtype=np.float64)

    # aPp
    M.row[0:Ny*Nx] = np.arange(0, Ny*Nx)
    M.col[0:Ny*Nx] = np.arange(0, Ny*Nx)
    # M.data[0:Ny*Nx] = np.arange(1, Ny*Nx+1)

    for r in range(0, Ny):
        # aWp
        M.row [Ny*Nx + r*(Nx-1) : Ny*Nx + (r+1)*(Nx-1)] = np.arange(r*Nx+1,(r+1)*Nx)
        M.col [Ny*Nx + r*(Nx-1) : Ny*Nx + (r+1)*(Nx-1)] = np.arange(r*Nx,(r+1)*Nx-1)
        # M.data[Ny*Nx + r*(Nx-1) : Ny*Nx + (r+1)*(Nx-1)] = np.arange(r*Nx+2,(r+1)*Nx+1)

    # need another for loop because of memory
    for r in range(0, Ny):
        # aEp
        M.row [2*Ny*Nx-Ny + r*(Nx-1) : 2*Ny*Nx-Ny + (r+1)*(Nx-1)] = np.arange(r*Nx,(r+1)*Nx-1)
        M.col [2*Ny*Nx-Ny + r*(Nx-1) : 2*Ny*Nx-Ny + (r+1)*(Nx-1)] = np.arange(r*Nx+1,(r+1)*Nx)
        # M.data[2*Ny*Nx-Ny + r*(Nx-1) : 2*Ny*Nx-Ny + (r+1)*(Nx-1)] = np.arange(r*Nx+1,(r+1)*Nx)

    # aSp
    M.row[3*Ny*Nx-2*Ny:4*Ny*Nx-2*Ny-Nx] = np.arange(0, (Ny-1)*Nx)
    M.col[3*Ny*Nx-2*Ny:4*Ny*Nx-2*Ny-Nx] = np.arange(Nx, Ny*Nx)
    # M.data[3*Ny*Nx-2*Ny:4*Ny*Nx-2*Ny-Nx] = np.arange(1, Ny*Nx-Nx+1)

    # aNp
    M.row[4*Ny*Nx-2*Ny-Nx:5*Ny*Nx-2*Ny-2*Nx] = np.arange(Nx, Ny*Nx)
    M.col[4*Ny*Nx-2*Ny-Nx:5*Ny*Nx-2*Ny-2*Nx] = np.arange(0, (Ny-1)*Nx)
    # M.data[4*Ny*Nx-2*Ny-Nx:5*Ny*Nx-2*Ny-2*Nx] = np.arange(Nx+1, Ny*Nx+1)

    # connect with matrix in memory
    aP = M.data[0:Ny*Nx].reshape(Ny, Nx)
    aW = M.data[Ny*Nx:2*Ny*Nx-Ny].reshape(Ny, Nx-1)
    aE = M.data[2*Ny*Nx-Ny:3*Ny*Nx-2*Ny].reshape(Ny, Nx-1)
    aS = M.data[3*Ny*Nx-2*Ny:4*Ny*Nx-2*Ny-Nx].reshape(Ny-1, Nx)
    aN = M.data[4*Ny*Nx-2*Ny-Nx:5*Ny*Nx-2*Ny-2*Nx].reshape(Ny-1, Nx)

    return M, aP, aW, aE, aS, aN



iterNumber = 1
gridOpt = 0
relaxOpt = 2

mu = 17.218*10**-6
# rho = 1.2758
rho = 1
L = 2.0
mass_flow = 1
AA = 0.5
AE = 0.1
pA = 9.6
pE = 0

alphaUV, alphaP = SetRelaxation(relaxOpt)
pNx, pNy = SetGridPoints(gridOpt)

uNx = pNx - 1
uNy = pNy
vNx = pNx
vNy = pNy -1

# p0, p1, p2 = CreateFields(3, pNy, pNx)
p0, Ap = CreateFields(2, pNy, pNx)
# u0, u1, u2 = CreateFields(3, uNy, uNx)
u0, Au = CreateFields(2, uNy, uNx)
# v0, v1, v2 = CreateFields(3, vNy, vNx)
v0, Av = CreateFields(2, vNy, vNx)


for c in range(0, pNx):
    Ap[:,c] = (AA - (AA - AE)/(pNx - 1) * c)/pNy
for c in range(0, uNx):
    Au[:,c] = (Ap[:,c] + Ap[:,c +1])/2
dx = L/(vNx-1)
for c in range(1, vNx-1):
    Av[:,c] = dx
Av[:,0] = dx/2
Av[:,-1] = dx/2

u0 = mass_flow/(uNy * rho * Au)
for c in range(0, pNx):
    p0[:,c] = pA - (pA - pE)/(pNx - 1) * c

# u0[0,0] = 12
# u0[5,8] = 8



# u0[:,0] = 0.1
# u0[:,1] = 0.2
# u0[:,2] = 0.3
# u0[:,3] = 0.4

# u0[:,0] = -0.1
# u0[:,1] = -0.2
# u0[:,2] = -0.3
# u0[:,3] = -0.4

# v0[:,0] = 0.1
# v0[:,1] = 2
# v0[3,:] = -2
# v0[:,2] = 0.3
# v0[:,3] = 0.4
# v0[:,4] = 0.5

# v0[:,0] = -0.1
# v0[:,1] = -0.2
# v0[:,2] = -0.3
# v0[:,3] = -0.4
# v0[:,4] = -0.5

# print(u0)
# print(v0)

P, aPp, aWp, aEp, aSp, aNp = CreateMatrix(pNy, pNx)
U, aPu, aWu, aEu, aSu, aNu = CreateMatrix(uNy, uNx)
V, aPv, aWv, aEv, aSv, aNv = CreateMatrix(vNy, vNx)

Fxu = CreateFields(1, uNy, uNx+1)[0]
Fyu = CreateFields(1, uNy-1, uNx)[0]
Fxv = CreateFields(1, vNy, vNx+1)[0]
Fyv = CreateFields(1, vNy-1, vNx)[0]

aWup, aEup, aSup, aNup, bu = CreateFields(5, uNy, uNx)
Fwu, Feu, Fsu, Fnu, Spu, du = CreateFields(6, uNy, uNx)
bu_  =  bu.reshape(1, uNy*uNx)[0,:]

aWvp, aEvp, aSvp, aNvp, bv = CreateFields(5, vNy, vNx)
Fwv, Fev, Fsv, Fnv, Spv, dv = CreateFields(6, vNy, vNx)
bv_  =  bv.reshape(1, vNy*vNx)[0,:]

aWpp, aEpp, aSpp, aNpp, bp = CreateFields(5, pNy, pNx)
bp_  =  bp.reshape(1, pNy*(pNx))[0,:]

test = CreateFields(1, uNy, uNx)[0]


###############################################################################
for iteration in range(0, iterNumber):
    print("Current iteration:", iteration)

    for r in range(0, uNy):
        for c in range(1, uNx):
            Fxu[r,c] = (u0[r,c-1] + u0[r,c])/2 * rho * Ap[r,c]
        # Fw[0] is rho * u[0] * Au[0] multiplied with Au[0]/Ap[0]
        # Fxu[r,0] = u0[r,0] * rho * Au[r,0] * Au[r,0]/Ap[r,0]
        Fxu[r,-1] = u0[r,-1] * rho * Au[r,-1]

    for r in range(0, uNy-1):
        for c in range(0, uNx):
            Fyu[r,c] = (v0[r,c] + v0[r,c+1])/2 * rho * dx


    for r in range(0, uNy):
        for c in range(0, uNx-1):
            aWu[r,c] = - np.max([ Fxu[r,c+1], 0])
            # aWu[r,c] = - np.max([ Fxu[r,c], Dxuv+Fxu[r,c]/2, 0])
            aEu[r,c] = - np.max([-Fxu[r,c+1], 0])
            # aEu[r,c] = - np.max([-Fxu[r,c], Dxuv-Fxu[r,c]/2, 0])

    for r in range(0, uNy-1):
        for c in range(0, uNx):
            aSu[r,c] = - np.max([ Fyu[r,c], 0])
            # aSu[r,c] = - np.max([ Fyu[r,c], Dyuv+Fyu[r,c]/2, 0])
            aNu[r,c] = - np.max([-Fyu[r,c], 0])
            # aNu[r,c] = - np.max([-Fyu[r,c], Dyuv-Fyu[r,c]/2, 0])


    aWup[:,  1:] = -aWu.copy()
    aEup[:, :-1] = -aEu.copy()
    aSup[:-1, :] = -aSu.copy()
    aNup[ 1:, :] = -aNu.copy()
    Fwu = Fxu.copy()[:,:-1]
    Feu = Fxu.copy()[:,1:]
    Fsu[:-1, :] = Fyu.copy()
    Fnu[1:, :] = Fyu.copy()

    # Spu[:,0] = -Fxu.copy()[:,0]


    for r in range(0, uNy):
         for c in range(0, uNx):
            aPu[r,c] = ( aWup[r,c] + aEup[r,c] + aSup[r,c] + aNup[r,c]
                         - Fwu[r,c] + Feu[r,c] - Fsu[r,c] + Fnu[r,c] - Spu[r,c])
            aPu[r,c] = aPu[r,c]/alphaUV
            du[r,c] = Au[r,c]/aPu[r,c]
            bu[r,c] = (p0[r,c] - p0[r,c+1]) * Au[r,c] + (1 - alphaUV) * aPu[r,c] * u0[r,c]

    # bu[:,0] = bu[:,0] + (u0[:,0] * rho * Au[:,0] * Au[:,0]/Ap[:,0]) * u0[:,0]

    xu  = spslinalg.bicgstab(U, bu_, tol=1e-8)[0].reshape((uNy,uNx))
    # xu0 = spslinalg.gmres(U, bu_, tol=1e-8)[0].reshape((uNy,uNx))


################################################################################

    for r in range(0, vNy):
        for c in range(1, vNx):
            Fxv[r,c] = (u0[r+1,c-1] + u0[r,c-1])/2 * rho * (Au[r+1,c-1] + Au[r,c-1])/2
        Fxv[r,-1] = (u0[r+1,-1] + u0[r,-1])/2 * rho * (Au[r+1,-1] + Au[r,-1])/2

    # v0[0,:] = 0.1
    # v0[1,:] = 0.2
    # v0[2,:] = 0.3
    # v0[3,:] = 0.4

    for r in range(0, vNy-1):
        for c in range(0, vNx):
            Fyv[r,c] = (v0[r,c] + v0[r+1,c])/2 * rho * Av[r,c]



    for r in range(0, vNy):
        for c in range(0, vNx-1):
            aWv[r,c] = - np.max([ Fxv[r,c+1], 0])
            aEv[r,c] = - np.max([-Fxv[r,c+1], 0])

    for r in range(0, vNy-1):
        for c in range(0, vNx):
            aSv[r,c] = - np.max([ Fyv[r,c], 0])
            aNv[r,c] = - np.max([-Fyv[r,c], 0])

    aWvp[:,  1:] = -aWv.copy()
    aEvp[:, :-1] = -aEv.copy()
    aSvp[:-1, :] = -aSv.copy()
    aNvp[ 1:, :] = -aNv.copy()
    Fwv = Fxv.copy()[:,:-1]
    Fev = Fxv.copy()[:,1:]
    Fsv[:-1, :] = Fyv.copy()
    Fnv[1:, :] = Fyv.copy()


    for r in range(0, vNy):
        for c in range(0, vNx):
            aPv[r,c] = ( aWvp[r,c] + aEvp[r,c] + aSvp[r,c] + aNvp[r,c]
                         - Fwv[r,c] + Fev[r,c] - Fsv[r,c] + Fnv[r,c] - Spv[r,c])
            aPv[r,c] = aPv[r,c]/alphaUV
            dv[r,c] = Av[r,c]/aPv[r,c]
            bv[r,c] = (p0[r+1,c] - p0[r,c]) * Av[r,c] + (1 - alphaUV) * aPv[r,c] * v0[r,c]


    xv = spslinalg.bicgstab(V, bv_, tol=1e-8)[0].reshape((vNy,vNx))
    # xv0 = spslinalg.gmres(V, bv_, tol=1e-8)[0].reshape((vNy,vNx))

################################################################################

    for r in range(0, pNy):
        for c in range(0, pNx-2):
            aWp[r,c] = - du[r,c] * rho * Au[r,c]

    for r in range(0, pNy):
        for c in range(1, pNx-1):
            aEp[r,c] = - du[r,c] * rho * Au[r,c]

    for r in range(0, pNy-1):
        for c in range(1, pNx-1):
            aSp[r,c] = - dv[r,c] * rho * Av[r,c]
            aNp[r,c] = - dv[r,c] * rho * Av[r,c]


    aWpp[:,  1:] = -aWp.copy()
    aEpp[:, :-1] = -aEp.copy()
    aSpp[:-1, :] = -aSp.copy()
    aNpp[ 1:, :] = -aNp.copy()


    for r in range(0, pNy):
        for c in range(0, pNx):
            aPp[r,c] = aWpp[r,c] + aEpp[r,c] + aSpp[r,c] + aNpp[r,c]

    for r in range(0, pNy):
        for c in range(1, pNx-1):
            bp[r,c] = xu[r,c-1] * rho * Au[r,c-1] - xu[r,c] * rho *Au[r,c]

    # xv[0,:] = 0.1
    # xv[1,:] = 0.2
    # xv[2,:] = 0.3
    # xv[3,:] = 0.4

    # bp.fill(0)
    for r in range(1, pNy-1):
        for c in range(1, pNx-1):
            bp[r,c] = bp[r,c] + xv[r,c] * rho * Av[r,c] - xv[r-1,c] * rho *Av[r-1,c]
    bp[0,1:-1] = bp[0,1:-1] + xv[0,1:-1] * rho * Av[0,1:-1]
    bp[-1,1:-1] = bp[-1,1:-1] - xv[-1,1:-1] * rho * Av[-1,1:-1]

    xp = spslinalg.bicgstab(P, bp_, tol=1e-8)[0].reshape((pNy,pNx))
    # xp0 = spslinalg.gmres(P, bp_, tol=1e-8)[0].reshape((pNy,pNx))
    # xp1 = spslinalg.cg(P, bp_, tol=1e-8)[0].reshape((pNy,pNx))

    for r in range(0, pNy):
        for c in range(0, pNx):
            p0[r,c] = p0[r,c] + alphaP * xp[r,c]

    for r in range(0, uNy):
        for c in range(0, uNx):
            u0[r,c] = xu[r,c] + du[r,c] * (xp[r,c] - xp[r,c+1])

    for r in range(0, vNy):
        for c in range(0, vNx):
            v0[r,c] = xv[r,c] + dv[r,c] * (xp[r+1,c] - xp[r,c])


    diff = rho * test * Au - rho * u0 *Au
    if abs(np.max(diff)) < 10**-5:
        print("fertig")
        print(iteration)
        break
    test = u0.copy()
