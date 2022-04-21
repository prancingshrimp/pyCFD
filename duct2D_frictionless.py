#!/usr/bin/env python3.3

import numpy as np
import scipy.sparse as sps

from scipy import linalg
from scipy.sparse import linalg as spslinalg
import matplotlib.pyplot as plt


np.set_printoptions(linewidth=np.nan)

def SetGridPoints(options):
    if options == 0:
        return 3, 5
        # return 3, 3
    if options == 1:
        return 10, 15
    if options == 2:
        return 40, 20
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
    if options == 3:
        return 0.2, 0.2

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

def mplot():

    plt.subplot(3, 1, 1)
    plt.pcolor(u0)
    plt.colorbar()
    plt.subplot(3, 1, 2)
    plt.pcolor(v0)
    plt.colorbar()
    plt.subplot(3, 1, 3)
    plt.pcolor(p0)
    plt.colorbar()
    plt.show()
    wait = input("PRESS ENTER TO CONTINUE.")


################################################################################

iterNumber = 400
gridOpt = 2
relaxOpt = 2

mu = 17.218*10**-6
# rho = 1.2758
rho = 1
L = 2.0
h = 0.06
# mass_flow = 1
# AA = 0.5
# AE = 0.1
# pA = 9.6
# u_in = 10

alphaUV, alphaP = SetRelaxation(relaxOpt)
pNx, pNy = SetGridPoints(gridOpt)

uNx = pNx
uNy = pNy
vNx = pNx+1
vNy = pNy-1

# p0, p1, p2 = CreateFields(3, pNy, pNx)
p0, Axp, Ayp = CreateFields(3, pNy, pNx)
# u0, u1, u2 = CreateFields(3, uNy, uNx)
u0, Axu, Ayu = CreateFields(3, uNy, uNx)
# v0, v1, v2 = CreateFields(3, vNy, vNx)
v0, Axv, Ayv = CreateFields(3, vNy, vNx)

p_out, Axp_out, Ayp_out = CreateFields(3, pNy, 1)
u_in, Axu_in, Ayu_in = CreateFields(3, uNy, 1)



dx = L/(pNx+1)
dy = h/pNy

Axp.fill(dx)
Ayp.fill(dy)
Axp_out.fill(dx)
Ayp_out.fill(dy)
Axu_in.fill(dx)
Ayu_in.fill(dy)
# Axp[0,0] = 1
# Axp[0,1] = 2
# Axp[0,2] = 3
# Axp[1,0] = 1
# Axp[1,1] = 2
# Axp[1,2] = 3
# Axp[2,0] = 1
# Axp[2,1] = 2
# Axp[2,2] = 3
# Ayp = Axp.copy().transpose()

for r in range(0, uNy):
    for c in range(0, uNx-1):
        Axu[r,c] = (Axp[r,c] + Axp[r,c+1])/2
Axu[:,-1] = (Axp[:,-1] + Axp_out[:,0])/2

for r in range(0, uNy):
    for c in range(0, uNx-1):
        Ayu[r,c] = (Ayp[r,c] + Ayp[r,c+1])/2
Ayu[:,-1] = (Ayp[:,-1] + Ayp_out[:,0])/2


for r in range(0, vNy):
    for c in range(0, vNx-1):
        Axv[r,c] = (Axp[r,c] + Axp[r+1,c])/2
    Axv[r,-1] = (Axp_out[r,0] + Axp_out[r+1,0])/2

for r in range(0, vNy):
    for c in range(0, vNx-1):
        Ayv[r,c] = (Ayp[r,c] + Ayp[r+1,c])/2
    Ayv[r,-1] = (Ayp_out[r,0] + Ayp_out[r+1,0])/2



u_in.fill(10)
u0.fill(0.01)
# u0[:,0].fill(10)
# u0[:,1].fill(1)
# u0[:,2].fill(15)
# u0[:,-5].fill(0.1)

p_out.fill(1)

# pA = 6
# for c in range(0, pNx):
#     p0[:,c] = pA - (pA - p_out[:,0])/(pNx) * c
# p0[:,0].fill(5)
# p0[:,0].fill(4)
# p0[:,1].fill(-2)
# p0[:,9].fill(30)

# u0[:,0] = 0.1
# u0[:,1] = 0.2
# u0[:,2] = 0.3
# u0[:,3] = 0.4

# u0[:,0] = -0.1
# u0[:,1] = -0.2
# u0[:,2] = -0.3
# u0[:,3] = -0.4

# v0[:,0] = 0.1
# v0[:,1] = 0.2
# v0[2,2] = 0.0005
# v0[:,3] = 0.4

# v0[:,0] = -0.1
# v0[:,1] = -0.2
# v0[:,2] = -0.3
# v0[:,3] = -0.4
# v0[:,4] = -0.5

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
test2 = CreateFields(1, pNy, pNx)[0]

Dxuv = mu * dy / dx
Dyuv = mu * dx / dy

mplot()
###############################################################################
for iteration in range(0, iterNumber):
    print("Current iteration:", iteration)



    for r in range(0, uNy):
        for c in range(1, uNx):
            Fxu[r,c] = (u0[r,c-1] + u0[r,c])/2 * rho * Ayp[r,c]
        Fxu[r,0] = (u_in[r] + u0[r,0])/2 * rho * Ayp[r,0]
        Fxu[r,-1] = u0[r,-1] * rho * Ayp_out[r,0]

    for r in range(0, uNy-1):
        for c in range(0, uNx):
            Fyu[r,c] = (v0[r,c] + v0[r,c+1])/2 * rho * Axu[r,c]


    for r in range(0, uNy):
        for c in range(0, uNx-1):
            aWu[r,c] = - np.max([ Fxu[r,c+1], 0])
            # aWu[r,c] = - np.max([ Fxu[r,c+1], Dxuv+Fxu[r,c+1]/2, 0])
            aEu[r,c] = - np.max([-Fxu[r,c+1], 0])
            # aEu[r,c] = - np.max([-Fxu[r,c+1], Dxuv-Fxu[r,c+1]/2, 0])

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
    Feu = Fxu.copy()[:, 1:]
    Fsu[:-1, :] = Fyu.copy()
    Fnu[1:, :] = Fyu.copy()

    Spu[:,0] = -Fxu.copy()[:,0]


    for r in range(0, uNy):
        for c in range(0, uNx):
            aPu[r,c] = ( aWup[r,c] + aEup[r,c] + aSup[r,c] + aNup[r,c]
                         - Fwu[r,c] + Feu[r,c] - Fsu[r,c] + Fnu[r,c] - Spu[r,c])
            aPu[r,c] = aPu[r,c]/alphaUV
            du[r,c] = Ayu[r,c]/aPu[r,c]


    for r in range(0, uNy):
        for c in range(0, uNx-1):
            bu[r,c] = (p0[r,c] - p0[r,c+1]) * Ayu[r,c] + (1 - alphaUV) * aPu[r,c] * u0[r,c]
        bu[r,-1] = (p0[r,-1] - p_out[r,0]) * Ayu[r,-1] + (1 - alphaUV) * aPu[r,-1] * u0[r,-1]
        bu[r,0] = bu[r,0] + (u_in[r,0] * rho * Ayu_in[r,0]) * u_in[r,0]


    xu = spslinalg.bicgstab(U, bu_, tol=1e-8)[0].reshape((uNy,uNx))

################################################################################

    # u0[0,:] = 0.1
    # u0[1,:] = 0.2
    # u0[2,:] = 0.3
    # u0[3,:] = 0.4
    # u0[4,:] = 0.5





    for r in range(0, vNy):
        for c in range(1, vNx):
            Fxv[r,c] = (u0[r+1,c-1] + u0[r,c-1])/2 * rho * (Ayu[r+1,c-1] + Ayu[r,c-1])/2
        Fxv[r,-1] = (u0[r+1,-1] + u0[r,-1])/2 * rho * (Ayu[r+1,-1] + Ayu[r,-1])/2
        Fxv[r,0] = (u_in[r+1] + u_in[r])/2 * rho * (Ayu_in[r+1] + Ayu_in[r])/2

    # v0[0,:] = 0.1
    # v0[1,:] = 0.2
    # v0[2,:] = 0.3
    # v0[3,:] = 0.4

    for r in range(0, vNy-1):
        for c in range(0, vNx):
            Fyv[r,c] = (v0[r,c] + v0[r+1,c])/2 * rho * Axv[r,c]

    # Fxv[0,3] = 34


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

    Spv[:,0] = -Fxv.copy()[:,0]

    for r in range(0, vNy):
        for c in range(0, vNx):
            aPv[r,c] = ( aWvp[r,c] + aEvp[r,c] + aSvp[r,c] + aNvp[r,c]
                         - Fwv[r,c] + Fev[r,c] - Fsv[r,c] + Fnv[r,c] - Spv[r,c])
            aPv[r,c] = aPv[r,c]/alphaUV
            dv[r,c] = Axv[r,c]/aPv[r,c]

    for r in range(0, vNy):
        for c in range(0, vNx-1):
            bv[r,c] = (p0[r+1,c] - p0[r,c]) * Axv[r,c] + (1 - alphaUV) * aPv[r,c] * v0[r,c]
        bv[r,-1] = (p_out[r+1] - p_out[r]) * Axv[r,c] + (1 - alphaUV) * aPv[r,-1] * v0[r,-1]
        # bv[r,0] = bv[r,0] + Fxv[r,0] * u_in[r,0]


    xv = spslinalg.bicgstab(V, bv_, tol=1e-8)[0].reshape((vNy,vNx))

#################################################################################




    for r in range(0, pNy):
        for c in range(0, pNx-1):
            aWp[r,c] = - du[r,c] * rho * Ayu[r,c]

    for r in range(0, pNy):
        for c in range(0, pNx-1):
            aEp[r,c] = - du[r,c] * rho * Ayu[r,c]

    for r in range(0, pNy-1):
        for c in range(0, pNx):
            aSp[r,c] = - dv[r,c] * rho * Axv[r,c]
            aNp[r,c] = - dv[r,c] * rho * Axv[r,c]


    aWpp[:, 1:] = -aWp.copy()
    aEpp[:, :-1] = -aEp.copy()
    aEpp[:,-1] = du[:,-1] * rho * Ayu[:,-1]
    aSpp[:-1, :] = -aSp.copy()
    aNpp[ 1:, :] = -aNp.copy()


    for r in range(0, pNy):
        for c in range(0, pNx):
            aPp[r,c] = aWpp[r,c] + aEpp[r,c] + aSpp[r,c] + aNpp[r,c]

    for r in range(0, pNy):
        for c in range(1, pNx):
            bp[r,c] = xu[r,c-1] * rho * Ayu[r,c-1] - xu[r,c] * rho * Ayu[r,c]
        bp[r,0] = u_in[r,0] * rho * Ayu_in[r,0] - xu[r,0] * rho * Ayu[r,0]

    # bp.fill(0)
    for r in range(1, pNy-1):
        for c in range(0, pNx):
            bp[r,c] = bp[r,c] + xv[r,c] * rho * Axv[r,c] - xv[r-1,c] * rho *Axv[r-1,c]
    bp[0,:] = bp[0,:] + xv[0,:-1] * rho * Axv[0,:-1]
    bp[-1,:] = bp[-1,:] +- xv[-1,:-1] * rho * Axv[-1,:-1]

    xp = spslinalg.bicgstab(P, bp_, tol=1e-8)[0].reshape((pNy,pNx))


    for r in range(0, pNy):
        for c in range(0, pNx):
            p0[r,c] = p0[r,c] + alphaP * xp[r,c]

    for r in range(0, uNy):
        for c in range(0, uNx-1):
            u0[r,c] = xu[r,c] + du[r,c] * (xp[r,c] - xp[r,c+1])
        u0[r,-1] = xu[r,-1] + du[r,-1] * (xp[r,-1])


    for r in range(0, vNy):
        for c in range(0, vNx-1):
            v0[r,c] = xv[r,c] + dv[r,c] * (xp[r+1,c] - xp[r,c])
        v0[r,-1] = xv[r,-1] + dv[r,-1] * (p_out[r+1] - p_out[r])


    diff = rho * test * Ayu - rho * u0 *Ayu
    if abs(np.max(diff)) < 10**-10:
        diff2 = test2 - p0
        if abs(np.max(diff2)) < 10**-12:
            print("fertig")
            print(iteration)
            mplot()
            break
    test = u0.copy()
    test2 = p0.copy()



