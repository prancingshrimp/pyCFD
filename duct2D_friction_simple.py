#!/usr/bin/env python3.3

import numpy as np
import scipy.sparse as sps

from scipy import linalg
from scipy.sparse import linalg as spslinalg
import matplotlib.pyplot as plt


np.set_printoptions(linewidth=np.nan)

def SetGridPoints(options):
    if options == 0:
        return 4, 3
        # return 3, 3
    if options == 1:
        return 10, 15
    if options == 2:
        return 40, 21
    if options == 3:
        return 225, 61
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
        return 0.6, 0.7
    if options == 3:
        return 0.5, 0.5
    if options == 4:
        return 0.3, 0.3

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

def mplot(n):
    plt.subplot(3, 1, 1)
    plt.pcolor(u0)
    plt.colorbar()
    plt.subplot(3, 1, 2)
    plt.pcolor(v0)
    plt.colorbar()
    plt.subplot(3, 1, 3)
    plt.pcolor(p0)
    plt.colorbar()
    plt.savefig("fig_" + str(n) + ".png")
    # wait = input("PRESS ENTER TO CONTINUE.")
    plt.close()

def nplot(field):
    plt.pcolor(field)
    plt.colorbar()
    plt.savefig("fig.png")
    # wait = input("PRESS ENTER TO CONTINUE.")
    plt.close()


################################################################################

iterNumber = 5000
gridOpt = 3
relaxOpt = 2

mu = 17.218*10**-6
rho = 1
L = 2.0
h = 0.06

alphaUV, alphaP = SetRelaxation(relaxOpt)
pNx, pNy = SetGridPoints(gridOpt)

uNx = pNx
uNy = pNy
vNx = pNx
vNy = pNy+1


p0 = CreateFields(1, pNy, pNx)[0]
p0_ = p0.reshape(1, pNy*pNx)[0,:]
u0 = CreateFields(1, uNy, uNx)[0]
u0_ = u0.reshape(1, uNy*uNx)[0,:]
v0 = CreateFields(1, vNy, vNx)[0]
v0_ = v0.reshape(1, vNy*vNx)[0,:]

p_wall = CreateFields(1, 1, pNx)[0]
p_out = CreateFields(1, pNy, 1)[0]
u_in = CreateFields(1, uNy, 1)[0]

dx = L/(pNx+1)
dy = h/pNy

Ax = dx
Ay = dy

u_in.fill(0.05)
for c in range(0, uNx):
    u0[:,c] = u_in[:,0] - u_in[:,0]/(uNx+1) * (c+1)
u0.fill(0.05)

p_out.fill(0)
p_wall.fill(0)


# u0[:,0] = 0.1
# u0[:,1] = 0.2
# u0[:,2] = 0.3
# u0[:,3] = 0.4

# u0[:,0] = -0.1
# u0[:,1] = -0.2
# u0[:,2] = -0.3
# u0[:,3] = -0.4

# v0[0,:] = 0.1
# v0[1,:] = 0.2
# v0[2,:] = 0.3
# v0[3,:] = 0.4

# v0[:,0] = -0.1
# v0[:,1] = -0.2
# v0[:,2] = -0.3
# v0[:,3] = -0.4
# v0[:,4] = -0.5

P, aPp, aWp, aEp, aSp, aNp = CreateMatrix(pNy, pNx)
U, aPu, aWu, aEu, aSu, aNu = CreateMatrix(uNy, uNx)
V, aPv, aWv, aEv, aSv, aNv = CreateMatrix(vNy, vNx)

Fxu = CreateFields(1, uNy, uNx+1)[0]
Fyu = CreateFields(1, uNy+1, uNx)[0]
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

u0_t = CreateFields(1, uNy, uNx)[0]
v0_t = CreateFields(1, vNy, vNx)[0]
p0_t = CreateFields(1, pNy, pNx)[0]

Dxuv = mu * dy / dx
Dyuv = mu * dx / dy

switch = 0

np.save("u0_0", u0)
np.save("v0_0", v0)
np.save("p0_0", p0)


###############################################################################
for iteration in range(0, iterNumber):
    print("Current iteration:", iteration)
    # mplot(iteration)

    for r in range(0, uNy):
        for c in range(1, uNx):
            Fxu[r,c] = (u0[r,c-1] + u0[r,c])/2 * rho * Ay
        Fxu[r,0] = (u_in[r] + u0[r,0])/2 * rho * Ay
        Fxu[r,-1] = u0[r,-1] * rho * Ay

    for r in range(0, uNy+1):
        for c in range(0, uNx-1):
            Fyu[r,c] = (v0[r,c] + v0[r,c+1])/2 * rho * Ax
        Fyu[r,-1] = v0[r,-1] * rho * Ax

    for r in range(0, uNy):
        for c in range(0, uNx-1):
            # aWu[r,c] = - np.max([ Fxu[r,c+1], 0])
            aWu[r,c] = - np.max([ Fxu[r,c+1], Dxuv+Fxu[r,c+1]/2, 0])
            # aEu[r,c] = - np.max([-Fxu[r,c+1], 0])
            aEu[r,c] = - np.max([-Fxu[r,c+1], Dxuv-Fxu[r,c+1]/2, 0])

    for r in range(0, uNy-1):
        for c in range(0, uNx):
            # aSu[r,c] = - np.max([ Fyu[r,c], 0])
            aSu[r,c] = - np.max([ Fyu[r,c], Dyuv+Fyu[r,c]/2, 0])
            # aNu[r,c] = - np.max([-Fyu[r,c], 0])
            aNu[r,c] = - np.max([-Fyu[r,c], Dyuv-Fyu[r,c]/2, 0])

    aWup[:,  1:] = -aWu.copy()
    aEup[:, :-1] = -aEu.copy()
    aSup[:-1, :] = -aSu.copy()
    aNup[ 1:, :] = -aNu.copy()
    Fwu = Fxu.copy()[:,:-1]
    Feu = Fxu.copy()[:, 1:]
    Fsu[:,:] = Fyu.copy()[1:,:]
    Fnu[:,:] = Fyu.copy()[:-1,:]

    Spu[:,0] = -Fxu.copy()[:,0]
    Spu[0,:] = -mu/dy*dx
    Spu[-1,:] = -mu/dy*dx
    # Spu[:,-1] = -Fxu.copy()[:,-1]

    for r in range(0, uNy):
        for c in range(0, uNx):
            aPu[r,c] = ( aWup[r,c] + aEup[r,c] + aSup[r,c] + aNup[r,c]
                         - Fwu[r,c] + Feu[r,c] - Fsu[r,c] + Fnu[r,c] - Spu[r,c])
            aPu[r,c] = aPu[r,c]/alphaUV
            du[r,c] = Ay/aPu[r,c]

    for r in range(0, uNy):
        for c in range(0, uNx-1):
            bu[r,c] = (p0[r,c] - p0[r,c+1]) * Ay + (1 - alphaUV) * aPu[r,c] * u0[r,c]
        bu[r,-1] = (p0[r,-1] - p_out[r,0]) * Ay + (1 - alphaUV) * aPu[r,-1] * u0[r,-1]
        bu[r,0] = bu[r,0] + (u_in[r,0] * rho * Ay) * u_in[r,0]


    xu = spslinalg.bicgstab(U, bu_, x0=u0_, tol=1e-12)[0].reshape((uNy,uNx))


################################################################################

    # u0[:,0] = 0.1
    # u0[:,1] = 0.2
    # u0[:,2] = 0.3



    for r in range(1, vNy-1):
        for c in range(1, vNx+1):
            Fxv[r,c] = (u0[r-1,c-1] + u0[r,c-1])/2 * rho * Ay
        Fxv[r,0] = (u_in[r-1] + u_in[r])/2 * rho * Ay
    Fxv[0,0] = (u_in[0] + 0)/2 * rho * Ay
    Fxv[-1,0] = (u_in[-1] + 0)/2 * rho * Ay
    Fxv[0,1:] = u0[0,:]/2 * rho * Ay
    Fxv[-1,1:] = u0[-1,:]/2 * rho * Ay



    # v0[0,:] = 0.1
    # v0[1,:] = 0.2
    # v0[2,:] = 0.3
    # v0[3,:] = 0.4

    for r in range(0, vNy-1):
        for c in range(0, vNx):
            Fyv[r,c] = (v0[r,c] + v0[r+1,c])/2 * rho * Ax

#     # Fxv[0,2] = 34

    for r in range(0, vNy):
        for c in range(0, vNx-1):
            # aWv[r,c] = - np.max([ Fxv[r,c+1], 0])
            aWv[r,c] = - np.max([ Fxv[r,c+1], Dxuv+Fxv[r,c+1]/2, 0])
            # aEv[r,c] = - np.max([-Fxv[r,c+1], 0])
            aEv[r,c] = - np.max([-Fxv[r,c+1], Dxuv-Fxv[r,c+1]/2, 0])

    for r in range(0, vNy-1):
        for c in range(0, vNx):
            # aSv[r,c] = - np.max([ Fyv[r,c], 0])
            aSv[r,c] = - np.max([ Fyv[r,c], Dyuv+Fyv[r,c]/2, 0])
            # aNv[r,c] = - np.max([-Fyv[r,c], 0])
            aNv[r,c] = - np.max([-Fyv[r,c], Dyuv-Fyv[r,c]/2, 0])


    aWvp[:,  1:] = -aWv.copy()
    aEvp[:, :-1] = -aEv.copy()
    aSvp[:-1, :] = -aSv.copy()
    aNvp[ 1:, :] = -aNv.copy()
    Fwv = Fxv.copy()[:,:-1]
    Fev = Fxv.copy()[:,1:]
    Fsv[:-1, :] = Fyv.copy()
    Fnv[1:, :] = Fyv.copy()

    # aSvp[-1,:].fill(Dyuv)
    # aNvp[0,:].fill(Dyuv)

    Spv[:,0] = -Fxv.copy()[:,0]
    # Spv[:,-1] = -Fxv.copy()[:,-1]

    for r in range(0, vNy):
        for c in range(0, vNx):
            aPv[r,c] = ( aWvp[r,c] + aEvp[r,c] + aSvp[r,c] + aNvp[r,c]
                         - Fwv[r,c] + Fev[r,c] - Fsv[r,c] + Fnv[r,c] - Spv[r,c])
            aPv[r,c] = aPv[r,c]/alphaUV
            dv[r,c] = Ax/aPv[r,c]

    # p0[0,:] = 0.1
    # p0[1,:] = 0.2
    # p0[2,:] = 0.3
    # p0[3,:] = 0.4


    for r in range(1, vNy-1):
        for c in range(0, vNx):
            bv[r,c] = (p0[r,c] - p0[r-1,c]) * Ax + (1 - alphaUV) * aPv[r,c] * v0[r,c]
        # bv[r,0] = bv[r,0] + Fxv[r,0] * u_in[r,0]

    bv[0,:] = (p0[0,:] - p_wall[0,:]) * Ax + (1 - alphaUV) * aPv[0,:] * v0[0,:]
    bv[-1,:] = (p_wall[-1,:]- p0[-1,:]) * Ax + (1 - alphaUV) * aPv[-1,:] * v0[-1,:]



    xv = spslinalg.bicgstab(V, bv_, x0=v0_, tol=1e-12)[0].reshape((vNy,vNx))

################################################################################

    # dv[0,1] = 99

    for r in range(0, pNy):
        for c in range(0, pNx-1):
            aWp[r,c] = - du[r,c] * rho * Ay

    for r in range(0, pNy):
        for c in range(0, pNx-1):
            aEp[r,c] = - du[r,c] * rho * Ay

    for r in range(0, pNy-1):
        for c in range(0, pNx):
            aSp[r,c] = - dv[r+1,c] * rho * Ax
            aNp[r,c] = - dv[r+1,c] * rho * Ax


    aWpp[:, 1:] = -aWp.copy()
    aEpp[:, :-1] = -aEp.copy()
    aEpp[:,-1] = du[:,-1] * rho * Ay
    aSpp[:-1, :] = -aSp.copy()
    aSpp[-1, :] = dv[-1,:] * rho * Ax
    aNpp[ 1:, :] = -aNp.copy()
    aNpp[0, :] = dv[0,:] * rho * Ax


    for r in range(0, pNy):
        for c in range(0, pNx):
            aPp[r,c] = aWpp[r,c] + aEpp[r,c] + aSpp[r,c] + aNpp[r,c]

    for r in range(0, pNy):
        for c in range(1, pNx):
            bp[r,c] = xu[r,c-1] * rho * Ay - xu[r,c] * rho * Ay
        bp[r,0] = u_in[r,0] * rho * Ay - xu[r,0] * rho * Ay

    # bp.fill(0)
    for r in range(0, pNy):
        for c in range(0, pNx):
            bp[r,c] = bp[r,c] + xv[r+1,c] * rho * Ax - xv[r,c] * rho * Ax
    # bp[0,:] = bp[0,:] + xv[0,:] * rho * Axv[0,:]
    # bp[-1,:] = bp[-1,:] - xv[-1,:] * rho * Axv[-1,:]

    xp = spslinalg.bicgstab(P, bp_, x0=p0_, tol=1e-12)[0].reshape((pNy,pNx))


    for r in range(0, pNy):
        for c in range(0, pNx):
            p0[r,c] = p0[r,c] + alphaP * xp[r,c]

    for r in range(0, uNy):
        for c in range(0, uNx-1):
            u0[r,c] = xu[r,c] + du[r,c] * (xp[r,c] - xp[r,c+1])
        u0[r,-1] = xu[r,-1] + du[r,-1] * xp[r,-1]


    for r in range(1, vNy-1):
        for c in range(0, vNx):
            v0[r,c] = xv[r,c] + dv[r,c] * (xp[r,c] - xp[r-1,c])
    v0[0,:] = xv[0,:] + dv[0,:] * xp[0,:]
    v0[-1,:] = xv[-1,:] + dv[-1,:] * -xp[-1,:]


    # diff = rho * test * Ay - rho * u0 *Ay
    # if abs(np.max(diff)) < 10**-10:
    #     diff2 = test2 - p0
    #     if abs(np.max(diff2)) < 10**-12:
    #         print("fertig")
    #         print(iteration)
    #         # mplot()
    #         break
    # test = u0.copy()
    # test2 = p0.copy()

    if switch == 10:
        switch = 0
        np.save("u0_" + str(iteration), u0)
        np.save("v0_" + str(iteration), v0)
        np.save("p0_" + str(iteration), p0)
    switch = switch + 1

    # if abs(np.max(u0-u0_t)) < 10**-6:
    #     if abs(np.max(v0-v0_t)) < 10**-6:
    #         if abs(np.max(p0-p0_t)) < 10**-6:
    #             print("fertig")
    #             print(iteration)
    #             # mplot()
    #             break

    u0_t = u0.copy()
    v0_t = v0.copy()
    p0_t = p0.copy()


