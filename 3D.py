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
        return 40, 20
    if options == 3:
        return 200, 60
    if options == 4:
        return 400, 120

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

def SaveFields(n, fields, fieldnames, pathname):
    for item in range(0, len(fields)):
        np.save(pathname + "/" + fieldnames[item] + "_" + str(n), fields[item])



################################################################################
# def main():
iterNumber = 1
gridOpt = 0
relaxOpt = 1

mu = 17.218*10**-6
rho = 1.2758
L = 2.0
h = 0.06

alphaUV, alphaP = SetRelaxation(relaxOpt)
pNx, pNy = SetGridPoints(gridOpt)

uNx = pNx
uNy = pNy
vNx = pNx
vNy = pNy+1

p0, p1, p2 = CreateFields(3, pNy, pNx)
p0_ = p0.reshape(1, pNy*pNx)[0,:]
p1_ = p1.reshape(1, pNy*pNx)[0,:]
u0, u1, u2 = CreateFields(3, uNy, uNx)
u0_ = u0.reshape(1, uNy*uNx)[0,:]
v0, v1, v2 = CreateFields(3, vNy, vNx)
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


SaveFields(0, [u0, v0, p0], ["u0", "v0", "p0"], "./pisoTest")

aau = CreateFields(1, uNy, uNx)[0]
aav = CreateFields(1, vNy, vNx)[0]
bp1 = CreateFields(1, pNy, pNx)[0]


###############################################################################
for iteration in range(0, iterNumber):
    print("Current iteration:", iteration)
    # mplot(iteration)

    for r in range(0, uNy):
        for c in range(1, uNx):
    # for c in range(1, uNx):
    #     for r in range(0, uNy):
            Fxu[r,c] = (u0[r,c-1] + u0[r,c])/2 * rho * Ay
    Fxu[:,0] = (u_in[:,0] + u0[:,0])/2 * rho * Ay
    Fxu[:,-1] = u0[:,-1] * rho * Ay

    for r in range(0, uNy+1):
        for c in range(0, uNx-1):
    # for c in range(0, uNx-1):
    #     for r in range(0, uNy+1):
            Fyu[r,c] = (v0[r,c] + v0[r,c+1])/2 * rho * Ax
    Fyu[:,-1] = v0[:,-1] * rho * Ax


    # blaWu = np.maximum( Fxu, Dxuv+Fxu/2, testFxu)
    # blaEu = np.maximum(-Fxu, Dxuv-Fxu/2, testFxu)
    # blaSu = np.maximum( Fyu, Dyuv+Fyu/2, testFyu)
    # blaNu = np.maximum(-Fyu, Dyuv-Fyu/2, testFyu)

    for r in range(0, uNy):
        for c in range(0, uNx-1):
    # for c in range(0, uNx-1):
    #     for r in range(0, uNy):
            # aWu[r,c] = - blaWu[r,c+1]
            # aEu[r,c] = - blaEu[r,c+1]
            aWu[r,c] = 0
            if  Fxu[r,c+1] > Dxuv+Fxu[r,c+1]/2 and  Fxu[r,c+1] > 0:
                aWu[r,c] = - Fxu[r,c+1]
            if Dxuv+Fxu[r,c+1]/2 > Fxu[r,c+1] and Dxuv+Fxu[r,c+1]/2 > 0:
                aWu[r,c] = - (Dxuv+Fxu[r,c+1]/2)

            aEu[r,c] = 0
            if -Fxu[r,c+1] > Dxuv-Fxu[r,c+1]/2 and -Fxu[r,c+1] > 0:
                aEu[r,c] =   Fxu[r,c+1]
            if Dxuv-Fxu[r,c+1]/2 > -Fxu[r,c+1] and Dxuv-Fxu[r,c+1]/2 > 0:
                aEu[r,c] = - (Dxuv-Fxu[r,c+1]/2)

            # aWu[r,c] = - np.max([ Fxu[r,c+1], Dxuv+Fxu[r,c+1]/2, 0])
            # aEu[r,c] = - np.max([-Fxu[r,c+1], Dxuv-Fxu[r,c+1]/2, 0])



    for r in range(0, uNy-1):
        for c in range(0, uNx):
    # for c in range(0, uNx):
    #     for r in range(0, uNy-1):
            # aSu[r,c] = - blaSu[r,c]
            # aNu[r,c] = - blaNu[r,c]
            aSu[r,c] = 0
            if Fyu[r,c] > Dyuv+Fyu[r,c]/2 and Fyu[r,c] > 0:
                aSu[r,c] = - Fyu[r,c]
            if Dyuv+Fyu[r,c]/2 > Fyu[r,c] and Dyuv+Fyu[r,c]/2 > 0:
                aSu[r,c] = - (Dyuv+Fyu[r,c]/2)

            aNu[r,c] = 0
            if -Fyu[r,c] > Dyuv-Fyu[r,c]/2 and -Fyu[r,c] > 0:
                aNu[r,c] =  Fyu[r,c]
            if Dyuv-Fyu[r,c]/2 > -Fyu[r,c] and Dyuv-Fyu[r,c]/2 > 0:
                aNu[r,c] = - (Dyuv-Fyu[r,c]/2)

            # aSu[r,c] = - np.max([ Fyu[r,c], Dyuv+Fyu[r,c]/2, 0])
            # aNu[r,c] = - np.max([-Fyu[r,c], Dyuv-Fyu[r,c]/2, 0])



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
    # for c in range(0, uNx):
    #     for r in range(0, uNy):
            aPu[r,c] = ( aWup[r,c] + aEup[r,c] + aSup[r,c] + aNup[r,c]
                         - Fwu[r,c] + Feu[r,c] - Fsu[r,c] + Fnu[r,c] - Spu[r,c])
            aPu[r,c] = aPu[r,c]/alphaUV
            du[r,c] = Ay/aPu[r,c]

    for r in range(0, uNy):
        for c in range(0, uNx-1):
    # for c in range(0, uNx-1):
    #     for r in range(0, uNy):
            bu[r,c] = (p0[r,c] - p0[r,c+1]) * Ay + (1 - alphaUV) * aPu[r,c] * u0[r,c]
    bu[:,-1] = (p0[:,-1] - p_out[:,0]) * Ay + (1 - alphaUV) * aPu[:,-1] * u0[:,-1]
    bu[:,0] = bu[:,0] + (u_in[:,0] * rho * Ay) * u_in[:,0]

    # print(aPu[8,7]*u0[8,7]-aWup[8,7]*u0[7,7]-aEup[8,7]*u0[9,7]-aSup[8,7]*u0[8,8]-aNup[8,7]*u0[8,6]-Spu[8,7])

    xu = spslinalg.bicgstab(U, bu_, x0=u0_, tol=1e-8)[0].reshape((uNy,uNx))



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


    # blaWv = np.maximum( Fxv, Dxuv+Fxv/2, testFxv)
    # blaEv = np.maximum(-Fxv, Dxuv-Fxv/2, testFxv)
    # blaSv = np.maximum( Fyv, Dyuv+Fyv/2, testFyv)
    # blaNv = np.maximum(-Fyv, Dyuv-Fyv/2, testFyv)



    for r in range(0, vNy):
        for c in range(0, vNx-1):
            # aWv[r,c] = - blaWv[r,c+1]
            # aEv[r,c] = - blaEv[r,c+1]
            aWv[r,c] = 0
            if  Fxv[r,c+1] > Dxuv+Fxv[r,c+1]/2 and  Fxv[r,c+1] > 0:
                aWv[r,c] = - Fxv[r,c+1]
            if Dxuv+Fxv[r,c+1]/2 > Fxv[r,c+1] and Dxuv+Fxv[r,c+1]/2 > 0:
                aWv[r,c] = - (Dxuv+Fxv[r,c+1]/2)

            aEv[r,c] = 0
            if -Fxv[r,c+1] > Dxuv-Fxv[r,c+1]/2 and -Fxv[r,c+1] > 0:
                aEv[r,c] =   Fxv[r,c+1]
            if Dxuv-Fxv[r,c+1]/2 > -Fxv[r,c+1] and Dxuv-Fxv[r,c+1]/2 > 0:
                aEv[r,c] = - (Dxuv-Fxv[r,c+1]/2)

            # aWv[r,c] = - np.max([ Fxv[r,c+1], Dxuv+Fxv[r,c+1]/2, 0])
            # aEv[r,c] = - np.max([-Fxv[r,c+1], Dxuv-Fxv[r,c+1]/2, 0])


    for r in range(0, vNy-1):
        for c in range(0, vNx):
            # aSv[r,c] = - blaSv[r,c]
            # aNv[r,c] = - blaNv[r,c]
            aSv[r,c] = 0
            if Fyv[r,c] > Dyuv+Fyv[r,c]/2 and Fyv[r,c] > 0:
                aSv[r,c] = - Fyv[r,c]
            if Dyuv+Fyv[r,c]/2 > Fyv[r,c] and Dyuv+Fyv[r,c]/2 > 0:
                aSv[r,c] = - (Dyuv+Fyv[r,c]/2)

            aNv[r,c] = 0
            if -Fyv[r,c] > Dyuv-Fyv[r,c]/2 and -Fyv[r,c] > 0:
                aNv[r,c] =  Fyv[r,c]
            if Dyuv-Fyv[r,c]/2 > -Fyv[r,c] and Dyuv-Fyv[r,c]/2 > 0:
                aNv[r,c] = - (Dyuv-Fyv[r,c]/2)

            # aSv[r,c] = - np.max([ Fyv[r,c], Dyuv+Fyv[r,c]/2, 0])
            # aNv[r,c] = - np.max([-Fyv[r,c], Dyuv-Fyv[r,c]/2, 0])




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



    xv = spslinalg.bicgstab(V, bv_, x0=v0_, tol=1e-8)[0].reshape((vNy,vNx))

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

    xp = spslinalg.bicgstab(P, bp_, x0=p0_, tol=1e-8)[0].reshape((pNy,pNx))


    for r in range(0, pNy):
        for c in range(0, pNx):
            p1[r,c] = p0[r,c] + alphaP * xp[r,c]
            # p0[r,c] = p0[r,c] + alphaP * xp[r,c]

    for r in range(0, uNy):
        for c in range(0, uNx-1):
            u1[r,c] = xu[r,c] + du[r,c] * (xp[r,c] - xp[r,c+1])
            # u0[r,c] = xu[r,c] + du[r,c] * (xp[r,c] - xp[r,c+1])
        u1[r,-1] = xu[r,-1] + du[r,-1] * xp[r,-1]
        # u0[r,-1] = xu[r,-1] + du[r,-1] * xp[r,-1]


    for r in range(1, vNy-1):
        for c in range(0, vNx):
            v1[r,c] = xv[r,c] + dv[r,c] * (xp[r,c] - xp[r-1,c])
            # v0[r,c] = xv[r,c] + dv[r,c] * (xp[r,c] - xp[r-1,c])
    v1[0,:] = xv[0,:] + dv[0,:] * xp[0,:]
    # v0[0,:] = xv[0,:] + dv[0,:] * xp[0,:]
    v1[-1,:] = xv[-1,:] + dv[-1,:] * -xp[-1,:]
    # v0[-1,:] = xv[-1,:] + dv[-1,:] * -xp[-1,:]


################################################################################

    for r in range(1, uNy-1):
        for c in range(1, uNx-1):
            aau[r,c] = ( aWup[r,c]*(u1[r,c-1]-xu[r,c-1]) +
                         aEup[r,c]*(u1[r,c+1]-xu[r,c+1]) +
                         aSup[r,c]*(u1[r+1,c]-xu[r+1,c]) +
                         aNup[r,c]*(u1[r-1,c]-xu[r-1,c])
                       )

    aau[ 0, 0] = ( aEup[ 0, 0]*( u1[ 0, 1] - xu[ 0, 1]) +
                   aSup[ 0, 0]*( u1[ 1, 0] - xu[ 1, 0])
                 )
    aau[ 0,-1] = ( aWup[ 0,-1]*( u1[ 0,-2] - xu[ 0,-2]) +
                   aSup[ 0,-1]*( u1[ 1,-1] - xu[ 1,-1])
                 )
    aau[-1, 0] = ( aEup[ -1,0]*( u1[-1, 1] - xu[-1, 1]) +
                   aNup[ -1,0]*( u1[-2, 0] - xu[-2, 0])
                 )
    aau[-1,-1] = ( aWup[-1,-1]*( u1[-1,-2] - xu[-1,-2]) +
                   aNup[-1,-1]*( u1[-2,-1] - xu[-2,-1])
                 )

    for c in range(1, uNx-1):
        aau[0,c] = ( aWup[0,c]*(u1[0,c-1]-xu[0,c-1]) +
                     aEup[0,c]*(u1[0,c+1]-xu[0,c+1]) +
                     aSup[0,c]*(u1[1,c]-xu[1,c])
                   )
        aau[-1,c] = ( aWup[-1,c]*(u1[-1,c-1]-xu[-1,c-1]) +
                      aEup[-1,c]*(u1[-1,c+1]-xu[-1,c+1]) +
                      aNup[-1,c]*(u1[-2,c]-xu[-2,c])
                   )

    for r in range(1, uNy-1):
        aau[r,0] = ( aEup[r,0]*(u1[r,1]-xu[r,1]) +
                     aSup[r,0]*(u1[r+1,0]-xu[r+1,0]) +
                     aNup[r,0]*(u1[r-1,0]-xu[r-1,0])
                   )
        aau[r,-1] = ( aWup[r,-1]*(u1[r,-2]-xu[r,-2]) +
                      aSup[r,-1]*(u1[r+1,-1]-xu[r+1,-1]) +
                      aNup[r,-1]*(u1[r-1,-1]-xu[r-1,-1])
                   )



    for r in range(1, vNy-1):
        for c in range(1, vNx-1):
            aav[r,c] = ( aWvp[r,c]*(v1[r,c-1]-xv[r,c-1]) +
                         aEvp[r,c]*(v1[r,c+1]-xv[r,c+1]) +
                         aSvp[r,c]*(v1[r+1,c]-xv[r+1,c]) +
                         aNvp[r,c]*(v1[r-1,c]-xv[r-1,c])
                       )

    aav[0,0] =  ( aEvp[0,0]*( v1[0,1]-  xv[0,1]) +
                  aSvp[0,0]*( v1[1,0]-  xv[1,0])
                )

    aav[0,-1] = ( aWvp[0,-1]*( v1[0,-2]-  xv[0,-2]) +
                  aSvp[0,-1]*( v1[1,-1]-  xv[1,-1])
                )

    aav[-1,0] = ( aEvp[-1,0]*( v1[-1,1]-  xv[-1,1])  +
                  aNvp[-1,0]*( v1[-2,0]-  xv[-2,0])
                )

    aav[-1,-1] = ( aWvp[-1,-1]*( v1[-1,-2]-  xv[-1,-2])  +
                   # aEvp[-1,-1]*( v1[-1,-1]-  xv[-1,-1])  +
                   aNvp[-1,-1]*( v1[-2,-1]-  xv[-2,-1])
                 )

    for c in range(1, vNx-1):
        aav[0,c] = ( aWvp[0,c]*(v1[0,c-1]-xv[0,c-1]) +
                     aEvp[0,c]*(v1[0,c+1]-xv[0,c+1]) +
                     aSvp[0,c]*(v1[1,c]-xv[1,c])
                   )
        aav[-1,c] = ( aWvp[-1,c]*(v1[-1,c-1]-xv[-1,c-1]) +
                      aEvp[-1,c]*(v1[-1,c+1]-xv[-1,c+1]) +
                      aNvp[-1,c]*(v1[-2,c]-xv[-2,c])
                   )

    for r in range(1, vNy-1):
        aav[r,0] = ( aEvp[r,0]*(v1[r,1]-xv[r,1]) +
                     aSvp[r,0]*(v1[r+1,0]-xv[r+1,0]) +
                     aNvp[r,0]*(v1[r-1,0]-xv[r-1,0])
                   )
        aav[r,-1] = ( aWvp[r,-1]*(v1[r,-2]-xv[r,-2]) +
                      # aEvp[r,-1]*(v1[r,-1]-xv[r,-1]) +
                      aSvp[r,-1]*(v1[r+1,-1]-xv[r+1,-1]) +
                      aNvp[r,-1]*(v1[r-1,-1]-xv[r-1,-1])
                   )

    for r in range(0, pNy):
        for c in range(0, pNx):
            bp1[r,c] = (rho*Ax/aPv[r+1,c]*aav[r+1,c] - rho*Ax/aPv[r,c]*aav[r,c] )

    for r in range(1, pNy-1):
        for c in range(1, pNx):
            bp1[r,c] = (bp1[r,c] + rho*Ay/aPu[r,c-1]*aau[r,c-1] - rho*Ay/aPu[r,c]*aau[r,c] )

    bp1[:,0] = bp1[:,0] - rho*Ay/aPu[:,0]*aau[:,0]
    for c in range(1, pNx):
        bp1[0,c] = bp1[0,c] + ( rho*Ay/aPu[0,c-1]*aau[0,c-1] - rho*Ay/aPu[0,c]*aau[0,c] )
        bp1[-1,c] = bp1[-1,c] +( rho*Ay/aPu[-1,c-1]*aau[-1,c-1] - rho*Ay/aPu[-1,c]*aau[-1,c] )


    bp1_  =   bp1.reshape(1, pNy*pNx)[0,:]

    xpp = spslinalg.bicgstab(P, bp1_, x0=p1_, tol=1e-8)[0].reshape((pNy,pNx))


    for r in range(0, pNy):
        for c in range(0, pNx):
            p2[r,c] = p1[r,c] + alphaP * xpp[r,c]




    for r in range(0, uNy):
        for c in range(0, uNx-1):
            u2[r,c] = u1[r,c] + du[r,c] * (xpp[r,c] - xpp[r,c+1]) + aau[r,c]/aPu[r,c]
        u2[r,-1] = u1[r,-1] + du[r,-1] * (xpp[r,-1] - 0) + aau[r,-1]/aPu[r,-1]


    for r in range(1, vNy-1):
        for c in range(0, vNx):
            v2[r,c] = v1[r,c] + dv[r,c] * (xpp[r,c] - xpp[r-1,c]) + aav[r,c]/aPv[r,c]
    v2[0,:] = v1[0,:] + dv[0,:] * xpp[0,:] + aav[0,:]/aPv[0,:]
    v2[-1,:] = v1[-1,:] + dv[-1,:] * -xpp[-1,:] + aav[-1,:]/aPv[-1,:]


    u0 = u2.copy()
    v0 = v2.copy()
    p0 = p2.copy()
    # u0 = u1.copy()
    # v0 = v1.copy()
    # p0 = p1.copy()


    if switch == 10:
        switch = 0
        SaveFields(iteration, [u0, v0, p0], ["u0", "v0", "p0"], "./pisoTest")

    if abs(np.max(u0-u0_t)) < 10**-6:
        if abs(np.max(v0-v0_t)) < 10**-6:
            if abs(np.max(p0-p0_t)) < 10**-6:
                print("fertig")
                print(iteration)
                # mplot()
                SaveFields(iteration, [u0, v0, p0], ["u0", "v0", "p0"], "./pisoTest")
                break

    u0_t = u0.copy()
    v0_t = v0.copy()
    p0_t = p0.copy()

    switch = switch + 1
