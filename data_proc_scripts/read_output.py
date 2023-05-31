import numpy as np 
import struct
import matplotlib.pyplot as plt 
from matplotlib import cm
import sys
import glob
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
import array
import os

def read_grid():
    # read nx, xgrid-------------------------
    file_grid = open('./output/xgrid.dat', 'rb')

    nx = ((struct.unpack("N",file_grid.read(8)))[0])

    xgrid = np.zeros(nx)

    for i in range(nx):
        xgrid[i] = (struct.unpack("d",file_grid.read(8)))[0]

    file_grid.close()

    # read ny, ygrid-------------------------
    file_grid = open('./output/ygrid.dat', 'rb')

    ny = ((struct.unpack("N",file_grid.read(8)))[0])

    ygrid = np.zeros(ny)

    for i in range(ny):
        ygrid[i] = (struct.unpack("d",file_grid.read(8)))[0]

    file_grid.close()

    return xgrid, ygrid


def read_output(filename, nvar,nx,ny):
    arr = array.array('d')
    arr.fromfile(open(filename, 'rb'), nvar*nx*ny + 1)

    t = arr[0]
    uu = np.reshape(np.asarray(arr[1:len(arr)]), [nvar,nx,ny], order='C')
    return t, uu


# read grid
xgrid, ygrid = read_grid()
nx = len(xgrid)
ny = len(ygrid)

print('nx = ', nx, ', ny = ', ny)

XX, YY = np.meshgrid(xgrid,ygrid,indexing='ij')

files = sorted(glob.glob('./output/out*.dat'))

nout = len(files)
for nt in range(nout):
    t, uu = read_output(files[nt],8,nx,ny)
    print('Read file: t = {:.3f}'.format(t))

    Bx = uu[4,:,:]
    By = uu[5,:,:]
    Bz = uu[6,:,:]
    ux = uu[1,:,:]
    uy = uu[2,:,:]
    uz = uu[3,:,:]
    rho = uu[0,:,:]
    p = uu[7,:,:]

    # # 1D cut along x---------------------
    # fig = plt.figure(figsize=[10,10])

    # sub = fig.add_subplot(411)
    # sub.plot(xgrid,rho[:,int(ny/2)],label=r'$\rho$')
    # sub.legend(loc='lower left')
    # subt = sub.twinx()
    # subt.plot(xgrid,p[:,int(ny/2)],label=r'$P$',color='C1')
    # subt.spines['right'].set_color('C1')
    # subt.tick_params(axis='y',which='both',colors='C1')
    # subt.legend(loc='lower right')

    # sub = fig.add_subplot(412,sharex=sub)
    # sub.plot(xgrid,ux[:,int(ny/2)],label=r'$u_x$')
    # sub.plot(xgrid,Bx[:,int(ny/2)],label=r'$B_x$',color='C1')
    # sub.legend()

    # sub = fig.add_subplot(413,sharex=sub)
    # sub.plot(xgrid,uy[:,int(ny/2)],label=r'$u_y$')
    # sub.plot(xgrid,By[:,int(ny/2)],label=r'$B_y$',color='C1')
    # sub.legend()

    # sub = fig.add_subplot(414,sharex=sub)
    # sub.plot(xgrid,uz[:,int(ny/2)],label=r'$u_z$')
    # sub.plot(xgrid,Bz[:,int(ny/2)],label=r'$B_z$',color='C1')
    # sub.legend()

    # fig.subplots_adjust(hspace=0)
    # fig.suptitle(r'$t={:.2f}$'.format(t))
    # fig.savefig('./figure/1D_cut_along_x_{:04d}.png'.format(nt))
    # plt.close(fig)




    # # 1D cut along y---------------------
    # fig = plt.figure(figsize=[10,10])

    # sub = fig.add_subplot(411)
    # sub.plot(ygrid,rho[int(nx/2),:],label=r'$\rho$')
    # sub.legend(loc='lower left')
    # subt = sub.twinx()
    # subt.plot(ygrid,p[int(nx/2),:],label=r'$P$',color='C1')
    # subt.spines['right'].set_color('C1')
    # subt.tick_params(axis='y',which='both',colors='C1')
    # subt.legend(loc='lower right')

    # sub = fig.add_subplot(412,sharex=sub)
    # sub.plot(ygrid,ux[int(nx/2),:],label=r'$u_x$')
    # sub.plot(ygrid,Bx[int(nx/2),:],label=r'$B_x$',color='C1')
    # sub.legend()

    # sub = fig.add_subplot(413,sharex=sub)
    # sub.plot(ygrid,uy[int(nx/2),:],label=r'$u_y$')
    # sub.plot(ygrid,By[int(nx/2),:],label=r'$B_y$',color='C1')
    # sub.legend()

    # sub = fig.add_subplot(414,sharex=sub)
    # sub.plot(ygrid,uz[int(nx/2),:],label=r'$u_z$')
    # sub.plot(ygrid,Bz[int(nx/2),:],label=r'$B_z$',color='C1')
    # sub.legend()

    # fig.subplots_adjust(hspace=0)
    # fig.suptitle(r'$t={:.2f}$'.format(t))
    # fig.savefig('./figure/1D_cut_along_y_{:04d}.png'.format(nt))
    # plt.close(fig)



    figsize=[5,12]
    # Bx---
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    cont = sub.pcolormesh(XX,YY,Bx,shading='gouraud',
        cmap=cm.RdBu_r,norm=TwoSlopeNorm(vcenter=0))
    cb = fig.colorbar(cont,shrink=1,aspect=15)
    sub.set_xlabel(r'x',fontsize=15)
    sub.set_ylabel(r'y',fontsize=15)
    sub.set_aspect('equal')
    fig.suptitle(r'$t={:.3f}$'.format(t))
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig('./figure/Bx_{:04d}.png'.format(nt))
    plt.close(fig)

    # By---
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    cont = sub.pcolormesh(XX,YY,By,shading='gouraud',cmap=cm.RdBu_r)
    cb = fig.colorbar(cont,shrink=1,aspect=15)
    sub.set_xlabel(r'x',fontsize=15)
    sub.set_ylabel(r'y',fontsize=15)
    sub.set_aspect('equal')
    fig.suptitle(r'$t={:.3f}$'.format(t))
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig('./figure/By_{:04d}.png'.format(nt))
    plt.close(fig)

    # Bz---
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    cont = sub.pcolormesh(XX,YY,Bz,shading='gouraud',
        cmap=cm.RdBu_r,norm=TwoSlopeNorm(vcenter=0))
    cb = fig.colorbar(cont,shrink=1,aspect=15)
    sub.set_xlabel(r'x',fontsize=15)
    sub.set_ylabel(r'y',fontsize=15)
    sub.set_aspect('equal')
    fig.suptitle(r'$t={:.3f}$'.format(t))
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig('./figure/Bz_{:04d}.png'.format(nt))
    plt.close(fig)


    # Ux---
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    cont = sub.pcolormesh(XX,YY,ux,shading='gouraud',
        cmap=cm.RdBu_r,norm=TwoSlopeNorm(vcenter=0))
    cb = fig.colorbar(cont,shrink=1,aspect=15)
    sub.set_xlabel(r'x',fontsize=15)
    sub.set_ylabel(r'y',fontsize=15)
    sub.set_aspect('equal')
    fig.suptitle(r'$t={:.3f}$'.format(t))
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig('./figure/Ux_{:04d}.png'.format(nt))
    plt.close(fig)

    # Uy---
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    cont = sub.pcolormesh(XX,YY,uy,shading='gouraud',
        cmap=cm.RdBu_r,norm=TwoSlopeNorm(vcenter=0))
    cb = fig.colorbar(cont,shrink=1,aspect=15)
    sub.set_xlabel(r'x',fontsize=15)
    sub.set_ylabel(r'y',fontsize=15)
    sub.set_aspect('equal')
    fig.suptitle(r'$t={:.3f}$'.format(t))
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig('./figure/Uy_{:04d}.png'.format(nt))
    plt.close(fig)

    # Uz---
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    cont = sub.pcolormesh(XX,YY,uz,shading='gouraud',
        cmap=cm.RdBu_r,norm=TwoSlopeNorm(vcenter=0))
    cb = fig.colorbar(cont,shrink=1,aspect=15)
    sub.set_xlabel(r'x',fontsize=15)
    sub.set_ylabel(r'y',fontsize=15)
    sub.set_aspect('equal')
    fig.suptitle(r'$t={:.3f}$'.format(t))
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig('./figure/Uz_{:04d}.png'.format(nt))
    plt.close(fig)




    # rho---
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    cont = sub.pcolormesh(XX,YY,rho,shading='gouraud',cmap=cm.plasma)
    cb = fig.colorbar(cont,shrink=1,aspect=15)
    sub.set_xlabel(r'x',fontsize=15)
    sub.set_ylabel(r'y',fontsize=15)
    sub.set_aspect('equal')
    fig.suptitle(r'$t={:.3f}$'.format(t))
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig('./figure/rho_{:04d}.png'.format(nt))
    plt.close(fig)


    # P---
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    cont = sub.pcolormesh(XX,YY,p,shading='gouraud',cmap=cm.plasma)
    cb = fig.colorbar(cont,shrink=1,aspect=15)
    sub.set_xlabel(r'x',fontsize=15)
    sub.set_ylabel(r'y',fontsize=15)
    sub.set_aspect('equal')
    fig.suptitle(r'$t={:.3f}$'.format(t))
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig('./figure/P_{:04d}.png'.format(nt))
    plt.close(fig)

