import numpy as np 
import struct
import matplotlib.pyplot as plt 
from matplotlib import cm
import sys
import glob
from matplotlib.gridspec import GridSpec
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



xgrid, ygrid = read_grid()
nx = len(xgrid)
ny = len(ygrid)

print('nx = ', nx, ', ny = ', ny)

ky = np.fft.rfftfreq(ny,ygrid[1]-ygrid[0]) * 2*np.pi

files = sorted(glob.glob('./output/out*.dat'))

nout = len(files)

Bx_four_time = np.zeros([len(ky),nout])

time = np.zeros(nout)
for nt in range(nout):
    t, uu = read_output(files[nt],8,nx,ny)
    print('Read file: t = {:.3f}'.format(t))

    Bx_cut_along_y = uu[4,int(nx/2),:]

    Bx_fft = np.abs(np.fft.rfft(Bx_cut_along_y,norm='ortho'))

    Bx_four_time[:,nt] = Bx_fft[:]

    time[nt] = t


nk_plot = 40

growth_rate = np.zeros(nk_plot)

for i in range(nk_plot):
    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.plot(time,Bx_four_time[i+1,:])
    sub.set_ylabel(r'$|\hat{B}_x|$',fontsize=15)
    sub.set_xlabel(r'$t$',fontsize=15)
    sub.set_yscale('log',basey=np.e)
    fig.suptitle(r'$k={:.3f}$'.format(ky[i+1]))
    fig.savefig('./figure/Bx_fourier_time_mode_{:03d}.png'.format(i+1))
    plt.close(fig)

    fit_x = time[40:60]
    fit_y = np.log(Bx_four_time[i+1,40:60])

    fit_ = np.polyfit(fit_x, fit_y, 1)

    growth_rate[i] = fit_[0]


data_ = np.load('./growth_rate_S_1E3.npy')

fig = plt.figure()
sub = fig.add_subplot(111)
sub.plot(ky[1:nk_plot+1], growth_rate)
sub.plot(data_[0,:], data_[1,:],'o')
sub.set_xlabel(r'$k$',fontsize=15)
sub.set_ylabel(r'$\gamma$',fontsize=15)
fig.savefig('./growth_rate_k.png')
plt.close(fig)
