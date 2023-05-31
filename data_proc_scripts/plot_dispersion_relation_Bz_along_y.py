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

XX, YY = np.meshgrid(xgrid,ygrid,indexing='ij')

files = sorted(glob.glob('./output/out*.dat'))
nout = len(files)


time = np.zeros(nout)
Bz_cut = np.zeros([nout,ny])


for nt in range(nout):
    t, uu = read_output(files[nt],8,nx,ny)

    print('t = {:.3f}'.format(t))

    time[nt] = t 

    Bz_cut[nt,:] = uu[6,int(nx/2),:]


arr = np.copy(Bz_cut)

dt = time[1] - time[0]
dy = ygrid[1] - ygrid[0]

# fourier transform
nk = ny
no = nt 

arr_tk = np.zeros([nt,nk],dtype='complex')

for i in range(nt):
    temp = arr[i,:]
    temp_fft = np.fft.fftshift(np.fft.fft(temp))
    arr_tk[i,:] = temp_fft

arr_ok = np.zeros([no,nk],dtype='complex')
for i in range(nk):
    temp = arr_tk[:,i]
    temp_fft = np.fft.fftshift(np.fft.fft(temp))
    arr_ok[:,i] = temp_fft

arr_ko = np.transpose(arr_ok)

k_arr = np.fft.fftshift(np.fft.fftfreq(ny,dy))
o_arr = np.fft.fftshift(np.fft.fftfreq(nt,dt))


# crop the omega - k array
kmax=40
omax=0

if kmax==0:
    kmax = max(k_arr)
if omax==0:
    omax = max(o_arr)

ind_kmax = ((np.where(np.abs(k_arr+kmax) <  (k_arr[1]-k_arr[0]) ))[0])[0]
ind_k0 = ((np.where(np.abs(k_arr) <  1e-7 ))[0])[0]

ind_omax = ((np.where(np.abs(o_arr-omax) <  (o_arr[1]-o_arr[0]) ))[0])[0]
ind_o0 = ((np.where(np.abs(o_arr) <  1e-7 ))[0])[0]


k_arr_new = np.copy(np.flip(np.abs(k_arr[ind_kmax:ind_k0+1]),0))


o_arr_new = np.copy(o_arr[ind_o0:ind_omax+1])
arr_ko_new = np.copy(arr_ko[ind_kmax:ind_k0+1, ind_o0: ind_omax+1])

for i in range(len(o_arr_new)):
    arr_ko_new[:, i] = np.flip(arr_ko_new[:, i],0)

k_arr = np.copy(k_arr_new)
o_arr = np.copy(o_arr_new)
arr_ko = np.copy(arr_ko_new)

# theoretical predictions
def disp_rel_theory(k,di,Va):
    k_ = k * 2 * np.pi  #be careful about the k and omega

    ome_p = 0.5 * (np.sqrt( 4+ (k_*di)**2)+ k_ * di) * k_ * Va
    ome_m = 0.5 * (np.sqrt( 4+ (k_*di)**2)- k_ * di) * k_ * Va

    return ome_p/(2*np.pi), ome_m/(2*np.pi)



k_mesh, o_mesh = np.meshgrid(k_arr,o_arr,indexing='ij')
abs_arr_ko = np.abs(arr_ko)
ome_p, ome_m = disp_rel_theory(k_arr, 0.05, 1.0)


# levels = np.linspace(0,80,32)
fig = plt.figure()
sub = fig.add_subplot(111)
cont = sub.contourf(k_mesh, o_mesh, abs_arr_ko,np.linspace(0,0.3,32),cmap=cm.plasma)
sub.plot(k_arr, ome_p,'--',color='w')
sub.plot(k_arr, ome_m,'--',color='w')


cb = fig.colorbar(cont,aspect=15,shrink=1.0)
sub.set_ylim([np.min(o_arr),np.max(o_arr)])
sub.set_xlabel(r'$k$',fontsize=15)
sub.set_ylabel(r'$f$',fontsize=15)
fig.savefig('./dispersion_relation_Bz_along_y.png',dpi=300)
plt.close(fig)
