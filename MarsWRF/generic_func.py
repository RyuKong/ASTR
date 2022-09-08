import numpy as np
from netCDF4 import Dataset
import scipy.ndimage 
from scipy.interpolate import griddata
import dwell.fft as fft
import matplotlib
matplotlib.use('Agg')
#import matplotlib
import matplotlib.pylab as plt
import matplotlib.colors as colors
import matplotlib.tri as tri
import os
import glob
import sys
import tarfile

matplotlib.rcParams.update({'font.size': 9})

def name(**variable): 
    return [x for x in variable][0]

def load_temp(filename, data):
    t0 = 300. #data.variables['T00'][:] # base temperature
    p0 = 610. #data.variables['P00'][:] # base pressure
    r_d = 191.8366
    cp = 767.3466
    g = 3.727
    gamma = r_d/cp
    
    t = data.variables['T'][:] # perturbation potential temp
    
    p = data.variables['P'][:] # perturbation pressure
    pb = data.variables['PB'][:] # base state pressure
    
    ls = data.variables['L_S'][:] # solar long
    ph = data.variables['PH'][:] # perturbation geopot height
    phb = data.variables['PHB'][:] # base state geopot height
    
    u = data.variables['U'][:] # zonal wind
    
    pot_temp = t + t0 # potential temperature
    press = p + pb # pressure
    geo_height = ph + phb
    geo_height = geo_height.mean(axis = 3)
    
    temp = pot_temp*(press/p0)**gamma # temperature
    
    temp = temp.mean(axis = 3) # avg in long
    press = press.mean(axis = 3) # avg in long
    u = u.mean(axis = 3) # avg in long
    
    return ls, temp, press, geo_height, u

def load_misc(filename, data, var_name):
    var = var_name + '_PHY_AM'
    t_am = data.variables[var][:] # solar long
    
    var2 = var_name + '_PHY_PM'
    t_pm = data.variables[var2][:] # solar long

    t_d = t_pm.mean(axis=3) - t_am.mean(axis=3) 
    t_d_2pa = t_pm[:,42,:,:] - t_am[:,42,:,:]

    return (t_d)/2., t_d_2pa/2.

def load_misc3D(filename, data, var_name):
    temp = data.variables[var_name][:]
    ls = data.variables['L_S'][:]
    return temp, ls
    
def init_reduction(filedir):
    arg = sys.argv[1]
    def wrfout():
        filepath = filedir + '/wrfout_d01*'
        print (filepath)
        
        tar = tarfile.open(filedir+'/reduction.tar.gz', 'w:gz')
        if not os.path.exists(filedir+'/reduction'): 
            os.mkdir(filedir+'/reduction')
        	#    print filedir
        	    
        for num, i in enumerate(sorted(glob.glob(filepath))):
            print(i)
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
                
                ls, temp, press, geoH, u = load_temp(nc_file, data)
            else:
                nc_file = i
                data = Dataset(nc_file, mode='r')
        		    
                ls2, temp2, press2, geoH2, u2 = load_temp(nc_file, data)
        		    
                ls = np.concatenate((ls, ls2),axis=0)
                temp = np.concatenate((temp, temp2),axis=0)
                press = np.concatenate((press, press2),axis=0)
                geoH = np.concatenate((geoH, geoH2),axis=0)
                u = np.concatenate((u, u2),axis=0)
        	    
        filedir2 = i.replace('wrfout','reduction/wrfout')
        var_list = ['_ls','_temp','_press','_geoH', '_u']
        for num, i in enumerate([ls,temp,press,geoH,u]):
            print('Saving', var_list[num])
            np.save(filedir2 + var_list[num], i)
            filename = filedir2.replace(filedir, '')
            tar.add(filedir2 + var_list[num]+'.npy', arcname = filename + var_list[num]+ '.npy')
        tar.close()
    
    def auxhist9():
        print ('hi')  
        print ('cool')
        filepath2 = filedir + '/auxhist9*'
        print (filepath2)
        
        tar = tarfile.open(filedir+'/reduction.tar.gz', 'w:gz')
        for num, i in enumerate(sorted(glob.glob(filepath2))):
            print (i)
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
		
                t_d, t_d_2Pa = load_misc(nc_file, data, 'T' )
            else: 
                nc_file = i
                data = Dataset(nc_file, mode='r')
		
                t_d2, t_d2_2Pa = load_misc(nc_file, data, 'T')
                t_d = np.concatenate((t_d, t_d2),axis=0)
                t_d_2Pa = np.concatenate((t_d_2Pa, t_d2_2Pa),axis=0)
	
        filedir3 = i.replace('auxhist9','reduction/wrfout')        
        var_list = ['_t_d','_t_d_2Pa']
        for num, i in enumerate([t_d,t_d_2Pa]):
            np.save(filedir3 + var_list[num], i)
            filename = filedir3.replace(filedir, '')
            tar.add(filedir3 + var_list[num]+'.npy', arcname = filename + var_list[num]+ '.npy')
        tar.close()
 
    def auxhist5():    
        filepath2 = filedir+'/auxhist5*'
        
        tar = tarfile.open(filedir+'/reduction.tar.gz', 'w:gz')
        for num, i in enumerate(sorted(glob.glob(filepath2))):
            print (i)
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')

                psfc, ls_psfc = load_misc3D(nc_file, data, 'PSFC' )
            else:
                nc_file = i
                data = Dataset(nc_file, mode='r')

                psfc2, ls_psfc2 = load_misc3D(nc_file, data, 'PSFC')
                psfc = np.concatenate((psfc, psfc2),axis=0)
                ls_psfc = np.concatenate((ls_psfc, ls_psfc2),axis=0)

        filedir3 = i.replace('auxhist5','reduction/wrfout')
        var_list = ['_psfc','_ls_psfc']
        for num, i in enumerate([psfc,ls_psfc]):
            print(psfc.shape)
            np.save(filedir3 + var_list[num], i)
            filename = filedir3.replace(filedir, '')
            tar.add(filedir3 + var_list[num]+'.npy', arcname = filename + var_list[num]+ '.npy')
        tar.close()

    if arg == 'wrfout': return wrfout()
    if arg == 'auxhist9': return auxhist9()
    if arg == 'auxhist5': return auxhist5()
#init_reduction('./../planetWRF/WRFV3/run/new_wbm')
 
def find_ls_idx(ls_arr, ls):
    idx = (np.abs(ls_arr-ls)).argmin() # finding index corresponding to wanted solar long
    assert np.min(np.abs(ls_arr-ls)) < 1, 'Solar Longitude {} not in the file'.format(ls)
    return idx


def zonal_avg(filedir, month, ls2):
    print ('Looking at zonal avg')
    ### ls1, ls2 - solar longitude 1, 2 ###
    
    filepath = glob.glob(filedir + '*_temp.npy')[0]
    print (filepath)
    temp = np.load(filepath)
    
#    filepath = filedir + '*_u.npy'
#    print filepath
#    u = np.load(filepath)
    
    filepath = glob.glob(filedir + '*_press.npy')[0]
    print (filepath)
    p = np.load(filepath)
        
    filepath = glob.glob(filedir + '*_ls.npy')[0]
    ls = np.load(filepath)[::2]

    idx1 = np.where(ls == 360)[0][1] # only looking at the second year
    idx2 = np.where(ls == 360)[0][2]
    
    ls = ls[idx1:idx2]
    temp = temp[idx1:idx2]
#    u = u[idx:]
    p = p[idx1:idx2]
    print (ls)
    
    zonal_t = []
#    zonal_u = []
    zonal_p = []
    for i in np.arange(0, 12):
        print ('month', i+1)
        idx = np.where((ls>i*30)&(ls<(i+1)*30))
        zonal_t.append(temp[idx].mean(axis=0))
#        zonal_u.append(u[idx].mean(axis=0))
        zonal_p.append(p[idx].mean(axis=0))
    del temp, p
    zonal_t = np.array(zonal_t)
#    zonal_u = np.array(zonal_u)
    zonal_p = np.array(zonal_p)
    
    print ('Plotting some cool shit')
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
        press = zonal_p[i-1][4:]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        temp_press = np.linspace(1e-2, 900, 52)[4:]
        
        lat, temp_press = np.meshgrid(lat, temp_press)
        
        temp = zonal_t[i-1][4:]
        print (temp_press.shape, lat.shape, temp.shape)
#        temp = temp[4:]
        
        im = ax.contourf(lat, press, temp, 12, cmap='viridis')
        ax.contour(lat, press, temp, 12, colors='k')
        
        ax.set_title(r'Avg Zonal T$_d$ LS {}-{}'.format((i)*30, (i+1)*30))
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_ylim([900, 1e-2])
    
    print ('Saving 1st cool shit')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
    plt.savefig('zonal_temp.png',bbox_inches='tight', dpi=600)
    
    print ('Saving 2nd cool shit')
    
#    cb = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
#    plt.savefig('test2.png',bbox_inches='tight', dpi=600)

    
#    print 'Plotting zonal wind' 
#    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
#    for i, ax in enumerate(axes.flat):
#        u = zonal_u[i-1]
#        im = ax.contourf(lat, press, u, cmap='inferno')
#        ax.contour(lat, press, u, linewidths=0.5, colors='black')
#        ax.set_title('Avg Zonal Wind LS {}-{}'.format((i)*30, (i+1)*30))
#        ax.invert_yaxis()
#        ax.set_yscale('log')
        
#    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
#    plt.savefig('test2.png',bbox_inches='tight', dpi=600)

def zonal_diff(filedir, month, ls2):
    print ('Looking at thermal tides')
    
    filepath = glob.glob(filedir + '*_press.npy')[0]
    print (filepath)
    p = np.load(filepath)
        
    filepath = glob.glob(filedir + '*_t_d.npy')[0]
    t_d = np.load(filepath)
    t_d = t_d
    
    filepath = glob.glob(filedir + '*_t_d_2Pa.npy')[0]
    t_d_2Pa = np.load(filepath)
    t_d_2Pa = t_d_2Pa 
    
    filepath = glob.glob(filedir + '*_ls.npy')[0]
    ls = np.load(filepath)[::2]

    idx1 = np.where(ls == 360)[0][1] # only looking at the second year
    idx2 = np.where(ls == 360)[0][2]
    
    ls = ls[idx1:idx2]
    t_d = t_d[idx1:idx2]
    t_d_2Pa = t_d_2Pa[idx1:idx2]
    p = p[idx1:idx2]
    
    print (t_d_2Pa.shape)
    
#    test_diff = t_d_2Pa.reshape((223,3,36,72)).mean(axis=1)
    ampl, phase, axis = fft.spec1d(t_d_2Pa, 1/72., use_axes = 2)
    ampl = ampl
    phase = phase
    print (np.min(ampl))
#    ampl.hstack(axis=0)
#    phase.hstack(axis=0)
    
    test = np.zeros((669,36))
    test2 = np.zeros((669,36))
    for i in np.arange(0,669):
        test[i] = ampl[i,:,0]
        test2[i] = phase[i,:,0]
        
    lat = np.linspace(-90, 90, 36) 
    ls = np.linspace(0,360, 669)
    t_lat, t_ls = np.meshgrid(lat, ls)
    plt.figure()
    plt.contourf(t_ls, t_lat, test, cmap='viridis')
#    plt.contour(t_ls, t_lat, np.tan(test2), colors='k')
    plt.colorbar()
    
    zonal_p = []
    zonal_t_d = []
    zonal_t_2P = []
    for i in np.arange(0, 12):
        print ('month', i+1)
        idx = np.where((ls>i*30)&(ls<(i+1)*30))
        zonal_p.append(p[idx].mean(axis=0))
        zonal_t_d.append(t_d[idx].mean(axis=0))
        zonal_t_2P.append(t_d_2Pa[idx].mean(axis=0))
    zonal_t_d = np.array(zonal_t_d)
    zonal_t_2P = np.array(zonal_t_2P)
    zonal_p = np.array(zonal_p)
    
    print ('Plotting some cool shit')
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
        press = zonal_p[i-1][4:]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        temp_press = np.linspace(1e-2, 900, 52)[4:]
        
        #lat, press = np.meshgrid(lat, press2[5:])
        lat, temp_press = np.meshgrid(lat, temp_press)
        
        temp = zonal_t_d[i-1][4:]
        
        im = ax.contourf(lat, press, temp, 12, cmap='viridis')
        ax.contour(lat, press, temp, 12, linewidths=0.5, colors='black')
        
        ax.set_title(r'Avg Zonal T$_d$ LS {0}-{1}'.format((i)*30, (i+1)*30))
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_ylim([900, 1e-2])
        
    print ('Saving 1st cool shit')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
#    plt.savefig('zonal_tdiff.png',bbox_inches='tight', dpi=600)
    
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
        t_2P = zonal_t_2P[i-1]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        lon = np.linspace(0, 360, 72)
        
        lon, lat = np.meshgrid(lon, lat)

        im = ax.contourf(lon, lat, t_2P, cmap='viridis')
        ax.contour(lon, lat, t_2P, linewidths=0.5, colors='black')
        ax.set_title(r'Avg T$_d$ LS {0}-{1} at 2 Pa'.format((i)*30, (i+1)*30))
        
    print ('Saving 2nd cool shit')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
#    plt.savefig('zonal_tdiff_2Pa.png',bbox_inches='tight', dpi=600)

#zonal_avg('./test_data/reduction/',12,2)
#zonal_diff('./test_data/reduction/',12,2)

def fft_hovmoller(filedir):
    print ('Looking at spectral of surface presssure')
    
    filepath = glob.glob(filedir + '*0_psfc.npy')[0]
    print (filepath)
    psfc = np.load(filepath)
    
    filepath = glob.glob(filedir + '*_ls_psfc.npy')[0]
    print (filepath)
    ls = np.load(filepath)
    
    idx1 = np.where(ls == 360)[0][1] # only looking at the second year
    idx2 = np.where(ls == 360)[0][2]
    ls = ls[idx1:idx2]
    psfc = psfc[idx1:idx2]
    
    #avging for specific latitude
    lat = np.linspace(-90,90,36)
    idx = np.where((lat>-10)&(lat<10))[0]
    
    psfc = psfc[:,18,:]#.mean(axis=1)
    
    idx = np.where((ls>300)&(ls<320))[0]
    ls = ls[idx]
    psfc = psfc[idx] 
    psfc = psfc - psfc.mean(axis=0)
#    
    lon = np.linspace(0,360,72)
    lon, ls = np.meshgrid(lon, ls)
    
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    plt.contourf(lon, ls, psfc)
    plt.ylabel('Solar Longitude')
    plt.xlabel('Longitude')
    plt.title(r'Hovm$\mathrm{\"{o}}$ller Diagram of Surface Pressure')

    lat = np.linspace(-90,90,36)
    
    plt.subplot(1,2,2)
    ampl, phase, cycle, waven = fft.spec(psfc, 1./8, 1./72, axes=[0,1])
    plt.contourf(cycle, waven[:5], np.log10(ampl).T[:5])
    plt.ylabel('Wavenumber')
    plt.xlabel('Cycles/sol')
    plt.title(r'Amplitude of FFT of Hovm$\mathrm{\"{o}}$ller Diagram')
    
fft_hovmoller('./test_data/reduction/')

def zonal_temperature2(filename, ls1, ls2):
    nc_file = filename
    data = Dataset(nc_file, mode='r')
    
    t = data.variables['T'][:] # perturbation potential temp
    p = data.variables['P'][:] # perturbation pressure
    pb = data.variables['PB'][:] # base state pressure
    ls = data.variables['L_S'][:] # solar longitude
    
    idx1 = find_ls_idx(ls, ls1)
    idx2 = find_ls_idx(ls, ls2)

    t0 = 300. #data.variables['T00'][:] # base temperature
    p0 = 610. #data.variables['P00'][:] # base pressure
    r_d = 191.8366
    cp = 767.3466
    g = 3.727
    gamma = r_d/cp
    
    pot_temp = t + t0 # potential temperature
    press = p + pb # pressure
    
    temp = pot_temp*(press/p0)**gamma # temperature
    
    temp_avg = temp.mean(axis = 3)[idx1:idx2].mean(axis = 0) # avg in long and solar long
    press_avg = press.mean(axis = 3).mean(axis = 2)[idx1:idx2].mean(axis = 0) # avg in lat, long and solar long
    
    temp_total_avg = temp_avg.mean(axis=1).mean() # total avg temp
    h = (r_d*(temp_total_avg/g))*(np.log(p0/press_avg)/np.log(np.e)) # hypsometric equation?
    
    return temp_avg, press_avg, h/1000.

def zonal_wind(filename, ls1, ls2):
    nc_file = filename
    data = Dataset(nc_file, mode='r')
    
    u = data.variables['U'][:] # zonal wind speed in m/s
#    p = data.variables['P'][:] # perturbation pressure
#    pb = data.variables['PB'][:] # base state pressure
    ls = data.variables['L_S'][:] # solar longitude
    
    idx1 = find_ls_idx(ls, ls1)
    idx2 = find_ls_idx(ls, ls2)
    
    u_avg = u.mean(axis = 3)[idx1:idx2].mean(axis = 0) # avg in long and solar long
    
    return u_avg

def net_hr_aer(filename, ls1, ls2):
    nc_file = filename
    data = Dataset(nc_file, mode='r')
    
    hrvis = data.variables['HRAERVIS'][:] # heating rate in visible
    hrir = data.variables['HRAERIR'][:] # heating rate in infrared
    p = data.variables['P'][:] # perturbation pressure
    pb = data.variables['PB'][:] # base state pressure
    ls = data.variables['L_S'][:] # solar longitude
    
    idx1 = find_ls_idx(ls, ls1)
    idx2 = find_ls_idx(ls, ls2)
    
    press = p + pb # pressure

    hrvis_avg = hrvis.mean(axis = 3)[idx1:idx2].mean(axis = 0) # avg in long and solar long 
    hrir_avg = hrir.mean(axis = 3)[idx1:idx2].mean(axis = 0) # avg in long and solar long
    net_ir = hrvis_avg + hrir_avg # net heating rate due to aerosols
    
    press_avg = press.mean(axis = 3).mean(axis = 2)[idx1:idx2].mean(axis = 0) # avg in lat, long and solar long
    
    return net_ir, press_avg
    
