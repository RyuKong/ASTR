import numpy as np
from netCDF4 import Dataset
import scipy.ndimage 
import os
import glob
import sys
import tarfile
import gzip
from tqdm import tqdm

def name(**variable): 
    return [x for x in variable][0]

def load_zm(filename, data, varlist):
    t0 = 300. #data.variables['T00'][:] # base temperature
    p0 = 610. #data.variables['P00'][:] # base pressure
    r_d = 191.8366
    cp = 767.3466
    g = 3.727
    gamma = r_d/cp
    
    ls = data.variables['L_S'][:] # solar long
    
    tmp = []
    for var in varlist:
        if var == 'T':
            t = data.variables['T'][:][:,:,16:20,60:64] # perturbation potential temp
            
            p = data.variables['P'][:][:,:,16:20,60:64] # perturbation pressure
            pb = data.variables['PB'][:][:,:,16:20,60:64] # base state pressure
            
            pot_temp = t + t0 # potential temperature
            p = p + pb # pressure
            del t, pb
            
            tmp.append((pot_temp*(p/p0)**gamma)) # temperature
            tmp.append(p)
        elif var == 'PH':
            ph = data.variables['PH'][:]
            phb = data.variables['PHB'][:]
            
            tmp.append((ph+phb).mean(axis = 3))
        elif var == 'U':
            u = data.variables['U'][:]
            tmp.append(u[:,:,16:20,61:65])
        elif var == 'PSFC':
            u = data.variables['PSFC'][:][:,16:20,60:64]
            tmp.append(u)
        else: 
            tmp2 = data.variables[var][:]
            tmp.append(tmp2[:,:,16:20,60:64])
    return ls, tmp

#def load_zm(filename, data, varlist):
#    t0 = 300. #data.variables['T00'][:] # base temperature
#    p0 = 610. #data.variables['P00'][:] # base pressure
#    r_d = 191.8366
#    cp = 767.3466
#    g = 3.727
#    gamma = r_d/cp
#    
#    ls = data.variables['L_S'][:] # solar long
#    
#    tmp = []
#    for var in varlist:
#        if var == 'T':
#            t = data.variables['T'][:]#[:,:,16:] # perturbation potential temp
#            
#            p = data.variables['P'][:]#[:,:,16:] # perturbation pressure
#            pb = data.variables['PB'][:]#[:,:,16:] # base state pressure
#            
#            pot_temp = t + t0 # potential temperature
#            p = p + pb # pressure
#            del t, pb
#            
#            tmp.append((pot_temp*(p/p0)**gamma)) # temperature
#            tmp.append(p.mean(axis=3))
#        elif var == 'PH':
#            ph = data.variables['PH'][:][:,:,16:].mean(axis = 3)
#            phb = data.variables['PHB'][:][:,:,16:].mean(axis = 3)
#            
#            tmp.append((ph+phb)/(1000*g))
#        elif var == 'U':
#            u = data.variables['U'][:]
#            tmp.append(u[:,:,16:20,61:65])
#        elif var == 'PSFC':
#            u = data.variables['PSFC'][:][:,16:20,60:64]
#            tmp.append(u)
#        else: 
#            tmp2 = data.variables[var][:][:,:,16:]
#            tmp.append(tmp2)
#    return ls, tmp

def load_misc(filename, data, var_name):
    var = var_name + '_AM'
    t_am = data.variables[var][:] # solar long
    
    var2 = var_name + '_PM'
    t_pm = data.variables[var2][:] # solar long
    
    ls = data.variables['L_S'][:]

    t_d = t_pm - t_am
    t_a = t_pm + t_am 
    t_d_2pa = t_pm[:,42] - t_am[:,42]
    t_a_2pa = t_pm[:,42] + t_am[:,42]

    return t_a/2., t_d/2.,  t_a_2pa/2., t_d_2pa/2., ls

def load_TAU(filename, data, var_name):
    var = var_name + '_AM'
    t_am = data.variables[var][:] # solar long
    
    var2 = var_name + '_PM'
    t_pm = data.variables[var2][:] # solar long

#    t_d = t_pm - t_am
#    t_a = t_pm + t_am 

    return t_am, t_pm

def load_misc_zm(filename, data, var_name):
    temp = data.variables[var_name][:].mean(axis=3)
    return temp

def load_misc3D(filename, data, var_name):
    temp = data.variables[var_name][:]
    ls = data.variables['L_S'][:]
    return temp, ls

def load_misc4D(filename, data, var_name):
    temp = data.variables[var_name][:]
    return temp
    
def init_reduction(filedir):
    arg = sys.argv[1]
    tar_cond = int(sys.argv[2])
    def wrfout():
        filepath = filedir + '/wrfout_d01*'
        print (filepath)
        
        if not os.path.exists(filedir+'/reduction'): 
            os.mkdir(filedir+'/reduction')
        	#    print filedir 
        varlist = np.array(['T']) #, 'PH', 'TAU_CL', 'TAU_OD'
        for num, i in enumerate(tqdm(sorted(glob.glob(filepath)[:]))):
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
                
                lsd, tmp = load_zm(nc_file, data, varlist)
                print (np.array(tmp[0]).shape)
            else:
                nc_file = i
                data = Dataset(nc_file, mode='r')
        		    
                lsd2, tmp2 = load_zm(nc_file, data, varlist)
                
                lsd = np.concatenate((lsd, lsd2))
                tmp = tmp + tmp2
        
        def create3D_var(varnameList, units, data):   
            tmp2 = dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3],), zlib=True)
            tmp2.units = (units)
            tmp2[:] = data
        def create4D_var(varnameList, units, data):   
            tmp2 = dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3], varnameList[4],), zlib=True)
            tmp2.units = (units)
            tmp2[:] = data
        
        dataset = Dataset('./test2.nc', 'w')
        
        solar_long = dataset.createDimension('time', None)
        pressure = dataset.createDimension('bottom_top', None)
        height = dataset.createDimension('bottom_top_stag', None)
        latitude = dataset.createDimension('south_north', None)
        longitude = dataset.createDimension('west_east', None)
        
        ls = dataset.createVariable('LS', np.float32, ('time',))
        lat = dataset.createVariable('LAT', np.float32, ('south_north',))
        long = dataset.createVariable('LONG', np.float32, ('west_east',))
        
        
        varlen = varlist.size + 1
        
        create4D_var(['T', 'time', 'bottom_top', 'south_north', 'west_east'], 'K', np.vstack(tmp[0::varlen]))
        create3D_var(['P', 'time', 'bottom_top', 'south_north'], 'Pa', np.vstack(tmp[1::varlen]))
        #create3D_var(['PH', 'time', 'bottom_top_stag', 'south_north'], 'km', np.vstack(tmp[2::varlen]))
        #create4D_var(['TAU_CL', 'time', 'bottom_top', 'south_north', 'west_east'], 'unitless', np.vstack(tmp[3::varlen]))
        #create4D_var(['TAU_OD', 'time', 'bottom_top', 'south_north', 'west_east'], 'unitless', np.vstack(tmp[4::varlen]))
#        create4D_var(['U', 'time', 'bottom_top', 'south_north', 'west_east_stag'], 'm s-1', tmp[2])
#        create3D_var(['PSFC', 'time', 'south_north', 'west_east_stag'], 'Pa', tmp[3])
                                                
        #amp = dataset.createVariable('M0_ZM', np.float32, ('martian_months', 'bottom_top', 'south_north',))
        
        ls.units = ('Solar Longtitude')
        lat.units = ('Degree')
        long.units = ('Degree')
        b = np.linspace(-180,180,72)
        a = np.linspace(-90,90,36)
        ls[:] = lsd
        lat[:] = a[16:]
        long[:] = b
#        lat[:] = np.array([-7.71428571, -2.57142857,  2.57142857,  7.71428571])
#        long[:] = np.array([ 124.22535211,  129.29577465,  134.36619718,  139.43661972])
#        long_u[:] = np.array([ 125.,  130.,  135., 140.])
        
        dataset.close()
        
#        filedir2 = i.replace(i,'{}/reduction/wrfout_{}'.format(filedir, 'ls'))
#        np.save(filedir2, ls)
#        del ls
#        
##        varlist.insert(1, 'P')
#        for j, var in enumerate(varlist):
#            filedir2 = i.replace(i,'{}/reduction/wrfout_ {}'.format(filedir, var+'FULL2'))
#            print('Saving', var)
#            np.save(filedir2, tmp[j])
            
    def wrfout_ext():
        filepath = filedir + '/wrfout_d01*'
        print (filepath)
        
        if not os.path.exists(filedir+'/reduction'): 
            os.mkdir(filedir+'/reduction')
        	    
        for num, i in enumerate(sorted(glob.glob(filepath))):
            print(i)
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
                
                t = load_misc4D(nc_file, data, 'T' )[:,16] # t at 2 km
                psfc, ls_psfc = load_misc3D(nc_file, data, 'PSFC' )
                v = load_misc_zm(nc_file, data, 'V' )
            else:
                nc_file = i
                data = Dataset(nc_file, mode='r')
        		    
                t2 = load_misc4D(nc_file, data, 'T')[:,16]
                psfc2, ls_psfc2 = load_misc3D(nc_file, data, 'PSFC')
                v2 = load_misc_zm(nc_file, data, 'V' )
        		    
                t = np.concatenate((t, t2),axis=0)
                psfc = np.concatenate((psfc, psfc2),axis=0)
                v = np.concatenate((v, v2),axis=0)
                ls_psfc = np.concatenate((ls_psfc, ls_psfc2),axis=0)
        	    
        filedir2 = i.replace(i,'{}/reduction/wrfout'.format(filedir))
        varlist = ['_T2KM', '_PSFC', '_ls_PSFC', '_V']
        for num, i in enumerate([t, psfc, ls_psfc, v]):
            print('Saving', varlist[num])
            np.save(filedir2 + varlist[num], i)
    
    def auxhist9():
        filepath = filedir + '/auxhist9*'
        print (filepath)
        
        if not os.path.exists(filedir+'/reduction'): 
            os.mkdir(filedir+'/reduction')
        
        varlist = ['T_PHY']
        for var in varlist:
            #print (sorted(glob.glob(filepath)))
            for num, i in enumerate(tqdm(sorted(glob.glob(filepath)))):
                if num == 0:
                    nc_file = i
                    data = Dataset(nc_file, mode='r')
        		
                    avg, diff, avg2Pa, diff2Pa, ls = load_misc(nc_file, data, var)
                else: 
                    nc_file = i
                    data = Dataset(nc_file, mode='r')
        		
                    avg2, diff2, avg2Pa2, diff2Pa2, ls2 = load_misc(nc_file, data, var)
                    
                    avg = np.concatenate((avg, avg2),axis=0)
                    diff = np.concatenate((diff, diff2),axis=0)
                    
                    avg2Pa = np.concatenate((avg2Pa, avg2Pa2),axis=0)
                    diff2Pa = np.concatenate((diff2Pa, diff2Pa2),axis=0)
                    
                    ls = np.concatenate((ls, ls2), axis=0)
	
            filedir2 = i.replace(i,'{}/reduction/wrfout_'.format(filedir))        
            varlist2 = ['_AVG', '_DIFF', '_AVG2PA', '_DIFF2PA', '_ls_AUX9']
            for num, i in enumerate([avg, diff, avg2Pa, diff2Pa, ls]):
                np.save(filedir2 + var + varlist2[num], i)
                
        varlist = ['TAU_OD', 'TAU_CL']
        for var in varlist:
            print (var)
            for num, i in enumerate(tqdm(sorted(glob.glob(filepath)))):
                if num == 0:
                    nc_file = i
                    data = Dataset(nc_file, mode='r')
        		
                    am, pm = load_TAU(nc_file, data, var)
                else: 
                    nc_file = i
                    data = Dataset(nc_file, mode='r')
        		
                    am2, pm2 = load_TAU(nc_file, data, var)
                    
                    am = np.concatenate((am, am2),axis=0)
                    pm = np.concatenate((pm, pm2),axis=0)
	
            filedir2 = i.replace(i,'{}/reduction/wrfout_'.format(filedir))        
            varlist2 = ['_AM', '_PM']
            for num, i in enumerate([am, pm]):
                np.save(filedir2 + var + varlist2[num], i)
 
    def auxhist5():    
        filepath2 = filedir+'/auxhist5*'
        
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
                print (psfc2.shape, ls_psfc2.shape)
                psfc = np.concatenate((psfc, psfc2),axis=0)
                ls_psfc = np.concatenate((ls_psfc, ls_psfc2),axis=0)

        filedir3 = i.replace('auxhist5','reduction/wrfout')
        var_list = ['_psfc','_ls_psfc']
        for num, i in enumerate([psfc,ls_psfc]):
            print(psfc.shape)
            np.save(filedir3 + var_list[num], i)
            
    def misc():
        filepath2 = filedir+'/wrfout*'
        
        varlist = ['T', 'U', 'PSFC']
        for num, i in enumerate(tqdm(sorted(glob.glob(filepath2)))):
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
                
                lsd, tmp = load_zm(nc_file, data, varlist)
            else:
                nc_file = i
                data = Dataset(nc_file, mode='r')
        		    
                lsd2, tmp2 = load_zm(nc_file, data, varlist)
                
                lsd = np.concatenate((lsd, lsd2))
                tmp = tmp + tmp2
          
        def create3D_var(varnameList, units, data):   
            tmp2 = dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3],))
            tmp2.units = (units)
            tmp2[:] = data
        def create4D_var(varnameList, units, data):   
            tmp2 = dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3], varnameList[4],))
            tmp2.units = (units)
            tmp2[:] = data
        
        dataset = Dataset('./r14p1_galeCrater.nc', 'w')
        
        # create complex128 compound data type.
    #    complex64 = np.dtype([("real",np.float32),("imag",np.float32)])
    #    c64 = dataset.createCompoundType(complex64, "complex64")
        
        solar_long = dataset.createDimension('time', None)
        pressure = dataset.createDimension('bottom_top', None)
        latitude = dataset.createDimension('south_north', None)
        longitude = dataset.createDimension('west_east', None)
        longitude_u = dataset.createDimension('west_east_stag', None)
        
        ls = dataset.createVariable('LS', np.float32, ('time',))
        lat = dataset.createVariable('LAT', np.float32, ('south_north',))
        long = dataset.createVariable('LONG', np.float32, ('west_east',))
        long_u = dataset.createVariable('LONG_U', np.float32, ('west_east_stag',))
            
        create4D_var(['PHY_T', 'time', 'bottom_top', 'south_north', 'west_east'], 'K', np.vstack(tmp[0::4]))
        create4D_var(['P', 'time', 'bottom_top', 'south_north', 'west_east'], 'Pa', np.vstack(tmp[1::4]))
        create4D_var(['U', 'time', 'bottom_top', 'south_north', 'west_east_stag'], 'm s-1', np.vstack(tmp[2::4]))
        create3D_var(['PSFC', 'time', 'south_north', 'west_east'], 'Pa', np.vstack(tmp[3::4]))
                                                
        #amp = dataset.createVariable('M0_ZM', np.float32, ('martian_months', 'bottom_top', 'south_north',))
        
        ls.units = ('Solar Longtitude')
        lat.units = ('Degree')
        long.units = ('Degree')
        long_u.units = ('Degree')
        ls[:] = lsd
        lat[:] = np.array([-7.71428571, -2.57142857,  2.57142857,  7.71428571])
        long[:] = np.array([ 124.22535211,  129.29577465,  134.36619718,  139.43661972])
        long_u[:] = np.array([ 125.,  130.,  135., 140.])
        
        dataset.close()

    if arg == 'wrfout': wrfout()
    if arg == 'wrfout_ext': wrfout_ext()
    if arg == 'auxhist9': auxhist9()
    if arg == 'auxhist5': auxhist5()
    
    if arg == 'misc': misc()
    
    print (tar_cond)
    if tar_cond == 1:
        tar = tarfile.open(filedir+'/reduction.tar.gz', 'w:gz')
        for i in glob.glob(filedir+'/reduction/wrfout*'):
            print ('Tarring file,', i)            
            tar.add(i, arcname = i.replace(filedir, ''))
            tar.close()
init_reduction('./../r14p1dustL45')
 
