# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 09:40:05 2017

@author: visal
"""

from netCDF4 import Dataset
import numpy as np

from netCDF4 import Dataset

def create3D_var(varnameList, units, data):   
    tmp = dataset.createVariable(varnameList[0], c64, (varnameList[1], varnameList[2], varnameList[3],))
    tmp.units = ('K')
    tmp[:] = data
#    tmp['imag'] = data.imag
    
    
dataset = Dataset('./test.nc', 'w')

# create complex128 compound data type.
complex64 = np.dtype([("real",np.float32),("imag",np.float32)])
c64 = dataset.createCompoundType(complex64, "complex64")
 

month = dataset.createDimension('martian_months', 12)
pressure = dataset.createDimension('bottom_top', 52)
latitude = dataset.createDimension('south_north', 36)

months = dataset.createVariable('MONTH', np.int32, ('martian_months',))
latitudes = dataset.createVariable('LAT', np.float32, ('south_north',))

pressures = dataset.createVariable('P', np.float32, ('martian_months','bottom_top', 'south_north',))
create3D_var(['M0_ZM', 'martian_months', 'bottom_top', 'south_north'], 'K', np.zeros((12,52,36), complex64))
create3D_var(['M0_DIUR', 'martian_months', 'bottom_top', 'south_north'], 'K', np.zeros((12,52,36)))
create3D_var(['M0_SEMI', 'martian_months', 'bottom_top', 'south_north'], 'K', np.zeros((12,52,36)))
create3D_var(['M0_TER', 'martian_months', 'bottom_top', 'south_north'], 'K', np.zeros((12,52,36)))
#amp = dataset.createVariable('M0_ZM', np.float32, ('martian_months', 'bottom_top', 'south_north',))

months.units = ('440_solar_days_avg')
pressures.units = ('Pa')
latitudes.units = ('Degree')
#amp.units = ('K')

latitudes[:] = np.linspace(-90,90,36)