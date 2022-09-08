import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import matplotlib.cm as cmaps
from generic_func import *

matplotlib.rcParams.update({'font.size': 9})

ls1, ls2 = 355,5 # range of the two solar longitudes we want to look at
#temp_avg, press_avg, h = zonal_temperature('/scratch2/p/pen/sokvisal/planetWRF/WRFV3/run/wrfout_d01_0001-00601_00:00:00', ls1, ls2)
#u_avg = zonal_wind('/scratch2/p/pen/sokvisal/planetWRF/WRFV3/run/wrfout_d01_0001-00601_00:00:00', ls1, ls2)
temp_avg, press_avg, h = zonal_temperature('/scratch2/p/pen/sokvisal/mars/planetWRF-dev-new/WRFV3/test/em_global_mars/wrfout_d01_0001-00651_00:00:00', ls1, ls2)
u_avg = zonal_wind('/scratch2/p/pen/sokvisal/mars/planetWRF-dev-new/WRFV3/test/em_global_mars/wrfout_d01_0001-00651_00:00:00', ls1, ls2)


lat = np.linspace(-90,90,36) # latitude
xlong = np.linspace(0,360,72)
lat, press = np.meshgrid(lat, press_avg)

fig = plt.subplots(figsize=(12, 8))
plt.subplot(2,2,1)
ax1 = plt.contourf(lat, press, temp_avg, 11, cmap = cmaps.RdBu_r, origin='upper')
#plt.clabel(ax1, color='k')
ax2 = plt.contour(lat, press, temp_avg, 11, colors='black')
plt.clabel(ax2, color='k', inline = 1, fmt='%1i')
plt.xlabel('Latitude (degree)')
plt.ylabel('Pressure (Pa)')
plt.title(r'Avg Zonal Temperature LS {}-{}'.format(ls1, ls2))
plt.gca().invert_yaxis()
plt.yscale('log')

plt.subplot(2,2,2)
ax1 = plt.contourf(lat, press, u_avg, 11, cmap = cmaps.RdBu_r, origin='upper')
ax2 = plt.contour(lat, press, u_avg, 11, colors='black')
plt.clabel(ax2, color='k', inline = 1, fmt='%1i')
#plt.colorbar(ax1)
plt.xlabel('Latitude (degree)')
plt.ylabel('Pressure (Pa)')
plt.title(r'Avg Zonal Wind LS {}-{}'.format(ls1, ls2))
plt.gca().invert_yaxis()
plt.yscale('log')

ls1, ls2 = 85,95 # range of the two solar longitudes we want to look at
#temp_avg, press_avg, h = zonal_temperature('/scratch2/p/pen/sokvisal/planetWRF/WRFV3/run/wrfout_d01_0002-00132_00:00:00', ls1, ls2)
#u_avg = zonal_wind('/scratch2/p/pen/sokvisal/planetWRF/WRFV3/run/wrfout_d01_0002-00132_00:00:00', ls1, ls2)
temp_avg, press_avg, h = zonal_temperature('/scratch2/p/pen/sokvisal/mars/planetWRF-dev-new/WRFV3/test/em_global_mars/wrfout_d01_0002-00182_00:00:00', ls1, ls2)
u_avg = zonal_wind('/scratch2/p/pen/sokvisal/mars/planetWRF-dev-new/WRFV3/test/em_global_mars/wrfout_d01_0002-00182_00:00:00', ls1, ls2)

plt.subplot(2,2,3)
ax1 = plt.contourf(lat, press, temp_avg, 11, cmap = cmaps.RdBu_r, origin='upper')
#plt.clabel(ax1, color='k')
ax2 = plt.contour(lat, press, temp_avg, 11, colors='black')
plt.clabel(ax2, color='k', inline = 1, fmt='%1i')
plt.xlabel('Latitude (degree)')
plt.ylabel('Pressure (Pa)')
plt.title(r'Avg Zonal Temperature LS {}-{}'.format(ls1, ls2))
plt.gca().invert_yaxis()
plt.yscale('log')

plt.subplot(2,2,4)
ax1 = plt.contourf(lat, press, u_avg, 11, cmap = cmaps.RdBu_r, origin='upper')
ax2 = plt.contour(lat, press, u_avg, 11, colors='black')
plt.clabel(ax2, color='k', inline = 1, fmt='%1i')
#plt.colorbar(ax1)
plt.xlabel('Latitude (degree)')
plt.ylabel('Pressure (Pa)')
plt.title(r'Avg Zonal Wind LS {}-{}'.format(ls1, ls2))
plt.gca().invert_yaxis()
plt.yscale('log')
plt.savefig('old_wbm_ls_0_90.png')
plt.show()
