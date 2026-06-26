#ice microphysical cross sections
#import needed libraries
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pyart
from scipy.special import gamma
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm
import metpy.calc as mpcalc
from metpy.interpolate import cross_section
from matplotlib import colors as mcolors
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm#ice microphysical cross sections
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable

#specify date, hour, and half-hour
date="2023072312"
hour="011"
half_hour="00"
half_hour_selector=1 #0 or 1, time1 index

# Load data
hrdps25data = xr.open_dataset(r'/home/clintonmacadam/data/'+date+'_'+hour+'', engine="fstd")

# Choose nearest lat-lon line for cross-section (e.g., along constant latitude)
target_lat=-0.89
start_lon=619
end_lon=658

#start = (585, 0.2)
#end = (635, 0)
#cross = cross_section(hrdps25datatest, start, end).set_coords(('rlat1', 'rlon1'))
#print(cross)

# Set up colormap and white mask
cmap=plt.get_cmap('magma')
new_colors = cmap(np.linspace(0, 1, 256))
white = np.array([1, 1, 1, 1])  # RGBA for white

# Extract data along the specfied lat index across longitude and levels for variables of interest
reflectivity = hrdps25data['ZET']
reflectivityimg = hrdps25data['ZEC']
reflectivityimg = reflectivityimg[1,:,:]

reflectivity = reflectivity[half_hour_selector,:,:,start_lon:end_lon]
reflectivity = reflectivity.sel(level1=slice(0.15,0.9975))
reflectivity = reflectivity.sel(rlat1=target_lat, method='nearest')

icenumbermr = hrdps25data['NTI1']
icenumbermr = icenumbermr[0,:, :, start_lon:end_lon]
icenumbermr= icenumbermr.sel(level1=slice(0.15,0.9975))
icenumbermr = icenumbermr.sel(rlat1=target_lat, method='nearest')

specifichumidity=hrdps25data['HU']
specifichumidity=specifichumidity[half_hour_selector,:,:,start_lon:end_lon]
specifichumidity = specifichumidity.sel(level1=slice(0.15,0.9975))
specifichumidity = specifichumidity.sel(rlat1=target_lat, method='nearest')
mixingratio=specifichumidity*(1-0.622)

icemassmr = hrdps25data['QTI1']
icemassmr = icemassmr[half_hour_selector,:, :, :]

rimedicemassmr = hrdps25data['QMI1']
print(rimedicemassmr)
rimedicemassmr = rimedicemassmr[0,:,:,:]
rimedicefraction = rimedicemassmr/icemassmr
print(rimedicefraction)
rimedicefraction = rimedicefraction[:, :, start_lon:end_lon]
rimedicefraction = rimedicefraction.sel(level1=slice(0.15,0.9975))
rimedicefraction = rimedicefraction.sel(rlat1=target_lat, method='nearest')
print(rimedicefraction)

icemassmr = icemassmr[:,:,start_lon:end_lon]
icemassmr = icemassmr.sel(level1=slice(0.15,0.9975))
icemassmr = icemassmr.sel(rlat1=target_lat, method='nearest')

rainmassmr = hrdps25data['MPQR']
rainmassmr = rainmassmr[half_hour_selector,:, :, start_lon:end_lon]
rainmassmr = rainmassmr.sel(level1=slice(0.15,0.9975))
rainmassmr = rainmassmr.sel(rlat1=target_lat, method='nearest')

cloudmassmr = hrdps25data['MPQC']
cloudmassmr = cloudmassmr[half_hour_selector,:,:,start_lon:end_lon]
cloudmassmr = cloudmassmr.sel(level1=slice(0.15,0.9975))
cloudmassmr = cloudmassmr.sel(rlat1=target_lat, method='nearest')

liquidmassmr = rainmassmr + cloudmassmr

temperature = hrdps25data['TT']
temperature = temperature[half_hour_selector,:,:,start_lon:end_lon]
temperature = temperature.sel(level1=slice(0.15,0.9975))
temperature = temperature.sel(rlat1=target_lat, method='nearest')
tempcalc=temperature+273.15 #Convert temperature to Kelvin

dewpoint = hrdps25data['TD']
dewpoint = dewpoint[0,:,:,start_lon:end_lon]
dewpoint = dewpoint.sel(level1=slice(0.15,0.9975))
dewpoint = dewpoint.sel(rlat1=target_lat, method='nearest')
dewpointcalc = dewpoint+273.15 #Convert dewpoint to Kelvin

rimedicemassmr = hrdps25data['QMI1']
rimedicemassmr = rimedicemassmr[0,:, :, start_lon:end_lon]
rimedicemassmr = rimedicemassmr.sel(level1=slice(0.15,0.9975))
rimedicemassmr = rimedicemassmr.sel(rlat1=target_lat, method='nearest')

rimevolumemr = hrdps25data['BMI1']
rimevolumemr = rimevolumemr[0,:, :, start_lon:end_lon]
rimevolumemr = rimevolumemr.sel(level1=slice(0.15,0.9975))
rimevolumemr = rimevolumemr.sel(rlat1=target_lat, method='nearest')

rimedensity=rimedicemassmr/rimevolumemr

omega = hrdps25data['WW']
omega = omega[half_hour_selector,:,:, start_lon:end_lon]
omega = omega.sel(level5=slice(0.15,0.9975))
omega = omega.sel(rlat1=target_lat, method='nearest')

es=0.6113*np.exp(((2.5*10**6)/(461))*((1/273.15)-(1/tempcalc)))
e=0.6113*np.exp(((2.5*10**6)/(461))*((1/273.15)-(1/dewpointcalc)))
rh=e/es*100
wetbulbtemp=temperature*np.arctan(0.151977*(rh+8.313659)**(1/2))-4.686035+np.arctan(temperature+rh)-np.arctan(rh-1.676331)+(0.00391838*(rh)**(3/2)*np.arctan(0.023101*rh))

#Determine the level where the wet bulb temperature is closest to 0C/the approximate freezing level
abs_diff = np.abs(wetbulbtemp-0)
min_idx = abs_diff.argmin(dim='level1')
min_idx = min_idx.compute().to_numpy()
print(min_idx)

geopotentialheight = hrdps25data['GZ']
geopotentialheight = geopotentialheight[half_hour_selector,:,:,start_lon:end_lon]
geopotentialheight = geopotentialheight*10
ro=6356766 #radius of the earth, m
a=29.3 #m/K
rd=0.287053 #kPa K-1 m3 kg-1

#Convert geopotential height to geometric height
geometricheight = (ro*geopotentialheight)/(ro-geopotentialheight)
geometricheight = xr.concat([geometricheight.isel(level2=slice(1, 125, 2)), geometricheight.isel(level2=124)], dim="level2")
geometricheight = geometricheight.sel(level2=slice(0.15,0.9975))
geometricheight = geometricheight.sel(rlat1=target_lat, method='nearest')
elevation=geometricheight.sel(level2=1, method='nearest')
geometricheight = geometricheight-elevation #Convert geometric height in ASL to AGL

print(geometricheight[40,30].values)

freezinglevelheight = geometricheight[min_idx,:] #Calculate freezing level height (approximate)
sfcpressure = hrdps25data['P0']
sfcpressure = sfcpressure[half_hour_selector,:,start_lon:end_lon]
sfcpressure = sfcpressure.sel(rlat1=target_lat, method='nearest')

pressure=xr.zeros_like(temperature)
pressure[41,:]=sfcpressure

virtualtemp=tempcalc*(1+0.61*mixingratio-liquidmassmr-icemassmr) #Calculate virtual temperature

i=41 
while i>0: #Calculate pressure on model levels using the hypsometric equation
        dz = geometricheight[i-1, :] - geometricheight[i, :]
        pressure[i-1,:]= pressure[i,:]*np.exp(-dz/(a* (0.5 * (virtualtemp[i,:] + virtualtemp[i-1,:]))))
        i=i-1

pressure=pressure/10 #Convert pressure to kPa

airdensity=pressure/(rd*virtualtemp) #Calculate air density
omega = omega.rename({'level5': 'level1'})

g=9.81 #acceleration due to gravity
verticalvelocity=-omega/(airdensity*g)
print(verticalvelocity)

lat=reflectivity.lat_1.values
lon=reflectivity.lon_1.values
print(lat)
print(lon)
start_lat=lat[0]
end_lat=lat[len(lat)-1]
start_lon=lon[0]
end_lon=lon[len(lon)-1]

#Constants for maximum hail size calculation
rd = 287.15
rhosui = 60000/(rd*253.15)
rho = airdensity
inv_rho = 1/rho
rhofaci = (rhosui*inv_rho)**0.54
dD=0.001
Dmax_psd=0.15
Ncrit=5**-4
Rcrit=1**-3
ch=206.89
dh=0.6384
cx=(np.pi/6)*900
print(cx)
dx=3
c1=3.7
c2=0.3
c3=9.0
c4=6.5
c5=1.0
c6=6.5
nd=int(Dmax_psd/dD)
#Set up arrays for max hail size calculation
R_tail = xr.zeros_like(icemassmr)
maxHailSize = xr.zeros_like(icemassmr)
R_crit=xr.ones_like(icemassmr)
R_crit=R_crit*0.001
#Dmh=((rho*icemassmr)/(cx*icenumbermr))**(1/dx)
Dmh=((6*icemassmr)/(np.pi*rimedensity*icenumbermr))**(1/dx) #Calculate mean-mass diameter of hail
icedensity=rimedensity*rimedicefraction+(1-rimedicefraction)*917
mu=0
lam=((np.pi*icedensity*icenumbermr*gamma(mu+4))/((gamma(mu+1))*icemassmr))**(1/(mu+3))
n0=(icenumbermr)*(lam)**(mu+1)/(gamma(mu+1)) #Calculate n0 using relation used in P3

for i in range(nd,1,-1): #Loop to calculate max hail size
    Di  = i*dD
    V_h = rhofaci*(ch*Di**dh) #Calculate terminal fall speed
    R_tail = R_tail + V_h*n0*Di**mu*np.exp(-lam*Di)*dD
    maxHailSize = xr.where((R_tail > R_crit) & (Di > maxHailSize), Di, maxHailSize)

icemassmr=icemassmr*airdensity
icenumbermr=icenumbermr*airdensity
icemassmr=icemassmr*1000
maxHailSize=maxHailSize*1000 #Convert max hail size to mm
Dmh=Dmh*1000

masked_reflectivityimg = np.ma.masked_where(reflectivityimg <= 0, reflectivityimg)
masked_reflectivity = np.ma.masked_where(reflectivity <= 0, reflectivity)
masked_rimedicefraction = np.ma.masked_where(rimedicefraction <= 0, rimedicefraction)
masked_icenumbermr = np.ma.masked_where(icenumbermr <= 0.001, icenumbermr)
masked_icenumbermr = np.log(masked_icenumbermr)
masked_icenumbermr2 = np.ma.masked_where(icemassmr <= 0.01, masked_icenumbermr)
masked_icemassmr = np.ma.masked_where(icemassmr <= 0.01, icemassmr)
masked_rimedensity = np.ma.masked_where(rimedensity <= 0, rimedensity)
masked_maxhailsizea = np.ma.masked_where(maxHailSize <= 0, maxHailSize)
masked_maxhailsizeb = np.ma.masked_where(rimedicefraction < 0.7, masked_maxhailsizea)
masked_maxhailsizec = np.ma.masked_where(rimedensity < 700, masked_maxhailsizeb)
masked_dmh = np.ma.masked_where(Dmh <= 0, Dmh)

#Print lat and lon values to verify location of cross sections
print(reflectivity.lat_1.values)
print(reflectivity.lon_1.values)

# Get lon and level for the cross section
lon = np.linspace(0,150,39)  # 1D lon array at chosen latitude
geometricheight=geometricheight/1000
freezinglevelheight=freezinglevelheight/1000
lev = geometricheight[:,0]
print(lev)
#lev=lev.reindex(level1=lev.level1[::-1])
# Create 2D meshgrid for plotting
LON, LEV = np.meshgrid(lon, lev)

lev2= verticalvelocity[0:39,0]
print(lev2)
LON2, LEV2 = np.meshgrid(lon,lev2)

cmap = ListedColormap(new_colors, 256)
cmap.set_bad(color='white')  # Will be used for masked values
hail = mcolors.LinearSegmentedColormap.from_list('hail', [(0, "blue"), (0.125, "blue"), (0.125, "cyan"), (0.25, "cyan"), (0.25, "green"), (0.375, "green"), (0.375, "lime"), (0.5, "lime"), (0.5, "yellow"), (0.625, "yellow"), (0.625, "orange"), (0.75, "orange"), (0.75, "red"), (0.875, "red"), (0.875, "purple"), (1.0, "purple")])
num_colors = 10  # Number of discrete colors
magma = plt.cm.get_cmap('magma', num_colors)  # Get a discrete version of 'magma'
num_colors = 8
magma2= plt.cm.get_cmap('magma', num_colors)
num_colors = 6
magma3= plt.cm.get_cmap('magma', num_colors)
num_colors = 9
magma4= plt.cm.get_cmap('magma', num_colors)
num_colors = 8
vvel1= plt.cm.get_cmap('seismic', num_colors)

levels = [0, 0.5, 1, 2, 3, 4, 5, 6]

magma5 = plt.get_cmap('magma', len(levels)-1)

norm= BoundaryNorm(levels, magma5.N)

fig = plt.figure(figsize=(17, 11))

gs = gridspec.GridSpec(2, 4, figure=fig, height_ratios=[1, 1])

# Reflectivity
ax2 = fig.add_subplot(gs[0,0])
cf2 = ax2.pcolormesh(LON, LEV, masked_reflectivity, shading='auto', cmap='NWSRef', vmin=0, vmax=70)
ax2.set_title('Reflectivity', fontsize=22)
cbar2=fig.colorbar(cf2, ax=ax2, shrink=0.95)
ax2.plot(lon, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
#ax2.set_xticks([0,100,200],labels=['0','100','200'],fontsize=20)
ax2.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax2.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar2.set_ticks([0,10,20,30,40,50,60,70], labels=['0','10','20','30','40','50','60','70'], fontsize=20)
ax2.set_xlabel('x [km]', fontsize=20)
ax2.set_xlim(0,150)
ax2.set_ylim(0,12)
ax2.set_ylabel('Height (km AGL)', fontsize=20)
ax2.text(0.045, 0.965, 'a)', transform=ax2.transAxes, fontsize=22, va='top')
# Ice number mixing ratio
ax3 = fig.add_subplot(gs[0,1])
cf3 = ax3.pcolormesh(LON, LEV, masked_icenumbermr2, shading='auto', cmap=magma2, vmin=0, vmax=16)
ax3.set_title('Number Concentration', fontsize=22)
ax3.plot(lon, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
#ax3.set_xticks([0,100,200],labels=['0','100','200'],fontsize=20)
ax3.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax3.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
ax3.set_xlabel('x [km]', fontsize=20)
ax3.set_xlim(0,150)
ax3.set_ylim(0,12)
cbar3=fig.colorbar(cf3, ax=ax3, shrink=0.95)
cbar3.set_ticks([0,2,4,6,8,10,12,14,16], labels=['0','2','4','6','8','10','12','14','16'], fontsize=20)
ax3.text(0.045, 0.965, 'b)', transform=ax3.transAxes, fontsize=22, va='top')
# Ice mass mixing ratio
ax4 = fig.add_subplot(gs[0,2])
cf4 = ax4.pcolormesh(LON, LEV, masked_icemassmr, shading='auto', cmap=magma5, norm=norm)
ax4.set_title('Total Mass Content', fontsize=22)
cbar4=fig.colorbar(cf4, ax=ax4,shrink=0.95)
ax4.plot(lon, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
#ax4.set_xticks([0,100,200],labels=['0','100','200'],fontsize=20)
ax4.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax4.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar4.set_ticks([0,0.5,1,1.5,2,3,4,5,6], labels=['0','0.5','1','1.5','2','3','4','5','6'], fontsize=20)
ax4.set_xlabel('x [km]', fontsize=20)
ax4.set_xlim(0,150)
ax4.set_ylim(0,12)
ax4.text(0.045, 0.965, 'c)', transform=ax4.transAxes, fontsize=22, va='top')
# Rime Fraction
ax5 = fig.add_subplot(gs[0,3])
cf5 = ax5.pcolormesh(LON, LEV, masked_rimedicefraction, shading='auto', cmap=magma, vmin=0, vmax=1)
ax5.set_title('Rime Fraction', fontsize=22)
cbar5=fig.colorbar(cf5, ax=ax5, shrink=0.95)
ax5.plot(lon, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
#ax5.set_xticks([0,100,200],labels=['0','100','200'],fontsize=20)
ax5.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax5.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar5.set_ticks([0,0.2,0.4,0.6,0.8,1.0], labels=['0','0.2','0.4','0.6','0.8','1.0'], fontsize=20)
ax5.set_xlabel('x [km]', fontsize=20)
ax5.set_xlim(0,150)
ax5.set_ylim(0,12)
ax5.text(0.045, 0.965, 'd)', transform=ax5.transAxes, fontsize=22, va='top')
# Rime Ice Density
ax6 = fig.add_subplot(gs[1,0])
cf6 = ax6.pcolormesh(LON, LEV, masked_rimedensity, shading='auto', cmap=magma4, vmin=0, vmax=900)
ax6.set_title('Rime Density', fontsize=22)
cbar6=fig.colorbar(cf6, ax=ax6, shrink=0.95)
ax6.plot(lon, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
#ax6.set_xticks([0,100,200],labels=['0','100','200'],fontsize=20)
ax6.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax6.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar6.set_ticks([0,100,200,300,400,500,600,700,800,900],labels=['0','100','200','300','400','500','600','700','800','900'],fontsize=20)
ax6.set_xlabel('x [km]', fontsize=20)
ax6.set_xlim(0,150)
ax6.set_ylabel('Height (km AGL)', fontsize=20)
ax6.set_ylim(0,12)
ax6.text(0.045, 0.965, 'e)', transform=ax6.transAxes, fontsize=22, va='top')
# Maximum hail size
ax7 = fig.add_subplot(gs[1,1])
cf7 = ax7.pcolormesh(LON, LEV, masked_maxhailsizec, shading='auto', cmap=hail, vmin=0, vmax=16)
ax7.set_title('Maximum Hail Size', fontsize=22)
ax7.plot(lon, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
#ax7.set_xticks([0,100,200],labels=['0','100','200'],fontsize=20)
ax7.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax7.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
ax7.set_xlabel('x [km]', fontsize=20)
ax7.set_xlim(0,150)
ax7.set_ylim(0,12)
cbar7=fig.colorbar(cf7, ax=ax7, shrink=0.95)
cbar7.set_ticks([0,2,4,6,8,10,12,14,16], labels=['0','2','4','6','8','10','12','14','16'], fontsize=20)
ax7.text(0.045, 0.965, 'f)', transform=ax7.transAxes, fontsize=22, va='top')
# Mean-mass diameter of hail
ax8 = fig.add_subplot(gs[1, 2])
cf8 = ax8.pcolormesh(LON, LEV, masked_dmh, shading='auto', cmap=hail, vmin=0, vmax=4)
ax8.set_title('Mean-Mass Diameter', fontsize=22)
cbar8=fig.colorbar(cf8, ax=ax8, shrink=0.95)
ax8.plot(lon, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
#ax8.set_xticks([0,100,200],labels=['0','100','200'],fontsize=20)
ax8.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax8.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar8.set_ticks([0,0.5,1,1.5,2,2.5,3,3.5,4], labels=['0','0.5','1','1.5','2','2.5','3','3.5','4'], fontsize=20)
ax8.set_xlabel('x [km]', fontsize=20)
ax8.set_xlim(0,150)
ax8.set_ylim(0,12)
ax8.text(0.045, 0.965, 'g)', transform=ax8.transAxes, fontsize=22, va='top')
#Omega
ax9 = fig.add_subplot(gs[1,3])
cf9 = ax9.pcolormesh(LON, LEV, verticalvelocity, shading='auto', cmap='seismic', vmin=-20, vmax=20)
ax9.set_title('Vertical Velocity', fontsize=22)
cbar9=fig.colorbar(cf9, ax=ax9, shrink=0.95)
ax9.plot(lon, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
ax9.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=17)
#ax9.set_xticks([0,100,200],labels=['0','100','200'],fontsize=20)
ax9.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar9.set_ticks([-20,-15,-10,-5,0,5,10,15,20], labels=['-20','-15','-10','-5','0','5','10','15','20'], fontsize=20)
ax9.set_xlabel('x [km]', fontsize=20)
ax9.set_xlim(0,150)
ax9.set_ylim(0,12)
ax9.text(0.045, 0.965, 'h)', transform=ax9.transAxes, fontsize=22, va='top')
#Plot figure, save figure, and show figure
fig.tight_layout()
plt.savefig(''+date+'/hrdps25icempcrossection1_'+ date +'_'+ hour +'_'+ half_hour +'_test.png')
plt.show()
