#ice microphysical cross sections
#import needed libraries
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pyart
from matplotlib import colors as mcolors
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from scipy.special import gamma
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm
import matplotlib.gridspec as gridspec
import metpy.calc as mpcalc
from metpy.interpolate import cross_section

#specify date, hour, and half-hour
date="20240805"
timestep_selector=12
time=str(timestep_selector*5)

#Load data
cm1output = xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_2500m_mphrdps.nc', engine='netcdf4')

# Choose nearest lat-lon line for cross-section (e.g., along constant latitude)
target_y=13
start_x=-200
end_x=200
start_y=target_y
end_y=target_y

# Set up colormap and white mask
cmap=plt.get_cmap('magma')
new_colors = cmap(np.linspace(0, 1, 256))
white = np.array([1, 1, 1, 1])  # RGBA for white

# Extract data along the specfied lat index across longitude and levels for variables of interest
reflectivity = cm1output['dbz']
reflectivityimg = reflectivity[timestep_selector,:,:,:]
reflectivityimg = reflectivityimg.max(dim="zh")

reflectivity = reflectivity[timestep_selector,:,:,start_x:end_x]
reflectivity = reflectivity.sel(zh=slice(0,12))
reflectivity = reflectivity.sel(yh=target_y, method='nearest')

icenumbermr1 = cm1output['ni1']
icenumbermr1 = icenumbermr1[timestep_selector,:,:,:]
icenumbermr = icenumbermr1[:,:,start_x:end_x]
icenumbermr= icenumbermr.sel(zh=slice(0,12))
icenumbermr = icenumbermr.sel(yh=target_y, method='nearest')
icenumbermr1=icenumbermr1[0,:,:]

icemassmr1 = cm1output['qi1']
icemassmr1 = icemassmr1[timestep_selector,:,:,:]
icemassmr = icemassmr1[:,:,start_x:end_x]
icemassmr = icemassmr.sel(zh=slice(0,12))
icemassmr = icemassmr.sel(yh=target_y, method='nearest')
icemassmr1=icemassmr1[0,:,:]

rimedicemassmr1 = cm1output['ri1']
rimedicemassmr1 = rimedicemassmr1[timestep_selector,:, :, :]
rimedicemassmr = rimedicemassmr1[:,:,start_x:end_x]
rimedicemassmr = rimedicemassmr.sel(zh=slice(0,12))
rimedicemassmr = rimedicemassmr.sel(yh=target_y, method='nearest')
rimedicemassmr1 = rimedicemassmr1[0,:,:]

rimedicefraction = rimedicemassmr/icemassmr

rimedicevolumemr1 = cm1output['bi1']
rimedicevolumemr1 = rimedicevolumemr1[timestep_selector,:, :, :]
rimedicevolumemr1 = rimedicevolumemr1[0,:,:]

icedensity1=rimedicemassmr1/rimedicevolumemr1

icedensity = cm1output['p3_rhoi1']
icedensity = icedensity[timestep_selector,:,:,start_x:end_x]
icedensity = icedensity.sel(zh=slice(0,12))
icedensity = icedensity.sel(yh=target_y, method='nearest')

theta = cm1output['th']
theta = theta[timestep_selector,:,:,start_x:end_x]
theta = theta.sel(zh=slice(0,12))
theta = theta.sel(yh=target_y, method='nearest')

pressure = cm1output['prs']
pressure = pressure[timestep_selector,:,:,start_x:end_x]
pressure = pressure.sel(zh=slice(0,12))
pressure = pressure.sel(yh=target_y, method='nearest')
pressurecalc = pressure/1000

vertvelocity = cm1output['w']
vertvelocity = vertvelocity[timestep_selector,:,:,start_x:end_x]
vertvelocity = vertvelocity.sel(zf=slice(0,12))
vertvelocity = vertvelocity.sel(yh=target_y, method='nearest')

mixingratio = cm1output['qv']
mixingratio = mixingratio[timestep_selector,:,:,start_x:end_x]
mixingratio = mixingratio.sel(zh=slice(0,12))
mixingratio = mixingratio.sel(yh=target_y, method='nearest')

p0=100000
constant=0.28571
tempcalc=theta/((p0/pressure)**constant)
temperature=tempcalc-273.15

es=0.6113*np.exp(((2.5*10**6)/(461))*((1/273.15)-(1/tempcalc)))
e=(pressurecalc*mixingratio)/(0.622 + mixingratio)
rh=e/es*100
wetbulbtemp=temperature*np.arctan(0.151977*(rh+8.313659)**(1/2))-4.686035+np.arctan(temperature+rh)-np.arctan(rh-1.676331)+(0.00391838*(rh)**(3/2)*np.arctan(0.023101*rh))

#Determine the level where the temperature is closest to 0C/the approximate freezing level
abs_diff = np.abs(wetbulbtemp-0)
min_idx = abs_diff.argmin(dim='zh')
min_idx = min_idx.compute().to_numpy()

geometricheight=temperature.zh.values
freezinglevelheight = geometricheight[min_idx] #Calculate freezing level height (approximate)
print(freezinglevelheight)

airdensity1 = cm1output['rho']
airdensity1 = airdensity1[timestep_selector,:,:,:]
airdensity = airdensity1[:,:,start_x:end_x]
airdensity = airdensity.sel(zh=slice(0,12))
airdensity = airdensity.sel(yh=target_y, method='nearest')
airdensity1 = airdensity1[0,:,:]

cx=(np.pi/6)*900
dx=3

dmh=((6*icemassmr)/(np.pi*icedensity*icenumbermr))**(1/dx)
dmh1=((6*icemassmr1)/(np.pi*icedensity1*icenumbermr1))**(1/dx)

def calc_max_hail_size(icemassmr,icenumbermr,icedensity,dmh,airdensity):
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
    N_tot=rho*icenumbermr
    Dmh=((6*icemassmr)/(np.pi*icedensity*icenumbermr))**(1/dx) #Calculate mean-mass diameter of hail
    print(Dmh.values)
    icedensity=icedensity*rimedicefraction+(1-rimedicefraction)*917
    mu=0
    lam=((np.pi*icedensity*icenumbermr*gamma(mu+4))/((gamma(mu+1))*icemassmr))**(1/(mu+3))
    n0=(icenumbermr)*(lam)**(mu+1)/(gamma(mu+1)) #Calculate n0 using relation used in P3

    for i in range(nd,1,-1): #Loop to calculate max hail size
        Di  = i*dD
        V_h = rhofaci*(ch*Di**dh) #Calculate terminal fall speed
        R_tail = R_tail + V_h*n0*Di**mu*np.exp(-lam*Di)*dD
        maxHailSize = xr.where((R_tail > R_crit) & (Di > maxHailSize), Di, maxHailSize)
    return maxHailSize

maxHailSize=calc_max_hail_size(icemassmr,icenumbermr,icedensity,dmh,airdensity)
maxHailSize_img=calc_max_hail_size(icemassmr1,icenumbermr1,icedensity1,dmh1,airdensity1)

maxHailSize_img = maxHailSize_img[:,:,0]

#maxHailSize = cm1output['p3_dhmax1']
#maxHailSize = maxHailSize[timestep_selector,:,:,start_x:end_x]
#maxHailSize = maxHailSize.sel(zh=slice(0,12))
#maxHailSize = maxHailSize.sel(yh=target_y, method='nearest')
maxHailSize=maxHailSize*1000 #Convert max hail size to mm
maxHailSize_img=maxHailSize_img*1000
icenumbermr=icenumbermr/airdensity
icemassmr=icemassmr/airdensity
icemassmr=icemassmr*1000
icenumbermr=icenumbermr*1000
dmh=dmh*1000

masked_reflectivityimg = np.ma.masked_where(reflectivityimg <= 0, reflectivityimg)
masked_maxHailSize_img = np.ma.masked_where(maxHailSize_img <= 0, maxHailSize_img)
masked_reflectivity = np.ma.masked_where(reflectivity <= 0, reflectivity)
masked_icenumbermr = np.ma.masked_where(icenumbermr <= 0, icenumbermr)
masked_icenumbermr = np.log(masked_icenumbermr)
masked_icenumbermr2 = np.ma.masked_where(icemassmr <= 0.01, masked_icenumbermr)
masked_icemassmr = np.ma.masked_where(icemassmr <= 0.01, icemassmr)
masked_rimedicefraction = np.ma.masked_where(rimedicefraction <= 0, rimedicefraction)
masked_rimedensity = np.ma.masked_where(icedensity <= 0, icedensity)
masked_rimedensity2 = np.ma.masked_where(icemassmr <= 0.01, masked_rimedensity)
masked_maxhailsizea = np.ma.masked_where(maxHailSize <= 0, maxHailSize)
masked_maxhailsizeb = np.ma.masked_where(rimedicefraction < 0.7, masked_maxhailsizea)
masked_maxhailsizec = np.ma.masked_where(icedensity < 700, masked_maxhailsizeb)
masked_dmh = np.ma.masked_where(dmh <= 0, dmh)
masked_dmh2 = np.ma.masked_where(icemassmr <= 0.01, masked_dmh)

#Print lat and lon values to verify location of cross sections
print(reflectivity.xh.values)
print(reflectivity.yh.values)

x=reflectivity.xh.values+50
# Get lon and level for the cross section
lon = reflectivity['xh']+50
lev = geometricheight
print(geometricheight.shape)
#lev=lev.reindex(level1=lev.level1[::-1])
# Create 2D meshgrid for plotting
LON, LEV = np.meshgrid(lon, lev)

# Plotting
fig = plt.figure(figsize=(17, 16.5))

gs = gridspec.GridSpec(3, 4, figure=fig, height_ratios=[1.2,1, 1])

hail = mcolors.LinearSegmentedColormap.from_list('hail', [(0, "blue"), (0.125, "blue"), (0.125, "cyan"), (0.25, "cyan"), (0.25, "green"), (0.375, "green"), (0.375, "lime"), (0.5, "lime"), (0.5, "yellow"), (0.625, "yellow"), (0.625, "orange"), (0.75, "orange"), (0.75, "red"), (0.875, "red"), (0.875, "purple"), (1.0, "purple")])
num_colors = 10  # Number of discrete colors
magma = plt.cm.get_cmap('magma', num_colors)  # Get a discrete version of 'magma'
num_colors = 8
magma2= plt.cm.get_cmap('magma', num_colors)
num_colors = 6
magma3= plt.cm.get_cmap('magma', num_colors)
num_colors = 9
magma4= plt.cm.get_cmap('magma', num_colors)

levels = [0, 1, 2, 4, 6, 8, 10, 12]

magma5 = plt.get_cmap('magma', len(levels)-1)

norm= BoundaryNorm(levels, magma5.N)

# Reflectivity image                                                                                                    
ax0 = fig.add_subplot(gs[0, 0:2])
ax0.set_aspect('auto')
cf0 = ax0.pcolormesh(reflectivityimg['xh'], reflectivityimg['yh'], masked_reflectivityimg, cmap='NWSRef', vmin=0, vmax=70)
ax0.plot([start_x, end_x],[start_y, end_y], color='blue', linewidth=2)

cbar0=plt.colorbar(cf0,ax=ax0)
cbar0.set_label('Reflectivity (dBZ)', fontsize=20)
cbar0.set_ticks([0,10,20,30,40,50,60,70], labels=['0','10','20','30','40','50','60','70'], fontsize=19)

ax0.set_title('Column Maximum Reflectivity', fontsize=22)
ax0.set_xlabel('x [km]', fontsize=20)
ax0.set_ylabel('y [km]', fontsize=20)
ax0.set_xlim(-60,60)
ax0.set_ylim(-60,60)
ax0.set_xticks([-60,-30,0,30,60],labels=['-60','-30','0','30','60'],fontsize=20)
ax0.set_yticks([-60,-30,0,30,60],labels=['-60','-30','0','30','60'],fontsize=20)
ax0.text(0.025, 0.97, 'a)', transform=ax0.transAxes, fontsize=22, va='top')

ax1 = fig.add_subplot(gs[0, 2:4])
ax1.set_aspect('auto')
cf1 = ax1.pcolormesh(maxHailSize_img['xh'], maxHailSize_img['yh'], masked_maxHailSize_img, cmap=hail, vmin=0, vmax=16)
#ax1.plot([start_x, end_x],[start_y, end_y], color='blue', linewidth=2)

cbar1=plt.colorbar(cf1,ax=ax1)
cbar1.set_label('Maximum Hail Size (mm)', fontsize=20)
cbar1.set_ticks([0,2,4,6,8,10,12,14,16], labels=['0','2','4','6','8','10','12','14','16'], fontsize=20)

ax1.set_title('Surface Maximum Hail Size', fontsize=22)
ax1.set_xlabel('x [km]', fontsize=20)
ax1.set_ylabel('y [km]', fontsize=20)
ax1.set_xlim(-60,60)
ax1.set_ylim(-60,60)
ax1.set_xticks([-60,-30,0,30,60],labels=['-60','-30','0','30','60'],fontsize=20)
ax1.set_yticks([-60,-30,0,30,60],labels=['-60','-30','0','30','60'],fontsize=20)
ax1.text(0.025, 0.97, 'b)', transform=ax1.transAxes, fontsize=22, va='top')

#fig.suptitle('Ice Microphysical Cross Sections [New P3: 2-Moment], 20240805 Calgary Hailstorm Case', fontsize=18)
# Reflectivity
ax2 = fig.add_subplot(gs[1, 0])
cf2 = ax2.pcolormesh(LON, LEV, masked_reflectivity, shading='auto', cmap='NWSRef', vmin=0, vmax=70)
ax2.set_title('Reflectivity', fontsize=22)
cbar2=fig.colorbar(cf2, ax=ax2, shrink=0.95)
ax2.plot(x, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
ax2.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax2.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar2.set_ticks([0,10,20,30,40,50,60,70], labels=['0','10','20','30','40','50','60','70'], fontsize=20)
ax2.set_ylabel('Height (m AGL)', fontsize=20)
ax2.set_xlim(0,150)
ax2.text(0.045, 0.965, 'c)', transform=ax2.transAxes, fontsize=22, va='top')
# Ice number mixing ratio
ax3 = fig.add_subplot(gs[1, 1])
cf3 = ax3.pcolormesh(LON, LEV, masked_icenumbermr2, shading='auto', cmap=magma2, vmin=0, vmax=24)
ax3.set_title('Number Concentration', fontsize=22)
cbar3=fig.colorbar(cf3, ax=ax3, shrink=0.95)
ax3.plot(x, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
ax3.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax3.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar3.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=20)
ax3.set_xlim(0,150)
ax3.text(0.045, 0.965, 'd)', transform=ax3.transAxes, fontsize=22, va='top')
# Ice mass mixing ratio
ax4 = fig.add_subplot(gs[1, 2])
cf4 = ax4.pcolormesh(LON, LEV, masked_icemassmr, shading='auto', cmap=magma5, norm=norm)
ax4.set_title('Ice Mass Content', fontsize=22)
cbar4=fig.colorbar(cf4, ax=ax4,shrink=0.95)
ax4.plot(x, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
ax4.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax4.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar4.set_ticks([0,1,2,4,6,8,10,12], labels=['0','1','2','4','6','8','10','12'], fontsize=20)
ax4.set_xlim(0,150)
ax4.text(0.045, 0.965, 'e)', transform=ax4.transAxes, fontsize=22, va='top')
# Rime Fraction
ax5 = fig.add_subplot(gs[1,3])
cf5 = ax5.pcolormesh(LON, LEV, masked_rimedicefraction, shading='auto', cmap=magma, vmin=0, vmax=1)
ax5.set_title('Rime Fraction', fontsize=22)
cbar5=fig.colorbar(cf5, ax=ax5, shrink=0.95)
ax5.plot(x, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
ax5.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax5.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar5.set_ticks([0,0.2,0.4,0.6,0.8,1.0], labels=['0','0.2','0.4','0.6','0.8','1.0'], fontsize=20)
ax5.set_xlim(0,150)
ax5.text(0.045, 0.965, 'f)', transform=ax5.transAxes, fontsize=22, va='top')
# Rime Ice Density
ax6 = fig.add_subplot(gs[2, 0])
cf6 = ax6.pcolormesh(LON, LEV, masked_rimedensity2, shading='auto', cmap=magma4, vmin=0, vmax=900)
ax6.set_title('Rime Density', fontsize=22)
cbar6=fig.colorbar(cf6, ax=ax6, shrink=0.95)
ax6.plot(x, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
ax6.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax6.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
ax6.set_xlabel('x [km]', fontsize=20)
ax6.set_ylabel('Height (m AGL)', fontsize=20)
cbar6.set_ticks([0,100,200,300,400,500,600,700,800,900],labels=['0','100','200','300','400','500','600','700','800','900'],fontsize=20)
ax6.set_xlim(0,150)
ax6.text(0.045, 0.965, 'g)', transform=ax6.transAxes, fontsize=22, va='top')
# Maximum hail size
ax7 = fig.add_subplot(gs[2,1])
cf7 = ax7.pcolormesh(LON, LEV, masked_maxhailsizec, shading='auto', cmap=hail, vmin=0, vmax=16)
ax7.set_title('Maximum Hail Size', fontsize=22)
cbar7=fig.colorbar(cf7, ax=ax7, shrink=0.95)
ax7.plot(x, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
ax7.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax7.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar7.set_ticks([0,2,4,6,8,10,12,14,16], labels=['0','2','4','6','8','10','12','14','16'], fontsize=20)
ax7.set_xlabel('x [km]', fontsize=20)
ax7.set_xlim(0,150)
ax7.text(0.045, 0.965, 'h)', transform=ax7.transAxes, fontsize=22, va='top')
# Mean-mass diameter of hail
ax8 = fig.add_subplot(gs[2,2])
cf8 = ax8.pcolormesh(LON, LEV, masked_dmh2, shading='auto', cmap=hail, vmin=0, vmax=4)
ax8.set_title('Mean-Mass Diameter', fontsize=22)
cbar8=fig.colorbar(cf8, ax=ax8, shrink=0.95)
ax8.plot(x, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
ax8.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax8.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar8.set_ticks([0,0.5,1,1.5,2,2.5,3,3.5,4], labels=['0','0.5','1','1.5','2','2.5','3','3.5','4'], fontsize=20)
ax8.set_xlabel('x [km]', fontsize=20)
ax8.set_xlim(0,150)
ax8.text(0.045, 0.965, 'i)', transform=ax8.transAxes, fontsize=22, va='top')
#Omega
ax9 = fig.add_subplot(gs[2,3])
cf9 = ax9.pcolormesh(LON, LEV, vertvelocity, shading='auto', cmap='seismic', vmin=-20, vmax=20)
ax9.set_title('Vertical Velocity', fontsize=22)
cbar9=fig.colorbar(cf9, ax=ax9, shrink=0.95)
ax9.plot(x, freezinglevelheight, color='red', linewidth=3, label='0°C Level')
ax9.set_xticks([0,50,100,150],labels=['0','50','100','150'],fontsize=20)
ax9.set_yticks([0,2,4,6,8,10],labels=['0','2','4','6','8','10'],fontsize=20)
cbar9.set_ticks([-20,-15,-10,-5,0,5,10,15,20], labels=['-20','-15','-10','-5','0','5','10','15','20'], fontsize=20)
ax9.set_xlabel('x [km]', fontsize=20)
ax9.set_xlim(0,150)
ax9.text(0.045, 0.965, 'j)', transform=ax9.transAxes, fontsize=22, va='top')
#Plot figure, save figure, and show figure
fig.tight_layout()
plt.savefig(''+date+'/cm1icempcrossectionfigure_'+ date +'_'+ time +'_2500m_mphrdps.png')
plt.show()
