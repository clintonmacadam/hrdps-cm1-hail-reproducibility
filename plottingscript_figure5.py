#Calculates the maximum observable hail size from HRDPS 2.5km model output
#Import needed libraries
import xarray as xr
import h5py
from scipy.special import gamma
from math import tanh
import numpy as np
import matplotlib.pyplot as plt
import fstd2nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import colors as mcolors
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from PIL import Image, ImageDraw
from numba import njit, prange
import gc
import matplotlib.ticker as mticker

date="2024080512"
hour=['012','013','014','015','016']
half_hour="00"
half_hour_selector=1

testarray=hrdps25data=xr.open_dataset(r'/home/clintonmacadam/data/'+date+'_012', engine='fstd')
testarray=testarray['NTI1']
testarray=testarray.sel(level1=0.92, method='nearest')
testarray2=testarray[0,325:650,500:850]
maxHailSizeswath = xr.zeros_like(testarray2)

for i in range(5):
    hr=hour[i]
    date2="20240806"
    hr2=int(hr)-12
    if hr2<0:
        date2="20240805"
    else:
        date2=date2
    if hr2<0:
        hr2=hr2+24
    hr2=str(hr2)
    hrdps25data=xr.open_dataset(r'/home/clintonmacadam/data/'+date+'_'+hr+'', engine='fstd')
    #Select individual variables of interest from the model output file
    #Select individual variables of interest from the model output file
    nit=hrdps25data['NTI1']
    nit=nit[0,:,325:650,500:850]
    qx=hrdps25data['QTI1']
    qx=qx[half_hour_selector,:,325:650,500:850]
    rimedicemassmr = hrdps25data['QMI1']
    rimedicemassmr = rimedicemassmr[0,:,325:650,500:850]
    rimedicefraction = rimedicemassmr/qx

    rimevolumemr = hrdps25data['BMI1']
    rimevolumemr = rimevolumemr[0,:,325:650,500:850]

    rimedensity=rimedicemassmr/rimevolumemr
    
    specifichumidity=hrdps25data['HU']
    specifichumidity=specifichumidity[half_hour_selector,:,325:650,500:850]
    mixingratio=specifichumidity*(1-0.622)
    rainmassmr = hrdps25data['MPQR']
    rainmassmr = rainmassmr[half_hour_selector,:,325:650,500:850]
    cloudmassmr = hrdps25data['MPQC']
    cloudmassmr = cloudmassmr[half_hour_selector,:,325:650,500:850]
    liquidmassmr = rainmassmr+cloudmassmr
    ro=6356766 #radius of the earth, m
    a=29.3 #m/K
    rd=0.287053 #kPa K-1 m3 kg-1
    temperature = hrdps25data['TT']
    temperature = temperature[half_hour_selector,:,325:650,500:850]
    geopotentialheight = hrdps25data['GZ']
    geopotentialheight = geopotentialheight[half_hour_selector,:,325:650,500:850]
    geopotentialheight = geopotentialheight*10
    geometricheight = (ro*geopotentialheight)/(ro-geopotentialheight)
    geometricheight = xr.concat([geometricheight.isel(level2=slice(1, 125, 2)), geometricheight.isel(level2=124)], dim="level2")
    sfcpressure = hrdps25data['P0']
    sfcpressure = sfcpressure[half_hour_selector,325:650,500:850]

    pressure=xr.zeros_like(temperature)
    pressure[62,:]=sfcpressure

    tempcalc=temperature+273.15 #Convert temperature to Kelvin
    virtualtemp=tempcalc*(1+0.61*mixingratio-liquidmassmr-qx) #Calculate virtual temperature

    i=62
    while i>0: #Calculate pressure on model levels using the hypsometric equation
        dz = geometricheight[i-1, :] - geometricheight[i, :]
        pressure[i-1,:]= pressure[i,:]*np.exp(-dz/(a* (0.5 * (virtualtemp[i,:] + virtualtemp[i-1,:]))))
        i=i-1

    pressure=pressure/10 #Convert pressure to kPa
    airdensity=pressure/(rd*virtualtemp) #Calculate air density

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
    dx=3
    c1=3.7
    c2=0.3
    c3=9.0
    c4=6.5
    c5=1.0
    c6=6.5
    #nd=int(Dmax_psd/dD)
    #Set up arrays for max hail size calculation
    dD=0.001
    nd=int(Dmax_psd/dD)
    R_tail = xr.zeros_like(qx)
    maxHailSize = xr.zeros_like(qx)
    R_crit=xr.ones_like(qx)
    R_crit=R_crit*0.001
    #Dmh=((rho*qx)/(cx*nit))**(1/dx)
    Dmh=((6*qx)/(np.pi*rimedensity*nit))**(1/dx) #Calculate mean-mass diameter of hail
    icedensity=rimedensity*rimedicefraction+(1-rimedicefraction)*917
    mu=0
    lam=((np.pi*icedensity*nit*gamma(mu+4))/((gamma(mu+1))*qx))**(1/(mu+3))
    #mu1=c1*np.tanh(c2*(Dmh-c3))+c4
    #mu2=c5*Dmh-c6
    #mu=xr.where(Dmh>=0.008,mu2,mu1)
    #lam=(((gamma(1+dx+mu))/(gamma(1+mu)))*((cx*nit)/(rho*qx)))**(1/dx) #using relation from Milbrandt et al. 2005a
    n0=(nit)*(lam)**(mu+1)/(gamma(mu+1)) #Calculate n0 using relation used in P3
    
    #Integrate over the ice PSD and determine maximum hail size using Rcrit method
    for i in range(nd,1,-1):
        Di  = i*dD
        V_h = rhofaci*(ch*Di**dh)
        R_tail = R_tail + V_h*n0*Di**mu*np.exp(-lam*Di)*dD
        maxHailSize = xr.where((R_tail > R_crit) & (Di > maxHailSize), Di, maxHailSize)

    #Select a model level to plot the maximum hail size
    maxHailSize=maxHailSize.sel(level1=0.9975, method='nearest')
    maxHailSize=maxHailSize*1000
    maxHailSizeswath=xr.ufuncs.maximum(maxHailSizeswath,maxHailSize)

    target_lat=-1.08
    target_lon1=600
    target_lon2=657
    data = testarray[:, target_lon1:target_lon2]
    data = data.sel(rlat1=target_lat, method='nearest')
    lat=data.lat_1.values
    lon=data.lon_1.values
    start_lat=lat[0]
    end_lat=lat[len(lat)-1]
    start_lon=lon[0]
    end_lon=lon[len(lon)-1]

    cmap=plt.get_cmap('RdYlGn_r')
    new_colors = cmap(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 1])  # RGBA for white
    
    # Mask where data == 0
    masked_data = np.ma.masked_where(maxHailSizeswath <= 0, maxHailSizeswath)
    #print(masked_data)
    cmap = ListedColormap(new_colors, 256)
    cmap.set_bad(color='white')  # Will be used for masked values
    hail = mcolors.LinearSegmentedColormap.from_list('hail', [(0, "blue"), (0.125, "blue"), (0.125, "cyan"), (0.25, "cyan"), (0.25, "green"), (0.375, "green"), (0.375, "lime"), (0.5, "lime"), (0.5, "yellow"), (0.625, "yellow"), (0.625, "orange"), (0.75, "orange"), (0.75, "red"), (0.875, "red"), (0.875, "purple"), (1.0, "purple")])

    #Plot the maximum hail size across the domain of interest
    fig, ax = plt.subplots(figsize=(17, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    cf = ax.pcolormesh(maxHailSize['lon_1'], maxHailSize['lat_1'], masked_data, cmap=hail, transform=ccrs.PlateCarree(), vmin=0, vmax=16)
    #plt.plot([start_lon, end_lon],[start_lat, end_lat], color='blue', linewidth=2, transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,ax=ax, shrink=0.65)
    cbar.set_label('Maximum Hail Size (mm)', fontsize=22, labelpad=12)
    cbar.set_ticks([0,2,4,6,8,10,12,14,16], labels=['0','2','4','6','8','10','12','14','16'], fontsize=22)
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES)
    ax.set_extent([-115.4,-111.5,50.5,52.5])
    ax.set_aspect('auto')
    #fig.supxlabel('Longitude', fontsize=22,y=0.14,x=0.43)
    #fig.supylabel('Latitude', fontsize=22,y=0.5,x=0.02)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, alpha=0.5, linestyle='-', color='None')
    gl.xlabel_style = {'size': 22}
    gl.ylabel_style = {'size': 22}
    gl.top_labels=False   # suppress top labels
    gl.right_labels=False # suppress right labels
    gl.xlocator = mticker.FixedLocator([-115,-114,-113,-112])
    gl.ylocator = mticker.FixedLocator([50.5,51,51.5,52,52.5])
    ax.plot(-114.07, 51.04, marker='o', markersize=24, color='r', transform=ccrs.PlateCarree())
    ax.text(-113.99, 50.995, 'Calgary', fontsize=24, weight='bold', transform=ccrs.PlateCarree())
    ax.plot(-114.105, 51.79, marker='o', markersize=12, color='r', transform=ccrs.PlateCarree())
    ax.text(-114.035, 51.757, 'Olds', fontsize=24, weight='bold', transform=ccrs.PlateCarree())
    ax.plot(-113.81, 52.27, marker='o', markersize=18, color='r', transform=ccrs.PlateCarree())
    ax.text(-113.74, 52.241, 'Red Deer', fontsize=24, weight='bold', transform=ccrs.PlateCarree())
    ax.plot(-113.5, 53.54, marker='o', markersize=10, color='b', transform=ccrs.PlateCarree())
    ax.text(-113.43, 53.54, 'Edmonton', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    ax.plot(-114.4669, 51.2557, marker='o', markersize=24, color='b', transform=ccrs.PlateCarree())
    ax.text(-114.5169, 51.32, '27', fontsize=24, weight='bold', transform=ccrs.PlateCarree())
    ax.plot(-114.3282, 51.2541, marker='o', markersize=24, color='b', transform=ccrs.PlateCarree())
    ax.text(-114.3782, 51.32, '37', fontsize=24, weight='bold', transform=ccrs.PlateCarree())
    #ax0.plot(-114.3274, 51.2310, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax0.text(-114.3274, 51.23, '36', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax0.plot(-114.6692, 51.2845, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
    #ax0.text(-114.7192, 51.35, '25', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
    #ax0.plot(-114.1180, 51.2123, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
    #ax0.text(-114.1680, 51.28, '34', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
    ax.plot(-114.0486, 51.1869, marker='o', markersize=24, color='b', transform=ccrs.PlateCarree())
    ax.text(-114.0986, 51.255, '39', fontsize=24, weight='bold', transform=ccrs.PlateCarree())
    ax.plot(-113.8887, 51.1477, marker='o', markersize=24, color='b', transform=ccrs.PlateCarree())
    ax.text(-113.9387, 51.21, '33', fontsize=24, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-114.5961, 52.2316, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-114.5961, 52.3, '56', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-113.74265, 52.0264, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-113.74265, 52.1, '79', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-113.23995, 51.9405, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-113.23995, 52.01, '43', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-114.2874, 52.1809, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-114.2874, 52.25, '106', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-114.135, 49.883, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-114.135, 49.89, '17', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-110.386, 51.98, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-110.386, 52.04, '50', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-114.44, 53.34, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-114.44, 53.4, '10', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-114.85, 52.46, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-114.85, 52.53, '24', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-115.54, 52.49, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-115.54, 52.55, '15', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-114.15, 52.57, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-114.15, 52.64, '23', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-114.89, 52.57, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-114.89, 52.63, '53', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-114.94, 53.00, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-114.94, 53.06, '51', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-114.78, 52.95, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-114.78, 53.01, '38', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-115.88, 53.10, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-115.88, 53.16, '70', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-115.22, 53.27, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-115.22, 53.33, '38', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    #ax.plot(-115.15, 52.26, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
    #ax.text(-115.15, 52.32, '24', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
    plt.title('Surface Maximum Hail Size', fontsize=28)
    plt.tight_layout()
    plt.savefig(''+date+'/hrdps25maxhailsize_'+ date +'_'+ hr +'_'+ half_hour +'.png')

image_list=[''+ date +'/hrdps25maxhailsize_'+ date +'_012_00.png', ''+ date +'/hrdps25maxhailsize_'+ date +'_013_00.png', ''+ date +'/hrdps25maxhailsize_'+ date +'_014_00.png', ''+ date +'/hrdps25maxhailsize_' + date +'_015_00.png', ''+ date +'/hrdps25maxhailsize_'+ date +'_016_00.png']
images=[Image.open(image) for image in image_list]
images=[Image.open(image) for image in image_list]
images[0].save(''+ date +'/maxhailsizeloop_'+ date +'.gif', save_all=True, append_images=images[1:], duration=600, loop=0)
