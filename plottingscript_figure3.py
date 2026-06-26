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
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe

#specify date, hour, and half-hour
date="2024080512"
hour="014"
half_hour="00"
half_hour_selector=1 #0 or 1, time1 index

date2="2023072312"
hour2="011"

# Load data
hrdps25data = xr.open_dataset(r'/home/clintonmacadam/data/'+date+'_'+hour+'', engine="fstd")
hrdps25data2 = xr.open_dataset(r'/home/clintonmacadam/data/'+date2+'_'+hour2+'', engine="fstd")

# Choose nearest lat-lon line for cross-section (e.g., along constant latitude)
target_lat=-2.87
start_lon=622
end_lon=679

target_lat2=-0.89
start_lon2=619
end_lon2=658

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
reflectivity2 = hrdps25data2['ZET']
reflectivityimg = hrdps25data['ZEC']
reflectivityimg = reflectivityimg[1,:,:]
reflectivityimg2 = hrdps25data2['ZEC']
reflectivityimg2 = reflectivityimg2[1,:,:]

reflectivity = reflectivity[half_hour_selector,:,:,start_lon:end_lon]
reflectivity = reflectivity.sel(level1=slice(0.15,0.9975))
reflectivity = reflectivity.sel(rlat1=target_lat, method='nearest')

reflectivity2 = reflectivity2[half_hour_selector,:,:,start_lon2:end_lon2]
reflectivity2 = reflectivity2.sel(level1=slice(0.15,0.9975))
reflectivity2 = reflectivity2.sel(rlat1=target_lat2, method='nearest')

lat=reflectivity.lat_1.values
lon=reflectivity.lon_1.values
print(lat)
print(lon)
start_lat=lat[0]
end_lat=lat[len(lat)-1]
start_lon=lon[0]
end_lon=lon[len(lon)-1]

lat2=reflectivity2.lat_1.values
lon2=reflectivity2.lon_1.values
print(lat2)
print(lon2)
start_lat2=lat2[0]
end_lat2=lat2[len(lat2)-1]
start_lon2=lon2[0]
end_lon2=lon2[len(lon2)-1]

masked_reflectivityimg = np.ma.masked_where(reflectivityimg <= 0, reflectivityimg)
masked_reflectivityimg2 = np.ma.masked_where(reflectivityimg2 <= 0, reflectivityimg2)

#Print lat and lon values to verify location of cross sections
print(reflectivity.lat_1.values)
print(reflectivity.lon_1.values)

fig = plt.figure(figsize=(17, 13))

gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[1.01,1], height_ratios=[1,1])

ax0 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())

cf0 = ax0.pcolormesh(hrdps25data['lon_1'], hrdps25data['lat_1'], masked_reflectivityimg, cmap='NWSRef', transform=ccrs.PlateCarree(), vmin=0, vmax=70)
ax0.plot([start_lon, end_lon],[start_lat, end_lat], color='blue', linewidth=7, transform=ccrs.PlateCarree())

ax0.margins(x=0, y=0)
cbar0=plt.colorbar(cf0,ax=ax0, fraction=0.035, pad=0.02, shrink=1.15)
cbar0.set_ticks([0,10,20,30,40,50,60,70], labels=['0','10','20','30','40','50','60','70'], fontsize=20)

ax0.add_feature(cfeature.BORDERS)
ax0.add_feature(cfeature.STATES)
ax0.set_title('Column Maximum Reflectivity', fontsize=24)
ax0.set_extent([-115.5,-113,50.5,52.5])
ax0.set_aspect('auto')
ax0.set_xlabel('Longitude', fontsize='20')
ax0.set_ylabel('Latitude', fontsize='20')

gl = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, alpha=0.5, linestyle='-', color='None')
gl.top_labels=False # suppress top labels
gl.right_labels=False # suppress right labels
gl.xlabel_style = {'size': 20}
gl.ylabel_style = {'size': 20}
gl.xlocator = mticker.FixedLocator([-115,-114])
gl.ylocator = mticker.FixedLocator([50.5,51,51.5,52,52.5])

ax0.plot(-114.07, 51.04, marker='o', markersize=20, color='r', transform=ccrs.PlateCarree())
ax0.text(-113.985, 51.02, 'Calgary', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
ax0.plot(-114.105, 51.79, marker='o', markersize=20, color='r', transform=ccrs.PlateCarree())
ax0.text(-114.05, 51.76, 'Olds', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
ax0.plot(-113.81, 52.27, marker='o', markersize=20, color='r', transform=ccrs.PlateCarree())
ax0.text(-113.75, 52.24, 'Red Deer', fontsize=18, weight='bold', transform=ccrs.PlateCarree())

ax0.plot(-114.4669, 51.2557, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax0.text(-114.5169, 51.32, '27', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
ax0.plot(-114.3282, 51.2541, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax0.text(-114.3782, 51.32, '37', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
#ax0.plot(-114.3274, 51.2310, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
#ax0.text(-114.3274, 51.23, '36', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
#ax0.plot(-114.6692, 51.2845, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
#ax0.text(-114.7192, 51.35, '25', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
#ax0.plot(-114.1180, 51.2123, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
#ax0.text(-114.1680, 51.28, '34', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
ax0.plot(-114.0486, 51.1869, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax0.text(-114.0986, 51.255, '39', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
ax0.plot(-113.8887, 51.1477, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax0.text(-113.9387, 51.21, '33', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
ax0.text(0.01, 0.975, 'a)', transform=ax0.transAxes, fontsize=24, va='top')

ax1 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())

img = mpimg.imread(r'/home/clintonmacadam/nhpanalysis/'+date+'/calgary_radarimg.png')
ax1.imshow(
    img,
    origin="upper",
    extent=[-115.5, -113, 50.5, 52.5],  # [lon_min, lon_max, lat_min, lat_max]
    transform=ccrs.PlateCarree(),
    interpolation="nearest",
    zorder=1
)
ax1.set_title('Observed Radar Reflectivity', fontsize=24)
ax1.set_aspect('auto')

gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, alpha=0.5, linestyle='-', color='None')
gl1.top_labels=False # suppress top labels
gl1.right_labels=False # suppress right labels
gl1.xlabel_style = {'size': 20}
gl1.ylabel_style = {'size': 20}
gl1.xlocator = mticker.FixedLocator([-115,-114])
gl1.ylocator = mticker.FixedLocator([50.5,51,51.5,52,52.5])

ax1.plot(-114.07, 51.04, marker='o', markersize=20, color='r', transform=ccrs.PlateCarree())
ax1.text(-113.985, 51.02, 'Calgary', fontsize=18, color='white', weight='bold', transform=ccrs.PlateCarree())
ax1.plot(-114.105, 51.79, marker='o', markersize=20, color='r', transform=ccrs.PlateCarree())
ax1.text(-114.05, 51.76, 'Olds', fontsize=18, color='white', weight='bold', transform=ccrs.PlateCarree())
ax1.plot(-113.81, 52.27, marker='o', markersize=20, color='r', transform=ccrs.PlateCarree())
ax1.text(-113.75, 52.24, 'Red Deer', fontsize=18, color='white', weight='bold', transform=ccrs.PlateCarree())

ax1.plot(-114.4669, 51.2557, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax1.text(-114.5169, 51.32, '27', fontsize=18, weight='bold', color='white', transform=ccrs.PlateCarree())
ax1.plot(-114.3282, 51.2541, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax1.text(-114.3782, 51.32, '37', fontsize=18, weight='bold', color='white', transform=ccrs.PlateCarree())
#ax0.plot(-114.3274, 51.2310, marker='o', markersize=8, color='r', transform=ccrs.PlateCarree())
#ax0.text(-114.3274, 51.23, '36', fontsize=10, weight='bold', transform=ccrs.PlateCarree())
#ax0.plot(-114.6692, 51.2845, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
#ax0.text(-114.7192, 51.35, '25', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
#ax1.plot(-114.1180, 51.2123, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
#ax1.text(-114.1680, 51.28, '34', fontsize=18, weight='bold', color='white', transform=ccrs.PlateCarree())
ax1.plot(-114.0486, 51.1869, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax1.text(-114.0986, 51.255, '39', fontsize=18, weight='bold', color='white', transform=ccrs.PlateCarree())
ax1.plot(-113.8887, 51.1477, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax1.text(-113.9387, 51.21, '33', fontsize=18, weight='bold', color='white', transform=ccrs.PlateCarree())
ax1.text(0.01, 0.975, 'b)', transform=ax1.transAxes, color='white', fontsize=24, va='top')

ax2 = fig.add_subplot(gs[1, 0], projection=ccrs.PlateCarree())
cf2 = ax2.pcolormesh(hrdps25data2['lon_1'], hrdps25data2['lat_1'], masked_reflectivityimg2, cmap='NWSRef', transform=ccrs.PlateCarree(), vmin=0, vmax=70)
ax2.plot([start_lon2, end_lon2],[start_lat2, end_lat2], color='blue', linewidth=7, transform=ccrs.PlateCarree())

ax2.margins(x=0, y=0)
cbar2=plt.colorbar(cf2,ax=ax2, fraction=0.035, pad=0.02, shrink=1.15)
cbar2.set_ticks([0,10,20,30,40,50,60,70], labels=['0','10','20','30','40','50','60','70'], fontsize=20)

ax2.add_feature(cfeature.BORDERS)
ax2.add_feature(cfeature.STATES)
ax2.set_title('Column Maximum Reflectivity', fontsize=24)
ax2.set_extent([-118.5,-114,52.5,54.5])
ax2.set_aspect('auto')
ax2.set_xlabel('Longitude', fontsize='20')
ax2.set_ylabel('Latitude', fontsize='20')

gl2 = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, alpha=0.5, linestyle='-', color='None')
gl2.top_labels=False # suppress top labels
gl2.right_labels=False # suppress right labels
gl2.xlabel_style = {'size': 20}
gl2.ylabel_style = {'size': 20}

gl2.xlocator = mticker.FixedLocator([-118,-117,-116, -115, -114])
gl2.ylocator = mticker.FixedLocator([52.5,53,53.5,54,54.5])

ax2.plot(-115.88, 53.10, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax2.text(-115.925, 53.15, '49', fontsize=18, weight='bold', transform=ccrs.PlateCarree())
ax2.text(0.01, 0.975, 'c)', transform=ax2.transAxes, fontsize=24, va='top')

ax3 = fig.add_subplot(gs[1, 1], projection=ccrs.PlateCarree())
img2 = mpimg.imread(r'/home/clintonmacadam/nhpanalysis/'+date2+'/edmonton_radarimg.png')
ax3.imshow(
    img2,
    origin="upper",
    extent=[-118.5, -114, 52.5, 54.5],  # [lon_min, lon_max, lat_min, lat_max]
    transform=ccrs.PlateCarree(),
    interpolation="nearest",
    zorder=1
)
ax3.set_title('Observed Radar Reflectivity', fontsize=24)
ax3.set_aspect('auto')

ax3.plot(-115.88, 53.10, marker='o', markersize=20, color='b', transform=ccrs.PlateCarree())
ax3.text(-115.965, 53.18, '49', fontsize=18, color='white', weight='bold', transform=ccrs.PlateCarree())
ax3.text(0.01, 0.975, 'd)', transform=ax3.transAxes, color='white', fontsize=24, va='top')

gl3 = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, alpha=0.5, linestyle='-', color='None')
gl3.top_labels=False # suppress top labels
gl3.right_labels=False # suppress right labels
gl3.xlabel_style = {'size': 20}
gl3.ylabel_style = {'size': 20}
gl3.xlocator = mticker.FixedLocator([-118,-117,-116,-115])
gl3.ylocator = mticker.FixedLocator([52.5,53,53.5,54,54.5])

#Plot figure, save figure, and show figure
fig.tight_layout()
plt.savefig(''+date+'/stormlocationvalidationfigure_'+ date +'_'+ hour +'_'+ half_hour +'_test.png')
plt.show()
