import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pyart
from matplotlib.colors import ListedColormap, BoundaryNorm
from PIL import Image, ImageDraw

date='20230723'

cm1output1 = xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20230723_1250m_mp50_2.nc', engine='netcdf4')
cm1output2 = xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20230723_1250m_mp51_nolf_2.nc', engine='netcdf4')

for i in range(7,36):
    reflectivity1 = cm1output1['dbz']
    maxreflectivity1 = reflectivity1.max(dim="zh")
    maxreflectivity1 = maxreflectivity1[i,:,:]
    reflectivity2 = cm1output2['dbz']
    maxreflectivity2 = reflectivity2.max(dim="zh")
    maxreflectivity2 = maxreflectivity2[i,:,:]
    
    maximumreflectivity1=str(maxreflectivity1.max(dim=['xh','yh']).values)
    maximumreflectivity2=str(maxreflectivity2.max(dim=['xh','yh']).values)

    saveindex=str((i*5)+5)

    cmap=plt.get_cmap('RdYlGn_r')
    new_colors = cmap(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 1])  # RGBA for white

    # Mask where data == 0
    masked_data1 = np.ma.masked_where(maxreflectivity1 <= 0, maxreflectivity1)
    masked_data2 = np.ma.masked_where(maxreflectivity2 <= 0, maxreflectivity2)

    cmap = ListedColormap(new_colors, 256)
    cmap.set_bad(color='white')  # Will be used for masked values
    
    num_colors = 8
    jet= plt.cm.get_cmap('jet', num_colors)
    fig, axs = plt.subplots(1, 2, figsize=(17, 6.5))

    maxreflectivity1['xh']=maxreflectivity1['xh']-30
    maxreflectivity2['xh']=maxreflectivity2['xh']-30

    #Plotting code
    cf1 = axs[0].pcolormesh(maxreflectivity1['xh'], maxreflectivity1['yh'], masked_data1, cmap="NWSRef", vmin=0, vmax=70)
    axs[0].set_title('Column Maximum Reflectivity', fontsize=24)
    cbar1=plt.colorbar(cf1,ax=axs[0])
    cbar1.set_label('Reflectivity (dBZ)', fontsize=22, labelpad=8)
    cbar1.set_ticks([0,10,20,30,40,50,60,70], labels=['0','10','20','30','40','50','60','70'], fontsize=22)
    axs[0].set_xlabel('x [km]', fontsize=22)
    axs[0].set_ylabel('y [km]', fontsize=22)
    axs[0].set_xlim(-60,60)
    axs[0].set_ylim(-60,60)
    axs[0].set_xticks([-60,-30,0,30,60],labels=['-60','-30','0','30','60'],fontsize=22)
    axs[0].set_yticks([-60,-30,0,30,60],labels=['-60','-30','0','30','60'],fontsize=22)
    axs[0].text(0.02, 0.97, 'a)', transform=axs[0].transAxes, fontsize=22, va='top')
    cf2 = axs[1].pcolormesh(maxreflectivity2['xh'], maxreflectivity2['yh'], masked_data2, cmap="NWSRef", vmin=0, vmax=70)
    axs[1].set_title('Column Maximum Reflectivity', fontsize=24)
    cbar2=plt.colorbar(cf2,ax=axs[1])
    cbar2.set_label('Reflectivity (dBZ)', fontsize=22, labelpad=8)
    axs[1].set_xlim(-60,60)
    axs[1].set_ylim(-60,60)
    axs[1].set_xticks([-60,-30,0,30,60],labels=['-60','-30','0','30','60'],fontsize=22)
    axs[1].set_yticks([-60,-30,0,30,60],labels=['-60','-30','0','30','60'],fontsize=22)
    axs[1].set_xlabel('x [km]', fontsize=22)
    axs[1].set_ylabel('y [km]', fontsize=22)
    axs[1].text(0.02, 0.97, 'b)', transform=axs[1].transAxes, fontsize=22, va='top')
    cbar2.set_ticks([0,10,20,30,40,50,60,70], labels=['0','10','20','30','40','50','60','70'], fontsize=22)
    #plt.figlegend(['domain max='+maximumreflectivity1+''], bbox_to_anchor=(-0.205, 0.38, 0.5, 0.5), fontsize=14, handlelength=0, handletextpad=0)
    #plt.figlegend(['domain max='+maximumreflectivity2+''], bbox_to_anchor=(0.2175, 0.38, 0.5, 0.5), fontsize=14, handlelength=0, handletextpad=0)
    plt.figlegend([''+saveindex+ ' min'], bbox_to_anchor=(-0.34, -0.24, 0.5, 0.5), fontsize=24, handlelength=0, handletextpad=0)
    plt.figlegend([''+saveindex+ ' min'], bbox_to_anchor=(0.156, -0.24, 0.5, 0.5), fontsize=24, handlelength=0, handletextpad=0)
    plt.figlegend(['2MOM'], bbox_to_anchor=(-0.0935, 0.4425, 0.5, 0.5), fontsize=24, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM'], bbox_to_anchor=(0.401, 0.4423, 0.5, 0.5), fontsize=24, handlelength=0, handletextpad=0)
    #plt.figlegend(['a)'], bbox_to_anchor=(-0.3875, -0.24, 0.5, 0.5), fontsize=24, handlelength=0, handletextpad=0)
    #plt.figlegend(['b)'], bbox_to_anchor=(0.108, -0.24, 0.5, 0.5), fontsize=24, handlelength=0, handletextpad=0)
    plt.tight_layout()
    plt.savefig(''+ date +'/cm1columnmaxreflectivity_2panel_'+ saveindex +'_min.png')

image_list=[''+ date +'/cm1columnmaxreflectivity_2panel_40_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_45_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_50_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_55_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_60_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_65_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_70_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_75_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_80_min.png',  ''+ date +'/cm1columnmaxreflectivity_2panel_85_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_90_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_95_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_100_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_105_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_110_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_115_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_120_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_125_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_130_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_135_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_140_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_145_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_150_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_155_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_160_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_165_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_170_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_175_min.png', ''+ date +'/cm1columnmaxreflectivity_2panel_180_min.png']

images=[Image.open(image) for image in image_list]
images[0].save(''+ date +'/columnmaxreflectivityloop_2panel_'+ date +'_1250m_comparemp51nolf_mp51.gif', save_all=True, append_images=images[1:], duration=600, loop=0)
