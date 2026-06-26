import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.special import gamma
from PIL import Image, ImageDraw

date='20240805'
cm1output1=xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_1250m_mp51_nolf.nc', engine='netcdf4')
cm1output2=xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_1250m_mp51.nc', engine='netcdf4')
cm1output3=xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_1250m_mp52.nc', engine='netcdf4')
cm1output4=xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_1250m_mp53.nc', engine='netcdf4')
cm1output5=xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_1250m_mp54.nc', engine='netcdf4')

examplearray1=cm1output1['qi1']
examplearray1=examplearray1[:,0,:,:]
examplearray2=cm1output2['qi1']
examplearray2=examplearray2[:,0,:,:]
examplearray3=cm1output3['qi1']
examplearray3=examplearray3[:,0,:,:]
examplearray4=cm1output4['qi1']
examplearray4=examplearray4[:,0,:,:]
examplearray5=cm1output5['qi1']
examplearray5=examplearray5[:,0,:,:]

maximumHailSizegri1=xr.zeros_like(examplearray1)
maximumHailSizemaxgr1=xr.zeros_like(examplearray1)
maximumHailSizegrshifted1=xr.zeros_like(examplearray1)
maximumHailSizegri2=xr.zeros_like(examplearray2)
maximumHailSizemaxgr2=xr.zeros_like(examplearray2)
maximumHailSizegrshifted2=xr.zeros_like(examplearray2)
maximumHailSizegri3=xr.zeros_like(examplearray3)
maximumHailSizemaxgr3=xr.zeros_like(examplearray3)
maximumHailSizegrshifted3=xr.zeros_like(examplearray3)
maximumHailSizegri4=xr.zeros_like(examplearray4)
maximumHailSizemaxgr4=xr.zeros_like(examplearray4)
maximumHailSizegrshifted4=xr.zeros_like(examplearray4)
maximumHailSizegri5=xr.zeros_like(examplearray5)
maximumHailSizemaxgr5=xr.zeros_like(examplearray5)
maximumHailSizegrshifted5=xr.zeros_like(examplearray5)

for i in range(7,36):
    maxHailSize1 = cm1output1['p3_dhmax1']
    maxHailSize1 = maxHailSize1[i,0,:,:]
    maxHailSize2 = cm1output2['p3_dhmax1']
    maxHailSize2 = maxHailSize2[i,0,:,:]

    maxHailSize3 = cm1output3[['p3_dhmax1','p3_dhmax2']]
    maxHailSize3 = maxHailSize3.isel(time=i,zh=0)
    maxHailSize3 = maxHailSize3.to_array(dim="category").max(dim="category")

    maxHailSize4 = cm1output4[['p3_dhmax1','p3_dhmax2','p3_dhmax3']]
    maxHailSize4 = maxHailSize4.isel(time=i,zh=0)
    maxHailSize4 = maxHailSize4.to_array(dim="category").max(dim="category")

    maxHailSize5 = cm1output5[['p3_dhmax1','p3_dhmax2','p3_dhmax3','p3_dhmax4']]
    maxHailSize5 = maxHailSize5.isel(time=i,zh=0)
    maxHailSize5 = maxHailSize5.to_array(dim="category").max(dim="category")

    saveindex=str((i*5)+5)

    umove=5.5
    vmove=4.4
    dt=300
    dxgrid=1250
    dygrid=1250
    dx=umove*dt*(i-0.5)
    dy=vmove*dt*(i-0.5)
    shift_x = int(np.round(dx / dxgrid))
    shift_y = int(np.round(dy / dygrid))
    maximumHailSizetimestep1=maxHailSize1
    maximumHailSizetimestepgr1=maximumHailSizetimestep1.roll(xh=shift_x,yh=shift_y, roll_coords=False)
    maximumHailSizegri1[i,:,:]=maximumHailSizetimestepgr1
    maximumHailSizemaxgr1[i,:,:]=maximumHailSizegri1.max(dim='time')
    maximumHailSizegrshifted1[i,:,:]=maximumHailSizemaxgr1[i,:,:].roll(xh=-shift_x,yh=-shift_y, roll_coords=False)
    maximumHailSizetimestep2=maxHailSize2
    maximumHailSizetimestepgr2=maximumHailSizetimestep2.roll(xh=shift_x,yh=shift_y, roll_coords=False)
    maximumHailSizegri2[i,:,:]=maximumHailSizetimestepgr2
    maximumHailSizemaxgr2[i,:,:]=maximumHailSizegri2.max(dim='time')
    maximumHailSizegrshifted2[i,:,:]=maximumHailSizemaxgr2[i,:,:].roll(xh=-shift_x,yh=-shift_y, roll_coords=False)
    maximumHailSizetimestep3=maxHailSize3
    maximumHailSizetimestepgr3=maximumHailSizetimestep3.roll(xh=shift_x,yh=shift_y, roll_coords=False)
    maximumHailSizegri3[i,:,:]=maximumHailSizetimestepgr3
    maximumHailSizemaxgr3[i,:,:]=maximumHailSizegri3.max(dim='time')
    maximumHailSizegrshifted3[i,:,:]=maximumHailSizemaxgr3[i,:,:].roll(xh=-shift_x,yh=-shift_y, roll_coords=False)
    maximumHailSizetimestep4=maxHailSize4
    maximumHailSizetimestepgr4=maximumHailSizetimestep4.roll(xh=shift_x,yh=shift_y, roll_coords=False)
    maximumHailSizegri4[i,:,:]=maximumHailSizetimestepgr4
    maximumHailSizemaxgr4[i,:,:]=maximumHailSizegri4.max(dim='time')
    maximumHailSizegrshifted4[i,:,:]=maximumHailSizemaxgr4[i,:,:].roll(xh=-shift_x,yh=-shift_y, roll_coords=False)
    maximumHailSizetimestep5=maxHailSize5
    maximumHailSizetimestepgr5=maximumHailSizetimestep5.roll(xh=shift_x,yh=shift_y, roll_coords=False)
    maximumHailSizegri5[i,:,:]=maximumHailSizetimestepgr5
    maximumHailSizemaxgr5[i,:,:]=maximumHailSizegri5.max(dim='time')
    maximumHailSizegrshifted5[i,:,:]=maximumHailSizemaxgr5[i,:,:].roll(xh=-shift_x,yh=-shift_y, roll_coords=False)

    maximumHailSizegr1=maximumHailSizegrshifted1[i,:,:]
    maximumHailSizegr1=maximumHailSizegr1*1000
    maximumHailSizegr2=maximumHailSizegrshifted2[i,:,:]
    maximumHailSizegr2=maximumHailSizegr2*1000
    maximumHailSizegr3=maximumHailSizegrshifted3[i,:,:]
    maximumHailSizegr3=maximumHailSizegr3*1000
    maximumHailSizegr4=maximumHailSizegrshifted4[i,:,:]
    maximumHailSizegr4=maximumHailSizegr4*1000
    maximumHailSizegr5=maximumHailSizegrshifted5[i,:,:]
    maximumHailSizegr5=maximumHailSizegr5*1000
    
    maximumhailsize1=str(maximumHailSizegr1.max(dim=['xh','yh']).values)
    maximumhailsize2=str(maximumHailSizegr2.max(dim=['xh','yh']).values)
    maximumhailsize3=str(maximumHailSizegr3.max(dim=['xh','yh']).values)
    maximumhailsize4=str(maximumHailSizegr4.max(dim=['xh','yh']).values)
    maximumhailsize5=str(maximumHailSizegr5.max(dim=['xh','yh']).values)
    #print(maximumhailsize1)
    #print(maximumhailsize2)

    accumulatedhail1=cm1output1['out2d1']
    accumulatedhail2=cm1output2['out2d1']
    accumulatedhail3=cm1output3['out2d1']
    accumulatedhail4=cm1output4['out2d1']
    accumulatedhail5=cm1output5['out2d1']
    dxgrid1=1250
    dxgrid2=1250
    dxgrid3=1250
    dxgrid4=1250
    dxgrid5=1250

    def calc_accumulatedhailgr(accumulatedhail,dxgrid):
        dygrid=dxgrid
        umove=5.5
        vmove=4.4
        dt=300
        accumulatedhailgr = xr.zeros_like(accumulatedhail)
        accumulatedhailgrshifted = xr.zeros_like(accumulatedhail)
        for i in range(8, 36):
            dx=umove*dt*(i-0.5)
            dy=vmove*dt*(i-0.5)
            shift_x = int(np.round(dx / dxgrid))
            shift_y = int(np.round(dy / dygrid))
            accumulatedhailtimestep=accumulatedhail.isel(time=i)-accumulatedhail.isel(time=i-1)
            accumulatedhailtimestepgr=accumulatedhailtimestep.roll(xh=shift_x,yh=shift_y, roll_coords=False)
            if i==0:
                accumulatedhailgr[i,:,:]=0
            else:
                accumulatedhailgr[i,:,:]=accumulatedhailgr[i-1,:,:]+accumulatedhailtimestepgr
                accumulatedhailgrshifted[i,:,:]=accumulatedhailgr[i,:,:].roll(xh=-shift_x,yh=-shift_y, roll_coords=False)
        return accumulatedhailgrshifted

    accumulatedhailgr1=calc_accumulatedhailgr(accumulatedhail1,dxgrid1)
    accumulatedhailgr2=calc_accumulatedhailgr(accumulatedhail2,dxgrid2)
    accumulatedhailgr3=calc_accumulatedhailgr(accumulatedhail3,dxgrid3)
    accumulatedhailgr4=calc_accumulatedhailgr(accumulatedhail4,dxgrid4)
    accumulatedhailgr5=calc_accumulatedhailgr(accumulatedhail5,dxgrid5)

    accumulatedhailgr1=accumulatedhailgr1*10
    accumulatedhailgr2=accumulatedhailgr2*10
    accumulatedhailgr3=accumulatedhailgr3*10
    accumulatedhailgr4=accumulatedhailgr4*10
    accumulatedhailgr5=accumulatedhailgr5*10

    accumulatedhailgr1instant=accumulatedhailgr1[i,:,:]
    accumulatedhailgr2instant=accumulatedhailgr2[i,:,:]
    accumulatedhailgr3instant=accumulatedhailgr3[i,:,:]
    accumulatedhailgr4instant=accumulatedhailgr4[i,:,:]
    accumulatedhailgr5instant=accumulatedhailgr5[i,:,:]

    #accumulatedhailgr1=accumulatedhailgr1.max(dim=['xh','yh']).values
    #accumulatedhailgr2=accumulatedhailgr2.max(dim=['xh','yh']).values
    #accumulatedhailgr3=accumulatedhailgr3.max(dim=['xh','yh']).values
    #accumulatedhailgr4=accumulatedhailgr4.max(dim=['xh','yh']).values
    #accumulatedhailgr5=accumulatedhailgr5.max(dim=['xh','yh']).values

    cmap=plt.get_cmap('RdYlGn_r')
    new_colors = cmap(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 1])  # RGBA for white

    # Mask where data == 0
    masked_data1 = np.ma.masked_where(maximumHailSizegr1 <= 0, maximumHailSizegr1)
    masked_data2 = np.ma.masked_where(maximumHailSizegr2 <= 0, maximumHailSizegr2)
    masked_data3 = np.ma.masked_where(maximumHailSizegr3 <= 0, maximumHailSizegr3)
    masked_data4 = np.ma.masked_where(maximumHailSizegr4 <= 0, maximumHailSizegr4)
    masked_data5 = np.ma.masked_where(maximumHailSizegr5 <= 0, maximumHailSizegr5)
    masked_data6 = np.ma.masked_where(accumulatedhailgr1instant <= 0, accumulatedhailgr1instant)
    masked_data7 = np.ma.masked_where(accumulatedhailgr2instant <= 0, accumulatedhailgr2instant)
    masked_data8 = np.ma.masked_where(accumulatedhailgr3instant <= 0, accumulatedhailgr3instant)
    masked_data9 = np.ma.masked_where(accumulatedhailgr4instant <= 0, accumulatedhailgr4instant)
    masked_data10 = np.ma.masked_where(accumulatedhailgr5instant <= 0, accumulatedhailgr5instant)

    cmap = ListedColormap(new_colors, 256)
    cmap.set_bad(color='white')  # Will be used for masked values
    hail = mcolors.LinearSegmentedColormap.from_list('hail', [(0, "blue"), (0.125, "blue"), (0.125, "cyan"), (0.25, "cyan"), (0.25, "green"), (0.375, "green"), (0.375, "lime"), (0.5, "lime"), (0.5, "yellow"), (0.625, "yellow"), (0.625, "orange"), (0.75, "orange"), (0.75, "red"), (0.875, "red"), (0.875, "purple"), (1.0, "purple")])
    hailaccum = mcolors.LinearSegmentedColormap.from_list('hail', [(0, "blue"), (0.05, "blue"), (0.05, "cyan"), (0.1, "cyan"), (0.1, "green"), (0.15, "green"), (0.15, "lime"), (0.2, "lime"), (0.2, "yellow"), (0.4, "yellow"), (0.4, "orange"), (0.6, "orange"), (0.6, "red"), (0.8, "red"), (0.8, "purple"), (1.0, "purple")])
    fig, axs = plt.subplots(5, 2, figsize=(14, 20))

    # Reflectivity
    cf1 = axs[0,0].pcolormesh(maxHailSize1['xh'], maxHailSize1['yh'], masked_data1, cmap=hail, vmin=0, vmax=24)
    axs[0,0].set_title('Surface Maximum Hail Size', fontsize=24)
    cbar1=plt.colorbar(cf1,ax=axs[0,0], shrink=0.95)
    axs[0,0].set_ylabel('y [km]', fontsize=20)
    axs[0,0].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[0,0].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[0,0].set_xlim(-90,90)
    axs[0,0].set_ylim(-90,90)
    cbar1.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=20)
    axs[0,0].text(0.02, 0.97, 'a)', transform=axs[0,0].transAxes, fontsize=20, va='top')
    cf2 = axs[0,1].pcolormesh(accumulatedhail1['xh'], accumulatedhail1['yh'], masked_data6, cmap=hailaccum, vmin=0, vmax=1)
    axs[0,1].set_title('Accumulated Hail', fontsize=24)
    cbar2=plt.colorbar(cf2,ax=axs[0,1], shrink=0.95)
    axs[0,1].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[0,1].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[0,1].set_xlim(-90,90)
    axs[0,1].set_ylim(-90,90)
    cbar2.set_ticks([0,0.1,0.2,0.4,0.6,0.8,1.0], labels=['0','0.1','0.2','0.4','0.6','0.8','1.0'], fontsize=20)
    axs[0,1].text(0.02, 0.97, 'b)', transform=axs[0,1].transAxes, fontsize=20, va='top')
    cf3 = axs[1,0].pcolormesh(maxHailSize2['xh'], maxHailSize2['yh'], masked_data2, cmap=hail, vmin=0, vmax=24)
    axs[1,0].set_title('Surface Maximum Hail Size', fontsize=24)
    cbar3=plt.colorbar(cf3,ax=axs[1,0], shrink=0.95)
    axs[1,0].set_ylabel('y [km]', fontsize=20)
    axs[1,0].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[1,0].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    cbar3.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=20)
    axs[1,0].set_xlim(-90,90)
    axs[1,0].set_ylim(-90,90)
    axs[1,0].text(0.02, 0.97, 'c)', transform=axs[1,0].transAxes, fontsize=20, va='top')
    cf4 = axs[1,1].pcolormesh(accumulatedhail2['xh'], accumulatedhail2['yh'], masked_data7, cmap=hailaccum, vmin=0, vmax=1)
    axs[1,1].set_title('Accumulated Hail', fontsize=24)
    cbar4=plt.colorbar(cf4,ax=axs[1,1], shrink=0.95)
    axs[1,1].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[1,1].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    cbar4.set_ticks([0,0.1,0.2,0.4,0.6,0.8,1.0], labels=['0','0.1','0.2','0.4','0.6','0.8','1.0'], fontsize=20)
    axs[1,1].set_xlim(-90,90)
    axs[1,1].set_ylim(-90,90)
    axs[1,1].text(0.02, 0.97, 'd)', transform=axs[1,1].transAxes, fontsize=20, va='top')
    cf5 = axs[2,0].pcolormesh(maxHailSize3['xh'], maxHailSize3['yh'], masked_data3, cmap=hail, vmin=0, vmax=24)
    axs[2,0].set_title('Surface Maximum Hail Size', fontsize=24)
    cbar5=plt.colorbar(cf5,ax=axs[2,0], shrink=0.95)
    axs[2,0].set_ylabel('y [km]', fontsize=20)
    axs[2,0].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[2,0].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    cbar5.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=20)
    axs[2,0].set_xlim(-90,90)
    axs[2,0].set_ylim(-90,90)
    axs[2,0].text(0.02, 0.97, 'e)', transform=axs[2,0].transAxes, fontsize=20, va='top')
    cf6 = axs[2,1].pcolormesh(accumulatedhail3['xh'], accumulatedhail3['yh'], masked_data8, cmap=hailaccum, vmin=0, vmax=1)
    axs[2,1].set_title('Accumulated Hail', fontsize=24)
    cbar6=plt.colorbar(cf6,ax=axs[2,1], shrink=0.95)
    axs[2,1].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[2,1].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    cbar6.set_ticks([0,0.1,0.2,0.4,0.6,0.8,1.0], labels=['0','0.1','0.2','0.4','0.6','0.8','1.0'], fontsize=20)
    axs[2,1].set_xlim(-90,90)
    axs[2,1].set_ylim(-90,90)
    axs[2,1].text(0.02, 0.97, 'f)', transform=axs[2,1].transAxes, fontsize=20, va='top')
    cf7 = axs[3,0].pcolormesh(maxHailSize4['xh'], maxHailSize4['yh'], masked_data4, cmap=hail, vmin=0, vmax=24)
    axs[3,0].set_title('Surface Maximum Hail Size', fontsize=24)
    cbar7=plt.colorbar(cf7,ax=axs[3,0], shrink=0.95)
    axs[3,0].set_ylabel('y [km]', fontsize=20)
    axs[3,0].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[3,0].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    cbar7.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=20)
    axs[3,0].set_xlim(-90,90)
    axs[3,0].set_ylim(-90,90)
    axs[3,0].text(0.02, 0.97, 'g)', transform=axs[3,0].transAxes, fontsize=20, va='top')
    cf8 = axs[3,1].pcolormesh(accumulatedhail4['xh'], accumulatedhail4['yh'], masked_data9, cmap=hailaccum, vmin=0, vmax=1)
    axs[3,1].set_title('Accumulated Hail', fontsize=24)
    cbar8=plt.colorbar(cf8,ax=axs[3,1], shrink=0.95)
    axs[3,1].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[3,1].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    cbar8.set_ticks([0,0.1,0.2,0.4,0.6,0.8,1.0], labels=['0','0.1','0.2','0.4','0.6','0.8','1.0'], fontsize=20)
    axs[3,1].set_xlim(-90,90)
    axs[3,1].set_ylim(-90,90)
    axs[3,1].text(0.02, 0.97, 'h)', transform=axs[3,1].transAxes, fontsize=20, va='top')
    cf9 = axs[4,0].pcolormesh(maxHailSize5['xh'], maxHailSize5['yh'], masked_data5, cmap=hail, vmin=0, vmax=24)
    axs[4,0].set_title('Surface Maximum Hail Size', fontsize=24)
    cbar9=plt.colorbar(cf9,ax=axs[4,0], shrink=0.95)
    axs[4,0].set_xlabel('x [km]', fontsize=20)
    axs[4,0].set_ylabel('y [km]', fontsize=20)
    axs[4,0].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[4,0].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    cbar9.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=20)
    axs[4,0].set_xlim(-90,90)
    axs[4,0].set_ylim(-90,90)
    axs[4,0].text(0.02, 0.97, 'i)', transform=axs[4,0].transAxes, fontsize=20, va='top')
    cf10 = axs[4,1].pcolormesh(accumulatedhail5['xh'], accumulatedhail5['yh'], masked_data10, cmap=hailaccum, vmin=0, vmax=1)
    axs[4,1].set_title('Accumulated Hail', fontsize=24)
    cbar10=plt.colorbar(cf10,ax=axs[4,1], shrink=0.95)
    axs[4,1].set_xlabel('x [km]', fontsize=20)
    axs[4,1].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    axs[4,1].set_yticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=20)
    cbar10.set_ticks([0,0.1,0.2,0.4,0.6,0.8,1.0], labels=['0','0.1','0.2','0.4','0.6','0.8','1.0'], fontsize=20)
    axs[4,1].set_xlim(-90,90)
    axs[4,1].set_ylim(-90,90)
    axs[4,1].text(0.02, 0.97, 'j)', transform=axs[4,1].transAxes, fontsize=20, va='top')

    #plt.figlegend([''+date+''], bbox_to_anchor=(0.5, 0.505, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=21'], bbox_to_anchor=(-0.31, 0.3575, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=6.0'], bbox_to_anchor=(0.1725, 0.3575, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=0'], bbox_to_anchor=(-0.3225, 0.1625, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=1.7'], bbox_to_anchor=(0.1725, 0.1625, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=23'], bbox_to_anchor=(-0.31, -0.0325, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=0.1'], bbox_to_anchor=(0.1725, -0.0325, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=19'], bbox_to_anchor=(-0.31, -0.2275, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=0.2'], bbox_to_anchor=(0.1725, -0.2275, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=20'], bbox_to_anchor=(-0.31, -0.4225, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=0.2'], bbox_to_anchor=(0.1725, -0.4225, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM'], bbox_to_anchor=(-0.07, 0.48, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM'], bbox_to_anchor=(0.4075, 0.48, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM_LF'], bbox_to_anchor=(-0.07, 0.283, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM_LF'], bbox_to_anchor=(0.4075, 0.283, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM_2CAT_LF'], bbox_to_anchor=(-0.07, 0.091, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM_2CAT_LF'], bbox_to_anchor=(0.4075, 0.091, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM_3CAT_LF'], bbox_to_anchor=(-0.07, -0.104, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM_3CAT_LF'], bbox_to_anchor=(0.4075, -0.104, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM_4CAT_LF'], bbox_to_anchor=(-0.07, -0.299, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM_4CAT_LF'], bbox_to_anchor=(0.4075, -0.299, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    fig.tight_layout()
    plt.savefig(''+ date +'/cm1hailswathfigureoutput_2panel_'+ saveindex +'_min.png')

image_list=[''+ date +'/cm1hailswathfigureoutput_2panel_40_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_45_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_50_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_55_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_60_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_65_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_70_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_75_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_80_min.png',  ''+ date +'/cm1hailswathfigureoutput_2panel_85_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_90_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_95_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_100_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_105_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_110_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_115_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_120_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_125_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_130_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_135_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_140_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_145_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_150_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_155_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_160_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_165_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_170_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_175_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_180_min.png'] 

images=[Image.open(image) for image in image_list]
images[0].save(''+ date +'/hailswathloop_groundrelative_8panel_1250m_comparemp51mp52mp53mp54_'+date+'.gif', save_all=True, append_images=images[1:], duration=600, loop=0)
