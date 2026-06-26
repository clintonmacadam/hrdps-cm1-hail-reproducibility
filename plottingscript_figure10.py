import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.special import gamma
from PIL import Image, ImageDraw

date='20240805'
cm1output1=xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_2500m_mp50.nc', engine='netcdf4')
cm1output2=xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_2500m_mp51_nolf.nc', engine='netcdf4')

examplearray1=cm1output1['qi1']
examplearray1=examplearray1[:,0,:,:]
examplearray2=cm1output2['qi1']
examplearray2=examplearray2[:,0,:,:]
maximumHailSizegri1=xr.zeros_like(examplearray1)
maximumHailSizemaxgr1=xr.zeros_like(examplearray1)
maximumHailSizegrshifted1=xr.zeros_like(examplearray1)
maximumHailSizegri2=xr.zeros_like(examplearray2)
maximumHailSizemaxgr2=xr.zeros_like(examplearray2)
maximumHailSizegrshifted2=xr.zeros_like(examplearray2)

for i in range(7,36):
    icenumbermr1 = cm1output1['ni1']
    icenumbermr1 = icenumbermr1[i,0,:,:]
    icemassmr1 = cm1output1['qi1']
    icemassmr1 = icemassmr1[i,0,:,:]
    dmh1 = cm1output1['p3_dmi1']
    dmh1 = dmh1[i,0,:,:]
    airdensity1 = cm1output1['rho']
    airdensity1 = airdensity1[i,0,:,:]
    icedensity1 = cm1output1['p3_rhoi1']
    icedensity1 = icedensity1[i,0,:,:]
    rimedicemr1 = cm1output1['ri1']
    rimedicemr1 = rimedicemr1[i,0,:,:]
    rimedicefraction1=rimedicemr1/icemassmr1
    icenumbermr2 = cm1output2['ni1']
    icenumbermr2 = icenumbermr2[i,0,:,:]
    icemassmr2 = cm1output2['qi1']
    icemassmr2 = icemassmr2[i,0,:,:]
    dmh2 = cm1output2['p3_dmi1']
    dmh2 = dmh2[i,0,:,:]
    airdensity2 = cm1output2['rho']
    airdensity2 = airdensity2[i,0,:,:]
    #maxHailSize1 = cm1output1['p3_dhmax1']
    #maxHailSize1 = maxHailSize1[i,0,:,:]
    maxHailSize2 = cm1output2['p3_dhmax1']
    maxHailSize2 = maxHailSize2[i,0,:,:]
    accumulatedhail1=cm1output1['out2d1']
    accumulatedhail2=cm1output2['out2d1']
    accumulatedprecipa=cm1output1['rain2']
    accumulatedprecip1=accumulatedprecipa[i,:,:]
    accumulatedprecip1=accumulatedprecip1-accumulatedprecipa[7,:,:]
    accumulatedprecip1=accumulatedprecip1*10
    accumulatedprecipb=cm1output2['rain2']
    accumulatedprecip2=accumulatedprecipb[i,:,:]
    accumulatedprecip2=accumulatedprecip2-accumulatedprecipb[7,:,:]
    accumulatedprecip2=accumulatedprecip2*10
    saveindex=str((i*5)+5)

    def calc_accumulatedhailgr(accumulatedhail):
        umove=5.5
        vmove=4.4
        dt=300
        dxgrid=2500
        dygrid=2500
        accumulatedhailgr = xr.zeros_like(accumulatedhail)
        accumulatedhailgrshifted = xr.zeros_like(accumulatedhail)
        for i in range(8,36):
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
    
    def calc_max_hail_size(icenumbermr,icemassmr,airdensity,rimedicefraction,icedensity):
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
        mu=xr.zeros_like(icemassmr) #Set up DataArray for mu
        Dmh=((6*icemassmr)/(np.pi*icedensity*icenumbermr))**(1/dx) #Calculate mean-mass diameter of hail
        print(Dmh.values)
        icedensity=icedensity*rimedicefraction+(1-rimedicefraction)*917
        mu=0
        lam=((np.pi*icedensity*icenumbermr*gamma(mu+4))/((gamma(mu+1))*icemassmr))**(1/(mu+3))
        n0=(icenumbermr)*(lam)**(mu+1)/(gamma(mu+1)) #Calculate n0 using relation used in P3

        for j in range(nd,1,-1): #Loop to calculate max hail size
            Di  = j*dD
            V_h = rhofaci*(ch*Di**dh) #Calculate terminal fall speed
            R_tail = R_tail + V_h*n0*Di**mu*np.exp(-lam*Di)*dD
            maxHailSize = xr.where((R_tail > R_crit) & (Di > maxHailSize), Di, maxHailSize)
        return maxHailSize

    maxHailSize1=calc_max_hail_size(icenumbermr1,icemassmr1,airdensity1,rimedicefraction1,icedensity1)
    #maxHailSize2=calc_max_hail_size(icenumbermr2,icemassmr2,dmh2,airdensity2)
    
    accumulatedhailgr1=calc_accumulatedhailgr(accumulatedhail1)
    accumulatedhailgr2=calc_accumulatedhailgr(accumulatedhail2)
    accumulatedhailgr1=accumulatedhailgr1*10
    accumulatedhailgr2=accumulatedhailgr2*10

    accumulatedhailgr1instant=accumulatedhailgr1[i,:,:]
    accumulatedhailgr2instant=accumulatedhailgr2[i,:,:]

    accumulatedhailtotal1=accumulatedhailgr1instant.sum(dim=['xh','yh']).values 
    accumulatedhailtotal2=accumulatedhailgr2instant.sum(dim=['xh','yh']).values
    
    #print(accumulatedhailtotal1)
    #print(accumulatedhailtotal2)

    accumulatedhailmean1=accumulatedhailgr1instant.mean(dim=['xh','yh']).values
    accumulatedhailmean2=accumulatedhailgr2instant.mean(dim=['xh','yh']).values

    #print(accumulatedhailmean1*1000)
    #print(accumulatedhailmean2*1000)

    accumulatedpreciptotal1=accumulatedprecip1.sum(dim=['xh','yh']).values
    accumulatedpreciptotal2=accumulatedprecip2.sum(dim=['xh','yh']).values

    #print(accumulatedpreciptotal1)
    #print(accumulatedpreciptotal2)

    maxaccumulatedhailgr1=str(accumulatedhailgr1instant.max(dim=['xh','yh']).values)
    maxaccumulatedhailgr2=str(accumulatedhailgr2instant.max(dim=['xh','yh']).values)

    print(maxaccumulatedhailgr2)

    maxaccumulatedprecipgr1=str(accumulatedprecip1.max(dim=['xh','yh']).values)
    maxaccumulatedprecipgr2=str(accumulatedprecip2.max(dim=['xh','yh']).values)

    print(maxaccumulatedprecipgr2)

    accumulatedprecipmean1=accumulatedprecip1.mean(dim=['xh','yh']).values
    accumulatedprecipmean2=accumulatedprecip2.mean(dim=['xh','yh']).values

    #print(accumulatedprecipmean1*1000)
    #print(accumulatedprecipmean2*1000)

    umove=5.5
    vmove=4.4
    dt=300
    dxgrid1=2500
    dygrid1=2500
    dxgrid2=2500
    dygrid2=2500
    dx=umove*dt*(i-0.5)
    dy=vmove*dt*(i-0.5)
    shift_x1 = int(np.round(dx / dxgrid1))
    shift_y1 = int(np.round(dy / dygrid1))
    shift_x2 = int(np.round(dx / dxgrid2))
    shift_y2 = int(np.round(dy / dygrid2))
    maximumHailSizetimestep1=maxHailSize1
    maximumHailSizetimestepgr1=maximumHailSizetimestep1.roll(xh=shift_x1,yh=shift_y1, roll_coords=False)
    maximumHailSizegri1[i,:,:]=maximumHailSizetimestepgr1
    maximumHailSizemaxgr1[i,:,:]=maximumHailSizegri1.max(dim='time')
    maximumHailSizegrshifted1[i,:,:]=maximumHailSizemaxgr1[i,:,:].roll(xh=-shift_x1,yh=-shift_y1, roll_coords=False)
    maximumHailSizetimestep2=maxHailSize2
    maximumHailSizetimestepgr2=maximumHailSizetimestep2.roll(xh=shift_x2,yh=shift_y2, roll_coords=False)
    maximumHailSizegri2[i,:,:]=maximumHailSizetimestepgr2
    maximumHailSizemaxgr2[i,:,:]=maximumHailSizegri2.max(dim='time')
    maximumHailSizegrshifted2[i,:,:]=maximumHailSizemaxgr2[i,:,:].roll(xh=-shift_x2,yh=-shift_y2, roll_coords=False)

    maximumHailSizegr1=maximumHailSizegrshifted1[i,:,:]
    maximumHailSizegr1=maximumHailSizegr1*1000
    maximumHailSizegr2=maximumHailSizegrshifted2[i,:,:]
    maximumHailSizegr2=maximumHailSizegr2*1000
    
    maximumhailsize1=str(maximumHailSizegr1.max(dim=['xh','yh']).values)
    maximumhailsize2=str(maximumHailSizegr2.max(dim=['xh','yh']).values)

    #print(maximumhailsize1)
    #print(maximumhailsize2)

    maximumhailsizemean1=maximumHailSizegr1.mean(dim=['xh','yh']).values
    maximumhailsizemean2=maximumHailSizegr2.mean(dim=['xh','yh']).values

    #print(maximumhailsizemean1*1000)
    #print(maximumhailsizemean2*1000)

    cmap=plt.get_cmap('RdYlGn_r')
    new_colors = cmap(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 1])  # RGBA for white

    # Mask where data == 0
    masked_data1 = np.ma.masked_where(maximumHailSizegr1 <= 0, maximumHailSizegr1)
    masked_data2 = np.ma.masked_where(maximumHailSizegr2 <= 0, maximumHailSizegr2)
    masked_data3 = np.ma.masked_where(accumulatedhailgr1instant <= 0.001, accumulatedhailgr1instant)
    masked_data4 = np.ma.masked_where(accumulatedhailgr2instant <= 0.001, accumulatedhailgr2instant)
    masked_data5 = np.ma.masked_where(accumulatedprecip1 <= 1.0, accumulatedprecip1)
    masked_data6 = np.ma.masked_where(accumulatedprecip2 <= 1.0, accumulatedprecip2)

    cmap = ListedColormap(new_colors, 256)
    cmap.set_bad(color='white')  # Will be used for masked values
    hail = mcolors.LinearSegmentedColormap.from_list('hail', [(0, "blue"), (0.125, "blue"), (0.125, "cyan"), (0.25, "cyan"), (0.25, "green"), (0.375, "green"), (0.375, "lime"), (0.5, "lime"), (0.5, "yellow"), (0.625, "yellow"), (0.625, "orange"), (0.75, "orange"), (0.75, "red"), (0.875, "red"), (0.875, "purple"), (1.0, "purple")])
    hailaccum = mcolors.LinearSegmentedColormap.from_list('hail', [(0, "blue"), (0.025, "blue"), (0.025, "cyan"), (0.05, "cyan"), (0.05, "green"), (0.1, "green"), (0.1, "lime"), (0.2, "lime"), (0.2, "yellow"), (0.4, "yellow"), (0.4, "orange"), (0.6, "orange"), (0.6, "red"), (0.8, "red"), (0.8, "purple"), (1.0, "purple")])
    fig, axs = plt.subplots(3, 2, figsize=(17, 19.5))

    # Reflectivity
    cf1 = axs[0,0].pcolormesh(icemassmr1['xh'], icemassmr1['yh'], masked_data1, cmap=hail, vmin=0, vmax=24)
    axs[0,0].set_title('Surface Maximum Hail Size', fontsize=30)
    cbar1=plt.colorbar(cf1,ax=axs[0,0], shrink=0.98)
    axs[0,0].set_ylabel('y [km]', fontsize=28)
    axs[0,0].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=28)
    axs[0,0].set_yticks([-50,-25,0,25,50],labels=['-50','-25','0','25','50'],fontsize=28)
    axs[0,0].set_xlim(-75,75)
    axs[0,0].set_ylim(-75,75)
    cbar1.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=28)
    axs[0,0].text(0.02, 0.97, 'a)', transform=axs[0,0].transAxes, fontsize=28, va='top')
    cf2 = axs[0,1].pcolormesh(icemassmr2['xh'], icemassmr2['yh'], masked_data2, cmap=hail, vmin=0, vmax=24)
    axs[0,1].set_title('Surface Maximum Hail Size', fontsize=30)
    cbar2=plt.colorbar(cf2,ax=axs[0,1], shrink=0.98)
    axs[0,1].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=28)
    axs[0,1].set_yticks([-50,-25,0,25,50],labels=['-50','-25','0','25','50'],fontsize=28)
    axs[0,1].set_xlim(-75,75)
    axs[0,1].set_ylim(-75,75)
    axs[0,1].text(0.02, 0.97, 'b)', transform=axs[0,1].transAxes, fontsize=28, va='top')
    cbar2.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=28)
    cf3 = axs[1,0].pcolormesh(icemassmr1['xh'], icemassmr1['yh'], masked_data3, cmap=hailaccum, vmin=0, vmax=1)
    axs[1,0].set_title('Accumulated Hail', fontsize=30)
    cbar3=plt.colorbar(cf3,ax=axs[1,0], shrink=0.98)
    axs[1,0].set_ylabel('y [km]', fontsize=28)
    axs[1,0].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=28)
    axs[1,0].set_yticks([-50,-25,0,25,50],labels=['-50','-25','0','25','50'],fontsize=28)
    axs[1,0].set_xlim(-75,75)
    axs[1,0].set_ylim(-75,75)
    cbar3.set_ticks([0,0.1,0.2,0.4,0.6,0.8,1.0], labels=['0','0.1','0.2','0.4','0.6','0.8','1.0'], fontsize=28)
    axs[1,0].text(0.02, 0.97, 'c)', transform=axs[1,0].transAxes, fontsize=28, va='top')
    cf4 = axs[1,1].pcolormesh(icemassmr2['xh'], icemassmr2['yh'], masked_data4, cmap=hailaccum, vmin=0, vmax=1)
    axs[1,1].set_title('Accumulated Hail', fontsize=30)
    cbar4=plt.colorbar(cf4,ax=axs[1,1], shrink=0.98)
    axs[1,1].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=28)
    axs[1,1].set_yticks([-50,-25,0,25,50],labels=['-50','-25','0','25','50'],fontsize=28)
    axs[1,1].set_xlim(-75,75)
    axs[1,1].set_ylim(-75,75)
    cbar4.set_ticks([0,0.1,0.2,0.4,0.6,0.8,1.0], labels=['0','0.1','0.2','0.4','0.6','0.8','1.0'], fontsize=28)
    axs[1,1].text(0.02, 0.97, 'd)', transform=axs[1,1].transAxes, fontsize=28, va='top')
    cf5 = axs[2,0].pcolormesh(icemassmr1['xh'], icemassmr1['yh'], masked_data5, cmap=hail, vmin=0, vmax=24)
    axs[2,0].set_title('Accumulated Total Precipitation', fontsize=30)
    cbar5=plt.colorbar(cf5,ax=axs[2,0], shrink=0.98)
    axs[2,0].set_xlabel('x [km]', fontsize=28)
    axs[2,0].set_ylabel('y [km]', fontsize=28)
    axs[2,0].set_xlim(-75,75)
    axs[2,0].set_ylim(-75,75)
    axs[2,0].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=28)
    axs[2,0].set_yticks([-50,-25,0,25,50],labels=['-50','-25','0','25','50'],fontsize=28)
    cbar5.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=28)
    axs[2,0].text(0.02, 0.97, 'e)', transform=axs[2,0].transAxes, fontsize=28, va='top')
    cf6 = axs[2,1].pcolormesh(icemassmr2['xh'], icemassmr2['yh'], masked_data6, cmap=hail, vmin=0, vmax=24)
    axs[2,1].set_title('Accumulated Total Precipitation', fontsize=30)
    cbar6=plt.colorbar(cf6,ax=axs[2,1], shrink=0.98)
    axs[2,1].set_xlabel('x [km]', fontsize=28)
    axs[2,1].set_xticks([-75,-50,-25,0,25,50,75],labels=['-75','-50','-25','0','25','50','75'],fontsize=28)
    axs[2,1].set_yticks([-50,-25,0,25,50],labels=['-50','-25','0','25','50'],fontsize=28)
    axs[2,1].set_xlim(-75,75)
    axs[2,1].set_ylim(-75,75)
    cbar6.set_ticks([0,3,6,9,12,15,18,21,24], labels=['0','3','6','9','12','15','18','21','24'], fontsize=28)
    axs[2,1].text(0.02, 0.97, 'f)', transform=axs[2,1].transAxes, fontsize=28, va='top')
    #plt.figlegend([''+date+''], bbox_to_anchor=(0.5, 0.505, 0.5, 0.5), fontsize=20, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=18'], bbox_to_anchor=(-0.29, 0.2495, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=13'], bbox_to_anchor=(0.1935, 0.2495, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=0.9'], bbox_to_anchor=(-0.2835, -0.073, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=0.1'], bbox_to_anchor=(0.2, -0.073, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=28'], bbox_to_anchor=(-0.29, -0.397, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['MAX=24'], bbox_to_anchor=(0.1935, -0.397, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    #plt.figlegend(['TOTAL=196'], bbox_to_anchor=(-0.07, -0.073, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    #plt.figlegend(['TOTAL=93'], bbox_to_anchor=(0.4125, -0.073, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['2MOM'], bbox_to_anchor=(-0.0725, 0.475, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM'], bbox_to_anchor=(0.41, 0.475, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['2MOM'], bbox_to_anchor=(-0.0725, 0.153, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM'], bbox_to_anchor=(0.41, 0.153, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['2MOM'], bbox_to_anchor=(-0.0725, -0.17, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    plt.figlegend(['3MOM'], bbox_to_anchor=(0.41, -0.17, 0.5, 0.5), fontsize=28, handlelength=0, handletextpad=0)
    fig.tight_layout()
    plt.savefig(''+ date +'/cm1hailswathfigureoutput_2panel_'+ saveindex +'_min.png')

image_list=[''+ date +'/cm1hailswathfigureoutput_2panel_40_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_45_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_50_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_55_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_60_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_65_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_70_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_75_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_80_min.png',  ''+ date +'/cm1hailswathfigureoutput_2panel_85_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_90_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_95_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_100_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_105_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_110_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_115_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_120_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_125_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_130_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_135_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_140_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_145_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_150_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_155_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_160_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_165_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_170_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_175_min.png', ''+ date +'/cm1hailswathfigureoutput_2panel_180_min.png'] 

images=[Image.open(image) for image in image_list]
images[0].save(''+ date +'/hailswathloop_groundrelative_2panel_2500m_comparemp50mp51nolf_'+date+'.gif', save_all=True, append_images=images[1:], duration=600, loop=0)
