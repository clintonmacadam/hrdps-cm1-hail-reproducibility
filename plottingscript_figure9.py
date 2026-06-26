import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.special import gamma
from PIL import Image, ImageDraw

date='20240805'

cm1output1 = xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_2500m_mphrdps.nc', engine='netcdf4')
cm1output2 = xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_2500m_mp50.nc', engine='netcdf4')
#cm1output3 = xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_1250m_mp53.nc', engine='netcdf4')
#cm1output4 = xr.open_dataset(r'/home/cmacadam/projects/def-hanesiak/cmacadam/cm1r21.1/run/cm1out_20240805_1250m_mp54.nc', engine='netcdf4')

def calc_max_hail_size(cm1output):
    maximumhailsize=np.ones(36)
    for i in range(36):
        icenumbermr = cm1output['ni1']
        icenumbermr = icenumbermr[i,:,:,:]
        icemassmr = cm1output['qi1']
        icemassmr = icemassmr[i,:,:,:]
        icedensity = cm1output['p3_rhoi1']
        icedensity = icedensity[i,:,:,:]
        airdensity = cm1output['rho']
        airdensity = airdensity[i,:,:,:]
        rimedicemassmr = cm1output['ri1']
        rimedicemassmr = rimedicemassmr[i,:, :, :]

        rimedicefraction = rimedicemassmr/icemassmr
        saveindex=str(i*5)
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
        icedensity=icedensity*rimedicefraction+(1-rimedicefraction)*917
        mu=0
        lam=((np.pi*icedensity*icenumbermr*gamma(mu+4))/((gamma(mu+1))*icemassmr))**(1/(mu+3))
        n0=(icenumbermr)*(lam)**(mu+1)/(gamma(mu+1)) #Calculate n0 using relation used in P3

        for j in range(nd,1,-1): #Loop to calculate max hail size
            Di  = j*dD
            V_h = rhofaci*(ch*Di**dh) #Calculate terminal fall speed
            R_tail = R_tail + V_h*n0*Di**mu*np.exp(-lam*Di)*dD
            maxHailSize = xr.where((R_tail > R_crit) & (Di > maxHailSize), Di, maxHailSize)        

        maxHailSize=maxHailSize*1000 #Convert max hail size to mm
        print(maxHailSize)
        maximumhailsize[i]=maxHailSize.max(dim=['xh','yh','zh']).values
    return maximumhailsize

maximumhailsize1=calc_max_hail_size(cm1output1)
#maximumhailsize2=calc_max_hail_size(cm1output2)
#maximumhailsize3=calc_max_hail_size(cm1output3)
#maximumhailsize4=calc_max_hail_size(cm1output4)

maxHailSize2 = cm1output2['p3_dhmax1']
print(maxHailSize2)
maximumhailsize2=np.zeros(36)
maxHailSize2 = maxHailSize2*1000
for i in range(36):
    maximumhailsize2[i]=maxHailSize2[i,:,:].max(dim=['xh','yh','zh']).values

def calc_maxmeanmassicediameter(cm1output):
    maximumdmh=np.zeros(36)

    for i in range(36):
        qx = cm1output['qi1']
        qx = qx[i,0,:,:]
        nit = cm1output['ni1']
        nit = nit[i,0,:,:]
        icedensity = cm1output['p3_rhoi1']
        icedensity = icedensity[i,0,:,:]
        dmh=((6*qx)/(np.pi*icedensity*nit))**(1/3) #Calculate mean-mass diameter of hail
        dmh = dmh*1000
        maximumdmh[i]=dmh.max(dim=['xh','yh']).values
    return maximumdmh

maxmeanmassicediameter1=calc_maxmeanmassicediameter(cm1output1)
maxmeanmassicediameter2=calc_maxmeanmassicediameter(cm1output2)

#maxHailSize1 = cm1output1['p3_dhmax1']
#print(maxHailSize1)
#columnmaximumhailsize1=np.zeros(36)
#maxHailSize1 = maxHailSize1*1000

#for i in range(36):
#    columnmaximumhailsize1[i]=maxHailSize1[i,:,:,:].max(dim=['xh','yh','zh']).values

#maxHailSize2 = cm1output2[['p3_dhmax1','p3_dhmax2']]
#print(maxHailSize2)
#maxHailSize3 = cm1output3[['p3_dhmax1','p3_dhmax2','p3_dhmax3']]
#maxHailSize4 = cm1output4[['p3_dhmax1','p3_dhmax2','p3_dhmax3','p3_dhmax4']]

def calc_max_hail_size_2(maxHailSize):
    maximumhailsize=np.zeros(36)

    for i in range(36):
        maxhailsize = maxHailSize.isel(time=i)
        maxhailsize = maxhailsize.to_array(dim="category").max(dim="category")
        maxhailsize = maxhailsize*1000
        maximumhailsize[i]=maxhailsize.max(dim=['xh','yh','zh']).values
    return maximumhailsize

#columnmaximumhailsize2=calc_max_hail_size_2(maxHailSize2)
#columnmaximumhailsize3=calc_max_hail_size_2(maxHailSize3)
#columnmaximumhailsize4=calc_max_hail_size_2(maxHailSize4)

def calc_domainmaxreflectivity(cm1output):
    maxreflectivity=np.zeros(36)

    for i in range(36):
        reflectivity = cm1output['dbz']
        maximumreflectivity = reflectivity.max(dim="zh")
        maximumreflectivity = maximumreflectivity[i,:,:]
        maxreflectivity[i]=maximumreflectivity.max(dim=['xh','yh']).values
    return maxreflectivity

#maxreflectivity1=calc_domainmaxreflectivity(cm1output1)
#maxreflectivity2=calc_domainmaxreflectivity(cm1output2)
#maxreflectivity3=calc_domainmaxreflectivity(cm1output3)
#maxreflectivity4=calc_domainmaxreflectivity(cm1output4)

time=np.arange(5,185,5)

fig, axs = plt.subplots(1, 2, figsize=(17, 6.5))

axs[0].plot(time,maximumhailsize1, color='red',linewidth=4)
axs[0].plot(time,maximumhailsize2, color='blue',linewidth=4)
#axs[0].plot(time,maximumhailsize3)
#axs[0].plot(time,maximumhailsize4)
axs[0].legend(['v3 P3','v5 P3: 2MOM'], loc='lower right', fontsize=21)
axs[0].set_title('Column Max Hail Size', fontsize=29)
axs[0].set_xlabel('Simulation time (mins)', fontsize=26, labelpad=12)
axs[0].set_ylabel('Maximum Hail Size (mm)', fontsize=26, labelpad=12)
axs[0].set_xticks([40,60,80,100,120,140,160,180],labels=['40','60','80','100','120','140','160','180'],fontsize=26)
axs[0].set_xlim(40,180)
axs[0].set_yticks([5,10,15,20,25],labels=['5','10','15','20','25'],fontsize=26)
axs[0].set_ylim(0,30)
axs[0].text(0.02, 0.97, 'a)', transform=axs[0].transAxes, fontsize=26, va='top')

axs[1].plot(time,maxmeanmassicediameter1, color='red',linewidth=4)
axs[1].plot(time,maxmeanmassicediameter2, color='blue',linewidth=4)
#axs[1].plot(time,maxreflectivity3)
#axs[1].plot(time,maxreflectivity4)
axs[1].legend(['v3 P3','v5 P3: 2MOM'], loc='lower right', fontsize=21)
axs[1].set_title('Surface Max Mean-Mass Diameter', fontsize=29)
axs[1].set_xlabel('Simulation time (mins)', fontsize=26, labelpad=12)
axs[1].set_ylabel('Mean-Mass Diameter (mm)', fontsize=26, labelpad=12)
axs[1].set_xticks([40,60,80,100,120,140,160,180],labels=['40','60','80','100','120','140','160','180'],fontsize=26)
axs[1].set_xlim(40,180)
axs[1].set_yticks([2,4,6,8],labels=['2','4','6','8'],fontsize=26)
axs[1].set_ylim(0,10)
axs[1].text(0.02, 0.97, 'b)', transform=axs[1].transAxes, fontsize=26, va='top')

#axs[0].plot(time,columnmaximumhailsize1, color='red',linewidth=4)
#axs[0].plot(time,columnmaximumhailsize2, color='orange',linewidth=4)
#axs[0].plot(time,columnmaximumhailsize3, color='green', linewidth=4)
#axs[0].plot(time,columnmaximumhailsize4, color='blue', linewidth=4)
#axs[0].legend(['1 Ice Category','2 Ice Categories', '3 Ice Categories','4 Ice Categories'], ncol=1, loc='upper right', fontsize=21)
#axs[0].set_title('Column Maximum Hail Size', fontsize=29)
#axs[0].set_xlabel('Simulation time (mins)', fontsize=26, labelpad=12)
#axs[0].set_ylabel('Maximum Hail Size (mm)', fontsize=26, labelpad=12)
#axs[0].set_xticks([0,25,50,75,100,125,150,175],labels=['0','25','50','75','100','125','150','175'],fontsize=26)
#axs[0].set_xlim(0,180)
#axs[0].set_yticks([0,20,40,60,80,100,120,140],labels=['0','20','40','60','80','100','120','140'],fontsize=26)
#axs[0].set_ylim(0,160)
#axs[0].text(0.02, 0.97, 'a)', transform=axs[0].transAxes, fontsize=26, va='top')

#axs[1].plot(time,maxreflectivity1, color='red',linewidth=4)
#axs[1].plot(time,maxreflectivity2, color='orange',linewidth=4)
#axs[1].plot(time,maxreflectivity3, color='green', linewidth=4)
#axs[1].plot(time,maxreflectivity4, color='blue', linewidth=4)
#axs[1].legend(['1 Ice Category','2 Ice Categories','3 Ice Categories','4 Ice Categories'], ncol=1, loc='lower right', fontsize=21)
#axs[1].set_title('Domain Maximum Reflectivity', fontsize=29)
#axs[1].set_xlabel('Simulation time (mins)', fontsize=26, labelpad=12)
#axs[1].set_ylabel('Reflectivity (dBZ)', fontsize=26, labelpad=12)
#axs[1].set_xticks([0,25,50,75,100,125,150,175],labels=['0','25','50','75','100','125','150','175'],fontsize=26)
#axs[1].set_xlim(0,180)
#axs[1].set_yticks([20,30,40,50,60,70,80],labels=['20','30','40','50','60','70','80'],fontsize=26)
#axs[1].set_ylim(20,90)
#axs[1].text(0.02, 0.97, 'b)', transform=axs[1].transAxes, fontsize=26, va='top')

fig.tight_layout()
plt.savefig(''+ date +'/cm1timeseriesfigure2_'+date+'_comparemphrdpsmp50_2500m.png')
