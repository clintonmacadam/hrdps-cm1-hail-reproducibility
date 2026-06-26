import os
import pandas as pd
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean
from metpy.plots import SkewT, Hodograph
import metpy.calc as mpcalc
from metpy.units import units
import matplotlib.colors as mcolors
import matplotlib.patheffects as mpatheffects
import metpy.interpolate as mpinterpolate
from scipy.interpolate import interp1d
import warnings
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.thermo as thermo
import wrf

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# Folder containing all the profiles
folder_dir = r'C:\Users\clint\repos\data\input_sounding\\'

for file in os.listdir(folder_dir):
    if file.endswith('.csv'):  # Process only CSV files
        file_path = os.path.join(folder_dir, file)
        df = pd.read_csv(file_path)
    
        # Create the figure
        fig = plt.figure(figsize=(10, 10))
        skew = SkewT(fig, subplot=(1, 1, 1))
    
        # Plotting variables
        temps =['red']
        dews = ['green']
        tvs = ['fuchsia']
        capes = ['lightcoral']
        cins = ['cyan']
    
        # Extract data
        a = df.altitude
        p = df.pressure
        t = df.temp
        d = df.dpt

        # Calculate wind components and virtual temperature
        u, v, data_tv = [], [], []
        for i in range(len(df.ws)):
            u.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')), df.wd[i] * units.deg)[0]).magnitude)
            v.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')), df.wd[i] * units.deg)[1]).magnitude)
            data_tv.append(thermo.virtemp(p[i], t[i], d[i]))
    
        data_tv = np.array(data_tv)
    
        # Set altitude to AGL
        a = [x - a[0] for x in a]
    
        # Create profile
        wdir = mpcalc.wind_direction(u*units('m/s'),v*units('m/s'))
        wspd = mpcalc.wind_speed(u*units('m/s'),v*units('m/s'))
        
        #creating profile
        data_profset = profile.create_profile(profile='default', pres=p, hght=a,
                                              tmpc=t, dwpc=d, wdir=wdir, wspd=wspd)
    
        sfcpres = data_profset.pres[data_profset.sfc]
        # sfcpres = profset.pres[profset.sfc]
    
    
        #finding values related to the parcel
        parcel_sfcp = data_profset.pres[0]
        parcel_topp = interp.pres(data_profset, interp.to_msl(data_profset, 500.))
        parcel_depth = parcel_sfcp - parcel_topp
        
        # Define various parcels to lift
        data_sfcpcl = params.parcelx(data_profset, flag=1)  # Surface-based parcel
        data_fcstpcl = params.parcelx(data_profset, flag=2) # Forecast parcel (not sure what this means, though)
        data_mupcl = params.parcelx(data_profset, flag=3)   # Most-unstable parcel
        # alternative workaround
        mlpcl_vals = params.DefineParcel(data_profset, flag=4, pres=parcel_depth)
        data_mlpcl = params.parcelx(data_profset, lplvals=mlpcl_vals)

        #plotting data
        skew.plot(p,t,c='r',linewidth=1,ls='-',label='Temperature')
        skew.plot(p,data_tv, c='darkorange', linewidth=1,ls=':',label='Virtual Temperature')
        skew.plot(p,d,c='g',linewidth=1,ls='-',label='Dewpoint')
        # skew.plot(data_sfcpcl.ptrace,data_sfcpcl.ttrace,color=[0.5,0,0.5],linewidth=1, linestyle='--')
        skew.plot(data_mlpcl.ptrace,data_mlpcl.ttrace,color='magenta',linewidth=1, linestyle='-',label='ML Parcel Trace')
        skew.plot_barbs(p[0:30:5],u[0:30:5],v[0:30:5],color='k',sizes={'spacing':0.2},length=6.,xloc=1)
        skew.plot_barbs(p[30:78:2],u[30:78:2],v[30:78:2],color='k',sizes={'spacing':0.2},length=6.,xloc=1)
    
    
        parc_tvshade = wrf.interp1d(np.round(data_mlpcl.ttrace.data,4), np.round(data_mlpcl.ptrace.data,4), np.round(p,5))
        # parc_tdshade = wrf.interp1d(np.round(data_mlpcl.ptrace.data,5), np.round(data_mlpcl.ptrace.data,5), np.round(data_pres,5))
        parc_pshade = wrf.interp1d(np.round(data_mlpcl.ptrace.data,4), np.round(data_mlpcl.ptrace.data,4), np.round(p,5))
        
        parc_tvshade[0] = data_mlpcl.ttrace.data[0]
        parc_pshade[0] = data_mlpcl.ptrace.data[0]
    
            
        lcl_pres = data_mlpcl.lclpres
        lcl_temp = wrf.interp1d(np.array(data_tv), np.array(p), np.asarray([lcl_pres]))
        # print('Region:',region[fin], "ML LCL pressure:", lcl_pres)
        lfc_pres = data_mlpcl.lfcpres
        lfc_temp = wrf.interp1d(np.array(data_tv), np.array(p), np.asarray([lfc_pres]))
        el_pres = data_mlpcl.elpres
        el_temp = wrf.interp1d(np.array(data_tv), np.array(p), np.asarray([el_pres]))
        parc_eltemp = wrf.interp1d(parc_tvshade, parc_pshade, np.asarray([el_pres]))
    
        # skews[fin].plot(data_mlpcl.pres, data_mlpcl.tmpc, 'ro', markersize=2, markerfacecolor='red')
        
        cap = data_mlpcl.cap
        ccd = data_mlpcl.elhght-data_mlpcl.lfchght
        
        # sys.exit()
        
        if h == 0:
            align = 'left'
            value = 20.0
        else:
            align = 'right'
            value = 8
        

        #plotting ML LCL, LFC, EL
        skew.plot([lcl_pres, lcl_pres], [lcl_temp-1, lcl_temp+1], 'k-', lw=1)
        skew.ax.text(lcl_temp.data+value, lcl_pres, 'ML LCL',fontsize=8,verticalalignment='center',horizontalalignment=f'{align}',bbox=dict(facecolor='none',edgecolor='none',alpha=0.5,pad=0.1))
        skew.plot([lfc_pres, lfc_pres], [lfc_temp-1, lfc_temp+1], 'k-', lw=1)
        skew.ax.text(lfc_temp.data+value, lfc_pres, 'ML LFC',fontsize=8,verticalalignment='center',horizontalalignment=f'{align}',bbox=dict(facecolor='none',edgecolor='none',alpha=0.5,pad=0.1))
        skew.plot([el_pres, el_pres], [el_temp-1, el_temp+1], 'k-', lw=1)
        skew.ax.text(el_temp.data+value, el_pres, 'ML EL',fontsize=8,verticalalignment='center',horizontalalignment=f'{align}',bbox=dict(facecolor='none',edgecolor='none',alpha=0.5,pad=0.1))
        
        # skews[fin].shade_cin(data_pshade., data_tshade, data_mlpcl.ttrace.data, dewpoint=data_tdshade)
        # skews[fin].shade_cape(data_pshade.data, data_tshade, data_mlpcl.ttrace.data, dewpoint=data_tdshade)
    
        # skews[fin].shade_cin(data_mlpcl.ptrace, envtempv_shade, data_mlpcl.ttrace.data)
        # skews[fin].shade_cape(data_mlpcl.ptrace, envtempv_shade, data_mlpcl.ttrace.data)

        # find the level where environment Tv exceeds the parcel TV again.
        EL_lvlind = np.min(np.where(p < el_pres))-1
        
        data_pres_cinshade = np.append(p[0:EL_lvlind], el_pres)
        data_temp_cinshade = np.append(data_tv[0:EL_lvlind], el_temp)  # really the environmental virtual temperature
        parc_tv_cinshade = np.append(parc_tvshade[0:EL_lvlind], parc_eltemp)
    
        # skews[fin].shade_cin(data_pres[0:EL_lvlind], data_tv[0:EL_lvlind], parc_tvshade[0:EL_lvlind])
        skew.shade_cin(data_pres_cinshade, data_temp_cinshade, parc_tv_cinshade,color='cyan',alpha=0.3,label='CIN')
        skew.shade_cape(p, data_tv, parc_tvshade, alpha = 0.7,color='lightcoral',label='CAPE')
    
        skew.ax.set_xlim(-50,50)
        skew.ax.set_ylim(1000,100)
        skew.ax.set_xticks(np.arange(-50,50.01,10))
        skew.ax.set_yticks(np.arange(1000,99.99,-100))
        skew.ax.set_xticklabels(['-50','-40','-30','-20','-10','0','10','20','30','40','50'],fontsize=11)
        skew.ax.set_yticklabels(['1000','900','800','700','600','500','400','300','200','100'],fontsize=11)
        skew.ax.set_xlabel('Temperature ($\mathregular{^o}$C)',fontsize=15)
        skew.ax.set_ylabel('Pressure (hPa)', fontsize=15)
        #skew.ax.axvline(0,color='c',linestyle='--',linewidth=1.5)
        skew.plot_dry_adiabats(t0=(np.arange(-50,260,10))*units('degC'),color='gray',linewidth=0.5,linestyle='--', alpha = 0.3)
        skew.plot_moist_adiabats(color='gray',linewidth=0.5,linestyle='--', alpha = 0.3)
        skew.plot_mixing_lines(color='gray',linewidth=0.5,linestyle=':')
        skew.ax.axvline(-30, color='darkblue', linestyle='--', alpha=0.7)
        skew.ax.axvline(-10, color='darkblue', linestyle='--', alpha=0.7)
        lf=(1-(lcl_pres/1013.25)**0.190284)*145366.45 #coverting from hPa to feet
        lkm=(lf*0.3048) # from feet to m
        sfcp=(1-(p[0]/1013.25)**0.190284)*145366.45 #coverting from hPa to feet for the surface
        sfc=(sfcp*0.3048)# from feet to m
        skew.ax.axhline(y=p[0], color='blue', linestyle='--',lw=1, alpha=0.7,label='sfc = ' + f'{sfc:.1f} m')
        skew.ax.axhline(y=lcl_pres, color='k', linestyle='--', lw=1,alpha=0.7,label='ML_LCL = ' + f'{lkm:.1f} m')
        skew.ax.text(-21.5,300, 'Hail Growth\nLayer', rotation=60)


        skew.ax.legend(loc='upper right', bbox_to_anchor=(0.95,1))
        profile_name = os.path.basename(file).replace('.csv', '')
        skew.ax.set_title(f'{profile_name}',fontsize=20)
        folder_dire = r'C:\Users\ayitahe-INS\Desktop\NHP\Master_Dataset\John_data_recalculated\NSC_Profiles\Skt_Plots'
        
        # Save the plot as an image
        output_path = os.path.join(folder_dire, f'{profile_name}_SkewT.png')
        plt.savefig(output_path, dpi=300)
        plt.close(fig)  
