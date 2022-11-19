import os
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time


def hjd_calc(alpha,delta,cat,loc):
    obj_coo = SkyCoord(alpha, delta, unit=(u.deg, u.deg), frame='icrs')
    date_obs=fits.open(cat[:-3]+'new')[0].header['DATE-OBS']
    exptime=fits.open(cat[:-3]+'new')[0].header['EXPTIME']
    t = Time(date_obs, format='isot', scale='utc', location=loc)
    ltt_helio = t.light_travel_time(obj_coo, 'heliocentric')
    hjd=(t+ltt_helio+exptime/2/86400).jd
    return hjd

def phase_calc(coords,obj,hjd):
    ind=np.where(coords[:,0]==obj)[0][0]
    P=np.float(coords[ind,3])
    T0=np.float(coords[ind,4])
    ph=(hjd-T0)/P-int((hjd-T0)/P)
    if ph<0:
        ph=1+ph
    if ph>0.5:
        ph=ph-1
    return ph

def k_objfromcat(alpha0,delta0,cat,r_arcsec):
    k0=-1
    for k in range(len(cat)):
        alpha=cat[k,12]
        delta=cat[k,13]
        r=np.sqrt((alpha-alpha0)**2+(delta-delta0)**2)
        if r<=r_arcsec/3600:
            k0=k
    return k0

def Is_Obj_OnOfImage(alpha0,delta0,cat):   
    alpha_min=np.min(cat[:,12])
    alpha_max=np.max(cat[:,12])
    delta_min=np.min(cat[:,13])
    delta_max=np.max(cat[:,13])
    onImage=((alpha0>=alpha_min) and (alpha0<=alpha_max) and (delta0>=delta_min) and (delta0<=delta_max))
    return onImage

def limit_mag_estimation(cat,thresh_w, k_ref, mag_ref, errmag_ref):
    threshold=cat[0,-1]
    fwhms=cat[:,14]
    fwhm=np.median(fwhms)
    h=threshold*np.exp(4*np.log10(2)*thresh_w**2/fwhm**2)
    V=2*np.pi*h*(fwhm/2.35482)**2
    limit_mag_cat=-2.5*np.log10(V)
    limitmag=limit_mag_cat-cat[k_ref,6]+mag_ref
    errlimitmag=np.sqrt(cat[k_ref,8]**2+errmag_ref**2)
    return limitmag, errlimitmag

def rel_phot(alpha_obj, delta_obj, alpha_ref, delta_ref, mag_ref, errmag_ref, cat, r_arcsec):
    k_ref=k_objfromcat(alpha_ref, delta_ref, cat, r_arcsec); k_obj=k_objfromcat(alpha_obj, delta_obj, cat, r_arcsec)
    if k_ref!=-1 and k_obj!=-1: #the reference and the object are detected
        mag_obj=cat[k_obj,6]-cat[k_ref,6]+mag_ref
        errmag_obj=np.sqrt(cat[k_obj,8]**2+cat[k_ref,8]**2+errmag_ref**2)
        SNR_obj=cat[k_obj,5]/cat[k_obj,7]
    elif k_ref==-1: # the reference is not detected
        refOnImage=Is_Obj_OnOfImage(alpha_ref,delta_ref,cat) 
        if refOnImage:
            mag_obj=1000; errmag_obj=0; SNR_obj=0   # the reference is inside the frame
        else:
            mag_obj=10000; errmag_obj=0; SNR_obj=-1  # the reference is outside the frame
    elif k_obj==-1: # the object is not detected
        objOnImage=Is_Obj_OnOfImage(alpha_obj,delta_obj,cat)
        if objOnImage:
            mag_obj=None; errmag_obj=0; SNR_obj=0 # the object is inside the frame
        else:
            mag_obj=0; errmag_obj=0; SNR_obj=0 # the object is outside the frame
    return mag_obj, errmag_obj, SNR_obj

def writeTable_date(obj, date,Filter,res, code):
    res_table=open('result_'+obj+'_'+code+'_'+date+'_'+Filter+'.txt', 'w')
    h,w=np.shape(res)
    for i in range(h):   
        for j in range(w):     
            res_table.write(str(round(res[i,j],5))+'\t')
        res_table.write('\n')
    res_table.close()

obj='EX_Dra'; Filter='Green'

path='...'

r_arcsec=6
DETECT_MINAREA=12
thresh_w=np.sqrt(DETECT_MINAREA/np.pi) 
loc = EarthLocation(lon="41:25:33", lat="+43:38:59", height=2026)
d_o=np.loadtxt(path+'Dates_Objects.txt',dtype=np.str)  #table with columns: Date_c	Object	Green (number of frames)	Red(number of frames) 
coords=np.loadtxt(path+'Coords_Ephemeris.dat', dtype=np.str) #table with columns:obj	aplha 	delta	P_orb	T0
j=np.where(coords[:,0]==obj)[0]
alpha_obj=float(coords[j,1]); delta_obj=float(coords[j,2])
dates=np.sort(d_o[np.where(d_o[:,1]==obj)[0],0])

checks=np.loadtxt(path+'CataclChecks/'+obj+'_checks.txt') # table with columns: n	alpha	delta	V	errV	R	errR
# alpha_ref, delta_ref, mag_ref are parameters of reference star
alpha_ref=checks[0,1]; delta_ref=checks[0,2]
if Filter=='Green':
    mag_ref=checks[0,3]; errmag_ref=checks[0,4]
if Filter=='Red':
    mag_ref=checks[0,5]; errmag_ref=checks[0,6]

os.chdir(path+'/relPhotResults')
result_obj=open('result_'+obj+'_obj_'+Filter+'.txt', 'w')

dates_no0=[]
hjds=[]; phases=[]
res_obj=[]
limitmags=[]; limitmagerrs=[]

for date in dates:
    hjds_date=[]; phases_date=[]
    limitmags_date=[]; limitmagerrs_date=[]
    date_o=d_o[np.where(d_o[:,0]==date)[0], :]
    if Filter=='Green':
        n=2
    if Filter=='Red':
        n=3
    n_img=int(date_o[np.where(date_o[:,1]==obj)[0][0],n])
    if n_img==0:
        continue
    dates_no0.append(date)
    res_obj_date=np.zeros((n_img,3))
    os.chdir(path+date+'/'+obj)
    k=0
    cats=[]
    for f in np.sort(os.listdir()):
        if f.endswith('_WCS.cat'):
            os.chdir(path+date+'/'+obj)
            cats.append(f)
            filt=fits.open(f[:-3]+'new')[0].header['FILTER']
            if filt==Filter:
                cat=np.loadtxt(f)
                hjd=hjd_calc(alpha_obj,delta_obj,f,loc)
                hjds_date.append(hjd)
                phase=phase_calc(coords,obj,hjd)
                phases_date.append(phase)
                mag_obj, errmag_obj, SNR_obj=rel_phot(alpha_obj, delta_obj, alpha_ref, delta_ref, mag_ref, errmag_ref, cat, r_arcsec)
                res_obj_date[k,:]=[mag_obj, errmag_obj, SNR_obj]
                if mag_obj==None: # the object is inside the frame, but not detected
                    limitmag, limitmagerr=limit_mag_estimation(cat, thresh_w, k_objfromcat(alpha_ref, delta_ref, cat, r_arcsec), mag_ref, errmag_ref)
                    limitmags_date.append(limitmag); limitmagerrs_date.append(limitmagerr)
                    mag_obj=100
                else:
                    limitmag=0; limitmagerr=0
                    limitmags_date.append(0); limitmagerrs_date.append(0)
                result_obj.write(date+'\t'+f+'\t'+str(round(hjd,5))+'\t'+str(round(phase,5))+'\t'+str(round(mag_obj,3))+'\t'+str(round(errmag_obj,3))+'\t'+str(round(SNR_obj,3))+'\t')
                if limitmag!=0:
                    result_obj.write('# limit_mag= '+str(round(limitmag,3))+'\t'+str(round(limitmagerr,3)))
                result_obj.write('\n')
                k=k+1
            
    os.chdir(path+'/relPhotResults')
    writeTable_date(obj,date,Filter,res_obj_date, 'obj')
    result_obj.write('# '+date+'\t'+str(round(np.mean(hjds_date),5))+'\t'+str(round(np.mean(phases_date),5))+'\tmean mag= '+str(round(np.nanmean(res_obj_date[:,0]),3))+'\t'+str(round(np.nanstd(res_obj_date[:,0]),3))+'\t')
    if len(limitmags_date)>0:
        result_obj.write('limitmag= '+str(round(np.mean(limitmags_date),3))+'\t'+str(round(np.std(limitmags_date),3)))
    result_obj.write('\n')
    limitmags.append(limitmags_date); limitmagerrs.append(limitmagerrs_date)
    phases.append(phases_date)
    hjds.append(hjds_date)
    res_obj.append(res_obj_date)

os.chdir(path+'/relPhotResults') 
args=phases # or args=hjds

plt.figure(figsize=(6,6))
plt.title(obj+' '+Filter)
for i in range(len(dates_no0)):
    plt.scatter(args[i],res_obj[i][:,0],linestyle='solid', label=dates_no0[i][:-2], s=50) 
    plt.errorbar(args[i],res_obj[i][:,0],yerr=res_obj[i][:,1],fmt='o') 
    plt.scatter(args[i],limitmags[i], marker='v',label='limitmag_'+dates_no0[i][:-2])
plt.legend()
plt.gca().invert_yaxis()
plt.show()
