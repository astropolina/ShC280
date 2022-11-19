import os
import numpy as np
from astropy.io import  fits
import ccdproc as ccdp

def img_mbd(Obj_list,SImage,letter,list_name): # function for SBias or SDark substraction
    m_list=open(list_name+'_'+letter,'w')
    for ff in Obj_list:
        hdul = fits.open(ff)
        header = hdul[0].header
        fframe=fits.getdata(ff)
        fframe_i=fframe-SImage   
        fits.writeto(ff[:-4]+'_'+letter+'.fit',fframe_i,header=header,overwrite=True)
        m_list.write(ff[:-4]+'_'+letter+'.fit\n')
    m_list.close()

def flat_divide(f,SFlat,list_name): # function for division by Flat
    f_list=open(list_name+'f','w')
    if f.ndim==0:
        f=[f]
    for ff in f:
        hdul = fits.open(ff)
        header = hdul[0].header
        fframe=fits.getdata(ff)
        fframe_i=fframe/SFlat*np.median(SFlat)
        fits.writeto(ff[:-4]+'f'+'.fit',fframe_i,header=header,overwrite=True)
        f_list.write(ff[:-4]+'f'+'.fit\n')
    f_list.close()


date='20211029_c' # folder name with '_c' means folder with images with removed overscan
fway='.../'+date+'/'
fits_extension='.fit'
os.chdir(fway+'bdf')

bdf_list1=os.listdir()
bdf_list=[]
for bdff in bdf_list1:
    if bdff.endswith('.fit') and bdff.startswith('bdf_'):
        bdf_list.append(bdff)

bias_list=open('Bias','w')
dark_list=open('Dark','w')
dark_exptimes=[]
flatG_list=open('FlatG','w')
flatR_list=open('FlatR','w')

# bias, dark and flat in 2 bands (G and R) are divided into lists
for f in bdf_list:
    hdul=fits.open(f)
    imagetyp=hdul[0].header['IMAGETYP']
    if imagetyp=='Bias Frame':  
        bias_list.write(f+'\n')
    if imagetyp=='Dark Frame':
        dark_list.write(f+'\n')
        dark_exptime=hdul[0].header['EXPTIME']
        if dark_exptime not in dark_exptimes:
            dark_exptimes.append(dark_exptime)
    if f.endswith('FlatG.fit'):
        flatG_list.write(f+'\n')  
    if f.endswith('FlatR.fit'):
        flatR_list.write(f+'\n') 

bias_list.close()
dark_list.close()
flatG_list.close()
flatR_list.close()

#dark frames are sorted by exposure time
dark_list=np.loadtxt('Dark',dtype=np.str)
for dark_exptime in dark_exptimes:
    dark_exptime_list=open('Dark'+str(dark_exptime)+'s','w')
    for dfile in dark_list:
        hdul=fits.open(dfile)
        dfile_exptime=hdul[0].header['EXPTIME']
        if dark_exptime==dfile_exptime:
            dark_exptime_list.write(dfile+'\n')
    dark_exptime_list.close()

Bias_list=np.loadtxt('Bias',dtype=np.str)
SBias=ccdp.combine(Bias_list,'SBias.fit',method='median',unit='adu') # if we don't want to save SBias.fit: SBias=ccdp.combine(Bias_list,method='median',unit='adu')

for dark_exptime in dark_exptimes:
    strdark_exptime=np.str(dark_exptime)
    Dark_list=np.loadtxt('Dark'+strdark_exptime+'s',dtype=np.str)
    SDark=ccdp.combine(Dark_list,'SDark'+strdark_exptime+'s.fit',method='median',unit='adu')

SBias=fits.getdata('SBias.fit')
FlatG_list=np.loadtxt('FlatG',dtype=np.str)
FlatR_list=np.loadtxt('FlatR',dtype=np.str)
img_mbd(FlatG_list,SBias,'b','FlatG')
img_mbd(FlatR_list,SBias,'b','FlatR')

FlatG_b_list=np.loadtxt('FlatG_b',dtype=np.str)
SFlatG=ccdp.combine(FlatG_b_list,'SFlatG.fit',method='average',unit='adu')
FlatR_b_list=np.loadtxt('FlatR_b',dtype=np.str)
SFlatR=ccdp.combine(FlatR_b_list,'SFlatR.fit',method='average',unit='adu')

SFlatG=fits.getdata('SFlatG.fit')
SFlatR=fits.getdata('SFlatR.fit')

os.chdir(fway)
obj_list=os.listdir() #folder with date contains folders with objects fits and folder 'bdf' (bias, dark and flat frames)
del obj_list[obj_list.index('bdf')] 

for obj in obj_list:
    os.chdir(fway+obj)
    n=open(obj,'w')
    for n1 in os.listdir():
        if n1.endswith('s.fit') or n1.endswith('0.fit'):
            n.write(n1+'\n')
    n.close()
    n=np.loadtxt(obj,dtype=np.str,delimiter='\t')
    obj_d_list=open(obj+'_d','w')
    for ff in n:
        hdul = fits.open(ff)
        header = hdul[0].header
        ff_exptime=header['EXPTIME']
        fframe=fits.getdata(ff)
        if ff_exptime in dark_exptimes:
            SDark=fits.getdata(fway+'bdf/SDark'+np.str(ff_exptime)+'s.fit')
            fframe_d=fframe-SDark   
            ff_d_name=ff[:-4]+'_d.fit'
        else:
            fframe_d=fframe-SBias # if there are no SDark with necessary exptime, SBias is substracted from the frame
            ff_d_name=ff[:-4]+'_b.fit'
        if ' ' in ff_d_name:
            ff_d_name=ff_d_name.replace(' ','_')
        fits.writeto(ff_d_name,fframe_d,header=header,overwrite=True)
        obj_d_list.write(ff_d_name+'\n')
    obj_d_list.close()
    n_d=np.loadtxt(obj+'_d',dtype=np.str,delimiter='\t')
    n_df_list=open(obj+'_df','w')
    for f in n_d:
        hdul=fits.open(f)
        fframe=fits.getdata(f)
        header = hdul[0].header
        filter = header['FILTER']
        if filter=='Green':
            SFlat=SFlatG
        if filter=='Red':
            SFlat=SFlatR
        fframe_f=fframe/SFlat*np.median(SFlat)
        fits.writeto(f[:-4]+'f'+'.fit',fframe_f,header=header,overwrite=True)
        n_df_list.write(f[:-4]+'f'+'.fit\n')
    n_df_list.close()






    
