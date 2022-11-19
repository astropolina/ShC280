import os
import numpy as np
from shutil import copyfile

date='20211029_c' # date folder name, '_c' means that overscan is removed
path='.../'
fway=path+date+'/' # directory name

os.chdir(fway)
flist=os.listdir()
del flist[flist.index('bdf')]
polDegree=6

#configuration parameters for SExtractor
DETECT_MINAREA='12'
DETECT_THRESH='3.0'
ANALYSIS_THRESH='3.0'
FILTER='N'
DEBLEND_NTHRESH='32'
DEBLEND_MINCONT='0.1'
CLEAN='N'              # Clean spurious detections? (Y or N)?
CLEAN_PARAM='1.0'            # Cleaning efficiency
PHOT_APERTURES='9'              # MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS='2.5, 3.5'       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
SATUR_LEVEL='64000.0'        # level (in ADUs) at which arises saturation
#SATUR_KEY='SATURATE'       # keyword for saturation level (in ADUs)
SATUR_KEY='64000.0'       # keyword for saturation level (in ADUs)
MAG_ZEROPOINT='0'            # magnitude zero-point
GAIN='0.5'            # detector gain in e-/ADU
#GAIN_KEY='GAIN'           # keyword for detector gain in e-/ADU
GAIN_KEY='0.5'           # keyword for detector gain in e-/ADU
PIXEL_SCALE='0.0'           # size of pixel in arcsec (0=use FITS WCS info)
BACK_SIZE='160' # '64'             # Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE='3'             # Background filter: <size> or <width>,<height>
BACKPHOTO_TYPE='GLOBAL'         # can be GLOBAL or LOCAL
CHECKIMAGE_TYPEs_LIST=['BACKGROUND','OBJECTS','APERTURES']  #can be NONE, BACKGROUND, BACKGROUND_RMS,MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,FILTERED, OBJECTS, -OBJECTS, SEGMENTATION or APERTURES
CHECKIMAGE_TYPE=CHECKIMAGE_TYPEs_LIST[0]
for i in range(len(CHECKIMAGE_TYPEs_LIST)-1):
    CHECKIMAGE_TYPE=CHECKIMAGE_TYPE+','+CHECKIMAGE_TYPEs_LIST[i+1]

df='_df' # '_df' is in the end of fits-files names where dark and flat are taken into account in our work

flist=['EX_Dra'] # names of objects which fits-files we wish to have WCS solvation and be sextracted

#astrometry.net
for obj in flist:
    os.chdir(fway+obj) #folder with date contains folder with objects fits
    objfs=np.loadtxt(obj+df,dtype=np.str,delimiter='\t')
    objfs_WCS=open(obj+df+'_WCS','w')
    for objf in objfs:
        fileNameWCS=objf[:-4]+'_WCS.new'
        os.system('solve-field '+objf+' --out '+fileNameWCS+' --use-sextractor -t '+str(polDegree)+' --overwrite') #run Astrometry.net
        if os.path.isfile(fileNameWCS):
            objfs_WCS.write(fileNameWCS+'\n')
    objfs_WCS.close()

#SExtractor
for obj in flist:
    os.chdir(fway+obj)
    if os.path.isfile('default.param')==True:
        os.remove('default.param')
    copyfile(path+'default.param','default.param') # default.param is file with SExtractor's output parameters, it replace by the new one
    objfs_WCS=np.loadtxt(obj+df+'_WCS',dtype=np.str,delimiter='\t')
    cats=open(obj+df+'_WCS_cat','w') # file with names of catalogues created in SExtractor
    if len(objfs_WCS)==0: # if there are no fits with WCS solvation
        break
    for f in objfs_WCS:
        config_name='config_'+f[0:-3]+'sex'
        config=open(config_name, 'w') # configuration file creating
        catalog_name=f[0:-3]+'cat'
        cats.write(catalog_name+'\n')
        config.write('CATALOG_NAME\t'+catalog_name+'\n')
        config.write('CATALOG_TYPE\tASCII_HEAD\n')
        config.write('PARAMETERS_NAME\tdefault.param\n')
        config.write('DETECT_TYPE\tCCD\n')
        config.write('DETECT_MINAREA\t'+DETECT_MINAREA)
        config.write('\n')
        config.write('DETECT_THRESH\t'+DETECT_THRESH)
        config.write('\n')
        config.write('ANALYSIS_THRESH\t'+ANALYSIS_THRESH)
        config.write('\n')
        config.write('FILTER\t'+FILTER)
        config.write('\n')
        config.write('DEBLEND_NTHRESH\t'+DEBLEND_NTHRESH)
        config.write('\n')
        config.write('DEBLEND_MINCONT\t'+DEBLEND_MINCONT)
        config.write('\n')
        config.write('CLEAN\t'+CLEAN)
        config.write('\n')
        config.write('CLEAN_PARAM\t'+CLEAN_PARAM)
        config.write('\n')
        config.write('MASK_TYPE\tCORRECT')
        config.write('\n')
        config.write('PHOT_APERTURES\t'+PHOT_APERTURES)
        config.write('\n')
        config.write('PHOT_AUTOPARAMS\t'+PHOT_AUTOPARAMS)
        config.write('\n')
        config.write('SATUR_LEVEL\t'+SATUR_LEVEL)
        config.write('\n')
        config.write('SATUR_KEY\t'+SATUR_KEY)
        config.write('\n')
        config.write('MAG_ZEROPOINT\t'+MAG_ZEROPOINT)
        config.write('\n')
        config.write('GAIN\t'+GAIN)
        config.write('\n')
        config.write('GAIN_KEY\t'+GAIN_KEY)
        config.write('\n')
        config.write('PIXEL_SCALE\t'+PIXEL_SCALE)
        config.write('\n')
        config.write('BACK_SIZE\t'+BACK_SIZE)
        config.write('\n')
        config.write('BACK_FILTERSIZE\t'+BACK_FILTERSIZE)
        config.write('\n')
        config.write('BACKPHOTO_TYPE\t'+BACKPHOTO_TYPE)
        config.write('\n')
        config.write('CHECKIMAGE_TYPE\t'+CHECKIMAGE_TYPE)
        config.write('\n')
        config.write('CHECKIMAGE_NAME\t')
        CHECKIMAGE_NAME=CHECKIMAGE_TYPEs_LIST[0]+'_'+f[:-3]+'fit'
        for i in range(len(CHECKIMAGE_TYPEs_LIST)-1):
            CHECKIMAGE_NAME=CHECKIMAGE_NAME+','+CHECKIMAGE_TYPEs_LIST[i+1]+'_'+f[:-3]+'fit'
        config.write(CHECKIMAGE_NAME)
        config.close()
        os.system('sex '+f+' -c ' + config_name) # run SExtractor
        if os.path.isfile(catalog_name):
            cats.write(catalog_name+'\n')
    cats.close()