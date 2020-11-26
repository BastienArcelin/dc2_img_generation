#### Import librairies
import sys
import os
import time
import numpy as np
import pandas as pd
import warnings
import healpy as hp
from astropy.table import Table

sys.path.insert(0,'../scripts/tools_for_VAE/')
from tools_for_VAE import cutout_img_dc2, utils

import FoFCatalogMatching
import GCRCatalogs
import lsst.geom
import lsst.daf.persistence as dafPersist

tract = '4855'# str(sys.argv[1]) # test: 4855 # training: 5074 # validation: 4637
training_test_val = 'test_mag_26.5'#str(sys.argv[2]) # test, training, validation
N = 1000#int(sys.argv[3]) # Usually 10 000 images per file

# Read in the observed galaxy catalog data.
with warnings.catch_warnings():
    warnings.filterwarnings('ignore')
    gc_obs = GCRCatalogs.load_catalog('dc2_object_run2.2i_dr6_wfd')

# Read in the truth galaxy catalog data.
with warnings.catch_warnings():
    warnings.filterwarnings('ignore')
    gc = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image')

# Let's define a magnitude cut
mag_filters = [
    (np.isfinite, 'mag_r'),
    'mag_r < 26.5']

# Load ra and dec from object, using both of the filters we just defined.
object_data = gc_obs.get_quantities(['ra', 'dec', 'blendedness', 'snr_r_cModel',
                                    'ext_shapeHSM_HsmShapeRegauss_e1','ext_shapeHSM_HsmShapeRegauss_e2',
                                    'mag_r_cModel'],
                filters=(mag_filters), native_filters=['tract == '+str(tract)]) # test: 4855 # training: 4438 #5074 # validation: 4637

# Match the corresponding area for the truth catalog
max_ra = np.nanmax(object_data['ra'])
min_ra = np.nanmin(object_data['ra'])
max_dec = np.nanmax(object_data['dec'])
min_dec = np.nanmin(object_data['dec'])
vertices = hp.ang2vec(np.array([min_ra, max_ra, max_ra, min_ra]),
                      np.array([min_dec, min_dec, max_dec, max_dec]), lonlat=True)
ipix = hp.query_polygon(32, vertices, inclusive=True)
native_filter = f'(healpix_pixel == {ipix[0]})'
for ipx in ipix:
    native_filter=native_filter+f' | (healpix_pixel == {ipx})'
pos_filters=[f'ra >= {min_ra}',f'ra <={max_ra}', f'dec >= {min_dec}', f'dec <= {max_dec}']

# Define a mag cut for truth catalog 
truth_mag_filters = ['mag_r < 26.5']

# Load wanted quantities from truth catalog (https://github.com/LSSTDESC/gcr-catalogs/blob/master/GCRCatalogs/SCHEMA.md)
quantities = ['galaxy_id', 'ra', 'dec', 
              'redshift', 'redshift_true',
              'mag_r_lsst',
              'ellipticity_1_true', 'ellipticity_2_true',
              'convergence',
              'shear_1', 'shear_2']
truth_data = gc.get_quantities(quantities, filters=truth_mag_filters+pos_filters, 
                                      native_filters=native_filter)

# FoFCatalogMatching.match takes a dictionary of catalogs to match, a friends-of-friends linking length. 
# Because our "catalog" is not an astropy table or pandas dataframe, 
# `len(truth_coord)` won't give the actual length of the table
# so we need to specify `catalog_len_getter` so that the code knows how to get the length of the catalog.
results = FoFCatalogMatching.match(
    catalog_dict={'truth': truth_data, 'object': object_data},
    linking_lengths=0.1, # Linking length of 1 arcsecond, you can play around with the values!
    catalog_len_getter=lambda x: len(x['ra']),
)

## Create a list of (ra,dec) of the wanted galaxies
# first we need to know which rows are from the truth catalog and which are from the object
truth_mask = results['catalog_key'] == 'truth'
object_mask = ~truth_mask

# then np.bincount will give up the number of id occurrences (like historgram but with integer input)
n_groups = results['group_id'].max() + 1
n_truth = np.bincount(results['group_id'][truth_mask], minlength=n_groups)
n_object = np.bincount(results['group_id'][object_mask], minlength=n_groups)

one_to_one_group_mask = np.in1d(results['group_id'], np.flatnonzero((n_truth == 1) & (n_object == 1)))

truth_idx = results['row_index'][one_to_one_group_mask & truth_mask]
object_idx = results['row_index'][one_to_one_group_mask & object_mask]

id_ra_dec = Table.from_pandas(pd.DataFrame.from_dict(object_data))

# Create butlers catalog for accessing images
repo_grizy = '/sps/lssttest/dataproducts/desc/DC2/Run2.2i/v19.0.0-v1/rerun/run2.2i-coadd-wfd-dr6-v1-grizy'
repo_u = '/sps/lssttest/dataproducts/desc/DC2/Run2.2i/v19.0.0-v1/rerun/run2.2i-coadd-wfd-dr6-v1-u'

butler_grizy = dafPersist.Butler(repo_grizy)
butler_u = dafPersist.Butler(repo_u)

# Create the numpy arrays
t_3 = time.time()

img_sample = np.zeros((N,59,59,6))
psf_sample = np.zeros((N,59,59,6))


indices = np.random.choice(list(range(len(truth_idx))), size=N, replace=False)

def gen_function(indices):
    np.random.seed() # important for multiprocessing !
    e1 = []
    e2 = []
    hsm_e1 = []
    hsm_e2 = []
    shear1=[]
    shear2=[]
    redshift=[]
    idx = []
    blend = []
    redshift_true=[]
    convergence=[]
    snr = []
    mag_r_meas = []
    mag_r_true = []
    i = np.random.choice(indices)
    print(i)
    first = id_ra_dec[object_idx[i]]
    ra, dec = first['ra'], first['dec']

    img = np.zeros((59,59,6))
    psf = np.zeros((59,59,6))
    filters = ['u','g','r','i','z','y']
    for k, filter_k in enumerate (filters):
        if k == 0:
            cutout = cutout_img_dc2.cutout_coadd_ra_dec(butler_u, ra, dec, filter=filter_k)
        else:
            cutout = cutout_img_dc2.cutout_coadd_ra_dec(butler_grizy, ra, dec, filter=filter_k)
        
        radec = lsst.geom.SpherePoint(ra, dec, lsst.geom.degrees)
        xy = cutout.getWcs().skyToPixel(radec)  # returns a Point2D
        
        img[:,:,k]= cutout.image.array
        if (cutout.getPsf().computeKernelImage(xy).array.size != 3481):
            print('not taken into account')
            break
        else:
            psf[:,:,k]= cutout.getPsf().computeKernelImage(xy).array
        
    idx.append(truth_data['galaxy_id'][truth_idx[i]])
    e1.append(truth_data['ellipticity_1_true'][truth_idx[i]])
    e2.append(truth_data['ellipticity_2_true'][truth_idx[i]])
    hsm_e1.append(object_data['ext_shapeHSM_HsmShapeRegauss_e1'][object_idx[i]])
    hsm_e2.append(object_data['ext_shapeHSM_HsmShapeRegauss_e2'][object_idx[i]])
    shear1.append(truth_data['shear_1'][truth_idx[i]])
    shear2.append(truth_data['shear_2'][truth_idx[i]])
    redshift.append(truth_data['redshift'][truth_idx[i]])
    blend.append(object_data['blendedness'][object_idx[i]])
    snr.append(object_data['snr_r_cModel'][object_idx[i]])
    redshift_true.append(truth_data['redshift_true'][truth_idx[i]])
    convergence.append(truth_data['convergence'][truth_idx[i]])
    mag_r_meas.append(object_data['mag_r_cModel'][object_idx[i]])
    mag_r_true.append(truth_data['mag_r_lsst'][truth_idx[i]])
    
    return img, psf, idx, e1, e2, hsm_e1, hsm_e2, shear1, shear2, redshift, blend, snr, convergence, redshift_true, mag_r_meas, mag_r_true
img_sample = []
psf_sample = []

# Here we save data for all datasets
res = utils.apply_ntimes(gen_function, N, ([indices]))

e1 = []
e2 = []
hsm_e1 = []
hsm_e2 = []
shear1=[]
shear2=[]
redshift=[]
idx = []
blend = []
redshift_true=[]
convergence=[]
snr = []
mag_r_meas = []
mag_r_true = []
for i in range (N):
    img_temp, psf_temp, idx_temp, e1_temp, e2_temp, hsm_e1_temp, hsm_e2_temp, shear1_temp, shear2_temp, redshift_temp, blend_temp, snr_temp, convergence_temp,redshift_true_temp, mag_r_meas_temp, mag_r_true_temp = res[i]
    idx.append(idx_temp)
    e1.append(e1_temp)
    e2.append(e2_temp)
    hsm_e1.append(hsm_e1_temp)
    hsm_e2.append(hsm_e2_temp)
    shear1.append(shear1_temp)
    shear2.append(shear2_temp)
    redshift.append(redshift_temp)
    redshift_true.append(redshift_true_temp)
    blend.append(blend_temp)
    convergence.append(convergence_temp)
    snr.append(snr_temp)
    mag_r_meas.append(mag_r_meas_temp)
    mag_r_true.append(mag_r_true_temp)
    img_sample.append(img_temp)
    psf_sample.append(psf_temp)


t_4 = time.time()
print(np.array(idx).shape)
print(t_4-t_3)

df = pd.DataFrame()
df['id']=np.array(idx)[:,0] # Galaxy index in truth catalog
df['e1']=np.array(e1)[:,0] # True ellipticity before lensing
df['e2']=np.array(e2)[:,0] # True ellipticity before lensing
df['shear_1']=np.array(shear1)[:,0] # True shear applied
df['shear_2']=np.array(shear2)[:,0] # True shear applied
df['redshift']=np.array(redshift)[:,0] # Cosmological redshift with line-of-sight motion
df['redshift_true']=np.array(redshift_true)[:,0] # Cosmological redshift
df['blendedness']=np.array(blend)[:,0] # blendedness
df['snr_r']=np.array(snr)[:,0] # Signal to noise ratio in r-band
df['convergence']=np.array(convergence)[:,0] # Shear convergence
df['e1_hsm_regauss']=np.array(hsm_e1)[:,0] # e1 parameter measured with HSM REGAUSS method
df['e2_hsm_regauss']=np.array(hsm_e2)[:,0] # e2 parameter measured with HSM REGAUSS method
df['mag_r_meas']=np.array(mag_r_meas)[:,0] # Magnitude measured on image by LSST pipeline
df['mag_r_true']=np.array(mag_r_true)[:,0] # Magnitude in input of simulation

# Save the arrays
data_dir = str(os.environ.get('IMGEN_DC2_DATA'))
np.save(data_dir+str(training_test_val)+'/img_sample.npy', img_sample)
np.save(data_dir+str(training_test_val)+'/psf_sample.npy', psf_sample)
df.to_csv(data_dir+str(training_test_val)+'/img_data.csv', index=False)
