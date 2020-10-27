#### Import librairies
import sys
import os
import numpy as np
import pandas as pd
import warnings
import healpy as hp
from astropy.table import Table

sys.path.insert(0,'../scripts/tools_for_VAE/')
from tools_for_VAE import cutout_img_dc2

import FoFCatalogMatching
import GCRCatalogs
import lsst.geom
import lsst.daf.persistence as dafPersist

tract = str(sys.argv[1]) # test: 4855 # training: 5074 # validation: 4637
training_test_val = str(sys.argv[2]) # test, training, validation
N = int(sys.argv[3]) # Usually 10 000 images per file

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
object_data = gc_obs.get_quantities(['ra', 'dec', 'blendedness'],
                filters=(mag_filters), native_filters=['tract == '+str(tract)]) # test: 4855 # training: 5074 # validation: 4637

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
              'ellipticity_1_true', 'ellipticity_2_true', 
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
import time
t_3 = time.time()

img_sample = np.zeros((N,59,59,6))
psf_sample = np.zeros((N,59,59,6))

e1 = []
e2 = []
shear1=[]
shear2=[]
redshift=[]
idx = []
blend = []
redshift_true=[]

indices = np.random.choice(list(range(len(truth_idx))), size=N, replace=False)

print('beginning of for loop')
for z, i in enumerate (indices):
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
    
    img_sample[z]=img
    psf_sample[z]=psf
    
    idx.append(truth_data['galaxy_id'][truth_idx[i]])
    e1.append(truth_data['ellipticity_1_true'][truth_idx[i]])
    e2.append(truth_data['ellipticity_2_true'][truth_idx[i]])
    shear1.append(truth_data['shear_1'][truth_idx[i]])
    shear2.append(truth_data['shear_2'][truth_idx[i]])
    redshift.append(truth_data['redshift'][truth_idx[i]])
    blend.append(object_data['blendedness'][object_idx[i]])
    redshift_true.append(truth_data['redshift_true'][truth_idx[i]])
    
t_4 = time.time()

print(t_4-t_3)

df = pd.DataFrame()
df['id']=np.array(idx) # Galaxy index in truth catalog
df['e1']=np.array(e1) # True ellipticity before lensing
df['e2']=np.array(e2) # True ellipticity before lensing
df['shear_1']=np.array(shear1) # True shear applied
df['shear_2']=np.array(shear2) # True shear applied
df['redshift']=np.array(redshift) # Cosmological redshift with line-of-sight motion
df['redshift_true']=np.array(redshift_true) # Cosmological redshift
df['blendedness']=np.array(blend) 

# Save the arrays
np.save('/sps/lsst/users/barcelin/data/dc2_test/'+str(training_test_val)+'/img_sample_5.npy', img_sample)
np.save('/sps/lsst/users/barcelin/data/dc2_test/'+str(training_test_val)+'/psf_sample_5.npy', psf_sample)
df.to_csv('/sps/lsst/users/barcelin/data/dc2_test/'+str(training_test_val)+'/img_data_5.csv', index=False)
