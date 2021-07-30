import numpy as np
import os
import galsim



## Compute lensed ellipticities from shear and convergence
def calc_lensed_ellipticity_1(es1, es2, gamma1, gamma2, kappa):
    gamma = gamma1 + gamma2*1j # shear (as a complex number)
    es =  es1 + es2*1j # intrinsic ellipticity (as a complex number)
    g = gamma / (1.0 - kappa) # reduced shear
    e = (es + g) / (1.0 + g.conjugate()*es) # lensed ellipticity
    return np.real(e)

def calc_lensed_ellipticity_2(es1, es2, gamma1, gamma2, kappa):
    gamma = gamma1 + gamma2*1j # shear (as a complex number)
    es =   es1 + es2*1j # intrinsic ellipticity (as a complex number)
    g = gamma / (1.0 - kappa) # reduced shear
    e = (es + g) / (1.0 + g.conjugate()*es) # lensed ellipticity
    return np.imag(e)


# read in the LSST filters
filters = {}
lsst_filters_dir = '../share_galsim/bandpasses/'
filter_names_lsst = 'ugrizy'
for filter_name in filter_names_lsst:
    filter_filename = os.path.join(lsst_filters_dir, 'LSST_{0}.dat'.format(filter_name))
    filters[filter_name] = galsim.Bandpass(filter_filename, wave_type='nm')
    filters[filter_name] = filters[filter_name].thin(rel_err=1e-4)
    

# read in SEDs
datapath = '../share_galsim/'
SED_names = ['CWW_E_ext', 'CWW_Sbc_ext', 'CWW_Scd_ext', 'CWW_Im_ext']
SEDs = {}
for SED_name in SED_names:
    SED_filename = os.path.join(datapath, 'SEDs/{0}.sed'.format(SED_name))
    SED = galsim.SED(SED_filename, wave_type='Ang', flux_type='flambda')
    SEDs[SED_name] = SED.withFluxDensity(target_flux_density=1.0, wavelength=500)
    

# Function to generate noiseless galaxy from an id
def generate_noiseless_img_dc2(galaxy_id, psf_img, truth_i, truth_data, truth_idx, pixel_scale = 0.2, bulge_n = 4, disk_n = 1, xsize = 59, ysize = 59):
    idx = truth_i[galaxy_id]
    rng = galsim.BaseDeviate(0)

    # Depends on galaxy
    ## True shear from extragalactic catalog
    gal_g1 = -truth_data['shear_1'][truth_idx[idx]]
    gal_g2 = truth_data['shear_2'][truth_idx[idx]] 

    ## Disk and bulge part
    bulge_frac       = truth_data['bulge_to_total_ratio_i'][truth_idx[idx]]
    disk_frac        = 1 - bulge_frac
    knot_frac        = 0.
    smooth_disk_frac = disk_frac - knot_frac

    disk_e1 = truth_data['ellipticity_1_disk_true'][truth_idx[idx]]
    disk_e2 = truth_data['ellipticity_2_disk_true'][truth_idx[idx]]
    bulge_e1 = truth_data['ellipticity_1_bulge_true'][truth_idx[idx]]
    bulge_e2 = truth_data['ellipticity_2_bulge_true'][truth_idx[idx]]

    disk_hlr = truth_data['size_disk_true'][truth_idx[idx]]
    bulge_hlr = truth_data['size_bulge_true'][truth_idx[idx]]

    ## Create bulge + disk profiles
    bulge = galsim.Sersic(bulge_n, half_light_radius=bulge_hlr)
    disk = galsim.Sersic(disk_n, half_light_radius=disk_hlr)

    ## Compute ellipticities
    shear_1 = np.array(truth_data['shear_1'][truth_idx[idx]])
    shear_2 = np.array(truth_data['shear_2'][truth_idx[idx]])
    convergence = np.array(truth_data['convergence'][truth_idx[idx]])

    disk_1 = calc_lensed_ellipticity_1(-disk_e1, disk_e2, shear_1, shear_2, convergence)
    disk_2 = calc_lensed_ellipticity_2(-disk_e1, disk_e2, shear_1, shear_2, convergence)

    bulge_1 = calc_lensed_ellipticity_1(-bulge_e1, bulge_e2, shear_1, shear_2, convergence)
    bulge_2 = calc_lensed_ellipticity_2(-bulge_e1, bulge_e2, shear_1, shear_2, convergence)

    ## Add knots if necessary
    #knots = galsim.RandomKnots(n_knots, half_light_radius=disk_hlr, flux=knot_frac, rng=rng)
    
    # Shear bulge and disk
    bulge = bulge.shear(g1=bulge_1, g2=bulge_2)
    disk = disk.shear(g1=disk_1, g2=disk_2)
    
    ## Create the galaxy
    gal = bulge_frac * bulge + (1-bulge_frac) * disk

    ## Create output array of noiseless galaxy
    gal_noiseless = np.zeros((xsize,ysize,len(filters)))

    # Depends on filters
    for i, k in enumerate(filters):
        # As object_data['r_FLUXMAG0']=object_data['y_FLUXMAG0'] and is constant I supposed that the fluxmag0 is the same for each filter
        fluxmag0 = 6.30957344e+10 # object_data['r_FLUXMAG0']
        if k=='y':
            gal_flux = fluxmag0*10**(-(truth_data['mag_true_Y_lsst'][truth_idx[idx]])/2.5) # Scale flux as function of magnitude
        else:
            gal_flux = fluxmag0*10**(-(truth_data['mag_true_'+k+'_lsst'][truth_idx[idx]])/2.5) # Scale flux as function of magnitude

        gal = gal.withFlux(gal_flux)

        # convolve with the PSF
        psf_i = galsim.Image(psf_img[galaxy_id,:,:,i].copy())
        interp = 'lanczos15'#'lanczos15'
        
        psf_int = galsim.InterpolatedImage(psf_i,scale = 0.2)#x_interpolant= interp,

        # Create final image and store it
        final = galsim.Convolve([psf_int, gal])
        image = galsim.ImageF(xsize, ysize, scale=0.2)
        _ = final.drawImage(image=image)#, method = 'phot')
        gal_noiseless[:,:,i] = image.array.data
    
    return gal_noiseless