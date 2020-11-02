# DC2 images generation

This repo contains script enabling to generate galaxies and corresponding PSF images using the [LSST DESC DC2 catalogs](https://arxiv.org/pdf/2010.05926.pdf) as presented below. Images are only 59x59 pixels to match the most generic PSF image size that can be produced from lsst stack

**Important**: The script and method presented here is adapted to use ac [CC IN2P3](https://doc.cc.in2p3.fr/index.html) only. Links and installation procedure are specific to this usage. You need to be a member of LSST collaboration to have the required access to run this code.

This script is highly inspired from the work of Yao-Yuan Mao, Scott Daniel and Michael Wood-Vasey realised in these 3 notebooks ([1](https://github.com/LSSTDESC/gcr-catalogs/blob/master/examples/GCRCatalogs%20Demo.ipynb), [2](https://github.com/LSSTDESC/DC2-analysis/blob/master/tutorials/matching_fof.ipynb) and [3](https://github.com/LSSTDESC/DC2-analysis/blob/master/tutorials/dm_butler_postage_stamps.ipynb)).

## Catalogs used
I access catalog through [gcr-catalogs](https://github.com/LSSTDESC/gcr-catalogs) developped in the LSST DESC collaboration. This tools enables access to data stored at NERSC if needed.

I use the extragalactic catalog ``cosmoDC2_v1.1.4_image`` as true catalog, you can find the corresponding paper [here](https://arxiv.org/pdf/1907.06530.pdf), and the object catalog for Run 2.2i DR6 WFD ``dc2_object_run2.2i_dr6_wfd`` (created as described in section 7 [here](https://arxiv.org/pdf/2010.05926.pdf), see also [Bosh et al.](https://arxiv.org/pdf/1705.06766.pdf)).

Images at CC are located at the following addresses:
- for *u*-filter: ``/sps/lssttest/dataproducts/desc/DC2/Run2.2i/v19.0.0-v1/rerun/run2.2i-coadd-wfd-dr6-v1-u``
- for *grizy*-filters: ``/sps/lssttest/dataproducts/desc/DC2/Run2.2i/v19.0.0-v1/rerun/run2.2i-coadd-wfd-dr6-v1-grizy``

## Requirements
Several steps are necessary in order to be able to use this code:
1. Install the LSST Science Pipelines Software following [these instructions](https://pipelines.lsst.io/install/newinstall.html)

2. Source environnement
``` 
source /path/to/directory/lsst_stack/loadLSST.bash # for bash 
```

3. Install [gcr-Catalogs](https://github.com/LSSTDESC/gcr-catalogs):
```
git clone https://github.com/LSSTDESC/gcr-catalogs.git
cd gcr-catalogs/
python setup.py install
```

4. Install FoFMatching (developped by Yao-Yuan Mao [here](https://github.com/yymao/FoFCatalogMatching/))
```
pip install https://github.com/yymao/FoFCatalogMatching/archive/master.zip
```

## Before starting
1. Add a IMGEN_DATA environment variable, in the shell you are running, which points to the directory where you want your data to be stored.

Example, add to your ``.bashrc``:
```
export IMGEN_DC2_DATA='/path/to/img_gen/dc2_data'
```

2. Modify ``generation.sh`` with the correct link to the ``loadLSST.bash`` file. 

3. Install the rest of the required packages (see below).

## Launch generation
1. You can use directly the python file ``generate_dc2_img.py`` on an [interactive machine](https://doc.cc.in2p3.fr/fr/Computing/job-types/job-interactive.html) with:
```
python generate_dc2_img.py 4855 test 10000
```
to generate 10 000 images of galaxy and their corresponding PSF on the directory ``IMGEN_DC2_DATA/test/`` from the tract 4855. You can look at the tract position (a tract is a rectangular region of the sky with a common map projection) [see fig.14, p.27 here](https://arxiv.org/pdf/2010.05926.pdf).

2. Or the shell file ``generation.sh`` modifying it with the corresponding tract number, storing directory and number of images you want. Then run:
```
qsub generation.sh
```

## Required packages
- lsst science pipeline software
- gcr-catalogs
- FoFMatching
- numpy
- astropy
- matplotlib
- pandas
- healpy
- warnings

## Author
Bastien Arcelin - arcelin at apc.in2p3.fr
