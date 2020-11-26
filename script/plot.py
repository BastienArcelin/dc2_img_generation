import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

def plot_rgb(gal, bands=[5,6,7], ax=None, band_first=True, zoom=1.5, shifts=None, clip = False):
    if ax is None:
        ax = plt.subplot()
    if band_first:
        tr = [1,2,0]
    else:
        tr = [0,1,2]
    
    imsize = float(gal.shape[1]) / 2.
    if clip==True:
        ax.imshow(np.clip(gal[:,:,:].transpose(tr)[:,:,bands], a_min=0.0, a_max=1000000.), extent=(-imsize,imsize,-imsize,imsize), origin='lower left')
    else:
        ax.imshow(gal[:,:,:].transpose(tr)[:,:,bands], extent=(-imsize,imsize,-imsize,imsize), origin='lower left')
    if shifts is not None:
        k =1
        for (x,y) in shifts:
            if k == 1:
                ax.scatter(x, y,  marker='+', c='b')
            else:
                ax.scatter(x, y,  marker='+', c='r')
            k =2

    #ax.scatter(0., 0., marker='+', c='b')
    ax.set_xlim(-imsize/zoom,imsize/zoom)
    ax.set_ylim(-imsize/zoom,imsize/zoom)
    ax.axis('off')
