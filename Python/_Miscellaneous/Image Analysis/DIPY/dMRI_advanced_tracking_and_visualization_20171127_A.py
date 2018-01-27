#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 16:18:59 2017

@author: ptm
"""


# DIPY Advanced Tracking + Advanced Visualization Only


#https://dipy.herokuapp.com/documentation/0.12.0.dev/examples_built/tracking_quick_start/#example-tracking-quick-start


"""
====================
Tracking Quick Start
====================

This example shows how to perform fast fiber tracking using dipy_
[Garyfallidis12]_.

We will use Constrained Spherical Deconvolution (CSD) [Tournier07]_ for local
reconstruction and then generate deterministic streamlines using the fiber
directions (peaks) from CSD and fractional anisotropic (FA) from DTI as a
stopping criteria for the tracking.

Let's load the necessary modules.
"""

import numpy as np
from dipy.tracking.local import LocalTracking, ThresholdTissueClassifier
from dipy.tracking.utils import random_seeds_from_mask
from dipy.reconst.dti import TensorModel
from dipy.reconst.csdeconv import (ConstrainedSphericalDeconvModel,
                                   auto_response)
from dipy.direction import peaks_from_model
from dipy.data import fetch_stanford_hardi, read_stanford_hardi, get_sphere
from dipy.segment.mask import median_otsu
from dipy.viz import actor, window
from dipy.io.image import save_nifti
from nibabel.streamlines import save as save_trk
from nibabel.streamlines import Tractogram
from dipy.tracking.streamline import Streamlines

"""
Enables/disables interactive visualization
"""

#interactive = False
interactive = True

"""
Load one of the available datasets with 150 gradients on the sphere and 10 b0s
"""

fetch_stanford_hardi()
img, gtab = read_stanford_hardi()

data = img.get_data()

"""
Create a brain mask. This dataset is a bit difficult to segment with the
default ``median_otsu`` parameters (see :ref:`example_brain_extraction_dwi`)
therefore we use here more advanced options.
"""

maskdata, mask = median_otsu(data, 3, 1, False,
                             vol_idx=range(10, 50), dilate=2)

"""
For the Constrained Spherical Deconvolution we need to estimate the response
function (see :ref:`example_reconst_csd`) and create a model.
"""

response, ratio = auto_response(gtab, data, roi_radius=10, fa_thr=0.7)

csd_model = ConstrainedSphericalDeconvModel(gtab, response)

"""
Next, we use ``peaks_from_model`` to fit the data and calculated the fiber
directions in all voxels.
"""

sphere = get_sphere('symmetric724')

csd_peaks = peaks_from_model(model=csd_model,
                             data=data,
                             sphere=sphere,
                             mask=mask,
                             relative_peak_threshold=.5,
                             min_separation_angle=25,
                             parallel=True)

"""
For the tracking part, we will use the fiber directions from the ``csd_model``
but stop tracking in areas where fractional anisotropy is low (< 0.1).
To derive the FA, used here as a stopping criterion, we would need to fit a
tensor model first. Here, we fit the tensor using weighted least squares (WLS).
"""

tensor_model = TensorModel(gtab, fit_method='WLS')
tensor_fit = tensor_model.fit(data, mask)

fa = tensor_fit.fa

"""
In this simple example we can use FA to stop tracking. Here we stop tracking
when FA < 0.1.
"""

tissue_classifier = ThresholdTissueClassifier(fa, 0.1)

"""
Now, we need to set starting points for propagating each track. We call those
seeds. Using ``random_seeds_from_mask`` we can select a specific number of
seeds (``seeds_count``) in each voxel where the mask ``fa > 0.3`` is true.
"""

seeds = random_seeds_from_mask(fa > 0.3, seeds_count=1)

"""
For quality assurance we can also visualize a slice from the direction field
which we will use as the basis to perform the tracking.
"""

ren = window.Renderer()
ren.add(actor.peak_slicer(csd_peaks.peak_dirs,
                          csd_peaks.peak_values,
                          colors=None))

if interactive:
    window.show(ren, size=(900, 900))
else:
    window.record(ren, out_path='/Users/ptm/desktop/Current_working_directory/DIPY/csd_direction_field.png', size=(900, 900))

##
##
##
import cv2 as cv2
# Plot Image
#%matplotlib inline
image_name = "/Users/ptm/desktop/Current_working_directory/DIPY/csd_direction_field.png" #Image to be used as query
def plotImage(image_location):
    image = cv2.imread(image_name)
    image = cv2.resize(image, (512,512))
    plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB)); plt.axis('off')
    return
plotImage(image_name)
plt.show()
plt.clf()
##
##
##




"""
.. figure:: csd_direction_field.png
 :align: center

 **Direction Field (peaks)**

``EuDX`` [Garyfallidis12]_ is a fast algorithm that we use here to generate
streamlines. This algorithm is what is used here and the default option
when providing the output of peaks directly in LocalTracking.
"""

streamline_generator = LocalTracking(csd_peaks, tissue_classifier,
                                     seeds, affine=np.eye(4),
                                     step_size=0.5)

streamlines = Streamlines(streamline_generator)

"""
The total number of streamlines is shown below.
"""

print(len(streamlines))

"""
To increase the number of streamlines you can change the parameter
``seeds_count`` in ``random_seeds_from_mask``.

We can visualize the streamlines using ``actor.line`` or ``actor.streamtube``.
"""

ren.clear()
ren.add(actor.line(streamlines))

if interactive:
    window.show(ren, size=(900, 900))
else:
    print('Saving illustration as det_streamlines.png')
    window.record(ren, out_path='/Users/ptm/desktop/Current_working_directory/DIPY/det_streamlines.png', size=(900, 900))

##
##
##
import cv2 as cv2
# Plot Image
#%matplotlib inline
image_name = "/Users/ptm/desktop/Current_working_directory/DIPY/det_streamlines.png" #Image to be used as query
def plotImage(image_location):
    image = cv2.imread(image_name)
    image = cv2.resize(image, (512,512))
    plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB)); plt.axis('off')
    return
plotImage(image_name)
plt.show()
plt.clf()
##
##
##


"""
.. figure:: det_streamlines.png
 :align: center

 **Deterministic streamlines using EuDX (new framework)**

To learn more about this process you could start playing with the number of
seed points or, even better, specify seeds to be in specific regions of interest
in the brain.

Save the resulting streamlines in a Trackvis (.trk) format and FA as
Nifti1 (.nii.gz).
"""

save_trk(Tractogram(streamlines, affine_to_rasmm=img.affine),
         '/Users/ptm/desktop/Current_working_directory/DIPY/det_streamlines.trk')

save_nifti('/Users/ptm/desktop/Current_working_directory/DIPY/fa_map.nii.gz', fa, img.affine)



##
##
##
import nibabel as nib
#newImage = join(dname, 'HARDI193_S0.nii.gz')
newImage = '/Users/ptm/desktop/Current_working_directory/DIPY/fa_map.nii.gz'
img = nib.load(newImage)
data = img.get_data()
import matplotlib.pyplot as plt
axial_middle = data.shape[2] // 2
plt.figure('Showing the datasets')
plt.subplot(1, 2, 1).set_axis_off()
plt.imshow(data[:, :, axial_middle].T, cmap='gray', origin='lower')
#plt.imshow(data[:, :, axial_middle, 0].T, cmap='gray', origin='lower')
#plt.subplot(1, 2, 2).set_axis_off()
#plt.imshow(data[:, :, axial_middle, 10].T, cmap='gray', origin='lower')
plt.show()
#plt.savefig('tensor_fa.png', bbox_inches='tight')
plt.clf()
##
##
##


"""
In Windows if you get a runtime error about frozen executable please start
your script by adding your code above in a ``main`` function and use:

``
if __name__ == '__main__':
    import multiprocessing
    multiprocessing.freeze_support()
    main()
``

References
----------

.. [Garyfallidis12] Garyfallidis E., "Towards an accurate brain tractography",
   PhD thesis, University of Cambridge, 2012.

.. [Tournier07] J-D. Tournier, F. Calamante and A. Connelly, "Robust
   determination of the fibre orientation distribution in diffusion MRI:
   Non-negativity constrained super-resolved spherical deconvolution",
   Neuroimage, vol. 35, no. 4, pp. 1459-1472, 2007.


.. include:: ../links_names.inc

"""

