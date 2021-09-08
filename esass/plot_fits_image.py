from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import ndimage, misc

path='/Users/baotong/Desktop/period_LW/'
hdu = fits.open(path+'LW_e1.fits')[0]
image_data=hdu.data
fig = plt.figure(figsize=(10, 5))
ax1, ax2= fig.subplots(1, 2)
img_rot = ndimage.rotate(image_data, 18, reshape=False)
ax1.imshow(image_data, cmap='gray',norm=LogNorm())
ax2.imshow(img_rot, cmap='gray',norm=LogNorm())
print(image_data.shape)
print(img_rot.shape)
plt.show()