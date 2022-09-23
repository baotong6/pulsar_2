import matplotlib.pylab as plt
import numpy as np

# If you're not familiar with np.r_, don't worry too much about this. It's just
# a series with points from 0 to 1 spaced at 0.1, and 9 to 10 with the same spacing.
x = np.r_[0:1:0.1, 9:10:0.1]
y = np.sin(x)

fig,(ax,ax2) = plt.subplots(1, 2, sharey=True)
plt.title('sdadsdadasdasdasdasdasdasd')
# plot the same data on both axes
ax.plot(x, y, 'bo')
ax2.plot(x, y, 'bo')
# zoom-in / limit the view to different portions of the data
ax.set_xlim(0,1) # most of the data
ax2.set_xlim(9,10) # outliers only

# hide the spines between ax and ax2
ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax.yaxis.tick_left()
ax.tick_params(labeltop='off') # don't put tick labels at the top
ax2.yaxis.tick_right()

# Make the spacing between the two axes a bit smaller
plt.subplots_adjust(wspace=0.15)

# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1). Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.

d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1-d,1+d),(-d,+d), **kwargs) # top-left diagonal
ax.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-left diagonal

kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
ax2.plot((-d,d),(-d,+d), **kwargs) # top-right diagonal
ax2.plot((-d,d),(1-d,1+d), **kwargs) # bottom-right diagonal

# What's cool about this is that now if we vary the distance between
# ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
# the diagonal lines will move accordingly, and stay right at the tips
# of the spines they are 'breaking'
plt.show()

## wiritten by suzhao ##
fig = plt.figure(figsize=[4,3], dpi=160)
ax = fig.add_subplot(111)
ax.scatter(sample_a, e, s=0.2, c='k')
ax.set_xlim(0, 6.3)
ax.set_ylim(0, 1.)
ax.set_xlabel(r"$a~\rm [pc]$")
ax.set_ylabel(r"$e$")

display_to_fig = fig.transFigure.inverted().transform
ax_to_display = ax.transAxes.transform

left, bottom = display_to_fig(ax_to_display([0, 1]))
width, height = (display_to_fig(ax_to_display([1, 0.2])) - display_to_fig(ax_to_display([0, 0])))
ax2 = fig.add_axes([left, bottom, width, height], sharex=ax)
ax.tick_params(axis='x', which='major', top=True)
ax2.tick_params(axis='x', which='major', labelbottom=False)
ax2.hist(sample_a, density=True, color='k', fill=None, bins=20)

left, bottom = display_to_fig(ax_to_display([1, 0]))
width, height = (display_to_fig(ax_to_display([0.2, 1.0])) - display_to_fig(ax_to_display([0, 0])))
ax3 = fig.add_axes([left, bottom, width, height], sharey=ax)
ax.tick_params(axis='y', which='major', right=True)
ax3.tick_params(axis='y', which='major', labelleft=False)
ax3.hist(e, density=True, color='k', fill=None, bins=20, orientation='horizontal')