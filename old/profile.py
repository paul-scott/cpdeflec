import pattern
import camera
import pattern
from numpy import *
from numpy.linalg import norm
from scipy import optimize
import Image
import gc
import sys

sys.setrecursionlimit(8000)
gc.set_debug(gc.DEBUG_LEAK)

# Global parameters.
pixwid = 4288
pixhei = 2848
pixsize = 0.00554
# Get the principal values from photogrammetry calibration of camera.
prdist = 20.52
# Need to sort out coord system of principal point from photogrammetry.
# For now assume that prpoint is position where optical axis strikes array
# relative to centre of array. And that increase in direction related to px
# corresponds in increase in ax dir.
prpoint = array([0.1065,-0.2374], dtype=float32)
dotsize = 5.
camdistguess = 1500.
# For dists 0-1, 1-3, 3-1
dotsep = [1200., sqrt(2)*1200., 1200.]
# Specify rough locations of mirror corners in global coords.
# Might want to underestimate mirredge and overestimate mirrstart.
mirredge = 1175.
mirrstart = array([25.,25.,0.], dtype=float32)
# Width of repeating pattern segment.
segsize = 50
# Get pattern values from photogrammetry and a few calculations.
# Position is horizontal edge of horiz rotated pattern, and vertical edge of
# vert rotated pattern.
patternpos = array([2000., -500., -3000.], dtype=float32)
patterntrans = array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype=float32)
# Given number 0 -> 1, return a number 0 -> 1.
# For now using simple linear relation.
patternrela = poly1d([1,0])
# Starting depth of mirror edge. Closely linked to mirrstart[2] if starting
# pixel is at top left corner. 
startdepth = 0.

# System coordinates roughly copy coordinates of image.

# THINK ABOUT MAKING CLASS THAT GETS INITIALISED WITH CALIBRATION FILE....

def findprofile(filenames, initdots, outpath, calibfile):
	# Initialise camera and pattern. Likely use values from configuration file.
	if calibfile == '':
		print('No calibration file given... Using program defaults.')
	else:
		loadcalibfile(calibfile)

	cam = camera.Camera(pixwid, pixhei, pixsize, prdist, prpoint)
	patt = pattern.Pattern(segsize, patternpos, patterntrans, patternrela)

# Need to consider what happens when image has alpha channel.
	image = Image.open(filenames[0])
	print('Horizontal image loaded.')
	dots = []
	for initdot in initdots:
		# Dots go in as tuples and come out as arrays.
		dots.append(centroid(image, initdot,
			cam.getobjpix(dotsize, camdistguess)))

	print(dots)
# NOT WORKING YET.
# Could be fucking up the calculation of mcorn.
	cam.findposition(dots, dotsep, [camdistguess,camdistguess,
		camdistguess])
	print('Camera position determined.')
	print(cam.campos)

	# Determine global coordinates of mirror corners. Assumes that the mirror
	# has been lined up with test frame.
	mirrcorn = [mirrstart,mirrstart+array([0.,mirredge,0.]),
		mirrstart+array([mirredge,mirredge,0.]),
		mirrstart+array([mirredge,0.,0.])]
	print('Corners Found.')

	# Work out pixel coords of mirror corners. Give back a list of arrays.
	mcorn = []
	for corn in mirrcorn:
		mcorn.append(cam.getpix(corn))

	# Work out bounds to crop image to.
	# Note that x coordinate refers to global coords not camera coords.
	# DON'T NEED TO CHECK EVERY CORNER EACH TIME.
	xpixmax = max(mcorn[0][0],mcorn[1][0],mcorn[2][0],mcorn[3][0])
	xpixmin = min(mcorn[0][0],mcorn[1][0],mcorn[2][0],mcorn[3][0])
	ypixmax = max(mcorn[0][1],mcorn[1][1],mcorn[2][1],mcorn[3][1])
	ypixmin = min(mcorn[0][1],mcorn[1][1],mcorn[2][1],mcorn[3][1])
	print('Max and mins found.')

	# Convert corner pixels into indices of cropped new image.
	# This is only required for determining masking.
	for i in xrange(len(mcorn)):
		mcorn[i] = mcorn[i] - array([xpixmin,ypixmin])

# These are sometimes huge!
	xlen = xpixmax-xpixmin + 1
	ylen = ypixmax-ypixmin + 1
	print(xlen)
	print(ylen)
	print('Offset and lengths established.')

	# Crop image so that it doesn't take up as much room and convert to 
	# an array so it is easier to manipulate.
	box = (xpixmin, ypixmin, xpixmax+1, ypixmax+1)
	print(box)
	image = asarray(image.crop(box),dtype=float32)
	print(image.shape)
	# Need to swap rows and cols so that they correspond to x and y.
	# Rows are x, columns are y.
	image = image.transpose((1,0,2))
	print('Horizontal image cropped and converted.')

	# Define lines that bound the mirror surface to measure.
	slope = (mcorn[2][1]-mcorn[1][1])/float(mcorn[2][0]-mcorn[1][0])
	offset = mcorn[1][1] - slope*mcorn[1][0]
	bottline = poly1d([slope,offset])

	slope = (mcorn[0][1]-mcorn[3][1])/float(mcorn[0][0]-mcorn[3][0])
	offset = mcorn[3][1] - slope*mcorn[3][0]
	topline = poly1d([slope,offset])
	
	# These next two could potentially have an infinite slope.
	# The 4 loops up ahead won't run on left or right if this is the case.
	# Therefore we just set the polynomial to a constant.

	leftline = poly1d(0)
	if (mcorn[0][0] != mcorn[1][0]):
		slope = (mcorn[1][1]-mcorn[0][1])/float(mcorn[1][0]-mcorn[0][0])
		offset = mcorn[0][1] - slope*mcorn[0][0]
		leftline = poly1d([slope,offset])

	rightline = poly1d(0)
	if (mcorn[3][0] != mcorn[2][0]):
		slope = (mcorn[3][1]-mcorn[2][1])/float(mcorn[3][0]-mcorn[2][0])
		offset = mcorn[2][1] - slope*mcorn[2][0]
		rightline = poly1d([slope,offset])

	# ybound[0] is lower bound, ybound[1] is upper.
	ybound=zeros((xlen,2),dtype=int)

	# At most one of these loops will execute (other will have a range of []).
	for x in xrange(mcorn[0][0],mcorn[1][0]+1):
		ybound[x,0] = int(round(leftline(x)))
		ybound[x,1] = int(round(bottline(x)))
	for x in xrange(mcorn[1][0],mcorn[0][0]+1):
		ybound[x,0] = int(round(topline(x)))
		ybound[x,1] = int(round(leftline(x)))

	# At most one of these loops will execute (other will have a range of []).
	for x in xrange(mcorn[3][0],mcorn[2][0]+1):
		ybound[x,0] = int(round(topline(x)))
		ybound[x,1] = int(round(rightline(x)))
	for x in xrange(mcorn[2][0],mcorn[3][0]+1):
		ybound[x,0] = int(round(rightline(x)))
		ybound[x,1] = int(round(bottline(x)))

	# Set the middle part. Don't need +1 because it has already been set.
	for x in xrange(max(mcorn[0][0],mcorn[1][0]),min(mcorn[3][0],mcorn[2][0])):
		ybound[x,0] = int(round(topline(x)))
		ybound[x,1] = int(round(bottline(x)))

	print('y bounds determined as function of x')
	print(ybound)
	print((ybound >= ylen).sum())
	print((ybound < 0).sum())

	# Set up array of vectors to hold pattern positions.
	# Should rename this as it also holds normal vectors.
	# First axis in array is xpix, second is ypix.
	# Printing of array will give values with x as rows and y as cols.
	vecs = zeros((xlen, ylen, 3), dtype=float32)

	# Patt horizontal goes in vecs[,,0].
	# Patt vertical goes in vecs[,,1].
	# Patt z stays the same (0).
	
	#return 2
	# Starting in top left corner of mirror.
	calcpatterndist(patt, image, vecs, mcorn[0], ybound, segsize/2., 0)
	print('Information extracted from horizontal image.')

	image = Image.open(filenames[1])
	print('Vertical image loaded.')
	# Crop image so that it doesn't take up as much room and convert to 
	# an array so it is easier to manipulate.
	image = asarray(image.crop(box),dtype=int)
	# Need to swap rows and cols so that they correspond to x and y.
	# Rows are x, columns are y.
	image = image.transpose((1,0,2))
	print('Vertical image cropped and converted.')

	calcpatterndist(patt, image, vecs, mcorn[0], ybound, segsize/2., 1)
	print('Information extracted from vertical image.')

	# Clear up memory taken by image.
	del image

	# Pattern positions are converted into a vector in global coordinates.
	vecs = patt.getvector(vecs)

	# Start solving for positions and normals at top left corner.
	# This should be ok even though we will diverge as we approach centre.
	# First work towards centre then...

	poss = zeros((xlen, ylen, 3), dtype=float32)
	# Probably a good enough case to check if x == 0, since outside box.

	solvesurface(cam, vecs, poss, mcorn[0], array([xpixmin,ypixmin]), ybound,
			startdepth)
	print('Surface normals and positions calculated.')

	print('Curve fitted to surface.')

	print('Ideal curve fitted to surface.')

	print('Slope errors calculated relative to best fit curve.')

	print('Slope errors calculated relative to ideal curve.')

	return 1

# Think about writing output parameters to an ongoing log file instead of
# separate file for each panel or batch.

def loadcalibfile(file):
	"""Loads calibration file and sets parameters to file values."""
	pixwid = 'puthere'
	print('Calibration file: ' + file + ' loaded.')

def calcpatterndist(patt, image, vecs, ipix, yb, idist, orien):
	"""Calculates what part of pattern is reflected."""
	pdist = patt.getdist(image[ipix[0],ipix[1]], idist)
	vecs[ipix[0],ipix[1],orien] = pdist

	for y in xrange(ipix[1]-1,yb[ipix[0],0]-1,-1):
		pdist = patt.getdist(image[ipix[0],y], pdist)
		vecs[ipix[0],y,orien] = pdist

	pdist = vecs[ipix[0],ipix[1],orien]
	for y in xrange(ipix[1]+1,yb[ipix[0],1]+1):
		pdist = patt.getdist(image[ipix[0],y], pdist)
		vecs[ipix[0],y,orien] = pdist

	for x in xrange(ipix[0]-1,-1,-1):
		# Find a value between the bounds.
		ystart = min(max(ipix[1],yb[x,0]),yb[x,1])
		pdist = vecs[x+1,ystart,orien]
		# Search at and above the starting point.
		for y in xrange(ystart,yb[x,0]-1,-1):
			pdist = patt.getdist(image[x,y], pdist)
			vecs[x,y,orien] = pdist

		pdist = vecs[x,ystart,orien]
		# Search below the starting point.
		for y in xrange(ystart+1,yb[x,1]+1):
			pdist = patt.getdist(image[x,y], pdist)
			vecs[x,y,orien] = pdist

	for x in xrange(ipix[0]+1,image.shape[0]):
		# Find a value between the bounds.
		ystart = min(max(ipix[1],yb[x,0]),yb[x,1])
		pdist = vecs[x-1,ystart,orien]
		# Search at and above the starting point.
		for y in xrange(ystart,yb[x,0]-1,-1):
			pdist = patt.getdist(image[x,y], pdist)
			vecs[x,y,orien] = pdist

		pdist = vecs[x,ystart,orien]
		# Search below the starting point.
		for y in xrange(ystart+1,yb[x,1]+1):
			pdist = patt.getdist(image[x,y], pdist)
			vecs[x,y,orien] = pdist
	# Don't need to return array since its values have been manipulated.

def solvesurface(cam, vecs, poss, ipix, offset, yb, startdepth):
# Could solve norms as 2 angles instead of 3 distances.
	# Need to get pixel for camera since image has been cropped.
	campix = ipix + offset
	ch = cam.getpixdir(campix)
	# Determine start pos from guess and using pixel vector data.
	c = ch*(cam.campos[2]-startdepth)/ch[2]
	poss[ipix[0],ipix[1]] = cam.campos + c
	ph = vecs[ipix[0],ipix[1]] - poss[ipix[0],ipix[1]]
	ph = ph/norm(ph)
	vecs[ipix[0],ipix[1]] = surfnorm(ch, ph)

	for y in xrange(ipix[1]-1,yb[ipix[0],0]-1,-1):
		ch = cam.getpixdir(array([ipix[0],y])+offset)
		dist = extrapolate(poss[ipix[0],y+1], vecs[ipix[0],y+1], ch,
				cam.campos)
		poss[ipix[0],y] = dist*ch + cam.campos
		ph = vecs[ipix[0],y] - poss[ipix[0],y]
		ph = ph/norm(ph)
		vecs[ipix[0],y] = surfnorm(ch, ph)

	for y in xrange(ipix[1]+1,yb[ipix[0],1]+1):
		ch = cam.getpixdir(array([ipix[0],y])+offset)
		dist = extrapolate(poss[ipix[0],y-1], vecs[ipix[0],y-1], ch,
				cam.campos)
		poss[ipix[0],y] = dist*ch + cam.campos
		ph = vecs[ipix[0],y] - poss[ipix[0],y]
		ph = ph/norm(ph)
		vecs[ipix[0],y] = surfnorm(ch, ph)

	for x in xrange(ipix[0]-1,-1,-1):
		# Find a value between the bounds.
		ystart = min(max(ipix[1],yb[x,0]),yb[x,1])

		ch = cam.getpixdir(array([x,ystart])+offset)
		dist = extrapolate(poss[x+1,ystart], vecs[x+1,ystart], ch, cam.campos)
		poss[x,ystart] = dist*ch + cam.campos
		ph = vecs[x,ystart] - poss[x,ystart]
		ph = ph/norm(ph)
		vecs[x,ystart] = surfnorm(ch, ph)
		# Search above the starting point.
		for y in xrange(ystart-1,yb[x,0]-1,-1):
			ch = cam.getpixdir(array([x,y])+offset)
			num = 1
			dist = extrapolate(poss[x,y+1], vecs[x,y+1], ch, cam.campos)
			# Optional to add more accuracy.
			if (poss[x+1,y+1,0] != 0):
				num = num + 1
				dist = dist + extrapolate(poss[x+1,y+1], vecs[x+1,y+1], ch,
						cam.campos)
			if (poss[x+1,y,0] != 0):
				num = num + 1
				dist = dist + extrapolate(poss[x+1,y], vecs[x+1,y], ch,
						cam.campos)
			poss[x,y] = (dist/num)*ch + cam.campos
			ph = vecs[x,y] - poss[x,y]
			ph = ph/norm(ph)
			vecs[x,y] = surfnorm(ch, ph)

		# Search below the starting point.
		for y in xrange(ystart+1,yb[x,1]+1):
			ch = cam.getpixdir(array([x,y])+offset)
			num = 1
			dist = extrapolate(poss[x,y-1], vecs[x,y-1], ch, cam.campos)
			# Optional to add more accuracy.
			if (poss[x+1,y-1,0] != 0):
				num = num + 1
				dist = dist + extrapolate(poss[x+1,y-1], vecs[x+1,y-1], ch,
						cam.campos)
			if (poss[x+1,y,0] != 0):
				num = num + 1
				dist = dist + extrapolate(poss[x+1,y], vecs[x+1,y], ch,
						cam.campos)
			poss[x,y] = (dist/num)*ch + cam.campos
			ph = vecs[x,y] - poss[x,y]
			ph = ph/norm(ph)
			vecs[x,y] = surfnorm(ch, ph)

	for x in xrange(ipix[0]+1,vecs.shape[0]):
		# Find a value between the bounds.
		ystart = min(max(ipix[1],yb[x,0]),yb[x,1])

		ch = cam.getpixdir(array([x,ystart])+offset)
		dist = extrapolate(poss[x-1,ystart], vecs[x-1,ystart], ch, cam.campos)
		poss[x,ystart] = dist*ch + cam.campos
		ph = vecs[x,ystart] - poss[x,ystart]
		ph = ph/norm(ph)
		vecs[x,ystart] = surfnorm(ch, ph)
		# Search above the starting point.
		for y in xrange(ystart-1,yb[x,0]-1,-1):
			ch = cam.getpixdir(array([x,y])+offset)
			num = 1
			dist = extrapolate(poss[x,y+1], vecs[x,y+1], ch, cam.campos)
			# Optional to add more accuracy.
			if (poss[x-1,y+1,0] != 0):
				num = num + 1
				dist = dist + extrapolate(poss[x-1,y+1], vecs[x-1,y+1], ch,
						cam.campos)
			if (poss[x-1,y,0] != 0):
				num = num + 1
				dist = dist + extrapolate(poss[x-1,y], vecs[x-1,y], ch,
						cam.campos)
			poss[x,y] = (dist/num)*ch + cam.campos
			ph = vecs[x,y] - poss[x,y]
			ph = ph/norm(ph)
			vecs[x,y] = surfnorm(ch, ph)

		# Search below the starting point.
		for y in xrange(ystart+1,yb[x,1]+1):
			ch = cam.getpixdir(array([x,y])+offset)
			num = 1
			dist = extrapolate(poss[x,y-1], vecs[x,y-1], ch, cam.campos)
			# Optional to add more accuracy.
			if (poss[x-1,y-1,0] != 0):
				num = num + 1
				dist = dist + extrapolate(poss[x-1,y-1], vecs[x-1,y-1], ch,
						cam.campos)
			if (poss[x-1,y,0] != 0):
				num = num + 1
				dist = dist + extrapolate(poss[x-1,y], vecs[x-1,y], ch,
						cam.campos)
			poss[x,y] = (dist/num)*ch + cam.campos
			ph = vecs[x,y] - poss[x,y]
			ph = ph/norm(ph)
			vecs[x,y] = surfnorm(ch, ph)

def surfnorm(v1, v2):
	# v1 is incident, v2 exiting.
	vnor = v2 - v1
	return vnor/norm(vnor)

def extrapolate(pos, normal, camvec, campos):
	return dot((pos-campos),normal)/dot(camvec,normal)

def centroid(image, startdot, pixwidth):
	"""Find the centre of a dot using a weighted arithmetric mean."""
	# Make searching area 3 times width of dot either side.
	# startdot comes in as tuple.
	lowx = startdot[0] - 3*pixwidth
	highx = startdot[0] + 3 *pixwidth
	lowy = startdot[1] - 3*pixwidth
	highy = startdot[1] + 3 *pixwidth
	if lowx < 0:
		lowx = 0
	if highx > pixwid-1:
		highx = pixwid-1
	if lowy < 0:
		lowy = 0
	if highy > pixhei-1:
		highy = pixhei-1
	xs = arange(lowx, highx+1)
	ys = arange(lowy, highy+1)
	# Find weighting.
	weight = zeros((len(xs), len(ys)), dtype=float32)
	for x in xs:
		for y in ys:
			# getpixel will only except pure python ints.
			pix = image.getpixel((int(x),int(y)))
			# First convert to greyscale.
			# Each channel is weighted equally.
			greypix = (pix[0] + pix[1] + pix[2])/3.
			# Then Invert.
			greypix = 255 - greypix
			weight[x-lowx,y-lowy] = greypix
	# Subtracting background.
	backgrnd = (mean(weight[0,:]) + mean(weight[weight.shape[0]-1,:]) +
		mean(weight[:,0]) + mean(weight[:,weight.shape[1]-1]))/4.
	weight = weight - backgrnd
	# Need to get rid of negatives.....
	weight = weight*(weight >= 0)
	# Normalising weight.
	weight = weight/sum(weight)
	# Calculate centres using weighted arithmetic mean.
	xav = sum(xs*transpose(weight))
	yav = sum(ys*weight)
	# Round off and save as integer. 
# Might think about allowing sub pixel accuracy.
	return array([round(xav), round(yav)], dtype=int)
