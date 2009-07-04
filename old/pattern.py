import colorsys
from numpy import *

class Pattern:

	def __init__(self, segsize, pos, trans, relat):
		"""Note that might want to get pattern coords to work out trans here"""
		self.segsize = segsize
		self.pos = pos
		self.trans = trans
		# Polynomial coefficients.
		self.relat = relat
# Think about adding another parameter where read in colour of white area
# Make sure position takes into account both rotations and that transform
# accounts for any over or under rotation of pattern(maybe can't use simple
# matrix??).

	def getvector(self, coord):
		# Is now designed to be called on whole array of points.
		# Need to inner here because handles multidimensional case better.
		# Inner uses sum product over last indices of each.
		# output has dimensions a.shape[:-1] + b.shape[:-1]
		# For dot sum product over last and second last.
		# Could use tensordot(coords, self.trans, axes=(-1,-1)).
		return inner(coord, self.trans) + self.pos

	def getdist(self, colour, pdist=0.0):
		hue = self.gethue(colour)
		shift = self.shiftrela(hue)

		pshift = mod(pdist, self.segsize)
		if shift < (pshift-0.5*self.segsize):
			return (pdist - pshift + self.segsize + shift)
		elif shift > (pshift+0.5*self.segsize):
			return (pdist - pshift - self.segsize + shift)
		else:
			return (pdist - pshift + shift)

	def gethue(self, colour):
		h, l, s = colorsys.rgb_to_hls(colour[0]/255., colour[1]/255.,
			colour[2]/255.)
		return h
# Need to correct hue?
	
	def shiftrela(self, hue):
		return self.relat(hue)*self.segsize
