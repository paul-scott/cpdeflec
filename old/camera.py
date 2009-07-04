from numpy import *
from numpy.linalg import *
from scipy.optimize import fsolve

class Camera:

	def __init__(self, w, h, pxsize, prdist, prpnt):
		self.dim = array([w, h], dtype=int16)
		self.pxsize = pxsize
		self.prdist = prdist
		self.prpnt = prpnt
		self.campos = zeros(3, dtype=float32)
		self.camtrans = array([[1,0,0],[0,1,0],[0,0,1]], dtype=float32)

	def findposition(self, pixs, sep, distguess):
		"""Assumes first pnt is origin point"""
		# Also assume that first point is global y direction from origin.
		# And that global x lies in plane of points 1-2-3.
		# Need to pass points as arrays.
		dirs = list()
		for pix in pixs:
			dirs.append(self.getpixdircam(pix))
		print(dirs)
		dists = fsolve(Camera.system3, distguess, args=(dirs,sep))

		# Find unit vectors of global coordinate system described in camera
		# array system.
		gyh = (dists[1]*dirs[1]-dists[0]*dirs[0])
		gyh = gyh/norm(gyh)
		gzh = cross((dists[2]*dirs[2]-dists[0]*dirs[0]), gyh)
		gzh = gzh/norm(gzh)
		gxh = cross(gyh, gzh)
		self.camtrans = array([gxh,gyh,gzh])
		self.campos = dot(self.camtrans,-dists[0]*dirs[0])

	def getpix(self, pnt):
		"""Convert global position into a pixel value that intercepts."""
		# Convert direction of pnt to camera coordinate system.
# Need to check multiplication.
		pdir = dot(inv(self.camtrans),(pnt-self.campos))
		# Scale to right coords.
		pdir = self.prdist*pdir/pdir[2]
		pix = (-pdir[0:2]+self.prpnt)/self.pxsize+self.dim/2.-0.5
		return array(around(pix), dtype=int)

	def getpixdir(self, pix):
		"""Rotated into system coordinates"""
		return dot(self.camtrans, self.getpixdircam(pix))

	def getpixdircam(self, pix):
		"""In camera coordinates. Pix is array of int16"""
		# Define camera ax opposite to px and ay opposite to py. az is along
		# optical axis, therefore (x,y,z)=(0,0,0) optical centre. As position
		# of pixel in array increases, image pixel value increases.
		dir = zeros(3, dtype=float32)
		# Signs used here to get the direction right. Reverse signs and we
		# get position of pixel in the array from optical centre.
		dir[2] = self.prdist
		dir[0:2] = -(self.pxsize*(pix-self.dim/2.+0.5)+self.prpnt)
		return dir/norm(dir)
	
	def getobjpix(self, width, dist):
		"""Returns estimate of number of pixels across an object is."""
		return int(round((self.prdist*width/dist)/self.pxsize))

	def system4(p, ph, l):
# Need to form traingle otherwise not enough restrictions. Use system3.
		out = [p[0]**2 - 2*p[0]*p[1]*dot(ph[0],ph[1]) + p[1]**2 - l[0]**2]
		out.append(p[0]**2 - 2*p[0]*p[3]*dot(ph[0],ph[3]) + p[3]**2 - l[3]**2)
		out.append(p[1]**2 - 2*p[1]*p[2]*dot(ph[1],ph[2]) + p[2]**2 - l[1]**2)
		out.append(p[2]**2 - 2*p[2]*p[3]*dot(ph[2],ph[3]) + p[3]**2 - l[2]**2)
		return out
	system4 = staticmethod(system4)
		
# Don't need an instance of this function for each object. Could work out how
# to do this.
#	def system3(p, ph, l):
#		out = [(p[0]**2-2*p[0]*p[1]*dot(ph[0],ph[1])+p[1]**2-l[0]**2)**2]
#		out.append((p[1]**2-2*p[1]*p[2]*dot(ph[1],ph[2])+p[2]**2-l[1]**2)**2)
#		out.append((p[0]**2-2*p[0]*p[2]*dot(ph[0],ph[2])+p[2]**2-l[2]**2)**2)
#		print('err:  ')
#		print(out)
#		print(p)
#		return out
	def system3(p, ph, l):
		out = [(norm(p[0]*ph[0]-p[1]*ph[1])-l[0])**2]
		out.append((norm(p[1]*ph[1]-p[2]*ph[2])-l[1])**2)
		out.append((norm(p[0]*ph[0]-p[2]*ph[2])-l[2])**2)
		print('err:  ')
		print(out)
		print(p)
		return out

	system3 = staticmethod(system3)
