from scipy import *
from scipy.optimize import fsolve
from scipy.optimize import fmin

# Would use camera tilted on side since mirror vert length always largest.
# Might like to take into account optical axis relative to center of array.
pixelsize = 0.00554 # Pixel size in mm.
w = 2848. # Number of CCD array pixels in width.
f = 20.52 # Lens focal length.
m = 1175. # Mirror panel width.
R = 30000. # Mirror panel radius of curvature.

def mirror(x, R, sx, sy):
	return sy - sqrt(R**2 - (x - sx)**2)

def mirrornorm(x, R, sx):
	slope = (x - sx)/sqrt(R**2 - (x - sx)**2)
	return array([cos(arctan(slope)),sin(arctan(slope))])

def solvesurface(camdist, camang, pattdist, gap, pixrange, fn, direc='fwd'):
	cam = -array([-camdist*cos(camang*pi/180.),
		camdist*sin(camang*pi/180.)])
	pcam = array([pattdist*cos(camang*pi/180.),
		pattdist*sin(camang*pi/180.)])
	
	if direc == 'rev':
		pix = pixrange[1]
		vars = [m/2., camdist]
	else:
		pix = pixrange[0]
		vars = [-m/2., camdist]

	ch = getcamray(cam, pix)
	vars = fsolve(func, vars, args=(cam, ch, R, 0, R))
	l = array([vars[0],mirror(vars[0], R, 0, R)])
	n = mirrornorm(vars[0], R, 0)
	c = vars[1]*ch
	cd = (vars[1]+gap)*ch
	ld = -cam + cd
	nd, vsm = nextnorm(c, cd, n, l, pcam)
	arrld = resize(ld, (1,2))
	arrnd = resize(nd, (1,2))
	arrvsm = [vsm]
	
	ran = list()
	if direc == 'rev':
		ran = range(pixrange[0], pixrange[1])
		ran.reverse()
	else:
		ran = range(pixrange[0]+1, pixrange[1]+1)

	for pix in ran:
		ch = getcamray(cam, pix)
		vars = fsolve(func, vars, args=(cam, ch, R, 0, R))
		l = array([vars[0],mirror(vars[0], R, 0, R)])
		n = mirrornorm(vars[0], R, 0)
		c = vars[1]*ch
		cd = nextpoint(cd, nd, c, direc)
		ld = -cam + cd
		nd, vsm = nextnorm(c, cd, n, l, pcam)
		arrld = concatenate((arrld,resize(ld, (1,2))))
		arrnd = concatenate((arrnd,resize(nd, (1,2))))
		arrvsm.append(vsm)

	optimin = fmin(ftomin, [R,0,R], args=(arrld,), maxfun=1000)

	f = open(fn, 'w')
	f.write('camdist,camang,pattdist,gap,pixrange\n')
	f.write(str(camdist)+','+str(camang)+','+str(pattdist)+','+str(gap)+
			','+str(pixrange)+'\n')
	f.write('R,sx,sy\n')
	f.write(str(optimin[0])+','+str(optimin[1])+','+str(optimin[2])+'\n')
	f.write('ldx,ldy,slopeerror,vsm\n')

	for i in range(arrld.shape[0]):
		slopeerr  = arctan(dot(arrnd[i,:],mirrornorm(arrld[i,0], optimin[0],
					optimin[1])))
		f.write(str(float(arrld[i,0]))+','+str(float(arrld[i,1]))+','+
				str(slopeerr)+','+str(arrvsm[i])+'\n')

	f.close()

def getcamray(cam, pix):
	"""Returns a unit vector in the camera ray direction. Assumes w is even"""
	u = pixelsize*(pix - (w/2. + 0.5)) # Array position in mm.
	theta = arctan(u/f)
	c = rotatez(cam, theta)
	return c/sqrt(dot(c,c))

def func(vars, cam, ch, R, sx, sy):
	return -cam + vars[1]*ch - array([vars[0],mirror(vars[0], R, sx, sy)])

def nextpoint(cdp, ndp, c, direc):
	ch = c/sqrt(dot(c,c))
	sp = rotatez(ndp, -pi/2)
	if direc == 'rev': sp = -sp
	cdm = (sp[0]*cdp[1] - cdp[0]*sp[1])/(sp[0]*ch[1] - ch[0]*sp[1])
	return cdm*ch
	
def nextnorm(c, cd, n, l, pcam):
	ch = c/sqrt(dot(c,c))
	ph = -2*dot(n,ch)*n + ch
	vh = rotatez(pcam/sqrt(dot(pcam,pcam)), pi/2.)
	pm = (-l[1] + l[0]*vh[1]/vh[0] + pcam[1] - pcam[0]*vh[1]/vh[0])/(ph[1]
			- ph[0]*vh[1]/vh[0])
	p = pm*ph
	v = p + l - pcam
	pd = p - (cd - c)
	pdh = pd/sqrt(dot(pd,pd))
	return (pdh-ch)/sqrt(dot((pdh-ch),(pdh-ch))), sign(v[0])*sqrt(dot(v,v))


def ftomin(vars, arrld):
	err = sum((arrld[:,1] - mirror(arrld[:,0], vars[0], vars[1], vars[2]))**2)
	return err

def rotatez(vec, ang):
	return resize(dot(array([[cos(ang), -sin(ang)],[sin(ang), cos(ang)]]),
			resize(vec,(2,1))),(2))

