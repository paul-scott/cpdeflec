from scipy import *

def solvesurface(camdist, camang, pattdist, gap, calls, alpha=0.000125):
	# Will need input for: camera distance, camera angle, focal length, error in focal length, mirror size, pattern distance, error in position of first point, 
	c = -array([[-camdist*sin(camang*pi/180.)],
		[camdist*cos(camang*pi/180.)]])
	cd = -array([[-(camdist+gap)*sin(camang*pi/180.)],
		[(camdist+gap)*cos(camang*pi/180.)]])
	nd = nextnorm(c, cd, array([[0],[1]]), pattdist)
	cy = c[1]

	lst = list()
	lst.append((cd, nd))
	#lst.append(cd)

	for i in range(calls-1):
		c = rotatez(c, alpha)
		c = concatenate(([c[0]*cy/c[1]],[cy]), 0)
		cd = nextpoint(cd, nd, c)
		nd = nextnorm(c, cd, array([[0],[1]]), pattdist)
		lst.append((cd, nd))
		#lst.append(cd)

	return lst

def nextpoint(cdp, ndp, c):
	cm = sqrt(dot(transpose(c),c))
	ch = c/sqrt(dot(transpose(c),c))
	sp = rotatez(ndp, -pi/2)
	cdm = (sp[0]*cdp[1] - cdp[0]*sp[1])/(sp[0]*ch[1] - ch[0]*sp[1])
	return cdm*ch
	
def nextnorm(c, cd, n, l):
	ch = c/sqrt(dot(transpose(c),c))
	ph = -2*dot(transpose(n),ch)*n + ch
	p = l*ph
	pd = p - (cd - c)
	pdh = pd/sqrt(dot(transpose(pd),pd))
	return (pdh - ch)/sqrt(dot(transpose(pdh - ch),(pdh - ch)))

def rotatez(vec, ang):
	return dot(array([[cos(ang), -sin(ang)],[sin(ang), cos(ang)]]),vec)
