from scipy import *
from numpy.linalg import norm
from scipy.optimize import fsolve

pix1 = [-1.,-10.]
pix2 = [10.,-10.]
pix3 = [10.,10.]
pix4 = [-10.,10.]

f = 50.
L = 1000.

" can use matrix or array.  need dot(,) when multiplying arrays "

def direc(x):
	dir = array([x[0],x[1],f])
	return dir/norm(dir)


def func(x):
	out = [x[0]**2 - 2*x[0]*x[1]*dot(p1h,p2h) + x[1]**2 - L**2]
	out.append(x[0]**2 - 2*x[0]*x[2]*dot(p1h,p4h) + x[2]**2 - L**2)
	out.append(x[2]**2 - 2*x[2]*x[1]*dot(p4h,p2h) + x[1]**2 - 2*L**2)
	return out

p1h = direc(pix1)
p2h = direc(pix2)
p3h = direc(pix3)
p4h = direc(pix4)

print p1h, p2h, p3h, p4h

" given as p1, p2, p4 "
pout = fsolve(func, [1000.,1000.,1000.])
