import sys
sys.path.append('../python/.')
from common import Point

# List of points to be generated
#   - N.B.: mass must be a float

points = [
Point(mass=0.5,ctau=None,vv=1e-04),
Point(mass=1.0,ctau=None,vv=1e-04),
Point(mass=1.5,ctau=None,vv=1e-04),
Point(mass=2.0,ctau=None,vv=1e-04),
Point(mass=2.5,ctau=None,vv=1e-04),
Point(mass=3.0,ctau=None,vv=1e-04),

Point(mass=1.0,ctau=None,vv=5e-03),
Point(mass=1.0,ctau=None,vv=1e-03),
Point(mass=1.0,ctau=None,vv=5e-04),
Point(mass=1.0,ctau=None,vv=1e-04),
Point(mass=1.0,ctau=None,vv=5e-05),
Point(mass=1.0,ctau=None,vv=1e-05),
Point(mass=1.0,ctau=None,vv=5e-06),
Point(mass=1.0,ctau=None,vv=1e-06),
Point(mass=1.0,ctau=None,vv=5e-07),



]


#for p in points:
#  p.stamp()




