
import numpy as np

def getVV(mass=-99.,ctau=-99.):
   '''
   Implemnts theory relation between ctau and V^2,mass 
   Currently available from extrapolation of existing W-initiated samples
   ctau = k / (m^5 * V^2)
   '''
   
   ref_m = 1. # GeV
   ref_VV = 0.022360679775*0.022360679775
   ref_ctau = 1338.961 # mm
 
   k = ref_ctau * np.power(ref_m, 5) * ref_VV

   return k/(np.power(mass, 5) * ctau)

def getCtau(mass=-99,vv=-99):
   ''' inverse of getVV function ''' 
   ref_m = 1. # GeV
   ref_VV = 0.022360679775*0.022360679775
   ref_ctau = 1338.961 # mm
   k = ref_ctau * np.power(ref_m, 5) * ref_VV

   return k / (np.power(mass, 5) * vv) 


class Point(object):
  '''
  Class that contains information on mass,ctau,vv of a given signal point
  '''
  def __init__(self,mass,ctau=None,vv=None,isrw=False,orig_vv=None):
    self.mass = mass
    self.isrw = isrw
    if not vv: 
      self.ctau=ctau 
      self.vv=getVV(mass=self.mass, ctau=self.ctau)
    if not ctau:
      self.vv = vv
      self.ctau=getCtau(mass=self.mass, vv=self.vv)

    if self.isrw:
      self.orig_vv = orig_vv
      self.orig_ctau = getCtau(mass=self.mass,vv=orig_vv)
    else:
      self.orig_vv = self.vv 
      self.orig_ctau = self.ctau

  def stamp(self):
    attrs=[]
    for k,v in self.__dict__.items():
      attrs.append(' {}={} '.format(k,v))
    attrs=' '.join(attrs)
    print(attrs)

