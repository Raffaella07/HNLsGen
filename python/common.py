
import numpy as np
import math

# constants 
const_GF =  1.1663787e-05 # 1/(GeV*GeV)       # GF/(hbar c)^3 from http://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
const_Vud = 0.97370                           # from (12.7)  of http://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf 
const_fpi = 130.2 * 0.001 # GeV               # from (71.14) of http://pdg.lbl.gov/2020/reviews/rpp2020-rev-pseudoscalar-meson-decay-cons.pdf
const_pi = math.pi
const_hbar = 6.582119569e-22 * 1e-03 # GeV s  # from http://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
const_c = 299792458                           # m / s

def ctau_from_gamma(gamma):
    tau_natural = 1. / gamma                  # 1/GeV
    tau = tau_natural * const_hbar            # s
    ctau = tau * const_c * 1000               # mm
    return ctau

def gamma_partial(mass,vv):
    '''
    Partial width for N->mupi, from https://arxiv.org/abs/1607.04258, assumes massless mu,pi
    '''
    gamma_partial = const_GF*const_GF / (16 * const_pi)       * np.power(mass,3) * vv * const_Vud*const_Vud * const_fpi*const_fpi  # GeV
    return gamma_partial

def gamma_total(mass,vv):
    '''
    Total width for N, from https://arxiv.org/abs/1607.04258, various approximations 
    '''
    gamma_total =   const_GF*const_GF / (96 * np.power(pi,3)) * np.power(mass,5) * vv * 10.95              # GeV
    return gamma_total

def BR_HNLmupion(mass): # vv is irrelevant, as it cancels out in the ratio
    return gamma_partial(mass=mass,vv=1.)/gamma_total(mass=mass,vv=1.)

def getCtau(mass=-99,vv=-99):
    '''
    Helper function to go from vv,m -> ctau
    '''
    return ctau_from_gamma(gamma_total(mass=mass,vv=vv))

def getVV(mass=-99.,ctau=-99.):
    '''
    Helper function to go from ctau,m -> vv
    '''
    ref_m = 1. # GeV
    ref_vv = 1. 
    ref_ctau = ctau_from_gamma(gamma_total(mass=ref_m,vv=ref_vv))

    k = ref_ctau * np.power(ref_m,5) * ref_vv

    return k/(np.power(mass, 5) * ctau)
   


##################################
# obsolete
def OLD_getVV(mass=-99.,ctau=-99.):
   '''
   Implemnts theory relation between ctau and V^2,mass 
   from extrapolation of existing W-initiated samples
   ctau = k / (m^5 * V^2)
   '''
   
   ref_m = 1. # GeV
   ref_VV = 0.022360679775*0.022360679775
   ref_ctau = 1338.961 # mm
 
   k = ref_ctau * np.power(ref_m, 5) * ref_VV

   return k/(np.power(mass, 5) * ctau)

def OLD_getCtau(mass=-99,vv=-99):
   ''' inverse of getVV function ''' 
   ref_m = 1. # GeV
   ref_VV = 0.022360679775*0.022360679775
   ref_ctau = 1338.961 # mm
   k = ref_ctau * np.power(ref_m, 5) * ref_VV

   return k / (np.power(mass, 5) * vv) 
###################################

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

