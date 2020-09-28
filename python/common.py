
import numpy as np
import math
from scipy.stats import expon  

# constants 
const_GF =  1.1663787e-05 # 1/(GeV*GeV)       # GF/(hbar c)^3 from http://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
const_Vud = 0.97370                           # from (12.7)  of http://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf 
const_fpi = 130.2 * 0.001 # GeV               # from (71.14) of http://pdg.lbl.gov/2020/reviews/rpp2020-rev-pseudoscalar-meson-decay-cons.pdf
const_pi = math.pi
const_hbar = 6.582119569e-22 * 1e-03 # GeV s  # from http://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
const_c = 299792458. # m / s                  # from http://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
const_mumass = 105.6583745 * 1e-03 # GeV      # from http://pdg.lbl.gov/2020/listings/rpp2020-list-muon.pdf
const_pimass = 139.57039 * 1e-03 # GeV        # from http://pdg.lbl.gov/2020/listings/rpp2020-list-pi-plus-minus.pdf

def Lambda(a=None,b=None,c=None):
    if not a or not b or not c: 
      raise RuntimeError('mass correction factor cannot be calculated')
    return a*a + b*b + c*c - 2*a*b - 2*a*c - 2*b*c

def mass_correction_factor(mass):
    x_mu = const_mumass/mass
    x_pi = const_pimass/mass
    mcf = ( (1-x_mu*x_mu)*(1-x_mu*x_mu) - x_pi*x_pi*(1+x_mu*x_mu) ) * np.sqrt( Lambda(1,x_mu*x_mu,x_pi*x_pi))   # eq 3.5 of
    return mcf

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

    mass_correction_factor = 1. 

    return gamma_partial * mass_correction_factor

def gamma_total(mass,vv):
    '''
    Total width for N, from https://arxiv.org/abs/1607.04258, various approximations 
    '''
    gamma_total =   const_GF*const_GF / (96 * np.power(const_pi,3)) * np.power(mass,5) * vv * 10.95              # GeV
    return gamma_total

def BR_HNLmupion(mass): # vv is irrelevant, as it cancels out in the ratio
    return gamma_partial(mass=mass,vv=1.)/gamma_total(mass=mass,vv=1.)

def BR__B_D0_l_N(mass):
    '''
    Gamma(B+ -> D0 l N)/ Gamma_total_B, as extracted from figure 5 of https://link.springer.com/article/10.1007/JHEP11(2018)032
    l=electron was assumed
    '''
    if mass==0.:
      BR=0.0235  ## this equals  Gamma( B+ -> D0 l nu ) / Gamma_total_B, assuming mass N = 0 = mass nu
    elif mass==0.5:
      BR=0.0230
    elif mass==1.0:
      BR=0.0172
    elif mass==1.5:
      BR=0.0109
    elif mass==2.0:
      BR=0.00475
    elif mass==2.5:
      BR=0.00124
    elif mass==3.0:
      BR=0.0000962
    elif mass==3.5: 
      BR=0.
    elif mass==4.0:
      BR=0.
    elif mass==4.5:
      BR=0.
    elif mass==5.0:
      BR=0.
    else:
      raise RuntimeError('Do not have the BR for this mass value, please check')

    return BR

def BR__B_D0_l_nu():
    '''
    Gamma(B+ -> D0 l nu)/ Gamma_total_B
    valid for both electron and muon
    '''
    BR=0.0235  ## see above, extracted from http://pdglive.lbl.gov/BranchingRatio.action?desig=145&parCode=S041
    return BR

def BR__B_X_l_nu():
    '''
    Gamma(B+ -> X l nu) / Gamma_total_B
    valid for both electron and muon
    '''
    BR=0.1099 # from PDG  http://pdglive.lbl.gov/BranchingRatio.action?desig=220&parCode=S041
    return BR

def BR__D0_K_pi():
    '''
    Gamma(D0 -> K pi) / Gamma_total_D0
    '''
    BR=3.950E-2 # from PDG http://pdglive.lbl.gov/BranchingRatio.action?desig=1&parCode=S032
    return BR

def BR__B_l_nu():
    '''
    Gamma(B+ -> l nu) / Gamma_total_B
    l=muon was assumed
    '''
    BR = 1e-06  # 0.29 to 1.07 x 10-06 from http://pdglive.lbl.gov/BranchingRatio.action?desig=183&parCode=S041
    return BR
 
def BR__B_l_N(mass):
    '''
    Gamma(B+ -> l N) / Gamma_total_B, as extracted from figure 5 of https://link.springer.com/article/10.1007/JHEP11(2018)032
    l=electron was assumed
    '''
    if mass==0.:
      BR=5e-06  # see above 
    elif mass==0.5:
      BR=0.0000111
    elif mass==1.0:
      BR=0.0000385
    elif mass==1.5:
      BR=0.0000804
    elif mass==2.0:
      BR=0.000125
    elif mass==2.5:
      BR=0.0001601
    elif mass==3.0:
      BR=0.000177
    elif mass==3.5:
      BR=0.000161
    elif mass==4.0:
      BR=0.000118
    elif mass==4.5:
      BR=0.0000657
    elif mass==5.0:
      BR=0.0000106
    else:
      raise RuntimeError('Do not have the BR for this mass value, please check')
    return BR

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

  #def getExpMedian():
    rv = expon(scale=self.ctau) 
    self.median = rv.median()
    #return rv.mean(),rv.median()

  def stamp(self):
    attrs=[]
    for k,v in self.__dict__.items():
      attrs.append(' {}={} '.format(k,v))
    attrs=' '.join(attrs)
    print(attrs)

  def stamp_simpli(self):
    attrs=[]
    #attrs.append('{}'.format(self.mass))
    attrs.append('{}'.format(self.vv))
    attrs.append('{}'.format(self.ctau))
    attrs.append('{}'.format(self.median))
    attrs=' '.join(attrs)
    return attrs

if __name__ == "__main__":
  import matplotlib.pyplot as plt

  # test old vs new relation 
  mass = 1.0

  vvs = [5e-03, 1e-03, 5e-04, 1e-04, 5e-05, 1e-05, 5e-06, 1e-06, 5e-07]

  old_ctaus = [OLD_getCtau(mass=mass,vv=vv) for vv in vvs]
  new_ctaus = [getCtau(mass=mass,vv=vv) for vv in vvs]

  plt.figure()
  plt.plot(vvs, old_ctaus, label='old relation')
  plt.plot(vvs, new_ctaus, label='new relation')
  plt.xlabel('$|V|^2$')
  plt.ylabel('ctau [mm]')
  plt.grid(True) # means principal and minor 
  plt.yscale('log')
  plt.xscale('log')
  plt.legend()
  plt.savefig('./comparison_ctau.pdf')
  #plt.show()

  # test mass correction factor  
  masses = [0.5+i*0.1 for i in range(0,30) ] # GeV
  mass_correction_factors = [mass_correction_factor(mass=mass) for mass in masses]  
  plt.figure()
  plt.plot(masses, mass_correction_factors)
  plt.xlabel('HNL mass [GeV]')
  plt.ylabel('correction factor to partial width')
  plt.savefig('./mass_correction_factor.pdf')

  #### BR correction factor
  # test (BR(B->DmuN) + BR(B->Nmu)) / (  BR(B->Dmunu) + BR(B->munu) )
  masses = [0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
  
  average_Ntonu_BR_corrections = [ ( BR__B_l_N(mass=mass) + BR__B_D0_l_N(mass=mass) ) / ( BR__B_l_nu() + BR__B_D0_l_nu() ) for mass in masses]
  
  D0muN_Ntonu_BR_corrections = [ BR__B_D0_l_N(mass=mass) / BR__B_D0_l_nu() for mass in masses  ]
  muN_Ntonu_BR_corrections =   [ BR__B_l_N(mass=mass)    / BR__B_l_nu()    for mass in masses  ]

  plt.figure()
  plt.plot(masses, D0muN_Ntonu_BR_corrections, label='channel: $B^{+} -> \\bar{D}^{0} \mu^{+} N$')
  plt.plot(masses, muN_Ntonu_BR_corrections, label='channel: $B^{+} -> \mu^{+} N$')
  plt.plot(masses, average_Ntonu_BR_corrections, label='Average')
  plt.xlabel('mass (GeV)')
  plt.ylabel(' N -> $\\nu$ correction factor')
  plt.yscale('log')
  plt.grid(True) # means principal and minor 
  plt.legend()
  plt.ylim(top=10000)
  plt.savefig('./Ntonu_BR_corrections.pdf')



  # Calculate neutrino numbers
  nB_all = 9.6E9 # after purity...
  # fractions from http://cds.cern.ch/record/2704495/files/DP2019_043.pdf

  # Channel 1 
  fB = 0.4
  BR__B_munu = BR__B_l_nu()
  BR__B_Xlnu = BR__B_X_l_nu()
  n_nu = nB_all * fB * BR__B_munu / BR__B_Xlnu 
  print('Channel 1: expected N(nu)={:2.1e}'.format(n_nu))

  # Channel 2
  fB = 0.4 
  BR__B_D0lnu = BR__B_D0_l_nu()
  BR__B_Xlnu = BR__B_X_l_nu()
  BR__D0_Kpi = 0.0395 # http://pdglive.lbl.gov/BranchingRatio.action?desig=1&parCode=S032
  n_nu = nB_all * fB  * BR__B_D0lnu / BR__B_Xlnu * BR__D0_Kpi
  print('Channel 2: expected N(nu)={:2.1e}'.format(n_nu))

  # Channel 3
  fB0 = 0.4 
  BR__B0_Dmunu = 0.0231 # http://pdglive.lbl.gov/BranchingRatio.action?desig=37&parCode=S042
  BR__B0_Xlnu = 0.1033  # http://pdglive.lbl.gov/BranchingRatio.action?desig=94&parCode=S042
  BR__D_Kpipi = 0.0938  # Gamma49/Gamma http://pdglive.lbl.gov/BranchingRatio.action?desig=1&parCode=S031&expand=true
  n_nu = nB_all * fB0 * BR__B0_Dmunu / BR__B0_Xlnu * BR__D_Kpipi
  print('Channel 3: expected N(nu)={:2.1e}'.format(n_nu))

  # Channel 4 
  fB0S = 0.1
  BR__B0S_DSmunuX = 0.081 # NOTE: not fully exclusive  http://pdglive.lbl.gov/BranchingRatio.action?desig=4&parCode=S086# not fully exclusive
  BR__B0S_Xmunu = 0.102   # http://pdglive.lbl.gov/BranchingRatio.action?desig=99&parCode=S086
  BR__DS_KKpi = 0.0539    #http://pdglive.lbl.gov/BranchingRatio.action?desig=40&parCode=S034
  n_nu = nB_all * fB0S * BR__B0S_DSmunuX / BR__B0S_Xmunu * BR__DS_KKpi 
  print('Channel 4: expected N(nu)={:2.1e}'.format(n_nu))
  
  

  # table for benchmark points
  masses = [1,2,3]
  benchmark_vvs = [3E-06, 3E-06, 1E-05]

  for i,mass in enumerate(masses):
    print('\nMass={}GeV').format(mass)
    p = Point(mass=mass,vv=benchmark_vvs[i])
    print('{} '.format(1.0) + p.stamp_simpli())
    for mult in np.logspace(-1,3,50):
      p = Point(mass=mass,vv=benchmark_vvs[i]*mult)
      #print(mult)
      print('{} '.format(mult) + p.stamp_simpli())

#  points[1] = [
#    Point(mass=1,vv=3E-06),
#    Point(mass=1,vv=3E-06*2),
#    Point(mass=1,vv=3E-06*5),
#    Point(mass=1,vv=3E-06*10),
#    Point(mass=1,vv=3E-06*50),
#    Point(mass=1,vv=3E-06*100),
#  ]
#  points[2] = [
#    Point(mass=2,vv=3E-06),
#    Point(mass=2,vv=3E-06*2),
#    Point(mass=2,vv=3E-06*5),
#    Point(mass=2,vv=3E-06*10),
#    Point(mass=2,vv=3E-06*50),
#    Point(mass=2,vv=3E-06*100),
#  ]
#  points[3] = [
#    Point(mass=3,vv=1E-05),
#    Point(mass=3,vv=1E-05*2),
#    Point(mass=3,vv=1E-05*5),
#    Point(mass=3,vv=1E-05*10),
#    Point(mass=3,vv=1E-05*50),
#    Point(mass=3,vv=1E-05*100),
#  ]
#
#  for mass in [1,2,3]:
#    print('\nMass={}GeV'.format(mass))
#    for p in points[mass]:
#      p.stamp_nofield()
#
#
#
#  from scipy.stats import expon
#
