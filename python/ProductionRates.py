'''
This script collects the particle definition and creates the ProductionRates class, that lists
the production rates of interest, with respect to the HNL mass and the mixing angle
'''

from objects import Particle, Decay
import math


## CONSTANTS ##

# masses [in GeV] 
# mesons
m_B_pdg       = 5.27934 
m_B0_pdg      = 5.27965
m_B_sub_c_pdg = 6.2749
m_B_sub_s_pdg = 5.36688
m_D_pdg       = 1.86965  
m_D0_pdg      = 1.86483
m_D_sub_s_pdg = 1.96834
m_D0star_pdg  = 2.00685 # check
m_Dstar_pdg   = 2.001021
m_D_sub_sstar_pdg = 2.1122
m_K_pdg       = 0.493677
m_K0_pdg      = 0.497611
m_Kstar_pdg   = 0.89166 
m_K0star_pdg  = 0.89555 #0.824 #1.425
m_rho0_pdg    = 0.77526 #0.769 # 0.770 # 0.77526 # make sure this is the correct meson to consider # in GeV
m_rho_pdg     = 0.77511
m_pi0_pdg     = 0.1349768
m_pi_pdg      = 0.13957039

# leptons
m_el_pdg  = 0.510999 * 1e-3 
m_mu_pdg  = 0.105658 
m_tau_pdg = 1.77686


# meson decay constant --> to be checked
# in [GeV]
dc_B       = 0.1871 # Table8 Bondarenko et al. 0.19   # 176 MeV according to pdg
dc_B_sub_c = 0.434 # Table8 Bondarenko et al. # 0.48
dc_D       = 0.212 # Table8 Bondarenko et al. # 0.2226   
dc_D_sub_s = 0.249 # Table8 Bondarenko et al

# lifetime in [GeV-1]
sToGeVconv       = 1.5198e24
lifetime_B       = 2.49e12 
lifetime_B0      = 1519e-15  * sToGeVconv 
lifetime_B_sub_c = 0.51e-12  * sToGeVconv
lifetime_B_sub_s = 1.527e-12 * sToGeVconv
lifetime_D       = 1.04e-12  * sToGeVconv
lifetime_D_sub_s = 5.04e-13  * sToGeVconv
lifetime_D0      = 4.101e-13 * sToGeVconv

# B-fractions
fraction_B       = 0.4
fraction_B0      = 0.4
fraction_B_sub_s = 0.1
fraction_B_sub_c = 0.001

# CKM matrix elements
Vud_pdg = 0.97417 #0.97420
Vus_pdg = 0.2248
Vub_pdg = 0.00409
Vcb_pdg = 40.5e-3
Vcd_pdg = 0.220
Vcs_pdg = 0.995


## PARTICLES ##

# beauty mesons
B_meson           = Particle('B_meson'          , 'meson', m_B_pdg          , dc_B      , lifetime_B               , fraction=fraction_B)
B0_meson          = Particle('B0_meson'         , 'meson', m_B0_pdg                     , lifetime=lifetime_B0     , fraction=fraction_B0)
B_sub_c_meson     = Particle('B_sub_c_meson'    , 'meson', m_B_sub_c_pdg    , dc_B_sub_c, lifetime_B_sub_c         , fraction=fraction_B_sub_c)
B_sub_s_meson     = Particle('B_sub_s_meson'    , 'meson', m_B_sub_s_pdg                , lifetime=lifetime_B_sub_s, fraction=fraction_B_sub_s)

# charm mesons
D_meson           = Particle('D_meson'          , 'meson', m_D_pdg          , dc_D      , lifetime_D)
D0_meson          = Particle('D0_meson'         , 'meson', m_D0_pdg         , dc_D      , lifetime_D0)
D_sub_s_meson     = Particle('D_sub_s_meson'    , 'meson', m_D_sub_s_pdg    , dc_D_sub_s, lifetime_D_sub_s)
D0star_meson      = Particle('D0star_meson'     , 'meson', m_D0star_pdg)
Dstar_meson       = Particle('Dstar_meson'      , 'meson', m_Dstar_pdg)
D_sub_sstar_meson = Particle('D_sub_sstar_meson', 'meson', m_D_sub_sstar_pdg)

# strange mesons
K_meson           = Particle('K_meson'          , 'meson', m_K_pdg)
K0_meson          = Particle('K0_meson'         , 'meson', m_K0_pdg)
Kstar_meson       = Particle('Kstar_meson'      , 'meson', m_Kstar_pdg)
K0star_meson      = Particle('K0star_meson'     , 'meson', m_K0star_pdg)

# light mesons
rho0_meson        = Particle('rho0_meson'       , 'meson', m_rho0_pdg)
rho_meson         = Particle('rho_meson'        , 'meson', m_rho_pdg)
pi0_meson         = Particle('pi0_meson'        , 'meson', m_pi0_pdg) 
pi_meson          = Particle('pi_meson'         , 'meson', m_pi_pdg) 

# leptons
el                = Particle('el'               , 'lepton', m_el_pdg) 
mu                = Particle('mu'               , 'lepton', m_mu_pdg) 
tau               = Particle('tau'              , 'lepton', m_tau_pdg) 


## DECAYS ##

class ProductionRates(object):
  def __init__(self, mass, mixing_angle_square): # add model label? 
    self.mass = mass
    self.mixing_angle_square = mixing_angle_square

    # define the HNL
    hnl = Particle('hnl', 'lepton', self.mass)

    # define the HNL
    nu = Particle('nu', 'lepton', 0)
    
    # get the model
    V_el_square = self.mixing_angle_square
    V_mu_square = self.mixing_angle_square
    V_tau_square = self.mixing_angle_square

    # list of the decays of interest
    # leptonic
    B_to_eHNL  = Decay(B_meson,       el, hnl, V_el_square, Vub_pdg, 'leptonic').decay_rate
    Bc_to_eHNL = Decay(B_sub_c_meson, el, hnl, V_el_square, Vcb_pdg, 'leptonic').decay_rate 
    D_to_eHNL  = Decay(D_meson,       el, hnl, V_el_square, Vcd_pdg, 'leptonic').decay_rate
    Ds_to_eHNL = Decay(D_sub_s_meson, el, hnl, V_el_square, Vcs_pdg, 'leptonic').decay_rate
    
    B_to_uHNL  = Decay(B_meson,       mu, hnl, V_mu_square, Vub_pdg, 'leptonic').decay_rate if self.mass < 5.1 else 0
    Bc_to_uHNL = Decay(B_sub_c_meson, mu, hnl, V_mu_square, Vcb_pdg, 'leptonic').decay_rate 
    D_to_uHNL  = Decay(D_meson,       mu, hnl, V_mu_square, Vcd_pdg, 'leptonic').decay_rate if self.mass < 1.7 else 0
    Ds_to_uHNL = Decay(D_sub_s_meson, mu, hnl, V_mu_square, Vcs_pdg, 'leptonic').decay_rate if self.mass < 1.7 else 0

    B_to_tHNL  = Decay(B_meson,       tau, hnl, V_tau_square, Vub_pdg, 'leptonic').decay_rate if self.mass < 0.15 else 0
    Bc_to_tHNL = Decay(B_sub_c_meson, tau, hnl, V_tau_square, Vcb_pdg, 'leptonic').decay_rate if self.mass < 4.5 else 0 
    D_to_tHNL  = Decay(D_meson,       tau, hnl, V_tau_square, Vcd_pdg, 'leptonic').decay_rate if self.mass < 0.15 else 0
    Ds_to_tHNL = Decay(D_sub_s_meson, tau, hnl, V_tau_square, Vcs_pdg, 'leptonic').decay_rate if self.mass < 0.2 else 0

    # semileptonic into pseudoscalar meson
    B_to_D0eHNL   = Decay(B_meson      , [D0_meson, el]     , hnl, V_el_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.4 else 0
    B_to_pi0eHNL  = Decay(B_meson      , [pi0_meson, el]    , hnl, V_el_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi').decay_rate if self.mass < 5.1 else 0
    B0_to_pieHNL  = Decay(B0_meson     , [pi_meson, el]     , hnl, V_el_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi').decay_rate if self.mass < 5.1 else 0
    B0_to_DeHNL   = Decay(B0_meson     , [D_meson, el]      , hnl, V_el_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.4 else 0
    Bs_to_KeHNL   = Decay(B_sub_s_meson, [K_meson, el]      , hnl, V_el_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='Bs_to_K').decay_rate if self.mass < 4.8 else 0
    Bs_to_DseHNL  = Decay(B_sub_s_meson, [D_sub_s_meson, el], hnl, V_el_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.4 else 0
    D_to_K0eHNL   = Decay(D_meson      , [K0_meson, el]     , hnl, V_el_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K' ).decay_rate  if self.mass < 1.3 else 0
    D_to_pi0eHNL  = Decay(D_meson      , [pi0_meson, el]    , hnl, V_el_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi').decay_rate if self.mass < 1.7 else 0
    D0_to_pieHNL  = Decay(D0_meson     , [pi_meson, el]     , hnl, V_el_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi').decay_rate if self.mass < 1.7 else 0
    D0_to_KeHNL   = Decay(D0_meson     , [K_meson, el]      , hnl, V_el_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K' ).decay_rate  if self.mass < 1.3 else 0

    B_to_D0uHNL   = Decay(B_meson      , [D0_meson, mu]     , hnl, V_mu_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.3 else 0
    B_to_pi0uHNL  = Decay(B_meson      , [pi0_meson, mu]    , hnl, V_mu_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi').decay_rate if self.mass < 5.0 else 0
    B0_to_piuHNL  = Decay(B0_meson     , [pi_meson, mu]     , hnl, V_mu_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi').decay_rate if self.mass < 5.0 else 0
    B0_to_DuHNL   = Decay(B0_meson     , [D_meson, mu]      , hnl, V_mu_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.3 else 0
    Bs_to_KuHNL   = Decay(B_sub_s_meson, [K_meson, mu]      , hnl, V_mu_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='Bs_to_K').decay_rate if self.mass < 4.7 else 0
    Bs_to_DsuHNL  = Decay(B_sub_s_meson, [D_sub_s_meson, mu], hnl, V_mu_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.3 else 0
    D_to_K0uHNL   = Decay(D_meson      , [K0_meson, mu]     , hnl, V_mu_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K' ).decay_rate  if self.mass < 1.2 else 0
    D_to_pi0uHNL  = Decay(D_meson      , [pi0_meson, mu]    , hnl, V_mu_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi').decay_rate if self.mass < 1.6 else 0
    D0_to_piuHNL  = Decay(D0_meson     , [pi_meson, mu]     , hnl, V_mu_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi').decay_rate if self.mass < 1.6 else 0
    D0_to_KuHNL   = Decay(D0_meson     , [K_meson, mu]      , hnl, V_mu_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K' ).decay_rate  if self.mass < 1.2 else 0

    B_to_D0tHNL   = Decay(B_meson      , [D0_meson, tau]     , hnl, V_tau_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 1.6 else 0
    B_to_pi0tHNL  = Decay(B_meson      , [pi0_meson, tau]    , hnl, V_tau_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi').decay_rate if self.mass < 3.4 else 0
    B0_to_pitHNL  = Decay(B0_meson     , [pi_meson, tau]     , hnl, V_tau_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi').decay_rate if self.mass < 3.4 else 0
    B0_to_DtHNL   = Decay(B0_meson     , [D_meson, tau]      , hnl, V_tau_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 1.6 else 0
    Bs_to_KtHNL   = Decay(B_sub_s_meson, [K_meson, tau]      , hnl, V_tau_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='Bs_to_K').decay_rate if self.mass < 3.1 else 0
    Bs_to_DstHNL  = Decay(B_sub_s_meson, [D_sub_s_meson, tau], hnl, V_tau_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 1.6 else 0
    #D_to_K0tHNL   = Decay(D_meson      , [K0_meson, tau]     , hnl, V_tau_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K').decay_rate  if self.mass < 4 else 0
    #D_to_pi0tHNL  = Decay(D_meson      , [pi0_meson, tau]    , hnl, V_tau_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi').decay_rate if self.mass < 0.1 else 0
    #D0_to_pitHNL  = Decay(D0_meson     , [pi_meson, tau]     , hnl, V_tau_square, Vcd_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_pi').decay_rate if self.mass < 1.7 else 0
    #D0_to_KtHNL   = Decay(D0_meson     , [K_meson, tau]      , hnl, V_tau_square, Vcs_pdg, 'semileptonic_pseudoscalar', formFactorLabel='D_to_K').decay_rate  if self.mass < 1.4 else 0
    
    
    # semileptonic into vector meson
    B_to_rho0eHNL    = Decay(B_meson      , [rho0_meson, el]       , hnl, V_el_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ).decay_rate     if self.mass < 4.5 else 0
    B_to_D0stareHNL  = Decay(B_meson      , [D0star_meson, el]     , hnl, V_el_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ).decay_rate   if self.mass < 3.2 else 0
    B0_to_DstareHNL  = Decay(B0_meson     , [Dstar_meson, el]      , hnl, V_el_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ).decay_rate   if self.mass < 3.2 else 0
    B0_to_rhoeHNL    = Decay(B0_meson     , [rho_meson, el]        , hnl, V_el_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ).decay_rate     if self.mass < 4.5 else 0
    Bs_to_DsstareHNL = Decay(B_sub_s_meson, [D_sub_sstar_meson, el], hnl, V_el_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Dsstar').decay_rate if self.mass < 3.2 else 0 
    Bs_to_KstareHNL  = Decay(B_sub_s_meson, [Kstar_meson, el]      , hnl, V_el_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Kstar' ).decay_rate  if self.mass < 4.4 else 0
    D_to_K0stareHNL  = Decay(D_meson      , [K0star_meson, el]     , hnl, V_el_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar'  ).decay_rate   if self.mass < 0.98 else 0
    D0_to_KstareHNL  = Decay(D0_meson     , [Kstar_meson, el]      , hnl, V_el_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar'  ).decay_rate   if self.mass < 0.98 else 0
    
    B_to_rho0uHNL    = Decay(B_meson      , [rho0_meson, mu]       , hnl, V_mu_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ).decay_rate     if self.mass < 4.4 else 0
    B_to_D0staruHNL  = Decay(B_meson      , [D0star_meson, mu]     , hnl, V_mu_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ).decay_rate   if self.mass < 3.1 else 0
    B0_to_DstaruHNL  = Decay(B0_meson     , [Dstar_meson, mu]      , hnl, V_mu_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ).decay_rate   if self.mass < 3.1 else 0
    B0_to_rhouHNL    = Decay(B0_meson     , [rho_meson, mu]        , hnl, V_mu_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ).decay_rate     if self.mass < 4.4 else 0
    Bs_to_DsstaruHNL = Decay(B_sub_s_meson, [D_sub_sstar_meson, mu], hnl, V_mu_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Dsstar').decay_rate if self.mass < 3.1 else 0 
    Bs_to_KstaruHNL  = Decay(B_sub_s_meson, [Kstar_meson, mu]      , hnl, V_mu_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Kstar' ).decay_rate  if self.mass < 4.3 else 0
    D_to_K0staruHNL  = Decay(D_meson      , [K0star_meson, mu]     , hnl, V_mu_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar'  ).decay_rate   if self.mass < 0.9 else 0
    D0_to_KstaruHNL  = Decay(D0_meson     , [Kstar_meson, mu]      , hnl, V_mu_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar'  ).decay_rate   if self.mass < 0.9 else 0

    B_to_rho0tHNL    = Decay(B_meson      , [rho0_meson, tau]       , hnl, V_tau_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ).decay_rate     if self.mass < 2.7 else 0
    B_to_D0startHNL  = Decay(B_meson      , [D0star_meson, tau]     , hnl, V_tau_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ).decay_rate   if self.mass < 1.5 else 0
    B0_to_DstartHNL  = Decay(B0_meson     , [Dstar_meson, tau]      , hnl, V_tau_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ).decay_rate   if self.mass < 1.5 else 0
    B0_to_rhotHNL    = Decay(B0_meson     , [rho_meson, tau]        , hnl, V_tau_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ).decay_rate     if self.mass < 2.7 else 0
    Bs_to_DsstartHNL = Decay(B_sub_s_meson, [D_sub_sstar_meson, tau], hnl, V_tau_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Dsstar').decay_rate if self.mass < 1.4 else 0 
    Bs_to_KstartHNL  = Decay(B_sub_s_meson, [Kstar_meson, tau]      , hnl, V_tau_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Kstar' ).decay_rate  if self.mass < 2.7 else 0
    #D_to_K0startHNL  = Decay(D_meson      , [K0star_meson, tau]     , hnl, V_tau_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar').decay_rate   if self.mass < 0.98 else 0
    #D0_to_KstartHNL  = Decay(D0_meson     , [Kstar_meson, tau]      , hnl, V_tau_square, Vcs_pdg, 'semileptonic_vector', formFactorLabel='D_to_Kstar').decay_rate   if self.mass < 0.98 else 0
    
    # SM decays
    B_to_unu  = Decay(B_meson      ,  mu, nu, V_mu_square, Vub_pdg, 'leptonic').decay_rate if self.mass < 5.1 else 0
    Bc_to_unu = Decay(B_sub_c_meson, mu, nu, V_mu_square, Vcb_pdg, 'leptonic').decay_rate 

    B_to_D0unu   = Decay(B_meson      , [D0_meson, mu]     , nu, V_mu_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.3 else 0
    B_to_pi0unu  = Decay(B_meson      , [pi0_meson, mu]    , nu, V_mu_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi').decay_rate if self.mass < 5.0 else 0
    B0_to_piunu  = Decay(B0_meson     , [pi_meson, mu]     , nu, V_mu_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_pi').decay_rate if self.mass < 5.0 else 0
    B0_to_Dunu   = Decay(B0_meson     , [D_meson, mu]      , nu, V_mu_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.3 else 0
    Bs_to_Kunu   = Decay(B_sub_s_meson, [K_meson, mu]      , nu, V_mu_square, Vub_pdg, 'semileptonic_pseudoscalar', formFactorLabel='Bs_to_K').decay_rate if self.mass < 4.7 else 0
    Bs_to_Dsunu  = Decay(B_sub_s_meson, [D_sub_s_meson, mu], nu, V_mu_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.3 else 0
    
    B_to_rho0unu    = Decay(B_meson      , [rho0_meson, mu]       , nu, V_mu_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ).decay_rate     if self.mass < 4.4 else 0
    B_to_D0starunu  = Decay(B_meson      , [D0star_meson, mu]     , nu, V_mu_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ).decay_rate   if self.mass < 3.1 else 0
    B0_to_Dstarunu  = Decay(B0_meson     , [Dstar_meson, mu]      , nu, V_mu_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='B_to_Dstar'  ).decay_rate   if self.mass < 3.1 else 0
    B0_to_rhounu    = Decay(B0_meson     , [rho_meson, mu]        , nu, V_mu_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='B_to_rho'    ).decay_rate     if self.mass < 4.4 else 0
    Bs_to_Dsstarunu = Decay(B_sub_s_meson, [D_sub_sstar_meson, mu], nu, V_mu_square, Vcb_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Dsstar').decay_rate if self.mass < 3.1 else 0 
    Bs_to_Kstarunu  = Decay(B_sub_s_meson, [Kstar_meson, mu]      , nu, V_mu_square, Vub_pdg, 'semileptonic_vector', formFactorLabel='Bs_to_Kstar' ).decay_rate  if self.mass < 4.3 else 0
    
    B_to_enu  = Decay(B_meson,       el, nu, V_el_square, Vub_pdg, 'leptonic').decay_rate if self.mass < 5.1 else 0
    B_to_D0enu   = Decay(B_meson      , [D0_meson, el]     , nu, V_el_square, Vcb_pdg, 'semileptonic_pseudoscalar', formFactorLabel='B_to_D' ).decay_rate  if self.mass < 3.3 else 0

    
    # production rates

    # B-channel
    self.Gamma_Be = B_to_eHNL          + Bc_to_eHNL      + B_to_D0eHNL \
                    + B_to_pi0eHNL     + B0_to_pieHNL    + B0_to_DeHNL  \
                    + Bs_to_KeHNL      + Bs_to_DseHNL    + B_to_rho0eHNL \
                    + B_to_D0stareHNL  + B0_to_DstareHNL + B0_to_rhoeHNL \
                    + Bs_to_DsstareHNL + Bs_to_KstareHNL

    self.Gamma_Bu = B_to_uHNL          + Bc_to_uHNL      + B_to_D0uHNL \
                    + B_to_pi0uHNL     + B0_to_piuHNL    + B0_to_DuHNL  \
                    + Bs_to_KuHNL      + Bs_to_DsuHNL    + B_to_rho0uHNL \
                    + B_to_D0staruHNL  + B0_to_DstaruHNL + B0_to_rhouHNL \
                    + Bs_to_DsstaruHNL + Bs_to_KstaruHNL
    
    #self.BR_Bu =    B_to_uHNL*B_meson.lifetime        + Bc_to_uHNL*B_sub_c_meson.lifetime      + B_to_D0uHNL*B_meson.lifetime \
    #                + B_to_pi0uHNL*B_meson.lifetime      + B0_to_piuHNL*B0_meson.lifetime    + B0_to_DuHNL*B0_meson.lifetime  \
    #                + Bs_to_KuHNL*B_sub_s_meson.lifetime      + Bs_to_DsuHNL*B_sub_s_meson.lifetime    + B_to_rho0uHNL*B_meson.lifetime \
    #                + B_to_D0staruHNL*B_meson.lifetime  + B0_to_DstaruHNL*B0_meson.lifetime + B0_to_rhouHNL*B0_meson.lifetime \
    #                + Bs_to_DsstaruHNL*B_sub_s_meson.lifetime + Bs_to_KstaruHNL*B_sub_s_meson.lifetime

    self.BR_Bu_HNL =   B_meson.lifetime  * B_meson.fraction  * (B_to_uHNL + B_to_D0uHNL + B_to_pi0uHNL + B_to_rho0uHNL + B_to_D0staruHNL) \
                 + B0_meson.lifetime * B0_meson.fraction * (B0_to_piuHNL + B0_to_DuHNL + B0_to_DstaruHNL + B0_to_rhouHNL) \
                 + B_sub_s_meson.lifetime * B_sub_s_meson.fraction * (Bs_to_KuHNL + Bs_to_DsuHNL + Bs_to_DsstaruHNL + Bs_to_KstaruHNL) \
                 + B_sub_c_meson.lifetime * B_sub_c_meson.fraction * Bc_to_uHNL

    #self.BR_Bu_HNL_small =  B_meson.lifetime  * B_meson.fraction  * (B_to_uHNL + B_to_D0uHNL) 
    self.BR_Bu_HNL_small =  B_meson.lifetime  *  (B_to_uHNL + B_to_D0uHNL) 
    self.BR_Be_HNL_small =  B_meson.lifetime  *  (B_to_eHNL + B_to_D0eHNL) 
    self.BR_Be_HNL_verysmall =  B_meson.lifetime  * B_to_eHNL 
    self.BR_Bu_HNL_medium =   B_meson.lifetime  * (B_to_uHNL + B_to_D0uHNL + B_to_pi0uHNL + B_to_rho0uHNL + B_to_D0staruHNL) 
    self.BR_Be_HNL_medium =   B_meson.lifetime  * (B_to_eHNL + B_to_D0eHNL + B_to_pi0eHNL + B_to_rho0eHNL + B_to_D0stareHNL) 
    self.BR_Bt_HNL_medium =   B_meson.lifetime  * (B_to_tHNL + B_to_D0tHNL + B_to_pi0tHNL + B_to_rho0tHNL + B_to_D0startHNL) 
    
    self.BR_Bu_nu =   B_meson.lifetime  * B_meson.fraction  * (B_to_unu + B_to_D0unu + B_to_pi0unu + B_to_rho0unu + B_to_D0starunu) \
                 + B0_meson.lifetime * B0_meson.fraction * (B0_to_piunu + B0_to_Dunu + B0_to_Dstarunu + B0_to_rhounu) \
                 + B_sub_s_meson.lifetime * B_sub_s_meson.fraction * (Bs_to_Kunu + Bs_to_Dsunu + Bs_to_Dsstarunu + Bs_to_Kstarunu) \
                 + B_sub_c_meson.lifetime * B_sub_c_meson.fraction * Bc_to_unu

    #self.BR_Bu_nu_small =   B_meson.lifetime  * B_meson.fraction  * (B_to_unu + B_to_D0unu) 
    self.BR_Bu_nu_small =   B_meson.lifetime * (B_to_unu + B_to_D0unu) 
    self.BR_Be_nu_small =   B_meson.lifetime * (B_to_enu + B_to_D0enu) 

    
    self.BR_Bl_test =   B_meson.lifetime  * B_meson.fraction  * (B_to_unu + B_to_D0unu + B_to_pi0unu + B_to_rho0unu + B_to_D0starunu) / (self.mixing_angle_square * 0.1099) \
                 + B0_meson.lifetime * B0_meson.fraction * (B0_to_piunu + B0_to_Dunu + B0_to_Dstarunu + B0_to_rhounu) / (self.mixing_angle_square * 0.1033) \
                 + B_sub_s_meson.lifetime * B_sub_s_meson.fraction * (Bs_to_Kunu + Bs_to_Dsunu + Bs_to_Dsstarunu + Bs_to_Kstarunu) / (self.mixing_angle_square * 0.096) #\
                 #+ B_sub_c_meson.lifetime * B_sub_c_meson.fraction * Bc_to_unu / (self.mixing_angle_square * 0.1099)


    self.Gamma_Bt = B_to_tHNL          + Bc_to_tHNL      + B_to_D0tHNL \
                    + B_to_pi0tHNL     + B0_to_pitHNL    + B0_to_DtHNL  \
                    + Bs_to_KtHNL      + Bs_to_DstHNL    + B_to_rho0tHNL \
                    + B_to_D0startHNL  + B0_to_DstartHNL + B0_to_rhotHNL \
                    + Bs_to_DsstartHNL + Bs_to_KstartHNL

    self.Gamma_Beut = self.Gamma_Be + self.Gamma_Bu + self.Gamma_Bt
    
    self.Gamma_Be_leptonic = B_to_eHNL + Bc_to_eHNL
    self.Gamma_Bu_leptonic = B_to_uHNL + Bc_to_uHNL
    self.Gamma_Bt_leptonic = B_to_tHNL + Bc_to_tHNL
    self.Gamma_Beut_leptonic = self.Gamma_Be_leptonic + self.Gamma_Bu_leptonic + self.Gamma_Bt_leptonic

    self.Gamma_Be_semileptonic_pseudoscalar = B_to_D0eHNL    + B_to_pi0eHNL \
                                              + B0_to_pieHNL + B0_to_DeHNL  \
                                              + Bs_to_KeHNL  + Bs_to_DseHNL 

    self.Gamma_Bu_semileptonic_pseudoscalar = B_to_D0uHNL    + B_to_pi0uHNL \
                                              + B0_to_piuHNL + B0_to_DuHNL  \
                                              + Bs_to_KuHNL  + Bs_to_DsuHNL 

    self.Gamma_Be_semileptonic_vector = B_to_rho0eHNL      + B_to_D0stareHNL \
                                        + B0_to_DstareHNL  + B0_to_rhoeHNL \
                                        + Bs_to_DsstareHNL + Bs_to_KstareHNL

    self.Gamma_Bu_semileptonic_vector = B_to_rho0uHNL      + B_to_D0staruHNL \
                                        + B0_to_DstaruHNL  + B0_to_rhouHNL \
                                        + Bs_to_DsstaruHNL + Bs_to_KstaruHNL




    # D-channel
    self.Gamma_De_leptonic = D_to_eHNL + Ds_to_eHNL
    self.Gamma_Du_leptonic = D_to_uHNL + Ds_to_uHNL
    self.Gamma_De_semileptonic_pseudoscalar = D_to_K0eHNL + D_to_pi0eHNL + D0_to_pieHNL + D0_to_KeHNL
    self.Gamma_Du_semileptonic_pseudoscalar = D_to_K0uHNL + D_to_pi0uHNL + D0_to_piuHNL + D0_to_KuHNL
    self.Gamma_De_semileptonic_vector = D_to_K0stareHNL + D0_to_KstareHNL 
    self.Gamma_Du_semileptonic_vector = D_to_K0staruHNL + D0_to_KstaruHNL 
    self.Gamma_De = self.Gamma_De_leptonic + self.Gamma_De_semileptonic_pseudoscalar + self.Gamma_De_semileptonic_vector
    self.Gamma_Du = self.Gamma_Du_leptonic + self.Gamma_Du_semileptonic_pseudoscalar + self.Gamma_Du_semileptonic_vector

    # Inclusive B-D 
    self.Gamma_BDe = self.Gamma_De + self.Gamma_Be
    self.Gamma_BDu = self.Gamma_Du + self.Gamma_Bu


