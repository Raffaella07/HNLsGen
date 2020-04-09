# http://home.thep.lu.se/~torbjorn/talks/fnal04lha.pdf

import sys
import ROOT
import numpy as np
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from scipy.constants import c as speed_of_light


def isAncestor(a, p):
    if a == p :
        return True
    for i in xrange(0,p.numberOfMothers()):
        if isAncestor(a,p.mother(i)):
            return True
    return False

def weight_to_new_ctau(old_ctau, old_v2, new_v2, ct):
    '''
    Returns an event weight based on the ratio of the normalised lifetime distributions.
    old_ctau: reference ctau
    old_v2  : reference coupling squared
    new_v2  : target coupling squared
    ct      : heavy neutrino lifetime in the specific event
    '''
    new_ctau = old_ctau * old_v2 / new_v2
    weight = old_ctau/new_ctau * np.exp( (1./old_ctau - 1./new_ctau) * ct )
    return weight, new_ctau

def scale_to_new_xs(old_v2, new_v2):
    '''
    At equal mass, the xs is simply proportional to v2 
    '''
    return new_v2 / old_v2


branches = [
    'run',  
    'lumi', 
    'event',

    'b_pt',
    'b_eta',
    'b_phi',
    'b_mass',
    'b_q',
    'b_pdgid',
    
    'hnl_pt',
    'hnl_eta',
    'hnl_phi',
    'hnl_mass',
    'hnl_q',
    'hnl_pdgid',
    #'hnl_ct_lhe', # from LHE information
    #'hnl_ct_reco', # from Lxyz = ct\gamma\beta --> ct = Lxyz/(\gamma\beta)
    'hnl_beta',  # Lorentz
    'hnl_gamma', # Lorentz

    'pi_pt',
    'pi_eta',
    'pi_phi',
    'pi_mass',
    'pi_q',
    'pi_pdgid',
    
    'k_pt',
    'k_eta',
    'k_phi',
    'k_mass',
    'k_q',
    'k_pdgid',
    
    'l0_pt',
    'l0_eta',
    'l0_phi',
    'l0_mass',
    'l0_q',
    'l0_pdgid',

    'l1_pt',
    'l1_eta',
    'l1_phi',
    'l1_mass',
    'l1_q',
    'l1_pdgid',

    'pi1_pt',
    'pi1_eta',
    'pi1_phi',
    'pi1_mass',
    'pi1_q',
    'pi1_pdgid',
    
    #'Lxy', # 2D transverse displacement
    #'Lxyz',  # 3D displacement
    #'Lxy_cos', # cosine of the pointing angle in the transverse plane
]

# couplings to be tested, for which the reweight is run
new_v2s = [
    1e-10, 
    5e-10, 
    1e-9, 
    5e-9, 
    1e-8, 
    5e-8, 
    1e-7, 
    5e-7, 
    1e-6, 
    5e-6, 
    6e-06, 
    8e-06, 
    1e-5, 
    2e-5, 
    3e-5, 
    4e-5, 
    5e-5, 
    7e-05, 
    0.0001, 
    0.0002, 
    0.00025, 
    0.0003, 
    0.0005, 
    0.0012,
]

# ctau weights, each coupling gets its own weight
weights = OrderedDict(zip(new_v2s, np.ones(len(new_v2s))))

# add weight branches
for vv in new_v2s:
    branches.append('weight_%s'      %(str(vv).replace('-', 'm')))
    branches.append('ctau_%s'        %(str(vv).replace('-', 'm')))
    branches.append('xs_scale_to_%s' %(str(vv).replace('-', 'm')))

handles = OrderedDict()
handles['genP'] = ('genParticles' , Handle('std::vector<reco::GenParticle>')     )
#handles['genp_packed'] = ('packedGenParticles' , Handle('std::vector<pat::PackedGenParticle>'))
#handles['lhe']         = ('externalLHEProducer', Handle('LHEEventProduct'))

# output file and tree gymnastics
outfile = ROOT.TFile.Open('lifetimes_test.root', 'recreate')
ntuple  = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       



# get to the real thing
print 'loading the file ...'
#events = Events('samples/BPH-test_numEvent50.root')
events = Events('root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/TEST0/BPH-test_numEvent10000.root')
print '... done!'

for i, event in enumerate(events):
  # access the handles
  for k, v in handles.iteritems():
      event.getByLabel(v[0], v[1])
      setattr(event, k, v[1].product())

  if i%1000==0:
      percentage = float(i)/events.size()*100.
      print '\t===> processing %d / %d event \t completed %.1f%s' %(i, events.size(), percentage, '%')

  #print '\n Event {a}'.format(a=i)

  #for the moment, comment everything that has to do with lhe   
  #hepup = event.lhe.hepeup()

  # for each particle at LHE level, its lifetime is saved. If the particle decays immediately, it is set to 0.
  # http://home.thep.lu.se/~torbjorn/talks/fnal04lha.pdf
  #ctaus = [hepup.VTIMUP.at(ictau) for ictau in xrange(hepup.VTIMUP.size()) if hepup.VTIMUP.at(ictau)>0.]

  # there should be only one particle at LHE level with lifetime > 0, that is the HNL
  #event.hnl_ct_lhe = ctaus[0]
  #####

  # all gen particles
  event.genp = [ip for ip in event.genP] #+ [ip for ip in event.genp_packed]


  # find the B
  # we sort the Bs from highest pt to lowest
  #the_bs = sorted([ip for ip in event.genP if ip.isLastCopy() and ip.statusFlags().isPrompt() and abs(ip.pdgId())==521], key = lambda x : x.pt(), reverse=True)
  #if len(the_bs):
      #print '1. found a B'
  #    event.the_b = the_bs[0]
  #else:
  #    event.the_b = None
  
  # find the mu
  # we sort the muons from highest pt to lowest
  #the_mus = sorted([ip for ip in event.genP if ip.isLastCopy() and ip.statusFlags().isPrompt() and abs(ip.pdgId())==13], key = lambda x : x.pt(), reverse=True)
  #if len(the_mus):
  #    #print 'found a mu'
  #    event.the_mu = the_mus[0]
  #else:
  #    event.the_mu = None
  
  # D0s daughters
  # find the pion
  # we sort the pions from highest pt to lowest
  the_pis = sorted([ip for ip in event.genP if ip.isLastCopy() and ip.statusFlags().isPrompt() and abs(ip.pdgId())==211], key = lambda x : x.pt(), reverse=True)
  if len(the_pis):
      #print '2. found a pion'
      event.the_pi = the_pis[0]
  else:
      event.the_pi = None
      
  # find the kaon
  # we sort the kaons from highest pt to lowest
  the_ks = sorted([ip for ip in event.genP if ip.isLastCopy() and ip.statusFlags().isPrompt() and abs(ip.pdgId())==321], key = lambda x : x.pt(), reverse=True)
  if len(the_ks):
      #print '3. found a kaon'
      event.the_k = the_ks[0]
  else:
      event.the_k = None
      
  # get the heavy neutrino
  the_hns = [ip for ip in event.genP if abs(ip.pdgId()) in [9900015, -9990015] and ip.isLastCopy()] # 9900012 is Majorana, 9990012 is Dirac. Dirac comes in two species, particle and anti-particle!
  if len(the_hns):
     event.the_hn = the_hns[0] # one per event
     #print '4. found hnls of pdgId {a}'.format(a=event.the_hn.pdgId())
      
  # find the B mother
  # we sort the Bs from highest pt to lowest
  #the_bs = sorted([ip for ip in event.genP if ip.isLastCopy() and ip.statusFlags().isPrompt() and abs(ip.pdgId())==521], key = lambda x : x.pt(), reverse=True)
  the_bs = sorted([ip for ip in event.genP if ip.isLastCopy() and abs(ip.pdgId())==521 and isAncestor(event.the_hn, ip)], key = lambda x : x.pt(), reverse=True)
  if len(the_bs):
      #print '1. found a B'
      event.the_b = the_bs[0]
  else:
      event.the_b = None
  
  # prompt lepton
  #event.the_pl = [ip for ip in event.genP if abs(ip.pdgId()) in [11,13,15] and ip.isPromptFinalState() and not isAncestor(event.the_hn, ip)]#[0]
  the_pls = [ip for ip in event.genP if abs(ip.pdgId()) in [11, 13, 15] and ip.isLastCopy() and not isAncestor(event.the_hn, ip)] 
  #the_pls = [ip for ip in event.genP if abs(ip.pdgId()) in [11,13,15] and ip.isPromptFinalState() and not isAncestor(event.the_hn, ip)]
  if len(the_pls):
     #print '5. found a prompt lepton'
     event.the_pl = the_pls[0]
  else:
     event.the_pl = None
  
  
  # get the immediate daughters of the heavy neutrino decay
  event.the_hn.initialdaus = [event.the_hn.daughter(jj) for jj in range(event.the_hn.numberOfDaughters())]

  #event.the_hn.lep = max([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt())
  the_lep_daughters = sorted([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt())
  if len(the_lep_daughters):
     #print '6. found a daughter lepton'
     event.the_hn.lep = the_lep_daughters[0]
  else: 
     event.the_hn.lep = None
  
  the_pi_daughters = sorted([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId())==211], key = lambda x : x.pt())
  if len(the_pi_daughters):
     #print '7. found a daughter pion'
     event.the_hn.pi = the_pi_daughters[0]
  else:
     event.the_hn.pi = None

  #event.the_dilepton = event.the_hn.lep1.p4() + event.the_hn.lep2.p4()
 
  # identify the primary vertex
  # for that, needs the prompt lepton
  #if len(the_pls):
  #  event.the_hn.the_pv = event.the_pl.vertex()
  
  # identify the secondary vertex
  #if len(the_lep_daughters):
  #  event.the_hn.the_sv = event.the_hn.lep.vertex()
  
  # 2D transverse and 3D displacement, Pythagoras
  #event.Lxy  = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
  #                     (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2)

  #event.Lxyz = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
  #                     (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2 + \
  #                     (event.the_hn.the_pv.z() - event.the_hn.the_sv.z())**2)

  # per event ct, as derived from the flight distance and Lorentz boost  
  event.the_hn.beta  = event.the_hn.p4().Beta()
  event.the_hn.gamma = event.the_hn.p4().Gamma()
  #event.the_hn.ct_reco = event.Lxyz / (event.the_hn.beta * event.the_hn.gamma)

  # pointing angle    
  #hn_pt_vect = ROOT.math.XYZVector(event.the_dilepton.px(),
  #                                 event.the_dilepton.py(),
  #                                 0.)

  #Lxy_vect = ROOT.GlobalPoint(-1*((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())), 
  #                            -1*((event.the_hn.the_pv.y() - event.the_hn.the_sv.y())),
  #                             0)

  #vperptau = ROOT.math.XYZVector(Lxy_vect.x(), Lxy_vect.y(), 0.)

  #event.cos_pointing = vperptau.Dot(hn_pt_vect)/(vperptau.R()*hn_pt_vect.R())
  
  # reset before filling
  for k, v in tofill.iteritems(): tofill[k] = -99. # initialise before filling

  tofill['run'        ] = event.eventAuxiliary().run()
  tofill['lumi'       ] = event.eventAuxiliary().luminosityBlock()
  tofill['event'      ] = event.eventAuxiliary().event()
   
  if event.the_b:
      tofill['b_pt'   ] = event.the_b.pt()     
      tofill['b_eta'  ] = event.the_b.eta()    
      tofill['b_phi'  ] = event.the_b.phi()    
      tofill['b_mass' ] = event.the_b.mass()   
      tofill['b_q'    ] = event.the_b.charge()   
      tofill['b_pdgid'] = event.the_b.pdgId()   
   
  tofill['hnl_pt'     ] = event.the_hn.pt()     
  tofill['hnl_eta'    ] = event.the_hn.eta()    
  tofill['hnl_phi'    ] = event.the_hn.phi()    
  tofill['hnl_mass'   ] = event.the_hn.mass()   
  tofill['hnl_q'      ] = event.the_hn.charge()   
  #tofill['hnl_ct_lhe' ] = event.hnl_ct_lhe 
  #tofill['hnl_ct_reco'] = event.the_hn.ct_reco * 10. # convert cm to mm 
  tofill['hnl_beta'   ] = event.the_hn.beta  
  tofill['hnl_gamma'  ] = event.the_hn.gamma  
  tofill['hnl_pdgid'  ] = event.the_hn.pdgId()  

  if event.the_k:
      tofill['k_pt'   ] = event.the_k.pt()     
      tofill['k_eta'  ] = event.the_k.eta()    
      tofill['k_phi'  ] = event.the_k.phi()    
      tofill['k_mass' ] = event.the_k.mass()   
      tofill['k_q'    ] = event.the_k.charge()   
      tofill['k_pdgid'] = event.the_k.pdgId()   
   
  if event.the_pi:
      tofill['pi_pt'   ] = event.the_pi.pt()     
      tofill['pi_eta'  ] = event.the_pi.eta()    
      tofill['pi_phi'  ] = event.the_pi.phi()    
      tofill['pi_mass' ] = event.the_pi.mass()   
      tofill['pi_q'    ] = event.the_pi.charge()   
      tofill['pi_pdgid'] = event.the_pi.pdgId()   

  if event.the_pl:
     tofill['l0_pt'      ] = event.the_pl.pt()
     tofill['l0_eta'     ] = event.the_pl.eta()
     tofill['l0_phi'     ] = event.the_pl.phi()
     tofill['l0_mass'    ] = event.the_pl.mass()
     tofill['l0_q'       ] = event.the_pl.charge()
     tofill['l0_pdgid'   ] = event.the_pl.pdgId()

  if event.the_hn.lep:
     tofill['l1_pt'      ] = event.the_hn.lep.pt()
     tofill['l1_eta'     ] = event.the_hn.lep.eta()
     tofill['l1_phi'     ] = event.the_hn.lep.phi()
     tofill['l1_mass'    ] = event.the_hn.lep.mass()
     tofill['l1_q'       ] = event.the_hn.lep.charge()
     tofill['l1_pdgid'   ] = event.the_hn.lep.pdgId()

  if event.the_hn.pi:   
     tofill['pi1_pt'      ] = event.the_hn.pi.pt()
     tofill['pi1_eta'     ] = event.the_hn.pi.eta()
     tofill['pi1_phi'     ] = event.the_hn.pi.phi()
     tofill['pi1_mass'    ] = event.the_hn.pi.mass()
     tofill['pi1_q'       ] = event.the_hn.pi.charge()
     tofill['pi1_pdgid'   ] = event.the_hn.pi.pdgId()

  #tofill['Lxy'        ] = event.Lxy
  #tofill['Lxyz'       ] = event.Lxyz
  #tofill['Lxy_cos'    ] = event.cos_pointing

  #for iv2 in new_v2s:
      #tofill['weight_%s'      %(str(iv2).replace('-', 'm'))] = weight_to_new_ctau(mysample.ctau, mysample.v2, iv2, event.hnl_ct_lhe)[0]
      #tofill['ctau_%s'        %(str(iv2).replace('-', 'm'))] = weight_to_new_ctau(mysample.ctau, mysample.v2, iv2, event.hnl_ct_lhe)[1]
      #tofill['xs_scale_to_%s' %(str(iv2).replace('-', 'm'))] = iv2 / mysample.v2
  
  ntuple.Fill(array('f',tofill.values()))


outfile.cd()
ntuple.Write()
outfile.Close()
