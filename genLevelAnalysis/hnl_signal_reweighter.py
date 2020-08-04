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


class Vertex(object):
   def __init__(self, index=0, x=0, y=0, z=0, pT=0):    
      self.index = index
      self.x     = x
      self.y     = y
      self.z     = z
      self.pT    = pT

def getpT(input):
   return input.pT

branches = [
    'run',  
    'lumi', 
    'event',
     
    # the mother
    'b_pt',
    'b_eta',
    'b_phi',
    'b_mass',
    'b_q',
    'b_pdgid',
    'b_ct_reco',
    
    # daughters of the B
    # # the HNL
    'hnl_pt',
    'hnl_eta',
    'hnl_phi',
    'hnl_mass',
    'hnl_q',
    'hnl_pdgid',
    #'hnl_ct_lhe', # from LHE information
    'hnl_ct_reco', # from Lxyz = ct\gamma\beta --> ct = Lxyz/(\gamma\beta)
    'hnl_beta',  # Lorentz
    'hnl_gamma', # Lorentz

    # # the D meson
    'd_pt',
    'd_eta',
    'd_phi',
    'd_mass',
    'd_q',
    'd_pdgid',

    # # the triggered lepton
    'l0_pt',
    'l0_eta',
    'l0_phi',
    'l0_mass',
    'l0_q',
    'l0_pdgid',

    # daughters of the D meson
    # # the pion
    'pi_pt',
    'pi_eta',
    'pi_phi',
    'pi_mass',
    'pi_q',
    'pi_pdgid',
    
    # # the kaon
    'k_pt',
    'k_eta',
    'k_phi',
    'k_mass',
    'k_q',
    'k_pdgid',
    
    # daughters of the HNL
    # # the lepton
    'l1_pt',
    'l1_eta',
    'l1_phi',
    'l1_mass',
    'l1_q',
    'l1_pdgid',

    # # the pion
    'pi1_pt',
    'pi1_eta',
    'pi1_phi',
    'pi1_mass',
    'pi1_q',
    'pi1_pdgid',
   
    # invariant masses
    'lep_pi_invmass',
    'k_pi_invmass',
    'hn_d_pl_invmass',

    'Lxy', # 2D transverse displacement for the HNL
    'Lxyz',  # 3D displacement for the HNL
    'Lxy_cos', # cosine of the pointing angle in the transverse plane

    'Lxyz_b', #3D displacement of the B wrt to primary vertex
    'Lxyz_l0' #3D displacement of the triggered lepton wrt to B vertex
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
handles['genP'] = ('genParticles' , Handle('std::vector<reco::GenParticle>'))
#handles['genp_packed'] = ('packedGenParticles' , Handle('std::vector<pat::PackedGenParticle>'))
#handles['lhe']         = ('externalLHEProducer', Handle('LHEEventProduct'))

# output file and tree gymnastics
outfile = ROOT.TFile.Open('genNTuples_ctau100mm_fullEvtGen_10k.root', 'recreate')
#outfile = ROOT.TFile.Open('genNTuples_10k.root', 'recreate')
ntuple  = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       



# get to the real thing
print 'loading the file ...'
#events = Events('samples/BPH-test_numEvent50.root')
#events = Events('samples/BPH-test_numEvent100_Bmod.root')
#events = Events('root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/TEST0/BPH-test_numEvent10000.root')
#events = Events('root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/TEST_ctau100mm/BPH-test_numEvent10000.root')
#events = Events('root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/TEST_ctau1_m1/BPH-test_numEvent10000.root')
events = Events('root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/TEST_ctau1_m1_fullEvtGen/BPH-test_numEvent10000.root')
#events = Events('root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/TEST_mupt5_ctau100mm/BPH-test_numEvent80000.root')
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
  #event.genp = [ip for ip in event.genP] #+ [ip for ip in event.genp_packed]
      
  # get the heavy neutrino
  the_hns = [ip for ip in event.genP if abs(ip.pdgId())==9900015 and ip.isLastCopy()] 
  if len(the_hns):
     event.the_hn = the_hns[0] # one per event
     #print '0. found hnl of pdgId {a}'.format(a=event.the_hn.pdgId())
  
  # find the B mother  
  event.the_hn.mothers = [event.the_hn.mother(jj) for jj in range(event.the_hn.numberOfMothers())]
  the_b_mothers = sorted([ii for ii in event.the_hn.mothers if abs(ii.pdgId())==521], key = lambda x : x.pt(), reverse=True)
  if len(the_b_mothers):
    event.the_b_mother = the_b_mothers[0]
    #print '1. found a B of pdgId {a}'.format(a=event.the_b_mother.pdgId())
  else:
    event.the_b_mother = None

  # get the other two daughters of the B meson
  event.the_b_mother.daughters = [event.the_b_mother.daughter(jj) for jj in range(event.the_b_mother.numberOfDaughters())]

  # # first the D0 meson
  the_ds = sorted([ii for ii in event.the_b_mother.daughters if abs(ii.pdgId())==421], key = lambda x : x.pt(), reverse=True)
  if len(the_ds):
    event.the_d = the_ds[0]
    #print '2. found a D of pdgId {a}'.format(a=event.the_d.pdgId())
  else:
    event.the_d = None

  # # then the triggered lepton
  the_pls = sorted([ii for ii in event.the_b_mother.daughters if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt(), reverse=True)
  if len(the_pls):
    event.the_pl = the_pls[0]
    #print '3. found a triggered lepton of pdgId {a}'.format(a=event.the_pl.pdgId())
  else:
    event.the_pl = None
  
  # D0s daughters
  event.the_d.daughters = [event.the_d.daughter(jj) for jj in range(event.the_d.numberOfDaughters())]
  
  # #find the pion
  the_pis = sorted([ii for ii in event.the_d.daughters if abs(ii.pdgId())==211], key = lambda x : x.pt(), reverse=True)
  if len(the_pis):
      event.the_pi = the_pis[0]
      #print '4. found a pion with pdgId {a}'.format(a=event.the_pi.pdgId())
  else:
      event.the_pi = None
      
  # find the kaon
  the_ks = sorted([ii for ii in event.the_d.daughters if abs(ii.pdgId())==321], key = lambda x : x.pt(), reverse=True)
  if len(the_ks):
      event.the_k = the_ks[0]
      #print '5. found a kaon with pdgId {a}'.format(a=event.the_k.pdgId())
  else:
      event.the_k = None
 

  # HNL daughters
  event.the_hn.initialdaus = [event.the_hn.daughter(jj) for jj in range(event.the_hn.numberOfDaughters())]

  # # the lepton
  the_lep_daughters = sorted([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt(), reverse=True)
  if len(the_lep_daughters):
     #print '6. found a daughter lepton'
     event.the_hn.lep = the_lep_daughters[0]
  else: 
     event.the_hn.lep = None
 
  # # the pion
  the_pi_daughters = sorted([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId())==211], key = lambda x : x.pt(), reverse=True)
  if len(the_pi_daughters):
     #print '7. found a daughter pion'
     event.the_hn.pi = the_pi_daughters[0]
  else:
     event.the_hn.pi = None
 
  # invariant masses   
  # # to get the invariant mass of the HNL daughters
  if len(the_pi_daughters) and len(the_lep_daughters):   
    event.the_hnldaughters = event.the_hn.lep.p4() + event.the_hn.pi.p4()
  
  # # to get the invariant mass of the D0 daughters
  if len(the_ks) and len(the_pis):   
    event.the_d0daughters = event.the_k.p4() + event.the_pi.p4()

  # # to get the invariant mass of the B daughters
  if len(the_ds) and len(the_pls) and len(the_hns):
     event.the_bdaughters = event.the_hn.p4() + event.the_d.p4() + event.the_pl.p4()
     #print event.the_bdaughters.mass()
    
  # identify the primary vertex
  # for that, needs the triggered lepton
  if len(the_pls):
    event.the_hn.the_pv = event.the_pl.vertex()
  
  # identify the secondary vertex
  if len(the_lep_daughters):
    event.the_hn.the_sv = event.the_hn.lep.vertex()
  
  # 2D transverse and 3D displacement, Pythagoras
  if len(the_pls) and len(the_lep_daughters):
    #print 'pl vx: {a} l1 vx: {b} hnl vx: {c}'.format(a=event.the_pl.vx(), b=event.the_hn.lep.vx(), c=event.the_hn.vx()) 
    event.Lxy  = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
                         (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2)

    event.Lxyz = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
                         (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2 + \
                         (event.the_hn.the_pv.z() - event.the_hn.the_sv.z())**2)

    #event.Lxyz = np.sqrt((event.the_hn.vx() - event.the_hn.lep.vx())**2 + \
    #                  (event.the_hn.vy() - event.the_hn.lep.vy())**2 + \
    #                  (event.the_hn.vz() - event.the_hn.lep.vz())**2)

    #print 'x vertex: triggered lep | hnl | hnl daughter lep | diff(hnl,hnl daughter lep)'
    #print '{a} | {b} | {c} | {d}'.format(a=event.the_pl.vx(), b=event.the_hn.vx(), c=event.the_hn.pi.vx(), d=event.the_hn.vx()-event.the_hn.pi.vx())

    #print 'y vertex: triggered lep | hnl | hnl daughter lep | diff(hnl,hnl daughter lep)'
    #print '{a} | {b} | {c} | {d}'.format(a=event.the_pl.vy(), b=event.the_hn.vy(), c=event.the_hn.lep.vy(), d=event.the_hn.vy()-event.the_hn.lep.vy())


    #print 'z vertex: triggered lep | hnl | hnl daughter lep | diff(hnl,hnl daughter lep)'
    #print '{a} | {b} | {c} | {d}'.format(a=event.the_pl.vz(), b=event.the_hn.vz(), c=event.the_hn.lep.vz(), d=event.the_hn.vz()-event.the_hn.lep.vz())
    
    #print 'pvx: {a} svx: {b} pvy: {c} svy: {d} Lxy: {e}'.format(a=event.the_hn.the_pv.x(), b= event.the_hn.the_sv.x(), c=event.the_hn.the_pv.y(), d=event.the_hn.the_sv.y(), e=event.Lxy)
    #print 'pvz: {a} svz: {b} Lxyz: {c}'.format(a=event.the_hn.the_pv.z(), b= event.the_hn.the_sv.z(), c=event.Lxyz)
  # per event ct, as derived from the flight distance and Lorentz boost  
  event.the_hn.beta  = event.the_hn.p4().Beta()
  event.the_hn.gamma = event.the_hn.p4().Gamma()
  
  # we get the lifetime from the kinematics 
  if len(the_pls) and len(the_lep_daughters):
    event.the_hn.ct_reco = event.Lxyz / (event.the_hn.beta * event.the_hn.gamma)
    #print 'hnl ct reco: {a}'.format(a=event.the_hn.ct_reco)
  # pointing angle    
  #if len(the_pls) and len(the_lep_daughters):
  #  hn_pt_vect = ROOT.math.XYZVector(event.the_hnldaughters.px(),
  #                                   event.the_hnldaughters.py(),
  #                                   0.)

  #  Lxy_vect = ROOT.GlobalPoint(-1*((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())), 
  #                              -1*((event.the_hn.the_pv.y() - event.the_hn.the_sv.y())),
  #                               0)

  #  vperptau = ROOT.math.XYZVector(Lxy_vect.x(), Lxy_vect.y(), 0.)
 
    #event.cos_pointing = vperptau.Dot(hn_pt_vect)/(vperptau.R()*hn_pt_vect.R())


  '''
  # get the primary vertex (of the event)
  # all gen particles
  the_genps = [ip for ip in event.genP] #+ [ip for ip in event.genp_packed]
  #print 'len tot: {a}'.format(a=len(the_genps))

  allvertices = []
  for genParticle in the_genps:
    tpl = genParticle.vx(), genParticle.vy(), genParticle.vz()
    allvertices.append(tpl)
  
  #print 'len all vertices: {a}'.format(a=len(allvertices))
  #print 'len vertices: {a}'.format(a=len(set(allvertices)))

  vertices = []
  idx = 0
  for vertex in set(allvertices):
     vertices.append(Vertex(idx, vertex[0], vertex[1], vertex[2], 0))
     idx = idx + 1

  for vertex in vertices:
    the_pt = 0
    for genParticle in the_genps:
      if genParticle.vx() == vertex.x and genParticle.vy() == vertex.y and genParticle.vz() == vertex.z:
        the_pt = the_pt + genParticle.pt()
    vertex.pT = the_pt


  #print 'before'  
  #for vertex in vertices:
  #  print '{a} {b} {c} {d} {e}'.format(a=vertex.index, b=vertex.x, c=vertex.y, d=vertex.z, e=vertex.pT)

  vertices.sort(key=getpT, reverse=True)  
  
  #print '\n after' 
  #for vertex in vertices:
  #  print '{a} {b} {c} {d} {e}'.format(a=vertex.index, b=vertex.x, c=vertex.y, d=vertex.z, e=vertex.pT)

  event.theprimaryvertex = vertices[0]

  #event.theBvertex = Vertex(0, event.the_b

  event.Lxy_bmother = np.sqrt((event.theprimaryvertex.x - event.the_b_mother.vx())**2 + \
                      (event.theprimaryvertex.y - event.the_b_mother.vy())**2)

  event.Lxyz_bmother = np.sqrt((event.theprimaryvertex.x - event.the_b_mother.vx())**2 + \
                      (event.theprimaryvertex.y - event.the_b_mother.vy())**2 + \
                      (event.theprimaryvertex.z - event.the_b_mother.vz())**2)
  #print event.Lxyz_bmother
  '''

  event.Lxyz_pl = np.sqrt((event.the_b_mother.vx() - event.the_pl.vx())**2 + \
                      (event.the_b_mother.vy() - event.the_pl.vy())**2 + \
                      (event.the_b_mother.vz() - event.the_pl.vz())**2)
  
  #print event.Lxyz_pl
 
  # get the lifetime of the B
  event.the_b_mother.beta  = event.the_b_mother.p4().Beta()
  event.the_b_mother.gamma = event.the_b_mother.p4().Gamma()
 
  event.the_b_mother.ct_reco = event.Lxyz_pl / (event.the_b_mother.beta * event.the_b_mother.gamma)
   
  # reset before filling
  for k, v in tofill.iteritems(): tofill[k] = -99. # initialise before filling

  tofill['run'        ] = event.eventAuxiliary().run()
  tofill['lumi'       ] = event.eventAuxiliary().luminosityBlock()
  tofill['event'      ] = event.eventAuxiliary().event()
   
  if event.the_b_mother:
      tofill['b_pt'   ] = event.the_b_mother.pt()     
      tofill['b_eta'  ] = event.the_b_mother.eta()    
      tofill['b_phi'  ] = event.the_b_mother.phi()    
      tofill['b_mass' ] = event.the_b_mother.mass()   
      tofill['b_q'    ] = event.the_b_mother.charge()   
      tofill['b_pdgid'] = event.the_b_mother.pdgId()   
      tofill['b_ct_reco'] = event.the_b_mother.ct_reco   
   
  tofill['hnl_pt'     ] = event.the_hn.pt()     
  tofill['hnl_eta'    ] = event.the_hn.eta()    
  tofill['hnl_phi'    ] = event.the_hn.phi()    
  tofill['hnl_mass'   ] = event.the_hn.mass()   
  tofill['hnl_q'      ] = event.the_hn.charge()   
  #tofill['hnl_ct_lhe' ] = event.hnl_ct_lhe 
  if len(the_lep_daughters):
    tofill['hnl_ct_reco'] = event.the_hn.ct_reco #* 10. # convert cm to mm 
  tofill['hnl_beta'   ] = event.the_hn.beta  
  tofill['hnl_gamma'  ] = event.the_hn.gamma  
  tofill['hnl_pdgid'  ] = event.the_hn.pdgId()  

  if event.the_d:
      tofill['d_pt'   ] = event.the_d.pt()     
      tofill['d_eta'  ] = event.the_d.eta()    
      tofill['d_phi'  ] = event.the_d.phi()    
      tofill['d_mass' ] = event.the_d.mass()   
      tofill['d_q'    ] = event.the_d.charge()   
      tofill['d_pdgid'] = event.the_d.pdgId()   

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

  # invariant mass
     tofill['lep_pi_invmass' ] = event.the_hnldaughters.mass()
  tofill['k_pi_invmass' ] = event.the_d0daughters.mass()
  tofill['hn_d_pl_invmass'] = event.the_bdaughters.mass()
  
  
  if len(the_pls) and len(the_lep_daughters):
    tofill['Lxy'        ] = event.Lxy
    tofill['Lxyz'       ] = event.Lxyz
    #tofill['Lxy_cos'    ] = event.cos_pointing

  #tofill['Lxyz_b']  = event.Lxyz_bmother
  tofill['Lxyz_l0'] = event.Lxyz_pl

  #for iv2 in new_v2s:
      #tofill['weight_%s'      %(str(iv2).replace('-', 'm'))] = weight_to_new_ctau(mysample.ctau, mysample.v2, iv2, event.hnl_ct_lhe)[0]
      #tofill['ctau_%s'        %(str(iv2).replace('-', 'm'))] = weight_to_new_ctau(mysample.ctau, mysample.v2, iv2, event.hnl_ct_lhe)[1]
      #tofill['xs_scale_to_%s' %(str(iv2).replace('-', 'm'))] = iv2 / mysample.v2
  
  ntuple.Fill(array('f',tofill.values()))


outfile.cd()
ntuple.Write()
outfile.Close()
