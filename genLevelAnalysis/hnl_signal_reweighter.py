# http://home.thep.lu.se/~torbjorn/talks/fnal04lha.pdf

import ROOT
import numpy as np
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from scipy.constants import c as speed_of_light

from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p00244948974278_mu_Dirac_massiveAndCKM_LO
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p00282842712475_mu_Dirac_massiveAndCKM_LO
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p00316227766017_mu_Dirac_massiveAndCKM_LO
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p004472135955_mu_Dirac_massiveAndCKM_LO  
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p00547722557505_mu_Dirac_massiveAndCKM_LO
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p00707106781187_mu_Dirac_massiveAndCKM_LO
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p00836660026534_mu_Dirac_massiveAndCKM_LO
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p01_mu_Dirac_massiveAndCKM_LO            
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p0141421356237_mu_Dirac_massiveAndCKM_LO 
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p0157162336455_mu_Dirac_massiveAndCKM_LO 
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p0173205080757_mu_Dirac_massiveAndCKM_LO 
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p022360679775_mu_Dirac_massiveAndCKM_LO  
from CMGTools.HNL.samples.signals_2016 import HN3L_M_2_V_0p0350713558335_mu_Dirac_massiveAndCKM_LO 

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

samples = [
#     HN3L_M_2_V_0p00244948974278_mu_Dirac_massiveAndCKM_LO,
#     HN3L_M_2_V_0p00282842712475_mu_Dirac_massiveAndCKM_LO,
#     HN3L_M_2_V_0p00316227766017_mu_Dirac_massiveAndCKM_LO,
#     HN3L_M_2_V_0p004472135955_mu_Dirac_massiveAndCKM_LO  ,
#     HN3L_M_2_V_0p00547722557505_mu_Dirac_massiveAndCKM_LO,
    HN3L_M_2_V_0p00707106781187_mu_Dirac_massiveAndCKM_LO,
#     HN3L_M_2_V_0p00836660026534_mu_Dirac_massiveAndCKM_LO,
#     HN3L_M_2_V_0p01_mu_Dirac_massiveAndCKM_LO            ,
#     HN3L_M_2_V_0p0141421356237_mu_Dirac_massiveAndCKM_LO ,
#     HN3L_M_2_V_0p0157162336455_mu_Dirac_massiveAndCKM_LO ,
#     HN3L_M_2_V_0p0173205080757_mu_Dirac_massiveAndCKM_LO ,
#     HN3L_M_2_V_0p022360679775_mu_Dirac_massiveAndCKM_LO  ,
#     HN3L_M_2_V_0p0350713558335_mu_Dirac_massiveAndCKM_LO ,
]

# mysample = HN3L_M_2_V_0p00707106781187_mu_Dirac_massiveAndCKM_LO

branches = [
    'run',  
    'lumi', 
    'event',

    'w_pt',
    'w_eta',
    'w_phi',
    'w_mass',
    'w_q',
    'w_pdgid',

    'hnl_pt',
    'hnl_eta',
    'hnl_phi',
    'hnl_mass',
    'hnl_ct_lhe', # from LHE information
    'hnl_ct_reco', # from Lxyz = ct\gamma\beta --> ct = Lxyz/(\gamma\beta)
    'hnl_beta',  # Lorentz
    'hnl_gamma', # Lorentz

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

    'l2_pt',
    'l2_eta',
    'l2_phi',
    'l2_mass',
    'l2_q',
    'l2_pdgid',

    'hn_pt',
    'hn_eta',
    'hn_phi',
    'hn_mass',

    'n_pt',
    'n_eta',
    'n_phi',
    'n_mass',
    'n_pdgid',
    
    'Lxy', # 2D transverse displacement
    'Lxyz',  # 3D displacement
    'Lxy_cos', # cosine of the pointing angle in the transverse plane
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
handles['genp_pruned'] = ('prunedGenParticles' , Handle('std::vector<reco::GenParticle>')     )
handles['genp_packed'] = ('packedGenParticles' , Handle('std::vector<pat::PackedGenParticle>'))
handles['lhe']         = ('externalLHEProducer', Handle('LHEEventProduct'))

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


for mysample in samples:

    print 'processing', mysample.name

    # output file and tree gymnastics
    outfile = ROOT.TFile.Open('lifetimes_%s.root' %mysample.name, 'recreate')
    ntuple  = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
    tofill = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       

    # get to the real thing
    print 'loading files ...'
    events = Events(mysample.files[:75])
    print '... done!'

    for i, event in enumerate(events):
    
        # access the handles
        for k, v in handles.iteritems():
            event.getByLabel(v[0], v[1])
            setattr(event, k, v[1].product())
    
        if i%1000==0:
            percentage = float(i)/events.size()*100.
            print '\t===> processing %d / %d event \t completed %.1f%s' %(i, events.size(), percentage, '%')

        hepup = event.lhe.hepeup()
    
        # for each particle at LHE level, its lifetime is saved. If the particle decays immediately, it is set to 0.
        # http://home.thep.lu.se/~torbjorn/talks/fnal04lha.pdf
        ctaus = [hepup.VTIMUP.at(ictau) for ictau in xrange(hepup.VTIMUP.size()) if hepup.VTIMUP.at(ictau)>0.]
    
        # there should be only one particle at LHE level with lifetime > 0, that is the HNL
        event.hnl_ct_lhe = ctaus[0]
    
        # all gen particles
        event.genp = [ip for ip in event.genp_pruned] + [ip for ip in event.genp_packed]
    
        # find the W
        the_ws = sorted([ip for ip in event.genp_pruned if ip.isLastCopy() and ip.statusFlags().isPrompt() and abs(ip.pdgId())==24], key = lambda x : x.pt(), reverse=True)
        if len(the_ws):
            event.the_w = the_ws[0]
        else:
            event.the_w = None
    
        # get the heavy neutrino
        the_hns = [ip for ip in event.genp_pruned if abs(ip.pdgId()) in [9900012, 9990012] and ip.isLastCopy()] # 9900012 is Majorana, 9990012 is Dirac. Dirac comes in two species, particle and anti-particle!
        event.the_hn = the_hns[0] # one per event

        # prompt lepton
        event.the_pl = [ip for ip in event.genp_pruned if abs(ip.pdgId()) in [11,13,15] and ip.isPromptFinalState() and not isAncestor(event.the_hn, ip)][0]      

        # get the immediate daughters of the heavy neutrino decay
        event.the_hn.initialdaus = [event.the_hn.daughter(jj) for jj in range(event.the_hn.numberOfDaughters())]

        event.the_hn.lep1 = max([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt())
        event.the_hn.lep2 = min([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt())
        event.the_hn.neu  =     [ii for ii in event.the_hn.initialdaus if abs(ii.pdgId()) in [12, 14, 16]][0] # there can be only one

        event.the_dilepton = event.the_hn.lep1.p4() + event.the_hn.lep2.p4()

        # identify the primary vertex
        event.the_hn.the_pv = event.the_pl.vertex()

        # identify the secondary vertex
        event.the_hn.the_sv = event.the_hn.lep1.vertex()

        # 2D transverse and 3D displacement, Pythagoras
        event.Lxy  = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
                             (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2)

        event.Lxyz = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
                             (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2 + \
                             (event.the_hn.the_pv.z() - event.the_hn.the_sv.z())**2)
    
        # per event ct, as derived from the flight distance and Lorentz boost  
        event.the_hn.beta  = event.the_hn.p4().Beta()
        event.the_hn.gamma = event.the_hn.p4().Gamma()
        event.the_hn.ct_reco = event.Lxyz / (event.the_hn.beta * event.the_hn.gamma)


        # pointing angle    
        hn_pt_vect = ROOT.math.XYZVector(event.the_dilepton.px(),
                                         event.the_dilepton.py(),
                                         0.)
  
        Lxy_vect = ROOT.GlobalPoint(-1*((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())), 
                                    -1*((event.the_hn.the_pv.y() - event.the_hn.the_sv.y())),
                                     0)
  
        vperptau = ROOT.math.XYZVector(Lxy_vect.x(), Lxy_vect.y(), 0.)
  
        event.cos_pointing = vperptau.Dot(hn_pt_vect)/(vperptau.R()*hn_pt_vect.R())

        # reset before filling
        for k, v in tofill.iteritems(): tofill[k] = -99. # initialise before filling    
    
        tofill['run'        ] = event.eventAuxiliary().run()
        tofill['lumi'       ] = event.eventAuxiliary().luminosityBlock()
        tofill['event'      ] = event.eventAuxiliary().event()

        if event.the_w:
            tofill['w_pt'   ] = event.the_w.pt()     
            tofill['w_eta'  ] = event.the_w.eta()    
            tofill['w_phi'  ] = event.the_w.phi()    
            tofill['w_mass' ] = event.the_w.mass()   
            tofill['w_q'    ] = event.the_w.charge()   
            tofill['w_pdgid'] = event.the_w.pdgId()   
         
        tofill['hnl_pt'     ] = event.the_hn.pt()     
        tofill['hnl_eta'    ] = event.the_hn.eta()    
        tofill['hnl_phi'    ] = event.the_hn.phi()    
        tofill['hnl_mass'   ] = event.the_hn.mass()   
        tofill['hnl_ct_lhe' ] = event.hnl_ct_lhe 
        tofill['hnl_ct_reco'] = event.the_hn.ct_reco * 10. # convert cm to mm 
        tofill['hnl_beta'   ] = event.the_hn.beta  
        tofill['hnl_gamma'  ] = event.the_hn.gamma  

        tofill['l0_pt'      ] = event.the_pl.pt()
        tofill['l0_eta'     ] = event.the_pl.eta()
        tofill['l0_phi'     ] = event.the_pl.phi()
        tofill['l0_mass'    ] = event.the_pl.mass()
        tofill['l0_q'       ] = event.the_pl.charge()
        tofill['l0_pdgid'   ] = event.the_pl.pdgId()

        tofill['l1_pt'      ] = event.the_hn.lep1.pt()
        tofill['l1_eta'     ] = event.the_hn.lep1.eta()
        tofill['l1_phi'     ] = event.the_hn.lep1.phi()
        tofill['l1_mass'    ] = event.the_hn.lep1.mass()
        tofill['l1_q'       ] = event.the_hn.lep1.charge()
        tofill['l1_pdgid'   ] = event.the_hn.lep1.pdgId()

        tofill['l2_pt'      ] = event.the_hn.lep2.pt()
        tofill['l2_eta'     ] = event.the_hn.lep2.eta()
        tofill['l2_phi'     ] = event.the_hn.lep2.phi()
        tofill['l2_mass'    ] = event.the_hn.lep2.mass()
        tofill['l2_q'       ] = event.the_hn.lep2.charge()
        tofill['l2_pdgid'   ] = event.the_hn.lep2.pdgId()

        tofill['hn_pt'      ] = event.the_dilepton.pt()     
        tofill['hn_eta'     ] = event.the_dilepton.eta()    
        tofill['hn_phi'     ] = event.the_dilepton.phi()    
        tofill['hn_mass'    ] = event.the_dilepton.mass()   

        tofill['n_pt'       ] = event.the_hn.neu.pt()
        tofill['n_eta'      ] = event.the_hn.neu.eta()
        tofill['n_phi'      ] = event.the_hn.neu.phi()
        tofill['n_pdgid'    ] = event.the_hn.neu.pdgId()

        tofill['Lxy'        ] = event.Lxy
        tofill['Lxyz'       ] = event.Lxyz
        tofill['Lxy_cos'    ] = event.cos_pointing

        for iv2 in new_v2s:
            tofill['weight_%s'      %(str(iv2).replace('-', 'm'))] = weight_to_new_ctau(mysample.ctau, mysample.v2, iv2, event.hnl_ct_lhe)[0]
            tofill['ctau_%s'        %(str(iv2).replace('-', 'm'))] = weight_to_new_ctau(mysample.ctau, mysample.v2, iv2, event.hnl_ct_lhe)[1]
            tofill['xs_scale_to_%s' %(str(iv2).replace('-', 'm'))] = iv2 / mysample.v2
        
        ntuple.Fill(array('f',tofill.values()))


    outfile.cd()
    ntuple.Write()
    outfile.Close()
