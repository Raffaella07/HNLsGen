import os
import sys
import numpy as np
from ROOT import TTree, TFile, TH1F, TEfficiency, TGraph, TCanvas, gROOT, TAxis, TMath, TLegend, gStyle, gPad, TLine, TGraphAsymmErrors, TGraphErrors, TChain
import ROOT
sys.path.append('/work/mratti/plotting/myplotting')
from spares import *
from glob import glob
import array

# couplings to be tested, for which the reweight is run

#from genTreeProducer import new_vvs
#small_new_vvs = new_vvs[::4]  # only select one every three elements
#small_new_vvs.reverse()
#new_vvs.reverse()
#small_new_vvs = new_vvs[0:10]

from python.common import *

global graph_saver
graph_saver=[]


def getOverflowedHisto(h):
  htemp = h.Clone()
  nbins = htemp.GetNbinsX()
  last_plus_overflow = htemp.GetBinContent(nbins) + htemp.GetBinContent(nbins+1)
  last_plus_overflow_error = TMath.Sqrt( htemp.GetBinError(nbins)*htemp.GetBinError(nbins)  + htemp.GetBinError(nbins+1)*htemp.GetBinError(nbins+1))
  htemp.SetBinContent(nbins,last_plus_overflow )
  htemp.SetBinError(nbins,last_plus_overflow_error)
  return htemp


class SampleQuantity(object):
  def __init__(self, name='vv', axis='x', title='|V|^{2}', log=True, Range=(0.,1), forceRange=False, err=False):
    self.name = name
    self.axis = axis
    self.title = title
    self.log = log
    self.Range = Range
    self.forceRange = forceRange
    self.err = err


class PlotOpt(object):
  def __init__(self, what, binning, xtitle, ytitle, logX, logY):
    self.what = what
    self.binning = binning
    self.xtitle = xtitle
    self.ytitle = ytitle
    self.logX = logX
    self.logY = logY

class Sample(object):
  '''
  Handles information of a single sample, i.e. mass,ctau , w/ or w/o reweighting
  '''
  def __init__(self,mass=-99, ctau=-99, vv=-99, infileNames=None, isrw=False, orig_vv=None, label='V00'):
    self.mass = mass
    self.ctau = ctau
    self.vv = vv
    if vv==-99: 
      self.vv=getVV(mass=self.mass, ctau=self.ctau)
    if ctau==-99:
      self.ctau=getCtau(mass=self.mass, vv=self.vv)
    self.treeName='tree'
    self.infileNames=infileNames
    #if not os.path.isfile(self.infileName): raise RuntimeError('file %s not found' % self.infileName)
    self.isrw=isrw
    if self.isrw:
      self.orig_vv = orig_vv
      self.orig_ctau = getCtau(mass=self.mass,vv=orig_vv)
    else:
      self.orig_vv = self.vv
      self.orig_ctau = self.ctau
    self.label=label
    self.ngenevts = int(label.split('_n')[1])
    self.njobs = int(label.split('_njt')[1])
    #self.isMajorana
    self.name='bhnl_mass{m}_ctau{ctau}'.format(m=self.mass, ctau=self.ctau)
    self.legname='{:3}: m={:.1f}GeV |V|^{{2}}={:.1e} c#tau={:.1f}mm {}'.format('rw' if self.isrw else 'gen',self.mass,self.vv,self.ctau,'- orig |V|^{{2}}={:.1e}'.format(self.orig_vv) if self.isrw else '')
    if self.isrw:
      self.evt_w = '(weight_{vv})'.format(vv=str(self.vv).replace('-', 'm')) # reweight the tree (created with orig vv) to the vv of this sample
    else:
      self.evt_w = '(1)'
    self.acc = -99
    self.acc_errup = -99
    self.acc_errdn = -99
    self.effFilter = -99
    self.effFilter_errup = -99
    self.effFilter_errdn = -99
  
   
    self.histoDefs = {
    # b mother
    'b_pt'          : PlotOpt('b_pt', '(60,0,60)', 'B meson p_{T} [GeV]', 'a.u.', False, True),
    'b_signid'      : PlotOpt('b_pdgid/abs(b_pdgid)', '(2,-1.1,1.1)', 'B meson pdg id sign', 'a.u.', False, True),   
    'b_eta'         : PlotOpt('b_eta', '(60,-6,6)', 'B meson #eta', 'a.u.', False, True),
    'b_ct'          : PlotOpt('b_ct_reco', '(50,0,100)', 'B meson ct [mm]', 'a.u.', False, True),
    #'b_ct_large'    : PlotOpt('b_ct_reco', '(100,0,1000)', 'B meson ct [mm]', 'a.u.', False, True),     
    # daughters of the B
    ## the HNL
    'hnl_pt'        : PlotOpt('hnl_pt', '(60,0,30)', 'HNL p_{T} [GeV]', 'a.u.', False, True),   
    'hnl_signid'    : PlotOpt('hnl_pdgid/abs(hnl_pdgid)', '(2,-1.1,1.1)', 'HNL pdg id sign', 'a.u.', False, True),   
    'hnl_eta'       : PlotOpt('hnl_eta', '(60,-6,6)', 'HNL #eta', 'a.u.', False, True),      
    'hnl_ct'        : PlotOpt('hnl_ct_reco', '(50,0,1000)', 'HNL ct [mm]', 'a.u.', False, True),
    'hnl_ct_large'  : PlotOpt('hnl_ct_reco', '(100,0,10000)', 'HNL ct [mm]', 'a.u.', False, True),    
    'hnl_Lxy'       : PlotOpt('Lxy', '(50,0,1000)', 'L_{xy} [mm]', 'a.u.', False, True),
    'hnl_Lxy_large' : PlotOpt('Lxy', '(100,0,10000)', 'L_{xy} [mm]', 'a.u.', False, True),    
    'hnl_Lxyz'       : PlotOpt('Lxyz', '(50,0,1000)', 'L_{xyz} [mm]', 'a.u.', False, True),
    'hnl_Lxyz_large' : PlotOpt('Lxyz', '(100,0,10000)', 'L_{xyz} [mm]', 'a.u.', False, True),    

    ### the D meson
#    'd_pt'          : PlotOpt('d_pt', '(60,0,30)', 'D meson p_{T} [GeV]', 'a.u.', False, True),   
#    'd_eta'         : PlotOpt('d_eta', '(60,-6,6)', 'D meson #eta', 'a.u.', False, True),      

    ### the trigger lepton
    'mutrig_pt'     : PlotOpt('l0_pt', '(60,0,30)', '#mu^{trig} p_{T} [GeV]', 'a.u.', False, True),   
    'mutrig_eta'    : PlotOpt('l0_eta', '(60,-6,6)', '#mu^{trig} #eta', 'a.u.', False, True),      

    #### daughters of the D meson
    ###### the pion
#    'piD_pt'        : PlotOpt('pi_pt', '(60,0,30)', '#pi (from D) p_{T} [GeV]', 'a.u.', False, True),  
#    'piD_eta'       : PlotOpt('pi_eta', '(60,-6,6)', '#pi (from D) #eta', 'a.u.', False, True),      
    
    ###### the kaon
#    'k_pt'          : PlotOpt('k_pt', '(60,0,30)', 'K (from D) p_{T} [GeV]', 'a.u.', False, True),  
#    'k_eta'         : PlotOpt('k_eta', '(60,-6,6)', 'K (from D) #eta', 'a.u.', False, True),      
    
    #### daughters of the HNL
    ###### the lepton
    'mu_pt'         : PlotOpt('l1_pt', '(60,0,30)', '#mu (from HNL) p_{T} [GeV]', 'a.u.', False, True),  
    'mu_eta'        : PlotOpt('l1_eta', '(60,-6,6)', '#mu (from HNL) #eta', 'a.u.', False, True),      

    ###### the pion
    'pi_pt'         : PlotOpt('pi1_pt', '(60,0,30)', '#pi (from HNL) p_{T} [GeV]', 'a.u.', False, True),   
    'pi_eta'        : PlotOpt('pi1_eta', '(60,-6,6)', '#pi (from HNL) #eta', 'a.u.', False, True),      
   
    # invariant masses
    'hnl_invmass'   : PlotOpt('lep_pi_invmass', '(50,0,5)', 'HNL invariant mass, m(#mu,#pi) [GeV]', 'a.u.', False, False),     
#    'd_invmass'     : PlotOpt('k_pi_invmass', '(50,0,5)', 'D meson invariant mass, m(K,#pi) [GeV]', 'a.u.', False, False),      
    #'b_invmass'     : PlotOpt('b_invmass', '(50,2,7)', 'B invariant mass [GeV]', 'a.u.', False, False),     
    'bpartial_invmass'     : PlotOpt('bpartial_invmass', '(50,2,7)', 'm(HNL,#mu^{trig}) [GeV]', 'a.u.', False, False),     
    'bpartial_invmass_log'     : PlotOpt('bpartial_invmass', '(50,2,7)', 'm(HNL,#mu^{trig}) [GeV]', 'a.u.', False, True),     
    #'Lxy_cos', # cosine of the pointing angle in the transverse plane
    #'Lxyz_b', #3D displacement of the B wrt to primary vertex
    #'Lxyz_l0' #3D displacement of the prompt lepton wrt to B vertex
    }
    self.histos = {}
   
  def stamp(self):
    '''
    This should print the basic information of the sample
    '''
    if debug: print('Info in Sample.stamp()')
    print('mass={m}GeV, ctau={ctau}mm VV={vv}, isrw={isrw}, orig_VV={ovv}, acc={acc}, au={au}, ad={ad}, num={num}, den={den}, evt_w={ew} effFilter={ef}'.format( \
            m=self.mass,ctau=self.ctau,vv=self.vv,isrw=self.isrw,ovv=self.orig_vv,acc=self.acc,au=self.acc_errup,ad=self.acc_errdn,num=self.num,den=self.den,ew=self.evt_w,ef=self.effFilter))
  
  def fillHistos(self,sel='(1)',sellabel='noSel'):
    '''
    This is to fill the histograms
    '''
    if debug: print('Info in Sample.fillHistos()')
    chain = TChain(self.treeName)
    for fn in self.infileNames:
      chain.Add(fn)
    #f = TFile.Open(self.infileName)
    #t = f.Get(self.treeName)
    #if not t:
    #  print 'ERROR: no tree in file %s' % fname
    #  return False

    for name,spec in self.histoDefs.items():
      chain.Draw('%s>>%s%s' % (spec.what, name, spec.binning), sel+'*'+self.evt_w, 'goff') 
      #chain.Project('{}'.format(name),'{}'.format(spec.what))
      h = chain.GetHistogram().Clone('%s_%s_%s' % (self.name, name, sellabel))
      h.SetDirectory(0)
      h.SetTitle(h.GetName())
      h.GetXaxis().SetTitle(spec.xtitle)
      h.GetXaxis().SetNdivisions(505)
      h.GetYaxis().SetTitle(spec.ytitle)
      h.GetYaxis().SetNdivisions(505)
      #h.GetYaxis().SetTitleSize(0.05)
      #h.GetXaxis().SetTitleSize(0.05)

      #newh = getOverflowedHisto(h)
      #newh.SetDirectory(0)
      #self.histos[name] = newh
      h.SetDirectory(0)
      self.histos[name] = h

    #f.IsA().Destructor(f)
    return True

  def saveHistos(self, norm):
    '''
    This is to save the histograms of the sample in a plot
    '''
    if debug: print('Info in Sample.saveHistos()')
    c = TCanvas(self.name, self.name, 800, 600)
    c.DivideSquare(len(self.histos.keys()))

    for i, what in enumerate(self.histos.keys()):
      h = self.histos[what]
      pad = c.cd(i + 1)
      pad.SetLogx(self.histoDefs[what].logX)
      pad.SetLogy(self.histoDefs[what].logY)
      if norm and h.Integral() != 0:
        hdrawn = h.DrawNormalized('LPE')
      else:
        hdrawn = h.Draw('PLE')

      norm_suffix='_norm' if norm else ''

    c.SaveAs('./plots/' + self.label + suffix + '/' + c.GetName() + norm_suffix + '.png')
    c.SaveAs('./plots/' + self.label + suffix + '/' + c.GetName() + norm_suffix + '.pdf')
    c.SaveAs('./plots/' + self.label + suffix + '/' + c.GetName() + norm_suffix + '.C')

  def fillAcceptance(self):
    '''
    This is to calculate the acceptance of the sample
    '''
    if debug: print('Info in Sample.fillAcceptance()')
    #f = TFile.Open(self.infileName)
    #t = f.Get(self.treeName)
    #if not t:
    #  raise RuntimeError( 'ERROR: no tree in file %s' % self.infileName)
    chain = TChain(self.treeName)
    for fn in self.infileNames:
      chain.Add(fn)

    self.effnum = ROOT.TH1F('effnum', 'effnum', 1, 0, 13000) #dict([(k, ROOT.TH1F('effnum_%s' % k, 'effnum_%s' % k, 1, 0, 1)) for k in self.settings])
    self.effden = ROOT.TH1F('effden', 'effden', 1, 0, 13000)

    if doInclusive:
      cutsnum = '(l0_pt>{mp} && abs(l0_eta)<1.5'.format(mp=muTrigPt)
    else: 
      cutsnum = '(l0_pt>{mp} && abs(l0_eta)<1.5 && k_pt>1 && abs(k_eta)<2.5 && pi_pt>1 && abs(pi_eta)<2.5'.format(mp=muTrigPt) 
    if doSkipDispl:
      cutsnum += ''
    else:
      cutsnum += '&& Lxy < 1000'
    if doSkipHNLptEta:
      cutsnum += ''
    else:
      cutsnum += '&& l1_pt>3 && abs(l1_eta)<2.5 && pi1_pt>0.8 && abs(pi1_eta)<2.5'
    
    cutsnum += ')'

    cutsden = '(l0_pt>{mp} && abs(l0_eta)<1.5)'.format(mp=muTrigPt)
    

    chain.Draw('hnl_pt>>effnum', cutsnum+'*'+self.evt_w, 'goff')
    chain.Draw('hnl_pt>>effden', cutsden+'*'+self.evt_w, 'goff')
 
    if TEfficiency.CheckConsistency(self.effnum,self.effden): 
      peff = TEfficiency(self.effnum,self.effden)
      
      self.acc = peff.GetEfficiency(1)
      self.acc_errup = peff.GetEfficiencyErrorUp(1)
      self.acc_errdn = peff.GetEfficiencyErrorLow(1)

    # for debugging purposes
    self.num = self.effnum.GetEntries()
    if self.num==0: print('**** 0 entries for mass={}'.format(self.mass))
    self.den = self.effden.GetEntries()

  def fillExpNevts(self):
    if debug: print('Info in Sample.fillExpNevts()')
    N_B = 3.84E9 * 2. # 3.84E9, inclusive, = n(B+/-) * 2  
    # factor ~2 is for considering B0s as well, which are ~50% of the B-parking dataset

    if doInclusive:
      corr = ( BR__B_l_N(mass=self.mass) + BR__B_D0_l_N(mass=self.mass) ) / ( BR__B_l_nu() + BR__B_D0_l_nu() )
      N_HNL_VV1 = N_B * corr
    else:
      corr = BR__B_D0_l_N(mass=self.mass) / BR__B_D0_l_nu() 
      N_HNL_VV1 = N_B * BR__B_D0_l_nu() / BR__B_X_l_nu() * corr * BR__D0_K_pi() 
          
    self.expNevts = N_HNL_VV1 * BR_HNLmupion(mass=self.mass) * self.vv * self.acc
    self.expNevts_errup = N_HNL_VV1 * BR_HNLmupion(mass=self.mass) * self.vv * self.acc_errup
    self.expNevts_errdn = N_HNL_VV1 * BR_HNLmupion(mass=self.mass) * self.vv * self.acc_errdn

  def fillFilterEff(self,dostamp=True):
    '''
    To retrieve and save the filter efficiency - from the minigentree
    TODO: retrieve the cross-section => for that you would need to run without separate jobs
    '''
    ##if debug: print('Info in Sample.fillFilterEff()')

    efffnum = ROOT.TH1F('efffnum', 'efffnum', 1, 0, 13000)
    efffden = ROOT.TH1F('efffden', 'efffden', 1, 0, 13000)

    chain = TChain(self.treeName)
    for fn in self.infileNames:
      chain.Add(fn)
    
    efffnum.SetBinContent(1, chain.GetEntries())
    efffnum.SetBinError(1, ROOT.TMath.Sqrt(chain.GetEntries()))
    
    # denominator = number of events that were run in the first place # access storage element... 
    #self.ngenevts_succ_afterfilter=0
    self.njobs_succ=0
    path = '/pnfs/psi.ch/cms/trivcat/store/user/{u}/BHNLsGen/{pl}/mass{m}_ctau{ctau}/step1*root'.format(u=os.environ['USER'],pl=self.label,m=self.mass,ctau=self.orig_ctau)
    for fname in glob(path):
        #if debug: print 'fname=',fname
        f = TFile.Open(fname)
        if f.GetListOfKeys().Contains('Events'):
          self.njobs_succ += 1

    #njobs_succ = len(glob(path))
    self.ngenevts_succ = float(self.ngenevts) / float(self.njobs) * float(self.njobs_succ)
    
    efffden.SetBinContent(1,self.ngenevts_succ)
    efffden.SetBinError(1,ROOT.TMath.Sqrt(self.ngenevts_succ))

    if TEfficiency.CheckConsistency(efffnum,efffden): 
      geneff = TEfficiency(efffnum,efffden) 
      self.effFilter = geneff.GetEfficiency(1)
      self.effFilter_errup = geneff.GetEfficiencyErrorUp(1)
      self.effFilter_errdn = geneff.GetEfficiencyErrorLow(1)

    # stamp basic info about filter eff
    if dostamp: print('mass={m}GeV, VV={vv:.1e}, nTot={nt}, nSucc={ns}, nSuccAfterFilter={nsf}, effFilter={ef:.3f}%, errup={eu:.3f}%, errdn={ed:.3f}% '.format( \
                m=self.mass,vv=self.vv,nt=self.ngenevts,ns=self.ngenevts_succ,nsf=efffnum.GetBinContent(1),ef=self.effFilter*100,eu=self.effFilter_errup*100,ed=self.effFilter_errdn*100))

class SampleList(object):
  '''
  Handles the plotting of several samples, both the histograms, and the general quantities (graphs)
  '''
  def __init__(self, name, samples, label=None):
    self.name = name
    self.samples = samples  # a list of objects of type Sample
    self.label = samples[0].label if label==None else label
    self.colors = [ ROOT.kBlack,  ROOT.kOrange+1, ROOT.kRed, ROOT.kMagenta+2, ROOT.kViolet+8, ROOT.kAzure-8, ROOT.kAzure+6 ,
                      ROOT.kGreen+1, ROOT.kSpring+4, ROOT.kYellow -5, ROOT.kYellow -3, ROOT.kYellow, ROOT.kOrange
                  ]
    self.styles = [  1,1,1,1,1,1,1,1,1,1,1,1,1,1, ]
    self.graphs={}

    # x,y quantities for graphs
    self.quantities={
      #'acc' : SampleQuantity(name='acc', axis='y', title='Acceptance', log=True, Range=(0.,1), forceRange=True, err=True),
      'acc' : SampleQuantity(name='acc', axis='y', title='Acceptance', log=True, Range=(0.00001,2), forceRange=True, err=True),
      #'acc' : SampleQuantity(name='acc', axis='y', title='Acceptance', log=False, Range=(0.,1), forceRange=True, err=True),
      'expNevts': SampleQuantity(name='expNevts', axis='y', title='N_{#nu} x BR(HN#rightarrow#mu#pi) x Acc x V^{2}', log=True, err=True),
      'num' : SampleQuantity(name='num', axis='y', title='Numerator', log=True, Range=(1, 3000), forceRange=True, err=False),
      'den' : SampleQuantity(name='den', axis='y', title='Denominator', log=True, Range=(1, 3000), forceRange=True, err=False),
      #'filterEff' : 
      # 'xsec'
      'vv'  : SampleQuantity(name='vv', axis='x', title='|V|^{2}', log=True),
      'ctau': SampleQuantity(name='ctau', axis='x', title='c#tau [mm]', log=False),
      'mass': SampleQuantity(name='mass', axis='x', title='mass [GeV]', log=False),
    } 

    if len(self.samples) > len(self.colors):
      raise RuntimeError('need more colors (%d vs %d)' % (len(self.colors), len(self.samples)))
    self.checkConsistency()
  def add(self, sample):
    self.samples.insert(sample)
    if len(self.samples) > len(self.colors):
      raise RuntimeError('need more colors (%d vs %d)' % (len(self.colors), len(self.samples)))
    self.checkConsistency()
  def remove(self, sample):
    self.samples.remove(sample)
  def checkConsistency(self):
    '''
    Makes sure histograms of samples have same binning and histograms
    '''
    for sample1 in self.samples:
      for sample2 in self.samples:
        if len(sample1.histoDefs) != len(sample2.histoDefs):
          raise RuntimeError('number of PlotOpts differs for %s and %s' % (sample1.name, sample2.name))
        if len(sample1.histos) != len(sample2.histos):
          raise RuntimeError('number of histos differs for %s and %s' % (sample1.name, sample2.name))

  def plotHistos(self, norm=False, sameCanvas=False):
    ''' 
    Superimpose histograms for all samples, on the same canvas or on separate canvases
    '''
    if sameCanvas:
      c = TCanvas(self.name, self.name, 6400, 4800)
      tosave = { c.GetName() : c }
    else:
      c = None
      tosave = {}

    legends = []

    hMax = {}
    for j, what in enumerate(self.samples[0].histos.keys()):
      hMax[what] = []
      leg=defaultLegend(x1=0.15,y1=0.7,x2=0.95,y2=0.95,mult=1.2)
      legends.append(leg)
    # do the actual plotting
    for i, sample in enumerate(self.samples):
      if i == 0:
        opt = 'histE'
        if sameCanvas: c.DivideSquare(len(sample.histos.keys()))
      else: opt = 'histEsame'
      for j, what in enumerate(sample.histos.keys()):
        h = sample.histos[what]
        graph_saver.append(h)

        if sameCanvas: pad = c.cd(j + 1)
        else:
          cname = '%s_%s' % (self.name, what)
          if cname not in tosave.keys():
            tosave[cname] = TCanvas(cname, cname, 700,600)
          pad = tosave[cname].cd()

        if i == 0:
          pad.SetLogx(sample.histoDefs[what].logX)
          pad.SetLogy(sample.histoDefs[what].logY)
        if norm and h.Integral():
          hdrawn = h.DrawNormalized(opt)
          hMax[what].append(hdrawn)
          hdrawn.SetLineColor(self.colors[i])
          hdrawn.SetMarkerColor(self.colors[i])
          hdrawn.SetMarkerStyle(self.styles[i])
          hdrawn.SetMarkerSize(3)
          hdrawn.SetLineWidth(2)
          legends[j].AddEntry(hdrawn,sample.legname, 'LP')

        else:
          h.Draw(opt)
          hMax[what].append(h)
          h.SetLineColor(self.colors[i])
          h.SetMarkerColor(self.colors[i])
          h.SetMarkerStyle(self.styles[i])
          h.SetMarkerSize(4)
          h.SetLineWidth(2)
          legends[j].AddEntry(hdrawn,sample.legname, 'LP')

        if i == len(self.samples) - 1:
          legends[j].Draw('same')

    # compute the max of norm histo and set the range accordingly
    for j, what in enumerate(self.samples[0].histos.keys()):
      if len(hMax[what])>0:
        max = 0
        for ih in hMax[what]:
          if(max < ih.GetMaximum()):
            max=ih.GetMaximum()
        if (sample.histoDefs[what].logY==True):
          #hMax[name][0].GetYaxis().SetRangeUser(0.01, max * 7)
          hMax[what][0].SetMaximum(max*30)
          #print hMax[what][0].GetMaximum()
        else:
          #hMax[name][0].GetYaxis().SetRangeUser(0.01, max * 1.5)
          hMax[what][0].SetMaximum(max * 1.7)

    norm_suffix='_norm' if norm else ''

    for cname, canv in tosave.items():
      canv.SaveAs('./plots/' + self.label + suffix + '/' + cname + norm_suffix + '.png')
      canv.SaveAs('./plots/' + self.label + suffix + '/' + cname + norm_suffix + '.pdf')
      canv.SaveAs('./plots/' + self.label + suffix + '/' + cname + norm_suffix + '.C')

  def plotGraph(self, x='vv', y='acc'): 
    '''
    Plot a graph with specified quantities on x and y axes , and saves the graph
    '''

    if (x not in self.quantities.keys()) or (y not in self.quantities.keys()):
      raise RuntimeError('selected quantities not available, available quantities are: \n{}'.format(self.quantities.keys()))

    xq = self.quantities[x]
    yq = self.quantities[y]

    #graph = TGraphAsymmErrors()
    #graph = TGraph()
    graph = TGraphErrors()
    for i,s in enumerate(self.samples):
      graph.SetPoint(i,getattr(s, xq.name), getattr(s, yq.name) )
      #if xq.err: 
      #  graph.SetPointEXhigh(i, getattr(s, xq.name+'_errup'))   # errup errdn
      #  graph.SetPointEXlow (i, getattr(s, xq.name+'_errdn'))
      if yq.err: 
        graph.SetPointError(i, 0, getattr(s, yq.name+'_errup'))
      #  graph.SetPointEYhigh(i, getattr(s, yq.name+'_errup'))  
      #  graph.SetPointEYlow (i, getattr(s, yq.name+'_errdn'))

    c = TCanvas()
    graph.SetLineWidth(2)
    graph.SetMarkerStyle(22)
    graph.SetTitle(';{x};{y}'.format(y=yq.title,x=xq.title))
    graph.Draw('APLE')

    if yq.forceRange:
      graph.SetMinimum(yq.Range[0])
      graph.SetMaximum(yq.Range[1])

    gPad.Modified()
    gPad.Update()
    if yq.name=='expNevts':
      line = TLine(gPad.GetUxmin(),3,gPad.GetUxmax(),3)
      line.SetLineColor(ROOT.kBlue)
      line.Draw('same')
      #graph.SetMinimum(0.01)
      #graph.SetMaximum(1E06)

    if xq.log: c.SetLogx()
    if yq.log: c.SetLogy()
    c.SetGridx()
    c.SetGridy()
    c.SaveAs('./plots/{}{}/{}_{}VS{}.pdf'.format(self.label,suffix,self.name,yq.name,xq.name))
    c.SaveAs('./plots/{}{}/{}_{}VS{}.C'.format(self.label,suffix,self.name,yq.name,xq.name))
    c.SaveAs('./plots/{}{}/{}_{}VS{}.png'.format(self.label,suffix,self.name,yq.name,xq.name))

    self.graphs['{}VS{}'.format(yq.name,xq.name)] = graph
    # add the graph container for memory?
    graph_saver.append(graph)

def doAnalysis(path,pl,points,name,path2=None,pl2=None,points2=None,path3=None,pl3=None,points3=None):
  '''
  Perform plotting of samples lists
  Optionally allows to create sample list out of 2 different productions
  TODO: quantities to plot in graph as an option...
  '''
  print('\n*************************************************')
  print('  => Going to do plotting for production={}, name={}'.format(pl,name))
  samples = []
  for p in points:
    fns = glob(path.format(m=p.mass,ctau=p.orig_ctau))
    s = Sample(mass=p.mass, ctau=p.ctau, vv=p.vv, infileNames=fns, isrw=p.isrw, orig_vv=p.orig_vv, label=pl)
    s.fillHistos() # no selection
    #s.fillHistos(sel='(l0_pt>5 && abs(l0_eta)<1.6)', sellabel='pt5_eta1p6') # filter selection
    s.fillAcceptance()
    s.fillExpNevts()
    s.fillFilterEff()
    s.stamp()
    samples.append(s)

  if path2 and pl2 and points2:
    print('\n*************************************************')
    print('  => Going to do plotting for production={}, name={}'.format(pl2,name))
    for p in points2:
      fns = glob(path2.format(m=p.mass,ctau=p.orig_ctau))
      s = Sample(mass=p.mass, ctau=p.ctau, vv=p.vv, infileNames=fns, isrw=p.isrw, orig_vv=p.orig_vv, label=pl2)
      s.fillHistos()
      #s.fillHistos(sel='(l0_pt>5 && abs(l0_eta)<1.6)', sellabel='pt5_eta1p6') # filter selection
      s.fillAcceptance()
      s.fillExpNevts()
      s.fillFilterEff()
      s.stamp()
      samples.append(s)

  if path3 and pl3 and points3:
    print('\n*************************************************')
    print('  => Going to do plotting for production={}, name={}'.format(pl3,name))
    for p in points3:
      fns = glob(path3.format(m=p.mass,ctau=p.orig_ctau))
      s = Sample(mass=p.mass, ctau=p.ctau, vv=p.vv, infileNames=fns, isrw=p.isrw, orig_vv=p.orig_vv, label=pl3)
      s.fillHistos()
      #s.fillHistos(sel='(l0_pt>5 && abs(l0_eta)<1.6)', sellabel='pt5_eta1p6') # filter selection
      s.fillAcceptance()
      s.fillExpNevts()
      s.fillFilterEff()
      s.stamp()
      samples.append(s)


  label = pl if not pl2 else ('{}_VS_{}'.format(pl,pl2) if not pl3 else '{}_VS_{}_VS_{}'.format(pl,pl2,pl3))

  slist = SampleList(name=name, samples=samples, label=label)
  slist.plotHistos(norm=True, sameCanvas=False)
  if 'fixedM' in name or 'closure' in name: 
    slist.plotGraph(x='vv',y='acc')
    slist.plotGraph(x='vv',y='expNevts')
  elif 'fixedVV' in name:
    slist.plotGraph(x='mass',y='acc')
    slist.plotGraph(x='mass',y='num')
    slist.plotGraph(x='mass',y='den')
    ###slist.plotGraph(x='mass',y='expNevts')

  return slist


def doGraphComparison(list1,list2,what):

  c = TCanvas()
  list1.graphs[what].Draw('APL')
  list2.graphs[what].Draw('same')
  c.SaveAs(what+'.pdf')
  c.SaveAs(what+'.png')

def doLimitGraph(slists,what, label):

  c = TCanvas()
  graph = TGraph()
  for ilist,slist in enumerate(slists):
    # get the inverted graph
    inverted_graph = TGraph()
    xs = slist.graphs[what].GetX()
    ys = slist.graphs[what].GetY()
    N =  slist.graphs[what].GetN()
    xs_arr = np.ndarray(N,'d',xs) 
    ys_arr = np.ndarray(N,'d',ys) 
    inverted_graph = TGraph(N,ys,xs)
    # use Eval
    limit = inverted_graph.Eval(3.) # linear interpolation between closest points
    graph.SetPoint(ilist,slist.samples[0].mass,limit)
  graph.SetLineWidth(2)
  if doInclusive: graph.SetMarkerStyle(22)
  else: graph.SetMarkerStyle(24)
  graph.SetTitle(';mass [GeV]; exp. limit on V^{2}')
  graph.SetMaximum(1E-02)
  graph.SetMinimum(1E-08)
  graph.Draw('APL')
  graph.GetXaxis().SetLimits(0.1, 10)
  graph.Draw('APL')
  c.Update()

  c.SetLogx() 
  c.SetLogy()
  #c.SetGridx()
  c.SetGridy()
  c.SaveAs('./plots/{}{}/limit_vv.pdf'.format(label,suffix))

def checkFiles(path,points):
  '''
  Check existence of files for analysis, return existing points
  '''
  for p in points:
    expr = path.format(m=p.mass,ctau=p.orig_ctau) # check the file where the tree is stored
    nfiles = 0
    for fn in glob(expr):
      if os.path.isfile(fn):
        nfiles+=1
    if nfiles==0:
      p.missing=True
    else: 
      p.missing=False
 
  existing_points = [p for p in points if not p.missing]
  
  if len(existing_points)>0:
    return existing_points 
  else:
    raise RuntimeError('no files for analysis')
    return []


  print('  => Going to do analysis on following points:')
  print('  ')

def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='', add_help=True)
  parser.add_argument('--pl', type=str, dest='pl', help='production label of input', default='V_blabla')  
  parser.add_argument('--pl2', type=str, dest='pl2', help='production label of second input', default=None)  
  parser.add_argument('--pl3', type=str, dest='pl3', help='production label of second input', default=None)  
  return parser.parse_args()

if __name__ == "__main__":

  gROOT.SetBatch(True)
  ROOT.TH1.SetDefaultSumw2()
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  gStyle.SetTitleXOffset(1.1);
  gStyle.SetTitleYOffset(1.45);

  global debug
  global doInclusive
  global doSkipDispl
  global suffix
  global muTrigPt

  #### options
  debug = False
  doInclusive = True # 
  doSkipDispl = False #
  doSkipHNLptEta = False
  doCompareAnalysis = False #
  doTestAnalysis = False
  doFixedMassAnalysis = False
  doRwAnalysis = False
  doFixedVVAnalysis = True
  muTrigPt = 9 # 0 1 2 5 7 9
  ####

  suffix='_incl' if doInclusive else '_excl'
  if doSkipDispl: suffix += '_skipDispl'
  if doSkipHNLptEta: suffix += '_skipHNLptEta'
  suffix += '_muTrigPt{mp}'.format(mp=muTrigPt)
  opt = getOptions()
  
  complabel = '' if not opt.pl2 else ('_VS_{}'.format(opt.pl2) if not opt.pl3 else '_VS_{}_VS_{}'.format(opt.pl2,opt.pl3))
  os.system('mkdir -p ./plots/{}{}{}'.format(opt.pl,complabel,suffix))

  path = './outputfiles/' + opt.pl + '/mass{m}_ctau{ctau}_miniGenTree*.root'
  #path = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGenDump/' + opt.pl + "/mass{m}_ctau{ctau}"  + '/mass{m}_ctau{ctau}_miniGenTree*.root'
  if opt.pl2: 
    path2 = './outputfiles/' + opt.pl2 + '/mass{m}_ctau{ctau}_miniGenTree*.root'
    #path2 = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGenDump/' + opt.pl2 + "/mass{m}_ctau{ctau}"  + '/mass{m}_ctau{ctau}_miniGenTree*.root'
  if opt.pl3: 
    path3 = './outputfiles/' + opt.pl3 + '/mass{m}_ctau{ctau}_miniGenTree*.root'
    #path3 = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGenDump/' + opt.pl3 + "/mass{m}_ctau{ctau}"  + '/mass{m}_ctau{ctau}_miniGenTree*.root'
  

  if doCompareAnalysis:
    points = [Point(mass=1.5,ctau=None,vv=1e-03,isrw=False)]
    points2 = [Point(mass=1.5,ctau=None,vv=1e-03,isrw=False)]
    points3 = [Point(mass=1.5,ctau=None,vv=1e-03,isrw=False)]  

    for p in points:
      p.stamp()
    existing_points=checkFiles(path=path,points=points)
    for p in points2:
      p.stamp()
    existing_points2=checkFiles(path=path2,points=points2)
    if not opt.pl3: 
      doAnalysis(path=path,pl=opt.pl,points=existing_points,name='comp_norw',path2=path2,pl2=opt.pl2,points2=existing_points2)
    else:
      for p in points3:
        p.stamp()
      existing_points3=checkFiles(path=path3,points=points3)
      doAnalysis(path=path,pl=opt.pl,points=existing_points,name='comp_norw',path2=path2,pl2=opt.pl2,points2=existing_points2,path3=path3,pl3=opt.pl3,points3=existing_points3)

  if doTestAnalysis:
    #points = [Point(mass=1.0,ctau=None,vv=5e-04,isrw=False)]
    points = [Point(mass=1.5,ctau=None,vv=1e-03,isrw=False)]
    #points = [Point(mass=3.0,ctau=None,vv=5e-05,isrw=False)]
    for p in points:
      p.stamp()
    existing_points=checkFiles(path=path,points=points)
    doAnalysis(path=path,pl=opt.pl,points=existing_points,name='testpoint_norw')

  if doFixedMassAnalysis:
   
    slists_fixedMass = []

    ################
    points = [
      Point(mass=3.0,ctau=None,vv=5e-03,isrw=False),
      Point(mass=3.0,ctau=None,vv=1e-03,isrw=False),
      Point(mass=3.0,ctau=None,vv=5e-04,isrw=False),
      Point(mass=3.0,ctau=None,vv=1e-04,isrw=False),
      Point(mass=3.0,ctau=None,vv=5e-05,isrw=False),
      Point(mass=3.0,ctau=None,vv=1e-05,isrw=False),
      Point(mass=3.0,ctau=None,vv=5e-06,isrw=False),
      Point(mass=3.0,ctau=None,vv=1e-06,isrw=False),
      Point(mass=3.0,ctau=None,vv=5e-07,isrw=False),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass3.0_norw'))
    ################
    points = [
      Point(mass=2.5,ctau=None,vv=5e-03,isrw=False),
      Point(mass=2.5,ctau=None,vv=1e-03,isrw=False),
      Point(mass=2.5,ctau=None,vv=5e-04,isrw=False),
      Point(mass=2.5,ctau=None,vv=1e-04,isrw=False),
      Point(mass=2.5,ctau=None,vv=5e-05,isrw=False),
      Point(mass=2.5,ctau=None,vv=1e-05,isrw=False),
      Point(mass=2.5,ctau=None,vv=5e-06,isrw=False),
      Point(mass=2.5,ctau=None,vv=1e-06,isrw=False),
      Point(mass=2.5,ctau=None,vv=5e-07,isrw=False),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass2.5_norw'))
    ################
    points = [
      Point(mass=2.0,ctau=None,vv=5e-03,isrw=False),
      Point(mass=2.0,ctau=None,vv=1e-03,isrw=False),
      Point(mass=2.0,ctau=None,vv=5e-04,isrw=False),
      Point(mass=2.0,ctau=None,vv=1e-04,isrw=False),
      Point(mass=2.0,ctau=None,vv=5e-05,isrw=False),
      Point(mass=2.0,ctau=None,vv=1e-05,isrw=False),
      Point(mass=2.0,ctau=None,vv=5e-06,isrw=False),
      Point(mass=2.0,ctau=None,vv=1e-06,isrw=False),
      Point(mass=2.0,ctau=None,vv=5e-07,isrw=False),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass2.0_norw'))
    ################
    points = [
      Point(mass=1.5,ctau=None,vv=5e-03,isrw=False),
      Point(mass=1.5,ctau=None,vv=1e-03,isrw=False),
      Point(mass=1.5,ctau=None,vv=5e-04,isrw=False),
      Point(mass=1.5,ctau=None,vv=1e-04,isrw=False),
      Point(mass=1.5,ctau=None,vv=5e-05,isrw=False),
      Point(mass=1.5,ctau=None,vv=1e-05,isrw=False),
      Point(mass=1.5,ctau=None,vv=5e-06,isrw=False),
      Point(mass=1.5,ctau=None,vv=1e-06,isrw=False),
      Point(mass=1.5,ctau=None,vv=5e-07,isrw=False),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass1.5_norw'))
    ################
    points = [
      Point(mass=1.0,ctau=None,vv=5e-03,isrw=False),
      Point(mass=1.0,ctau=None,vv=1e-03,isrw=False),
      Point(mass=1.0,ctau=None,vv=5e-04,isrw=False),
      Point(mass=1.0,ctau=None,vv=1e-04,isrw=False),
      Point(mass=1.0,ctau=None,vv=5e-05,isrw=True,orig_vv=1e-04),
      Point(mass=1.0,ctau=None,vv=1e-05,isrw=False),
      Point(mass=1.0,ctau=None,vv=5e-06,isrw=False),
      Point(mass=1.0,ctau=None,vv=1e-06,isrw=False),
      Point(mass=1.0,ctau=None,vv=5e-07,isrw=False),
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slist_norw = doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass1.0_norw')
    slists_fixedMass.append(slist_norw)
    ################
    points = [
      Point(mass=0.5,ctau=None,vv=5e-03,isrw=False),
      Point(mass=0.5,ctau=None,vv=1e-03,isrw=False),
      Point(mass=0.5,ctau=None,vv=5e-04,isrw=False),
      Point(mass=0.5,ctau=None,vv=1e-04,isrw=False),
      Point(mass=0.5,ctau=None,vv=5e-05,isrw=False),
      Point(mass=0.5,ctau=None,vv=1e-05,isrw=False),
      Point(mass=0.5,ctau=None,vv=5e-06,isrw=False),
      Point(mass=0.5,ctau=None,vv=1e-06,isrw=False),
      Point(mass=0.5,ctau=None,vv=5e-07,isrw=False),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass0.5_norw'))

    doLimitGraph(slists_fixedMass,'expNevtsVSvv', label=opt.pl)

    if doRwAnalysis:
      ################
      points = [
        Point(mass=1.0,ctau=None,vv=5e-03,isrw=False),
        Point(mass=1.0,ctau=None,vv=1e-03,isrw=True,orig_vv=5e-03),
        Point(mass=1.0,ctau=None,vv=5e-04,isrw=True,orig_vv=5e-03),
        Point(mass=1.0,ctau=None,vv=1e-04,isrw=True,orig_vv=5e-03),
        Point(mass=1.0,ctau=None,vv=5e-05,isrw=True,orig_vv=5e-03),
        Point(mass=1.0,ctau=None,vv=1e-05,isrw=True,orig_vv=5e-03),
        Point(mass=1.0,ctau=None,vv=5e-06,isrw=True,orig_vv=5e-03),
      ]
      for p in points:
       p.stamp()
      existing_points=checkFiles(path=path,points=points)
      slist_rw = doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass1.0_rwFrom5em03')
    
      ## compare the graphs w/ and w/o reweighting
      #doGraphComparison(slist_norw,slist_rw,what='accVSvv')
      #doGraphComparison(slist_norw,slist_rw,what='expNevtsVSvv')
  
      ###############
      points = [
        Point(mass=1.0,ctau=None,vv=1e-03, isrw=False),
        Point(mass=1.0,ctau=None,vv=1e-03, isrw=True, orig_vv=5e-03),
        Point(mass=1.0,ctau=None,vv=5e-03, isrw=False),
      ]
      existing_points = checkFiles(path,points)
      doAnalysis(path=path,pl=opt.pl,points=existing_points,name='closureRw_Mass1.0')


  if doFixedVVAnalysis:

    ################
    points = [
#      Point(mass=0.5,ctau=None,vv=1e-04,isrw=False),
#      Point(mass=1.0,ctau=None,vv=1e-04,isrw=False),
#      Point(mass=1.5,ctau=None,vv=1e-04,isrw=False),
#      Point(mass=2.0,ctau=None,vv=1e-04,isrw=False),
#      Point(mass=2.5,ctau=None,vv=1e-04,isrw=False),
#      Point(mass=3.0,ctau=None,vv=1e-04,isrw=False),
      Point(mass=1.0,ctau=None,vv=1e-05,isrw=False),
      Point(mass=2.0,ctau=None,vv=1e-05,isrw=False),
      Point(mass=3.0,ctau=None,vv=1e-05,isrw=False),
      #Point(mass=3.0,ctau=None,vv=1e-05,orig_vv=5e-05,isrw=True),
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedVV1em05')
  
