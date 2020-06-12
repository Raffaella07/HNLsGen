import os
import sys
import numpy as np
from ROOT import TTree, TFile, TH1F, TEfficiency, TGraph, TCanvas, gROOT, TAxis, TMath, TLegend, gStyle, gPad, TLine
import ROOT
sys.path.append('/work/mratti/plotting/myplotting')
from spares import *


# couplings to be tested, for which the reweight is run

from genTreeProducer import new_vvs
#small_new_vvs = new_vvs[::4]  # only select one every three elements
#small_new_vvs.reverse()
new_vvs.reverse()
small_new_vvs = new_vvs[0:10]

sys.path.append('../python/') 
from common import getVV,getCtau

def getOverflowedHisto(h):
  htemp = h.Clone()
  nbins = htemp.GetNbinsX()
  last_plus_overflow = htemp.GetBinContent(nbins) + htemp.GetBinContent(nbins+1)
  last_plus_overflow_error = TMath.Sqrt( htemp.GetBinError(nbins)*htemp.GetBinError(nbins)  + htemp.GetBinError(nbins+1)*htemp.GetBinError(nbins+1))
  htemp.SetBinContent(nbins,last_plus_overflow )
  htemp.SetBinError(nbins,last_plus_overflow_error)
  return htemp

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
  def __init__(self,mass=-99, ctau=-99, vv=-99, infileName=None, isrw=False, mass_rw, ctau_rw):
    self.mass = mass
    self.ctau = ctau
    if not vv: 
      self.ctau=ctau 
      self.vv=getVV(mass=self.mass, ctau=self.ctau)
    if not ctau:
      self.vv = vv
      self.ctau=getCtau(mass=self.mass, vv=self.vv)
    self.treeName='tree'
    self.infileName=infileName
    if os.path.isfile(self.infileName): raise RuntimeError('file %s not found' % s.infileName)
    self.isrw=isrw
    #self.isMajorana
    self.nickname='bhnl_mass{m}_ctau{ctau}'.format(m=self.mass, ctau=self.ctau)
    self.legname='{}: m={:.1f}GeV |V|^{{2}}={:.1e} ctau={:.1f}mm'.format('rw' if self.isrw else 'gen' ,self.m,self.vv, self.ctau)
    if isrw:
      evt_w = '(weight_{vv})'.format(vv=str(self.vv).replace('-', 'm'))
    else:
      evt_w = '(1)'


   
    self.histoDefs = {
    # b mother
    'b_pt'          : PlotOpt('b_pt', '(30,0,30)', 'B meson p_{T} [GeV]', 'a.u.', False, True),
    'b_eta'         : PlotOpt('b_eta', '(30,-6,6)', 'B meson #eta', 'a.u.', False, True),
    'b_ct'          : PlotOpt('b_ct_reco', '(50,0,1000)', 'B meson ct [mm]', False, True),
    'b_ct_large'    : PlotOpt('b_ct_reco', '(100,0,10000)', 'B meson ct [mm]', False, True),     
    # daughters of the B
    ## the HNL
    'hnl_pt'        : PlotOpt('hnl_pt', '(30,0,30)', 'HNL p_{T} [GeV]', 'a.u.', False, True),   
    'hnl_eta'       : PlotOpt('hnl_eta', '(30,-6,6)', 'HNL #eta', 'a.u.', False, True),      
    'hnl_ct'        : PlotOpt('hnl_ct_reco', '(50,0,1000)', 'HNL ct [mm]', False, True),
    'hnl_ct_large'  : PlotOpt('hnl_ct_reco' '(100,0,10000)', 'HNL ct [mm]', False, True),    
    'hnl_Lxy'       : PlotOpt('Lxy', '(50,0,1000)', 'L_{xy} [mm]', False, True),
    'hnl_Lxy_large' : PlotOpt('Lxy' '(100,0,10000)', 'L_{xy} [mm]', False, True),    
    'hnl_Lxyz'       : PlotOpt('Lxyz', '(50,0,1000)', 'L_{xyz} [mm]', False, True),
    'hnl_Lxyz_large' : PlotOpt('Lxyz' '(100,0,10000)', 'L_{xyz} [mm]', False, True),    

    ### the D meson
    'd_pt'          : PlotOpt('d_pt', '(30,0,30)', 'D meson p_{T} [GeV]', 'a.u.', False, True),   
    'd_eta'         : PlotOpt('d_eta', '(30,-6,6)', 'D meson #eta', 'a.u.', False, True),      

    ### the trigger lepton
    'mutrig_pt'     : PlotOpt('l0_pt', '(30,0,30)', '#mu^{trig} p_{T} [GeV]', 'a.u.', False, True),   
    'mutrig_eta'    : PlotOpt('l0_eta', '(30,-6,6)', '#mu^{trig} #eta', 'a.u.', False, True),      

    #### daughters of the D meson
    ###### the pion
    'piD_pt'        : PlotOpt('pi_pt', '(30,0,30)', '#pi (from D) p_{T} [GeV]', 'a.u.', False, True),  
    'piD_eta'       : PlotOpt('pi_eta', '(30,-6,6)', '#pi (from D) #eta', 'a.u.', False, True),      
    
    ###### the kaon
    'k_pt'          : PlotOpt('k_pt', '(30,0,30)', 'K (from D) p_{T} [GeV]', 'a.u.', False, True),  
    'k_eta'         : PlotOpt('k_eta', '(30,-6,6)', 'K (from D) #eta', 'a.u.', False, True),      
    
    #### daughters of the HNL
    ###### the lepton
    'mu_pt'         : PlotOpt('hnl_pt', '(30,0,30)', '#mu (from HNL) p_{T} [GeV]', 'a.u.', False, True),  
    'mu_eta'        : PlotOpt('hnl_eta', '(30,-6,6)', '#mu (from HNL) #eta', 'a.u.', False, True),      

    ###### the pion
    'pi_pt'         : PlotOpt('pi_pt', '(30,0,30)', '#pi (from HNL) p_{T} [GeV]', 'a.u.', False, True),   
    'pi_eta'        : PlotOpt('pi_eta', '(30,-6,6)', '#pi (from HNL) #eta', 'a.u.', False, True),      
   
    # invariant masses
    'hnl_invmass'   : PlotOpt('lep_pi_invmass', '(20,0,5)', 'HNL invariant mass, m(#mu,#pi) [GeV]', 'a.u.', False, False),     
    'd_invmass'     : PlotOpt('k_pi_invmass', '(20,0,5)', 'D meson invariant mass, m(K,#pi) [GeV]', 'a.u.', False, False),      
    'b_invmass'     : PlotOpt('hn_d_pl_invmass', '(20,0,5)', 'B meson invariant mass, m(HNL,D,#mu^{trig}) [GeV]', 'a.u.', False, False),     
    #'Lxy_cos', # cosine of the pointing angle in the transverse plane
    #'Lxyz_b', #3D displacement of the B wrt to primary vertex
    #'Lxyz_l0' #3D displacement of the prompt lepton wrt to B vertex
    }
    self.histos = {}
  
  def fillMetaData(self):
    '''
    This should contain code to retrieve mainly the filter efficiency, perhaps the cross-section from the log file
    '''
    pass
    
  def stamp(self):
    '''
    This should print the basic information of the sample
    '''
    pass
  
  def fillHistos(self):
    '''
    This is to fill the histograms
    '''
    f = TFile.Open(self.infileName)
    t = f.Get(self.treeName)
    if not t:
      print 'ERROR: no tree in file %s' % fname
      return False

    for name,spec in self.histoDefs.items():
      t.Draw('%s>>%s%s' % (spec.what, name, spec.binning), self.evt_w, 'goff') 
      h = t.GetHistogram().Clone('%s_%s' % (self.name, name))
      h.SetDirectory(0)
      h.SetTitle(h.GetName())
      h.GetXaxis().SetTitle(spec.xtitle)
      h.GetXaxis().SetNdivisions(505)
      h.GetYaxis().SetTitle(spec.ytitle)
      h.GetYaxis().SetNdivisions(505)
      #h.GetYaxis().SetTitleSize(0.05)
      #h.GetXaxis().SetTitleSize(0.05)

      newh = getOverflowedHisto(h)
      newh.SetDirectory(0)
      self.histos[name] = newh

    f.IsA().Destructor(f)
    return True

  def saveHistos(self, norm):
    '''
    This is to save the histograms of the sample in a plot
    '''
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

    c.SaveAs('./plots/all_plots/' + c.GetName() + norm_suffix + '.png')
    c.SaveAs('./plots/all_plots/' + c.GetName() + norm_suffix + '.pdf')

  def fillAcceptance(self):
    '''
    This is to calculate the acceptance of the sample
    '''
    f = TFile.Open(self.infileName)
    t = f.Get(self.treeName)

    self.effnum = ROOT.TH1F('effnum', 'effnum', 1, 0, 13000) #dict([(k, ROOT.TH1F('effnum_%s' % k, 'effnum_%s' % k, 1, 0, 1)) for k in self.settings])
    self.effden = ROOT.TH1F('effden', 'effden', 1, 0, 13000)

    cutsnum = '(l0_pt>7 && abs(l0_eta)<1.5 && l1_pt>3 && abs(l1_eta)<2.5 && pi_pt>1 && abs(pi_eta)<2.5 && k_pt>1 && abs(k_eta)<2.5 && pi1_pt>1 && abs(pi1_eta)<2.5 && Lxy < 1000)'
    cutsden = '(l0_pt>7 && abs(l0_eta)<1.5)'
    

    tree.Draw('hnl_pt>>effnum', cutsnum+'*'+self.evt_w, 'goff')
    tree.Draw('hnl_pt>>effden', cutsden+'*'+self.evt_w, 'goff')
 
    if TEfficiency.CheckConsistency(effnum,effden): 
      peff = TEfficiency(effnum,effden)
      
      self.acc = peff.GetEfficiency(1)
      self.acc_errup = peff.GetEfficiencyErrorUp(1)
      self.acc_errdn = peff.GetEfficiencyErrorLow(1)

    else:
      self.acc = -99
      self.acc_errup = -99
      self.acc_errdn = -99

class SampleList(object):
  '''
  Handles the plotting of several samples
  '''
  def __init__(self, name, samples):
    self.name = name
    self.samples = samples  # a list of objects of type Sample
    self.colors = [   ROOT.kOrange+1, ROOT.kRed, ROOT.kMagenta+2, ROOT.kViolet+8, ROOT.kAzure-8, ROOT.kAzure+6 ,
                      ROOT.kGreen+1, ROOT.kSpring+4, ROOT.kYellow -5, ROOT.kYellow -3, ROOT.kYellow, ROOT.kOrange
                  ]
    self.styles = [  1,1,1,1,1,1,1,1,1,1,1,1,1,1, ]
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
    Makes sure samples have same binning and histograms
    '''
    for sample1 in self.samples:
      for sample2 in self.samples:
        if len(sample1.histoDefs) != len(sample2.histoDefs):
          raise RuntimeError('number of PlotOpts differs for %s and %s' % (sample1.name, sample2.name))
        if len(sample1.histos) != len(sample2.histos):
          raise RuntimeError('number of histos differs for %s and %s' % (sample1.name, sample2.name))
  def saveAll(self, norm=False, sameCanvas=False):
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
    graph_saver = []

    hMax = {}
    for j, what in enumerate(self.samples[0].histos.keys()):
      hMax[what] = []
      leg=defaultLegend(x1=0.3,y1=0.7,x2=0.95,y2=0.95,mult=1.2)
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
            tosave[cname] = TCanvas(cname, cname, 800, 600)
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
          hMax[what][0].SetMaximum(max*7)
          #print hMax[what][0].GetMaximum()
        else:
          #hMax[name][0].GetYaxis().SetRangeUser(0.01, max * 1.5)
          hMax[what][0].SetMaximum(max * 1.5)

    norm_suffix='_norm' if norm else ''

    for cname, canv in tosave.items():
      canv.SaveAs('./plots/all_plots/' + cname + norm_suffix +'.png')
      canv.SaveAs('./plots/all_plots/' + cname + norm_suffix +'.pdf')
      canv.SaveAs('./plots/all_plots/' + cname + norm_suffix +'.C')


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='', add_help=True)
  parser.add_argument('--pl', type=str, dest='pl', help='production label of input', default='V_blabla')  
  return parser.parse_args()


def getLxyPlot(s):
  f = TFile.Open(s.infileName)
  if f:
    tree = f.Get(s.treeName)
    if not tree: raise RuntimeError('no tree found')
  else: raise RuntimeError('file %s not found' % infileName)  

  c=TCanvas() 

  holdtemp = TH1F('holdtemp', 'holdtemp', 100, 0., 10000.) 
  tree.Draw('Lxy>>holdtemp', '', 'goff')
  hold = getOverflowedHisto(holdtemp)
  hold.SetLineColor(ROOT.kBlack)
  hold.SetMarkerColor(ROOT.kBlack)
  hold.SetMarkerStyle(1)
  hold.GetXaxis().SetTitle('L_{xy} [mm]') 
  hold.GetYaxis().SetTitle('a.u.') 
  hold.DrawNormalized('histE')

  leg = TLegend() 
  leg=defaultLegend(x1=0.3,y1=0.7,x2=0.95,y2=0.90,mult=1.2)
  leg.AddEntry(hold, 'gen: |V|^{{2}}={:.1e} ctau={:.1f}mm'.format(s.vv, s.ctau), 'PEL')   # also here, property of the sample needs to be recognised...

  hnewtemps=[]
  hnews=[]

  for i,vv in enumerate(small_new_vvs):
    # only do 

    hnewtemps.append( TH1F('hnewtemp_{}'.format(i), 'hnewtemp_{}'.format(i), 100,0.,10000.))
    weight = '(weight_{vv})'.format(vv=str(vv).replace('-', 'm'))
    tree.Draw('Lxy>>hnewtemp_{}'.format(i), weight, 'goff')

    hnews.append(getOverflowedHisto(hnewtemps[i]))
        
    hnews[i].SetLineColor(colors[i])
    hnews[i].SetMarkerColor(colors[i])
    hnews[i].SetMarkerStyle(1)
    
    hnews[i].DrawNormalized('histEsame')

    leg.AddEntry(hnews[i], 'new: |V|^{{2}}={:.1e} ctau={:.1f}mm'.format(vv, getCtau(s.mass,vv)), 'PEL')   

  leg.Draw('same')
  c.SaveAs('plots/lxy_reweight.pdf')



def getAcceptance(s):

  f = TFile.Open(s.infileName)
  if f:
    tree = f.Get(s.treeName)
    if not tree: raise RuntimeError('no tree found')
  else: raise RuntimeError('file %s not found' % s.infileName)  

  cutsnum = '(l0_pt>7 && abs(l0_eta)<1.5 && l1_pt>3 && abs(l1_eta)<2.5 && pi_pt>1 && abs(pi_eta)<2.5 && k_pt>1 && abs(k_eta)<2.5 && pi1_pt>1 && abs(pi1_eta)<2.5 && Lxy < 1000)'
  cutsden = '(l0_pt>7 && abs(l0_eta)<1.5)'

  graph = TGraph()
  graph_acc = TGraph()

  for i,vv in enumerate(new_vvs): 
    weight='*(weight_{vv})'.format(vv=str(vv).replace('-', 'm'))
    effnum = TH1F('effnum', 'effnum', 1, 0, 13000)
    effden = TH1F('effden', 'effden', 1, 0, 13000)

    tree.Draw('hnl_pt>>effnum', cutsnum+weight, 'goff')
    tree.Draw('hnl_pt>>effden', cutsden+weight, 'goff')

 
    if TEfficiency.CheckConsistency(effnum,effden): 
      peff = TEfficiency(effnum,effden)
      
      N_nu = 3.2E7
      BR = 19.7 / 100.
      eff = peff.GetEfficiency(1)
      eff_errup = peff.GetEfficiencyErrorUp(1)
      eff_errdn = peff.GetEfficiencyErrorLow(1)
      ## xsec_rescale = vv / 0.0013

      print 'VV={:.2e} origVV={:.2e}'.format(vv,s.vv)
      print 'eff={:.2f}+{:.2f}-{:.2f}'.format(eff,eff_errup,eff_errdn)
      print 'N_nu * BR * eff x VV = {:.2e} '.format(N_nu*BR*eff*vv)

      graph.SetPoint(i,vv,N_nu*BR*eff*vv)
      graph_acc.SetPoint(i,vv,eff)
    
  c = TCanvas()
  graph.SetLineWidth(2)
  graph.SetMarkerStyle(22)
  graph.SetTitle(';|V|^{2};N_{#nu} x BR(HN#rightarrow#mu#pi) x Acc x V^{2}')
  graph.Draw('APL')

  gPad.Modified()
  gPad.Update()
  line = TLine(gPad.GetUxmin(),3,gPad.GetUxmax(),3)
  line.SetLineColor(ROOT.kBlue)
  line.Draw('same')

  c.SetLogx()
  c.SetLogy()
  c.SetGridx()
  c.SetGridy()
  c.SaveAs('./plots/expected_nevts.pdf')
  c.SaveAs('./plots/expected_nevts.C')
  #graph.SetLineColor(k)

  c_acc = TCanvas()
  c_acc.SetLogx()
  c_acc.SetLogy()
  
  
  graph_acc.SetLineWidth(2)
  graph_acc.SetLineColor(kRed)
  graph_acc.SetTitle(';|V|^{2};Acceptance')
  graph_acc.Draw('APL')
  graph_acc.SetMaximum(0.1)
  #graph_acc.SetMinimum(0.025)
  c_acc.SaveAs('./plots/expected_acc.pdf')
  c_acc.SaveAs('./plots/expected_acc.C')


if __name__ == "__main__":

  gROOT.SetBatch(True)
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  gStyle.SetTitleXOffset(1.1);
  gStyle.SetTitleYOffset(1.45);


  opt = getOptions()

  infileName = './outputfiles/{}_miniGenTree.root'.format(opt.pl) 
  
  tags = opt.pl.split('_')
  this_mass = float([it.split('mass')[1] for it in tags if 'mass' in it][0])
  this_ctau = float([it.split('ctau')[1] for it in tags if 'ctau' in it][0])
  #this_vv = getVV(mass=this_mass, ctau=this_ctau)

  this_sample = Sample(mass=this_mass, ctau=this_ctau, treeName='tree', infileName=infileName)

  getAcceptance(s=this_sample)
  #getAcceptance(inFileName='./outputfiles/{}_miniGenTree.root'.format(opt.pl), treeName='tree')
  getLxyPlot(s=this_sample)
  #getLxyPlot(inFileName='./outputfiles/{}_miniGenTree.root'.format(opt.pl), treeName='tree')


