import os
import sys
import numpy as np
from ROOT import TTree, TFile, TH1F, TEfficiency, TGraph, TCanvas, gROOT, TAxis, TMath, TLegend, gStyle, gPad, TLine, TGraphAsymmErrors
import ROOT
sys.path.append('/work/mratti/plotting/myplotting')
from spares import *
from glob import glob

# couplings to be tested, for which the reweight is run

from genTreeProducer import new_vvs
#small_new_vvs = new_vvs[::4]  # only select one every three elements
#small_new_vvs.reverse()
new_vvs.reverse()
small_new_vvs = new_vvs[0:10]

sys.path.append('../python/') 
from common import getVV,getCtau,Point

def getOverflowedHisto(h):
  htemp = h.Clone()
  nbins = htemp.GetNbinsX()
  last_plus_overflow = htemp.GetBinContent(nbins) + htemp.GetBinContent(nbins+1)
  last_plus_overflow_error = TMath.Sqrt( htemp.GetBinError(nbins)*htemp.GetBinError(nbins)  + htemp.GetBinError(nbins+1)*htemp.GetBinError(nbins+1))
  htemp.SetBinContent(nbins,last_plus_overflow )
  htemp.SetBinError(nbins,last_plus_overflow_error)
  return htemp


class SampleQuantity(object):
  def __init__(self, name='vv', axis='x', title='|V|^{2}', log=True, range=(0.,1), forceRange=False, err=False):
    self.name = name
    self.axis = axis
    self.title = title
    self.log = log
    self.range = range
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
  def __init__(self,mass=-99, ctau=-99, vv=-99, infileName=None, isrw=False, orig_vv=None, label='V00'):
    self.mass = mass
    self.ctau = ctau
    self.vv = vv
    if vv==-99: 
      self.vv=getVV(mass=self.mass, ctau=self.ctau)
    if ctau==-99:
      self.ctau=getCtau(mass=self.mass, vv=self.vv)
    self.treeName='tree'
    self.infileName=infileName
    if not os.path.isfile(self.infileName): raise RuntimeError('file %s not found' % self.infileName)
    self.isrw=isrw
    if self.isrw:
      self.orig_vv = orig_vv
      self.orig_ctau = getCtau(mass=self.mass,vv=orig_vv)
    else:
      self.orig_vv = self.vv
      self.orig_ctau = self.ctau
    self.label=label
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
  
   
    self.histoDefs = {
    # b mother
    'b_pt'          : PlotOpt('b_pt', '(30,0,30)', 'B meson p_{T} [GeV]', 'a.u.', False, True),
    'b_eta'         : PlotOpt('b_eta', '(30,-6,6)', 'B meson #eta', 'a.u.', False, True),
    'b_ct'          : PlotOpt('b_ct_reco', '(50,0,100)', 'B meson ct [mm]', 'a.u.', False, True),
    #'b_ct_large'    : PlotOpt('b_ct_reco', '(100,0,1000)', 'B meson ct [mm]', 'a.u.', False, True),     
    # daughters of the B
    ## the HNL
    'hnl_pt'        : PlotOpt('hnl_pt', '(30,0,30)', 'HNL p_{T} [GeV]', 'a.u.', False, True),   
    'hnl_eta'       : PlotOpt('hnl_eta', '(30,-6,6)', 'HNL #eta', 'a.u.', False, True),      
    'hnl_ct'        : PlotOpt('hnl_ct_reco', '(50,0,1000)', 'HNL ct [mm]', 'a.u.', False, True),
    'hnl_ct_large'  : PlotOpt('hnl_ct_reco', '(100,0,10000)', 'HNL ct [mm]', 'a.u.', False, True),    
    'hnl_Lxy'       : PlotOpt('Lxy', '(50,0,1000)', 'L_{xy} [mm]', 'a.u.', False, True),
    'hnl_Lxy_large' : PlotOpt('Lxy', '(100,0,10000)', 'L_{xy} [mm]', 'a.u.', False, True),    
    'hnl_Lxyz'       : PlotOpt('Lxyz', '(50,0,1000)', 'L_{xyz} [mm]', 'a.u.', False, True),
    'hnl_Lxyz_large' : PlotOpt('Lxyz', '(100,0,10000)', 'L_{xyz} [mm]', 'a.u.', False, True),    

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
    'hnl_invmass'   : PlotOpt('lep_pi_invmass', '(50,0,5)', 'HNL invariant mass, m(#mu,#pi) [GeV]', 'a.u.', False, False),     
    'd_invmass'     : PlotOpt('k_pi_invmass', '(50,0,5)', 'D meson invariant mass, m(K,#pi) [GeV]', 'a.u.', False, False),      
    'b_invmass'     : PlotOpt('hn_d_pl_invmass', '(50,2,7)', 'B meson invariant mass, m(HNL,D,#mu^{trig}) [GeV]', 'a.u.', False, False),     
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
    print('mass={m}GeV, ctau={ctau}mm VV={vv}, isrw={isrw}, orig_VV={ovv}, acc={acc}, evt_w={ew}'.format( \
            m=self.mass,ctau=self.ctau,vv=self.vv,isrw=self.isrw,ovv=self.orig_vv,acc=self.acc,ew=self.evt_w))
  
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

      #newh = getOverflowedHisto(h)
      #newh.SetDirectory(0)
      #self.histos[name] = newh
      h.SetDirectory(0)
      self.histos[name] = h

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

    c.SaveAs('./plots/' + self.label + '/' + c.GetName() + norm_suffix + '.png')
    c.SaveAs('./plots/' + self.label + '/' + c.GetName() + norm_suffix + '.pdf')

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
    

    t.Draw('hnl_pt>>effnum', cutsnum+'*'+self.evt_w, 'goff')
    t.Draw('hnl_pt>>effden', cutsden+'*'+self.evt_w, 'goff')
 
    if TEfficiency.CheckConsistency(self.effnum,self.effden): 
      peff = TEfficiency(self.effnum,self.effden)
      
      self.acc = peff.GetEfficiency(1)
      self.acc_errup = peff.GetEfficiencyErrorUp(1)
      self.acc_errdn = peff.GetEfficiencyErrorLow(1)

  def fillExpNevts(self):
    N_nu = 3.2E7
    BR = 19.7 / 100.
    self.expNevts = self.acc * N_nu * BR * self.vv


class SampleList(object):
  '''
  Handles the plotting of several samples, both the histograms, and the general quantities (graphs)
  '''
  def __init__(self, name, samples, label=None):
    self.name = name
    self.samples = samples  # a list of objects of type Sample
    self.label = samples[0].label if label==None else label
    self.colors = [   ROOT.kOrange+1, ROOT.kRed, ROOT.kMagenta+2, ROOT.kViolet+8, ROOT.kAzure-8, ROOT.kAzure+6 ,
                      ROOT.kGreen+1, ROOT.kSpring+4, ROOT.kYellow -5, ROOT.kYellow -3, ROOT.kYellow, ROOT.kOrange
                  ]
    self.styles = [  1,1,1,1,1,1,1,1,1,1,1,1,1,1, ]

    # x,y quantities for graphs
    self.quantities={
      #'acc' : SampleQuantity(name='acc', axis='y', title='Acceptance', log=True, range=(0.,1), forceRange=True, err=True),
      'acc' : SampleQuantity(name='acc', axis='y', title='Acceptance', log=False, range=(0.,1), forceRange=False, err=True),
      'expNevts': SampleQuantity(name='expNevts', axis='y', title='N_{#nu} x BR(HN#rightarrow#mu#pi) x Acc x V^{2}', log=True, err=False),
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
    graph_saver = []

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
      canv.SaveAs('./plots/' + self.label + '/' + cname + norm_suffix + '.png')
      canv.SaveAs('./plots/' + self.label + '/' + cname + norm_suffix + '.pdf')

  def plotGraph(self, x='vv', y='acc'): 
    '''
    Plot a graph with specified quantities on x and y axes 
    '''

    if (x not in self.quantities.keys()) or (y not in self.quantities.keys()):
      raise RuntimeError('selected quantities not available, available quantities are: \n{}'.format(self.quantities.keys()))

    xq = self.quantities[x]
    yq = self.quantities[y]

    graph = TGraphAsymmErrors()
    for i,s in enumerate(self.samples):
      graph.SetPoint(i,getattr(s, xq.name), getattr(s, yq.name) )
      if xq.err: 
        graph.SetPointEXhigh(i, getattr(s, xq.name+'_errup'))   # errup errdn
        graph.SetPointEXlow (i, getattr(s, xq.name+'_errdn'))
      if yq.err: 
        graph.SetPointEYhigh(i, getattr(s, yq.name+'_errup'))  
        graph.SetPointEYlow (i, getattr(s, yq.name+'_errdn'))

    c = TCanvas()
    graph.SetLineWidth(2)
    graph.SetMarkerStyle(22)
    graph.SetTitle(';{x};{y}'.format(y=yq.title,x=xq.title))
    graph.Draw('APL')

    gPad.Modified()
    gPad.Update()
    if yq.name=='expNevts':
      line = TLine(gPad.GetUxmin(),3,gPad.GetUxmax(),3)
      line.SetLineColor(ROOT.kBlue)
      line.Draw('same')

    if xq.log: c.SetLogx()
    if yq.log: c.SetLogy()
    c.SetGridx()
    c.SetGridy()
    c.SaveAs('./plots/{}/{}_{}VS{}.pdf'.format(self.label,self.name,yq.name,xq.name))
    c.SaveAs('./plots/{}/{}_{}VS{}.png'.format(self.label,self.name,yq.name,xq.name))


def doAnalysis(path,points,name):
  '''
  Perform plotting of samples lists
  '''
  label = path.split('/')[2]
  print('  => Going to do closure of reweighting analysis for production={}, name={}'.format(label,name))
  samples = []
  for p in points:
    fn = path.format(m=p.mass,ctau=p.orig_ctau)
    s = Sample(mass=p.mass, ctau=p.ctau, infileName=fn, isrw=p.isrw, orig_vv=p.orig_vv, label=label)
    s.fillHistos()
    s.fillAcceptance()
    s.fillExpNevts()
    s.stamp()
    samples.append(s)

  slist = SampleList(name=name, samples=samples, label=label)
  slist.plotHistos(norm=True, sameCanvas=False)
  if 'fixedM' in name or 'closure' in name: 
    slist.plotGraph(x='vv',y='acc')
    slist.plotGraph(x='vv',y='expNevts')
  elif 'fixedVV' in name:
    slist.plotGraph(x='mass',y='acc')
    slist.plotGraph(x='mass',y='expNevts')

def checkFiles(path,points):
  '''
  Check existence of files for analysis, return existing points
  '''
  for p in points:
    fn = path.format(m=p.mass,ctau=p.orig_ctau) # check the file where the tree is stored
    if not os.path.isfile(fn):
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
  return parser.parse_args()

if __name__ == "__main__":

  gROOT.SetBatch(True)
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  gStyle.SetTitleXOffset(1.1);
  gStyle.SetTitleYOffset(1.45);

  global debug
  debug = False

  opt = getOptions()

  os.system('mkdir -p ./plots/{}'.format(opt.pl))

  path = './outputfiles/' + opt.pl + '/mass{m}_ctau{ctau}_miniGenTree.root'

  ################
  points_for_fixedMass_analysis = [
    Point(mass=1.5,ctau=None,vv=5e-03,isrw=False),
    Point(mass=1.5,ctau=None,vv=1e-03,isrw=False),
    Point(mass=1.5,ctau=None,vv=5e-04,isrw=False),
    Point(mass=1.5,ctau=None,vv=1e-04,isrw=False),
    Point(mass=1.5,ctau=None,vv=5e-05,isrw=False),
  ]
  for p in points_for_fixedMass_analysis:
   p.stamp()
  existing_points=checkFiles(path=path,points=points_for_fixedMass_analysis)
  doAnalysis(path=path,points=existing_points,name='fixedMass1.5')
  ################
  points_for_fixedMass_analysis = [
    Point(mass=0.5,ctau=None,vv=5e-03,isrw=False),
    Point(mass=0.5,ctau=None,vv=1e-03,isrw=False),
    Point(mass=0.5,ctau=None,vv=5e-04,isrw=False),
    Point(mass=0.5,ctau=None,vv=1e-04,isrw=False),
    Point(mass=0.5,ctau=None,vv=5e-05,isrw=False),
  ]
  for p in points_for_fixedMass_analysis:
   p.stamp()
  existing_points=checkFiles(path=path,points=points_for_fixedMass_analysis)
  doAnalysis(path=path,points=existing_points,name='fixedMass0.5')
  ################
  points_for_fixedMass_analysis = [
    Point(mass=1.0,ctau=None,vv=5e-03,isrw=False),
    Point(mass=1.0,ctau=None,vv=1e-03,isrw=False),
    Point(mass=1.0,ctau=None,vv=5e-04,isrw=False),
    Point(mass=1.0,ctau=None,vv=1e-04,isrw=False),
    Point(mass=1.0,ctau=None,vv=5e-05,isrw=False),
  ]
  for p in points_for_fixedMass_analysis:
   p.stamp()
  existing_points=checkFiles(path=path,points=points_for_fixedMass_analysis)
  doAnalysis(path=path,points=existing_points,name='fixedMass1.0')

  ###############
  points = [
    Point(mass=1.0,ctau=None,vv=1e-03, isrw=False),
    Point(mass=1.0,ctau=None,vv=1e-03, isrw=True, orig_vv=5e-03),
    Point(mass=1.0,ctau=None,vv=5e-03, isrw=False),
  ]
  existing_points = checkFiles(path,points)
  doAnalysis(path=path,points=existing_points,name='closureRw_Mass1.0')

  ################
  points_for_fixedMass_analysis = [
    Point(mass=0.5,ctau=None,vv=5e-04,isrw=False),
    Point(mass=1.0,ctau=None,vv=5e-04,isrw=False),
    Point(mass=1.5,ctau=None,vv=5e-04,isrw=False),
  ]
  for p in points_for_fixedMass_analysis:
   p.stamp()
  existing_points=checkFiles(path=path,points=points_for_fixedMass_analysis)
  doAnalysis(path=path,points=existing_points,name='fixedVV5em04')

