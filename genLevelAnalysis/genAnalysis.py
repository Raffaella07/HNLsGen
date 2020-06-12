import os
from ROOT import TTree, TFile, TH1F, TEfficiency, TGraph, TCanvas, gROOT, TAxis
from ROOT import kRed
# couplings to be tested, for which the reweight is run
# from Riccardo
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


def getOptions():

  from argparse import ArgumentParser

  parser = ArgumentParser(description='', add_help=True)

  parser.add_argument('--pl', type=str, dest='pl', help='production label of input', default='V_blabla')  

  return parser.parse_args()


def getAcceptance(inFileName, treeName):

  f = TFile.Open(inFileName)
  if f:
    tree = f.Get(treeName)
    if not tree: raise RuntimeError('no tree found')
  else: raise RuntimeError('file %s not found' % inFileName)  

  cutsnum = '(l0_pt>7 && abs(l0_eta)<1.5 && l1_pt>3 && abs(l1_eta)<2.5 && pi_pt>1 && abs(pi_eta)<2.5 && k_pt>1 && abs(k_eta)<2.5 && pi1_pt>1 && abs(pi1_eta)<2.5 && Lxy < 1000)'
  cutsden = '(l0_pt>7 && abs(l0_eta)<1.5)'


  graph = TGraph()
  graph_acc = TGraph()

  for i,iv2 in enumerate(new_v2s): 
    weight='*(weight_{vv})'.format(vv=str(iv2).replace('-', 'm'))
    effnum = TH1F('effnum', 'effnum', 1, 0, 13000)
    effden = TH1F('effden', 'effden', 1, 0, 13000)

    tree.Draw('b_pt>>effnum', cutsnum+weight, 'goff')
    tree.Draw('b_pt>>effden', cutsden+weight, 'goff')
    #tree.Draw('b_pt>>effnum', cutsnum, 'goff')
    #tree.Draw('b_pt>>effden', cutsden, 'goff')

 
    if TEfficiency.CheckConsistency(effnum,effden): 
      peff = TEfficiency(effnum,effden)
      
      N_nu = 3.6E7
      BR = 19.7 / 100.
      eff = peff.GetEfficiency(1)
      eff_errup = peff.GetEfficiencyErrorUp(1)
      eff_errdn = peff.GetEfficiencyErrorLow(1)
      xsec_rescale = iv2 / 0.0013

      print 'Old_v2 = 1.3E-3, new_v2={:2e}'.format(iv2)
      print 'eff = {:2f} + {:2f} - {:2f}'.format(eff,eff_errup,eff_errdn)
      print 'N_nu * BR * xsec_new/xsec_old * eff = {:2e} '.format(N_nu*BR*eff*xsec_rescale)

      graph.SetPoint(i,iv2,N_nu*BR*eff*xsec_rescale)
      graph_acc.SetPoint(i,iv2,eff)
    
  c = TCanvas()
  graph.SetLineWidth(2)
  graph.SetMarkerStyle(22)
  #graph.SetTitle(';target |V|^{2};N_{#nu}(B^{\pm}->D0(->K#pi)#mu#nu_{#mu} * BR(HN->#mu#pi) * eff');
  graph.SetTitle(';|V|^{2}_{target};N_{#nu} x BR(HN#rightarrow#mu#pi) x |V|^{2}_{target}/|V|^{2}_{orig}) x Acc')
  #graph.GetXaxis().SetTitle('target |V|^2')
  #graph.GetYaxis().SetTitle('N_nu(B^{\pm}->D0(->K#pi)#mu#nu_{#mu} * BR(HN->#mu#pi) * eff')
  #graph.SetMaximum(800000)  
  #graph.SetMinimum(100000)  
  graph.Draw('APL')
  
  c.SetLogx()
  c.SetLogy()
  c.SetGridx()
  c.SetGridy()
  c.SaveAs('expected_nevts_reweight.pdf')
  c.SaveAs('expected_nevts_reweight.C')
  #graph.SetLineColor(k)

  c_acc = TCanvas()
  c_acc.SetLogx()
  
  graph_acc.SetLineWidth(2)
  graph_acc.SetLineColor(kRed)
  graph_acc.SetTitle(';|V|^{2}_{target};Acceptance')
  graph_acc.Draw('APL')
  graph_acc.SetMaximum(0.11)
  graph_acc.SetMinimum(0.10)
  c_acc.SaveAs('expected_acc_reweight.pdf')
  c_acc.SaveAs('expected_acc_reweight.C')


if __name__ == "__main__":

  gROOT.SetBatch(True)
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')

  opt = getOptions()

  getAcceptance(inFileName='./{}_miniGenTree.root'.format(opt.pl), treeName='tree')
