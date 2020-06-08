import os
from ROOT import TTree, TFile, TH1F, TEfficiency


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

  cutsnum = 'l0_pt>7 && abs(l0_eta)<1.5 && l1_pt>5 && abs(l1_eta)<2.5 && pi_pt>5 && abs(pi_eta)<2.5 && k_pt>5 && abs(k_eta)<2.5 && pi1_pt>5 && abs(pi1_eta)<2.5'
  cutsden = 'l0_pt>7 && abs(l0_eta)<1.5'

  effnum = TH1F('effnum', 'effnum', 1, 0, 13000)
  effden = TH1F('effden', 'effden', 1, 0, 13000)

  tree.Draw('b_pt>>effnum', cutsnum, 'goff')
  tree.Draw('b_pt>>effden', cutsden, 'goff')
  
  if TEfficiency.CheckConsistency(effnum,effden): 
    eff = TEfficiency(effnum,effden)

    print 'eff = {:3f} + {:3f} - {:3f}'.format(eff.GetEfficiency(1), eff.GetEfficiencyErrorUp(1), eff.GetEfficiencyErrorLow(1))
 

if __name__ == "__main__":

  opt = getOptions()

  getAcceptance(inFileName='./{}_miniGenTree.root'.format(opt.pl), treeName='tree')
