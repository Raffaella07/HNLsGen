'''
Production helper for generating private samples for B-initiated HNLs
=> for the moment up to gen
=> would like to add a script that checks the status of the production
'''

import sys
import os
import subprocess



def makeTemplate():

  template = [
    '#!/bin/bash',
    '',
    '#SBATCH -J step1_m{m}_ctau{ctau}',
    '#SBATCH -o logs/step1_%a.log', 
    '#SBATCH -e logs/step1_%a.log',
    #'#SBATCH -o logs/step1_%A_%a.log', 
    #'#SBATCH -e logs/step1_%A_%a.log',
    '#SBATCH -p wn',
    '#SBATCH -t {hh}:00:00',
    '#SBATCH --array={arr}',
    '#SBATCH --ntasks=1',
    '#SBATCH --account=t3',
    '',
    'DIRNAME="{pl}"',
    'STARTDIR=/work/mratti/GEN_HNL_newPythia/CMSSW_10_2_3/src/HNLsGen/',
    'TOPWORKDIR="/scratch/{user}/"',
    'JOBDIR="gen_$SLURM_JOB_ID_$SLURM_ARRAY_TASK_ID"',
    'WORKDIR=$TOPWORKDIR/$JOBDIR',
    'INSEPREFIX="root://t3dcachedb.psi.ch:1094/"',
    'OUTSEPREFIX="root://t3dcachedb.psi.ch:1094/"',
    'SERESULTDIR="/pnfs/psi.ch/cms/trivcat/store/user/{user}/BHNLsGen/"$DIRNAME'
    '',
    'source $VO_CMS_SW_DIR/cmsset_default.sh',
    'shopt -s expand_aliases',
    'echo ""',
    'echo "Going to set up cms environment"',
    'cd $STARTDIR',
    'cmsenv',
    'echo ""',
    '',
    'echo "Going to create work dir"',
    'mkdir -p $WORKDIR',
    'echo "workdir: "',
    'echo $WORKDIR',
    'echo ""',
    '',
    'echo "Going to create the output dir"',
    'echo "May give an error if the directory already exists, which can be ignored"', # once python bindings to interact with SE are available, should be easier...
    'xrdfs $T3CACHE mkdir $SERESULTDIR',
    'echo ""',
    '',
    'echo "Going to copy cmsdriver to work dir"',
    'cp $STARTDIR/cmsDrivers/{jop} $WORKDIR/. ',
    'echo ""',
    '',
    'cd $WORKDIR',
    'pwd',
    'echo "Going to run"',
    'DATE_START=`date +%s`',
    'cmsRun {jop} maxEvents={nevtsjob} nThr={nthr} outputFile=BPH-test.root seedOffset=$SLURM_ARRAY_TASK_ID',
    'DATE_END=`date +%s`',
    'echo "Finished running"',
    'echo "Content of current directory"',
    'ls -al',
    'echo ""',
    '',
    'echo "Going to copy output to result directory"',
    'xrdcp -f $WORKDIR/BPH-test_numEvent{nevtsjob}.root $OUTSEPREFIX/$SERESULTDIR/step1_nj$SLURM_ARRAY_TASK_ID".root"',
    '',
    'echo ""',
    'echo "Cleaning up $WORKDIR"',
    'rm -rf $WORKDIR',
    'RUNTIME=$((DATE_END-DATE_START))',
    'echo "Wallclock running time: $RUNTIME s"',
    'cd $STARTDIR',
  ]

  template = '\n'.join(template)
  return template


def getOptions():

  # convention: no capital letters

  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for B-initiated HNL signals', add_help=True)

  parser.add_argument('-v','--ver', type=str, dest='ver', help='version of production, e.g. V00_v00', default='V00_v00')
  parser.add_argument('-n','--nevts', type=int, dest='nevts', help='total number of events to be generated', default=10)
  #parser.add_argument('--time', type=str, dest='time', help='allowed time for each job of each step in hours, including dumper, example 01,02,01,01,04 for 1h for step1, 2h for step2, 1h for step3, 1h for step4, 4h for dumper', default='02,02,05,01,05')
  parser.add_argument('--time', type=str, dest='time', help='allowed time for each job', default='02')
  parser.add_argument('--njobs', type=int, dest='njobs', help='number of parallel jobs to submit', default=10)
  parser.add_argument('--mass', type=int, dest='mass', help='mass to generate (GeV)', default=1)
  parser.add_argument('--ctau', type=int, dest='ctau', help='ctau to generate (mm)', default=100)
  parser.add_argument('--domultithread', dest='domultithread', help='run multithreaded', action='store_true', default=False)
  parser.add_argument('--domultijob', dest='domultijob', help='run several separate jobs', action='store_true', default=False)
  parser.add_argument('--dosubmit', dest='dosubmit', help='submit to slurm', action='store_true', default=False)

  return parser.parse_args()

if __name__ == "__main__":

  opt = getOptions()

  ##############################
  # job configurations
  #############################
  if opt.domultijob and opt.njobs <= 1: raise RuntimeError('when running multiple jobs, the number of parallel jobs should be larger than 1')
  if opt.domultijob and opt.nevts % opt.njobs != 0 and not opt.dorecofromeos: raise RuntimeError('cannot split events in njobs evenly, please change njobs / nevts')
  if opt.domultijob and opt.domultithread: raise RuntimeError('either multijob or multithread, choose, otherwise seed for generation will repeat')
  njobs = opt.njobs if opt.domultijob else 1
  nevtsjob = opt.nevts if not opt.domultijob else opt.nevts/opt.njobs
  prodLabel = '{v}_m{m}_ctau{ct}_n{n}_njt{nj}'.format(v=opt.ver,n=opt.nevts,m=opt.mass,ct=opt.ctau,nj=njobs)
  nthr = 8 if opt.domultithread else 1
  user = os.environ["USER"]

  ##############################
  # create local production directory 
  ############################# 
  if not os.path.isdir(prodLabel):
    os.system('mkdir -p ./{}'.format(prodLabel))
  # otherwise will overwrite 
  if not os.path.isdir(prodLabel+'/logs'):
    os.system('mkdir {}/logs'.format(prodLabel))

  ############################
  # template for batch submission 
  ############################ 
  template = makeTemplate()
  template = template.format(m=opt.mass, ctau=opt.ctau, hh=opt.time, arr='1-{}'.format(njobs),pl=prodLabel,user=user,jop='BPH_mod_cfg.py',nevtsjob=nevtsjob,nthr=nthr)

  launcherFile = '{}/slurm_step1.sh'.format(prodLabel)
  with open(launcherFile, 'w') as f:
    f.write(template)


  #
  print('==> Created directory for submission\n    {}'.format(prodLabel))
  if opt.dosubmit:
    os.chdir(prodLabel)
    os.system('sbatch slurm_step1.sh')
    os.chdir("../")

  ##  
  with open('{}/cfg.txt'.format(prodLabel), 'w') as f:
    f.write('Run genHelper.py with following configurations\n')
    for k,v in sorted(vars(opt).items()):
      f.write('{:15s}: {:10s}\n'.format(str(k),str(v)))
    
  



