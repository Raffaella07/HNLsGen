'''
Production helper for generating private samples for B-initiated HNLs
=> for the moment up to gen
=> would like to add a script that checks the status of the production
'''

import sys
import os
import subprocess

from python.common import Point

class Job(object):
  def __init__(self,opt):

    self.opt = opt
    for k,v in sorted(vars(opt).items()):
      setattr(self,k,v)

    if not os.path.isfile(self.pointFile): raise RuntimeError('Provided file for points to scan does not exist, {}'.format(self.pointFile))
    ps = __import__(self.pointFile.split('.py')[0])
    self.points = ps.points
    if self.domultijob and self.njobs <= 1: raise RuntimeError('when running multiple jobs, the number of parallel jobs should be larger than 1')
    if self.domultijob and self.nevts % self.njobs != 0: raise RuntimeError('cannot split events in njobs evenly, please change njobs / nevts')
    if self.domultijob and self.domultithread: raise RuntimeError('either multijob or multithread, choose, otherwise seed for generation will repeat')
    self.njobs = self.njobs if self.domultijob else 1
    self.nevtsjob = self.nevts if not self.domultijob else self.nevts/self.njobs
    #self.prodLabel = '{v}_mass{m}_ctau{ct}_n{n}_njt{nj}'.format(v=self.ver,n=self.nevts,m=self.mass,ct=self.ctau,nj=self.njobs)
    self.prodLabel = '{v}_n{n}_njt{nj}'.format(v=self.ver,n=self.nevts,nj=self.njobs)
    self.nthr = 8 if self.domultithread else 1
    self.user = os.environ["USER"]
    self.jop = 'BPH_mod_cfg.py'
   
  def makeProdDir(self):
    if not os.path.isdir(self.prodLabel):
      os.system('mkdir -p ./{}'.format(self.prodLabel))
    # otherwise will overwrite 
    if not os.path.isdir(self.prodLabel+'/logs'):
      os.system('mkdir {}/logs'.format(self.prodLabel)) 
    print('===> Created directory for submission {}\n'.format(self.prodLabel))

    print('===> Points to be run')
    for p in self.points:
      p.stamp()
    print('')

  def makeEvtGenData(self):
    for p in self.points:      
      hnl_lines = 'add  p Particle  hnl                          9900015  {:.7e}  0.0000000e+00  0.0000000e+00     0     1  {:.7e}    9900015\nadd  p Particle  anti_hnl                    -9900015  {:.7e}  0.0000000e+00  0.0000000e+00     0     1  {:.7e}   -9900015\n'.format(p.mass,p.ctau,p.mass,p.ctau)

      with open('../evtGenData/evt_2014.pdl', 'r') as fin:
        contents = fin.readlines()
        contents.insert(4, hnl_lines)
        contents = ''.join(contents)
      with open('../evtGenData/evt_2014_mass{m}_ctau{ctau}.pdl'.format(m=p.mass, ctau=p.ctau), 'w') as fout:
        fout.write(contents)
    print('===> Created evtGen particle property files\n')

  def makeTemplates(self):
    for p in self.points:
      template = [
        '#!/bin/bash',
        '',
        '#SBATCH -J step1_m{m}_ctau{ctau}',
        '#SBATCH -o logs/step1_mass{m}_ctau{ctau}_%a.log', 
        '#SBATCH -e logs/step1_mass{m}_ctau{ctau}_%a.log',
        #'#SBATCH -o logs/step1_%A_%a.log', 
        #'#SBATCH -e logs/step1_%A_%a.log',
        '#SBATCH -p wn',
        '#SBATCH -t {hh}:00:00',
        '#SBATCH --array={arr}',
        '#SBATCH --ntasks=1',
        '#SBATCH --account=t3',
        '',
        'DIRNAME="{pl}"/mass{m}_ctau{ctau}/',
        'STARTDIR=$CMSSW_BASE/src/HNLsGen/', # Will take the cmssw version used at submissio time
        'TOPWORKDIR="/scratch/{user}/"',
        'JOBDIR="gen_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"', # MIND THE PARENTHESIS
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
        'xrdfs $T3CACHE mkdir -p $SERESULTDIR',
        'echo ""',
        '',
        'echo "Going to copy cmsdriver to work dir"', # copy from here...
        'cd -', # going back to submission directory
        'cp {jop} $WORKDIR/. ',
        'echo ""',
        '',
        #'echo "Going to copy the evtGen particle data file"',  # currently not needed, cannot use local copy of file
        #'mkdir -p $WORKDIR/HNLsGen/evtGenData/',
        #'cp $STARTDIR/evtGenData/evt_2014_mass{m}_ctau{ctau}.pdl $WORKDIR/HNLsGen/evtGenData/.',
        #'echo ""',
        #'',
        'cd $WORKDIR',
        'pwd',
        'echo "Going to run"',
        'DATE_START=`date +%s`',
        'cmsRun {jop} maxEvents={nevtsjob} nThr={nthr} mass={m} ctau={ctau} outputFile=BPH-test.root seedOffset=$SLURM_ARRAY_TASK_ID',
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
      template = template.format(m=p.mass, ctau=p.ctau, hh=self.time, arr='1-{}'.format(self.njobs),pl=self.prodLabel,user=self.user,jop=self.jop,nevtsjob=self.nevtsjob,nthr=self.nthr)
      launcherFile = '{pl}/slurm_mass{m}_ctau{ctau}_step1.sh'.format(pl=self.prodLabel,m=p.mass,ctau=p.ctau)
      with open(launcherFile, 'w') as f:
        f.write(template)
    
    print('===> Created templates for batch submission\n')
    

  def writeCfg(self):
    with open('{}/cfg.txt'.format(self.prodLabel), 'w') as f:
      f.write('Run genHelper.py with following options\n')
      for k,v in sorted(vars(self.opt).items()):
        f.write('{:15s}: {:10s}\n'.format(str(k),str(v)))
    os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop,pl=self.prodLabel))

  def submit(self):
    os.chdir(self.prodLabel)
    for p in self.points:
      os.system('sbatch slurm_mass{m}_ctau{ctau}_step1.sh'.format(m=p.mass,ctau=p.ctau))
    os.chdir('../')
    print('')
    print('===> Submitted {n} job arrays for {pl}\n'.format(n=len(self.points),pl=self.prodLabel))
  

def getOptions():

  # convention: no capital letters

  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for B-initiated HNL signals', add_help=True)

  parser.add_argument('-v','--ver', type=str, dest='ver', help='version of production, e.g. V00_v00', default='V00_v00')
  parser.add_argument('-n','--nevts', type=int, dest='nevts', help='total number of events to be generated', default=10)
  #parser.add_argument('--time', type=str, dest='time', help='allowed time for each job of each step in hours, including dumper, example 01,02,01,01,04 for 1h for step1, 2h for step2, 1h for step3, 1h for step4, 4h for dumper', default='02,02,05,01,05')
  parser.add_argument('--time', type=str, dest='time', help='allowed time for each job', default='02')
  parser.add_argument('--njobs', type=int, dest='njobs', help='number of parallel jobs to submit', default=10)
  parser.add_argument('--points', type=str, dest='pointFile', help='name of file contaning information on scan to be run', default='points.py')
  #parser.add_argument('--mass', type=int, dest='mass', help='mass to generate (GeV)', default=1)
  #parser.add_argument('--ctau', type=int, dest='ctau', help='ctau to generate (mm)', default=100)
  parser.add_argument('--domultithread', dest='domultithread', help='run multithreaded', action='store_true', default=False)
  parser.add_argument('--domultijob', dest='domultijob', help='run several separate jobs', action='store_true', default=False)
  parser.add_argument('--dosubmit', dest='dosubmit', help='submit to slurm', action='store_true', default=False)

  return parser.parse_args()

if __name__ == "__main__":

  opt = getOptions()

  job = Job(opt)
  #for k,v in job.__dict__.items():
  #  print k,v
  job.makeProdDir()

  job.makeEvtGenData()

  job.makeTemplates()

  job.writeCfg()   

  if opt.dosubmit:
    job.submit()

