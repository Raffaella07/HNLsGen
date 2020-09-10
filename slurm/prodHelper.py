'''
Production helper for generating private samples for B-initiated HNLs
'''

import sys
import os
import subprocess

sys.path.append('../python/') 
from common import Point


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
    self.prodLabel = '{v}_n{n}_njt{nj}'.format(v=self.ver,n=self.nevts,nj=self.njobs)
    self.nthr = 8 if self.domultithread else 1
    self.npremixfiles = 20 # the number of events per file being 1200
    self.user = os.environ["USER"]
    self.jop1 = 'step1.py'
    self.jop2 = 'step2.py'
    self.jop3 = 'step3.py'
    self.jop4 = 'step4.py'
   

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


  '''
  def getStartDir(self):
    if self.user == 'mratti':
      startdir = '/work/mratti/GEN_HNL_newPythia/CMSSW_10_2_3/src/HNLsGen/'
    elif self.user == 'anlyon':
      startdir = '/t3home/anlyon/BHNL/CMSSW_10_2_3/src/HNLsGen'
    else:
      startdir = '`pwd`'
    return startdir
  '''

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


  def appendTemplate(self, jopa, jopb, nthr, nevtsjob, npremixfiles=0):
    if self.dogenonly:
      addlines = ''
    else:
      labela = jopa[0:len(jopa)-3]
      labelb = jopb[0:len(jopb)-3]
      if jopa == 'step2.py':
        command = 'cmsRun {jopa} randomizePremix=1 nPremixFiles={nprx} nThr={nthr} inputFile=BPH-{lblb}_numEvent{nevtsjob}.root outputFile=BPH-{lbla}.root seedOffset=$SLURM_ARRAY_TASK_ID'.format(jopa=jopa, nprx=npremixfiles, nthr=nthr, lblb=labelb, nevtsjob=nevtsjob, lbla=labela)
      else:
        command = 'cmsRun {jopa} nThr={nthr} inputFile=BPH-{lblb}.root outputFile=BPH-{lbla}.root seedOffset=$SLURM_ARRAY_TASK_ID'.format(jopa=jopa, nthr=1, lblb=labelb, lbla=labela)

      addlines = [
        '',
        '### {lbla} ###',
        'echo "Going to copy cmsdriver to work dir"',
        'cp $STARTDIR/cmsDrivers/{jopa} $WORKDIR/. ',
        'cp $STARTDIR/slurm/{lbldir}/{jopa} $WORKDIR/. ',
        'echo "Going to run {lbla}"',
        'DATE_START_{lbla}=`date +%s`',
        '{command}',
        'DATE_END_{lbla}=`date +%s`',
        'echo "Finished running {lbla}"',
        'echo "Content of current directory"',
        'ls -al',
        'echo ""',
        'RUNTIME_{lbla}=$((DATE_END_{lbla}-DATE_START_{lbla}))',
        'echo "Intermediate allclock running time {lbla}: $RUNTIME_{lbla} s"',
        '',
      ]
      if jopa == 'step4.py':
        addlines += [
        '',
        'echo "Going to copy output to result directory"',
        'xrdcp -f $WORKDIR/BPH-{lbla}.root $OUTSEPREFIX/$SERESULTDIR/{lbla}_nj$SLURM_ARRAY_TASK_ID".root"',
        '',
        ]

      addlines = '\n'.join(addlines)
      addlines = addlines.format(jopa=jopa, nthr=nthr, lblb=labelb, nevtsjob=nevtsjob, lbldir=self.prodLabel, lbla=labela, command=command)

    return addlines


  def makeTimeStamp(self):
    if self.dogenonly:
      timestamp=[
        'RUNTIME_step1=$((DATE_END_step1-DATE_START_step1))',
        'echo "Wallclock running time: $RUNTIME s"'
      ]
    else:
      timestamp=[
        'RUNTIME_step1=$((DATE_END_step1-DATE_START_step1))',
        'RUNTIME_step2=$((DATE_END_step2-DATE_START_step2))',
        'RUNTIME_step3=$((DATE_END_step3-DATE_START_step3))',
        'RUNTIME_step4=$((DATE_END_step4-DATE_START_step4))',
        'RUNTIME_tot=$((DATE_END_step4-DATE_START_step1))',
        'echo "Wallclock running time step1: $RUNTIME_step1 s"',
        'echo "                       step2: $RUNTIME_step2 s"',
        'echo "                       step3: $RUNTIME_step3 s"',
        'echo "                       step4: $RUNTIME_step4 s"',
        'echo "                         tot: $RUNTIME_tot s"'
      ]
    return '\n'.join(timestamp)


  def makeTemplates(self):
    for p in self.points:
      template = [
        '#!/bin/bash',
        '',
        '#SBATCH -J prod_m{m}_ctau{ctau}',
        '#SBATCH -o logs/prod_mass{m}_ctau{ctau}_%a.log', 
        '#SBATCH -e logs/prod_mass{m}_ctau{ctau}_%a.log',
        '#SBATCH -p wn',
        '#SBATCH -t {hh}:00:00',
        '#SBATCH --mem {mem}',
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
        'echo "Going to copy cmsdriver to work dir"',
        'cp $STARTDIR/slurm/{lbldir}/{jop1} $WORKDIR/. ',
        'echo ""',
        '',
        #'echo "Going to copy the evtGen particle data file"',  # currently not needed, cannot use local copy of file
        #'mkdir -p $WORKDIR/HNLsGen/evtGenData/',
        #'cp ../../evtGenData/evt_2014_mass{m}_ctau{ctau}.pdl $WORKDIR/HNLsGen/evtGenData/.',
        #'echo ""',
        #'',
        'cd $WORKDIR',
        'pwd',
        'echo "Going to run step1"',
        'DATE_START_step1=`date +%s`',
        'cmsRun {jop1} maxEvents={nevtsjob} nThr={nthr} mass={m} ctau={ctau} outputFile=BPH-step1.root seedOffset=$SLURM_ARRAY_TASK_ID',
        'DATE_END_step1=`date +%s`',
        'echo "Finished running step1"',
        'echo "Content of current directory"',
        'ls -al',
        'echo ""',
        '',
        'echo "Going to copy output to result directory"',
        'xrdcp -f $WORKDIR/BPH-step1_numEvent{nevtsjob}.root $OUTSEPREFIX/$SERESULTDIR/step1_nj$SLURM_ARRAY_TASK_ID".root"',
        '',
        '{addstep2}',
        '{addstep3}',
        '{addstep4}',
        'echo ""',
        'echo "Cleaning up $WORKDIR"',
        'rm -rf $WORKDIR',
        '{timestamp}',
        'cd $STARTDIR',
      ]
      template = '\n'.join(template)
      template = template.format(
          m=p.mass,
          ctau=p.ctau,
          hh=self.time,
          mem=self.mem,
          lbldir=self.prodLabel,
          arr='1-{}'.format(self.njobs),
          pl=self.prodLabel,
          user=self.user,
          jop1=self.jop1,
          nevtsjob=self.nevtsjob,
          nthr=self.nthr,
          jop2=self.jop2,
          jop3=self.jop3,
          jop4=self.jop4,
          addstep2=self.appendTemplate(self.jop2,self.jop1,self.nthr,self.nevtsjob,self.npremixfiles),
          addstep3=self.appendTemplate(self.jop3,self.jop2,self.nthr,self.nevtsjob),
          addstep4=self.appendTemplate(self.jop4,self.jop3,self.nthr,self.nevtsjob),
          timestamp=self.makeTimeStamp()
          )
      launcherFile = '{pl}/slurm_mass{m}_ctau{ctau}_prod.sh'.format(pl=self.prodLabel,m=p.mass,ctau=p.ctau)
      with open(launcherFile, 'w') as f:
        f.write(template)
    
    print('===> Created templates for batch submission\n')
    
  def writeCfg(self):
    with open('{}/cfg.txt'.format(self.prodLabel), 'w') as f:
      f.write('Run genHelper.py with following options\n')
      for k,v in sorted(vars(self.opt).items()):
        f.write('{:15s}: {:10s}\n'.format(str(k),str(v)))
    os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop1,pl=self.prodLabel))
    if not self.dogenonly:
      os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop2,pl=self.prodLabel))
      os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop3,pl=self.prodLabel))
      os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop4,pl=self.prodLabel))
      

  def submit(self):
    os.chdir(self.prodLabel)
    for p in self.points:
      os.system('sbatch slurm_mass{m}_ctau{ctau}_prod.sh'.format(m=p.mass,ctau=p.ctau))
    os.chdir('../')
    print('')
    print('===> Submitted {n} job arrays for {pl}\n'.format(n=len(self.points),pl=self.prodLabel))
  

def getOptions():

  # convention: no capital letters

  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for B-initiated HNL signals', add_help=True)

  parser.add_argument('-v','--ver', type=str, dest='ver', help='version of production, e.g. V00_v00', default='V00_v00')
  parser.add_argument('-n','--nevts', type=int, dest='nevts', help='total number of events to be generated', default=10)
  parser.add_argument('--time', type=str, dest='time', help='allowed time for each job', default='08')
  parser.add_argument('--mem', type=str, dest='mem', help='allowed memory for each job in [MB]', default='2500')
  parser.add_argument('--njobs', type=int, dest='njobs', help='number of parallel jobs to submit', default=10)
  parser.add_argument('--points', type=str, dest='pointFile', help='name of file contaning information on scan to be run', default='points.py')
  parser.add_argument('--domultithread', dest='domultithread', help='run multithreaded', action='store_true', default=False)
  parser.add_argument('--domultijob', dest='domultijob', help='run several separate jobs', action='store_true', default=False)
  parser.add_argument('--dosubmit', dest='dosubmit', help='submit to slurm', action='store_true', default=False)
  parser.add_argument('--dogenonly', dest='dogenonly', help='produce sample until gen', action='store_true', default=False)

  return parser.parse_args()

if __name__ == "__main__":

  opt = getOptions()

  job = Job(opt)

  job.makeProdDir()

  job.makeEvtGenData()

  job.makeTemplates()

  job.writeCfg()   

  if opt.dosubmit:
    job.submit()

