from genAnalysis import *
# import all functions from gen analysis and simply use them in a new main...


if __name__ == "__main__":

  gROOT.SetBatch(True)
  ROOT.TH1.SetDefaultSumw2()

  opt = getOptions()

  path = './outputfiles/' + opt.pl + '/mass{m}_ctau{ctau}_miniGenTree.root'
  label = path.split('/')[2]

  ################
  points = [
#   Point(mass=1.0,ctau=None,vv=5e-03),
#   Point(mass=1.0,ctau=None,vv=1e-03),
#   Point(mass=1.0,ctau=None,vv=5e-04),
#   Point(mass=1.0,ctau=None,vv=1e-04),
#   Point(mass=1.0,ctau=None,vv=5e-05),
#   Point(mass=1.0,ctau=None,vv=1e-05),
#   Point(mass=1.0,ctau=None,vv=5e-06),
#   Point(mass=1.0,ctau=None,vv=1e-06),
#   Point(mass=1.0,ctau=None,vv=5e-07),
#   
#   Point(mass=1.5,ctau=None,vv=5e-03),
#   Point(mass=1.5,ctau=None,vv=1e-03),
#   Point(mass=1.5,ctau=None,vv=5e-04),
#   Point(mass=1.5,ctau=None,vv=1e-04),
#   Point(mass=1.5,ctau=None,vv=5e-05),
#   Point(mass=1.5,ctau=None,vv=1e-05),
#   Point(mass=1.5,ctau=None,vv=5e-06),
#   Point(mass=1.5,ctau=None,vv=1e-06),
#   Point(mass=1.5,ctau=None,vv=5e-07),
#   
#   Point(mass=2.0,ctau=None,vv=5e-03),
#   Point(mass=2.0,ctau=None,vv=1e-03),
#   Point(mass=2.0,ctau=None,vv=5e-04),
#   Point(mass=2.0,ctau=None,vv=1e-04),
#   Point(mass=2.0,ctau=None,vv=5e-05),
#   Point(mass=2.0,ctau=None,vv=1e-05),
#   Point(mass=2.0,ctau=None,vv=5e-06),
#   Point(mass=2.0,ctau=None,vv=1e-06),
#   Point(mass=2.0,ctau=None,vv=5e-07),    
#      Point(mass=0.5,ctau=None,vv=1e-04,isrw=False),
#      Point(mass=1.0,ctau=None,vv=1e-04,isrw=False),
      #Point(mass=1.5,ctau=None,vv=1e-03,isrw=False),
#      Point(mass=2.0,ctau=None,vv=1e-04,isrw=False),
#      Point(mass=2.5,ctau=None,vv=1e-04,isrw=False),
#      Point(mass=3.0,ctau=None,vv=1e-04,isrw=False),

      Point(mass=1.0,ctau=None,vv=1e-05,isrw=False),
      Point(mass=2.0,ctau=None,vv=1e-05,isrw=False),
      Point(mass=3.0,ctau=None,vv=1e-05,orig_vv=5e-05,isrw=True),


  ]

  for p in points:
    #fn = path.format(m=p.mass,ctau=p.orig_ctau)
    fns = glob(path.format(m=p.mass,ctau=p.orig_ctau))
    s = Sample(mass=p.mass, ctau=p.ctau, vv=p.vv, infileNames=fns, isrw=p.isrw, orig_vv=p.orig_vv, label=label)
    s.fillFilterEff(dostamp=True)
   
   
