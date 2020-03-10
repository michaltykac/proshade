import sys
sys.path.append("/Users/mysak/BioCEV/proshade/exp2/install/python3")
import proshade

pSet = proshade.ProSHADE_settings ()

pSet.addStructure ( "/Users/mysak/LMB/proshade/exp/demo/C2.pdb" )
pSet.addStructure ( "/Users/mysak/LMB/proshade/exp/demo/C2.pdb" )
pSet.setResolution ( 10 )
pSet.task = proshade.Distances

pRun = proshade.ProSHADE_run ( pSet )

del pRun
del pSet
