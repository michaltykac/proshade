##############################################
##############################################
#   \file simpleAccess_distances.py
#   \brief This code demonstrates the usage of the ProSHADE tool in the simple mode for distances functionality.
#
#   This code demonstrates how the ProSHADE python module can be loaded,
#   the settings object created and filled (including all optional values),
#   supplied to the run object and how to three descriptor values can be
#   obtained programatically.
#
#   Copyright by Michal Tykac and individual contributors. All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#   1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#   2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#   3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
#   This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In     no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data     or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility     of such damage.
#
#   \author    Michal Tykac
#   \author    Garib N. Murshudov
#   \version   0.7.5.0
#   \date      DEC 2020
##############################################
##############################################


##############################################
### Import system modules
import sys
import numpy

### Import ProSHADE from non-system folder (local installation assumed)
sys.path.append                               ( "/Users/mysak/BioCEV/proshade/build" )
import pyproshade as proshade

ps = proshade.ProSHADE_settings ( proshade.Symmetry )

ps.addStructure("/Users/mysak/LMB/proshade/exp/demo/C2.pdb")
ps.requestedResolution = 8.0
ps.verbose = 2

rn = proshade.ProSHADE_run ( ps )

symType = rn.getSymmetryType()
symFold = rn.getSymmetryFold()
symAxis = rn.getSymmetryAxis(0)
allCs   = rn.getAllCSyms()
