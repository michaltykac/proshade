/*! \file proshadeBinary.cpp
    \brief This code is the same as the main function for the executable, but it uses the dynamic library linking instead.
 
    This is a test file, which was created to test how ProSHADE could be used in the case where its
    dynamic library is linked to the code instead of it being compiled directly from the code.
 
    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.6.2
    \date      DEC 2021
 */

//==================================================== ProSHADE
#include "../../src/proshade/ProSHADE.hpp"

//==================================================== Main
int main ( int argc, char **argv )
{
    //================================================ Create the settings object and parse the command line arguments
    ProSHADE_settings* settings                       = new ProSHADE_settings ( );
    settings->getCommandLineParams                    ( argc, argv );
    
    //================================================ Execute
    ProSHADE_run *run                                 = new ProSHADE_run ( settings );

    //================================================ Release the settings object
    delete settings;
    
    //================================================ Release the executive object
    delete run;
    
    //================================================ DONE
    return                                            ( EXIT_SUCCESS );
}
