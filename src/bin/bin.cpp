/*! \file bin.cpp
    \brief This code is the main function for the executable.
 
    In general, this file contains the documentation start page code. It also
    provides all the required code for running ProSHADE binary, that is it
    creates the ProSHADE_settings object, it reads the command line settings
    into it and it then passes this filled object to the ProSHADE_run constructor,
    which does all the computations required by the settings object. The ProSHADE_run
    object then also holds the results.
 
    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.3
    \date      JUN 2020
 */

//==================================================== DOxygen main page specifications

/*! \mainpage ProSHADE Documentation
 *
 * \section intro Introduction
 *
 * ProSHADE is a C++ language library and an associated tool providing functionalities for working with structural biology molecular structures. The library implements functions for computing shape-wise structural
 * distances between pairs of molecules, detecting symmetry over the centre of mass of a single structure, map re-sizing as well as matching density maps and PDB coordinate files into one another.
 * The executable implemented in the bin.cpp file then allows easy access to these functionalities without the need for library linking, while the python modules provide easy access to the functionality from
 * the python language. For help on how the executable should be used, refer to the -h option of it. For more details about the functionalities, see below.
 *
 * \section index Index
 *
 * 1) \ref intro
 *
 * 2) \ref index
 *
 * 3) \ref install
 *
 * \section install Installation
 *
 * The installation of the ProSHADE software should be done using the CMake system.
 *
 *
 */

//==================================================== ProSHADE
#include "../proshade/ProSHADE.hpp"

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
    return                                            ( 0 );
}
