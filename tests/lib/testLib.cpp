/*! \file testLib.cpp
 \brief This code is the same as the main function for the executable, but it uses the dynamic library linking instead.
 
 This is a test file, which was created to test how ProSHADE could be used in the case where its
  dynamic library is linked to the code instead of it being compiled directly from the code.
 
 This file is part of the ProSHADE library for calculating
 shape descriptors and symmetry operators of protein structures.
 This is a prototype code, which is by no means complete or fully
 tested. Its use is at your own risk only. There is no quarantee
 that the results are correct.
 
 \author    Michal Tykac
 \author    Garib N. Murshudov
 \version   0.7.3
 \date      MAY 2020
 */

//============================================ ProSHADE
#include "../../src/proshade/ProSHADE.hpp"

//============================================ Main
int main ( int argc, char **argv )
{
    //======================================== Create the settings object and parse the command line arguments
    ProSHADE_settings* settings               = new ProSHADE_settings ( );
    settings->getCommandLineParams            ( argc, argv );
    
    //======================================== Execute
    ProSHADE_run *run                             = new ProSHADE_run ( settings );

    //======================================== Release the settings object
    delete settings;
    
    //======================================== Release the executive object
    delete run;
    
    //======================================== DONE
    return 0;
}
