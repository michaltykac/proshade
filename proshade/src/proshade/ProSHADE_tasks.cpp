/*! \file ProSHADE_tasks.cpp
    \brief This source file contains the task functions, which drive the computation of a specific task.
 
    The funcions in this source file are of two types, firstly, there are the task functions, which are responsible for executing a particular task, that is executing a set of functions in the required order so that the
    task computation is achieved. The second type of functions in this source file are the sanity funnctions, which test that all the information required for a particular task were supplied by the user.
 
    Copyright by Michal Tykac and individual contributors. All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
 
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.5.4
    \date      MAR 2021
 */

//==================================================== ProSHADE
#include "ProSHADE_tasks.hpp"

/*! \brief The re-boxing task driver function.
 
    This function is called to proceed with the map re-boxing task according to the information placed in
    the settings object passed as the first argument.
 
    \param[in] settings ProSHADE_settings object specifying the details of how re-boxing should be done.
    \param[in] originalBounds Vector to which the original map boundaries of each re-boxed map will be saved into.
    \param[in] reboxedBounds Vector to which the re-boxed map boundaries of each re-boxed map will be saved into.
    \param[in] manipulatedMaps Vector to which the map values of each re-boxed map will be saved into.
 */
void ProSHADE_internal_tasks::MapManipulationTask ( ProSHADE_settings* settings, std::vector < proshade_signed* >* originalBounds, std::vector < proshade_signed* >* reboxedBounds, std::vector < proshade_double* >* manipulatedMaps )
{
    //================================================ Check the settings are complete and meaningful
    checkMapManipulationSettings                      ( settings );
    
    //================================================ For all inputted structures
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( settings->inputFiles.size() ); iter++ )
    {
        //============================================ Create a data object
        ProSHADE_internal_data::ProSHADE_data* strToRebox = new ProSHADE_internal_data::ProSHADE_data ( settings );
        
        //============================================ Read in the file
        strToRebox->readInStructure                   ( settings->inputFiles.at(iter), iter, settings );
        
        //============================================ Save the original boundaries
        ProSHADE_internal_misc::deepCopyBoundsSigPtrVector ( originalBounds, strToRebox->getXFromPtr(), strToRebox->getXToPtr(), strToRebox->getYFromPtr(), strToRebox->getYToPtr(), strToRebox->getZFromPtr(), strToRebox->getZToPtr() );
        
        //============================================ Internal data processing  (COM, norm, mask, extra space)
        strToRebox->processInternalMap                ( settings );
        
        //============================================ Create new structure for re-boxing
        ProSHADE_internal_data::ProSHADE_data* reBoxStr = new ProSHADE_internal_data::ProSHADE_data ( settings );
        
        //============================================ Re-box map, if need be
        if ( settings->reBoxMap )
        {
            //======================================== Find non-zero bounds
            proshade_signed* nonZeroBounds            = new proshade_signed[6];
            strToRebox->getReBoxBoundaries            ( settings, nonZeroBounds );
            
            //============================================ Create new structure from the bounds
            strToRebox->createNewMapFromBounds        ( settings, reBoxStr, nonZeroBounds );
            
            //======================================== Release memory
            delete[] nonZeroBounds;
        }
        
        //============================================ Save the modified structure
        std::stringstream ss;
        ss << settings->outName << "_" << iter << ".map";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Saving the re-boxed map into " + ss.str() );
        if ( settings->reBoxMap )  { reBoxStr->writeMap ( ss.str() ); }
        else { strToRebox->writeMap ( ss.str() ); }
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 2, "Structure saved." );
        
        //============================================ Save the re-boxed boundaries
        ProSHADE_internal_misc::deepCopyBoundsSigPtrVector ( reboxedBounds, reBoxStr->getXFromPtr(), reBoxStr->getXToPtr(), reBoxStr->getYFromPtr(),
                                                             reBoxStr->getYToPtr(), reBoxStr->getZFromPtr(), reBoxStr->getZToPtr() );
        
        //============================================ Save the map
        proshade_double* mapCopy                      = NULL;
        reBoxStr->deepCopyMap                         ( mapCopy, settings->verbose );
        ProSHADE_internal_misc::addToDblPtrVector     ( manipulatedMaps, mapCopy );
        
        //============================================ Release memory
        delete strToRebox;
        delete reBoxStr;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief The re-boxing settings checks.
 
    This function is called to check the settings object for having all the required information for
    the Re-Boxing task to proceed.
 
    \param[in] settings ProSHADE_settings object specifying the details of how re-boxing should be done.
 */
void ProSHADE_internal_tasks::checkMapManipulationSettings ( ProSHADE_settings* settings )
{
    //================================================ Is there a single file for processing?
    if ( settings->inputFiles.size () == 0 )
    {
        throw ProSHADE_exception ( "There is no input structure for map manipulation.", "EB00002", __FILE__, __LINE__, __func__, "The ProSHADE_settings object does not contain any\n                    : structure that could be manipulated. Please supply exactly\n                    : one structure using the addStructure() function." );
    }
    
    //================================================ Is the file type MAP? Warning if not
    if ( ProSHADE_internal_io::isFilePDB ( settings->inputFiles.at(0) ) )
    {
        ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! The input file is not of the MAP (MRC) format. Will output re-boxed map, but beware that this is simple PDB->MAP conversion and REFMAC5 should be used to compute more appropriate maps.", "WB00004" );
        
        //============================================ No resolution for PDB? Problem...
        if ( settings->requestedResolution == 0.0 )
        {
            throw ProSHADE_exception ( "No resolution given for PDB file re-boxing.", "EB00011", __FILE__, __LINE__, __func__, "The ProSHADE_settings object does not contain any\n                    : resolution value. However, resolution is required when\n                    : re-boxing structures read from PDB files. Please supply\n                    : the resolution value using the setResolution() function." );
        }
    }
    
    //================================================ Is there output file name?
    if ( settings->outName == "" )
    {
        throw ProSHADE_exception ( "No output file name.", "EB00016", __FILE__, __LINE__, __func__, "There is no output file name set in the settings object.\n                    : Please supply the file name to where the re-boxed map\n                    : should be saved using the setOutputFilename() function." );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief The distances computation task driver function.
 
    This function is called to proceed with the distances computation task according to the information placed in
    the settings object passed as the first argument.
 
    \param[in] settings ProSHADE_settings object specifying the details of how distances computation should be done.
    \param[in] enLevs Pointer to vector where all energy levels distances are to be saved into.
    \param[in] trSigm Pointer to vector where all trace sigma distances are to be saved into.
    \param[in] rotFun Pointer to vector where all rotation function distances are to be saved into.
 */
void ProSHADE_internal_tasks::DistancesComputationTask ( ProSHADE_settings* settings, std::vector< proshade_double >* enLevs, std::vector< proshade_double >* trSigm, std::vector< proshade_double >* rotFun )
{
    //================================================ Check the settings are complete and meaningful
    checkDistancesSettings                            ( settings );
    
    //================================================ Create a data object
    ProSHADE_internal_data::ProSHADE_data* compareAgainst  = new ProSHADE_internal_data::ProSHADE_data ( settings );
    
    //================================================ Read in the structure all others will be compared to
    compareAgainst->readInStructure                   ( settings->inputFiles.at(0), 0, settings );
    
    //================================================ Internal data processing  (COM, norm, mask, extra space)
    compareAgainst->processInternalMap                ( settings );
    
    //================================================ Map to sphere
    compareAgainst->mapToSpheres                      ( settings );
    
    //================================================ Get spherical harmonics
    compareAgainst->computeSphericalHarmonics         ( settings );
    
    //================================================ Now, for each other structure
    for ( proshade_unsign iter = 1; iter < static_cast<proshade_unsign> ( settings->inputFiles.size() ); iter++ )
    {
        //============================================ Create a data object
        ProSHADE_internal_data::ProSHADE_data* compareChanging = new ProSHADE_internal_data::ProSHADE_data ( settings );

        //============================================ Read in the compared structure
        compareChanging->readInStructure              ( settings->inputFiles.at(iter), iter, settings );

        //============================================ Internal data processing  (COM, norm, mask, extra space)
        compareChanging->processInternalMap           ( settings );

        //============================================ Map to sphere
        compareChanging->mapToSpheres                 ( settings );
        
        //============================================ Get spherical harmonics
        compareChanging->computeSphericalHarmonics    ( settings );
        
        //============================================ Get distances
        proshade_double enLevDist                     = 0.0;
        if ( settings->computeEnergyLevelsDesc ) { enLevDist  = ProSHADE_internal_distances::computeEnergyLevelsDescriptor ( compareAgainst, compareChanging, settings ); }
        else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Energy levels distance computation not required." ); }
        
        proshade_double trSigmDist                    = 0.0;
        if ( settings->computeTraceSigmaDesc   ) { trSigmDist = ProSHADE_internal_distances::computeTraceSigmaDescriptor ( compareAgainst, compareChanging, settings ); }
        else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Trace sigma distance computation not required." ); }
        
        proshade_double rotFunDist                    = 0.0;
        if ( settings->computeRotationFuncDesc ) { rotFunDist = ProSHADE_internal_distances::computeRotationunctionDescriptor ( compareAgainst, compareChanging, settings ); }
        else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Rotation function distance computation not required." ); }
        
        //============================================ Save results to the run object
        ProSHADE_internal_misc::addToDoubleVector     ( enLevs, enLevDist  );
        ProSHADE_internal_misc::addToDoubleVector     ( trSigm, trSigmDist );
        ProSHADE_internal_misc::addToDoubleVector     ( rotFun, rotFunDist );
        
        //============================================ Report results
        ReportDistancesResults                        ( settings, settings->inputFiles.at(0), settings->inputFiles.at(iter), enLevDist, trSigmDist, rotFunDist );
        
        //============================================ Release the memory
        delete compareChanging;
    }
    

    //================================================ Release memory
    delete compareAgainst;
    
    //================================================ Done
    return ;
    
}

/*! \brief Simple function for reporting the distances computation results.
 
    \param[in] settings ProSHADE_settings object specifying the details of how distances computation should be done.
    \param[in] str1 The name of the structure to which all other structures are to be compared to.
    \param[in] str2 The name of the structure which is compared to str1.
    \param[in] enLevDist The value of the energy levels descriptor for the two structures.
    \param[in] trSimDist The value of the trace sigma descriptor for the two structures.
    \param[in] rotFunDist The value of the roation function descriptor for the two structures.
 */
void ProSHADE_internal_tasks::ReportDistancesResults ( ProSHADE_settings* settings, std::string str1, std::string str2, proshade_double enLevDist, proshade_double trSigmDist, proshade_double rotFunDist )
{
    std::stringstream hlpSS;
    hlpSS << "Distances between " << str1 << " and " << str2;
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, hlpSS.str() );
    
    std::stringstream hlpSSE;
    hlpSSE << "Energy levels distance    : " << enLevDist;
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, hlpSSE.str() );
    
    std::stringstream hlpSSS;
    hlpSSS << "Trace sigma distance      : " << trSigmDist;
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, hlpSSS.str() );
    
    std::stringstream hlpSSR;
    hlpSSR << "Rotation function distance: " << rotFunDist;
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, hlpSSR.str() );
    
    //================================================ Done
    return ;
    
}

/*! \brief The distances computation settings checks.
 
    This function is called to check the settings object for having all the required information for
    the distances computation task to proceed.
 
    \param[in] settings ProSHADE_settings object specifying the details of how distances computation should be done.
 */
void ProSHADE_internal_tasks::checkDistancesSettings ( ProSHADE_settings* settings )
{
    //================================================ Are there at least two structures?
    if ( settings->inputFiles.size () < 2 )
    {
        throw ProSHADE_exception ( "There are not enough structures for distance computation.", "ED00012", __FILE__, __LINE__, __func__, "There needs to be at least two structures between which\n                    : distances are computed. The ProSHADE_settings object\n                    : contains less than two structures and therefore cannot\n                    : proceed. Please supply at least two structures by\n                    : repeatedly using the addStructure() function." );
    }
    
    //================================================ Is there resolution value set?
    if ( settings->requestedResolution == -1.0 )
    {
        throw ProSHADE_exception ( "Resolution value not set.", "ED00013", __FILE__, __LINE__, __func__, "The resolution value was not set. Please set the\n                    : resolution value for the distance computation by using\n                    : the setResolution() function." );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief The symmetry detection task driver function.
 
    This function is called to run the detect symmetries task according to the information placed in
    the settings object passed as the first argument.
 
    \param[in] settings ProSHADE_settings object specifying the details of how distances computation should be done.
    \param[in] axes A pointer to a vector to which all the axes of the recommended symmetry (if any) will be saved.
    \param[in] allCs A pointer to a vector to which all the detected cyclic symmetries will be saved into.
 */
void ProSHADE_internal_tasks::SymmetryDetectionTask ( ProSHADE_settings* settings, std::vector< proshade_double* >* axes, std::vector < std::vector< proshade_double > >* allCs )
{
    //================================================ Check the settings are complete and meaningful
    checkSymmetrySettings                             ( settings );
    
    //================================================ Now, for each other structure
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( settings->inputFiles.size() ); iter++ )
    {
        //============================================ Create a data object
        ProSHADE_internal_data::ProSHADE_data* symmetryStructure = new ProSHADE_internal_data::ProSHADE_data ( settings );
        
        //============================================ Read in the compared structure
        symmetryStructure->readInStructure            ( settings->inputFiles.at(iter), iter, settings );
        
        //============================================ Internal data processing  (COM, norm, mask, extra space)
        symmetryStructure->processInternalMap         ( settings );
        
        //============================================ Map to sphere
        symmetryStructure->mapToSpheres               ( settings );
        
        //============================================ Get spherical harmonics
        symmetryStructure->computeSphericalHarmonics  ( settings );
        
        //============================================ Compute auto-rotation map
        symmetryStructure->computeRotationFunction    ( settings );
        
        if ( settings->usePeakSearchInRotationFunctionSpace )
        {            
            //======================================== Detect point groups in the angle-axis space
            symmetryStructure->detectSymmetryFromAngleAxisSpace ( settings, axes, allCs );
        }
        else
        {
            //======================================== Detect symmetry using the peak detection in rotation function space
            symmetryStructure->detectSymmetryInStructure ( settings, axes, allCs );
        }
        
        //============================================ Report results
        symmetryStructure->reportSymmetryResults      ( settings );
        
        //============================================ Release memory
        delete symmetryStructure;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief The symmetry computation settings checks.
 
    This function is called to check the settings object for having all the required information for
    the symmetry computation task to proceed.
 
    \param[in] settings ProSHADE_settings object specifying the details of how symmetry detection should be done.
 */
void ProSHADE_internal_tasks::checkSymmetrySettings ( ProSHADE_settings* settings )
{
    //================================================ Are the any structures?
    if ( settings->inputFiles.size () < 1 )
    {
        throw ProSHADE_exception ( "There are not enough structures for symmetry detection.", "ES00028", __FILE__, __LINE__, __func__, "There needs to be at least one structure for which\n                    : symmetry is to be detected. Please supply at least one\n                    : structure by using the addStructure() function." );
    }
    
    //================================================ Is the axis tolerance set properly?
    if ( settings->axisErrTolerance < 0.0 )
    {
        throw ProSHADE_exception ( "Symmetry axis detection tolerance set to negative value.", "ES00053", __FILE__, __LINE__, __func__, "The symmetry axis detection tolerance was manually set to\n                    : negative value. This makes no sense, please supply\n                    : value >= 0.0." );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief The symmetry detection task driver function.
 
    This function is called to run the detect symmetries task according to the information placed in
    the settings object passed as the first argument.
 
    \param[in] settings ProSHADE_settings object specifying the details of how distances computation should be done.
    \param[in] rotationCentre Pointer to vector for saving the position of the centre of rotation about which the rotation is to be done.
    \param[in] eulerAngles Pointer to vector where the three Euler angles will be saved into.
    \param[in] finalTranslation Pointer to a vector where the translation required to move structure from origin to optimal overlay with static structure will be saved into.
 */
void ProSHADE_internal_tasks::MapOverlayTask ( ProSHADE_settings* settings, std::vector < proshade_double >* rotationCentre, std::vector < proshade_double >* eulerAngles, std::vector < proshade_double >* finalTranslation )
{
    //================================================ Check the settings are complete and meaningful
    checkOverlaySettings                              ( settings );
    
    //================================================ Initialise variables
    proshade_double eulA, eulB, eulG, trsX, trsY, trsZ;
    
    //================================================ Create the data objects initially (this time without phase)
    ProSHADE_internal_data::ProSHADE_data* staticStructure = new ProSHADE_internal_data::ProSHADE_data ( settings );
    ProSHADE_internal_data::ProSHADE_data* movingStructure = new ProSHADE_internal_data::ProSHADE_data ( settings );
    
    //================================================ First, run without phase and find best rotation angles
    settings->usePhase                                = false;
    ProSHADE_internal_overlay::getOptimalRotation     ( settings, staticStructure, movingStructure, &eulA, &eulB, &eulG );

    //================================================ Release memory
    delete staticStructure;
    delete movingStructure;
    
    //================================================ Create the data objects again (this time with phase)
    staticStructure                                   = new ProSHADE_internal_data::ProSHADE_data ( settings );
    movingStructure                                   = new ProSHADE_internal_data::ProSHADE_data ( settings );

    //================================================ Now, run with phase and find optimal translation
    settings->usePhase                                = true;
    settings->changeMapResolution                     = true;
    ProSHADE_internal_overlay::getOptimalTranslation  ( settings, staticStructure, movingStructure, &trsX, &trsY, &trsZ, eulA, eulB, eulG );
    
    //================================================ Compute the proper translations using the translation function output
    ProSHADE_internal_misc::addToDoubleVector         ( rotationCentre, movingStructure->originalPdbRotCenX );
    ProSHADE_internal_misc::addToDoubleVector         ( rotationCentre, movingStructure->originalPdbRotCenY );
    ProSHADE_internal_misc::addToDoubleVector         ( rotationCentre, movingStructure->originalPdbRotCenZ );
    ProSHADE_internal_misc::addToDoubleVector         ( finalTranslation, movingStructure->originalPdbTransX );
    ProSHADE_internal_misc::addToDoubleVector         ( finalTranslation, movingStructure->originalPdbTransY );
    ProSHADE_internal_misc::addToDoubleVector         ( finalTranslation, movingStructure->originalPdbTransZ );
    
    //================================================ Write out everything
    movingStructure->writeOutOverlayFiles              ( settings, eulA, eulB, eulG, rotationCentre, finalTranslation );
    
    //================================================ Save the rotation and rest of translations
    ProSHADE_internal_misc::addToDoubleVector         ( eulerAngles, eulA );
    ProSHADE_internal_misc::addToDoubleVector         ( eulerAngles, eulB );
    ProSHADE_internal_misc::addToDoubleVector         ( eulerAngles, eulG );
    
    //================================================ Report results to user
    movingStructure->reportOverlayResults             ( settings, rotationCentre, eulerAngles, finalTranslation );
    
    //================================================ Release memory
    delete staticStructure;
    delete movingStructure;
    
    //================================================ Done
    return ;
    
}

/*! \brief The map overlay computation settings checks.
 
    This function is called to check the settings object for having all the required information for
    the map overlay task to proceed.
 
    \param[in] settings ProSHADE_settings object specifying the details of how map overlay should be done.
 */
void ProSHADE_internal_tasks::checkOverlaySettings ( ProSHADE_settings* settings )
{
    //================================================ Are the any structures?
    if ( settings->inputFiles.size () != 2 )
    {
        throw ProSHADE_exception ( "There are not enough structures for map overlay\n                    : computation.", "EO00033", __FILE__, __LINE__, __func__, "There needs to be exactly two structures for map overlay\n                    : mode to work; the first structure is the static and the\n                    : second is the moving structure." );
    }
    
    //================================================ Done
    return ;
    
}
