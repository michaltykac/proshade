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
    \version   0.7.6.1
    \date      AUG 2021
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
        ProSHADE_internal_data::ProSHADE_data* strToRebox = new ProSHADE_internal_data::ProSHADE_data ( );
        
        //============================================ Read in the file
        strToRebox->readInStructure                   ( settings->inputFiles.at(iter), iter, settings );
        
        //============================================ Save the original boundaries
        ProSHADE_internal_misc::deepCopyBoundsSigPtrVector ( originalBounds, strToRebox->getXFromPtr(), strToRebox->getXToPtr(), strToRebox->getYFromPtr(), strToRebox->getYToPtr(), strToRebox->getZFromPtr(), strToRebox->getZToPtr() );
        
        //============================================ Internal data processing  (COM, norm, mask, extra space)
        strToRebox->processInternalMap                ( settings );
        
        //============================================ Create new structure for re-boxing
        ProSHADE_internal_data::ProSHADE_data* reBoxStr = new ProSHADE_internal_data::ProSHADE_data ( );
        
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
        proshade_double* mapCopy                      = nullptr;
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
        if ( settings->requestedResolution == 0.0f )
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
    ProSHADE_internal_data::ProSHADE_data* compareAgainst  = new ProSHADE_internal_data::ProSHADE_data ( );
    
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
        ProSHADE_internal_data::ProSHADE_data* compareChanging = new ProSHADE_internal_data::ProSHADE_data ( );

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
    const FloatingPoint< proshade_single > lhs ( settings->requestedResolution ), rhs ( -1.0f );
    if ( lhs.AlmostEquals ( rhs ) )
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
void ProSHADE_internal_tasks::SymmetryDetectionTask ( ProSHADE_settings* settings, std::vector< proshade_double* >* axes, std::vector < std::vector< proshade_double > >* allCs, std::vector< proshade_double >* mapCOMShift )
{
    //================================================ Check the settings are complete and meaningful
    checkSymmetrySettings                             ( settings );
    
    //================================================ Now, for each other structure
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( settings->inputFiles.size() ); iter++ )
    {
        //============================================ Create a data object
        ProSHADE_internal_data::ProSHADE_data* symmetryStructure = new ProSHADE_internal_data::ProSHADE_data ( );
        
        //============================================ Read in the compared structure
        symmetryStructure->readInStructure            ( settings->inputFiles.at(iter), iter, settings );
        
        if ( settings->findSymCentre )
        {
            //======================================== ...
            std::cout << "@@@ Attempting to find the symmetry centre using phase-less detection." << std::endl;
            
            //======================================== Detect symmetry without using the phase information
            settings->usePhase                        = false;
            std::string origReqSymType                = settings->requestedSymmetryType;
            settings->requestedSymmetryType           = "onlyC";
            settings->moveToCOM                       = false;
            symmetryStructure->processInternalMap               ( settings );
            symmetryStructure->mapToSpheres                     ( settings );
            symmetryStructure->computeSphericalHarmonics        ( settings );
            symmetryStructure->computeRotationFunction          ( settings );
            symmetryStructure->detectSymmetryFromAngleAxisSpace ( settings, axes, allCs );
            
            //== Find reliable phase-less Cs
            std::vector< proshade_unsign > reliableCs;
            for ( proshade_unsign axIt = 0; axIt < static_cast< proshade_unsign > ( allCs->size() ); axIt++ ) { std::cout << allCs->at(axIt)[0] << " | " << allCs->at(axIt)[1] << " x " << allCs->at(axIt)[2] << " x " << allCs->at(axIt)[3] << " | " << allCs->at(axIt)[5] << " || " << allCs->at(axIt)[6] << std::endl; }
            for ( proshade_unsign axIt = 0; axIt < static_cast< proshade_unsign > ( allCs->size() ); axIt++ ) { if ( ( allCs->at(axIt)[5] > 0.85 ) && ( allCs->at(axIt)[6] > 0.85 ) ) { ProSHADE_internal_misc::addToUnsignVector ( &reliableCs, axIt ); } }
            
            //== If none found, report
            if ( reliableCs.size() == 0 )
            {
                std::cout << "!!! Pruser - zadna osa !!!" << std::endl;
            }
            
            //== Are any orthogonal?
            proshade_double dotProduct, maxOrtSum = 0.0, curOrtSum = 0.0;
            proshade_unsign maxOrtAx1 = 0, maxOrtAx2 = 0;
            for ( size_t relAx1 = 0; relAx1 < reliableCs.size(); relAx1++ )
            {
                for ( size_t relAx2 = 1; relAx2 < reliableCs.size(); relAx2++ )
                {
                    //================================ Ignore same axes
                    if ( relAx1 >= relAx2 ) { continue; }
                    
                    //================================ Are the two axes orthogonal?
                    dotProduct                        = ProSHADE_internal_maths::computeDotProduct ( &allCs->at(reliableCs.at(relAx1))[1], &allCs->at(reliableCs.at(relAx1))[2],
                                                                                                     &allCs->at(reliableCs.at(relAx1))[3], &allCs->at(reliableCs.at(relAx2))[1],
                                                                                                     &allCs->at(reliableCs.at(relAx2))[2], &allCs->at(reliableCs.at(relAx2))[3] );
                    
                    //================================ If close to zero, these two axes are perpendicular
                    if ( std::abs( dotProduct ) < settings->axisErrTolerance )
                    {
                        //============================ Find sum
                        curOrtSum                     = allCs->at(reliableCs.at(relAx1))[5] + allCs->at(reliableCs.at(relAx1))[6] + allCs->at(reliableCs.at(relAx2))[5] + allCs->at(reliableCs.at(relAx2))[6];
                        
                        //============================ If best, save it
                        if ( curOrtSum > maxOrtSum )
                        {
                            maxOrtSum                     = curOrtSum;
                            maxOrtAx1                     = static_cast< proshade_unsign > ( relAx1 );
                            maxOrtAx2                     = static_cast< proshade_unsign > ( relAx2 );
                        }
                    }
                }
            }
            
            //== If no orthogonal axes, deal with this case (i.e. single axis centre detection)
            if ( maxOrtAx2 == 0 )
            {
                //== ...
            }
            
            //== Optimise the orthogonal pair
            std::vector < std::vector< proshade_double > > ortPair;
            std::vector< proshade_double > hlpVec;
            hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx1))[0] ); hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx1))[1] ); hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx1))[2] );
            hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx1))[3] ); hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx1))[4] ); hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx1))[5] );
            hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx1))[6] );
            ortPair.push_back ( hlpVec ); hlpVec.clear ( );
            hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx2))[0] ); hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx2))[1] ); hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx2))[2] );
            hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx2))[3] ); hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx2))[4] ); hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx2))[5] );
            hlpVec.push_back ( allCs->at(reliableCs.at(maxOrtAx2))[6] );
            ortPair.push_back ( hlpVec );
            
            ProSHADE_internal_symmetry::optimiseDGroupAngleFromAxesHeights ( &ortPair, symmetryStructure, settings );
            
            //== Assuming orthogonal set was found
            std::cout << "The best orthogonal set is: " << ortPair.at(0).at(0) << " | " << ortPair.at(0).at(1) << " x " << ortPair.at(0).at(2) << " x " << ortPair.at(0).at(3) << " | " << ortPair.at(0).at(5) << " | " << ortPair.at(0).at(6) << std::endl;
            std::cout << " and                      : " << ortPair.at(1).at(0) << " | " << ortPair.at(1).at(1) << " x " << ortPair.at(1).at(2) << " x " << ortPair.at(1).at(3) << " | " << ortPair.at(1).at(5) << " | " << ortPair.at(1).at(6) << std::endl;

            //== Generate the symmetry elements for the detected axes
            std::vector< proshade_unsign > axList;
            ProSHADE_internal_misc::addToUnsignVector ( &axList, maxOrtAx1 );
            ProSHADE_internal_misc::addToUnsignVector ( &axList, maxOrtAx2 );
            std::vector< std::vector < proshade_double > > symElems = symmetryStructure->getAllGroupElements ( settings, axList, "D", settings->axisErrTolerance + 0.1 );
            
            std::cout << "Found total of " << symElems.size() << " elements, the first being:" << std::endl;
            std::cout << symElems.at(0).at(0) << " | " << symElems.at(0).at(1) << " | " << symElems.at(0).at(2) << std::endl;
            std::cout << symElems.at(0).at(3) << " | " << symElems.at(0).at(4) << " | " << symElems.at(0).at(5) << std::endl;
            std::cout << symElems.at(0).at(6) << " | " << symElems.at(0).at(7) << " | " << symElems.at(0).at(8) << std::endl;
            
            //== Re-read the map, this time with phases
            delete symmetryStructure;
            settings->usePhase                        = true;
            settings->moveToCOM                       = false;
            settings->requestedSymmetryType           = origReqSymType;
            settings->addExtraSpace                   = 0.0;
            symmetryStructure                         = new ProSHADE_internal_data::ProSHADE_data ( );
            symmetryStructure->readInStructure        ( settings->inputFiles.at(iter), iter, settings );
            symmetryStructure->processInternalMap     ( settings );
            
            //== Convert original map to Fourier space
            fftw_complex *origMap                     = new fftw_complex [symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim()];
            fftw_complex *origCoeffs                  = new fftw_complex [symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim()];
            fftw_complex *rotMapComplex               = new fftw_complex [symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim()];
            fftw_complex *rotCoeffs                   = new fftw_complex [symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim()];
            fftw_complex *trFunc                      = new fftw_complex [symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim()];
            fftw_complex *trFuncCoeffs                = new fftw_complex [symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim()];
            ProSHADE_internal_misc::checkMemoryAllocation ( origMap,       __FILE__, __LINE__, __func__ );
            ProSHADE_internal_misc::checkMemoryAllocation ( origCoeffs,    __FILE__, __LINE__, __func__ );
            ProSHADE_internal_misc::checkMemoryAllocation ( rotMapComplex, __FILE__, __LINE__, __func__ );
            ProSHADE_internal_misc::checkMemoryAllocation ( rotCoeffs,     __FILE__, __LINE__, __func__ );
            ProSHADE_internal_misc::checkMemoryAllocation ( trFunc,        __FILE__, __LINE__, __func__ );
            ProSHADE_internal_misc::checkMemoryAllocation ( trFuncCoeffs,  __FILE__, __LINE__, __func__ );
            for ( size_t it = 0; it < static_cast< size_t > ( symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim() ); it++ ) { origMap[it][0] = symmetryStructure->getMapValue( it ); origMap[it][1] = 0.0; }
            fftw_plan planForwardFourier              = fftw_plan_dft_3d ( static_cast< int > ( symmetryStructure->getXDim() ), static_cast< int > ( symmetryStructure->getYDim() ), static_cast< int > ( symmetryStructure->getZDim() ), origMap, origCoeffs, FFTW_FORWARD,        FFTW_ESTIMATE );
            fftw_plan planForwardFourierRot           = fftw_plan_dft_3d ( static_cast< int > ( symmetryStructure->getXDim() ), static_cast< int > ( symmetryStructure->getYDim() ), static_cast< int > ( symmetryStructure->getZDim() ), rotMapComplex, rotCoeffs, FFTW_FORWARD,   FFTW_ESTIMATE );
            fftw_plan planReverseFourierComb          = fftw_plan_dft_3d ( static_cast< int > ( symmetryStructure->getXDim() ), static_cast< int > ( symmetryStructure->getYDim() ), static_cast< int > ( symmetryStructure->getZDim() ), trFuncCoeffs, trFunc, FFTW_BACKWARD,  FFTW_ESTIMATE );
            
            fftw_execute                              ( planForwardFourier );
            
            
            //== For each group element except the identity one
            proshade_double avgRX = 0.0, avgRY = 0.0, avgRZ = 0.0;
            proshade_double axX, axY, axZ, axAng;
            proshade_double *rMat                     = new proshade_double[9];
            ProSHADE_internal_misc::checkMemoryAllocation ( rMat, __FILE__, __LINE__, __func__ );
            for ( size_t grEl = 1; grEl < symElems.size(); grEl++ )
            {
                std::cout << std::endl;
                std::cout << symElems.at(grEl).at(0) << ", " << symElems.at(grEl).at(1) << ", " << symElems.at(grEl).at(2) << std::endl;
                std::cout << symElems.at(grEl).at(3) << ", " << symElems.at(grEl).at(4) << ", " << symElems.at(grEl).at(5) << std::endl;
                std::cout << symElems.at(grEl).at(6) << ", " << symElems.at(grEl).at(7) << ", " << symElems.at(grEl).at(8) << std::endl;
                std::cout << std::endl;
                
                //== Rotate the map by the rotation matrix
                proshade_double *rotMap;
                for ( size_t mIt = 0; mIt < 9; mIt++ ) { rMat[mIt] = symElems.at(grEl).at(mIt); }
                ProSHADE_internal_maths::getAxisAngleFromRotationMatrix ( rMat, &axX, &axY, &axZ, &axAng );
                symmetryStructure->rotateMapRealSpace ( axX, axY, axZ, axAng, rotMap );
                
                std::stringstream hh2;
                hh2 << "rotMap" << grEl << ".map";
                proshade_double* hlpMap2 = new proshade_double[symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim()];
                for ( int i = 0; i < symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim(); i++ ) { hlpMap2[i] = symmetryStructure->getInternalMap()[i]; }
                for ( int i = 0; i < symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim(); i++ ) { symmetryStructure->getInternalMap()[i] = rotMap[i]; }
                symmetryStructure->writeMap ( hh2.str() );
                for ( int i = 0; i < symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim(); i++ ) { symmetryStructure->getInternalMap()[i] = hlpMap2[i]; }
                
                //== Convert to Fourier space
                for ( size_t it = 0; it < static_cast< size_t > ( symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim() ); it++ ) { rotMapComplex[it][0] = rotMap[it]; rotMapComplex[it][1] = 0.0; }
                fftw_execute                          ( planForwardFourierRot );
                
                //== Combine coeffs for translation function
                ProSHADE_internal_overlay::combineFourierForTranslation ( origCoeffs, rotCoeffs, trFuncCoeffs, symmetryStructure->getXDim(), symmetryStructure->getYDim(), symmetryStructure->getZDim() );
                
                //== Compute translation function
                fftw_execute                          ( planReverseFourierComb );
                
                //== Find peak
                proshade_double mapPeak = 0.0, trsX, trsY, trsZ;
                ProSHADE_internal_overlay::findHighestValueInMap  ( trFunc, symmetryStructure->getXDim(), symmetryStructure->getYDim(), symmetryStructure->getZDim(), &trsX, &trsY, &trsZ, &mapPeak );
                
                //== Not over half
                if ( trsX > ( static_cast< proshade_double > ( symmetryStructure->getXDim() ) / 2.0 ) ) { trsX = trsX - static_cast< proshade_double > ( symmetryStructure->getXDim() ); }
                if ( trsY > ( static_cast< proshade_double > ( symmetryStructure->getYDim() ) / 2.0 ) ) { trsY = trsY - static_cast< proshade_double > ( symmetryStructure->getYDim() ); }
                if ( trsZ > ( static_cast< proshade_double > ( symmetryStructure->getZDim() ) / 2.0 ) ) { trsZ = trsZ - static_cast< proshade_double > ( symmetryStructure->getZDim() ); }
                
                trsX *= ( symmetryStructure->getXDimSize() / symmetryStructure->getXDim() );
                trsY *= ( symmetryStructure->getYDimSize() / symmetryStructure->getYDim() );
                trsZ *= ( symmetryStructure->getZDimSize() / symmetryStructure->getZDim() );
                std::cout << " ### Found best translation to be " << trsX << ", " << trsY << ", " << trsZ << " with peak " << mapPeak << std::endl;
                
                //== Optimise the translation
                
                //== Find Moore-Penrose pseudo inverse of I-Ri
                int dim                                           = 3;
                char job                                          = 'A';                                   // Save computation of parts of U and V matrices, they are not needed here
                double* singularValues                            = new double[dim];                       // The array of singular values
                double *rotMatU                                   = new double [dim*dim];                  // The U matrix space
                double *rotMatV                                   = new double [dim*dim];                  // The V^T matrix space
                double *work                                      = new double [static_cast< proshade_unsign >( ( 3 * dim ) + pow( dim, 2 ) * dim)]; // Workspace, minimum required is 4*dim^2 + 7*dim, using more for performance
                int workDim                                       = static_cast< int > ( 2 * ( ( 4 * dim * dim ) + ( 7 * dim ) ) ); // Formalism stating just that
                double* rwork                                     = new double[static_cast<proshade_unsign>((5 * dim) + 5 * pow(dim,2))]; // Required by LAPACK
                int* iwork                                        = new int[(8 * dim)];                    // Required by LAPACK
                int returnValue                                   = 0;                                     // This will tell if operation succeeded
                ProSHADE_internal_misc::checkMemoryAllocation     ( singularValues, __FILE__, __LINE__, __func__ );
                ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatU,        __FILE__, __LINE__, __func__ );
                ProSHADE_internal_misc::checkMemoryAllocation     ( rotMatV,        __FILE__, __LINE__, __func__ );
                ProSHADE_internal_misc::checkMemoryAllocation     ( work,           __FILE__, __LINE__, __func__ );
                ProSHADE_internal_misc::checkMemoryAllocation     ( rwork,          __FILE__, __LINE__, __func__ );
                ProSHADE_internal_misc::checkMemoryAllocation     ( iwork,          __FILE__, __LINE__, __func__ );
                
                //================================================ Load input data into array in column-major order
                double *matrixToDecompose                         = new double[dim*dim];
                ProSHADE_internal_misc::checkMemoryAllocation     ( matrixToDecompose, __FILE__, __LINE__, __func__ );
                for ( int rowIt = 0; rowIt < dim; rowIt++ )
                {
                    for ( int colIt = 0; colIt < dim; colIt++ )
                    {
                        if ( rowIt == colIt ) { matrixToDecompose[(colIt*dim)+rowIt] = 1.0 - rMat[(rowIt*dim)+colIt]; }
                        else                  { matrixToDecompose[(colIt*dim)+rowIt] = 0.0 - rMat[(rowIt*dim)+colIt]; }
                    }
                }
                
                //================================================ Run LAPACK ZGESDD
                dgesdd_                                           ( &job, &dim, &dim, matrixToDecompose, &dim, singularValues, rotMatU, &dim, rotMatV, &dim,
                                                                    work, &workDim, rwork, iwork, &returnValue );
                if ( returnValue != 0 )
                {
                    ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! SVD algorithm did not converge.", "WS00069" );
                }
                
                //== Determine positivity
                bool anyPositive = false;
                std::vector< bool > positivityTest;
                for ( proshade_unsign it = 0; it < static_cast< proshade_unsign > ( dim ); it++ )
                {
                    positivityTest.push_back ( singularValues[it] > 0.001 );
                    if ( positivityTest.at(it) ) { anyPositive = true; }
                }
                
                proshade_double* pseudoInverseMat;
                if ( anyPositive )
                {
                    //== Set all non-positive
                    if ( !positivityTest.at(0) )
                    {
                        singularValues[0]             = 0.0;
                        rotMatU[0]                    = 0.0;
                        rotMatU[1]                    = 0.0;
                        rotMatU[2]                    = 0.0;
                        rotMatV[0]                    = 0.0;
                        rotMatV[3]                    = 0.0;
                        rotMatV[6]                    = 0.0;
                    }
                    else { singularValues[0] = 1.0 / singularValues[0]; }

                    if ( !positivityTest.at(1) )
                    {
                        singularValues[1]             = 0.0;
                        rotMatU[3]                    = 0.0;
                        rotMatU[4]                    = 0.0;
                        rotMatU[5]                    = 0.0;
                        rotMatV[1]                    = 0.0;
                        rotMatV[4]                    = 0.0;
                        rotMatV[7]                    = 0.0;
                    }
                    else { singularValues[1] = 1.0 / singularValues[1]; }

                    //== The last singular value (they are in order) must be zero as Ri is a rotation matrix with at least one eigenvalue 1 and therefore I - Ri must have at least one eigenvalue 0.
                    singularValues[2]                 = 0.0;
                    rotMatU[6]                        = 0.0;
                    rotMatU[7]                        = 0.0;
                    rotMatU[8]                        = 0.0;
                    rotMatV[2]                        = 0.0;
                    rotMatV[5]                        = 0.0;
                    rotMatV[8]                        = 0.0;
                    
                    //== All positive values formula
                    proshade_double* diagMat          = ProSHADE_internal_maths::build3x3MatrixFromDiag ( singularValues );
                    proshade_double* hlpMat           = ProSHADE_internal_maths::compute3x3MatrixMultiplication ( diagMat, rotMatU );
                    pseudoInverseMat                  = ProSHADE_internal_maths::compute3x3MatrixMultiplication ( rotMatV, hlpMat );
                }
                else
                {
                    //== No axis in matrix
                    pseudoInverseMat                         = new proshade_double[9];
                    ProSHADE_internal_misc::checkMemoryAllocation     ( pseudoInverseMat, __FILE__, __LINE__, __func__ );
                    
                    for ( size_t mIt = 0; mIt < 9; mIt++ ) { pseudoInverseMat[mIt] = 0.0; }
                }
                
                //================================================ Free memory
                delete[] work;
                delete[] rwork;
                delete[] iwork;
                delete[] matrixToDecompose;
                delete[] singularValues;
                delete[] rotMatU;
                delete[] rotMatV;
                
                
                
                std::stringstream hh;
                hh << "trsRotMap" << grEl << ".map";
                proshade_double* hlpMap = new proshade_double[symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim()];
                for ( int i = 0; i < symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim(); i++ ) { hlpMap[i] = symmetryStructure->getInternalMap()[i]; }
                for ( int i = 0; i < symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim(); i++ ) { symmetryStructure->getInternalMap()[i] = rotMap[i]; }
                
                ProSHADE_internal_mapManip::moveMapByFourier ( symmetryStructure->getInternalMap(), trsX, trsY, trsZ,
                                                               symmetryStructure->getXDimSize(), symmetryStructure->getYDimSize(), symmetryStructure->getZDimSize(),
                                                               static_cast< proshade_signed > ( symmetryStructure->getXDim() ), static_cast< proshade_signed > ( symmetryStructure->getYDim() ),
                                                               static_cast< proshade_signed > ( symmetryStructure->getZDim() ) );
                
                symmetryStructure->writeMap ( hh.str() );
                for ( int i = 0; i < symmetryStructure->getXDim() * symmetryStructure->getYDim() * symmetryStructure->getZDim(); i++ ) { symmetryStructure->getInternalMap()[i] = hlpMap[i]; }
                
                
                //== Multiple translation with the inverted I- Ri to get centre of rotation
                proshade_double* rotCen               = ProSHADE_internal_maths::compute3x3MatrixVectorMultiplication ( pseudoInverseMat, trsX, trsY, trsZ );
                
                
                delete[] pseudoInverseMat;
                
                //== Sum the translations
                avgRX += rotCen[0];
                avgRY += rotCen[1];
                avgRZ += rotCen[2];
                
                std::cout << " ### Centre of rotation according to this group element is: " << rotCen[0] << " x " << rotCen[1] << " x " << rotCen[2] << std::endl;
                
                delete[] rotMap;
            }
            
            //== Average the translation sum
            avgRX /= static_cast< proshade_double > ( symElems.size() );
            avgRY /= static_cast< proshade_double > ( symElems.size() );
            avgRZ /= static_cast< proshade_double > ( symElems.size() );
            
            std::cout << "### Averaged rotation centre position is: " << avgRX << " x " << avgRY << " x " << avgRZ << std::endl;
            
            delete[] rMat;
            fftw_destroy_plan                         ( planReverseFourierComb );
            fftw_destroy_plan                         ( planForwardFourier );
            fftw_destroy_plan                         ( planForwardFourierRot );
            delete[] origMap;
            delete[] origCoeffs;
            delete[] rotMapComplex;
            delete[] rotCoeffs;
            delete[] trFunc;
            delete[] trFuncCoeffs;
            
            //== Translate
//            symmetryStructure->writeMap ( "mapNotCentred.map" );
//            ProSHADE_internal_mapManip::moveMapByFourier ( symmetryStructure->getInternalMap(), avgRX, avgRY, avgRZ,
//                                                           symmetryStructure->getXDimSize(), symmetryStructure->getYDimSize(), symmetryStructure->getZDimSize(),
//                                                           static_cast< proshade_signed > ( symmetryStructure->getXDim() ), static_cast< proshade_signed > ( symmetryStructure->getYDim() ),
//                                                           static_cast< proshade_signed > ( symmetryStructure->getZDim() ) );
//            symmetryStructure->writeMap ( "mapCentred.map" );
            
            exit(0);
        }
        
        //============================================ Internal data processing  (COM, norm, mask, extra space)
        symmetryStructure->processInternalMap         ( settings );
        
        //============================================ Map to sphere
        symmetryStructure->mapToSpheres               ( settings );
        
        //============================================ Get spherical harmonics
        symmetryStructure->computeSphericalHarmonics  ( settings );
        
        //============================================ Compute auto-rotation map
        symmetryStructure->computeRotationFunction    ( settings );
        
        //======================================== Detect point groups in the angle-axis space
        symmetryStructure->detectSymmetryFromAngleAxisSpace ( settings, axes, allCs );
        
        //============================================ Report results
        symmetryStructure->reportSymmetryResults      ( settings );
        
        //============================================ Save internal map shift to run object,
        ProSHADE_internal_misc::addToDoubleVector     ( mapCOMShift, symmetryStructure->mapCOMProcessChangeX );
        ProSHADE_internal_misc::addToDoubleVector     ( mapCOMShift, symmetryStructure->mapCOMProcessChangeY );
        ProSHADE_internal_misc::addToDoubleVector     ( mapCOMShift, symmetryStructure->mapCOMProcessChangeZ );
        
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
    ProSHADE_internal_data::ProSHADE_data* staticStructure = new ProSHADE_internal_data::ProSHADE_data ( );
    ProSHADE_internal_data::ProSHADE_data* movingStructure = new ProSHADE_internal_data::ProSHADE_data ( );

    //================================================ First, run without phase and find best rotation angles
    settings->usePhase                                = false;
    ProSHADE_internal_overlay::getOptimalRotation     ( settings, staticStructure, movingStructure, &eulA, &eulB, &eulG );

    //================================================ Release memory
    delete staticStructure;
    delete movingStructure;
    
    //================================================ Create the data objects again (this time with phase)
    staticStructure                                   = new ProSHADE_internal_data::ProSHADE_data ( );
    movingStructure                                   = new ProSHADE_internal_data::ProSHADE_data ( );

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
    
    //================================================ If centring is on, turn it off and report warning.
    if ( settings->moveToCOM )
    {
        ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Map centring was requested, but makes no sense for overlay mode. Turning it off.", "WO00066" );
        settings->moveToCOM                           = false;
    }
    
    //================================================ Done
    return ;
    
}
