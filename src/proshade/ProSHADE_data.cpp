/*! \file ProSHADE_data.cpp
    \brief This is the source file containing internal data representation and manipulation structures and functions.
 
    This source file contains the ProSHADE_data class definitions as well as the code for simple manipulations with the data
    (more complex manipulations are done in dedicated source files) and caller functions for the more complex manipulations.
    The class described here is how ProSHADE stores the structural data internally; however, the user should not need to access
    any of this code manually, as changes to this structure may have large consequences unforseen by the user.
   
    Copyright by Michal Tykac and individual contributors. All rights reserved.
     
    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
     
    This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In     no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data     or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility     of such damage.
     
    \author    Michal Tykac
    \author    Garib N. Murshudov
    \version   0.7.5.0
    \date      DEC 2020
 */

//==================================================== ProSHADE
#include "ProSHADE_data.hpp"

//==================================================== Gemmi PDB output - this cannot be with the rest of includes for some stb_sprintf library related reasons ...
#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_pdb.hpp>

/*! \brief Constructor for getting empty ProSHADE_data class.
 
    This constructor creates an empty data structure which can later be filled with data and used to process
    these data further.
 
    \param[in] settings ProSHADE_settings object specifying what should be done.
    \param[out] X Empty data object with deault values.
 */
ProSHADE_internal_data::ProSHADE_data::ProSHADE_data ( ProSHADE_settings* settings )
{
    //================================================ Initialise variables
    // ... Variables regarding input file
    this->fileName                                    = "";
    this->fileType                                    = ProSHADE_internal_io::UNKNOWN;
    
    // ... Variables regarding map
    this->internalMap                                 = NULL;
    
    // ... Variables regarding map information
    this->xDimSize                                    = 0.0;
    this->yDimSize                                    = 0.0;
    this->zDimSize                                    = 0.0;
    this->aAngle                                      = 0.0;
    this->bAngle                                      = 0.0;
    this->cAngle                                      = 0.0;
    this->xDimIndices                                 = 0;
    this->yDimIndices                                 = 0;
    this->zDimIndices                                 = 0;
    this->xGridIndices                                = 0;
    this->yGridIndices                                = 0;
    this->zGridIndices                                = 0;
    this->xAxisOrder                                  = 1;
    this->yAxisOrder                                  = 2;
    this->zAxisOrder                                  = 3;
    this->xAxisOrigin                                 = 0;
    this->yAxisOrigin                                 = 0;
    this->zAxisOrigin                                 = 0;
    this->comMovX                                     = 0.0;
    this->comMovY                                     = 0.0;
    this->comMovZ                                     = 0.0;
    this->xCom                                        = 0.0;
    this->yCom                                        = 0.0;
    this->zCom                                        = 0.0;
    
    // ... Variables regarding original input values (i.e. these do not change with ProSHADE manipulations)
    this->xDimSizeOriginal                            = 0.0;
    this->yDimSizeOriginal                            = 0.0;
    this->zDimSizeOriginal                            = 0.0;
    this->xDimIndicesOriginal                         = 0;
    this->yDimIndicesOriginal                         = 0;
    this->zDimIndicesOriginal                         = 0;
    this->xAxisOriginOriginal                         = 0;
    this->yAxisOriginOriginal                         = 0;
    this->zAxisOriginOriginal                         = 0;
    this->originalMapXCom                             = 0.0;
    this->originalMapYCom                             = 0.0;
    this->originalMapZCom                             = 0.0;
    this->mapPostRotXCom                              = 0.0;
    this->mapPostRotYCom                              = 0.0;
    this->mapPostRotZCom                              = 0.0;
    
    // ... Variables regarding rotation and translation of original input files
    this->originalPdbRotCenX                          = 0.0;
    this->originalPdbRotCenY                          = 0.0;
    this->originalPdbRotCenZ                          = 0.0;
    this->originalPdbTransX                           = 0.0;
    this->originalPdbTransY                           = 0.0;
    this->originalPdbTransZ                           = 0.0;
    
    // ... Variables regarding iterator positions
    this->xFrom                                       = 0;
    this->yFrom                                       = 0;
    this->zFrom                                       = 0;
    this->xTo                                         = 0;
    this->yTo                                         = 0;
    this->zTo                                         = 0;
    
    // ... Variables regarding SH mapping spheres
    this->spherePos                                   = std::vector<proshade_single> ( );
    this->noSpheres                                   = 0;
    this->spheres                                     = NULL;
    this->sphericalHarmonics                          = NULL;
    this->rotSphericalHarmonics                       = NULL;
    this->maxShellBand                                = 0;
    
    // ... Variables regarding shape distance computations
    this->rrpMatrices                                 = NULL;
    this->eMatrices                                   = NULL;
    this->so3Coeffs                                   = NULL;
    this->so3CoeffsInverse                            = NULL;
    this->wignerMatrices                              = NULL;
    this->integrationWeight                           = 0.0;
    this->maxCompBand                                 = 0;
    this->translationMap                              = NULL;
    
    
    // ... Control variables
    this->isEmpty                                     = true;
    
    //================================================ Done
    
}

/*! \brief Constructor for creating ProSHADE_data structure with data.
 
    This constructor creates a data structure with all the map information, so that maps obtained from other software could
    be loeaded and used. This function makes a lot of assumptions (all angles are 90 degrees, axis grids are equal to indices,
    axis order is XYZ and axis origin is the first index in all dimensions). If any of these are not true, the user is required
    to change the appropriate internal values after this function has returned the object.
 
    \param[in] settings ProSHADE_settings object specifying what should be done.
    \param[in] strName The name of the structure for reference.
    \param[in] mapVals A pointer to array where all the map data are.
    \param[in] len The length of this map values array.
    \param[in] xDmSz The size of the x-axis dimension in Angstroms.
    \param[in] yDmSz The size of the y-axis dimension in Angstroms.
    \param[in] zDmSz The size of the z-axis dimension in Angstroms.
    \param[in] xDmInd The size of the x-axis dimension in indices.
    \param[in] yDmInd The size of the y-axis dimension in indices.
    \param[in] zDmInd The size of the z-axis dimension in indices.
    \param[in] xFr The first index statting position along the x-axis.
    \param[in] yFr The first index statting position along the y-axis.
    \param[in] zFr The first index statting position along the z-axis.
    \param[in] xT The last index end position along the x-axis.
    \param[in] yT The last index end position along the y-axis.
    \param[in] zT The last index end position along the z-axis.
    \param[in] inputO The input order for this structure.
    \param[out] X Empty data object with filled in values and map.
 */
ProSHADE_internal_data::ProSHADE_data::ProSHADE_data ( ProSHADE_settings* settings, std::string strName, double *mapVals, int len, proshade_single xDmSz, proshade_single yDmSz, proshade_single zDmSz, proshade_unsign xDmInd, proshade_unsign yDmInd, proshade_unsign zDmInd, proshade_signed xFr, proshade_signed yFr, proshade_signed zFr, proshade_signed xT, proshade_signed yT, proshade_signed zT, proshade_unsign inputO )
{
    //================================================ Initialise variables
    // ... Variables regarding input file
    this->fileName                                    = strName;
    this->fileType                                    = ProSHADE_internal_io::MAP;
    
    // ... Variables regarding map
    this->internalMap                                 = NULL;
    
    // ... Variables regarding map information
    this->xDimSize                                    = xDmSz;
    this->yDimSize                                    = yDmSz;
    this->zDimSize                                    = zDmSz;
    this->aAngle                                      = 90.0;
    this->bAngle                                      = 90.0;
    this->cAngle                                      = 90.0;
    this->xDimIndices                                 = xDmInd;
    this->yDimIndices                                 = yDmInd;
    this->zDimIndices                                 = zDmInd;
    this->xGridIndices                                = xDmInd;
    this->yGridIndices                                = yDmInd;
    this->zGridIndices                                = zDmInd;
    this->xAxisOrder                                  = 1;
    this->yAxisOrder                                  = 2;
    this->zAxisOrder                                  = 3;
    this->xAxisOrigin                                 = xFr;
    this->yAxisOrigin                                 = yFr;
    this->zAxisOrigin                                 = zFr;
    this->comMovX                                     = 0.0;
    this->comMovY                                     = 0.0;
    this->comMovZ                                     = 0.0;
    this->xCom                                        = 0.0;
    this->yCom                                        = 0.0;
    this->zCom                                        = 0.0;
    
    // ... Variables regarding original input values (i.e. these do not change with ProSHADE manipulations)
    this->xDimSizeOriginal                            = 0.0;
    this->yDimSizeOriginal                            = 0.0;
    this->zDimSizeOriginal                            = 0.0;
    this->xDimIndicesOriginal                         = 0;
    this->yDimIndicesOriginal                         = 0;
    this->zDimIndicesOriginal                         = 0;
    this->xAxisOriginOriginal                         = 0;
    this->yAxisOriginOriginal                         = 0;
    this->zAxisOriginOriginal                         = 0;
    this->originalMapXCom                             = 0.0;
    this->originalMapYCom                             = 0.0;
    this->originalMapZCom                             = 0.0;
    this->mapPostRotXCom                              = 0.0;
    this->mapPostRotYCom                              = 0.0;
    this->mapPostRotZCom                              = 0.0;
    
    // ... Variables regarding rotation and translation of original input files
    this->originalPdbRotCenX                          = 0.0;
    this->originalPdbRotCenY                          = 0.0;
    this->originalPdbRotCenZ                          = 0.0;
    this->originalPdbTransX                           = 0.0;
    this->originalPdbTransY                           = 0.0;
    this->originalPdbTransZ                           = 0.0;
    
    // ... Variables regarding iterator positions
    this->xFrom                                       = xFr;
    this->yFrom                                       = yFr;
    this->zFrom                                       = zFr;
    this->xTo                                         = xT;
    this->yTo                                         = yT;
    this->zTo                                         = zT;
    
    // ... Variables regarding SH mapping spheres
    this->spherePos                                   = std::vector<proshade_single> ( );
    this->noSpheres                                   = 0;
    this->spheres                                     = NULL;
    this->sphericalHarmonics                          = NULL;
    this->rotSphericalHarmonics                       = NULL;
    this->maxShellBand                                = 0;
    
    // ... Variables regarding shape distance computations
    this->rrpMatrices                                 = NULL;
    this->eMatrices                                   = NULL;
    this->so3Coeffs                                   = NULL;
    this->so3CoeffsInverse                            = NULL;
    this->wignerMatrices                              = NULL;
    this->integrationWeight                           = 0.0;
    this->maxCompBand                                 = 0;
    this->translationMap                              = NULL;
        
    // ... Control variables
    this->isEmpty                                     = false;
    this->inputOrder                                  = inputO;
    
    //================================================ Sanity checks
    if ( static_cast<proshade_unsign> ( len ) != ( xDmInd * yDmInd * zDmInd ) )
    {
        throw ProSHADE_exception ( "Structure class input map has wrong dimensions.", "EP00044", __FILE__, __LINE__, __func__, "The supplied map array size has different dimensions to\n                    : the required map dimensions." );
    }
    
    if ( ( static_cast<proshade_signed> ( xT - xFr ) != static_cast<proshade_signed> ( xDmInd - 1 ) ) ||
         ( static_cast<proshade_signed> ( yT - yFr ) != static_cast<proshade_signed> ( yDmInd - 1 ) ) ||
         ( static_cast<proshade_signed> ( zT - zFr ) != static_cast<proshade_signed> ( zDmInd - 1 ) ) )
    {
        throw ProSHADE_exception ( "Structure class input dimensions not in line with map\n                    : to/from indices.", "EP00045", __FILE__, __LINE__, __func__, "The supplied map information does not add up. The\n                    : dimensions are not in line with the indexing start/stop\n                    : position distances and therefore proper map indexing\n                    : cannot be done. Please check the input values." );
    }
    
    //================================================ Allocate the map memory
    this->internalMap                                 = new proshade_double [this->xDimIndices * this->yDimIndices * this->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->internalMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy the values into the map
    proshade_unsign arrPos                            = 0;
    for ( proshade_unsign xIt = 0; xIt < this->xDimIndices; xIt++ )
    {
        for ( proshade_unsign yIt = 0; yIt < this->yDimIndices; yIt++ )
        {
            for ( proshade_unsign zIt = 0; zIt < this->zDimIndices; zIt++ )
            {
                arrPos                                = zIt + this->zDimIndices * ( yIt + this->yDimIndices * xIt );
                this->internalMap[arrPos]             = static_cast<proshade_double> ( mapVals[arrPos] );
            }
        }
    }
    
    //================================================ Done
    
}

/*! \brief Destructor for the ProSHADE_data class.
 
    This destructor is responsible for releasing all memory used by the data storing object
 
    \param[out] X N/A.
 */
ProSHADE_internal_data::ProSHADE_data::~ProSHADE_data ( )
{
    //================================================ Release the internal map
    if ( this->internalMap != NULL )
    {
        delete[] this->internalMap;
    }
    
    //================================================ Release the sphere mapping
    if ( this->spheres != NULL )
    {
        for ( proshade_unsign iter = 0; iter < this->noSpheres; iter++ )
        {
            if ( this->spheres[iter] != NULL )
            {
                delete this->spheres[iter];
                this->spheres[iter]                   = NULL;
            }
        }
        delete[] this->spheres;
    }
    
    //================================================ Release the spherical harmonics
    if ( this->sphericalHarmonics != NULL )
    {
        for ( proshade_unsign iter = 0; iter < this->noSpheres; iter++ )
        {
            if ( this->sphericalHarmonics[iter] != NULL )
            {
                delete[] this->sphericalHarmonics[iter];
                this->sphericalHarmonics[iter]        = NULL;
            }
        }
        delete[] this->sphericalHarmonics;
    }
    
    //================================================ Release the rotated spherical harmonics
    if ( this->rotSphericalHarmonics != NULL )
    {
        for ( proshade_unsign iter = 0; iter < this->noSpheres; iter++ )
        {
            if ( this->rotSphericalHarmonics[iter] != NULL )
            {
                delete[] this->rotSphericalHarmonics[iter];
                this->rotSphericalHarmonics[iter]     = NULL;
            }
        }
        delete[] this->rotSphericalHarmonics;
    }
    
    //================================================ Release the RRP matrices (pre-computation for the energy levels descriptor)
    if ( this->rrpMatrices != NULL )
    {
        for ( proshade_unsign bwIt = 0; bwIt < this->maxShellBand; bwIt++ )
        {
            if ( this->rrpMatrices[bwIt] != NULL )
            {
                for ( proshade_unsign shIt = 0; shIt < this->noSpheres; shIt++ )
                {
                    if ( this->rrpMatrices[bwIt][shIt] != NULL )
                    {
                        delete[] this->rrpMatrices[bwIt][shIt];
                    }
                }
                
                delete[] this->rrpMatrices[bwIt];
            }
        }
        
        delete[] this->rrpMatrices;
    }
    
    //================================================ Release the E matrices
    if ( this->eMatrices != NULL )
    {
        for ( proshade_unsign bandIter = 0; bandIter < this->maxCompBand; bandIter++ )
        {
            if ( this->eMatrices[bandIter] != NULL )
            {
                for ( proshade_unsign band2Iter = 0; band2Iter < static_cast<proshade_unsign> ( ( bandIter * 2 ) + 1 ); band2Iter++ )
                {
                    if ( this->eMatrices[bandIter][band2Iter] != NULL )
                    {
                        delete[] this->eMatrices[bandIter][band2Iter];
                    }
                }
                
                delete[] this->eMatrices[bandIter];
            }
        }
        
        delete[] this->eMatrices;
    }
    
    //================================================ Release SOFT and inverse SOFT coefficients
    if ( this->so3Coeffs != NULL )
    {
        delete[] this->so3Coeffs;
    }
    if ( this->so3CoeffsInverse != NULL )
    {
        delete[] this->so3CoeffsInverse;
    }
    
    //================================================ Release Wigner matrices
    if ( this->wignerMatrices != NULL )
    {
        for ( proshade_unsign bandIter = 1; bandIter < this->maxCompBand; bandIter++ )
        {
            if ( this->wignerMatrices[bandIter] != NULL )
            {
                for ( proshade_unsign order1Iter = 0; order1Iter < ( (bandIter * 2) + 1 ); order1Iter++ )
                {
                    if ( this->wignerMatrices[bandIter][order1Iter] != NULL )
                    {
                        delete[] this->wignerMatrices[bandIter][order1Iter];
                    }
                }
                delete[] this->wignerMatrices[bandIter];
            }
        }
        delete[] wignerMatrices;
    }
    
    //================================================ Release translation map
    if ( this->translationMap != NULL )
    {
        delete[] this->translationMap;
    }
    
    //================================================ Release the angle-axis space rotation function
    if ( this->sphereMappedRotFun.size() > 0 )
    {
        for ( proshade_unsign spIt = 0; spIt < static_cast<proshade_unsign> ( this->sphereMappedRotFun.size() ); spIt++ )
        {
            delete this->sphereMappedRotFun.at(spIt);
        }
    }
    
    //================================================ Done
    
}

/*! \brief This function initialises the basic ProSHADE_data variables and reads in a single structure.
 
    This function is basically the constructor for the ProSHADE_data class. It reads in a structure (independent of the structure type) and
    fills in all the appropriate variables of the class.
 
    \param[in] fName The file name of the file which should be loaded.
    \param[in] inputO The order of this structure in this run's input.
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
 */
void ProSHADE_internal_data::ProSHADE_data::readInStructure ( std::string fName, proshade_unsign inputO, ProSHADE_settings* settings )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting to read the structure: " + fName );
    
    //================================================ Check if instance is empty
    if ( !this->isEmpty )
    {
        throw ProSHADE_exception ( "Structure data class not empty.", "E000005", __FILE__, __LINE__, __func__, "Attempted to read in structure into a ProSHADE_data\n                    : object which already does have structure read in\n                    : i.e. " + this->fileName );
    }
    
    //================================================ Save the filename
    this->fileName                                    = fName;
    
    //================================================ Check what is the input format
    this->fileType                                    = ProSHADE_internal_io::figureDataType ( this->fileName );
    
    //================================================ Save input order
    this->inputOrder                                  = inputO;
    
    //================================================ Decide how to proceed
    switch ( this->fileType )
    {
        case ProSHADE_internal_io::UNKNOWN:
            throw ProSHADE_exception ( "Unknown file type.", "E000006", __FILE__, __LINE__, __func__, "When attempting to read the file\n                    : " + this->fileName + "\n                    : the file extension was determined as unknown. This could\n                    : mean either that the file does not exist, or that it is\n                    : not one of the supported extensions." );
            break;
        
        case ProSHADE_internal_io::PDB:
            this->readInPDB                           ( settings );
            break;
        
        case ProSHADE_internal_io::MAP:
            this->readInMAP                           ( settings );
            break;
    }
    
    //================================================ This structure is now full
    this->isEmpty                                     = false;
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Structure read in successfully." );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for reading map data using gemmi library.
 
    This function reads in the map data using the information from the settings object and saves all the results into the
    structure calling it.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
 */
void ProSHADE_internal_data::ProSHADE_data::readInMAP ( ProSHADE_settings* settings )
{
    //================================================ Open the file
    gemmi::Ccp4<float> map;
    map.read_ccp4                                     ( gemmi::MaybeGzipped ( this->fileName.c_str() ) );
    
    //================================================ Convert to XYZ and create complete map, if need be
    map.setup                                         ( gemmi::GridSetup::ReorderOnly, NAN );
    
    //================================================ Read in the rest of the map file header
    ProSHADE_internal_io::readInMapHeader             ( &map,
                                                        &this->xDimIndices,  &this->yDimIndices,  &this->zDimIndices,
                                                        &this->xDimSize,     &this->yDimSize,     &this->zDimSize,
                                                        &this->aAngle,       &this->bAngle,       &this->cAngle,
                                                        &this->xFrom,        &this->yFrom,        &this->zFrom,
                                                        &this->xAxisOrigin,  &this->yAxisOrigin,  &this->zAxisOrigin,
                                                        &this->xAxisOrder,   &this->yAxisOrder,   &this->zAxisOrder,
                                                        &this->xGridIndices, &this->yGridIndices, &this->zGridIndices );
    
    //================================================ Save the map density to ProSHADE variable
    ProSHADE_internal_io::readInMapData               ( &map, this->internalMap, this->xDimIndices, this->yDimIndices, this->zDimIndices, this->xAxisOrder, this->yAxisOrder, this->zAxisOrder );
    
    //================================================ Set resolution if need be
    if ( settings->requestedResolution < 0.0 )
    {
        settings->setResolution                       ( std::min ( static_cast<proshade_double> ( this->xDimSize ) / static_cast<proshade_double> ( this->xDimIndices ),
                                                        std::min ( static_cast<proshade_double> ( this->yDimSize ) / static_cast<proshade_double> ( this->yDimIndices ),
                                                                   static_cast<proshade_double> ( this->zDimSize ) / static_cast<proshade_double> ( this->zDimIndices ) ) ) * 2.0 );
    }
    
    //================================================ Set iterators from and to
    this->figureIndexStartStop                        ( );
    
    //================================================ If specific resolution is requested, make sure the map has it
    if ( settings->changeMapResolution || settings->changeMapResolutionTriLinear )
    {
        this->reSampleMap                             ( settings );
    }
    
    //================================================ Done
    
}

/*! \brief Function for reading pdb data.
 
    This function reads in the pdb data using the information from the settings object, converts the  co-ordinates onto
    a theoretical map and and saves all the results into the structure calling it.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
 
    \warning For multiple models, this function works, but the map is not perfectly fitted to the PDB file.
 */
void ProSHADE_internal_data::ProSHADE_data::readInPDB ( ProSHADE_settings* settings )
{
    //================================================ Set resolution if need be
    if ( settings->requestedResolution < 0.0 )
    {
        settings->setResolution                       ( 8.0 );
    }
    
    //================================================ Open PDB file for reading
    gemmi::Structure pdbFile                          = gemmi::read_structure ( gemmi::MaybeGzipped ( this->fileName ) );
    
    //================================================ Change B-factors if need be
    if ( settings->pdbBFactorNewVal >= 0.0 )
    {
        ProSHADE_internal_mapManip::changePDBBFactors ( &pdbFile, settings->pdbBFactorNewVal, settings->firstModelOnly );
    }
    
    //================================================ Remove waters if required
    if ( settings->removeWaters )
    {
        ProSHADE_internal_mapManip::removeWaters      ( &pdbFile, settings->firstModelOnly );
    }
    
    //================================================ Get PDB COM values
    proshade_double xCOMPdb, yCOMPdb, zCOMPdb;
    ProSHADE_internal_mapManip::findPDBCOMValues      ( pdbFile, &xCOMPdb, &yCOMPdb, &zCOMPdb, settings->firstModelOnly );
    
    //================================================ Find the ranges
    proshade_single xF, xT, yF, yT, zF, zT;
    ProSHADE_internal_mapManip::determinePDBRanges    ( pdbFile, &xF, &xT, &yF, &yT, &zF, &zT, settings->firstModelOnly );

    //================================================ Move ranges to have all FROM values 20
    proshade_single xMov                              = 20.0 - xF;
    proshade_single yMov                              = 20.0 - yF;
    proshade_single zMov                              = 20.0 - zF;
    ProSHADE_internal_mapManip::movePDBForMapCalc     ( &pdbFile, xMov, yMov, zMov, settings->firstModelOnly );

    //================================================ Set the angstrom sizes
    this->xDimSize                                    = xT - xF + 40.0;
    this->yDimSize                                    = yT - yF + 40.0;
    this->zDimSize                                    = zT - zF + 40.0;

    //================================================ Generate map from nicely placed atoms (cell size will be range + 40)
    ProSHADE_internal_mapManip::generateMapFromPDB    ( pdbFile, this->internalMap, settings->requestedResolution, this->xDimSize, this->yDimSize, this->zDimSize, &this->xTo, &this->yTo, &this->zTo, settings->forceP1, settings->firstModelOnly );
    
    //================================================ Set the internal variables to correct values
    this->setPDBMapValues                             ( );
    
    //================================================ Compute reverse movement based on COMs. If there is more than 1 models, simply moving back the xyzMov is not enough.
    proshade_double xCOMMap, yCOMMap, zCOMMap;
    ProSHADE_internal_mapManip::findMAPCOMValues      ( this->internalMap, &xCOMMap, &yCOMMap, &zCOMMap,
                                                        this->xDimSize, this->yDimSize, this->zDimSize,
                                                        this->xFrom, this->xTo, this->yFrom, this->yTo, this->zFrom, this->zTo );
    
    if ( pdbFile.models.size() > 1 )
    {
        xMov                                          = xCOMMap - xCOMPdb;
        yMov                                          = yCOMMap - yCOMPdb;
        zMov                                          = zCOMMap - zCOMPdb;
    }
    
    //================================================ Move map back to the original PDB location
    ProSHADE_internal_mapManip::moveMapByIndices      ( &xMov, &yMov, &zMov, this->xDimSize, this->yDimSize, this->zDimSize,
                                                        &this->xFrom, &this->xTo, &this->yFrom, &this->yTo, &this->zFrom, &this->zTo,
                                                        &this->xAxisOrigin, &this->yAxisOrigin, &this->zAxisOrigin );
    ProSHADE_internal_mapManip::moveMapByFourier      ( this->internalMap, xMov, yMov, zMov, this->xDimSize, this->yDimSize, this->zDimSize,
                                                        this->xDimIndices, this->yDimIndices, this->zDimIndices );
   
    //================================================ If specific resolution is requested, make sure the map has it
    if ( settings->changeMapResolution || settings->changeMapResolutionTriLinear )
    {
        this->reSampleMap                             ( settings );
    }
    
    //================================================ Done
    return;
    
}

/*! \brief Function for determining iterator start and stop positions.
 
    This function is called to set the xFrom, yFrom, ..., yTo and zTo iterator values for easier further calculations.
 */
void ProSHADE_internal_data::ProSHADE_data::setPDBMapValues ( void )
{
    //================================================ Set starts to 0 
    this->xFrom                                       = 0;
    this->yFrom                                       = 0;
    this->zFrom                                       = 0;
    
    //================================================ Set angles to 90 degrees
    this->aAngle                                      = 90.0;
    this->bAngle                                      = 90.0;
    this->cAngle                                      = 90.0;
    
    //================================================ Set dimension sizes in indices
    this->xDimIndices                                 = this->xTo;
    this->yDimIndices                                 = this->yTo;
    this->zDimIndices                                 = this->zTo;
    
    //================================================ Set grid indexing to cell indexing
    this->xGridIndices                                = this->xDimIndices;
    this->yGridIndices                                = this->yDimIndices;
    this->zGridIndices                                = this->zDimIndices;
    
    //================================================ Set axis order
    this->xAxisOrder                                  = 1;
    this->yAxisOrder                                  = 2;
    this->zAxisOrder                                  = 3;
    
    //================================================ Set origin to the first index
    this->xAxisOrigin                                 = this->xFrom;
    this->yAxisOrigin                                 = this->yFrom;
    this->zAxisOrigin                                 = this->zFrom;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for determining iterator start and stop positions.
 
    This function is called to set the xFrom, yFrom, ..., yTo and zTo iterator values for easier further calculations. It assumes that gemmi has read in the xFrom, yFrom and zFrom values already.
 */
void ProSHADE_internal_data::ProSHADE_data::figureIndexStartStop ( void )
{
    //================================================ Set ends to origin + size - 1
    this->xTo                                         = this->xFrom + this->xDimIndices - 1;
    this->yTo                                         = this->yFrom + this->yDimIndices - 1;
    this->zTo                                         = this->zFrom + this->zDimIndices - 1;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for writing out the internal structure representation in MRC MAP format.
 
    This function takes all the internal map representation information from the calling object and proceeds to write all this information in MRC MAP format for
    visualisation and possibly further processing by other software. This function will write out axis order XYZ and spacegroup P1 irrespective of the input
    axis order and spacegroup.
 
    \param[in] fName The filename (including path) to where the output MAP file should be saved.
    \param[in] title String with the map title to be written into the header - default value is "Created by ProSHADE and written by GEMMI"
    \param[in] mode The type of the data, leave at default 2 (mean float type) unless you specifically required other types.
 */
void ProSHADE_internal_data::ProSHADE_data::writeMap ( std::string fName, std::string title, int mode )
{
    //================================================ Create and prepare new Grid gemmi object
    gemmi::Grid<float> mapData;
    mapData.set_unit_cell                             ( this->xDimSize, this->yDimSize, this->zDimSize, this->aAngle, this->bAngle, this->cAngle );
    mapData.set_size_without_checking                 ( this->xDimIndices, this->yDimIndices, this->zDimIndices );
    mapData.axis_order                                = gemmi::AxisOrder::XYZ;
    mapData.spacegroup                                = &gemmi::get_spacegroup_p1();

    //================================================ Create and prepare new Ccp4 gemmi object
    gemmi::Ccp4<float> map;
    map.grid                                          = mapData;
    map.update_ccp4_header                           ( mode );
    
    //================================================ Fill in the header
    ProSHADE_internal_io::writeOutMapHeader           ( &map,
                                                        this->xDimIndices,  this->yDimIndices,  this->zDimIndices,
                                                        this->xDimSize,     this->yDimSize,     this->zDimSize,
                                                        this->aAngle,       this->bAngle,       this->cAngle,
                                                        this->xFrom,        this->yFrom,        this->zFrom,
                                                        this->xAxisOrigin,  this->yAxisOrigin,  this->zAxisOrigin,
                                                        this->xAxisOrder,   this->yAxisOrder,   this->zAxisOrder,
                                                        this->xGridIndices, this->yGridIndices, this->zGridIndices,
                                                        title, mode );
    
    //================================================ Copy internal map to grid
    proshade_unsign arrPos                            = 0;
    for ( proshade_unsign uIt = 0; uIt < this->xDimIndices; uIt++ )
    {
        for ( proshade_unsign vIt = 0; vIt < this->yDimIndices; vIt++ )
        {
            for ( proshade_unsign wIt = 0; wIt < this->zDimIndices; wIt++ )
            {
                arrPos                                = wIt + this->zDimIndices * ( vIt + this->yDimIndices * uIt );
                map.grid.set_value                    ( uIt, vIt, wIt, static_cast<float> ( this->internalMap[arrPos] ) );
            }
        }
    }
    
    //================================================ Update the statistics in the header
    map.update_ccp4_header                            ( mode, true );
    
    //================================================ Write out the map
    map.write_ccp4_map                                ( fName );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function writes out the PDB formatted file coresponding to the structure.

    This function first checks if this internal structure originated from co-ordinate file (only if co-ordinates are provided can they be written out). If so,
    it will proceed to

    \param[in] fName The filename (including path) to where the output PDB file should be saved.
    \param[in] euA The Euler angle alpha by which the co-ordinates should be rotated (leave empty if no rotation is required).
    \param[in] euB The Euler angle beta by which the co-ordinates should be rotated (leave empty if no rotation is required).
    \param[in] euG The Euler angle gamma by which the co-ordinates should be rotated (leave empty if no rotation is required).
    \param[in] transX The translation to be done along the X-axis in Angstroms.
    \param[in] transY The translation to be done along the Y-axis in Angstroms.
    \param[in] transZ The translation to be done along the Z-axis in Angstroms.
    \param[in] firstModel Should only the first model, or rather all of them be used?
*/
void ProSHADE_internal_data::ProSHADE_data::writePdb ( std::string fName, proshade_double euA, proshade_double euB, proshade_double euG, proshade_double transX, proshade_double transY, proshade_double transZ, bool firstModel )
{
    //================================================ Check for co-ordinate origin
    if ( !ProSHADE_internal_io::isFilePDB ( this->fileName ) )
    {
        throw ProSHADE_exception ( "Cannot write co-ordinate file if the input file did not contain co-ordinates.", "EP00047", __FILE__, __LINE__, __func__, "You have called the WritePDB function on structure which\n                    : was created by reading in a map. This is not allowed as\n                    : ProSHADE cannot create co-ordinates from map file." );
    }
    
    //================================================ Open PDB file for reading
    gemmi::Structure pdbFile                          = gemmi::read_structure ( gemmi::MaybeGzipped ( this->fileName ) );
    
    //================================================ If the map was rotated, do the same for the co-ordinates, making sure we take into account the rotation centre of the map
    if ( ( euA != 0.0 ) || ( euB != 0.0 ) || ( euG != 0.0 ) )
    {
        //============================================ Save original PDB COM position
        proshade_double xCOMOriginal = 0.0, yCOMOriginal = 0.0, zCOMOriginal = 0.0;
        ProSHADE_internal_mapManip::findPDBCOMValues  ( pdbFile, &xCOMOriginal, &yCOMOriginal, &zCOMOriginal, firstModel );
        
        //============================================ Compute the rotation centre for the co-ordinates
        proshade_double xRotPos                       = ( ( static_cast<proshade_double> ( this->xDimIndicesOriginal / 2 ) - this->xAxisOriginOriginal ) *
                                                          ( static_cast<proshade_double> ( this->xDimIndicesOriginal - 1 ) / this->xDimSizeOriginal ) ) -
                                                          (this->originalMapXCom - xCOMOriginal);
        proshade_double yRotPos                       = ( ( static_cast<proshade_double> ( this->yDimIndicesOriginal / 2 ) - this->yAxisOriginOriginal ) *
                                                          ( static_cast<proshade_double> ( this->yDimIndicesOriginal - 1 ) / this->yDimSizeOriginal ) ) -
                                                          (this->originalMapYCom - yCOMOriginal);
        proshade_double zRotPos                       = ( ( static_cast<proshade_double> ( this->zDimIndicesOriginal / 2 ) - this->zAxisOriginOriginal ) *
                                                          ( static_cast<proshade_double> ( this->zDimIndicesOriginal - 1 ) / this->zDimSizeOriginal ) ) -
                                                          (this->originalMapZCom - zCOMOriginal);
        
        //============================================ Save the values
        this->originalPdbRotCenX                      = xRotPos;
        this->originalPdbRotCenY                      = yRotPos;
        this->originalPdbRotCenZ                      = zRotPos;
        
        //============================================ Rotate the co-ordinates
        ProSHADE_internal_mapManip::rotatePDBCoordinates ( &pdbFile, euA, euB, euG, xRotPos, yRotPos, zRotPos, firstModel );

        //============================================ Compute the after rotation PDB COM position
        proshade_double xCOMRotated = 0.0, yCOMRotated = 0.0, zCOMRotated = 0.0;
        ProSHADE_internal_mapManip::findPDBCOMValues  ( pdbFile, &xCOMRotated, &yCOMRotated, &zCOMRotated, firstModel );

        //============================================ Compute the after rotation position correction
        proshade_double xPDBTrans                     = ( xCOMRotated - xCOMOriginal ) + ( this->mapPostRotXCom - this->originalMapXCom );
        proshade_double yPDBTrans                     = ( yCOMRotated - yCOMOriginal ) + ( this->mapPostRotYCom - this->originalMapYCom );
        proshade_double zPDBTrans                     = ( zCOMRotated - zCOMOriginal ) + ( this->mapPostRotZCom - this->originalMapZCom );

        //============================================ Correct the co-ordinate position after rotation
        ProSHADE_internal_mapManip::translatePDBCoordinates ( &pdbFile, xPDBTrans, yPDBTrans, zPDBTrans, firstModel );
    }

    //================================================ Save the values
    this->originalPdbTransX                           = this->comMovX + transX + this->originalPdbRotCenX;
    this->originalPdbTransY                           = this->comMovY + transY + this->originalPdbRotCenY;
    this->originalPdbTransZ                           = this->comMovZ + transZ + this->originalPdbRotCenZ;
    
    //================================================ Translate by required translation and the map centering (if applied)
    ProSHADE_internal_mapManip::translatePDBCoordinates ( &pdbFile, this->comMovX + transX, this->comMovY + transY, this->comMovZ + transZ, firstModel );

    //================================================ Write the PDB file
    std::ofstream outCoOrdFile;
    outCoOrdFile.open                                 ( fName.c_str() );
    
    if ( outCoOrdFile.is_open() )
    {
        gemmi::PdbWriteOptions opt;
        write_pdb                                     ( pdbFile, outCoOrdFile, opt );
    }
    else
    {
        std::stringstream hlpMessage;
        hlpMessage << "Failed to open the PDB file " << fName << " for output.";
        throw ProSHADE_exception ( hlpMessage.str().c_str(), "EP00048", __FILE__, __LINE__, __func__, "ProSHADE has failed to open the PDB output file. This is\n                    : likely caused by either not having the write privileges\n                    : to the required output path, or by making a mistake in\n                    : the path." );
    }
    
    outCoOrdFile.close                                ( );
    
    //================================================ Done
    return ;
}

/*! \brief Function for writing out a mask in MRC MAP format.
 
    This function takes a mask map and the filename and proceeds to write out the mask into the requested file name in th
    MRC MAP format. It assumes that the mask has the same dimmensions as the map.
 
    \param[in] fileName The filename (including path) to where the output should be saved.
    \param[in] mask Pointer to the mask map array.
 */
void ProSHADE_internal_data::ProSHADE_data::writeMask ( std::string fName, proshade_double* mask )
{
    //================================================ Allocate the memory
    proshade_double* hlpMap                           = new proshade_double[this->xDimIndices * this->yDimIndices * this->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation ( hlpMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy original map and over-write with the mask
    for ( proshade_unsign iter = 0; iter < ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        hlpMap[iter]                                  = this->internalMap[iter];
        this->internalMap[iter]                       = mask[iter];
    }
    
    //================================================ Write out the mask
    this->writeMap                                    ( fName );
    
    //================================================ Copy the original map values back
    for ( proshade_unsign iter = 0; iter < ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        this->internalMap[iter]                       = hlpMap[iter];
    }
    
    //================================================ Release memory
    delete[] hlpMap;
    
    //================================================ Done
    return ;
}

/*! \brief Function for inverting the map to its mirror image.
 
    This function switches all index values along the three axes from 0 ... max to max ... 0. This should not
    normally be done, but in the case where the wrong hand has been used in the map re-construction process, this
    may be helpful.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
 */
void ProSHADE_internal_data::ProSHADE_data::invertMirrorMap ( ProSHADE_settings* settings )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Map inversion." );
    
    //================================================ Initialise variables
    proshade_signed arrayPos, invPos;
    
    //================================================ Create helper map
    proshade_double* hlpMap                           = new proshade_double [this->xDimIndices * this->yDimIndices * this->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation     ( hlpMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Save map values to the helper map
    for ( proshade_unsign iter = 0; iter < ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        hlpMap[iter]                                  = this->internalMap[iter];
    }
    
    //================================================ Invert the values
    for ( proshade_signed xIt = 0; xIt < static_cast<proshade_signed> ( this->xDimIndices ); xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < static_cast<proshade_signed> ( this->yDimIndices ); yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < static_cast<proshade_signed> ( this->zDimIndices ); zIt++ )
            {
                //==================================== Var init
                arrayPos                              = zIt + this->zDimIndices * ( yIt + this->yDimIndices * xIt );
                invPos                                = ( (this->zDimIndices-1) - zIt ) + this->zDimIndices * ( ( (this->yDimIndices-1) - yIt ) + this->yDimIndices * ( (this->xDimIndices-1) - xIt ) );
                
                //==================================== And save
                this->internalMap[invPos]             = hlpMap[arrayPos];
            }
        }
    }
    
    //================================================ Release memory
    delete[] hlpMap;
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Map inversion completed." );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for normalising the map to mean 0 and sd 1..
 
    This function takes the map and changes its value to have mean 0 and standard deviation of 1. This should make
    wo maps with very different density levels more comparable, but it remains to be seen if this causes any trouble.
    Can be turned off using the settings options.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
 */
void ProSHADE_internal_data::ProSHADE_data::normaliseMap ( ProSHADE_settings* settings )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Map normalisation." );
    
    //================================================ Initialise vector of map values
    std::vector<proshade_double> mapVals              ( this->xDimIndices * this->yDimIndices * this->zDimIndices, 0.0 );
    
    //================================================ Get all map values
    for ( proshade_unsign iter = 0; iter < ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        mapVals.at(iter)                              = this->internalMap[iter];
    }
    
    //================================================ Get mean and sd
    proshade_double* meanSD                           = new proshade_double[2];
    ProSHADE_internal_maths::vectorMeanAndSD          ( &mapVals, meanSD );
    
    //================================================ Normalise the values
    for ( proshade_unsign iter = 0; iter < ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        this->internalMap[iter]                       = ( this->internalMap[iter] - meanSD[0] ) / meanSD[1];
    }
    
    //================================================ Clear the vector
    mapVals.clear                                     ( );
    
    //================================================ Release memory
    delete[] meanSD;
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Map normalisation completed." );
    
    //================================================ Done
    return ;
    
}


/*! \brief Function for computing the map mask using blurring and X IQRs from median.
 
    This function takes all the internal map representation information from the calling object and the internal map itself
    and proceeds to write all this information in MRC MAP format for visualisation and further processing by other software.
    It is dependent on the internal information being correct.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
 */
void ProSHADE_internal_data::ProSHADE_data::maskMap ( ProSHADE_settings* settings )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Computing mask." );
    
    //================================================ Initialise the blurred map
    proshade_double* blurredMap                       = new proshade_double[this->xDimIndices * this->yDimIndices * this->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation ( blurredMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Compute blurred map
    ProSHADE_internal_mapManip::blurSharpenMap        ( this->internalMap, blurredMap, this->xDimIndices, this->yDimIndices, this->zDimIndices,
                                                        this->xDimSize, this->yDimSize, this->zDimSize, settings->blurFactor );
    
    //================================================ Compute mask from blurred map and save it into the original map
    ProSHADE_internal_mapManip::getMaskFromBlurr      ( blurredMap, this->internalMap, this->xDimIndices, this->yDimIndices, this->zDimIndices, settings->maskingThresholdIQRs );
    
    //================================================ Print the mask if need be
    if ( settings->saveMask ) {  if ( settings->maskFileName == "" ) { this->writeMask ( "proshade_mask.map", blurredMap ); } else { std::stringstream ss; ss << settings->maskFileName << "_" << this->inputOrder << ".map"; this->writeMask ( ss.str(), blurredMap ); } }
    
    //================================================ Release memory
    delete[] blurredMap;
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Mask computed." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the boundaries enclosing positive map values and adds some extra space.
 
    This function firstly finds the boundaries which enclose the positive map values and then it proceeds to add a
    given amount of space to all dimensions (positive and negative) to make sure the map does not end exactly at the
    bounds. It returns the new boundaries in the ret variable if they are smaller than the original bounds, or just
    the original bounds in case decrease was not achieved.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
    \param[in] ret A pointer to proshade_signed array of 6 storing the results - (0 = minX; 1 = maxX; 2 = minY; 3 = maxY; 4 - minZ; 5 = maxZ).
 */
void ProSHADE_internal_data::ProSHADE_data::getReBoxBoundaries ( ProSHADE_settings* settings, proshade_signed*& ret )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Finding new boundaries." );
    
    //================================================ If same bounds as first one are required, test if possible and return these instead
    if ( settings->useSameBounds && ( this->inputOrder != 0 ) )
    {
        for ( proshade_unsign iter = 0; iter < 6; iter++ ) { ret[iter] = settings->forceBounds[iter]; }
    }
    //================================================ In this case, bounds need to be found de novo
    else
    {
        //============================================ Find the non-zero bounds
        ProSHADE_internal_mapManip::getNonZeroBounds  ( this->internalMap, this->xDimIndices, this->yDimIndices, this->zDimIndices,
                                                        this->xDimSize, this->yDimSize, this->zDimSize, ret );
        
        //============================================ Add the extra space
        ProSHADE_internal_mapManip::addExtraBoundSpace ( this->xDimIndices, this->yDimIndices, this->zDimIndices,
                                                         this->xDimSize, this->yDimSize, this->zDimSize, ret, settings->boundsExtraSpace );
        
        //============================================ Beautify boundaries
        ProSHADE_internal_mapManip::beautifyBoundaries ( ret, this->xDimIndices, this->yDimIndices, this->zDimIndices, settings->boundsSimilarityThreshold, settings->verbose );
        
        //============================================ Report function results
        std::stringstream ssHlp;
        ssHlp << "New boundaries are: " << ret[1] - ret[0] + 1 << " x " << ret[3] - ret[2] + 1 << " x " << ret[5] - ret[4] + 1;
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 3, ssHlp.str() );
        
        //============================================ If need be, save boundaries to be used for all other structure
        if ( settings->useSameBounds && ( this->inputOrder == 0 ) )
        {
            for ( proshade_unsign iter = 0; iter < 6; iter++ ) { settings->forceBounds[iter] = ret[iter]; }
        }
    }
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "New boundaries determined." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the boundaries enclosing positive map values and adds some extra space, returning values in Python friendly format.
 
    This function firstly finds the boundaries which enclose the positive map values and then it proceeds to add a
    given amount of space to all dimensions (positive and negative) to make sure the map does not end exactly at the
    bounds. It returns the new boundaries in the ret variable if they are smaller than the original bounds, or just
    the original bounds in case decrease was not achieved. The return format is Python friendly, otherwise it is the same
    function as the getReBoxBoundaries() function.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
    \param[in] ret A pointer to proshade_signed array of 6 storing the results - (0 = minX; 1 = maxX; 2 = minY; 3 = maxY; 4 - minZ; 5 = maxZ).
    \param[in] len The length of the ret array - must always be 6, but Swig/Numpy requires this number to be passed.
 */
void ProSHADE_internal_data::ProSHADE_data::getReBoxBoundariesPy ( ProSHADE_settings* settings, int* reBoxBounds, int len )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Finding new boundaries." );
    
    //================================================ Initialise internal variables
    proshade_signed* ret                              = NULL;
    ret                                               = new proshade_signed[6];
    ProSHADE_internal_misc::checkMemoryAllocation     ( ret, __FILE__, __LINE__, __func__ );
    
    //================================================ If same bounds as first one are required, test if possible and return these instead
    if ( settings->useSameBounds && ( this->inputOrder != 0 ) )
    {
        for ( proshade_unsign iter = 0; iter < 6; iter++ ) { ret[iter] = settings->forceBounds[iter]; }
    }
    //================================================ In this case, bounds need to be found de novo
    else
    {
        //============================================ Find the non-zero bounds
        ProSHADE_internal_mapManip::getNonZeroBounds ( this->internalMap, this->xDimIndices, this->yDimIndices, this->zDimIndices, this->xDimSize, this->yDimSize, this->zDimSize, ret );
        
        //============================================ Add the extra space
        ProSHADE_internal_mapManip::addExtraBoundSpace ( this->xDimIndices, this->yDimIndices, this->zDimIndices, this->xDimSize, this->yDimSize, this->zDimSize, ret, settings->boundsExtraSpace );
        
        //============================================ Beautify boundaries
        ProSHADE_internal_mapManip::beautifyBoundaries ( ret, this->xDimIndices, this->yDimIndices, this->zDimIndices, settings->boundsSimilarityThreshold, settings->verbose );
        
        //============================================ Report function results
        std::stringstream ssHlp;
        ssHlp << "New boundaries are: " << ret[1] - ret[0] + 1 << " x " << ret[3] - ret[2] + 1 << " x " << ret[5] - ret[4] + 1;
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 3, ssHlp.str() );
        
        //============================================ If need be, save boundaries to be used for all other structure
        if ( settings->useSameBounds && ( this->inputOrder == 0 ) )
        {
            for ( proshade_unsign iter = 0; iter < 6; iter++ ) { settings->forceBounds[iter] = ret[iter]; }
        }
    }
    
    //================================================ Copy internal to output array
    reBoxBounds[0]                                    = static_cast<int> ( ret[0] );
    reBoxBounds[1]                                    = static_cast<int> ( ret[1] );
    reBoxBounds[2]                                    = static_cast<int> ( ret[2] );
    reBoxBounds[3]                                    = static_cast<int> ( ret[3] );
    reBoxBounds[4]                                    = static_cast<int> ( ret[4] );
    reBoxBounds[5]                                    = static_cast<int> ( ret[5] );
    
    //================================================ Release memory
    delete[] ret;
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "New boundaries determined." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function creates a new structure from the calling structure and new bounds values.
 
    This function takes a pointer to uninitialised structure and fills it with the calling structure values adjusted for
    the new bounds. If the bounds are the same, the two structures should be identical except the file (the new structure
    does not have an input file associated) and the type (no type for the new structure). It can deal with both larger and
    smaller bounds than the original values.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
    \param[in] newStr A pointer reference to a new structure class which has all the same values except for the new bounds and adequately changed map.
    \param[in] newBounds A pointer to proshade_signed array of 6 storing the results - (0 = minX; 1 = maxX; 2 = minY; 3 = maxY; 4 - minZ; 5 = maxZ).
 */
void ProSHADE_internal_data::ProSHADE_data::createNewMapFromBounds ( ProSHADE_settings* settings, ProSHADE_data*& newStr, proshade_signed* newBounds )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Creating new structure according to the new  bounds." );
    
    //================================================ Fill in basic info
    newStr->fileName                                  = "N/A";
    newStr->fileType                                  = ProSHADE_internal_io::MAP;
    
    //================================================ Fill in new structure values
    newStr->xDimIndices                               = static_cast<proshade_signed> ( newBounds[1] ) - static_cast<proshade_signed> ( newBounds[0] ) + 1;
    newStr->yDimIndices                               = static_cast<proshade_signed> ( newBounds[3] ) - static_cast<proshade_signed> ( newBounds[2] ) + 1;
    newStr->zDimIndices                               = static_cast<proshade_signed> ( newBounds[5] ) - static_cast<proshade_signed> ( newBounds[4] ) + 1;
            
    newStr->aAngle                                    = this->aAngle;
    newStr->bAngle                                    = this->aAngle;
    newStr->cAngle                                    = this->aAngle;
            
    newStr->xDimSize                                  = static_cast<proshade_single> ( newStr->xDimIndices ) * ( this->xDimSize / static_cast<proshade_single> ( this->xDimIndices ) );
    newStr->yDimSize                                  = static_cast<proshade_single> ( newStr->yDimIndices ) * ( this->yDimSize / static_cast<proshade_single> ( this->yDimIndices ) );
    newStr->zDimSize                                  = static_cast<proshade_single> ( newStr->zDimIndices ) * ( this->zDimSize / static_cast<proshade_single> ( this->zDimIndices ) );
            
    newStr->xGridIndices                              = newStr->xDimIndices;
    newStr->yGridIndices                              = newStr->yDimIndices;
    newStr->zGridIndices                              = newStr->zDimIndices;
            
    newStr->xAxisOrder                                = this->xAxisOrder;
    newStr->yAxisOrder                                = this->yAxisOrder;
    newStr->zAxisOrder                                = this->zAxisOrder;
            
    newStr->xAxisOrigin                               = this->xAxisOrigin + newBounds[0];
    newStr->yAxisOrigin                               = this->yAxisOrigin + newBounds[2];
    newStr->zAxisOrigin                               = this->zAxisOrigin + newBounds[4];
            
    newStr->xFrom                                     = this->xFrom + newBounds[0];
    newStr->yFrom                                     = this->yFrom + newBounds[2];
    newStr->zFrom                                     = this->zFrom + newBounds[4];
            
    newStr->xTo                                       = this->xTo - ( (this->xDimIndices-1) - newBounds[1] );
    newStr->yTo                                       = this->yTo - ( (this->yDimIndices-1) - newBounds[3] );
    newStr->zTo                                       = this->zTo - ( (this->zDimIndices-1) - newBounds[5] );
    
    //================================================ Allocate new structure map
    newStr->internalMap                               = new proshade_double[newStr->xDimIndices * newStr->yDimIndices * newStr->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation     ( newStr->internalMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy the map
    ProSHADE_internal_mapManip::copyMapByBounds       ( newStr->xFrom, newStr->xTo, newStr->yFrom, newStr->yTo, newStr->zFrom, newStr->zTo,
                                                        this->xFrom, this->yFrom, this->zFrom, newStr->yDimIndices, newStr->zDimIndices,
                                                        this->xDimIndices, this->yDimIndices, this->zDimIndices, newStr->internalMap, this->internalMap );
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "New structure created." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function creates a new structure from the calling structure and new bounds values in Python friendly format.
 
    This function takes a pointer to uninitialised structure and fills it with the calling structure values adjusted for
    the new bounds. If the bounds are the same, the two structures should be identical except the file (the new structure
    does not have an input file associated) and the type (no type for the new structure). It can deal with both larger and
    smaller bounds than the original values. This version of the function is identical to the createNewMapFromBounds()
    function, but it takes the input boundaries in Python friendly fassion.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
    \param[in] newStr A pointer reference to a new structure class which has all the same values except for the new bounds and adequately changed map.
    \param[in] newBounds A pointer to proshade_signed array of 6 storing the results - (0 = minX; 1 = maxX; 2 = minY; 3 = maxY; 4 - minZ; 5 = maxZ).
    \param[in] len The length of the newBounds array - this must always be 6, but Swig/Numpy requires this number to be passed.
 */
void ProSHADE_internal_data::ProSHADE_data::createNewMapFromBoundsPy ( ProSHADE_settings* settings, ProSHADE_data* newStr, int* newBounds, int len )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Creating new structure according to the new  bounds." );
    
    //================================================ Fill in basic info
    newStr->fileName                                  = "N/A";
    newStr->fileType                                  = ProSHADE_internal_io::MAP;
    
    //================================================ Fill in new structure values
    newStr->xDimIndices                               = static_cast<proshade_signed> ( newBounds[1] ) - static_cast<proshade_signed> ( newBounds[0] ) + 1;
    newStr->yDimIndices                               = static_cast<proshade_signed> ( newBounds[3] ) - static_cast<proshade_signed> ( newBounds[2] ) + 1;
    newStr->zDimIndices                               = static_cast<proshade_signed> ( newBounds[5] ) - static_cast<proshade_signed> ( newBounds[4] ) + 1;
            
    newStr->aAngle                                    = this->aAngle;
    newStr->bAngle                                    = this->aAngle;
    newStr->cAngle                                    = this->aAngle;
            
    newStr->xDimSize                                  = static_cast<proshade_single> ( newStr->xDimIndices ) * ( this->xDimSize / static_cast<proshade_single> ( this->xDimIndices ) );
    newStr->yDimSize                                  = static_cast<proshade_single> ( newStr->yDimIndices ) * ( this->yDimSize / static_cast<proshade_single> ( this->yDimIndices ) );
    newStr->zDimSize                                  = static_cast<proshade_single> ( newStr->zDimIndices ) * ( this->zDimSize / static_cast<proshade_single> ( this->zDimIndices ) );
            
    newStr->xGridIndices                              = newStr->xDimIndices;
    newStr->yGridIndices                              = newStr->yDimIndices;
    newStr->zGridIndices                              = newStr->zDimIndices;
            
    newStr->xAxisOrder                                = this->xAxisOrder;
    newStr->yAxisOrder                                = this->yAxisOrder;
    newStr->zAxisOrder                                = this->zAxisOrder;
            
    newStr->xAxisOrigin                               = this->xAxisOrigin + newBounds[0];
    newStr->yAxisOrigin                               = this->yAxisOrigin + newBounds[2];
    newStr->zAxisOrigin                               = this->zAxisOrigin + newBounds[4];
            
    newStr->xFrom                                     = this->xFrom + newBounds[0];
    newStr->yFrom                                     = this->yFrom + newBounds[2];
    newStr->zFrom                                     = this->zFrom + newBounds[4];
            
    newStr->xTo                                       = this->xTo - ( (this->xDimIndices-1) - newBounds[1] );
    newStr->yTo                                       = this->yTo - ( (this->yDimIndices-1) - newBounds[3] );
    newStr->zTo                                       = this->zTo - ( (this->zDimIndices-1) - newBounds[5] );
    
    //================================================ Allocate new structure map
    newStr->internalMap                               = new proshade_double[newStr->xDimIndices * newStr->yDimIndices * newStr->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation     ( newStr->internalMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy the map
    ProSHADE_internal_mapManip::copyMapByBounds       ( newStr->xFrom, newStr->xTo, newStr->yFrom, newStr->yTo, newStr->zFrom, newStr->zTo,
                                                        this->xFrom, this->yFrom, this->zFrom, newStr->yDimIndices, newStr->zDimIndices,
                                                        this->xDimIndices, this->yDimIndices, this->zDimIndices, newStr->internalMap, this->internalMap );
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "New structure created." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function changes the internal map sampling to conform to particular resolution value.
 
    This function will take the requested resolution value from the settings object and will proceed to change the internal
    map sampling to conform to requested resolution / 2 and therefore to the requested resolution map.
 
    \param[in] settings A pointer to settings class containing all the information required for reading in the map.
 */
void ProSHADE_internal_data::ProSHADE_data::reSampleMap ( ProSHADE_settings* settings )
{
    //================================================ Initialise the return variable
    proshade_single* changeVals                       = new proshade_single[6];
    
    //================================================ Now re-sample the map
    if ( settings->changeMapResolution )
    {
        ProSHADE_internal_mapManip::reSampleMapToResolutionFourier ( this->internalMap, settings->requestedResolution, this->xDimIndices, this->yDimIndices, this->zDimIndices,
                                                                     this->xDimSize, this->yDimSize, this->zDimSize, changeVals );
        
        if ( settings->changeMapResolutionTriLinear )
        {
            ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Requested both Fourier-space and real-space map re-sampling. Defaulting to only Fourier space re-samplling.", "WM00049" );
        }
    }
    if ( settings->changeMapResolutionTriLinear && !settings->changeMapResolution )
    {
        ProSHADE_internal_mapManip::reSampleMapToResolutionTrilinear ( this->internalMap, settings->requestedResolution, this->xDimIndices, this->yDimIndices, this->zDimIndices,
                                                                       this->xDimSize, this->yDimSize, this->zDimSize, changeVals );

    }
    
    //================================================ Set the internal values to reflect the new map size
    this->xDimIndices                                += static_cast<proshade_unsign>  ( changeVals[0] );
    this->yDimIndices                                += static_cast<proshade_unsign>  ( changeVals[1] );
    this->zDimIndices                                += static_cast<proshade_unsign>  ( changeVals[2] );
            
    this->xGridIndices                                = this->xDimIndices;
    this->yGridIndices                                = this->yDimIndices;
    this->zGridIndices                                = this->zDimIndices;
            
    this->xTo                                        += static_cast<proshade_unsign>  ( changeVals[0] );
    this->yTo                                        += static_cast<proshade_unsign>  ( changeVals[1] );
    this->zTo                                        += static_cast<proshade_unsign>  ( changeVals[2] );
            
    this->xDimSize                                    = changeVals[3];
    this->yDimSize                                    = changeVals[4];
    this->zDimSize                                    = changeVals[5];

    //================================================ Figure how much the new map moved
    proshade_single xMov                              = -( ( this->xFrom * ( this->xDimSize / static_cast<proshade_single> ( this->xDimIndices - changeVals[0] ) ) ) -
                                                           ( this->xFrom * ( this->xDimSize / static_cast<proshade_single> ( this->xDimIndices ) ) ) );
    proshade_single yMov                              = -( ( this->yFrom * ( this->yDimSize / static_cast<proshade_single> ( this->yDimIndices - changeVals[1] ) ) ) -
                                                           ( this->yFrom * ( this->yDimSize / static_cast<proshade_single> ( this->yDimIndices ) ) ) );
    proshade_single zMov                              = -( ( this->zFrom * ( this->zDimSize / static_cast<proshade_single> ( this->zDimIndices - changeVals[2] ) ) ) -
                                                           ( this->zFrom * ( this->zDimSize / static_cast<proshade_single> ( this->zDimIndices ) ) ) );

    //================================================ Move by indices (this should be sufficient)
    ProSHADE_internal_mapManip::moveMapByIndices      ( &xMov, &yMov, &zMov, this->xDimSize, this->yDimSize, this->zDimSize, &this->xFrom, &this->xTo,
                                                        &this->yFrom, &this->yTo, &this->zFrom, &this->zTo, &this->xAxisOrigin, &this->yAxisOrigin, &this->zAxisOrigin );
    
    ProSHADE_internal_mapManip::moveMapByFourier      ( this->internalMap, xMov, yMov, zMov, this->xDimSize, this->yDimSize, this->zDimSize,
                                                        this->xDimIndices, this->yDimIndices, this->zDimIndices );
    
    //================================================ Release memory
    delete[] changeVals;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function shits the map so that its COM is in the centre of the map.
 
    This function finds the Centre Of Mass (COM) for the internal map and proceeds to use Fourier to shift
    the COM to the centre of the map. There is an assumption that the COM and centre of map are close, as if they
    were far apart, the shift could move some part of the map through the boundaries and cause the map to become
    messy.
 
    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
 */
void ProSHADE_internal_data::ProSHADE_data::centreMapOnCOM ( ProSHADE_settings* settings )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Centering map onto its COM." );
    
    //================================================ Initialise local variables
    proshade_unsign arrPos                            = 0;
    proshade_single xCOM                              = 0.0;
    proshade_single yCOM                              = 0.0;
    proshade_single zCOM                              = 0.0;
    proshade_single totDens                           = 0.0;
    
    //================================================ Find the COM location
    for ( proshade_unsign xIt = 0; xIt < this->xDimIndices; xIt++ )
    {
        for ( proshade_unsign yIt = 0; yIt < this->yDimIndices; yIt++ )
        {
            for ( proshade_unsign zIt = 0; zIt < this->zDimIndices; zIt++ )
            {
                //==================================== Get index
                arrPos                                = zIt + this->zDimIndices * ( yIt + this->yDimIndices * xIt );
                
                //==================================== Get COM
                if ( this->internalMap[arrPos] > 0.0 )
                {
                    xCOM                             += static_cast<proshade_single> ( this->internalMap[arrPos] * xIt );
                    yCOM                             += static_cast<proshade_single> ( this->internalMap[arrPos] * yIt );
                    zCOM                             += static_cast<proshade_single> ( this->internalMap[arrPos] * zIt );
                    totDens                          += static_cast<proshade_single> ( this->internalMap[arrPos] );
                }
            }
        }
    }
    xCOM                                             /= totDens;
    yCOM                                             /= totDens;
    zCOM                                             /= totDens;
    
    //================================================ Find distance from COM to map centre in Angstroms
    proshade_single xDist                             = ( static_cast<proshade_single> ( this->xDimIndices / 2.0 ) - xCOM ) * static_cast<proshade_single> ( this->xDimSize / this->xDimIndices );
    proshade_single yDist                             = ( static_cast<proshade_single> ( this->yDimIndices / 2.0 ) - yCOM ) * static_cast<proshade_single> ( this->yDimSize / this->yDimIndices );
    proshade_single zDist                             = ( static_cast<proshade_single> ( this->zDimIndices / 2.0 ) - zCOM ) * static_cast<proshade_single> ( this->zDimSize / this->zDimIndices );
    
    //================================================ Save COM movement
    this->comMovX                                     = static_cast<proshade_double> ( xDist );
    this->comMovY                                     = static_cast<proshade_double> ( yDist );
    this->comMovZ                                     = static_cast<proshade_double> ( zDist );
    
    //================================================ Move the map within the box
    ProSHADE_internal_mapManip::moveMapByFourier      ( this->internalMap, xDist, yDist, zDist, this->xDimSize, this->yDimSize, this->zDimSize, this->xDimIndices, this->yDimIndices, this->zDimIndices );
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Map centered." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function increases the size of the map so that it can add empty space around it.
 
    This function adds a given number of angstroms (as given in the settings object) to the internal structure map. This requires all the internal
    variables to be adjusted for the extra space at the begginning and at the end, while also copying the map into a larger one with appropriate
    extra space.
 
    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
 */
void ProSHADE_internal_data::ProSHADE_data::addExtraSpace ( ProSHADE_settings* settings )
{
    //================================================ Report function start
    std::stringstream hlpSS;
    hlpSS << "Adding extra " << settings->addExtraSpace << " angstroms.";
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1,  hlpSS.str() );
    
    //================================================ Figure how much indices need to change
    proshade_unsign xAddIndices                       = ProSHADE_internal_mapManip::myRound ( settings->addExtraSpace / static_cast<proshade_single> ( this->xDimSize / this->xDimIndices ) );
    proshade_unsign yAddIndices                       = ProSHADE_internal_mapManip::myRound ( settings->addExtraSpace / static_cast<proshade_single> ( this->yDimSize / this->yDimIndices ) );
    proshade_unsign zAddIndices                       = ProSHADE_internal_mapManip::myRound ( settings->addExtraSpace / static_cast<proshade_single> ( this->zDimSize / this->zDimIndices ) );
    
    //================================================ Update internal data variables
    this->xDimSize                                   += static_cast<proshade_single> ( 2 * xAddIndices ) * static_cast<proshade_single> ( this->xDimSize / this->xDimIndices );
    this->yDimSize                                   += static_cast<proshade_single> ( 2 * yAddIndices ) * static_cast<proshade_single> ( this->yDimSize / this->yDimIndices );
    this->zDimSize                                   += static_cast<proshade_single> ( 2 * zAddIndices ) * static_cast<proshade_single> ( this->zDimSize / this->zDimIndices );
            
    this->xDimIndices                                += 2 * xAddIndices;
    this->yDimIndices                                += 2 * yAddIndices;
    this->zDimIndices                                += 2 * zAddIndices;
            
    this->xGridIndices                                = this->xDimIndices;
    this->yGridIndices                                = this->yDimIndices;
    this->zGridIndices                                = this->zDimIndices;
            
    this->xAxisOrigin                                -= xAddIndices;
    this->yAxisOrigin                                -= yAddIndices;
    this->zAxisOrigin                                -= zAddIndices;
            
    this->xFrom                                      -= xAddIndices;
    this->yFrom                                      -= yAddIndices;
    this->zFrom                                      -= zAddIndices;
            
    this->xTo                                        += xAddIndices;
    this->yTo                                        += yAddIndices;
    this->zTo                                        += zAddIndices;
    
    //================================================ Allocate new map
    proshade_double* newMap                           = new proshade_double[this->xDimIndices * this->yDimIndices * this->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation     ( newMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Set new map to zeroes
    for ( proshade_unsign iter = 0; iter < ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        newMap[iter]                                  = 0.0;
    }
    
    //================================================ Update the map
    proshade_unsign newMapIndex, oldMapIndex;
    for ( proshade_unsign xIt = 0; xIt < (this->xDimIndices - xAddIndices); xIt++ )
    {
        //============================================ Check if point is applicable
        if ( xIt < xAddIndices ) { continue; }
        
        for ( proshade_unsign yIt = 0; yIt < (this->yDimIndices - yAddIndices); yIt++ )
        {
            //======================================== Check if point is applicable
            if ( yIt < yAddIndices ) { continue; }
            
            for ( proshade_unsign zIt = 0; zIt < (this->zDimIndices - zAddIndices); zIt++ )
            {
                //==================================== Check if point is applicable
                if ( zIt < zAddIndices ) { continue; }
                
                //==================================== Var init
                newMapIndex                           = zIt + this->zDimIndices * ( yIt + this->yDimIndices * xIt );
                oldMapIndex                           = (zIt - zAddIndices) + (this->zDimIndices - ( 2 * zAddIndices ) ) * ( (yIt - yAddIndices) + (this->yDimIndices - ( 2 * yAddIndices ) ) * (xIt - xAddIndices) );
                
                newMap[newMapIndex]                   = this->internalMap[oldMapIndex];
            }
        }
    }
    
    //================================================ Copy new to old
    delete[] this->internalMap;
    
    this->internalMap                                 = new proshade_double[this->xDimIndices * this->yDimIndices * this->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->internalMap, __FILE__, __LINE__, __func__ );
    
    for ( proshade_unsign iter = 0; iter < ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        this->internalMap[iter]                       = newMap[iter];
    }
    
    //================================================ Release memory
    delete[] newMap;
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Extra space added." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function simply clusters several other functions which should be called together.
 
    This function serves to cluster the map normalisation, map masking, map centering and map extra space addition
    into a single function. This allows for simpler code and does not take any control away, as all the decisions
    are ultimately driven by the settings.
 
    This function also does some internal value saving and auto-determination of any parameters that the user did not
    supply. This, however, means, that this function MUST be called for every structure that is to be processed by
    ProSHADE. This is of importance to people whe want to use only a perticular functions.
 
    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
 
    \warning This function MUST be called on any structure that is to be processed by ProSHADE.
 */
void ProSHADE_internal_data::ProSHADE_data::processInternalMap ( ProSHADE_settings* settings )
{
    //================================================ Invert map
    if ( settings->invertMap ) { this->invertMirrorMap ( settings ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Map inversion (mirror image) not requested." ); }
 
    //================================================ Normalise map
    if ( settings->normaliseMap ) { this->normaliseMap ( settings ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Map normalisation not requested." ); }

    //================================================ Compute mask
    //if ( settings->maskMap ) { if ( settings->useCorrelationMasking ) { this->maskMapCorrelation ( settings ); } else { this->maskMap ( settings ); } }
    if ( settings->maskMap ) { this->maskMap ( settings ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Masking not requested." ); }

    //================================================ Centre map
    if ( settings->moveToCOM ) { this->centreMapOnCOM ( settings ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Map centering not requested." ); }
    
    //================================================ Add extra space
    if ( settings->addExtraSpace != 0.0 ) { this->addExtraSpace ( settings ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Extra space not requested." ); }
    
    //================================================ Remove phase, if required
    if ( !settings->usePhase ) { this->removePhaseInormation ( settings ); ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Phase information removed from the data." ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Phase information retained in the data." ); }
    
    //================================================ Compute and save the original map values
    this->setOriginalMapValues                        ( );
    
    //================================================ Set settings values which were left on AUTO by user and will not be set later
    settings->setVariablesLeftOnAuto                  ( );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function determines the sphere positions (radii) for sphere mapping.

    This function determines the radii of the concentric spheres (as measured from the centre of the map). This is
    done by checking if these values have already been as and if not, then the radii are placed between points of the
    map starting between the centre point and its neighbours and then adding spheres until the most outlying diagonal
    point is covered.

    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
 */
void ProSHADE_internal_data::ProSHADE_data::getSpherePositions ( ProSHADE_settings* settings )
{
    //================================================ Check the current settings value is set to auto
    if ( this->spherePos.size() != 0 )
    {
        std::stringstream hlpSS;
        hlpSS << "The sphere distances were determined as " << this->spherePos.at(0);
        for ( proshade_unsign iter = 1; iter < static_cast<proshade_unsign> ( this->spherePos.size() ); iter++ ) { hlpSS << "; " << this->spherePos.at(iter); }
        hlpSS << " Angstroms.";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 3, hlpSS.str() );
        return ;
    }
    
    //================================================ Find maximum diagonal
    proshade_unsign maxDim                            = std::max ( this->xDimSize, std::max ( this->yDimSize, this->zDimSize ) );
    proshade_unsign minDim                            = std::min ( this->xDimSize, std::min ( this->yDimSize, this->zDimSize ) );
    proshade_unsign midDim                            = 0;
    if      ( ( this->xDimSize < maxDim ) && ( this->xDimSize > minDim ) ) { midDim = this->xDimSize; }
    else if ( ( this->yDimSize < maxDim ) && ( this->yDimSize > minDim ) ) { midDim = this->yDimSize; }
    else                                                                   { midDim = this->zDimSize; }
    
    proshade_single maxDiag                           = std::sqrt ( std::pow ( static_cast<proshade_single> ( maxDim ), 2.0 ) +
                                                                    std::pow ( static_cast<proshade_single> ( midDim ), 2.0 ) );
    
    //================================================ Set between the points
    for ( proshade_single iter = 0.5; ( iter * settings->maxSphereDists ) < ( maxDiag / 2.0 ); iter += 1.0 )
    {
        ProSHADE_internal_misc::addToSingleVector     ( &this->spherePos, ( iter * settings->maxSphereDists ) );
    }
    
    //================================================ Save the number of spheres
    this->noSpheres                                   = static_cast<proshade_unsign> ( this->spherePos.size() );
    
    //================================================ Report progress
    std::stringstream hlpSS;
    hlpSS << "The sphere distances were determined as " << this->spherePos.at(0);
    for ( proshade_unsign iter = 1; iter < static_cast<proshade_unsign> ( this->spherePos.size() ); iter++ ) { hlpSS << "; " << this->spherePos.at(iter); }
    hlpSS << " Angstroms.";
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 3, hlpSS.str() );
    
    //================================================ Done
    return ;

}


/*! \brief This function simply clusters several other functions which should be called together.
 
    This function serves to cluster the map normalisation, map masking, map centering and map extra space addition
    into a single function. This allows for simpler code and does not take any control away, as all the decisions
    are ultimately driven by the settings.
 
    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
 */
void ProSHADE_internal_data::ProSHADE_data::mapToSpheres ( ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting sphere mapping procedure." );
    
    //================================================ Determine spherical harmonics variables
    settings->determineAllSHValues                    ( this->xDimIndices, this->yDimIndices, this->zDimIndices );
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Sphere settings determined." );
    
    //================================================ Find number of spheres supported
    this->getSpherePositions                          ( settings );
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 2, "Sphere positions obtained." );
    
    //================================================ Create sphere objects and map the density
    this->spheres                                     = new ProSHADE_internal_spheres::ProSHADE_sphere* [ this->noSpheres ];
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( this->spherePos.size() ); iter++ )
    {
        std::stringstream ss;
        ss << "Now mapping sphere " << iter << " .";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 4, ss.str() );
        
        this->spheres[iter]                           = new ProSHADE_internal_spheres::ProSHADE_sphere ( this->xDimIndices, this->yDimIndices, this->zDimIndices,
                                                                                                         this->xDimSize, this->yDimSize, this->zDimSize, iter,
                                                                                                        &this->spherePos, settings->progressiveSphereMapping, settings->maxBandwidth,
                                                                                                         this->internalMap, &this->maxShellBand );
    }
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Sphere mapping procedure completed." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the spherical harmonics decomposition for the whole structure.
 
    This function is called to compute the spherical harmonics decomposition of the mapped data on every available
    sphere. This is done sphere-wise and there is some sub-optimal memory management stemming from different shells
    having different resolutions.
 
    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
 */
void ProSHADE_internal_data::ProSHADE_data::computeSphericalHarmonics ( ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting spherical harmonics decomposition." );
    
    //================================================ Initialise memory
    this->sphericalHarmonics                          = new proshade_complex* [this->noSpheres];
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->sphericalHarmonics, __FILE__, __LINE__, __func__ );
    for ( proshade_unsign iter = 0; iter < this->noSpheres; iter++ )
    {
        this->sphericalHarmonics[iter]                = new proshade_complex [(this->spheres[iter]->getLocalBandwidth() * 2) * (this->spheres[iter]->getLocalBandwidth() * 2)];
        ProSHADE_internal_misc::checkMemoryAllocation ( this->sphericalHarmonics[iter], __FILE__, __LINE__, __func__ );
    }
    
    //================================================ Compute the spherical harmonics
    for ( proshade_unsign iter = 0; iter < this->noSpheres; iter++ )
    {
        //============================================ Report progress
        std::stringstream ss;
        ss << "Now decomposing sphere " << iter << ". " << "( Band is: " << this->spheres[iter]->getLocalBandwidth() << ").";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 4, ss.str() );
        
        //============================================ Compute
        ProSHADE_internal_sphericalHarmonics::computeSphericalHarmonics ( this->spheres[iter]->getLocalBandwidth(), this->spheres[iter]->getMappedData(), this->sphericalHarmonics[iter] );
    }
    
    //================================================ Report completion
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 2, "Spherical harmonics decomposition complete." );
    
    //======================================== Done
    return ;
    
}

/*! \brief This function runs the symmetry detection algorithms on this structure and saves the results in the settings object.
 
    This function is the symmetry detection starting point. It decides whether a specific symmetry is being sought after, or whether
    a general symmetry search is required. Consequently, it calls the appropriate functions and ends up with saving the resulting
    predictions into the settings object supplied. It also saves all the detected symmetry groups to the settings object for further
    processing and programmatical access.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] axes A vector to which all the axes of the recommended symmetry (if any) will be saved.
    \param[in] allCs A vector to which all the detected cyclic symmetries will be saved into.
 */
void ProSHADE_internal_data::ProSHADE_data::detectSymmetryInStructure ( ProSHADE_settings* settings, std::vector< proshade_double* >* axes, std::vector < std::vector< proshade_double > >* allCs )
{
    //================================================ Initialise variables
    std::vector< proshade_double* > CSyms             = this->getCyclicSymmetriesList ( settings );
    
    //================================================ Was any particular symmetry requested?
    if ( settings->requestedSymmetryType == "" )
    {
        //============================================ Run the symmetry detection functions for C, D, T, O and I symmetries
        std::vector< proshade_double* > DSyms         = this->getDihedralSymmetriesList ( settings, &CSyms );
        std::vector< proshade_double* > ISyms         = this->getIcosahedralSymmetriesList ( settings, &CSyms );
        std::vector< proshade_double* > OSyms; std::vector< proshade_double* > TSyms;
        if ( ISyms.size() < 31 ) {  OSyms = this->getOctahedralSymmetriesList ( settings, &CSyms ); if ( OSyms.size() < 13 ) { TSyms = this->getTetrahedralSymmetriesList ( settings, &CSyms ); } }
        
        //============================================ Decide on recommended symmetry
        this->saveRecommendedSymmetry                 ( settings, &CSyms, &DSyms, &TSyms, &OSyms, &ISyms, axes );
    }
    
    if ( settings->requestedSymmetryType == "C" )
    {
        //============================================ Run only the C symmetry detection and search for requested fold
        this->saveRequestedSymmetryC                  ( settings, &CSyms, axes );
    }
    
    if ( settings->requestedSymmetryType == "D" )
    {
        //============================================ Run only the D symmetry detection and search for requested fold
        std::vector< proshade_double* > DSyms         = this->getDihedralSymmetriesList ( settings, &CSyms );
        this->saveRequestedSymmetryD                  ( settings, &DSyms, axes );
    }
    
    if ( settings->requestedSymmetryType == "T" )
    {
        //============================================ Run only the T symmetry detection and search for requested fold
        std::vector< proshade_double* > TSyms         = this->getTetrahedralSymmetriesList ( settings, &CSyms );
        settings->setRecommendedFold                  ( 0 );
        if ( TSyms.size() == 7 )              { settings->setRecommendedSymmetry ( "T" ); for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( TSyms.size() ); it++ ) { settings->setDetectedSymmetry ( TSyms.at(it) ); ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, TSyms.at(it) ); } }
        else                                  { settings->setRecommendedSymmetry ( "" ); }
    }
    
    if ( settings->requestedSymmetryType == "O" )
    {
        //============================================ Run only the O symmetry detection and search for requested fold
        std::vector< proshade_double* > OSyms         = this->getOctahedralSymmetriesList ( settings, &CSyms );
        settings->setRecommendedFold                  ( 0 );
        if ( OSyms.size() == 13 )             { settings->setRecommendedSymmetry ( "O" ); for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( OSyms.size() ); it++ ) { settings->setDetectedSymmetry ( OSyms.at(it) ); ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, OSyms.at(it) ); } }
        else                                  { settings->setRecommendedSymmetry ( "" ); }
    }
    
    if ( settings->requestedSymmetryType == "I" )
    {
        //============================================ Run only the T symmetry detection and search for requested fold
        std::vector< proshade_double* > ISyms         = this->getIcosahedralSymmetriesList ( settings, &CSyms );
        settings->setRecommendedFold                  ( 0 );
        if ( ISyms.size() == 31 )             { settings->setRecommendedSymmetry ( "I" ); for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( ISyms.size() ); it++ ) { settings->setDetectedSymmetry ( ISyms.at(it) ); ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, ISyms.at(it) ); } }
        else                                  { settings->setRecommendedSymmetry ( "" ); }
    }
    
    if ( ( settings->requestedSymmetryType != "" ) && ( settings->requestedSymmetryType != "C" ) && ( settings->requestedSymmetryType != "D" ) && ( settings->requestedSymmetryType != "T" ) && ( settings->requestedSymmetryType != "O" ) && ( settings->requestedSymmetryType != "I" ) )
    {
        throw ProSHADE_exception ( "Requested symmetry supplied, but not recognised.", "ES00032", __FILE__, __LINE__, __func__, "There are only the following value allowed for the\n                    : symmetry type request: \"C\", \"D\", \"T\", \"O\" and \"I\". Any\n                    : other value will result in this error." );
    }
    
    //================================================ Save C symmetries to argument and if different from settings, to the settings as well
    bool isArgSameAsSettings                          = true;
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSyms.size() ); cIt++ )
    {
        std::vector< proshade_double > nextSym;
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms.at(cIt)[0] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms.at(cIt)[1] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms.at(cIt)[2] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms.at(cIt)[3] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms.at(cIt)[4] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms.at(cIt)[5] );
        ProSHADE_internal_misc::addToDoubleVectorVector ( allCs, nextSym );

        if ( ( cIt == 0 ) && ( settings->allDetectedCAxes.size() == 0 ) ) { isArgSameAsSettings = false; }
        if ( !isArgSameAsSettings ) { ProSHADE_internal_misc::addToDoubleVectorVector ( &settings->allDetectedCAxes, nextSym ); }

        nextSym.clear                                 ( );
    }

    //================================================ Release memory
    for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( CSyms.size() ); it++ ) { delete[] CSyms.at(it); }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function runs the symmetry detection algorithms on this structure saving the axes in the settings object only.
 
    This function runs the detectSymmetryInStructure() function without requiring the vector of double pointers arguments, so that
    it is simply callable from Python. The axes are saved in the settings object for later retrieval.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
 */
void ProSHADE_internal_data::ProSHADE_data::detectSymmetryInStructurePython ( ProSHADE_settings* settings )
{
    //================================================ Run the algorithm
    this->detectSymmetryInStructure                   ( settings, &settings->detectedSymmetry, &settings->allDetectedCAxes );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function runs the symmetry detection algorithms on this structure using the angle-axis space and saving the results in the settings object.
 
    This function firstly decides whether specific C symmetry was requested or not. This decision is important as knowing the required fold allows for a rather
    simplified algorithm to be applied. Thus, if specific fold is known, simplified algorithm will be used. Otherwise, this function will do a general search by firstly
    finding all cyclic point groups and then applying the dihedral, tetrahedral, octahedral and icosahedral searches.
 
    Once complete, the function will save both, the vector of ProSHADE formatted array pointers as well as vector of vectors of doubles with the same information
    containing all detected cyclic point groups into the supplied vector pointers.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] axes A vector to which all the axes of the recommended symmetry (if any) will be saved.
    \param[in] allCs A vector to which all the detected cyclic symmetries will be saved into.
 */
void ProSHADE_internal_data::ProSHADE_data::detectSymmetryFromAngleAxisSpace ( ProSHADE_settings* settings, std::vector< proshade_double* >* axes, std::vector < std::vector< proshade_double > >* allCs )
{
    //================================================ Modify axis tolerance by sampling, if required by user
    if ( settings->axisErrToleranceDefault )
    {
        settings->axisErrTolerance                    = std::max ( 0.01, ( 2.0 * M_PI ) / this->maxShellBand );
    }
    
    //================================================  If C was requested, we will do it immediately - this allows for a significant speed-up.
    if ( settings->requestedSymmetryType == "C" )
    {
        //============================================ Report progress
        std::stringstream hlpSS;
        hlpSS << "Starting detection of cyclic point group C" << settings->requestedSymmetryFold;
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, hlpSS.str() );
        
        //============================================ Do simplified search only in the applicable data
        proshade_double symThres                      = 0.0;
        std::vector< proshade_double* > CSyms         = this->findRequestedCSymmetryFromAngleAxis ( settings, settings->requestedSymmetryFold, &symThres );
        
        //============================================ Save the best axis as the recommended one
        if ( settings->detectedSymmetry.size() == 0 ) { if ( CSyms.size() > 0 ) { settings->setDetectedSymmetry ( CSyms.at(0) ); } }
        if ( CSyms.size() > 0 )
        {
            settings->setRecommendedSymmetry          ( "C" );
            settings->setRecommendedFold              ( settings->requestedSymmetryFold );
            
            ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSyms.at(0) );
            this->saveDetectedSymmetries              ( settings, &CSyms, allCs );
        }
        else
        {
            settings->setRecommendedSymmetry          ( "" );
            settings->setRecommendedFold              ( 0 );
        }
        
        //============================================ Done
        return ;
    }
    
    //============================================ Report progress
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Starting C symmetry detection." );

    //================================================ Initialise variables
    std::vector< proshade_double* > CSyms             = getCyclicSymmetriesListFromAngleAxis ( settings );
    
    //================================================ Sanity check - was the rotation function mapped properly?
    if ( this->sphereMappedRotFun.size() < 1 )
    {
        throw ProSHADE_exception ( "Rotation function was not converted into angle-axis space.", "ES00062", __FILE__, __LINE__, __func__, "It seems that the convertRotationFunction() function was\n                    : not yet called. Therefore, there are no data to detect the\n                    : symmetry from; please call the convertRotationFunction()\n                    : function before the detectSymmetryFromAngleAxisSpace()\n                    : function." );
    }
    
    //================================================ Sanity check - was any symmetry requested?
    if ( ( settings->requestedSymmetryType != "" ) && ( settings->requestedSymmetryType != "C" ) && ( settings->requestedSymmetryType != "D" ) && ( settings->requestedSymmetryType != "T" ) && ( settings->requestedSymmetryType != "O" ) && ( settings->requestedSymmetryType != "I" ) )
    {
        throw ProSHADE_exception ( "Requested symmetry supplied, but not recognised.", "ES00032", __FILE__, __LINE__, __func__, "There are only the following value allowed for the\n                    : symmetry type request: \"C\", \"D\", \"T\", \"O\" and \"I\". Any\n                    : other value will result in this error." );
    }
    
    //================================================ Are we doing general search?
    if ( settings->requestedSymmetryType == "" )
    {
        //============================================ Run the symmetry detection functions for C, D, T, O and I symmetries
        std::vector< proshade_double* > DSyms         = this->getDihedralSymmetriesList ( settings, &CSyms );
        std::vector< proshade_double* > ISyms         = this->getPredictedIcosahedralSymmetriesList ( settings, &CSyms );
//        std::vector< proshade_double* > ISyms         = this->getIcosahedralSymmetriesList ( settings, &CSyms );
        std::vector< proshade_double* > OSyms; std::vector< proshade_double* > TSyms;
        if ( ISyms.size() < 31 ) {  OSyms = this->getOctahedralSymmetriesList ( settings, &CSyms ); if ( OSyms.size() < 13 ) { TSyms = this->getTetrahedralSymmetriesList ( settings, &CSyms ); } }
        
        //============================================ Decide on recommended symmetry
        this->saveRecommendedSymmetry                 ( settings, &CSyms, &DSyms, &TSyms, &OSyms, &ISyms, axes );
    }
    
    if ( settings->requestedSymmetryType == "D" )
    {
        //============================================ Run only the D symmetry detection and search for requested fold
        std::vector< proshade_double* > DSyms         = this->getDihedralSymmetriesList ( settings, &CSyms );
        this->saveRequestedSymmetryD                  ( settings, &DSyms, axes );
    }
    
    if ( settings->requestedSymmetryType == "T" )
    {
        //============================================ Run only the T symmetry detection and search for requested fold
        std::cerr << "Sadly, this functionality is not yet implemented. Please use the -z option to use the original peak searching symmetry detection algorithm." << std::endl;
//        std::vector< proshade_double* > TSyms         = this->getTetrahedralSymmetriesList ( settings, &CSyms );
//        settings->setRecommendedFold                  ( 0 );
//        if ( TSyms.size() == 7 ) { settings->setRecommendedSymmetry ( "T" ); for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( TSyms.size() ); it++ ) { settings->setDetectedSymmetry ( TSyms.at(it) ); ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, TSyms.at(it) ); } }
//        else                     { settings->setRecommendedSymmetry ( "" ); }
    }
    
    if ( settings->requestedSymmetryType == "O" )
    {
        std::cerr << "Sadly, this functionality is not yet implemented. Please use the -z option to use the original peak searching symmetry detection algorithm." << std::endl;
    }
    
    if ( settings->requestedSymmetryType == "I" )
    {
        std::cerr << "Sadly, this functionality is not yet implemented. Please use the -z option to use the original peak searching symmetry detection algorithm." << std::endl;
    }
    
    //================================================ Save C symmetries to argument and if different from settings, to the settings as well
    this->saveDetectedSymmetries                      ( settings, &CSyms, allCs );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the results of point group searches and saves then into the output variables.
 
    This function takes the CSyms as they are returned by the findRequestedCSymmetryFromAngleAxis() or the
    getCyclicSymmetriesListFromAngleAxis() functions and re-saves then to the output variables of the detectSymmetryFromAngleAxisSpace()
    function. It also releases the memory of the CSyms argument.
 
    \warning This function releases the memory of the CSyms argument.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] CSyms A pointer to vector
    \param[in] axes A pointer to a vector to which all the axes of the recommended symmetry (if any) will be saved.
    \param[in] allCs A pointer to a vector to which all the detected cyclic symmetries will be saved into.
 */
void ProSHADE_internal_data::ProSHADE_data::saveDetectedSymmetries ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSyms, std::vector < std::vector< proshade_double > >* allCs )
{
    //================================================ Initialise variables
    bool isArgSameAsSettings                          = true;
    
    //================================================ For each detected point group
    for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( CSyms->size() ); cIt++ )
    {
        //============================================ Create vector to replace the pointer
        std::vector< proshade_double > nextSym;
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms->at(cIt)[0] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms->at(cIt)[1] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms->at(cIt)[2] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms->at(cIt)[3] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms->at(cIt)[4] );
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms->at(cIt)[5] );
        ProSHADE_internal_misc::addToDoubleVectorVector ( allCs, nextSym );

        //============================================ Copy the vector to output variable and if different, then also to settings object
        if ( ( cIt == 0 ) && ( settings->allDetectedCAxes.size() == 0 ) ) { isArgSameAsSettings = false; }
        if ( !isArgSameAsSettings ) { ProSHADE_internal_misc::addToDoubleVectorVector ( &settings->allDetectedCAxes, nextSym ); }

        //============================================ Release memory
        nextSym.clear                                 ( );
        delete[] CSyms->at(cIt);
    }
    
    //================================================ Done
    return ;
    
}


/*! \brief This function locates the best scoring C symmetry axis, returning the score and best symmetry index.
 
    This function takes the list of detected C symmetries and decides which of them is the best, taking into account the folds and the average heights. This
    is not the best approach, I would look into MLE of symmetry presence givent the height and fold, but for now this should do.
 
    \param[in] CSym This is the complete list of the ProSHADE detected C axes in the ProSHADE format.
    \param[in] symInd A pointer to variable where the best symmetry axis index will be stored.
    \param[out] ret The score of the best scoring C axis.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::findBestCScore ( std::vector< proshade_double* >* CSym, proshade_unsign* symInd )
{
    //================================================ Sanity check
    if ( CSym->size() == 0 ) { *symInd = 0; return ( 0.0 ); }
    
    //================================================ Initalise variables
    proshade_double ret                               = CSym->at(0)[5];
   *symInd                                            = 0;
    proshade_double frac                              = 0.0;
    
    //================================================ Check all other axes
// THIS NEEDS TO BE IMPROVED USING THE MAXIMUM LIKELIHOOD FOR THIS FOLD
    for ( proshade_unsign ind = 1; ind < static_cast<proshade_unsign>( CSym->size() ); ind++ )
    {
        //============================================ If higher fold than already leading one (do not care for lower fold and lower average height axes)
        if ( CSym->at(ind)[0] > CSym->at(*symInd)[0] )
        {
            //======================================== How much higher fold is it? Also, adding some protection against large syms supported only by a subset and a minimum requirement.
            frac                                      = ( std::abs( CSym->at(ind)[5]- 0.5 ) / std::abs( CSym->at(*symInd)[5] - 0.5 ) ) / ( CSym->at(*symInd)[0] / CSym->at(ind)[0] );
 
            //======================================== Check if the new is "better" according to this criteria.
            if ( frac >= 1.0 && ( ( CSym->at(*symInd)[5] * 0.85 ) < CSym->at(ind)[5] ) )
            {
                //==================================== And it is! Save and try next one.
               *symInd                                = ind;
                ret                                   = CSym->at(ind)[5];
            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function locates the best scoring D symmetry axis, returning the score and best symmetry index.
 
    This function takes the list of detected D symmetries and decides which of them is the best, taking into account the folds and the average heights. This
    is not the best approach, I would look into MLE of symmetry presence givent the height and folds, but for now this should do.
 
    \param[in] DSym This is the complete list of the ProSHADE detected D axes in the ProSHADE format.
    \param[in] symInd A pointer to variable where the best symmetry axis index will be stored.
    \param[out] ret The score of the best scoring D axis.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::findBestDScore ( std::vector< proshade_double* >* DSym, proshade_unsign* symInd )
{
    //================================================ Sort the vector
    std::sort                                         ( DSym->begin(), DSym->end(), ProSHADE_internal_misc::sortDSymHlpInv );
    
    //================================================ Initalise variables
    proshade_double ret                               = 0.0;
    proshade_double frac                              = 0.0;
    if ( DSym->size() > 0 )
    {
        ret                                           = ( ( DSym->at(0)[0] * DSym->at(0)[5] ) + ( DSym->at(0)[6] * DSym->at(0)[11] ) ) / ( DSym->at(0)[0] + DSym->at(0)[6] );
       *symInd                                        = 0;
    }
    else { return ( ret ); }

    //================================================ Check all other axes
// THIS NEEDS TO BE IMPROVED USING THE MAXIMUM LIKELIHOOD FOR THIS FOLD
    for ( proshade_unsign ind = 1; ind < static_cast<proshade_unsign>( DSym->size() ); ind++ )
    {
        //============================================ If higher fold than already leading one (do not care for lower fold and lower average height axes)
        if ( ( DSym->at(ind)[0] + DSym->at(ind)[6] ) > ( DSym->at(*symInd)[0] + DSym->at(*symInd)[6] ) )
        {
            //======================================== How much higher fold is it? Also, adding some protection against large syms supported only by a subset and a minimum requirement.
            frac                                      = std::max ( std::min ( ( ( DSym->at(*symInd)[0] + DSym->at(*symInd)[6] ) / ( DSym->at(ind)[0] + DSym->at(ind)[6] ) ) * 1.5, 0.9 ), 0.6 );

            //======================================== Check if the new is "better" according to this criteria.
            if ( ( ( ( DSym->at(*symInd)[0] * DSym->at(*symInd)[5] ) + ( DSym->at(*symInd)[6] * DSym->at(*symInd)[11] ) ) / ( DSym->at(*symInd)[0] + DSym->at(*symInd)[6] ) * frac ) < ( ( DSym->at(ind)[0] * DSym->at(ind)[5] ) + ( DSym->at(ind)[6] * DSym->at(ind)[11] ) ) / ( DSym->at(ind)[0] + DSym->at(ind)[6] ) )
            {
                //==================================== And it is! Save and try next one.
               *symInd                                = ind;
                ret                                   = ( ( DSym->at(ind)[0] * DSym->at(ind)[5] ) + ( DSym->at(ind)[6] * DSym->at(ind)[11] ) ) / ( DSym->at(ind)[0] + DSym->at(ind)[6] );
            }
        }
    }

    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function takes the list of tetrahedral axes and returns a score for deciding whether T symmetry should be recommended.
 
    This function simply checks if the complete T symmetry is present (returning 0.0 if not). If present, the function will compute the
    fold weighted average axis height for the whole symmetry and return this number.
 
    \param[in] TSym This is the complete list of the ProSHADE detected T symmetry axes in the ProSHADE format.
    \param[out] ret The score of the T symmetry.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::findTScore ( std::vector< proshade_double* >* TSym )
{
    //================================================ Initialise variables
    proshade_double ret                               = 0.0;
    proshade_double foldSum                           = 0.0;

    //================================================ Check the T symmetry for being complete
    if ( TSym->size() == 7 )
    {
        //============================================ Compute the weighted fold
        for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( TSym->size() ); cIt++ )
        {
            ret                                      += TSym->at(cIt)[0] * TSym->at(cIt)[5];
            foldSum                                  += TSym->at(cIt)[0];
        }
        
        //============================================ Weight
        ret                                          /= foldSum;
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function takes the list of octahedral axes and returns a score for deciding whether O symmetry should be recommended.
 
    This function simply checks if the complete O symmetry is present (returning 0.0 if not). If present, the function will compute the
    fold weighted average axis height for the whole symmetry and return this number.
 
    \param[in] OSym This is the complete list of the ProSHADE detected O symmetry axes in the ProSHADE format.
    \param[out] ret The score of the O symmetry.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::findOScore ( std::vector< proshade_double* >* OSym )
{
    //================================================ Initialise variables
    proshade_double ret                               = 0.0;
    proshade_double foldSum                           = 0.0;
    
    //================================================ Check the O symmetry for being complete
    if ( OSym->size() == 13 )
    {
        //============================================ Compute the weighted fold
        for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( OSym->size() ); cIt++ )
        {
            ret                                      += OSym->at(cIt)[0] * OSym->at(cIt)[5];
            foldSum                                  += OSym->at(cIt)[0];
        }
        
        //============================================ Weight
        ret                                          /= foldSum;
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function takes the list of icosahedral axes and returns a score for deciding whether I symmetry should be recommended.
 
    This function simply checks if the complete I symmetry is present (returning 0.0 if not). If present, the function will compute the
    fold weighted average axis height for the whole symmetry and return this number.
 
    \param[in] ISym This is the complete list of the ProSHADE detected I symmetry axes in the ProSHADE format.
    \param[out] ret The score of the I symmetry.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::findIScore ( std::vector< proshade_double* >* ISym )
{
    //================================================ Initialise variables
    proshade_double ret                               = 0.0;
    proshade_double foldSum                           = 0.0;
    
    //================================================ Check the T symmetry for being complete
    if ( ISym->size() == 31 )
    {
        //============================================ Compute the weighted fold
        for ( proshade_unsign cIt = 0; cIt < static_cast<proshade_unsign> ( ISym->size() ); cIt++ )
        {
            ret                                      += ISym->at(cIt)[0] * ISym->at(cIt)[5];
            foldSum                                  += ISym->at(cIt)[0];
        }
        
        //============================================ Weight
        ret                                          /= foldSum;
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function takes all the detected symmetry results and decides on which are to be recommended for this structure.
 
    This function starts by obtaining the scores (fold weighted height averages) for each of the detectable symmetry types. Then, it proceeds
    to compute which of these should be recommended by ProSHADE based on a little shaky combination of axes number and score. This part needs to
    be improved by using ML estimation, when I get the time.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] CSym A vector of pointers to double arrays, each array being a single Cyclic symmetry entry.
    \param[in] DSym A vector of pointers to double arrays, each array being a single Dihedral symmetry entry.
    \param[in] TSym A vector of pointers to double arrays, all of which together form the axes of tetrahedral symmetry.
    \param[in] OSym A vector of pointers to double arrays, all of which together form the axes of octahedral symmetry.
    \param[in] ISym A vector of pointers to double arrays, all of which together form the axes of icosahedral symmetry.
    \param[in] axes A vector to which all the axes of the recommended symmetry (if any) will be saved.
 */
void ProSHADE_internal_data::ProSHADE_data::saveRecommendedSymmetry ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSym, std::vector< proshade_double* >* DSym, std::vector< proshade_double* >* TSym, std::vector< proshade_double* >* OSym, std::vector< proshade_double* >* ISym, std::vector< proshade_double* >* axes )
{
    //================================================ Initialise variables
    proshade_double cScore = 0.0, dScore = 0.0, tScore = 0.0, oScore = 0.0, iScore = 0.0;
    proshade_unsign bestCIndex, bestDIndex;
    
    //================================================ Find a score for each input symmetry type.
    cScore                                            = this->findBestCScore ( CSym, &bestCIndex );
    dScore                                            = this->findBestDScore ( DSym, &bestDIndex );
    tScore                                            = this->findTScore     ( TSym );
    oScore                                            = this->findOScore     ( OSym );
    iScore                                            = this->findIScore     ( ISym );

    //================================================ Find the best available score
    proshade_double bestWeightedScore                 = std::max ( cScore, std::max ( dScore * 1.1, std::max ( tScore * 3.0, std::max ( oScore * 4.0, iScore * 5.0 ) ) ) );
    
    //================================================ No score? Well, no symmetry.
    if ( bestWeightedScore < 0.05 ) { settings->setRecommendedSymmetry ( "" ); return; }
    
    if ( bestWeightedScore == cScore )
    {
        settings->setRecommendedSymmetry              ( "C" );
        settings->setRecommendedFold                  ( CSym->at(bestCIndex)[0] );
        ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSym->at(bestCIndex) );
        if ( settings->detectedSymmetry.size() == 0 ) { settings->setDetectedSymmetry ( CSym->at(bestCIndex) ); }
        
        //============================================ Warn if resolution does not really support this fold
        if ( ( ( 360.0 / static_cast<double> ( CSym->at(bestCIndex)[0] ) ) - ( 360.0 / static_cast<double> ( CSym->at(bestCIndex)[0] + 1 ) ) ) <
             ( 360.0 / static_cast<double> ( settings->maxBandwidth * 4.0 ) ) )
        {
            std::stringstream hlpSS;
            hlpSS << "!!! ProSHADE WARNING !!! Reporting symmetry C" << CSym->at(bestCIndex)[0] << ", however, the grid sampling does not provide reasonable accuracy for symmetry with such high fold and therefore ProSHADE cannot responsibly claim this symmetry to be correct. It is suggested that the grid sampling is increased for more accurate symmetry detection. (Set higher resolution using -r).";
            ProSHADE_internal_messages::printWarningMessage ( settings->verbose, hlpSS.str(), "WS00054" );
        }
    }
    if ( bestWeightedScore == dScore * 1.1 )
    {
        settings->setRecommendedSymmetry              ( "D" );
        settings->setRecommendedFold                  ( std::max ( DSym->at(bestDIndex)[0], DSym->at(bestDIndex)[6] ) );
        ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, DSym->at(bestDIndex) );
        ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, &DSym->at(bestDIndex)[6] );
        if ( settings->detectedSymmetry.size() == 0 )
        {
            settings->setDetectedSymmetry             ( DSym->at(bestDIndex) );
            settings->setDetectedSymmetry             ( &DSym->at(bestDIndex)[6] );
        }
        
        //============================================ Warn if resolution does not really support this fold
        if ( ( ( 360.0 / static_cast<double> ( std::max ( DSym->at(bestDIndex)[0], DSym->at(bestDIndex)[6] ) ) ) - ( 360.0 / static_cast<double> ( std::max ( DSym->at(bestDIndex)[0], DSym->at(bestDIndex)[6] ) + 1 ) ) ) <
             ( 360.0 / static_cast<double> ( settings->maxBandwidth * 4.0 ) ) )
        {
            std::stringstream hlpSS;
            hlpSS << "!!! ProSHADE WARNING !!! Reporting symmetry D" << std::max ( DSym->at(bestDIndex)[0], DSym->at(bestDIndex)[6] ) << ", however, the grid sampling does not provide reasonable accuracy for symmetry with such high fold and therefore ProSHADE cannot responsibly claim this symmetry to be correct. It is suggested that the grid sampling is increased for more accurate symmetry detection. (Set higher resolution using -r).";
            ProSHADE_internal_messages::printWarningMessage ( settings->verbose, hlpSS.str(), "WS00054" );
        }
    }
    if ( bestWeightedScore == tScore * 3.0 )
    {
        settings->setRecommendedSymmetry              ( "T" );
        settings->setRecommendedFold                  ( 0 );
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( TSym->size() ); it++ ) { ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, TSym->at(it) ); }
        if ( settings->detectedSymmetry.size() == 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( TSym->size() ); it++ ) { settings->setDetectedSymmetry ( TSym->at(it) ); } }
    }
    if ( bestWeightedScore == oScore * 4.0 )
    {
        settings->setRecommendedSymmetry              ( "O" );
        settings->setRecommendedFold                  ( 0 );
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( OSym->size() ); it++ ) { ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, OSym->at(it) ); }
        if ( settings->detectedSymmetry.size() == 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( OSym->size() ); it++ ) { settings->setDetectedSymmetry ( OSym->at(it) ); } }
    }
    if ( bestWeightedScore == iScore * 5.0 )
    {
        settings->setRecommendedSymmetry              ( "I" );
        settings->setRecommendedFold                  ( 0 );
        for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( ISym->size() ); it++ ) { ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, ISym->at(it) ); }
        if ( settings->detectedSymmetry.size() == 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( ISym->size() ); it++ ) { settings->setDetectedSymmetry ( ISym->at(it) ); } }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the C symmetries and searched for the requested symmetry.
 
    This is a simple search function, which searches the symmetry results for the requested symmetry fold, and if more such
    symmetries are found, takes the one with the highest average peak height. If the requested fold was found, it will save
    it to the settings object, while it will set the object to fold 0 if the requested symmetry was not found (although there
    may be other symmetries present).
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] CSym A vector of pointers to double arrays, each array being a single Cyclic symmetry entry.
    \param[in] axes A vector to which all the axes of the requested symmetry (if any) will be saved.
 */
void ProSHADE_internal_data::ProSHADE_data::saveRequestedSymmetryC ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSym, std::vector< proshade_double* >* axes )
{
    //================================================ Initialise variables
    proshade_unsign bestIndex                         = 0;
    proshade_double highestSym                        = 0.0;
    
    //================================================ Search for best fold
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( CSym->size() ); iter++ )
    {
        //============================================ Check if it is tbe correct fold
        if ( CSym->at(iter)[0] != settings->requestedSymmetryFold ) { continue; }
        
        //============================================ If correct, is it the highest found?
        if ( CSym->at(iter)[5] > highestSym )
        {
            highestSym                                = CSym->at(iter)[5];
            bestIndex                                 = iter;
        }
    }
    
    //================================================ Found?
    if ( highestSym  > 0.0 )
    {
        settings->setRecommendedSymmetry              ( "C" );
        settings->setRecommendedFold                  ( CSym->at(bestIndex)[0] );
        ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSym->at(bestIndex) );
        
        if ( settings->detectedSymmetry.size() == 0 ) { settings->setDetectedSymmetry ( CSym->at(bestIndex) ); }
    }
    else
    {
        settings->setRecommendedSymmetry              ( "" );
        settings->setRecommendedFold                  ( 0 );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes the D symmetries and searched for the requested symmetry.
 
    This is a simple search function, which searches the symmetry results for the requested symmetry fold, and if more such
    symmetries are found, takes the one with the highest average peak height sum. If the requested fold was found, it will
    save it to the settings object, while it will set the object to fold 0 if the requested symmetry was not found (albeit
    there may be other symmetries present).
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] DSym A vector of pointers to double arrays, each array being a single Dihedral symmetry entry.
    \param[in] axes A vector to which all the axes of the requested symmetry (if any) will be saved.
 */
void ProSHADE_internal_data::ProSHADE_data::saveRequestedSymmetryD ( ProSHADE_settings* settings, std::vector< proshade_double* >* DSym, std::vector< proshade_double* >* axes )
{
    //================================================ Initialise variables
    proshade_unsign bestIndex                         = 0;
    proshade_double highestSym                        = 0.0;
    
    //================================================ Search for best fold
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( DSym->size() ); iter++ )
    {
        //============================================ Check if it is tbe correct fold
        if ( std::max ( DSym->at(iter)[0], DSym->at(iter)[6] ) != settings->requestedSymmetryFold ) { continue; }

        //============================================ If correct, is it the highest found?
        if ( ( DSym->at(iter)[5] + DSym->at(iter)[11] ) > highestSym )
        {
            highestSym                                = ( DSym->at(iter)[5] + DSym->at(iter)[11] );
            bestIndex                                 = iter;
        }
    }
    
    //================================================ Found?
    if ( highestSym  > 0.0 )
    {
        settings->setRecommendedSymmetry              ( "D" );
        settings->setRecommendedFold                  ( std::max ( DSym->at(bestIndex)[0], DSym->at(bestIndex)[6] ) );
        ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes,  DSym->at(bestIndex) );
        ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, &DSym->at(bestIndex)[6] );
        
        if ( settings->detectedSymmetry.size() == 0 )
        {
            settings->setDetectedSymmetry             ( DSym->at(bestIndex) );
            settings->setDetectedSymmetry             ( &DSym->at(bestIndex)[6] );
        }
    }
    else
    {
        settings->setRecommendedSymmetry              ( "" );
        settings->setRecommendedFold                  ( 0 );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes the group elements as rotation matrices (except for the identity element) for any detected cyclic point group.
 
    \param[in] allCSyms A vector of vectors of doubles, each array being a single Cyclic symmetry entry in a vector of all detected Cyclic symmetries.
    \param[in] grPosition An index of the C symmetry group which should have its group elements computed and returned.
    \param[out] val A vector containing vectors of 9 (rotation matrix) for each group element for the requested group, except for the identity element.
 */
std::vector<std::vector< proshade_double > > ProSHADE_internal_data::ProSHADE_data::computeGroupElementsForGroup ( std::vector<std::vector< proshade_double > >* allCSyms, proshade_unsign grPosition )
{
    //================================================ Sanity check
    if ( grPosition >= static_cast<proshade_unsign> ( allCSyms->size() ) )
    {
        std::stringstream hlpSS;
        hlpSS << "The request for group elements of group " << grPosition << " cannot be\n                    : processed, as the list of all groups does not have\n                    : group with this index.";
        throw ProSHADE_exception ( "Requested group elements for group which does not exist.", "ES00057", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Initialise variables
    std::vector<std::vector< proshade_double > > ret;
    proshade_double groupAngle                        = ( 2 * M_PI ) / static_cast<proshade_double> ( allCSyms->at(grPosition).at(0) );
    proshade_double* rotMat                           = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    
    //================================================ Generate Cn elements
    
    for ( proshade_unsign elIt = 1; elIt < static_cast<proshade_unsign> ( allCSyms->at(grPosition).at(0) ); elIt++ )
    {
        //============================================ Find the element angle
        proshade_double thisElementAngle              = static_cast<proshade_double> ( elIt ) * groupAngle;
        
        //============================================ Combine it with the group axis and get rotation matrix
        ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMat,
                                                                  allCSyms->at(grPosition).at(1),
                                                                  allCSyms->at(grPosition).at(2),
                                                                  allCSyms->at(grPosition).at(3),
                                                                  thisElementAngle );
        
        //============================================ Save the element rotation matrix to the return vector
        std::vector<proshade_double> retEl;
        for ( unsigned int matIt = 0; matIt < 9; matIt++ )
        {
            ProSHADE_internal_misc::addToDoubleVector ( &retEl, rotMat[matIt] );
        }
        ProSHADE_internal_misc::addToDoubleVectorVector ( &ret, retEl );
    }
    
    //================================================ Release memory
    delete[] rotMat;
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function computes the group elements as rotation matrices (except for the identity element) for any detected cyclic point group.
 
    \param[in] allCSyms A vector of double pointers, each array being a single Cyclic symmetry entry in a vector of all detected Cyclic symmetries.
    \param[in] grPosition An index of the C symmetry group which should have its group elements computed and returned.
    \param[out] val A vector containing vectors of 9 (rotation matrix) for each group element for the requested group, except for the identity element.
 */
std::vector<std::vector< proshade_double > > ProSHADE_internal_data::ProSHADE_data::computeGroupElementsForGroup ( std::vector< proshade_double* >* allCSyms, proshade_unsign grPosition )
{
    //================================================ Sanity check
    if ( grPosition >= static_cast<proshade_unsign> ( allCSyms->size() ) )
    {
        std::stringstream hlpSS;
        hlpSS << "The request for group elements of group " << grPosition << " cannot be\n                    : processed, as the list of all groups does not have\n                    : group with this index.";
        throw ProSHADE_exception ( "Requested group elements for group which does not exist.", "ES00057", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Initialise variables
    std::vector<std::vector< proshade_double > > ret;
    proshade_double groupAngle                        = ( 2 * M_PI ) / static_cast<proshade_double> ( allCSyms->at(grPosition)[0] );
    proshade_double* rotMat                           = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    
    //================================================ Generate Cn elements
    
    for ( proshade_unsign elIt = 1; elIt < static_cast<proshade_unsign> ( allCSyms->at(grPosition)[0] ); elIt++ )
    {
        //============================================ Find the element angle
        proshade_double thisElementAngle              = static_cast<proshade_double> ( elIt ) * groupAngle;
        
        //============================================ Combine it with the group axis and get rotation matrix
        ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMat,
                                                                  allCSyms->at(grPosition)[1],
                                                                  allCSyms->at(grPosition)[2],
                                                                  allCSyms->at(grPosition)[3],
                                                                  thisElementAngle );
        
        //============================================ Save the element rotation matrix to the return vector
        std::vector<proshade_double> retEl;
        for ( unsigned int matIt = 0; matIt < 9; matIt++ )
        {
            ProSHADE_internal_misc::addToDoubleVector ( &retEl, rotMat[matIt] );
        }
        ProSHADE_internal_misc::addToDoubleVectorVector ( &ret, retEl );
    }
    
    //================================================ Release memory
    delete[] rotMat;
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function checks that the required and obtained numbers of axes are correct, printing error if they are not.
 
    \param[in] requiredAxes Number of axes that are required by the symmetry type.
    \param[in] obtainedAxes Number of axes given.
    \param[in] groupType A string specifying for which symmetry type the group elements are to be computed.
 */
void axesToGroupTypeSanityCheck ( proshade_unsign requiredAxes, proshade_unsign obtainedAxes, std::string groupType )
{
    //================================================ Sanity check
    if ( obtainedAxes != requiredAxes )
    {
        std::stringstream hlpSS;
        hlpSS << "The supplied number of axes for group element\n                    : detection ( >" << obtainedAxes << "< ) does not match the group type ( >" << groupType << "< ).";
        throw ProSHADE_exception ( "Mismatch between supplied number of axes and\n                    : symmetry type.", "ES00059", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function adds identity matrix as the first element of the vector of vectors of doubles.
 
    \param[in] vecToPrepend Vector to which the identity element should be prepended to.
 */
void prependIdentity ( std::vector<std::vector< proshade_double > >* vecToPrepend )
{
    //================================================ Create the identity element
    std::vector< proshade_double > identity;
    ProSHADE_internal_misc::addToDoubleVector         ( &identity, 1.0 );
    ProSHADE_internal_misc::addToDoubleVector         ( &identity, 0.0 );
    ProSHADE_internal_misc::addToDoubleVector         ( &identity, 0.0 );
    ProSHADE_internal_misc::addToDoubleVector         ( &identity, 0.0 );
    ProSHADE_internal_misc::addToDoubleVector         ( &identity, 1.0 );
    ProSHADE_internal_misc::addToDoubleVector         ( &identity, 0.0 );
    ProSHADE_internal_misc::addToDoubleVector         ( &identity, 0.0 );
    ProSHADE_internal_misc::addToDoubleVector         ( &identity, 0.0 );
    ProSHADE_internal_misc::addToDoubleVector         ( &identity, 1.0 );
    
    //================================================ Prepend identity as first element
    vecToPrepend->insert                              ( vecToPrepend->begin() , identity );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function checks if the element list already contains a given matrix.
 
    \param[in] elements Vector containing all group elements.
    \param[in] elem A single element which should already be in the list.
    \param[out] elementFound A boolean value stating if the element was found int the elements list or not.
 */
bool checkElementAlreadyExists ( std::vector<std::vector< proshade_double > >* elements, std::vector< proshade_double >* elem )
{
    //================================================ Initialise variables
    bool elementFound                                 = false;
    proshade_double allowedError                      = 0.1;
    
    //================================================ For each existing element
    for ( proshade_unsign elIt = 0; elIt < static_cast<proshade_unsign> ( elements->size() ); elIt++ )
    {
        if ( ProSHADE_internal_maths::rotationMatrixSimilarity ( &elements->at(elIt), elem, allowedError ) )
        {
            elementFound                              = true;
            break;
        }
    }
    
    //================================================ Done
    return                                            ( elementFound );
    
}

/*! \brief This function checks if all group element products produce another group element.
 
    \param[in] elements Vector containing all group elements.
    \param[out] isGroup A boolean value stating if all group element products for another group element.
 */
bool checkElementsFormGroup ( std::vector<std::vector< proshade_double > >* elements )
{
    //================================================ Initialise variables
    bool isGroup                                      = true;
    
    //================================================ Multiply all group element pairs
    for ( proshade_unsign gr1 = 0; gr1 < static_cast<proshade_unsign> ( elements->size() ); gr1++ )
    {
        for ( proshade_unsign gr2 = 1; gr2 < static_cast<proshade_unsign> ( elements->size() ); gr2++ )
        {
            //======================================== Use unique pairs only
            if ( gr1 >= gr2 ) { continue; }
            
            //======================================== Multiply the two rotation matrices
            std::vector< proshade_double > product    = ProSHADE_internal_maths::multiplyGroupElementMatrices ( &elements->at(gr1), &elements->at(gr2) );
        
            //======================================== Check the group already contains the produces as an element
            if ( !checkElementAlreadyExists ( elements, &product ) )
            {
                isGroup                               = false;
                break;
            }
        }
        
        //============================================ Stop if problem was found
        if ( !isGroup ) { break; }
    }
    
    //================================================ Done
    return                                            ( isGroup );
    
}

/*! \brief This function joins two group element lists using only unique elements.
 
    \param[in] first Vector of group elements.
    \param[in] second Vector of group elements.
    \param[in] combine Should the element combinations be added as well?
    \param[out] ret A vector of group elements containing all unique elements from both input element groups.
 */
std::vector<std::vector< proshade_double > > ProSHADE_internal_data::joinElementsFromDifferentGroups ( std::vector<std::vector< proshade_double > >* first, std::vector<std::vector< proshade_double > >* second, bool combine )
{
    //================================================ Initialise variables
    std::vector< std::vector< proshade_double > > ret;
    
    //================================================ Add the first list to ret, checking for uniqueness
    for ( proshade_unsign elIt = 0; elIt < static_cast<proshade_unsign> ( first->size() ); elIt++ )
    {
        if ( !checkElementAlreadyExists( &ret, &first->at(elIt) ) )
        {
            ProSHADE_internal_misc::addToDoubleVectorVector ( &ret, first->at(elIt) );
        }
    }
    
    //================================================ Add the second list to ret, checking for uniqueness
    for ( proshade_unsign elIt = 0; elIt < static_cast<proshade_unsign> ( second->size() ); elIt++ )
    {
        if ( !checkElementAlreadyExists( &ret, &second->at(elIt) ) )
        {
            ProSHADE_internal_misc::addToDoubleVectorVector ( &ret, second->at(elIt) );
        }
    }
    
    //================================================ Multiply all combinations of first and second and check for uniqueness
    if ( combine )
    {
        for ( proshade_unsign gr1 = 0; gr1 < static_cast<proshade_unsign> ( first->size() ); gr1++ )
        {
            for ( proshade_unsign gr2 = 0; gr2 < static_cast<proshade_unsign> ( second->size() ); gr2++ )
            {
                //==================================== Multiply the two rotation matrices
                std::vector< proshade_double > product = ProSHADE_internal_maths::multiplyGroupElementMatrices ( &first->at(gr1), &second->at(gr2) );

                //==================================== Add
                if ( !checkElementAlreadyExists( &ret, &product ) )
                {
                    ProSHADE_internal_misc::addToDoubleVectorVector ( &ret, product );
                }

            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function returns the group elements as rotation matrices or any defined point group.
 
    This function generates a list of all point group elements for any group defined by a set of cyclic point groups. The set is supplied using the second
    parameter, where these need to be detected by ProSHADE first and then their index in the ProSHADE cyclic group detected list can be given here.
    
    This function can generate appropriate elementes for all ProSHADE supported point group types (i.e. C, D, T, O and I) as well as for any supplied set
    of cyclic point groups (use the groupType value of "X").
 
    Please note that the final set of point group elements will be checked for being a point group, i.e. for the fact that a product of any two members will
    be another already present member. If this condition is not met, error will be thrown. This poses some isses when the point group axes are slightly off,
    as this can lead to the point group check failing...
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] axesList A vector of ints specifying which C axes from the full list are members of the group.
    \param[in] groupType An optional string specifying for which symmetry type the group elements are to be computed. Leave empty if you want to use the supplied axes without any questions being asked.
    \param[out] val A vector containing a vector of 9 doubles (rotation matrix) for each group element for the requested group.
 */
std::vector<std::vector< proshade_double > > ProSHADE_internal_data::ProSHADE_data::getAllGroupElements ( ProSHADE_settings* settings, std::vector< proshade_unsign > axesList, std::string groupType )
{
    //================================================ Initialise variables
    std::vector<std::vector< proshade_double > > ret;
    
    //================================================ Select which symmetry type are we computing for
    if ( groupType == "C" )
    {
        //============================================ Sanity check
        axesToGroupTypeSanityCheck                    ( 1, static_cast<proshade_unsign> ( axesList.size() ), groupType );
        
        //============================================ Generate elements
        ret                                           = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(0) );
        
        //============================================ Prepend identity element
        prependIdentity                               ( &ret );

        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret ) ) { return ( ret ); }
        else
        {
            throw ProSHADE_exception ( "Computed point group elements do not form a group.", "ES00060", __FILE__, __LINE__, __func__, "The supplied cyclic groups list does not form a group and\n                    : therefore such group's elements cannot be obtained. Please\n                    : check the cyclic groups list supplied to the\n                    : getAllGroupElements() function." );
        }
    }
    else if ( groupType == "D" )
    {
        //============================================ Sanity check
        axesToGroupTypeSanityCheck                    ( 2, static_cast<proshade_unsign> ( axesList.size() ), groupType );
        
        //============================================ Generate elements for both axes
        std::vector<std::vector< proshade_double > > first  = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(0) );
        std::vector<std::vector< proshade_double > > second = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(1) );
        
        //============================================ Join the element lists
        ret                                           = joinElementsFromDifferentGroups ( &first, &second, true );
        
        //============================================ Prepend identity element
        prependIdentity                               ( &ret );
        
        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret ) ) { return ( ret ); }
        else
        {
            throw ProSHADE_exception ( "Computed point group elements do not form a group.", "ES00060", __FILE__, __LINE__, __func__, "The supplied cyclic groups list does not form a group and\n                    : therefore such group's elements cannot be obtained. Please\n                    : check the cyclic groups list supplied to the\n                    : getAllGroupElements() function." );
        }
    }
    else if ( groupType == "T" )
    {
        //============================================ Sanity check
        axesToGroupTypeSanityCheck                    ( 7, static_cast<proshade_unsign> ( axesList.size() ), groupType );
        
        //============================================ Generate elements for all four C3 axes first
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            if ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) == 3 )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(grIt) );
                
                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, false );
            }
        }
        
        //============================================ Generate elements for all three C2 axes second
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            if ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) == 2 )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(grIt) );
                
                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, false );
            }
        }
        
        //============================================ Prepend identity element
        prependIdentity                               ( &ret );
        
        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret ) ) { return ( ret ); }
        else
        {
            throw ProSHADE_exception ( "Computed point group elements do not form a group.", "ES00060", __FILE__, __LINE__, __func__, "The supplied cyclic groups list does not form a group and\n                    : therefore such group's elements cannot be obtained. Please\n                    : check the cyclic groups list supplied to the\n                    : getAllGroupElements() function." );
        }
    }
    else if ( groupType == "O" )
    {
        //============================================ Sanity check
        axesToGroupTypeSanityCheck                    ( 13, static_cast<proshade_unsign> ( axesList.size() ), groupType );
        
        //============================================ Generate elements for all three C4 axes first
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            if ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) == 4 )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(grIt) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, false );
            }
        }

        //============================================ Generate elements for all four C3 axes first
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            if ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) == 3 )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(grIt) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, false );
            }
        }
        
        //============================================ Generate elements for all six C2 axes next
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            if ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) == 2 )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(grIt) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, false );
            }
        }

        //============================================ Prepend identity element
        prependIdentity                               ( &ret );

        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret ) ) { return ( ret ); }
        else
        {
            throw ProSHADE_exception ( "Computed point group elements do not form a group.", "ES00060", __FILE__, __LINE__, __func__, "The supplied cyclic groups list does not form a group and\n                    : therefore such group's elements cannot be obtained. Please\n                    : check the cyclic groups list supplied to the\n                    : getAllGroupElements() function." );
        }
    }
    else if ( groupType == "I" )
    {
        //============================================ Sanity check
        axesToGroupTypeSanityCheck                    ( 31, static_cast<proshade_unsign> ( axesList.size() ), groupType );
        
        //============================================ Generate elements for all six C5 axes first
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C5 axis
            if ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) == 5 )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(grIt) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, false );
            }
        }
        
        //============================================ Generate elements for all ten C3 axes next
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            if ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) == 3 )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(grIt) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, false );
            }
        }
        
        //============================================ Generate elements for all fifteen C2 axes lastly
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            if ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) == 2 )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(grIt) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, false );
            }
        }
        
        //============================================ Prepend identity element
        prependIdentity                               ( &ret );

        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret ) ) { return ( ret ); }
        else
        {
            throw ProSHADE_exception ( "Computed point group elements do not form a group.", "ES00060", __FILE__, __LINE__, __func__, "The supplied cyclic groups list does not form a group and\n                    : therefore such group's elements cannot be obtained. Please\n                    : check the cyclic groups list supplied to the\n                    : getAllGroupElements() function." );
        }
    }
    else if ( groupType == "X" )
    {
        //============================================ User forced no checking for unspecified symmetry
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== Compute group elements
            std::vector<std::vector< proshade_double > > els = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, axesList.at(grIt) );
            
            //======================================== Join the elements to any already found
            ret                                       = joinElementsFromDifferentGroups ( &els, &ret, true );
        }
        
        //============================================ Prepend identity element
        prependIdentity                               ( &ret );
        
        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret ) ) { return ( ret ); }
        else
        {
            throw ProSHADE_exception ( "Computed point group elements do not form a group.", "ES00060", __FILE__, __LINE__, __func__, "The supplied cyclic groups list does not form a group and\n                    : therefore such group's elements cannot be obtained. Please\n                    : check the cyclic groups list supplied to the\n                    : getAllGroupElements() function." );
        }
    }
    else
    {
        std::stringstream hlpSS;
        hlpSS << "Unknown symmetry type: >" << groupType << "<";
        throw ProSHADE_exception ( hlpSS.str().c_str(), "ES00058", __FILE__, __LINE__, __func__, "Function getAllGroupElements was called with symmetry type\n                    : value outside of the allowed values C, D, T, O, I\n                    : or empty for using all supplied axes." );
    }
        
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function returns the length of 1D array that could hold all group elements rotation matrices for any group.
 
    Note: This is required for passing the values to python, otherwise the function has no usage.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] grIndices An array of indices to the all detected cyclic groups list, specifying which cyclic groups will form the point group for which the elements will be obtained.
    \param[in] len The length of the grIndices array.
    \param[in] groupType A string specifying which group type the elements will be computed for - allowed values are C, D, T, O, I and X.
    \param[out] val The minimal length of a 1D array that can hold all the group elements rotation matrices.
 */
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getAllGroupElementsLength ( ProSHADE_settings* settings, int* grIndices, int len, std::string groupType )
{
    //================================================ Convert array to vector
    std::vector< proshade_unsign > axes;
    for ( int iter = 0; iter < len; iter++ )
    {
        ProSHADE_internal_misc::addToUnsignVector     ( &axes, grIndices[iter] );
    }
    
    //================================================ Get the elements
    std::vector<std::vector< proshade_double > > groupElements = this->getAllGroupElements ( settings, axes, groupType );
    
    //================================================ Return their size
    return                                            ( static_cast<proshade_unsign> ( groupElements.size() * 9 ) );
    
}

/*! \brief This function computes all point group elements for a group formed by any number of cyclic point groups.
 
    This function is a Python wrapper function for the getAllGroupElements() function. It takes the indices of the cyclic point groups which form the point group for which
    the group elements are required as a Python list and fills a Python compatible array with the results.
 
    \warning This function has specific signature for SWIG processing into proshade Python module, please use the getAllGroupElements() funtion
    for C++ access.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] grIndices A list of indices of cyclic point groups which should be combined to create the point group whose group elements are required.
    \param[in] len The length of the grIndices array.
    \param[in] groupType The type of the group. Allowed values are C, D, T, O, I and X.
    \param[in] allGroupElement This is the array to which the results will be saved into.
    \param[in] ln2 The lenght of the allGroupElement array.
 */
void ProSHADE_internal_data::ProSHADE_data::getAllGroupElementsPython ( ProSHADE_settings* settings, int* grIndices, int len, std::string groupType, double* allGroupElement, int ln2 )
{
    //================================================ Convert array to vector
    std::vector< proshade_unsign > axes;
    for ( int iter = 0; iter < len; iter++ )
    {
        ProSHADE_internal_misc::addToUnsignVector     ( &axes, grIndices[iter] );
    }
    
    //================================================ Get the elements in C++ format
    std::vector<std::vector< proshade_double > > groupElements = this->getAllGroupElements ( settings, axes, groupType );
    
    //================================================ Re-save them in Python format
    for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( groupElements.size() ); grIt++ )
    {
        for ( proshade_unsign matIt = 0; matIt < 9; matIt++ )
        {
            allGroupElement[(grIt*9)+matIt]           = groupElements.at(grIt).at(matIt);
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function returns the length of 1D array that could hold all the cyclic group elements rotation matrices.
 
    Note: This is required for passing the values to python, otherwise the function has no usage.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] grPosition An index of the C symmetry group which should have its group elements computed and returned.
    \param[out] val The minimal length of a 1D array that can hold all the group elements rotation matrices.
 */
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getCGroupElementsLength ( ProSHADE_settings* settings, proshade_unsign grPosition )
{
    //================================================ Sanity check
    if ( grPosition >= static_cast<proshade_unsign> ( settings->allDetectedCAxes.size() ) )
    {
        std::stringstream hlpSS;
        hlpSS << "The request for group elements of group " << grPosition << " cannot be\n                    : processed, as the list of all groups does not have\n                    : group with this index.";
        throw ProSHADE_exception ( "Requested group elements for group which does not exist.", "ES00057", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Done
    return                                            ( static_cast<proshade_unsign> ( settings->allDetectedCAxes.at(grPosition).at(0) - 1 ) * 9 );
    
}

/*! \brief This function computes the group elements rotation matrices (except for the identity element) for requested group and fills the supplied 1D array with them.
 
    \warning The identity element is ignored by this function.
    \warning This function has specific signature for SWIG processing into proshade Python module, please use the computeGroupElementsForGroup() funtion
    for C++ access.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] grPosition An index of the C symmetry group which should have its group elements computed and returned.
 */
void ProSHADE_internal_data::ProSHADE_data::getCGroupElementsPython ( ProSHADE_settings* settings, double* groupElements, int len, proshade_unsign grPosition )
{
    //================================================ Sanity check
    if ( grPosition >= static_cast<proshade_unsign> ( settings->allDetectedCAxes.size() ) )
    {
        std::stringstream hlpSS;
        hlpSS << "The request for group elements of group " << grPosition << " cannot be\n                    : processed, as the list of all groups does not have\n                    : group with this index.";
        throw ProSHADE_exception ( "Requested group elements for group which does not exist.", "ES00057", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
        
    //================================================ Get the matrices
    std::vector<std::vector< proshade_double > > grElements = this->computeGroupElementsForGroup ( &settings->allDetectedCAxes, grPosition );
    
    //================================================ Copy to Python array
    for ( proshade_unsign elIt = 0; elIt < static_cast<proshade_unsign> ( grElements.size() ); elIt++ )
    {
        for ( proshade_unsign matIt = 0; matIt < 9; matIt++ )
        {
            groupElements[(elIt*9)+matIt]             = grElements.at(elIt).at(matIt);
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function copies the internal map into the supplied pointer, which it also allocates.
 
    This function is provided so that the user can provide a pointer and have it allocated and filled with the map values.
 
    \param[in] saveTo A pointer where the internal map should be deep copied into.
    \param[in] verbose How loud the run should be?
 */
void ProSHADE_internal_data::ProSHADE_data::deepCopyMap ( proshade_double*& saveTo, proshade_unsign verbose )
{
    //================================================ Sanity check
    if ( saveTo != NULL )
    {
        ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! The deep copy pointer is not set to NULL. Cannot proceed and returning unmodified pointer.", "WB00040" );
        return ;
    }
    
    //================================================ Allocate the memory
    saveTo                                            = new proshade_double[this->xDimIndices * this->yDimIndices * this->zDimIndices];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation ( saveTo, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy internal map to the new pointer
    for ( proshade_unsign iter = 0; iter < ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ )
    {
        saveTo[iter]                                  = this->internalMap[iter];
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes prints the report for symmetry detection.
 
    This is a very simple function which provides the basic textual output for the symmetry detection task.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection reporting.
 */
void ProSHADE_internal_data::ProSHADE_data::reportSymmetryResults ( ProSHADE_settings* settings )
{
    //================================================ Improve this!
    if ( settings->recommendedSymmetryType == "" )
    {
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, "Did not detect any symmetry!" );
    }
    else
    {
        std::stringstream ssHlp;
        ssHlp << std::endl << "Detected " << settings->recommendedSymmetryType << " symmetry with fold " << settings->recommendedSymmetryFold << " .";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        
        if ( settings->detectedSymmetry.size() > 0 )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << "  Fold       X           Y          Z           Angle        Height";
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( settings->detectedSymmetry.size() ); symIt++ )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << std::showpos << std::fixed << std::setprecision(0) << "   " << settings->detectedSymmetry.at(symIt)[0] << std::setprecision(5) << "     " << settings->detectedSymmetry.at(symIt)[1] << "   " << settings->detectedSymmetry.at(symIt)[2] << "   " << settings->detectedSymmetry.at(symIt)[3] << "     " << settings->detectedSymmetry.at(symIt)[4] << "      " << settings->detectedSymmetry.at(symIt)[5];
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        
        std::stringstream hlpSS3;
        ssHlp.clear(); ssHlp.str ( "" );
        hlpSS3 << std::endl << "However, since the selection of the recommended symmetry needs improvement, here is a list of all detected C symmetries:";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, hlpSS3.str() );
        
        if ( settings->allDetectedCAxes.size() > 0 )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << "  Fold       X           Y          Z           Angle        Height";
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( settings->allDetectedCAxes.size() ); symIt++ )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << std::showpos << std::fixed << std::setprecision(0) << "   " << settings->allDetectedCAxes.at(symIt)[0] << std::setprecision(5) << "     " << settings->allDetectedCAxes.at(symIt)[1] << "   " << settings->allDetectedCAxes.at(symIt)[2] << "   " << settings->allDetectedCAxes.at(symIt)[3] << "     " << settings->allDetectedCAxes.at(symIt)[4] << "      " << settings->allDetectedCAxes.at(symIt)[5];
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        
        hlpSS3.clear(); hlpSS3.str ( "" );
        hlpSS3 << std::endl << "Also, for the same reason, here is a list of all detected D symmetries:";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, hlpSS3.str() );
        
        if ( settings->allDetectedDAxes.size() > 0 )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << "  Fold       X           Y          Z           Angle        Height";
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( settings->allDetectedDAxes.size() ); symIt++ )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << std::showpos << std::fixed << std::setprecision(0) << "   " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(0))[0] << std::setprecision(5) << "     " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(0))[1] << "   " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(0))[2] << "   " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(0))[3] << "     " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(0))[4] << "      " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(0))[5];
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
            
            for ( proshade_unsign axIt = 1; axIt < static_cast<proshade_unsign> ( settings->allDetectedDAxes.at(symIt).size() ); axIt++ )
            {
                ssHlp.clear(); ssHlp.str ( "" );
                ssHlp << std::showpos << std::fixed << std::setprecision(0) << "   " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(axIt))[0] << std::setprecision(5) << "     " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(axIt))[1] << "   " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(axIt))[2] << "   " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(axIt))[3] << "     " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(axIt))[4] << "      " << settings->allDetectedCAxes.at(settings->allDetectedDAxes.at(symIt).at(axIt))[5];
                ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
            }
            
            ssHlp.clear(); ssHlp.str ( "" );
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function sets the original structure representation value to a set of dedicated variables.
 
 
    This function is called once a structure hes been read in and initially processed by ProSHADE. Therefore, the values it saves are already after axis inversion, COM centering, etc., but before any map rotation or map padding. The reason why these
    values need to be saved is because in some circumstances (such as PDB file wiritng) these are required and the normal variables holding these may be over-written by the map manipulation process.
*/
void ProSHADE_internal_data::ProSHADE_data::setOriginalMapValues ()
{
    //================================================ Compute and save the COM
    this->findMapCOM                                  ( );
    this->originalMapXCom                             = this->xCom;
    this->originalMapYCom                             = this->yCom;
    this->originalMapZCom                             = this->zCom;
    
    //================================================ Save the original sizes
    this->xDimSizeOriginal                            = this->xDimSize;
    this->yDimSizeOriginal                            = this->yDimSize;
    this->zDimSizeOriginal                            = this->zDimSize;
    
    //================================================ Save the original index counts
    this->xDimIndicesOriginal                         = this->xDimIndices;
    this->yDimIndicesOriginal                         = this->yDimIndices;
    this->zDimIndicesOriginal                         = this->zDimIndices;
    
    //================================================ Save the original axis origins
    this->xAxisOriginOriginal                         = this->xAxisOrigin;
    this->yAxisOriginOriginal                         = this->yAxisOrigin;
    this->zAxisOriginOriginal                         = this->zAxisOrigin;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the centre of mass of the internal map representation.

    This function simply computes the centre of mass for the given ProSHADE_data object map in the "real space" (i.e. the space that counts Angstroms from the bottom left further corner). These are then saved into the ProSHADE_data object.
*/
void ProSHADE_internal_data::ProSHADE_data::findMapCOM ( )
{
    //================================================ Initialise variables
    this->xCom                                        = 0.0;
    this->yCom                                        = 0.0;
    this->zCom                                        = 0.0;
    proshade_double totNonZeroPoints                  = 0.0;
    proshade_signed mapIt                             = 0;
    
    //================================================ Compute COM from 0 ; 0 ; 0
    for ( proshade_signed xIt = 0; xIt < static_cast<proshade_signed> ( this->xDimIndices ); xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < static_cast<proshade_signed> ( this->yDimIndices ); yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < static_cast<proshade_signed> ( this->zDimIndices ); zIt++ )
            {
                //==================================== Find map index
                mapIt                                 = zIt  + this->zDimIndices * ( yIt  + this->yDimIndices * xIt  );
                
                //==================================== Use only positive density
                if ( this->internalMap[mapIt] <= 0.0 ) { continue; }
                
                //==================================== Compute Index COM
                this->xCom                           += this->internalMap[mapIt] * static_cast<proshade_double> ( xIt + 1 );
                this->yCom                           += this->internalMap[mapIt] * static_cast<proshade_double> ( yIt + 1 );
                this->zCom                           += this->internalMap[mapIt] * static_cast<proshade_double> ( zIt + 1 );
                totNonZeroPoints                     += this->internalMap[mapIt];
            }
        }
    }
    
    this->xCom                                       /= totNonZeroPoints;
    this->yCom                                       /= totNonZeroPoints;
    this->zCom                                       /= totNonZeroPoints;
    
    //================================================ Convert to real world
    this->xCom                                        = ( (this->xCom-1) - this->xAxisOrigin ) * ( static_cast<proshade_double> ( this->xDimIndices - 1 ) / this->xDimSize );
    this->yCom                                        = ( (this->yCom-1) - this->yAxisOrigin ) * ( static_cast<proshade_double> ( this->yDimIndices - 1 ) / this->yDimSize );
    this->zCom                                        = ( (this->zCom-1) - this->zAxisOrigin ) * ( static_cast<proshade_double> ( this->zDimIndices - 1 ) / this->zDimSize );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function returns the number of spheres which contain the whole object.
 
    \param[out] X The total number of spheres to which the structure is mapped.
 */
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getMaxSpheres ( )
{
    //================================================ Return the value
    return                                            ( this->noSpheres );
}

/*! \brief This function returns the internal map representation value of a particular array position.
 
    \param[in] pos The position in the map array, of which the value should be returned.
    \param[out] X The internal map representation value at position pos.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::getMapValue ( proshade_unsign pos )
{
    //================================================ Return the value
    return                                            ( this->internalMap[pos] );
}

/*! \brief This function returns the maximum band value for the object.
 
    \param[out] X The largest number of bands used in any shell of the object.
 */
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getMaxBand ( )
{
    //================================================ Return the value
    return                                            ( this->maxShellBand );
}

/*! \brief This function allows access to the priva internal RRP matrices.
 
    \param[out] X The value of the internal private RRP matrix for the given indices.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::getRRPValue ( proshade_unsign band, proshade_unsign sh1, proshade_unsign sh2 )
{
    //================================================ Return the value
    return                                            ( this->rrpMatrices[band][sh1][sh2] );
}

/*! \brief This function checks if particular shell has a  particular band.
 
    This function is useful for the progressive shell mapping, where it may not be clear in one part of the code
    whether a particular shell does or does not have a particular band value. Therefore, this function allows simple
    check.
 
    \param[in] shell The index (number) of the shell for which the check should be done.
    \param[in] bandVal The band value which should be sought for the shell.
    \param[out] X True if the shell has the band, false otherwise.
 */
bool ProSHADE_internal_data::ProSHADE_data::shellBandExists ( proshade_unsign shell, proshade_unsign bandVal )
{
    if ( this->spheres[shell]->getLocalBandwidth( ) >= bandVal )
    {
        return                                        ( true );
    }
    else
    {
        return                                        ( false );
    }
}

/*! \brief This function removes phase from the map.
 
    This function is called when the phase information needs to be removed from the internal map representation. It
    does the forward Fourier transform, removes the phase from the Fourier coefficients and then the inverse Fourier
    transform, thus resulting with the Patterson map. It does write over the original map.
 
    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
 */
void ProSHADE_internal_data::ProSHADE_data::removePhaseInormation ( ProSHADE_settings* settings )
{
    //================================================ Report function start
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Centering map onto its COM." );
    
    //================================================ Copy map for processing
    fftw_complex* mapCoeffs                           = new fftw_complex[this->xDimIndices * this->yDimIndices * this->zDimIndices];
    fftw_complex* pattersonMap                        = new fftw_complex[this->xDimIndices * this->yDimIndices * this->zDimIndices];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( mapCoeffs, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( pattersonMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy data to mask
    for ( proshade_unsign iter = 0; iter < (this->xDimIndices * this->yDimIndices * this->zDimIndices); iter++ )
    {
        pattersonMap[iter][0]                         = this->internalMap[iter];
        pattersonMap[iter][1]                         = 0.0;
    }
    
    //================================================ Prepare FFTW plans
    fftw_plan forward                                 = fftw_plan_dft_3d ( this->xDimIndices, this->yDimIndices, this->zDimIndices,
                                                                           pattersonMap, mapCoeffs, FFTW_FORWARD,  FFTW_ESTIMATE );
    fftw_plan inverse                                 = fftw_plan_dft_3d ( this->xDimIndices, this->yDimIndices, this->zDimIndices,
                                                                           mapCoeffs, pattersonMap, FFTW_BACKWARD, FFTW_ESTIMATE );
    
    //================================================ Run forward Fourier
    fftw_execute                                      ( forward );
    
    //================================================ Remove the phase
    ProSHADE_internal_mapManip::removeMapPhase        ( mapCoeffs, this->xDimIndices, this->yDimIndices, this->zDimIndices );
    
    //================================================ Run inverse Fourier
    fftw_execute                                      ( inverse );
    
    //================================================ Save the results
    proshade_signed mapIt, patIt, patX, patY, patZ;
    for ( proshade_signed xIt = 0; xIt < static_cast<proshade_signed> ( this->xDimIndices ); xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < static_cast<proshade_signed> ( this->yDimIndices ); yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < static_cast<proshade_signed> ( this->zDimIndices ); zIt++ )
            {
                //==================================== Centre patterson map
                patX = xIt - ( static_cast<proshade_signed> ( this->xDimIndices ) / 2 ); if ( patX < 0 ) { patX += this->xDimIndices; }
                patY = yIt - ( static_cast<proshade_signed> ( this->yDimIndices ) / 2 ); if ( patY < 0 ) { patY += this->yDimIndices; }
                patZ = zIt - ( static_cast<proshade_signed> ( this->zDimIndices ) / 2 ); if ( patZ < 0 ) { patZ += this->zDimIndices; }
                
                //==================================== Find indices
                mapIt                                 = zIt  + this->zDimIndices * ( yIt  + this->yDimIndices * xIt  );
                patIt                                 = patZ + this->zDimIndices * ( patY + this->yDimIndices * patX );
                
                //==================================== Copy
                this->internalMap[mapIt]              = pattersonMap[patIt][0];
            }
        }
    }
    
    //================================================ Release memory
    delete[] pattersonMap;
    delete[] mapCoeffs;
    
    //================================================ Delete FFTW plans
    fftw_destroy_plan                                 ( forward );
    fftw_destroy_plan                                 ( inverse );
    
    //================================================ Report function completion
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, "Phase information removed." );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allows access to the private internal real spherical harmonics values.
 
    \param[out] X Pointer to the value of the internal private spherical harmonics real value of the given index.
 */
proshade_double* ProSHADE_internal_data::ProSHADE_data::getRealSphHarmValue ( proshade_unsign band, proshade_unsign order, proshade_unsign shell )
{
    //================================================ Done
    return                                            ( &this->sphericalHarmonics[shell][seanindex ( static_cast<proshade_signed> ( order ) - static_cast<proshade_signed> ( band ),
                                                                                                     band,
                                                                                                     this->spheres[shell]->getLocalBandwidth() )][0] );
    
}

/*! \brief This function allows access to the private internal imaginary spherical harmonics values.
 
    \param[out] X Pointer to the value of the internal private spherical harmonics imaginary value of the given index.
 */
proshade_double* ProSHADE_internal_data::ProSHADE_data::getImagSphHarmValue ( proshade_unsign band, proshade_unsign order, proshade_unsign shell )
{
    //================================================ Done
    return                                            ( &this->sphericalHarmonics[shell][seanindex ( static_cast<proshade_signed> ( order ) - static_cast<proshade_signed> ( band ),
                                                                                                     band,
                                                                                                     this->spheres[shell]->getLocalBandwidth() )][1] );
    
}

/*! \brief This function allows access to the radius of any particular sphere.
 
    \param[out] X The distance of the requested sphere to the centre of the coordinates.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::getAnySphereRadius ( proshade_unsign shell )
{
    //================================================ Done
    return                                            ( this->spheres[shell]->getShellRadius() );
    
}

/*! \brief This function allows access to the integration weight for the object.
 
    \param[out] X The integration weight for the object or 0.0 if not yet computed.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::getIntegrationWeight  ( void )
{
    //================================================ Done
    return                                            ( this->integrationWeight );
    
}

/*! \brief This function allows access to the bandwidth of a particular shell.
 
    \param[in] shell The index of the shell for which the bandwidth is required.
    \param[out] X The bandwidth of the requested shell.
 */
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getShellBandwidth ( proshade_unsign shell )
{
    //================================================ Done
    return                                            ( this->spheres[shell]->getLocalBandwidth ( ) );
    
}

/*! \brief This function allows access to sphere positions.
 
    \param[in] shell The index of the sphere for which the position (radius) is to be obtained.
    \param[out] X The radius of the sphere with index shell.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::getSpherePosValue ( proshade_unsign shell )
{
    //================================================ Done
    return                                            ( this->spherePos.at(shell) );
    
}

/*! \brief This function allows access to E matrix for a particular band.
 
    \param[in] band The band for which the E matrix subset order * order should be returned.
    \param[out] X Pointer to pointer of complex matrix with dimensions order * order of the E matrices.
 */
proshade_complex** ProSHADE_internal_data::ProSHADE_data::getEMatrixByBand ( proshade_unsign band )
{
    //================================================ Done
    return                                            ( this->eMatrices[band] );
    
}

/*! \brief This function allows access to E matrix by knowing the band, order1 and order2 indices.
 
    \param[in] band The band for which the E matrix value should be returned.
    \param[in] order1 The first order for which the E matrix value should be returned.
    \param[in] order2 The second order for which the E matrix value should be returned.
    \param[in] valueReal The proshade_double number pointer to where the real part of the value will be saved.
    \param[in] valueImag The proshade_double number pointer to where the imaginary part of the value will be saved.
 */
void ProSHADE_internal_data::ProSHADE_data::getEMatrixValue ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_double* valueReal, proshade_double* valueImag )
{
    //================================================ Set pointer
   *valueReal                                         = this->eMatrices[band][order1][order2][0];
   *valueImag                                         = this->eMatrices[band][order1][order2][1];
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allows access to the inverse SO(3) coefficients array.
 
    \param[out] X The internal inverse SO(3) coefficients array variable.
 */
proshade_complex* ProSHADE_internal_data::ProSHADE_data::getInvSO3Coeffs ( )
{
    //================================================ Done
    return                                            ( this->so3CoeffsInverse );
    
}

/*! \brief This function allows access to the SO(3) coefficients array.
 
    \param[out] X The internal SO(3) coefficients array variable.
 */
proshade_complex* ProSHADE_internal_data::ProSHADE_data::getSO3Coeffs ( )
{
    //================================================ Done
    return                                            ( this->so3Coeffs );
    
}

/*! \brief This function allows access to the maximum band for the comparison.
 
    \param[out] X The bandwidth used for this comparison.
 */
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getComparisonBand ( )
{
    //================================================ Done
    return                                            ( this->maxCompBand );
    
}

/*! \brief This function allows access to the Wigner D matrix by knowing the band, order1 and order2 indices.
 
    \param[in] band The band for which the Wigner D matrix value should be returned.
    \param[in] order1 The first order for which the Wigner D matrix value should be returned.
    \param[in] order2 The second order for which the Wigner D matrix value should be returned.
    \param[in] valueReal The proshade_double number pointer to where the real part of the value will be saved.
    \param[in] valueImag The proshade_double number pointer to where the imaginary part of the value will be saved.
 */
void ProSHADE_internal_data::ProSHADE_data::getWignerMatrixValue ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_double* valueReal, proshade_double* valueImag )
{
    //================================================ Set pointer
    *valueReal                                        = this->wignerMatrices[band][order1][order2][0];
    *valueImag                                        = this->wignerMatrices[band][order1][order2][1];
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allows access to the map size in angstroms along the X axis.
 
    \param[out] xDimSize The size of the internal map in angstroms along the X axis.
 */
proshade_single ProSHADE_internal_data::ProSHADE_data::getXDimSize ( )
{
    //================================================ Return the requested value
    return                                            ( this->xDimSize );
}

/*! \brief This function allows access to the map size in angstroms along the Y axis.
 
    \param[out] xDimSize The size of the internal map in angstroms along the Y axis.
 */
proshade_single ProSHADE_internal_data::ProSHADE_data::getYDimSize ( )
{
    //================================================ Return the requested value
    return                                            ( this->yDimSize );
}

/*! \brief This function allows access to the map size in angstroms along the Z axis.
 
    \param[out] xDimSize The size of the internal map in angstroms along the Z axis.
 */
proshade_single ProSHADE_internal_data::ProSHADE_data::getZDimSize ( )
{
    //================================================ Return the requested value
    return                                            ( this->zDimSize );
}

/*! \brief This function allows access to the map size in indices along the X axis.
 
    \param[out] xDimSize The size of the internal map in indices along the X axis.
 */
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getXDim ( )
{
    //================================================ Return the requested value
    return                                            ( this->xDimIndices );
}

/*! \brief This function allows access to the map size in indices along the Y axis.
 
    \param[out] xDimSize The size of the internal map in indices along the Y axis.
 */
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getYDim ( )
{
    //================================================ Return the requested value
    return                                            ( this->yDimIndices );
}

/*! \brief This function allows access to the map size in indices along the Z axis.
 
    \param[out] xDimSize The size of the internal map in indices along the Z axis.
 */
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getZDim ( )
{
    //================================================ Return the requested value
    return                                            ( this->zDimIndices );
}

/*! \brief This function allows access to the map start along the X axis.
 
    \param[out] xFrom Pointer to the starting index along the X axis.
 */
proshade_signed* ProSHADE_internal_data::ProSHADE_data::getXFromPtr ( )
{
    //================================================ Return the requested value
    return                                            ( &this->xFrom );
}

/*! \brief This function allows access to the map start along the Y axis.
 
    \param[out] yFrom Pointer to the starting index along the Y axis.
 */
proshade_signed* ProSHADE_internal_data::ProSHADE_data::getYFromPtr ( )
{
    //================================================ Return the requested value
    return                                            ( &this->yFrom );
}

/*! \brief This function allows access to the map start along the Z axis.
 
    \param[out] zFrom Pointer to the starting index along the Z axis.
 */
proshade_signed* ProSHADE_internal_data::ProSHADE_data::getZFromPtr ( )
{
    //================================================ Return the requested value
    return                                            ( &this->zFrom );
}

/*! \brief This function allows access to the map last position along the X axis.
 
    \param[out] xFrom Pointer to the final index along the X axis.
 */
proshade_signed* ProSHADE_internal_data::ProSHADE_data::getXToPtr ( )
{
    //================================================ Return the requested value
    return                                            ( &this->xTo );
}

/*! \brief This function allows access to the map last position along the Y axis.
 
    \param[out] yFrom Pointer to the final index along the Y axis.
 */
proshade_signed* ProSHADE_internal_data::ProSHADE_data::getYToPtr ( )
{
    //================================================ Return the requested value
    return                                            ( &this->yTo );
}

/*! \brief This function allows access to the map last position along the Z axis.
 
    \param[out] zFrom Pointer to the final index along the Z axis.
 */
proshade_signed* ProSHADE_internal_data::ProSHADE_data::getZToPtr ( )
{
    //================================================ Return the requested value
    return                                            ( &this->zTo );
}

/*! \brief This function allows access to the map X axis origin value.
 
    \param[out] xAxisOrigin The value of X axis origin for the map.
 */
proshade_signed* ProSHADE_internal_data::ProSHADE_data::getXAxisOrigin ( )
{
    //================================================ Return the requested value
    return                                            ( &this->xAxisOrigin );
}

/*! \brief This function allows access to the map Y axis origin value.
 
    \param[out] yAxisOrigin The value of Y axis origin for the map.
 */
proshade_signed* ProSHADE_internal_data::ProSHADE_data::getYAxisOrigin ( )
{
    //================================================ Return the requested value
    return                                            ( &this->yAxisOrigin );
}

/*! \brief This function allows access to the map Z axis origin value.
 
    \param[out] zAxisOrigin The value of Z axis origin for the map.
 */
proshade_signed* ProSHADE_internal_data::ProSHADE_data::getZAxisOrigin ( )
{
    //================================================ Return the requested value
    return                                            ( &this->zAxisOrigin );
}

/*! \brief This function allows access to the first map array value address.

    \param[out] internalMap Pointer to the first position in the internal map array.
*/
proshade_double*& ProSHADE_internal_data::ProSHADE_data::getInternalMap ( )
{
    //================================================ Return the requested value
    return                                            ( this->internalMap );
}

/*! \brief This function allows access to the translation function through a pointer.

    \param[out] translationMap Pointer to the first position in the translation function map array.
*/
proshade_complex* ProSHADE_internal_data::ProSHADE_data::getTranslationFnPointer ( void )
{
    //================================================ Return the requested value
    return                                            ( this->translationMap );
}

/*! \brief This function allows setting the integration weight for the object.
 
    \param[in] intW The integration weight to be set for this object.
 */
void ProSHADE_internal_data::ProSHADE_data::setIntegrationWeight ( proshade_double intW )
{
    //================================================ Mutate
    this->integrationWeight                           = intW;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allows setting the cumulative integration weight for the object.
 
    \param[in] intW The integration weight to be added to the current value for this object.
 */
void ProSHADE_internal_data::ProSHADE_data::setIntegrationWeightCumul ( proshade_double intW )
{
    //================================================ Mutate
    this->integrationWeight                          += intW;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allows setting the E matrix value.
 
    \param[in] band The band indice of the E matrix to which the value should be assigned.
    \param[in] order1 The first order indice of the E matrix to which the value should be assigned.
    \param[in] order2 The second order indice of the E matrix to which the value should be assigned.
    \param[in] val The value which should be saved.
 */
void ProSHADE_internal_data::ProSHADE_data::setEMatrixValue ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_complex val )
{
    //================================================ Mutate
    this->eMatrices[band][order1][order2][0]          = val[0];
    this->eMatrices[band][order1][order2][1]          = val[1];
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allows normalising the E matrix value.
 
    \param[in] band The band indice of the E matrix to which the value should be assigned.
    \param[in] order1 The first order indice of the E matrix to which the value should be assigned.
    \param[in] order2 The second order indice of the E matrix to which the value should be assigned.
    \param[in] normF The value by which the original E matrix value will be divided to normalise it.
 */
void ProSHADE_internal_data::ProSHADE_data::normaliseEMatrixValue ( proshade_unsign band, proshade_unsign order1, proshade_unsign order2, proshade_double normF )
{
    //================================================ Mutate
    this->eMatrices[band][order1][order2][0]         /= normF;
    this->eMatrices[band][order1][order2][1]         /= normF;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allows setting the SOFT coefficient values using array position and value.
 
    \param[in] position The 1D array position at which the new value should be saved.
    \param[in] val Complex value to be saved into the array.
 */
void ProSHADE_internal_data::ProSHADE_data::setSO3CoeffValue ( proshade_unsign position, proshade_complex val )
{
    //================================================ Mutate
    this->so3Coeffs[position][0]                      = val[0];
    this->so3Coeffs[position][1]                      = val[1];
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allows setting the Wigner D matrix value by its band, order1 and order2 co-ordinate.
 
    \param[in] val proshade_complex value of the Wigner D matrix at position band, order1, order2.
    \param[in] band The band of the Wigner D matrix value.
    \param[in] order1 The first order of the Wigner D matrix value.
    \param[in] order2 The second order of the Wigner D matrix value.
 */
void ProSHADE_internal_data::ProSHADE_data::setWignerMatrixValue ( proshade_complex val, proshade_unsign band, proshade_unsign order1, proshade_unsign order2 )
{
    //================================================ Mutate
    this->wignerMatrices[band][order1][order2][0]     = val[0];
    this->wignerMatrices[band][order1][order2][1]     = val[1];
    
    //================================================ Done
    return ;
    
}

/*! \brief This function simplyfies Python access to maps - it returns how long the map array is, so that the getMap function take take this value as input.
 
    \param[out] val The size of the map array.
 */
int ProSHADE_internal_data::ProSHADE_data::getMapArraySizePython ( void )
{
    //================================================ Done
    return                                            ( static_cast<int> ( this->xDimIndices * this->yDimIndices * this->zDimIndices ) );
    
}

/*! \brief This function takes an array and fills it with the internal map density data. This is necessary for Python Numpy arrays output.

    \param[in] mapArrayPython The array to which the map data will be saved into.
    \param[in] len The length of the array, expected to be given by the getMapArraySizePython() function.
*/
void ProSHADE_internal_data::ProSHADE_data::getMapPython ( double *mapArrayPython, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        mapArrayPython[iter]                          = static_cast<double> ( this->internalMap[iter] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes an array and re-writes the internal map value with the argument map values.

    \param[in] mapChangedInPython The array from which the map data will be read.
    \param[in] len The length of the array, expected to be given by the getMapArraySizePython() function.
*/
void ProSHADE_internal_data::ProSHADE_data::setMapPython ( double *mapChangedInPython, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        this->internalMap[iter]                       = static_cast<proshade_double> ( mapChangedInPython[iter] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes an array, deletes old map, allocates a new map and fills the internal map values with the argument map values.

    \param[in] mapChangedInPython The array from which the map data will be read.
    \param[in] len The length of the array, expected to be given by the getMapArraySizePython() function.
*/
void ProSHADE_internal_data::ProSHADE_data::setNewMapPython ( double *mapChangedInPython, int len )
{
    //================================================ Delete the old map data
    delete[] this->internalMap;
    this->internalMap                                 = NULL;
    
    //================================================ Allocate new map space
    this->internalMap                                 = new proshade_double [this->xDimIndices * this->yDimIndices * this->zDimIndices];
    ProSHADE_internal_misc::checkMemoryAllocation     ( this->internalMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        this->internalMap[iter]                       = static_cast<proshade_double> ( mapChangedInPython[iter] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills a given array with the real parts of the spherical harmonics values for a given shell.

    \param[in] shellNo The index of the shell for which the spherical harmonics are being retried.
    \param[in] verbose How loud this run should be.
    \param[in] sphericalHarmsReal The array to which the data will be loaded into.
    \param[in] len The length of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getRealSphericalHarmonicsForShell ( proshade_unsign shellNo, proshade_signed verbose, double *sphericalHarmsReal, int len )
{
    //================================================ Sanity check
    if ( shellNo >= this->noSpheres )
    {
        ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! Attempted to access shell index outside of the shell range.", "WP00043" );
        return ;
    }
    
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        sphericalHarmsReal[iter]                      = static_cast<double> ( this->sphericalHarmonics[shellNo][iter][0] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills a given array with the imaginary parts of the spherical harmonics values for a given shell.

    \param[in] shellNo The index of the shell for which the spherical harmonics are being retried.
    \param[in] verbose How loud this run should be.
    \param[in] sphericalHarmonics The array to which the data will be loaded into.
    \param[in] len The length of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getImagSphericalHarmonicsForShell ( proshade_unsign shellNo, proshade_signed verbose, double *sphericalHarmsImag, int len )
{
    //================================================ Sanity check
    if ( shellNo >= this->noSpheres )
    {
        ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! Attempted to access shell index outside of the shell range.", "WP00043" );
        return ;
    }
    
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        sphericalHarmsImag[iter]                      = static_cast<double> ( this->sphericalHarmonics[shellNo][iter][1] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function gets the spherical harmonics array length for a particular shell.

    \param[in] shellNo The index of the shell for which the spherical harmonics are being retried.
    \param[in] verbose How loud this run should be.
    \param[out] val Length of spherical harmonics array or zero if error occured.
*/
int ProSHADE_internal_data::ProSHADE_data::getSphericalHarmonicsLenForShell  ( proshade_unsign shellNo, proshade_signed verbose )
{
    //================================================ Sanity check
    if ( shellNo >= this->noSpheres )
    {
        ProSHADE_internal_messages::printWarningMessage ( verbose, "!!! ProSHADE WARNING !!! Attempted to access shell index outside of the shell range.", "WP00043" );
        return                                        ( 0 );
    }
    
    //================================================ Return the value
    return                                            ( static_cast<int> ( (this->spheres[shellNo]->getLocalBandwidth() * 2) * (this->spheres[shellNo]->getLocalBandwidth() * 2) ) );
}

/*! \brief This function gets the spherical harmonics array index for a particular spherical harmonics band and order.

    \param[in] order The order for which the spherical harmonics value index is requested.
    \param[in] band The band for which the spherical harmonics value index is requested.
    \param[in] shell The shell for which the spherical harmonics value index is requested.
    \param[out] val Index value of the spherical harmonics value.
*/
int ProSHADE_internal_data::ProSHADE_data::sphericalHarmonicsIndex ( proshade_signed order, proshade_signed band, proshade_signed shell )
{
    //================================================ Return the value
    return                                            ( static_cast<int> ( seanindex ( order, band, this->spheres[shell]->getLocalBandwidth() ) ) );
}

/*! \brief This function fills the input array with the real E matrix values for particular band and order1 (l as opposed to l').

    \param[in] band The band for which the real E matrix values are requested.
    \param[in] order The order for which the real E matrix values are requested.
    \param[in] eMatsLMReal The array to which the values will be written into.
    \param[in] len The lenght of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getRealEMatrixValuesForLM ( proshade_signed band, proshade_signed order1, double *eMatsLMReal, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        eMatsLMReal[iter]                             = static_cast<double> ( this->eMatrices[band][order1][iter][0] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills the input array with the imaginary E matrix values for particular band and order1 (l as opposed to l').

    \param[in] band The band for which the imaginary E matrix values are requested.
    \param[in] order The order for which the imaginary E matrix values are requested.
    \param[in] eMatsLMImag The array to which the values will be written into.
    \param[in] len The lenght of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getImagEMatrixValuesForLM ( proshade_signed band, proshade_signed order1, double *eMatsLMImag, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        eMatsLMImag[iter]                             = static_cast<double> ( this->eMatrices[band][order1][iter][1] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills the input array with the real SO(3) coefficient values.

    \param[in] so3CoefsReal The array to which the values will be written into.
    \param[in] len The lenght of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getRealSO3Coeffs ( double *so3CoefsReal, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        so3CoefsReal[iter]                            = static_cast<double> ( this->so3Coeffs[iter][0] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills the input array with the imaginary SO(3) coefficient values.

    \param[in] so3CoefsImag The array to which the values will be written into.
    \param[in] len The lenght of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getImagSO3Coeffs ( double *so3CoefsImag, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        so3CoefsImag[iter]                            = static_cast<double> ( this->so3Coeffs[iter][1] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function gets the SO(3) coefficients array index for a particular so(3) band, order1 and order2 position.
 
    It should be noted that this function assumes that the orders are in the format -l < 0 < l and NOT 0 < 2l + 1.

    \param[in] order1 The first order for which the SO(3) value index is requested.
    \param[in] order2 The second order for which the SO(3) value index is requested.
    \param[in] band The band for which the SO(3) value index is requested.
    \param[out] val Index position of the SO(3) value.
*/
int ProSHADE_internal_data::ProSHADE_data::so3CoeffsArrayIndex ( proshade_signed order1, proshade_signed order2, proshade_signed band )
{
    //================================================ Return the value
    return                                            ( static_cast<int> ( so3CoefLoc ( order1, order2, band, this->getMaxBand() ) ) );
}

/*! \brief This function fills the input array with the real rotation function values.

    \param[in] rotFunReal The array to which the values will be written into.
    \param[in] len The lenght of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getRealRotFunction ( double *rotFunReal, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        rotFunReal[iter]                              = static_cast<double> ( this->so3CoeffsInverse[iter][0] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills the input array with the imaginary rotation function values.

    \param[in] rotFunImag The array to which the values will be written into.
    \param[in] len The lenght of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getImagRotFunction ( double *rotFunImag, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        rotFunImag[iter]                              = static_cast<double> ( this->so3CoeffsInverse[iter][1] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills the input array with the real translation function values.

    \param[in] trsFunReal The array to which the values will be written into.
    \param[in] len The lenght of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getRealTranslationFunction ( double *trsFunReal, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        trsFunReal[iter]                              = static_cast<double> ( this->translationMap[iter][0] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function fills the input array with the imaginary translation function values.

    \param[in] trsFunImag The array to which the values will be written into.
    \param[in] len The lenght of the array.
*/
void ProSHADE_internal_data::ProSHADE_data::getImagTranslationFunction ( double *trsFunImag, int len )
{
    //================================================ Save the data into the output array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        trsFunImag[iter]                              = static_cast<double> ( this->translationMap[iter][1] );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes rotation function indices, converts them to Euler angles and these to rotation matrix, which it then returns.

    \param[in] aI The index along the Euler alpha dimension.
    \param[in] bI The index along the Euler beta dimension.
    \param[in] gI The index along the Euler gamma dimension.
    \param[in] rotMat The array to which the rotation matrix will be written into.
    \param[in] len The lenght of the array (must be 9).
*/
void ProSHADE_internal_data::ProSHADE_data::getRotMatrixFromRotFunInds ( proshade_signed aI, proshade_signed bI, proshade_signed gI, double *rotMat, int len )
{
    //================================================ Get Euler angles
    proshade_double eA, eB, eG;
    ProSHADE_internal_maths::getEulerZXZFromSOFTPosition ( this->getMaxBand(), aI, bI, gI, &eA, &eB, &eG );
    
    //================================================ Prepare internal rotation matrix memory
    proshade_double* rMat                             = NULL;
    rMat                                              = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rMat, __FILE__, __LINE__, __func__ );
    
    //================================================ Convert to rotation matrix
    ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( eA, eB, eG, rMat );
    
    //================================================ Copy to output
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( len ); iter++ )
    {
        rotMat[iter]                                  = static_cast<double> ( rMat[iter] );
    }
    
    //================================================ Release internal memory
    delete[] rMat;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function simply returns the detected recommended symmetry type.

    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
*/
std::string ProSHADE_internal_data::ProSHADE_data::getRecommendedSymmetryType ( ProSHADE_settings* settings )
{
    //================================================ Return the value
    return                                            ( settings->recommendedSymmetryType );
    
}

/*! \brief This function simply returns the detected recommended symmetry fold.

    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
*/
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getRecommendedSymmetryFold ( ProSHADE_settings* settings )
{
    //================================================ Return the value
    return                                            ( settings->recommendedSymmetryFold );
    
}

/*! \brief This function returns the number of detected recommended symmetry axes.

    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
    \param[out] val The length of the recommended symmetry axes vector.
*/
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getNoRecommendedSymmetryAxes ( ProSHADE_settings* settings )
{
    //================================================ Return the value
    return                                            ( static_cast<proshade_unsign> ( settings->detectedSymmetry.size() ) );
}

/*! \brief This function returns the length of 1D array which would contain all detected axes info.

    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
    \param[out] val The required length of 1D array to hold all detected axes info.
*/
proshade_unsign ProSHADE_internal_data::ProSHADE_data::getAllSymsOneArrayLength ( ProSHADE_settings* settings )
{
    //================================================ Return the value
    return                                            ( static_cast<proshade_unsign> ( settings->allDetectedCAxes.size() * 6 ) );
}

/*! \brief This function returns a single symmetry axis as a vector of strings from the recommended symmetry axes list.

    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] axisNo The index of the axis to be returned.
    \param[out] val A vector of strings containing the symmetry axis fold, x, y, z axis element, angle and peak height in this order.
*/
std::vector< std::string > ProSHADE_internal_data::ProSHADE_data::getSymmetryAxis ( ProSHADE_settings* settings, proshade_unsign axisNo )
{
    //================================================ Sanity checks
    if ( static_cast<proshade_unsign> ( settings->detectedSymmetry.size() ) <= axisNo )
    {
        ProSHADE_internal_messages::printWarningMessage ( settings->verbose, "!!! ProSHADE WARNING !!! Requested symmetry index does not exist. Returning empty vector.", "WS00039" );
        return                                        ( std::vector< std::string > ( ) );
    }
    
    //================================================ Initialise local variables
    std::vector< std::string > ret;
    
    //================================================ Input the axis data as strings
    std::stringstream ssHlp;
    ssHlp << settings->detectedSymmetry.at(axisNo)[0];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
    
    ssHlp << settings->detectedSymmetry.at(axisNo)[1];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
    
    ssHlp << settings->detectedSymmetry.at(axisNo)[2];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
    
    ssHlp << settings->detectedSymmetry.at(axisNo)[3];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
    
    ssHlp << settings->detectedSymmetry.at(axisNo)[4];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
    
    ssHlp << settings->detectedSymmetry.at(axisNo)[5];
    ProSHADE_internal_misc::addToStringVector         ( &ret, ssHlp.str() );
    ssHlp.str                                         ( "" );
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function writes out the rotated map, co-ordinates and transformation JSON file.

    This function takes basically all the results of the overlay mode and appropriately applies them to write out the
    moved density map, if possible the moved co-ordinates and also the overlay operations listing JSON file.
 
    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
    \param[in] trsX The optimal x-axis translation value.
    \param[in] trsY The optimal y-axis translation value.
    \param[in] trsZ The optimal z-axis translation value.
    \param[in] eulA The Euler alpha angle value, by which the moving structure is to be rotated by.
    \param[in] eulB The Euler beta angle value, by which the moving structure is to be rotated by.
    \param[in] eulG The Euler gamma angle value, by which the moving structure is to be rotated by.
    \param[in] rotCentre The rotation centre position as determined by the computeOverlayTranslations function.
    \param[in] ultimateTranslation The final translation as determined by the computeOverlayTranslations function.
*/
void ProSHADE_internal_data::ProSHADE_data::writeOutOverlayFiles ( ProSHADE_settings* settings, proshade_double trsX, proshade_double trsY, proshade_double trsZ, proshade_double eulA, proshade_double eulB, proshade_double eulG, std::vector< proshade_double >* rotCentre, std::vector< proshade_double >* ultimateTranslation )
{
    //================================================ Write out rotated map
    std::stringstream fNameHlp;
    fNameHlp << settings->overlayStructureName << ".map";
    this->writeMap                                    ( fNameHlp.str() );
    
    //================================================ Write out rotated co-ordinates if possible
    if ( ProSHADE_internal_io::isFilePDB ( this->fileName ) )
    {
        fNameHlp.str("");
        fNameHlp << settings->overlayStructureName << ".pdb";
        this->writePdb                                ( fNameHlp.str(), eulA, eulB, eulG, trsX, trsY, trsZ, settings->firstModelOnly );
    }
    
    //================================================ Write out the json file with the results
    ProSHADE_internal_io::writeRotationTranslationJSON ( -rotCentre->at(0), -rotCentre->at(1), -rotCentre->at(2),
                                                         eulA, eulB, eulG,
                                                         ultimateTranslation->at(0), ultimateTranslation->at(1), ultimateTranslation->at(2),
                                                         this->comMovX, this->comMovY, this->comMovZ, settings->rotTrsJSONFile );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function sets the correct translation values for the overlay mode.

    This function takes the translation as computed by the translation function (in the last three input variables) and proceeds  to
    compute and save the initial centre of rotation as well as the ultimate translation to be done after rotation to obtain the optimal
    overlay of the moving structure over the static structure.
    
    \param[in] rcX Pointer to where to save the rotation centre position along the X-axis in Angstroms.
    \param[in] rcY Pointer to where to save the rotation centre position along the Y-axis in Angstroms.
    \param[in] rcZ Pointer to where to save the rotation centre position along the Z-axis in Angstroms.
    \param[in] transX Pointer to where to save the translation to be done along the X-axis in Angstroms. This variable should already have the computed translation from the translation map.
    \param[in] transY Pointer to where to save the translation to be done along the Y-axis in Angstroms. This variable should already have the computed translation from the translation map.
    \param[in] transZ Pointer to where to save the translation to be done along the Z-axis in Angstroms. This variable should already have the computed translation from the translation map.
*/
void ProSHADE_internal_data::ProSHADE_data::computeOverlayTranslations ( proshade_double* rcX, proshade_double* rcY, proshade_double* rcZ, proshade_double* transX, proshade_double* transY, proshade_double* transZ )
{
    //================================================ Write out the json file with the results
    if ( ProSHADE_internal_io::isFilePDB ( this->fileName ) )
    {
        //============================================ If PDB, we already have these
       *rcX                                           = this->originalPdbRotCenX;
       *rcY                                           = this->originalPdbRotCenY;
       *rcZ                                           = this->originalPdbRotCenZ;
         
       *transX                                        = *transX + this->originalPdbRotCenX;
       *transY                                        = *transY + this->originalPdbRotCenY;
       *transZ                                        = *transZ + this->originalPdbRotCenZ;
    }
    else
    {
        //============================================ Compute the rotation centre for the co-ordinates
        proshade_double xRotPos                       = ( ( static_cast<proshade_double> ( this->xDimIndicesOriginal / 2 ) - this->xAxisOriginOriginal ) *
                                                          ( static_cast<proshade_double> ( this->xDimIndicesOriginal - 1 ) / this->xDimSizeOriginal ) ) -
                                                          (this->comMovX);
        proshade_double yRotPos                       = ( ( static_cast<proshade_double> ( this->yDimIndicesOriginal / 2 ) - this->yAxisOriginOriginal ) *
                                                          ( static_cast<proshade_double> ( this->yDimIndicesOriginal - 1 ) / this->yDimSizeOriginal ) ) -
                                                          (this->comMovY);
        proshade_double zRotPos                       = ( ( static_cast<proshade_double> ( this->zDimIndicesOriginal / 2 ) - this->zAxisOriginOriginal ) *
                                                          ( static_cast<proshade_double> ( this->zDimIndicesOriginal - 1 ) / this->zDimSizeOriginal ) ) -
                                                          (this->comMovZ);
        
        //============================================ And save
        *rcX                                           = xRotPos;
        *rcY                                           = yRotPos;
        *rcZ                                           = zRotPos;
          
        *transX                                        = *transX + xRotPos;
        *transY                                        = *transY + yRotPos;
        *transZ                                        = *transZ + zRotPos;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function reports the results of the overlay mode.

    \param[in] settings ProSHADE_settings object specifying the details of how the computations should be done.
    \param[in] rotationCentre Pointer to vector for saving the position of the centre of rotation about which the rotation is to be done.
    \param[in] mapBoxMovement Pointer to vector for saving the sum of all translations done internally by ProSHADE to this input map.
    \param[in] eulerAngles Pointer to vector where the three Euler angles will be saved into.
    \param[in] finalTranslation Pointer to a vector where the translation required to move structure from origin to optimal overlay with static structure will be saved into.
*/
void ProSHADE_internal_data::ProSHADE_data::reportOverlayResults ( ProSHADE_settings* settings, std::vector < proshade_double >* rotationCentre, std::vector< proshade_double >* mapBoxMovement, std::vector < proshade_double >* eulerAngles, std::vector < proshade_double >* finalTranslation )
{
    //================================================ Empty line
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, "" );
    
    //================================================ Write out rotation centre translation results
    std::stringstream rotCen; rotCen << std::setprecision (3) << std::showpos << "The rotation centre to origin translation vector is: " << -rotationCentre->at(0) << "     " << -rotationCentre->at(1) << "     " << -rotationCentre->at(2);
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, rotCen.str() );
    
    //================================================ Write out internal map translation results
    std::stringstream mapBox; mapBox << std::setprecision (3) << std::showpos << "The within box internal map translation vector is  : " << mapBoxMovement->at(0) << "     " << mapBoxMovement->at(1) << "     " << mapBoxMovement->at(2);
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, mapBox.str() );
    
    //================================================ Write out rotation matrix about origin
    proshade_double* rotMat                           = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( eulerAngles->at(0), eulerAngles->at(1), eulerAngles->at(2), rotMat );
    
    std::stringstream rotMatSS;
    rotMatSS << std::setprecision (3) << std::showpos << "The rotation matrix about origin is                : " << rotMat[0] << "     " << rotMat[1] << "     " << rotMat[2] << std::endl;
    rotMatSS << std::setprecision (3) << std::showpos << "                                                   : " << rotMat[3] << "     " << rotMat[4] << "     " << rotMat[5] << std::endl;
    rotMatSS << std::setprecision (3) << std::showpos << "                                                   : " << rotMat[6] << "     " << rotMat[7] << "     " << rotMat[8];
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, rotMatSS.str() );
    
    delete[] rotMat;
    
    //================================================ Write out origin to overlay translation results
    std::stringstream finTrs; finTrs << std::setprecision (3) << std::showpos << "The origin to overlay translation vector is        : " << finalTranslation->at(0) << "     " << finalTranslation->at(1) << "     " << finalTranslation->at(2);
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, finTrs.str() );
    
    //================================================ Done
    return ;
    
}
