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
    \version   0.7.5.5
    \date      MAY 2021
 */

//==================================================== ProSHADE
#include "ProSHADE_data.hpp"

//==================================================== Do not use the following flags for the included files - this causes a lot of warnings that have nothing to do with ProSHADE
#if defined ( __GNUC__ )
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wpedantic"
    #pragma GCC diagnostic ignored "-Wshadow"
    #pragma GCC diagnostic ignored "-Wall"
    #pragma GCC diagnostic ignored "-Wextra"
    #pragma GCC diagnostic ignored "-Wdouble-promotion"
    #pragma GCC diagnostic ignored "-Wconversion"
#endif

//==================================================== Do not use the following flags for the included files - this causes a lot of warnings that have nothing to do with ProSHADE
#if defined ( __clang__ )
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wpedantic"
    #pragma clang diagnostic ignored "-Wshadow"
    #pragma clang diagnostic ignored "-Wall"
    #pragma clang diagnostic ignored "-Wextra"
    #pragma clang diagnostic ignored "-Wdouble-promotion"
    #pragma clang diagnostic ignored "-Weverything"
#endif

//==================================================== Remove MSVC C4996 Warnings caused by Gemmi code
#if defined ( _MSC_VER )
    #pragma warning ( disable:4996 )
#endif

//==================================================== Gemmi PDB output - this cannot be with the rest of includes for some stb_sprintf library related reasons ...
#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_pdb.hpp>

//==================================================== Enable MSVC C4996 Warnings for the rest of the code
#if defined ( _MSC_VER )
    #pragma warning ( default:4996 )
#endif

//==================================================== Now the flags can be restored and used as per the CMakeLists.txt file.
#if defined ( __GNUC__ )
    #pragma GCC diagnostic pop
#endif

//==================================================== Now the flags can be restored and used as per the CMakeLists.txt file.
#if defined ( __clang__ )
    #pragma clang diagnostic pop
#endif

//==================================================== Local functions prototypes
void        axesToGroupTypeSanityCheck                ( proshade_unsign requiredAxes, proshade_unsign obtainedAxes, std::string groupType );
bool        checkElementAlreadyExists                 ( std::vector<std::vector< proshade_double > >* elements, std::vector< proshade_double >* elem, proshade_double matrixTolerance );
bool        checkElementsFormGroup                    ( std::vector<std::vector< proshade_double > >* elements, proshade_double matrixTolerance );
bool        sortProSHADESymmetryByFSC                 ( proshade_double* a, proshade_double* b );

/*! \brief Constructor for getting empty ProSHADE_data class.
 
    This constructor creates an empty data structure which can later be filled with data and used to process
    these data further.
 
    \param[in] settings ProSHADE_settings object specifying what should be done.
    \param[out] X Empty data object with deault values.
 */
ProSHADE_internal_data::ProSHADE_data::ProSHADE_data ( )
{
    //================================================ Initialise variables
    // ... Variables regarding input file
    this->fileName                                    = "";
    this->fileType                                    = ProSHADE_internal_io::UNKNOWN;
    
    // ... Variables regarding map
    this->internalMap                                 = nullptr;
    
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
    this->mapMovFromsChangeX                          = 0.0;
    this->mapMovFromsChangeY                          = 0.0;
    this->mapMovFromsChangeZ                          = 0.0;
    this->mapCOMProcessChangeX                        = 0.0;
    this->mapCOMProcessChangeY                        = 0.0;
    this->mapCOMProcessChangeZ                        = 0.0;
    
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
    this->spheres                                     = nullptr;
    this->sphericalHarmonics                          = nullptr;
    this->rotSphericalHarmonics                       = nullptr;
    this->maxShellBand                                = 0;
    
    // ... Variables regarding shape distance computations
    this->rrpMatrices                                 = nullptr;
    this->eMatrices                                   = nullptr;
    this->so3Coeffs                                   = nullptr;
    this->so3CoeffsInverse                            = nullptr;
    this->wignerMatrices                              = nullptr;
    this->integrationWeight                           = 0.0;
    this->maxCompBand                                 = 0;
    this->translationMap                              = nullptr;
    
    
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
ProSHADE_internal_data::ProSHADE_data::ProSHADE_data ( std::string strName, double *mapVals, int len, proshade_single xDmSz, proshade_single yDmSz, proshade_single zDmSz, proshade_unsign xDmInd, proshade_unsign yDmInd, proshade_unsign zDmInd, proshade_signed xFr, proshade_signed yFr, proshade_signed zFr, proshade_signed xT, proshade_signed yT, proshade_signed zT, proshade_unsign inputO )
{
    //================================================ Initialise variables
    // ... Variables regarding input file
    this->fileName                                    = strName;
    this->fileType                                    = ProSHADE_internal_io::MAP;
    
    // ... Variables regarding map
    this->internalMap                                 = nullptr;
    
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
    this->mapMovFromsChangeX                          = 0.0;
    this->mapMovFromsChangeY                          = 0.0;
    this->mapMovFromsChangeZ                          = 0.0;
    this->mapCOMProcessChangeX                        = 0.0;
    this->mapCOMProcessChangeY                        = 0.0;
    this->mapCOMProcessChangeZ                        = 0.0;
    
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
    this->spheres                                     = nullptr;
    this->sphericalHarmonics                          = nullptr;
    this->rotSphericalHarmonics                       = nullptr;
    this->maxShellBand                                = 0;
    
    // ... Variables regarding shape distance computations
    this->rrpMatrices                                 = nullptr;
    this->eMatrices                                   = nullptr;
    this->so3Coeffs                                   = nullptr;
    this->so3CoeffsInverse                            = nullptr;
    this->wignerMatrices                              = nullptr;
    this->integrationWeight                           = 0.0;
    this->maxCompBand                                 = 0;
    this->translationMap                              = nullptr;
        
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
    
    //================================================ Release memory (it was allocated by the PyBind11 lambda function and needs to be released)
    delete[] mapVals;
    
    //================================================ Done
    
}

/*! \brief Destructor for the ProSHADE_data class.
 
    This destructor is responsible for releasing all memory used by the data storing object
 
    \param[out] X N/A.
 */
ProSHADE_internal_data::ProSHADE_data::~ProSHADE_data ( )
{
    //================================================ Release the internal map
    if ( this->internalMap != nullptr )
    {
        delete[] this->internalMap;
    }
    
    //================================================ Release the sphere mapping
    if ( this->spheres != nullptr )
    {
        for ( proshade_unsign iter = 0; iter < this->noSpheres; iter++ )
        {
            if ( this->spheres[iter] != nullptr )
            {
                delete this->spheres[iter];
                this->spheres[iter]                   = nullptr;
            }
        }
        delete[] this->spheres;
    }
    
    //================================================ Release the spherical harmonics
    if ( this->sphericalHarmonics != nullptr )
    {
        for ( proshade_unsign iter = 0; iter < this->noSpheres; iter++ )
        {
            if ( this->sphericalHarmonics[iter] != nullptr )
            {
                delete[] this->sphericalHarmonics[iter];
                this->sphericalHarmonics[iter]        = nullptr;
            }
        }
        delete[] this->sphericalHarmonics;
    }
    
    //================================================ Release the rotated spherical harmonics
    if ( this->rotSphericalHarmonics != nullptr )
    {
        for ( proshade_unsign iter = 0; iter < this->noSpheres; iter++ )
        {
            if ( this->rotSphericalHarmonics[iter] != nullptr )
            {
                delete[] this->rotSphericalHarmonics[iter];
                this->rotSphericalHarmonics[iter]     = nullptr;
            }
        }
        delete[] this->rotSphericalHarmonics;
    }
    
    //================================================ Release the RRP matrices (pre-computation for the energy levels descriptor)
    if ( this->rrpMatrices != nullptr )
    {
        for ( proshade_unsign bwIt = 0; bwIt < this->maxShellBand; bwIt++ )
        {
            if ( this->rrpMatrices[bwIt] != nullptr )
            {
                for ( proshade_unsign shIt = 0; shIt < this->noSpheres; shIt++ )
                {
                    if ( this->rrpMatrices[bwIt][shIt] != nullptr )
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
    if ( this->eMatrices != nullptr )
    {
        for ( proshade_unsign bandIter = 0; bandIter < this->maxCompBand; bandIter++ )
        {
            if ( this->eMatrices[bandIter] != nullptr )
            {
                for ( proshade_unsign band2Iter = 0; band2Iter < static_cast<proshade_unsign> ( ( bandIter * 2 ) + 1 ); band2Iter++ )
                {
                    if ( this->eMatrices[bandIter][band2Iter] != nullptr )
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
    if ( this->so3Coeffs != nullptr )
    {
        delete[] this->so3Coeffs;
    }
    if ( this->so3CoeffsInverse != nullptr )
    {
        delete[] this->so3CoeffsInverse;
    }
    
    //================================================ Release Wigner matrices
    if ( this->wignerMatrices != nullptr )
    {
        for ( proshade_unsign bandIter = 1; bandIter < this->maxCompBand; bandIter++ )
        {
            if ( this->wignerMatrices[bandIter] != nullptr )
            {
                for ( proshade_unsign order1Iter = 0; order1Iter < ( (bandIter * 2) + 1 ); order1Iter++ )
                {
                    if ( this->wignerMatrices[bandIter][order1Iter] != nullptr )
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
    if ( this->translationMap != nullptr )
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
    
    //================================================ If mask is supplied and the correct task is used
    if ( ( settings->appliedMaskFileName != "" ) && ( ( settings->task == MapManip ) || ( settings->task == Symmetry ) ) )
    {
        //================================================ Report progress
        std::stringstream hlpSS;
        hlpSS << "Reading mask file " << settings->appliedMaskFileName;
        ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 2, hlpSS.str() );
        
        //============================================ Open the mask
        gemmi::Ccp4<float> mask;
        mask.read_ccp4                                ( gemmi::MaybeGzipped ( settings->appliedMaskFileName.c_str() ) );
        
        //============================================ Convert to XYZ and create complete mask, if need be
        mask.setup                                    ( gemmi::GridSetup::ReorderOnly, 0 );
        
        //============================================ Read in the rest of the mask file header
        proshade_unsign xDI, yDI, zDI, xAOR, yAOR, zAOR, xGI, yGI, zGI;
        proshade_single xDS, yDS, zDS, aA, bA, cA;
        proshade_signed xF, yF, zF, xAO, yAO, zAO;
        ProSHADE_internal_io::readInMapHeader         ( &mask,
                                                        &xDI,  &yDI,  &zDI,
                                                        &xDS,  &yDS,  &zDS,
                                                        &aA,   &bA,   &cA,
                                                        &xF,   &yF,   &zF,
                                                        &xAO,  &yAO,  &zAO,
                                                        &xAOR, &yAOR, &zAOR,
                                                        &xGI,  &yGI,  &zGI );
        
        //============================================ Sanity check
        if ( ( this->xDimIndices != xDI ) || ( this->yDimIndices != yDI ) || ( this->zDimIndices != zDI ) )
        {
            throw ProSHADE_exception                  ( "The supplied map mask has different dimensions than the\n                    : density map.", "EM00065", __FILE__, __LINE__, __func__, "Most likely the mask is not the correct mask for this map,\n                    : as it has different dimensions from the density map.\n                    : Please review that the supplied map and mask form a pair." );
        }
        
        //============================================ Save the mask values to ProSHADE variable
        proshade_double* internalMask = nullptr;
        ProSHADE_internal_io::readInMapData           ( &mask, internalMask, xDI, yDI, zDI, xAOR, yAOR, zAOR );
        
        //============================================ Apply the mask to the map
        for ( size_t iter = 0; iter < static_cast< size_t > ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ ) { this->internalMap[iter] *= internalMask[iter]; }
        
        //============================================ Release the memory
        delete[] internalMask;
        
        //============================================ Report progress
        ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 3, "Mask read in and applied successfully." );
    }
    
    //================================================ Remove negative values if so required
    if ( settings->removeNegativeDensity ) { for ( size_t iter = 0; iter < static_cast< size_t > ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ ) { if ( this->internalMap[iter] < 0.0 ) { this->internalMap[iter] = 0.0; } } }
    
    //================================================ Set resolution if need be
    if ( settings->requestedResolution < 0.0f )
    {
        settings->setResolution                       ( std::min ( static_cast<proshade_single> ( this->xDimSize ) / static_cast<proshade_single> ( this->xDimIndices ),
                                                        std::min ( static_cast<proshade_single> ( this->yDimSize ) / static_cast<proshade_single> ( this->yDimIndices ),
                                                                   static_cast<proshade_single> ( this->zDimSize ) / static_cast<proshade_single> ( this->zDimIndices ) ) ) * 2.0f );
    }
    
    //================================================ Set iterators from and to
    this->figureIndexStartStop                        ( );
    
    //================================================ If specific resolution is requested, make sure the map has it
    if ( settings->changeMapResolution || settings->changeMapResolutionTriLinear )
    {
        //============================================ Before re-sampling sampling rate
        proshade_single xSampRate                     = this->xDimSize / static_cast< proshade_single > ( this->xTo - this->xFrom );
        proshade_single ySampRate                     = this->yDimSize / static_cast< proshade_single > ( this->yTo - this->yFrom );
        proshade_single zSampRate                     = this->zDimSize / static_cast< proshade_single > ( this->zTo - this->zFrom );
        
        //============================================ Bofore re-sampling first index position
        proshade_single xStartPosBefore               = static_cast< proshade_single > ( this->xFrom ) * xSampRate;
        proshade_single yStartPosBefore               = static_cast< proshade_single > ( this->yFrom ) * ySampRate;
        proshade_single zStartPosBefore               = static_cast< proshade_single > ( this->zFrom ) * zSampRate;
        
        //============================================ Find COM before map re-sampling
        proshade_double xMapCOMPreReSampl = 0.0, yMapCOMPreReSampl = 0.0, zMapCOMPreReSampl = 0.0;
        ProSHADE_internal_mapManip::findMAPCOMValues  ( this->internalMap, &xMapCOMPreReSampl, &yMapCOMPreReSampl, &zMapCOMPreReSampl, this->xDimSize, this->yDimSize, this->zDimSize, this->xFrom, this->xTo, this->yFrom, this->yTo, this->zFrom, this->zTo );
        
        //============================================ Re-sample map
        this->reSampleMap                             ( settings );
    
        //============================================ After re-sampling sampling rate
        xSampRate                                     = this->xDimSize / static_cast< proshade_single > ( this->xTo - this->xFrom );
        ySampRate                                     = this->yDimSize / static_cast< proshade_single > ( this->yTo - this->yFrom );
        zSampRate                                     = this->zDimSize / static_cast< proshade_single > ( this->zTo - this->zFrom );
        
        //============================================ After re-sampling first index position
        proshade_single xStartPosAfter                = static_cast< proshade_single > ( this->xFrom ) * xSampRate;
        proshade_single yStartPosAfter                = static_cast< proshade_single > ( this->yFrom ) * ySampRate;
        proshade_single zStartPosAfter                = static_cast< proshade_single > ( this->zFrom ) * zSampRate;
        
        //============================================ Translate by change in corners to make the boxes as similarly placed as possible
        proshade_single xMov                          = static_cast< proshade_single > ( xStartPosAfter - xStartPosBefore );
        proshade_single yMov                          = static_cast< proshade_single > ( yStartPosAfter - yStartPosBefore );
        proshade_single zMov                          = static_cast< proshade_single > ( zStartPosAfter - zStartPosBefore );
        ProSHADE_internal_mapManip::moveMapByIndices  ( &xMov, &yMov, &zMov, this->xDimSize, this->yDimSize, this->zDimSize,
                                                        &this->xFrom, &this->xTo, &this->yFrom, &this->yTo, &this->zFrom, &this->zTo,
                                                        &this->xAxisOrigin, &this->yAxisOrigin, &this->zAxisOrigin );
        
        //============================================ Find COM after map re-sampling and corner move
        proshade_double xMapCOMPostReSampl = 0.0, yMapCOMPostReSampl = 0.0, zMapCOMPostReSampl = 0.0;
        ProSHADE_internal_mapManip::findMAPCOMValues  ( this->internalMap, &xMapCOMPostReSampl, &yMapCOMPostReSampl, &zMapCOMPostReSampl, this->xDimSize, this->yDimSize, this->zDimSize, this->xFrom, this->xTo, this->yFrom, this->yTo, this->zFrom, this->zTo );
        
        //============================================ Match the COMs to get as close position of the re-sampled structure to the original as possible
        ProSHADE_internal_mapManip::moveMapByFourier  ( this->internalMap,
                                                        static_cast< proshade_single > ( xMapCOMPreReSampl - xMapCOMPostReSampl ),
                                                        static_cast< proshade_single > ( yMapCOMPreReSampl - yMapCOMPostReSampl ),
                                                        static_cast< proshade_single > ( zMapCOMPreReSampl - zMapCOMPostReSampl ),
                                                        this->xDimSize, this->yDimSize, this->zDimSize,
                                                        static_cast< proshade_signed > ( this->xDimIndices ),
                                                        static_cast< proshade_signed > ( this->yDimIndices ),
                                                        static_cast< proshade_signed > ( this->zDimIndices ) );
    }
    
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
    
    //================================================ Compute and save the COM
    this->findMapCOM                                  ( );
    this->originalMapXCom                             = this->xCom;
    this->originalMapYCom                             = this->yCom;
    this->originalMapZCom                             = this->zCom;
    
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
    if ( settings->requestedResolution < 0.0f )
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
    proshade_single xMov                              = static_cast< proshade_single > ( 20.0f - xF );
    proshade_single yMov                              = static_cast< proshade_single > ( 20.0f - yF );
    proshade_single zMov                              = static_cast< proshade_single > ( 20.0f - zF );
    ProSHADE_internal_mapManip::movePDBForMapCalc     ( &pdbFile, xMov, yMov, zMov, settings->firstModelOnly );
    
    //================================================ Set the angstrom sizes
    this->xDimSize                                    = static_cast< proshade_single > ( xT - xF + 40.0f );
    this->yDimSize                                    = static_cast< proshade_single > ( yT - yF + 40.0f );
    this->zDimSize                                    = static_cast< proshade_single > ( zT - zF + 40.0f );

    //================================================ Generate map from nicely placed atoms (cell size will be range + 40)
    ProSHADE_internal_mapManip::generateMapFromPDB    ( pdbFile, this->internalMap, settings->requestedResolution, this->xDimSize, this->yDimSize, this->zDimSize, &this->xTo, &this->yTo, &this->zTo, settings->forceP1, settings->firstModelOnly );
    
    //================================================ Remove negative values if so required
    if ( settings->removeNegativeDensity ) { for ( size_t iter = 0; iter < static_cast< size_t > ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ ) { if ( this->internalMap[iter] < 0.0 ) { this->internalMap[iter] = 0.0; } } }
    
    //================================================ Set the internal variables to correct values
    this->setPDBMapValues                             ( );
    
    //================================================ Compute reverse movement based on COMs. If there is more than 1 models, simply moving back the xyzMov is not enough.
    proshade_double xCOMMap, yCOMMap, zCOMMap;
    ProSHADE_internal_mapManip::findMAPCOMValues      ( this->internalMap, &xCOMMap, &yCOMMap, &zCOMMap,
                                                        this->xDimSize, this->yDimSize, this->zDimSize,
                                                        this->xFrom, this->xTo, this->yFrom, this->yTo, this->zFrom, this->zTo );
    
    if ( pdbFile.models.size() > 1 )
    {
        xMov                                          = static_cast< proshade_single > ( xCOMMap - xCOMPdb );
        yMov                                          = static_cast< proshade_single > ( yCOMMap - yCOMPdb );
        zMov                                          = static_cast< proshade_single > ( zCOMMap - zCOMPdb );
    }
    
    //================================================ Move map back to the original PDB location
    ProSHADE_internal_mapManip::moveMapByIndices      ( &xMov, &yMov, &zMov, this->xDimSize, this->yDimSize, this->zDimSize,
                                                        &this->xFrom, &this->xTo, &this->yFrom, &this->yTo, &this->zFrom, &this->zTo,
                                                        &this->xAxisOrigin, &this->yAxisOrigin, &this->zAxisOrigin );
    ProSHADE_internal_mapManip::moveMapByFourier      ( this->internalMap, xMov, yMov, zMov, this->xDimSize, this->yDimSize, this->zDimSize,
                                                        static_cast< proshade_signed > ( this->xDimIndices ), static_cast< proshade_signed > ( this->yDimIndices ),
                                                        static_cast< proshade_signed > ( this->zDimIndices ) );
    
    //================================================ If specific resolution is requested, make sure the map has it
    if ( settings->changeMapResolution || settings->changeMapResolutionTriLinear )
    {
        //============================================ Before re-sampling sampling rate
        proshade_single xSampRate                     = this->xDimSize / static_cast< proshade_single > ( this->xTo - this->xFrom );
        proshade_single ySampRate                     = this->yDimSize / static_cast< proshade_single > ( this->yTo - this->yFrom );
        proshade_single zSampRate                     = this->zDimSize / static_cast< proshade_single > ( this->zTo - this->zFrom );
        
        //============================================ Bofore re-sampling first index position
        proshade_single xStartPosBefore               = static_cast< proshade_single > ( this->xFrom ) * xSampRate;
        proshade_single yStartPosBefore               = static_cast< proshade_single > ( this->yFrom ) * ySampRate;
        proshade_single zStartPosBefore               = static_cast< proshade_single > ( this->zFrom ) * zSampRate;
        
        //============================================ Find COM before map re-sampling
        proshade_double xMapCOMPreReSampl = 0.0, yMapCOMPreReSampl = 0.0, zMapCOMPreReSampl = 0.0;
        ProSHADE_internal_mapManip::findMAPCOMValues  ( this->internalMap, &xMapCOMPreReSampl, &yMapCOMPreReSampl, &zMapCOMPreReSampl, this->xDimSize, this->yDimSize, this->zDimSize, this->xFrom, this->xTo, this->yFrom, this->yTo, this->zFrom, this->zTo );
        
        //============================================ Re-sample map
        this->reSampleMap                             ( settings );
    
        //============================================ After re-sampling sampling rate
        xSampRate                                     = this->xDimSize / static_cast< proshade_single > ( this->xTo - this->xFrom );
        ySampRate                                     = this->yDimSize / static_cast< proshade_single > ( this->yTo - this->yFrom );
        zSampRate                                     = this->zDimSize / static_cast< proshade_single > ( this->zTo - this->zFrom );
        
        //============================================ After re-sampling first index position
        proshade_single xStartPosAfter                = static_cast< proshade_single > ( this->xFrom ) * xSampRate;
        proshade_single yStartPosAfter                = static_cast< proshade_single > ( this->yFrom ) * ySampRate;
        proshade_single zStartPosAfter                = static_cast< proshade_single > ( this->zFrom ) * zSampRate;
        
        //============================================ Translate by change in corners to make the boxes as similarly placed as possible
        proshade_single xMovHlp                       = static_cast< proshade_single > ( xStartPosAfter - xStartPosBefore );
        proshade_single yMovHlp                       = static_cast< proshade_single > ( yStartPosAfter - yStartPosBefore );
        proshade_single zMovHlp                       = static_cast< proshade_single > ( zStartPosAfter - zStartPosBefore );
        ProSHADE_internal_mapManip::moveMapByIndices  ( &xMovHlp, &yMovHlp, &zMovHlp, this->xDimSize, this->yDimSize, this->zDimSize,
                                                        &this->xFrom, &this->xTo, &this->yFrom, &this->yTo, &this->zFrom, &this->zTo,
                                                        &this->xAxisOrigin, &this->yAxisOrigin, &this->zAxisOrigin );
        
        //============================================ Find COM after map re-sampling and corner move
        proshade_double xMapCOMPostReSampl = 0.0, yMapCOMPostReSampl = 0.0, zMapCOMPostReSampl = 0.0;
        ProSHADE_internal_mapManip::findMAPCOMValues  ( this->internalMap, &xMapCOMPostReSampl, &yMapCOMPostReSampl, &zMapCOMPostReSampl, this->xDimSize, this->yDimSize, this->zDimSize, this->xFrom, this->xTo, this->yFrom, this->yTo, this->zFrom, this->zTo );
        
        //============================================ Match the COMs to get as close position of the re-sampled structure to the original as possible
        ProSHADE_internal_mapManip::moveMapByFourier  ( this->internalMap,
                                                        static_cast< proshade_single > ( xMapCOMPreReSampl - xMapCOMPostReSampl ),
                                                        static_cast< proshade_single > ( yMapCOMPreReSampl - yMapCOMPostReSampl ),
                                                        static_cast< proshade_single > ( zMapCOMPreReSampl - zMapCOMPostReSampl ),
                                                        this->xDimSize, this->yDimSize, this->zDimSize,
                                                        static_cast< proshade_signed > ( this->xDimIndices ),
                                                        static_cast< proshade_signed > ( this->yDimIndices ),
                                                        static_cast< proshade_signed > ( this->zDimIndices ) );
    }
    
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
    
    //================================================ Compute and save the COM
    this->findMapCOM                                  ( );
    this->originalMapXCom                             = this->xCom;
    this->originalMapYCom                             = this->yCom;
    this->originalMapZCom                             = this->zCom;
    
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
    this->xDimIndices                                 = static_cast< proshade_unsign > ( this->xTo );
    this->yDimIndices                                 = static_cast< proshade_unsign > ( this->yTo );
    this->zDimIndices                                 = static_cast< proshade_unsign > ( this->zTo );
    
    //================================================ Set the to indices properly
    this->xTo                                        -= 1;
    this->yTo                                        -= 1;
    this->zTo                                        -= 1;
    
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
    this->xTo                                         = this->xFrom + static_cast< proshade_signed > ( this->xDimIndices ) - 1;
    this->yTo                                         = this->yFrom + static_cast< proshade_signed > ( this->yDimIndices ) - 1;
    this->zTo                                         = this->zFrom + static_cast< proshade_signed > ( this->zDimIndices ) - 1;
    
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
    mapData.set_unit_cell                             ( static_cast< double > ( this->xDimSize ), static_cast< double > ( this->yDimSize ), static_cast< double > ( this->zDimSize ), static_cast< double > ( this->aAngle ), static_cast< double > ( this->bAngle ), static_cast< double > ( this->cAngle ) );
    mapData.set_size_without_checking                 ( static_cast< int > ( this->xDimIndices ), static_cast< int > ( this->yDimIndices ), static_cast< int > ( this->zDimIndices ) );
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
                map.grid.set_value                    ( static_cast< int > ( uIt ), static_cast< int > ( vIt ), static_cast< int > ( wIt ), static_cast<float> ( this->internalMap[arrPos] ) );
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

/*! \brief This function writes out the PDB formatted file coresponding to the structure so that its COM is at specific position.

    This function first checks if this internal structure originated from co-ordinate file (only if co-ordinates are provided can they be written out). If so,
    it will proceed to read in the original co-ordinates, rotate and translate them according to the arguments and then write the resulting co-ordinates
    into a new file.

    \param[in] fName The filename (including path) to where the output PDB file should be saved.
    \param[in] euA The Euler angle alpha by which the co-ordinates should be rotated (leave empty if no rotation is required).
    \param[in] euB The Euler angle beta by which the co-ordinates should be rotated (leave empty if no rotation is required).
    \param[in] euG The Euler angle gamma by which the co-ordinates should be rotated (leave empty if no rotation is required).
    \param[in] trsX The translation to be done along X-axis in Angstroms.
    \param[in] trsY The translation to be done along Y-axis in Angstroms.
    \param[in] trsZ The translation to be done along Z-axis in Angstroms.
    \param[in] firstModel Should only the first model, or rather all of them be used?
*/
void ProSHADE_internal_data::ProSHADE_data::writePdb ( std::string fName, proshade_double euA, proshade_double euB, proshade_double euG, proshade_double trsX, proshade_double trsY, proshade_double trsZ, bool firstModel )
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
        //============================================ Rotate the co-ordinates
        ProSHADE_internal_mapManip::rotatePDBCoordinates ( &pdbFile, euA, euB, euG, this->originalPdbRotCenX, this->originalPdbRotCenY, this->originalPdbRotCenZ, firstModel );
    }

    //================================================ Translate by required translation and the map centering (if applied)
    ProSHADE_internal_mapManip::translatePDBCoordinates ( &pdbFile, trsX, trsY, trsZ, firstModel );

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
                arrayPos                              = zIt + static_cast< proshade_signed > ( this->zDimIndices ) * ( yIt + static_cast< proshade_signed > ( this->yDimIndices ) * xIt );
                invPos                                = ( static_cast< proshade_signed > ( this->zDimIndices - 1 ) - zIt ) + static_cast< proshade_signed > ( this->zDimIndices ) * ( ( static_cast< proshade_signed > ( this->yDimIndices - 1 ) - yIt ) + static_cast< proshade_signed > ( this->yDimIndices ) * ( static_cast< proshade_signed > ( this->xDimIndices - 1 ) - xIt ) );
                
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

/*! \brief Function for normalising the map values to mean 0 and sd 1..
 
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
        ProSHADE_internal_mapManip::getNonZeroBounds  ( this->internalMap,
                                                        static_cast< proshade_signed > ( this->xDimIndices ),
                                                        static_cast< proshade_signed > ( this->yDimIndices ),
                                                        static_cast< proshade_signed > ( this->zDimIndices ),
                                                        ret );
        
        //============================================ Add the extra space
        ProSHADE_internal_mapManip::addExtraBoundSpace ( this->xDimIndices, this->yDimIndices, this->zDimIndices,
                                                         this->xDimSize, this->yDimSize, this->zDimSize, ret, settings->boundsExtraSpace );
        
        //============================================ Beautify boundaries
        ProSHADE_internal_mapManip::beautifyBoundaries ( ret, this->xDimIndices, this->yDimIndices, this->zDimIndices, settings->boundsSimilarityThreshold );
        
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
    newStr->xDimIndices                               = static_cast< proshade_unsign > ( newBounds[1] ) - static_cast< proshade_unsign > ( newBounds[0] ) + 1;
    newStr->yDimIndices                               = static_cast< proshade_unsign > ( newBounds[3] ) - static_cast< proshade_unsign > ( newBounds[2] ) + 1;
    newStr->zDimIndices                               = static_cast< proshade_unsign > ( newBounds[5] ) - static_cast< proshade_unsign > ( newBounds[4] ) + 1;
            
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
            
    newStr->xTo                                       = this->xTo - ( static_cast< proshade_signed > ( this->xDimIndices - 1 ) - newBounds[1] );
    newStr->yTo                                       = this->yTo - ( static_cast< proshade_signed > ( this->yDimIndices - 1 ) - newBounds[3] );
    newStr->zTo                                       = this->zTo - ( static_cast< proshade_signed > ( this->zDimIndices - 1 ) - newBounds[5] );
    
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
    proshade_single xMov                              = -( ( static_cast<proshade_single> ( this->xFrom ) * ( this->xDimSize / static_cast<proshade_single> ( this->xDimIndices ) - changeVals[0] ) ) -
                                                           ( static_cast<proshade_single> ( this->xFrom ) * ( this->xDimSize / static_cast<proshade_single> ( this->xDimIndices ) ) ) );
    proshade_single yMov                              = -( ( static_cast<proshade_single> ( this->yFrom ) * ( this->yDimSize / static_cast<proshade_single> ( this->yDimIndices ) - changeVals[1] ) ) -
                                                           ( static_cast<proshade_single> ( this->yFrom ) * ( this->yDimSize / static_cast<proshade_single> ( this->yDimIndices ) ) ) );
    proshade_single zMov                              = -( ( static_cast<proshade_single> ( this->zFrom ) * ( this->zDimSize / static_cast<proshade_single> ( this->zDimIndices ) - changeVals[2] ) ) -
                                                           ( static_cast<proshade_single> ( this->zFrom ) * ( this->zDimSize / static_cast<proshade_single> ( this->zDimIndices ) ) ) );

    //================================================ Move by indices (this should be sufficient)
    ProSHADE_internal_mapManip::moveMapByIndices      ( &xMov, &yMov, &zMov, this->xDimSize, this->yDimSize, this->zDimSize, &this->xFrom, &this->xTo,
                                                        &this->yFrom, &this->yTo, &this->zFrom, &this->zTo, &this->xAxisOrigin, &this->yAxisOrigin, &this->zAxisOrigin );
    
    ProSHADE_internal_mapManip::moveMapByFourier      ( this->internalMap, xMov, yMov, zMov, this->xDimSize, this->yDimSize, this->zDimSize,
                                                        static_cast< proshade_signed > ( this->xDimIndices ), static_cast< proshade_signed > ( this->yDimIndices ), static_cast< proshade_signed > ( this->zDimIndices ) );
    
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
    proshade_double xCOM                              = 0.0;
    proshade_double yCOM                              = 0.0;
    proshade_double zCOM                              = 0.0;
    
    //================================================ Find the COM location
    ProSHADE_internal_mapManip::findMAPCOMValues      ( this->internalMap, &xCOM, &yCOM, &zCOM,
                                                        this->xDimSize, this->yDimSize, this->xDimSize, this->xFrom,
                                                        this->xTo, this->yFrom, this->yTo, this->zFrom, this->zTo );
    
    //================================================ Find the sampling rates
    proshade_single xSampRate                         = static_cast< proshade_single > ( this->xDimSize ) / static_cast< proshade_single > ( this->xTo - this->xFrom );
    proshade_single ySampRate                         = static_cast< proshade_single > ( this->yDimSize ) / static_cast< proshade_single > ( this->yTo - this->yFrom );
    proshade_single zSampRate                         = static_cast< proshade_single > ( this->zDimSize ) / static_cast< proshade_single > ( this->zTo - this->zFrom );
    
    //================================================ Convert to position in indices starting from 0
    xCOM                                             /= static_cast< proshade_double > ( xSampRate );
    yCOM                                             /= static_cast< proshade_double > ( ySampRate );
    zCOM                                             /= static_cast< proshade_double > ( zSampRate );
         
    xCOM                                             -= static_cast< proshade_double > ( this->xFrom );
    yCOM                                             -= static_cast< proshade_double > ( this->yFrom );
    zCOM                                             -= static_cast< proshade_double > ( this->zFrom );
    
    //================================================ Find distance from COM to map centre in Angstroms
    proshade_double xDist                             = ( ( static_cast<proshade_double> ( this->xDimIndices ) / 2.0 ) - xCOM ) * static_cast<proshade_double> ( this->xDimSize ) / static_cast<proshade_double> ( this->xDimIndices );
    proshade_double yDist                             = ( ( static_cast<proshade_double> ( this->yDimIndices ) / 2.0 ) - yCOM ) * static_cast<proshade_double> ( this->yDimSize ) / static_cast<proshade_double> ( this->yDimIndices );
    proshade_double zDist                             = ( ( static_cast<proshade_double> ( this->zDimIndices ) / 2.0 ) - zCOM ) * static_cast<proshade_double> ( this->zDimSize ) / static_cast<proshade_double> ( this->zDimIndices );
    
    //================================================ Move the map within the box
    ProSHADE_internal_mapManip::moveMapByFourier      ( this->internalMap,
                                                        static_cast< proshade_single > ( xDist ),
                                                        static_cast< proshade_single > ( yDist ),
                                                        static_cast< proshade_single > ( zDist ),
                                                        this->xDimSize, this->yDimSize, this->zDimSize,
                                                        static_cast< proshade_signed > ( this->xDimIndices ),
                                                        static_cast< proshade_signed > ( this->yDimIndices ),
                                                        static_cast< proshade_signed > ( this->zDimIndices ) );
    
    //================================================ Note the change due to centering
    this->mapCOMProcessChangeX                       -= xDist;
    this->mapCOMProcessChangeY                       -= yDist;
    this->mapCOMProcessChangeZ                       -= zDist;
    
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
    proshade_unsign xAddIndices                       = static_cast< proshade_unsign > ( ProSHADE_internal_mapManip::myRound ( settings->addExtraSpace / this->xDimSize / static_cast<proshade_single> ( this->xDimIndices ) ) );
    proshade_unsign yAddIndices                       = static_cast< proshade_unsign > ( ProSHADE_internal_mapManip::myRound ( settings->addExtraSpace / this->yDimSize / static_cast<proshade_single> ( this->yDimIndices ) ) );
    proshade_unsign zAddIndices                       = static_cast< proshade_unsign > ( ProSHADE_internal_mapManip::myRound ( settings->addExtraSpace / this->zDimSize / static_cast<proshade_single> ( this->zDimIndices ) ) );
    
    //================================================ Update internal data variables
    this->xDimSize                                   += 2 * static_cast<proshade_single> ( xAddIndices ) * this->xDimSize / static_cast<proshade_single> ( this->xDimIndices );
    this->yDimSize                                   += 2 * static_cast<proshade_single> ( yAddIndices ) * this->yDimSize / static_cast<proshade_single> ( this->yDimIndices );
    this->zDimSize                                   += 2 * static_cast<proshade_single> ( zAddIndices ) * this->zDimSize / static_cast<proshade_single> ( this->zDimIndices );
            
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
    if ( settings->maskMap ) { this->maskMap ( settings ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Masking not requested." ); }

    //================================================ Centre map
    if ( settings->moveToCOM ) { this->centreMapOnCOM ( settings ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Map centering not requested." ); }
    
    //================================================ Add extra space
    if ( settings->addExtraSpace != 0.0f ) { this->addExtraSpace ( settings ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Extra space not requested." ); }
    
    //================================================ Remove phase, if required
    if ( !settings->usePhase ) { this->removePhaseInormation ( settings ); ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Phase information removed from the data." ); }
    else { ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Phase information retained in the data." ); }
    
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
    proshade_unsign maxDim                            = static_cast< proshade_unsign > ( std::max ( this->xDimSize, std::max ( this->yDimSize, this->zDimSize ) ) );
    proshade_unsign minDim                            = static_cast< proshade_unsign > ( std::min ( this->xDimSize, std::min ( this->yDimSize, this->zDimSize ) ) );
    proshade_unsign midDim                            = static_cast< proshade_unsign > ( 0 );
    if      ( ( this->xDimSize < static_cast< proshade_single > ( maxDim ) ) && ( this->xDimSize > static_cast< proshade_single > ( minDim ) ) ) { midDim = static_cast< proshade_unsign > ( this->xDimSize ); }
    else if ( ( this->yDimSize < static_cast< proshade_single > ( maxDim ) ) && ( this->yDimSize > static_cast< proshade_single > ( minDim ) ) ) { midDim = static_cast< proshade_unsign > ( this->yDimSize ); }
    else                                                                                                                                         { midDim = static_cast< proshade_unsign > ( this->zDimSize ); }
    
    proshade_single maxDiag                           = static_cast< proshade_single > ( std::sqrt ( std::pow ( static_cast<proshade_single> ( maxDim ), 2.0 ) +
                                                                                                     std::pow ( static_cast<proshade_single> ( midDim ), 2.0 ) ) );
    
    //================================================ Set between the points
    for ( proshade_single iter = 0.5f; ( iter * settings->maxSphereDists ) < ( maxDiag / 2.0f ); iter += 1.0f )
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


/*! \brief This function converts the internal map onto a set of concentric spheres.
 
    This function starts by determining the spherical harmonics values which were not supplied by the user, these may be
    bandwidth, taylor series cap, integration order, etc. It then proceeds to determine the optimal sphere distances, unless
    these were determined by the user.
 
    Finally, the function creates a new instance of the ProSHADE_sphere class for each of the already determined sphere
    positions. Note: The constructor of ProSHADE_sphere is where the mapping then happens.
 
    \param[in] settings A pointer to settings class containing all the information required for map manipulation.
 */
void ProSHADE_internal_data::ProSHADE_data::mapToSpheres ( ProSHADE_settings* settings )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 1, "Starting sphere mapping procedure." );
    
    //================================================ Determine spherical harmonics variables
    settings->determineAllSHValues                    ( this->xDimIndices, this->yDimIndices,
                                                        this->xDimSize,    this->yDimSize,    this->zDimSize );
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
    //================================================ Prepare FSC computation memory and variables
    fftw_complex *mapData, *origCoeffs, *fCoeffs;
    fftw_plan planForwardFourier;
    proshade_double **bindata;
    proshade_signed *binIndexing, *binCounts, noBins;
    this->prepareFSCFourierMemory                     ( mapData, origCoeffs, fCoeffs, binIndexing, &noBins, bindata, binCounts, &planForwardFourier );
    
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
        this->saveRecommendedSymmetry                 ( settings, &CSyms, &DSyms, &TSyms, &OSyms, &ISyms, axes, mapData, fCoeffs, origCoeffs, &planForwardFourier, noBins, binIndexing, bindata, binCounts );
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
        this->saveRequestedSymmetryD                  ( settings, &DSyms, axes, mapData, fCoeffs, origCoeffs, &planForwardFourier, noBins, binIndexing, bindata, binCounts );
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
    
    //================================================ Release memory after FSC computation
    delete[] mapData;
    delete[] origCoeffs;
    delete[] fCoeffs;
    fftw_destroy_plan                                 ( planForwardFourier );
    delete[] binIndexing;
    for (size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ ) { delete[] bindata[binIt]; }
    delete[] bindata;
    delete[] binCounts;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function allows using std::sort to sort vectors of ProSHADE symmetry format.
 
    \param[in] a Pointer to a ProSHADE symmetry formatted array.
    \param[in] b Pointer to a ProSHADE symmetry formatted array.
    \param[out] X  Boolean whether a is larger than b.
 */
bool sortProSHADESymmetryByFSC ( proshade_double* a, proshade_double* b)
{
    //================================================ Done
    return ( a[6] > b[6] );

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
    //================================================ Modify axis tolerance and matrix tolerance by sampling, if required by user
    if ( settings->axisErrToleranceDefault )
    {
        settings->axisErrTolerance                    = std::min ( std::max ( 0.01, ( 2.0 * M_PI ) / static_cast< proshade_double > ( this->maxShellBand ) ), 0.05 );
    }
    
    //================================================ Prepare FSC computation memory and variables
    fftw_complex *mapData, *origCoeffs, *fCoeffs;
    fftw_plan planForwardFourier;
    proshade_double **bindata;
    proshade_signed *binIndexing, *binCounts, noBins;
    this->prepareFSCFourierMemory                     ( mapData, origCoeffs, fCoeffs, binIndexing, &noBins, bindata, binCounts, &planForwardFourier );
    
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
        
        //============================================ Deal with the rotation function 0 1 0 PI issue
        bool possible010PIIssue                       = false;
        for ( size_t iter = 0; iter < CSyms.size(); iter++ ) { if ( ProSHADE_internal_maths::isAxisUnique ( &CSyms, 0.0, 1.0, 0.0, 2.0, 0.1 ) ) { possible010PIIssue = true; } }
        if ( possible010PIIssue )
        {
            proshade_double* addAxis                  = new proshade_double[7];
            addAxis[0] = 2.0; addAxis[1] = 0.0; addAxis[2] = 1.0; addAxis[3] = 0.0; addAxis[4] = M_PI; addAxis[5] = -999.9; addAxis[6] = -1.0;
            ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( &CSyms, addAxis );
            delete[] addAxis;
        }
        
        //============================================ Compute FSC for all possible axes
        for ( size_t cIt = 0; cIt < CSyms.size(); cIt++ ) { const FloatingPoint< proshade_double > lhs ( CSyms.at(cIt)[5] ), rhs ( -999.9 ); if ( ( CSyms.at(cIt)[5] > settings->peakThresholdMin ) || ( lhs.AlmostEquals ( rhs ) ) ) { this->computeFSC ( settings, &CSyms, cIt, mapData, fCoeffs, origCoeffs, &planForwardFourier, noBins, binIndexing, bindata, binCounts ); } }
        
        //============================================ Sort by FSC
        std::sort                                     ( CSyms.begin(), CSyms.end(), sortProSHADESymmetryByFSC );
        
        //============================================ Save the best axis as the recommended one
        if ( settings->detectedSymmetry.size() == 0 ) { if ( CSyms.size() > 0 ) { settings->setDetectedSymmetry ( CSyms.at(0) ); } }
        if ( CSyms.size() > 0 )
        {
            bool passedTests                          = false;
            for ( size_t cIt = 0; cIt < CSyms.size(); cIt++ )
            {
                if ( CSyms.at(0)[6] > settings->fscThreshold )
                {
                    settings->setRecommendedSymmetry  ( "C" );
                    settings->setRecommendedFold      ( settings->requestedSymmetryFold );
                    
                    ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSyms.at(0) );
                    this->saveDetectedSymmetries      ( settings, &CSyms, allCs );
                    
                    passedTests                       = true;
                    break;
                }
            }
            
            if ( !passedTests )
            {
                settings->setRecommendedSymmetry      ( "" );
                settings->setRecommendedFold          ( 0 );
            }
        }
        else
        {
            settings->setRecommendedSymmetry          ( "" );
            settings->setRecommendedFold              ( 0 );
        }
        
        //============================================ Release memory after FSC computation
        delete[] mapData;
        delete[] origCoeffs;
        delete[] fCoeffs;
        fftw_destroy_plan                             ( planForwardFourier );
        delete[] binIndexing;
        for (size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ ) { delete[] bindata[binIt]; }
        delete[] bindata;
        delete[] binCounts;
        
        //============================================ Done
        return ;
    }
    
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Starting C symmetry detection." );

    //================================================ Initialise variables
    std::vector< proshade_double* > CSyms             = getCyclicSymmetriesListFromAngleAxis ( settings );
    
    //============================================ Deal with the rotation function 0 1 0 PI issue
    bool possible010PIIssue                       = false;
    for ( size_t iter = 0; iter < CSyms.size(); iter++ ) { if ( ProSHADE_internal_maths::isAxisUnique ( &CSyms, 0.0, 1.0, 0.0, 2.0, 0.1 ) ) { possible010PIIssue = true; } }
    if ( possible010PIIssue )
    {
        proshade_double* addAxis                  = new proshade_double[7];
        addAxis[0] = 2.0; addAxis[1] = 0.0; addAxis[2] = 1.0; addAxis[3] = 0.0; addAxis[4] = M_PI; addAxis[5] = -999.9; addAxis[6] = -1.0;
        ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( &CSyms, addAxis );
        delete[] addAxis;
    }
    
    for ( size_t cIt = 0; cIt < CSyms.size(); cIt++ ) { const FloatingPoint< proshade_double > lhs ( CSyms.at(cIt)[5] ), rhs ( -999.9 ); if ( lhs.AlmostEquals ( rhs ) ) { this->computeFSC ( settings, &CSyms, cIt, mapData, fCoeffs, origCoeffs, &planForwardFourier, noBins, binIndexing, bindata, binCounts ); } }
    
    //================================================ Report progress
    std::stringstream ss;
    ss << "Detected " << CSyms.size() << " C symmetries.";
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 2, ss.str() );
    
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
        std::vector< proshade_double* > OSyms         = this->getPredictedOctahedralSymmetriesList ( settings, &CSyms );
        std::vector< proshade_double* > TSyms         = this->getPredictedTetrahedralSymmetriesList ( settings, &CSyms );
        
        //============================================ Decide on recommended symmetry
        this->saveRecommendedSymmetry                 ( settings, &CSyms, &DSyms, &TSyms, &OSyms, &ISyms, axes, mapData, fCoeffs, origCoeffs, &planForwardFourier, noBins, binIndexing, bindata, binCounts );
    }
    
    if ( settings->requestedSymmetryType == "D" )
    {
        //============================================ Run only the D symmetry detection and search for requested fold
        std::vector< proshade_double* > DSyms         = this->getDihedralSymmetriesList ( settings, &CSyms );
        this->saveRequestedSymmetryD                  ( settings, &DSyms, axes, mapData, fCoeffs, origCoeffs, &planForwardFourier, noBins, binIndexing, bindata, binCounts );
    }
    
    if ( settings->requestedSymmetryType == "T" )
    {
        //============================================ Run only the T symmetry detection
        std::vector< proshade_double* > TSyms         = this->getPredictedTetrahedralSymmetriesList ( settings, &CSyms );
        
        if ( TSyms.size() == 7 )
        {
            //======================================== Initialise decision vars
            proshade_double fscVal                    = 0.0;
            proshade_double fscValAvg                 = 0.0;
            
            //======================================== Check if axes have high enough FSC and peak height
            for ( size_t tIt = 0; tIt < 7; tIt++ ) { if ( CSyms.at(settings->allDetectedTAxes.at(tIt))[5] > settings->peakThresholdMin ) { fscVal = this->computeFSC ( settings, &CSyms, settings->allDetectedTAxes.at(tIt), mapData, fCoeffs, origCoeffs, &planForwardFourier, noBins, binIndexing, bindata, binCounts ); fscValAvg += fscVal; } }
            fscValAvg                                /= 7.0;
            
            //======================================== If C3 and C5 are found and have correct angle (must have if they are both in ISym)
            if ( fscValAvg >= ( settings->fscThreshold * 0.8 ) )
            {
                //==================================== The decision is T
                settings->setRecommendedSymmetry      ( "T" );
                settings->setRecommendedFold          ( 0 );
                for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedTAxes.size() ); it++ ) { ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSyms.at(settings->allDetectedTAxes.at(it)) ); }
                if ( settings->detectedSymmetry.size() == 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedTAxes.size() ); it++ ) { settings->setDetectedSymmetry ( CSyms.at(settings->allDetectedTAxes.at(it)) ); } }
            }
        }
    }
    
    if ( settings->requestedSymmetryType == "O" )
    {
        //============================================ Run only the O symmetry detection
        std::vector< proshade_double* > OSyms         = this->getPredictedOctahedralSymmetriesList ( settings, &CSyms );
        
        if ( OSyms.size() == 13 )
        {
            //======================================== Initialise decision vars
            proshade_double fscVal                    = 0.0;
            proshade_double fscValAvg                 = 0.0;
            
            //======================================== Check if at least one C5 and one C3 with the correct angle have high FSC and peak height
            for ( size_t oIt = 0; oIt < 13; oIt++ ) { if ( CSyms.at(settings->allDetectedOAxes.at(oIt))[5] > settings->peakThresholdMin ) { fscVal = this->computeFSC ( settings, &CSyms, settings->allDetectedOAxes.at(oIt), mapData, fCoeffs, origCoeffs, &planForwardFourier, noBins, binIndexing, bindata, binCounts ); fscValAvg += fscVal; } }
            fscValAvg                                /= 13.0;
            
            //======================================== If C3 and C5 are found and have correct angle (must have if they are both in ISym)
            if ( fscValAvg >= ( settings->fscThreshold * 0.8 ) )
            {
                //==================================== The decision is O
                settings->setRecommendedSymmetry      ( "O" );
                settings->setRecommendedFold          ( 0 );
                for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedOAxes.size() ); it++ ) { ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSyms.at(settings->allDetectedOAxes.at(it)) ); }
                if ( settings->detectedSymmetry.size() == 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedOAxes.size() ); it++ ) { settings->setDetectedSymmetry ( CSyms.at(settings->allDetectedOAxes.at(it)) ); } }
            }
        }
    }
    
    if ( settings->requestedSymmetryType == "I" )
    {
        //============================================ Run only the I symmetry detection
        std::vector< proshade_double* > ISyms         = this->getPredictedIcosahedralSymmetriesList ( settings, &CSyms );
        
        if ( ISyms.size() == 31 )
        {
            //======================================== Initialise decision vars
            proshade_double fscVal                    = 0.0;
            proshade_double fscValAvg                 = 0.0;
            
            //======================================== Check if at least one C5 and one C3 with the correct angle have high FSC and peak height
            for ( size_t iIt = 0; iIt < 31; iIt++ ) { if ( CSyms.at(settings->allDetectedIAxes.at(iIt))[5] > settings->peakThresholdMin ) { fscVal = this->computeFSC ( settings, &CSyms, settings->allDetectedIAxes.at(iIt), mapData, fCoeffs, origCoeffs, &planForwardFourier, noBins, binIndexing, bindata, binCounts ); fscValAvg += fscVal; } }
            fscValAvg                                /= 31.0;
            
            //======================================== If C3 and C5 are found and have correct angle (must have if they are both in ISym)
            if ( fscValAvg >= ( settings->fscThreshold * 0.8 ) )
            {
                //==================================== The decision is I
                settings->setRecommendedSymmetry      ( "I" );
                settings->setRecommendedFold          ( 0 );
                for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedIAxes.size() ); it++ ) { ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSyms.at(settings->allDetectedIAxes.at(it)) ); }
                if ( settings->detectedSymmetry.size() == 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedIAxes.size() ); it++ ) { settings->setDetectedSymmetry ( CSyms.at(settings->allDetectedIAxes.at(it)) ); } }
            }
        }
    }
    
    //================================================ Save C symmetries to argument and if different from settings, to the settings as well
    this->saveDetectedSymmetries                      ( settings, &CSyms, allCs );
    
    //================================================ Release memory after FSC computation
    delete[] mapData;
    delete[] origCoeffs;
    delete[] fCoeffs;
    fftw_destroy_plan                                 ( planForwardFourier );
    delete[] binIndexing;
    for (size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ ) { delete[] bindata[binIt]; }
    delete[] bindata;
    delete[] binCounts;
    
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
        ProSHADE_internal_misc::addToDoubleVector     ( &nextSym, CSyms->at(cIt)[6] );
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

/*! \brief This function finds the distinct group of axes with highest peak heights.
 
    This function starts by getting a vector of all peak heights detected in all C symmetries detected in the structure. It then proceeds to convert
    this vector into a histogram, which it then smoothens using Gaussian convolution according to the input parameters. Finally, it searches for
    peaks in the smoothened histogram function and reports the minimal value that average peak height must be in order for the axis to be
    considered part of this top group.
 
    \param[in] CSym A vector of pointers to double arrays, each array being a single Cyclic symmetry entry.
    \param[in] step The granulosity of the interval <0,1> using which the search should be done.
    \param[in] sigma The variance of the Gaussian used to smoothen the peak height histogram.
    \param[in] windowSize The width of the window over which smoothening is done.
    \param[out] threshold The minimum peak height that an axis needs to have to be considered a member of the distinct top group.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::findTopGroupSmooth ( std::vector< proshade_double* >* CSym, size_t peakPos, proshade_double step, proshade_double sigma, proshade_signed windowSize )
{
    //================================================ Initialise local variables
    proshade_double threshold                         = 0.0;
    proshade_signed totSize                           = static_cast< proshade_signed > ( ( 1.0 / step ) + 1 );
    std::vector< std::pair < proshade_double, proshade_unsign > > vals;
    std::vector< proshade_double > hist               ( static_cast< unsigned long int > ( totSize ), 0.0 );
    proshade_unsign histPos                           = 0;
    
    //================================================ Make sure window size is odd
    if ( windowSize % 2 == 0 )                        { windowSize += 1; }
    
    //================================================ Get vector of pairs of peak heights and indices in CSym array
    for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( CSym->size() ); symIt++ ) { vals.emplace_back ( std::pair < proshade_double, proshade_unsign > ( CSym->at(symIt)[peakPos], symIt ) ); }
    
    //================================================ Convert all found heights to histogram from 0.0 to 1.0 by step
    for ( proshade_double it = 0.0; it <= 1.0; it = it + step )
    {
        for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( vals.size() ); symIt++ )
        {
            //======================================== Is this height in the range?
            if ( ( vals.at(symIt).first > it ) && ( vals.at(symIt).first <= ( it + step ) ) ) { hist.at(histPos) += 1.0; }
        }
        
        //============================================ Update counter and continue
        histPos                                      += 1;
    }
    
    //================================================ Smoothen the distribution
    std::vector< proshade_double > smoothened         = ProSHADE_internal_maths::smoothen1D ( step, windowSize, sigma, hist );
    
    //================================================ Find peaks in smoothened data
    std::vector< proshade_signed > peaks              = ProSHADE_internal_peakSearch::findPeaks1D ( smoothened );
    
    //================================================ Take best peaks surroundings and produce a new set of "high" axes
    proshade_signed bestHistPos;
    if ( peaks.size() > 0 ) { bestHistPos = peaks.at(peaks.size()-1) + ( ( windowSize - 1 ) / 2 ); }
    else                    { bestHistPos = 0.0; }
    
    threshold                                         = ( static_cast< proshade_double > ( bestHistPos ) * step ) - ( windowSize * step );
    
    //================================================ Done
    return                                            ( threshold );
    
}

/*! \brief This function allocates the memory and makes all preparations required for FSC computation.
 
    \param[in] mapData The input array for Fourier transform.
    \param[in] origCoeffs The array for holding the Fourier coefficients of the original (non-rotated) density.
    \param[in] fCoeffs The array for holding the results of Fourier transform.
    \param[in] binIndexing A map that will be filled with binning indices for fast binning.
    \param[in] noBins The number of bins will be stored in this variable.
    \param[in] bindata An array to store the bin sums and other FSC computation temporary results.
    \param[in] binCounts An array that will be used to store the number of reflactions in each bin.
 */
void ProSHADE_internal_data::ProSHADE_data::prepareFSCFourierMemory ( fftw_complex*& mapData, fftw_complex*& origCoeffs, fftw_complex*& fCoeffs, proshade_signed*& binIndexing, proshade_signed* noBins, proshade_double**& bindata, proshade_signed*& binCounts, fftw_plan* planForwardFourier )
{
    //================================================ Decide number of bins and allocate which reflection belongs to which bin
    ProSHADE_internal_maths::binReciprocalSpaceReflections ( this->xDimIndices, this->yDimIndices, this->zDimIndices, noBins, binIndexing );
    
    //================================================ Allocate memory for FSC sums
    bindata                                           = new proshade_double*[*noBins];
    binCounts                                         = new proshade_signed [*noBins];
    
    //================================================ Allcate memory for bin sumation
    for ( size_t binIt = 0; binIt < static_cast< size_t > ( *noBins ); binIt++ )
    {
        bindata[binIt]                                = new proshade_double[12];
        ProSHADE_internal_misc::checkMemoryAllocation ( bindata[binIt], __FILE__, __LINE__, __func__ );
    }
    
    //================================================ Allocate memory for Fourier transform imputs and outputs
    mapData                                           = new fftw_complex [this->xDimIndices * this->yDimIndices * this->zDimIndices];
    origCoeffs                                        = new fftw_complex [this->xDimIndices * this->yDimIndices * this->zDimIndices];
    fCoeffs                                           = new fftw_complex [this->xDimIndices * this->yDimIndices * this->zDimIndices];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( mapData,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( origCoeffs,    __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( fCoeffs,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( bindata,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( binCounts,     __FILE__, __LINE__, __func__ );
    
    //================================================ Prepare memory for Fourier transform
   *planForwardFourier                                = fftw_plan_dft_3d ( static_cast< int > ( this->xDimIndices ), static_cast< int > ( this->yDimIndices ), static_cast< int > ( this->zDimIndices ), mapData, fCoeffs, FFTW_FORWARD,  FFTW_ESTIMATE );
    
    //================================================ Compute Fourier transform of the original map
    for ( size_t iter = 0; iter < static_cast< size_t > ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ ) { mapData[iter][0] = this->internalMap[iter]; mapData[iter][1] = 0.0; }
    fftw_execute                                      ( *planForwardFourier );
    for ( size_t iter = 0; iter < static_cast< size_t > ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ ) { origCoeffs[iter][0] = fCoeffs[iter][0]; origCoeffs[iter][1] = fCoeffs[iter][1]; }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function computes FSC for any given axis in the supplied CSym symmetry axes vector.
 
    This function drives the FSC computation for symmetry detection. It iterates over all rotations specified by the symmetry axis (except for
    the I rotation) and for each of these, computes the map rotation accordingly. Next, it computes the Fourier transform of such rotated map
    and procceds to bin the coefficients according to the supplied key. From these bins, the sums required for finding the FSC between the
    rotated map and the original map can be obtained. Finally, averaging the FSC'c for all rotations the final FSC for the whole symmetry axis
    is obtained and returned.
 
    \warning This function does not really work on its own, it makes plethora of assumptions, the main ones being: Fourier transform of
    original (non-rotated) map is already computed and supplied, the FFTW plan for the Fourier transform of the rotated map is prepared and
    supplied, the workspace is allocated and supplied, the binning mask and number of bins is already determined ans supplied and so on.
    It would not be hard to remove these, but then repeated calls of this function would become much slower as all of these steps would have
    to be done for each call separately.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] CSym A vector of pointers to double arrays, each array being a single Cyclic symmetry entry.
    \param[in] symIndex The index of the symmetry axis in the CSym vector for which FSC should be computed.
    \param[in] mapData FFTW complex array which will be the input for the Fourier transform done by the planForwardFourier plan.
    \param[in] fCoeffs FFTW complex array where the Fourier transform coefficients will be saved into by the planForwardFourier plan.
    \param[in] origCoeffs FFTW complex array holding already compute Fourier transform of the non-rotated map.
    \param[in] planForwardFourier A prepared FFTW3 plan for transforming the mapData onto fCoeffs.
    \param[in] noBins The number of bins as already pre-computed.
    \param[in] binIndexing A map of pre-computed bin assignments for each reflection in the format as outputted by FFTW.
    \param[in] bindata Pre-allocated array of dimensions noBins x 12 serving as workspace for the bin summation and FSC computation. This array is modified by the function in case the caller would be interested in these results.
    \param[in] binCounts Pre-allocated array of dimension noBins serving to store the bin sizes for FSC computation. This array is modified by the function in case the caller would be interested in these results.
    \param[out] fsc The FSC value found for the first (smallest) rotated map along the symmetry axis and the original map.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::computeFSC ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSym, size_t symIndex, fftw_complex* mapData, fftw_complex* fCoeffs, fftw_complex* origCoeffs, fftw_plan* planForwardFourier, proshade_signed noBins, proshade_signed *binIndexing, proshade_double**& bindata, proshade_signed*& binCounts )
{
    //================================================ Sanity check
    if ( symIndex >= CSym->size() )
    {
        std::cerr << "The supplied symmetry axes vector does not contain element number " << symIndex << ". Returning FSC 0.0." << std::endl;
        return                                        ( -2.0 );
    }
    
    //================================================ Ignore if already computed
    const FloatingPoint< proshade_double > lhs1 ( CSym->at(symIndex)[6] ), rhs1 ( -1.0 );
    if ( !lhs1.AlmostEquals ( rhs1 ) ) { return ( CSym->at(symIndex)[6] ); }
    
    //================================================ Report progress
    std::stringstream ss2;
    ss2 << "Computing FSC for symmetry C" << CSym->at(symIndex)[0] << " ( " << CSym->at(symIndex)[1] << " ; " << CSym->at(symIndex)[2] << " ; " << CSym->at(symIndex)[3] << " ) with peak height " << CSym->at(symIndex)[5];
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 4, ss2.str() );
    
    //================================================ Initialise local variables
    proshade_double *rotMap;
    
    //================================================ For each rotation along the axis
    proshade_double averageFSC                        = 0.0;
    for ( proshade_double rotIter = 1.0; rotIter < CSym->at(symIndex)[0]; rotIter += 1.0 )
    {
        //============================================ Get rotated map by the smallest fold angle along the symmetry axis
        this->rotateMapRealSpace                      ( CSym->at(symIndex)[1], CSym->at(symIndex)[2], CSym->at(symIndex)[3], ( ( 2.0 * M_PI ) / CSym->at(symIndex)[0] ) * rotIter, rotMap );
        
        //============================================ Get Fourier for the rotated map
        for ( size_t iter = 0; iter < static_cast< size_t > ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ ) { mapData[iter][0] = rotMap[iter]; mapData[iter][1] = 0.0; }
        fftw_execute                                  ( *planForwardFourier );
        
        //============================================ Clean FSC computation memory
        for ( size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ ) { for ( size_t valIt = 0; valIt < 12; valIt++ ) { bindata[binIt][valIt] = 0.0; } }
        for ( size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ ) { binCounts[binIt] = 0; }
        
        //============================================ Compute FSC
        averageFSC                                   += ProSHADE_internal_maths::computeFSC ( origCoeffs, fCoeffs, this->xDimIndices, this->yDimIndices, this->zDimIndices, noBins, binIndexing, bindata, binCounts );
        
        //============================================ Release memory
        delete[] rotMap;
    }
    
    //================================================ Convert sum to average
    averageFSC                                       /= ( CSym->at(symIndex)[0] - 1.0 );
    
    //================================================ Save result to the axis
    CSym->at(symIndex)[6]                             = averageFSC;
    
    //================================================ Report progress
    std::stringstream ss3;
    ss3 << "FSC value is " << averageFSC << " .";
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 5, ss3.str() );
    
    //================================================ Done
    return                                            ( averageFSC );
    
}

/*! \brief This function computes FSC for any given axis in the supplied CSym symmetry axes vector.
 
    This function drives the FSC computation for symmetry detection. It iterates over all rotations specified by the symmetry axis (except for
    the I rotation) and for each of these, computes the map rotation accordingly. Next, it computes the Fourier transform of such rotated map
    and procceds to bin the coefficients according to the supplied key. From these bins, the sums required for finding the FSC between the
    rotated map and the original map can be obtained. Finally, averaging the FSC'c for all rotations the final FSC for the whole symmetry axis
    is obtained and returned.
 
    \warning This function does not really work on its own, it makes plethora of assumptions, the main ones being: Fourier transform of
    original (non-rotated) map is already computed and supplied, the FFTW plan for the Fourier transform of the rotated map is prepared and
    supplied, the workspace is allocated and supplied, the binning mask and number of bins is already determined ans supplied and so on.
    It would not be hard to remove these, but then repeated calls of this function would become much slower as all of these steps would have
    to be done for each call separately.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] sym A single double array containing a single Cyclic symmetry entry in the ProSHADE format.
    \param[in] mapData FFTW complex array which will be the input for the Fourier transform done by the planForwardFourier plan.
    \param[in] fCoeffs FFTW complex array where the Fourier transform coefficients will be saved into by the planForwardFourier plan.
    \param[in] origCoeffs FFTW complex array holding already compute Fourier transform of the non-rotated map.
    \param[in] planForwardFourier A prepared FFTW3 plan for transforming the mapData onto fCoeffs.
    \param[in] noBins The number of bins as already pre-computed.
    \param[in] binIndexing A map of pre-computed bin assignments for each reflection in the format as outputted by FFTW.
    \param[in] bindata Pre-allocated array of dimensions noBins x 12 serving as workspace for the bin summation and FSC computation. This array is modified by the function in case the caller would be interested in these results.
    \param[in] binCounts Pre-allocated array of dimension noBins serving to store the bin sizes for FSC computation. This array is modified by the function in case the caller would be interested in these results.
    \param[out] fsc The FSC value found for the first (smallest) rotated map along the symmetry axis and the original map.
 */
proshade_double ProSHADE_internal_data::ProSHADE_data::computeFSC ( ProSHADE_settings* settings, proshade_double* sym, fftw_complex* mapData, fftw_complex* fCoeffs, fftw_complex* origCoeffs, fftw_plan* planForwardFourier, proshade_signed noBins, proshade_signed *binIndexing, proshade_double**& bindata, proshade_signed*& binCounts )
{
    //================================================ Ignore if already computed
    const FloatingPoint< proshade_double > lhs1 ( sym[6] ), rhs1 ( -1.0 );
    if ( !lhs1.AlmostEquals ( rhs1 ) ) { return ( sym[6] ); }
    
    //================================================ Report progress
    std::stringstream ss2;
    ss2 << "Computing FSC for symmetry C" << sym[0] << " ( " << sym[1] << " ; " << sym[2] << " ; " << sym[3] << " ) with peak height " << sym[5];
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 4, ss2.str() );
    
    //================================================ Initialise local variables
    proshade_double *rotMap;
    
    //================================================ For each rotation along the axis
    proshade_double averageFSC                        = 0.0;
    for ( proshade_double rotIter = 1.0; rotIter < sym[0]; rotIter += 1.0 )
    {
        //============================================ Get rotated map by the smallest fold angle along the symmetry axis
        this->rotateMapRealSpace                      ( sym[1], sym[2], sym[3], ( ( 2.0 * M_PI ) / sym[0] ) * rotIter, rotMap );
        
        //============================================ Get Fourier for the rotated map
        for ( size_t iter = 0; iter < static_cast< size_t > ( this->xDimIndices * this->yDimIndices * this->zDimIndices ); iter++ ) { mapData[iter][0] = rotMap[iter]; mapData[iter][1] = 0.0; }
        fftw_execute                                  ( *planForwardFourier );
        
        //============================================ Clean FSC computation memory
        for ( size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ ) { for ( size_t valIt = 0; valIt < 12; valIt++ ) { bindata[binIt][valIt] = 0.0; } }
        for ( size_t binIt = 0; binIt < static_cast< size_t > ( noBins ); binIt++ ) { binCounts[binIt] = 0; }
        
        //============================================ Compute FSC
        averageFSC                                   += ProSHADE_internal_maths::computeFSC ( origCoeffs, fCoeffs, this->xDimIndices, this->yDimIndices, this->zDimIndices, noBins, binIndexing, bindata, binCounts );
        
        //============================================ Release memory
        delete[] rotMap;
    }
    
    //================================================ Convert sum to average
    averageFSC                                       /= ( sym[0] - 1.0 );
    
    //================================================ Save result to the axis
    sym[6]                                            = averageFSC;
    
    //================================================ Report progress
    std::stringstream ss3;
    ss3 << "FSC value is " << averageFSC << " .";
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 5, ss3.str() );
    
    //================================================ Done
    return                                            ( averageFSC );
    
}

/*! \brief This function takes all the detected symmetry results and decides on which are to be recommended for this structure.
 
    This function is the brains of symmetry detection in the sense that it decides which symmetry group ProSHADE recommends as being
    detected. It starts by taking all C symmetries and building a histogram of their peak heights. From this histogram, it determines a threshold
    which contains only the most reliable axes.
 
    Next, the function tests for all axes being over this threshold for the polyhedral symmetries - I, O and T in this order. If all such
    axes (with appropriate folds) are found, their FSCs will be checked against the supplied (settings object) threshold (default: 0.80).
    If all axes pass the FSC test, then the corresponding polyhedral symmetry is determined as recommended.
 
    Should no polyhedral symmetries be found, the list of detected D symmetries will be tested next with very similar approach - both axes are
    required to pass the peak height threshold as well as the FSC threshold. Should multiple axes pairs pass, the one with the highest fold will
    be decided as the recommended one.
 
    Finally, if no dihedral symmetry is found, the C symmetries list will be searched, again with the peak height and FSC criteria. If multiple symmetry
    axes are found, the one with the highest fold will be determined as the recommended one, while if no symmetries axis passes both tests, then
    no symmetry will be returned as detected.
 
    \param[in] settings A pointer to settings class containing all the information required for map symmetry detection.
    \param[in] CSym A vector of pointers to double arrays, each array being a single Cyclic symmetry entry.
    \param[in] DSym A vector of pointers to double arrays, each array being a single Dihedral symmetry entry.
    \param[in] TSym A vector of pointers to double arrays, all of which together form the axes of tetrahedral symmetry.
    \param[in] OSym A vector of pointers to double arrays, all of which together form the axes of octahedral symmetry.
    \param[in] ISym A vector of pointers to double arrays, all of which together form the axes of icosahedral symmetry.
    \param[in] axes A vector to which all the axes of the recommended symmetry (if any) will be saved.
    \param[in] mapData FFTW complex array which will be the input for the Fourier transform done by the planForwardFourier plan.
    \param[in] fCoeffs FFTW complex array where the Fourier transform coefficients will be saved into by the planForwardFourier plan.
    \param[in] origCoeffs FFTW complex array holding already compute Fourier transform of the non-rotated map.
    \param[in] planForwardFourier A prepared FFTW3 plan for transforming the mapData onto fCoeffs.
    \param[in] noBins The number of bins as already pre-computed.
    \param[in] binIndexing A map of pre-computed bin assignments for each reflection in the format as outputted by FFTW.
    \param[in] bindata Pre-allocated array of dimensions noBins x 12 serving as workspace for the bin summation and FSC computation. This array is modified by the function in case the caller would be interested in these results.
    \param[in] binCounts Pre-allocated array of dimension noBins serving to store the bin sizes for FSC computation. This array is modified by the function in case the caller would be interested in these results.
 */
void ProSHADE_internal_data::ProSHADE_data::saveRecommendedSymmetry ( ProSHADE_settings* settings, std::vector< proshade_double* >* CSym, std::vector< proshade_double* >* DSym, std::vector< proshade_double* >* TSym, std::vector< proshade_double* >* OSym, std::vector< proshade_double* >* ISym, std::vector< proshade_double* >* axes, fftw_complex* mapData, fftw_complex* origCoeffs, fftw_complex* fCoeffs, fftw_plan* planForwardFourier, proshade_signed noBins, proshade_signed* binIndexing, proshade_double** bindata, proshade_signed* binCounts )
{
    //================================================ Report progress
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 1, "Starting recommended symmetry decision procedure." );
    
    //================================================ If no C symmetries, nothing to save...
    if ( CSym->size() == 0 )
    {
        settings->setRecommendedSymmetry              ( "" );
        settings->setRecommendedFold                  ( 0 );
        return;
    }
    
    //================================================ Find the top group minimum threshold using smoothened histogram
    proshade_double step                              = 0.01;
    proshade_double sigma                             = 0.03;
    proshade_signed windowSize                        = 9;
    proshade_double bestHistPeakStart                 = this->findTopGroupSmooth ( CSym, 5, step, sigma, windowSize );
    if ( bestHistPeakStart > settings->peakThresholdMin ) { bestHistPeakStart = settings->peakThresholdMin; }
    
    //================================================ Report progress
    proshade_unsign noPassed = 0; for ( size_t cIt = 0; cIt < CSym->size(); cIt++ ) { if ( CSym->at(cIt)[5] > bestHistPeakStart ) { noPassed += 1; } }
    std::stringstream ss;
    ss << "Smoothening has resolved in " << noPassed << " C symmetries.";
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 2, ss.str() );
    ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 2, "Starting FSC computation to confirm the C symmetries existence." );

    //================================================ Decide if I is the answer
    bool alreadyDecided                               = false;
    if ( ISym->size() == 31 )
    {
        //============================================ Initialise decision vars
        proshade_double fscVal                        = 0.0;
        proshade_double fscValAvg                     = 0.0;
        
        //============================================ Check if at least one C5 and one C3 with the correct angle have high FSC and peak height
        for ( size_t iIt = 0; iIt < 31; iIt++ ) { if ( CSym->at(settings->allDetectedIAxes.at(iIt))[5] > bestHistPeakStart ) { fscVal = this->computeFSC ( settings, CSym, settings->allDetectedIAxes.at(iIt), mapData, fCoeffs, origCoeffs, planForwardFourier, noBins, binIndexing, bindata, binCounts ); fscValAvg += fscVal; } }
        fscValAvg                                    /= 31.0;
        
        //============================================ If C3 and C5 are found and have correct angle (must have if they are both in ISym)
        if ( fscValAvg >= ( settings->fscThreshold * 0.85 ) )
        {
            //======================================== The decision is I
            settings->setRecommendedSymmetry          ( "I" );
            settings->setRecommendedFold              ( 0 );
            for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedIAxes.size() ); it++ ) { ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSym->at(settings->allDetectedIAxes.at(it)) ); }
            if ( settings->detectedSymmetry.size() == 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedIAxes.size() ); it++ ) { settings->setDetectedSymmetry ( CSym->at(settings->allDetectedIAxes.at(it)) ); } }
            
            //======================================== Done
            alreadyDecided                            = true;
        }
    }
    
    //================================================ Decide if O is the answer
    if ( ( OSym->size() == 13 ) && !alreadyDecided )
    {
        //============================================ Initialise decision vars
        proshade_double fscVal                        = 0.0;
        proshade_double fscValAvg                     = 0.0;
        
        //============================================ Check if at least one C5 and one C3 with the correct angle have high FSC and peak height
        for ( size_t oIt = 0; oIt < 13; oIt++ ) { if ( CSym->at(settings->allDetectedOAxes.at(oIt))[5] > bestHistPeakStart ) { fscVal = this->computeFSC ( settings, CSym, settings->allDetectedOAxes.at(oIt), mapData, fCoeffs, origCoeffs, planForwardFourier, noBins, binIndexing, bindata, binCounts ); fscValAvg += fscVal; } }
        fscValAvg                                    /= 13.0;
        
        //============================================ If C3 and C5 are found and have correct angle (must have if they are both in ISym)
        if ( fscValAvg >= ( settings->fscThreshold * 0.85 ) )
        {
            //======================================== The decision is O
            settings->setRecommendedSymmetry          ( "O" );
            settings->setRecommendedFold              ( 0 );
            for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedOAxes.size() ); it++ ) { ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSym->at(settings->allDetectedOAxes.at(it)) ); }
            if ( settings->detectedSymmetry.size() == 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedOAxes.size() ); it++ ) { settings->setDetectedSymmetry ( CSym->at(settings->allDetectedOAxes.at(it)) ); } }
            
            //======================================== Done
            alreadyDecided                            = true;
        }
    }
    
    //================================================ Decide if T is the answer
    if ( ( TSym->size() == 7 ) && !alreadyDecided )
    {
        //============================================ Initialise decision vars
        proshade_double fscVal                        = 0.0;
        proshade_double fscValAvg                     = 0.0;
        
        //============================================ Check if at least one C5 and one C3 with the correct angle have high FSC and peak height
        for ( size_t tIt = 0; tIt < 7; tIt++ ) { if ( CSym->at(settings->allDetectedTAxes.at(tIt))[5] > bestHistPeakStart ) { fscVal = this->computeFSC ( settings, CSym, settings->allDetectedTAxes.at(tIt), mapData, fCoeffs, origCoeffs, planForwardFourier, noBins, binIndexing, bindata, binCounts ); fscValAvg += fscVal; } }
        fscValAvg                                    /= 7.0;
        
        //============================================ If C3 and C5 are found and have correct angle (must have if they are both in ISym)
        if ( fscValAvg >= ( settings->fscThreshold * 0.85 ) )
        {
            //======================================== The decision is T
            settings->setRecommendedSymmetry          ( "T" );
            settings->setRecommendedFold              ( 0 );
            for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedTAxes.size() ); it++ ) { ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSym->at(settings->allDetectedTAxes.at(it)) ); }
            if ( settings->detectedSymmetry.size() == 0 ) { for ( proshade_unsign it = 0; it < static_cast<proshade_unsign> ( settings->allDetectedTAxes.size() ); it++ ) { settings->setDetectedSymmetry ( CSym->at(settings->allDetectedTAxes.at(it)) ); } }
            
            //======================================== Done
            alreadyDecided                            = true;
        }
    }

    //================================================ Decide if D is the answer
    if ( ( settings->allDetectedDAxes.size() > 0 ) && ( DSym->size() > 0 ) && !alreadyDecided )
    {
        //============================================ Initialise decision vars
        proshade_signed bestD                         = -1;
        proshade_unsign bestFold                      = 0;
        
        //============================================ Find FSCs
        for ( size_t dIt = 0; dIt < settings->allDetectedDAxes.size(); dIt++ )
        {
            //======================================== Do not consider more than top 20, takes time and is unlikely to produce anything...
            if ( dIt > 20 ) { continue; }
            
            //======================================== Check the peak heights
            const FloatingPoint< proshade_double > lhs999a ( CSym->at(settings->allDetectedDAxes.at(dIt).at(0))[5] ), lhs999b ( CSym->at(settings->allDetectedDAxes.at(dIt).at(1))[5] ), rhs999 ( static_cast< proshade_double > ( -999.9 ) );
            if ( ( CSym->at(settings->allDetectedDAxes.at(dIt).at(0))[5] < bestHistPeakStart ) && !( lhs999a.AlmostEquals( rhs999 ) ) ) { continue; }
            if ( ( CSym->at(settings->allDetectedDAxes.at(dIt).at(1))[5] < bestHistPeakStart ) && !( lhs999b.AlmostEquals( rhs999 ) ) ) { continue; }
            
            //======================================== Find FSCs
            this->computeFSC                          ( settings, CSym, settings->allDetectedDAxes.at(dIt).at(0), mapData, fCoeffs, origCoeffs, planForwardFourier, noBins, binIndexing, bindata, binCounts );
            this->computeFSC                          ( settings, CSym, settings->allDetectedDAxes.at(dIt).at(1), mapData, fCoeffs, origCoeffs, planForwardFourier, noBins, binIndexing, bindata, binCounts );
        }
        
        //============================================ Find FSC top group threshold
        proshade_double bestHistFSCStart              = this->findTopGroupSmooth ( CSym, 6, step, sigma, windowSize );
        
        //============================================ Check if both C symmetries are reliable
        for ( size_t dIt = 0; dIt < settings->allDetectedDAxes.size(); dIt++ )
        {
            //======================================== Check the peak heights
            const FloatingPoint< proshade_double > lhs999a2 ( CSym->at(settings->allDetectedDAxes.at(dIt).at(0))[5] ), lhs999b2 ( CSym->at(settings->allDetectedDAxes.at(dIt).at(1))[5] ), rhs999 ( static_cast< proshade_double > ( -999.9 ) );
            if ( ( CSym->at(settings->allDetectedDAxes.at(dIt).at(0))[5] < bestHistPeakStart ) && !( lhs999a2.AlmostEquals( rhs999 ) ) ) { continue; }
            if ( ( CSym->at(settings->allDetectedDAxes.at(dIt).at(1))[5] < bestHistPeakStart ) && !( lhs999b2.AlmostEquals( rhs999 ) ) ) { continue; }
            
            //======================================== Does this improve the best fold?
            if ( ( CSym->at(settings->allDetectedDAxes.at(dIt).at(0))[0] > static_cast< proshade_double > ( bestFold ) ) || ( CSym->at(settings->allDetectedDAxes.at(dIt).at(1))[0] > static_cast< proshade_double > ( bestFold ) ) )
            {
                //==================================== Check the FSC vals
                if ( CSym->at(settings->allDetectedDAxes.at(dIt).at(0))[6] < settings->fscThreshold ) { continue; }
                if ( CSym->at(settings->allDetectedDAxes.at(dIt).at(1))[6] < settings->fscThreshold ) { continue; }
                if ( std::max ( CSym->at(settings->allDetectedDAxes.at(dIt).at(0))[6], CSym->at(settings->allDetectedDAxes.at(dIt).at(1))[6] ) < bestHistFSCStart ) { continue; }
                
                //==================================== All good!
                bestFold                              = static_cast< proshade_unsign > ( std::max ( CSym->at(settings->allDetectedDAxes.at(dIt).at(0))[0], CSym->at(settings->allDetectedDAxes.at(dIt).at(1))[0] ) );
                bestD                                 = static_cast< proshade_signed > ( dIt );
            }
        }
        
        //============================================ Anything?
        if ( bestD != -1 )
        {
            //======================================== The decision is D
            settings->setRecommendedSymmetry          ( "D" );
            settings->setRecommendedFold              ( static_cast< proshade_unsign > ( std::max ( CSym->at(settings->allDetectedDAxes.at( static_cast< size_t > ( bestD ) ).at(0))[0], CSym->at(settings->allDetectedDAxes.at( static_cast< size_t > ( bestD ) ).at(1))[0] ) ) );
            ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSym->at(settings->allDetectedDAxes.at( static_cast< size_t > ( bestD ) ).at(0)) );
            ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSym->at(settings->allDetectedDAxes.at( static_cast< size_t > ( bestD ) ).at(1)) );
            if ( settings->detectedSymmetry.size() == 0 )
            {
                settings->setDetectedSymmetry         ( CSym->at(settings->allDetectedDAxes.at( static_cast< size_t > ( bestD ) ).at(0)) );
                settings->setDetectedSymmetry         ( CSym->at(settings->allDetectedDAxes.at( static_cast< size_t > ( bestD ) ).at(1)) );
            }
            
            //======================================== Done
            alreadyDecided                            = true;
        }
    }
    
    //================================================ Decide if C is the answer
    if ( ( CSym->size() > 0 ) && !alreadyDecided )
    {
        //============================================ Initialise decision vars
        proshade_signed bestC                         = -1;
        proshade_unsign bestFold                      = 0;
        
        //============================================ Find FSCs for C syms
        for ( size_t cIt = 0; cIt < CSym->size(); cIt++ )
        {
            //======================================== Do not consider more than top 20, takes time and is unlikely to produce anything...
            if ( cIt > 20 ) { continue; }
            
            //======================================== Check the peak height
            if ( CSym->at(cIt)[5]  < bestHistPeakStart ) { continue; }
            
            //======================================== Compute FSC
            this->computeFSC                          ( settings, CSym, cIt, mapData, fCoeffs, origCoeffs, planForwardFourier, noBins, binIndexing, bindata, binCounts );
        }
        
        //============================================ Find FSC top group threshold
        proshade_double bestHistFSCStart              = this->findTopGroupSmooth ( CSym, 6, step, sigma, windowSize );
        
        //============================================ Find reliable C syms
        for ( size_t cIt = 0; cIt < CSym->size(); cIt++ )
        {
            //======================================== Check if this improves the best already found fold
            if ( CSym->at(cIt)[0] > static_cast< proshade_double > ( bestFold ) )
            {
                //==================================== If FSC passes
                if ( ( CSym->at(cIt)[6] > settings->fscThreshold ) && ( CSym->at(cIt)[6] >= bestHistFSCStart ) )
                {
                    bestFold                          = static_cast< proshade_unsign > ( CSym->at(cIt)[0] );
                    bestC                             = static_cast< proshade_signed > ( cIt );
                }
            }
        }
        
        //============================================ Anything?
        if ( bestC != -1 )
        {
            //======================================== The decision is C
            settings->setRecommendedSymmetry              ( "C" );
            settings->setRecommendedFold                  ( static_cast< proshade_unsign > ( CSym->at( static_cast< size_t > ( bestC ) )[0] ) );
            ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, CSym->at( static_cast< size_t > ( bestC ) ) );
            if ( settings->detectedSymmetry.size() == 0 ) { settings->setDetectedSymmetry ( CSym->at( static_cast< size_t > ( bestC ) ) ); }
            
            //======================================== Done
            alreadyDecided                            = true;
        }
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
        const FloatingPoint< proshade_double > lhs1 ( CSym->at(iter)[0] ), rhs1 ( static_cast< proshade_double > ( settings->requestedSymmetryFold ) );
        if ( !lhs1.AlmostEquals ( rhs1 ) ) { continue; }
        
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
        settings->setRecommendedFold                  ( static_cast< proshade_unsign > ( CSym->at(bestIndex)[0] ) );
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
void ProSHADE_internal_data::ProSHADE_data::saveRequestedSymmetryD ( ProSHADE_settings* settings, std::vector< proshade_double* >* DSym, std::vector< proshade_double* >* axes, fftw_complex* mapData, fftw_complex* origCoeffs, fftw_complex* fCoeffs, fftw_plan* planForwardFourier, proshade_signed noBins, proshade_signed* binIndexing, proshade_double** bindata, proshade_signed* binCounts )
{
    //================================================ Initialise variables
    proshade_unsign bestIndex                         = 0;
    proshade_double highestSym                        = 0.0;
    
    //================================================ Search for best fold
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( DSym->size() ); iter++ )
    {
        //============================================ Check if it is tbe correct fold
        const FloatingPoint< proshade_double > lhs1 ( std::max ( DSym->at(iter)[0], DSym->at(iter)[7] ) ), rhs1 ( static_cast< proshade_double > ( settings->requestedSymmetryFold ) );
        if ( !lhs1.AlmostEquals ( rhs1 ) ) { continue; }

        //============================================ Check if peak height is decent
        const FloatingPoint< proshade_double > lhs999a ( DSym->at(iter)[5] ), lhs999b ( DSym->at(iter)[12] ), rhs999 ( static_cast< proshade_double > ( -999.9 ) );
        if ( ( DSym->at(iter)[5]  < settings->peakThresholdMin ) && !( lhs999a.AlmostEquals( rhs999 ) ) ) { continue; }
        if ( ( DSym->at(iter)[12] < settings->peakThresholdMin ) && !( lhs999b.AlmostEquals( rhs999 ) ) ) { continue; }
        
        //============================================ If correct, compute FSC
        this->computeFSC                              ( settings, &DSym->at(iter)[0], mapData, fCoeffs, origCoeffs, planForwardFourier, noBins, binIndexing, bindata, binCounts );
        this->computeFSC                              ( settings, &DSym->at(iter)[7], mapData, fCoeffs, origCoeffs, planForwardFourier, noBins, binIndexing, bindata, binCounts );
        
        //============================================ If best, store it
        if ( ( DSym->at(iter)[6] + DSym->at(iter)[13] ) > highestSym )
        {
            highestSym                                = ( DSym->at(iter)[6] + DSym->at(iter)[13] );
            bestIndex                                 = iter;
        }
    }
    
    //================================================ Found?
    if ( highestSym  > 0.0 )
    {
        settings->setRecommendedSymmetry              ( "D" );
        settings->setRecommendedFold                  ( static_cast< proshade_unsign > ( std::max ( DSym->at(bestIndex)[0], DSym->at(bestIndex)[7] ) ) );
        ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, &DSym->at(bestIndex)[0] );
        ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( axes, &DSym->at(bestIndex)[7] );
        
        if ( settings->detectedSymmetry.size() == 0 )
        {
            settings->setDetectedSymmetry             ( &DSym->at(bestIndex)[0] );
            settings->setDetectedSymmetry             ( &DSym->at(bestIndex)[7] );
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

/*! \brief This function computes the group elements as rotation matrices (except for the identity element) for any cyclic point group given as axis and fold.
 
    \param[in] xAx The x-axis element of the axis vector of the point group axis.
    \param[in] yAx The y-axis element of the axis vector of the point group axis.
    \param[in] zAx The z-axis element of the axis vector of the point group axis.
    \param[in] fold The fold of the point group.
    \param[out] val A vector containing vectors of 9 (rotation matrix) for each group element for the requested group, except for the identity element.
 */
std::vector<std::vector< proshade_double > > ProSHADE_internal_data::computeGroupElementsForGroup ( proshade_double xAx, proshade_double yAx, proshade_double zAx, proshade_signed fold )
{
    //================================================ Initialise variables
    std::vector< proshade_double > angList;
    std::vector<std::vector< proshade_double > > ret;
    
    //================================================ Allocate memory
    proshade_double* rotMat                           = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    
    
    //================================================ Normalise the axis to have magnitude of 1.0
    proshade_double normF                             = std::sqrt( std::pow ( xAx, 2.0 ) + std::pow ( yAx, 2.0 ) + std::pow ( zAx, 2.0 ) );
    xAx                                              /= normF;
    yAx                                              /= normF;
    zAx                                              /= normF;
    
    //================================================ Determine the list of angles
    if ( fold % 2 == 0 )
    {
        //============================================ If fold is even, add the negative angles
        for ( proshade_double iter = static_cast < proshade_double > ( -( ( fold / 2 ) - 1 ) ); iter <= static_cast < proshade_double > ( fold / 2 ); iter++ )
        {
            ProSHADE_internal_misc::addToDoubleVector ( &angList, ( ( 2.0 * M_PI ) / static_cast<proshade_double> ( fold ) ) * iter );
        }
    }
    else
    {
        //============================================ If fold is odd, do the same as for even, but start one index earlier
        for ( proshade_double iter = static_cast < proshade_double > ( -fold / 2 ); iter <= static_cast < proshade_double > ( fold / 2 ); iter++ )
        {
            ProSHADE_internal_misc::addToDoubleVector ( &angList, ( ( 2.0 * M_PI ) / static_cast<proshade_double> ( fold ) ) * iter );
        }
    }
    
    //================================================ For each detected angle
    for ( proshade_unsign iter = 0; iter < static_cast < proshade_unsign > ( angList.size() ); iter++ )
    {
        //============================================ Compute the rotation matrix
        ProSHADE_internal_maths::getRotationMatrixFromAngleAxis ( rotMat, xAx, yAx, zAx, angList.at(iter) );
        
        //============================================ Convert to vector of vectors of doubles and save to ret
        std::vector < proshade_double > retEl;
        for ( proshade_unsign matIt = 0; matIt < 9; matIt++ )
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

/*! \brief This function checks if the element list already contains a given matrix.
 
    \param[in] elements Vector containing all group elements.
    \param[in] elem A single element which should already be in the list.
    \param[in] matrixTolerance The maximum allowed trace error on rotation matrices to be still considered the same.
    \param[out] elementFound A boolean value stating if the element was found int the elements list or not.
 */
bool checkElementAlreadyExists ( std::vector<std::vector< proshade_double > >* elements, std::vector< proshade_double >* elem, proshade_double matrixTolerance )
{
    //================================================ Initialise variables
    bool elementFound                                 = false;
    
    //================================================ For each existing element
    for ( proshade_unsign elIt = 0; elIt < static_cast<proshade_unsign> ( elements->size() ); elIt++ )
    {
        if ( ProSHADE_internal_maths::rotationMatrixSimilarity ( &elements->at(elIt), elem, matrixTolerance ) )
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
    \param[in] matrixTolerance The maximum trace error for the matrices to be still considered the same.
    \param[out] isGroup A boolean value stating if all group element products for another group element.
 */
bool checkElementsFormGroup ( std::vector<std::vector< proshade_double > >* elements, proshade_double matrixTolerance )
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
            if ( !checkElementAlreadyExists ( elements, &product, matrixTolerance ) )
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
    \param[in] matrixTolerance The maximum trace error for rotation matrices to be still considered the same.
    \param[in] combine Should the element combinations be added as well?
    \param[out] ret A vector of group elements containing all unique elements from both input element groups.
 */
std::vector<std::vector< proshade_double > > ProSHADE_internal_data::joinElementsFromDifferentGroups ( std::vector<std::vector< proshade_double > >* first, std::vector<std::vector< proshade_double > >* second, proshade_double matrixTolerance, bool combine )
{
    //================================================ Initialise variables
    std::vector< std::vector< proshade_double > > ret;
    
    //================================================ Add the first list to ret, checking for uniqueness
    for ( proshade_unsign elIt = 0; elIt < static_cast<proshade_unsign> ( first->size() ); elIt++ )
    {
        if ( !checkElementAlreadyExists( &ret, &first->at(elIt), matrixTolerance ) )
        {
            ProSHADE_internal_misc::addToDoubleVectorVector ( &ret, first->at(elIt) );
        }
    }
    
    //================================================ Add the second list to ret, checking for uniqueness
    for ( proshade_unsign elIt = 0; elIt < static_cast<proshade_unsign> ( second->size() ); elIt++ )
    {
        if ( !checkElementAlreadyExists( &ret, &second->at(elIt), matrixTolerance ) )
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
                if ( !checkElementAlreadyExists( &ret, &product, matrixTolerance ) )
                {
                    ProSHADE_internal_misc::addToDoubleVectorVector ( &ret, product );
                }

            }
        }
    }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function returns the group elements as rotation matrices of any defined point group.
 
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
    \param[in] matrixTolerance The maximum allowed trace difference for two matrices to still be considered the same.
    \param[out] val A vector containing a vector of 9 doubles (rotation matrix) for each group element for the requested group.
 */
std::vector<std::vector< proshade_double > > ProSHADE_internal_data::ProSHADE_data::getAllGroupElements ( ProSHADE_settings* settings, std::vector< proshade_unsign > axesList, std::string groupType, proshade_double matrixTolerance )
{
    //================================================ Initialise variables
    std::vector<std::vector< proshade_double > > ret;
    
    //================================================ Select which symmetry type are we computing for
    if ( groupType == "C" )
    {
        //============================================ Sanity check
        axesToGroupTypeSanityCheck                    ( 1, static_cast< proshade_unsign > ( axesList.size() ), groupType );
        
        //============================================ Generate elements
        ret                                           = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(0)).at(1),
                                                                                       settings->allDetectedCAxes.at(axesList.at(0)).at(2),
                                                                                       settings->allDetectedCAxes.at(axesList.at(0)).at(3),
                                                                                       static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(0)).at(0) ) );

        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret, matrixTolerance ) ) { return ( ret ); }
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
        std::vector<std::vector< proshade_double > > first  = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(0)).at(1),
                                                                                             settings->allDetectedCAxes.at(axesList.at(0)).at(2),
                                                                                             settings->allDetectedCAxes.at(axesList.at(0)).at(3),
                                                                                             static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(0)).at(0) ) );
        std::vector<std::vector< proshade_double > > second = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(1)).at(1),
                                                                                             settings->allDetectedCAxes.at(axesList.at(1)).at(2),
                                                                                             settings->allDetectedCAxes.at(axesList.at(1)).at(3),
                                                                                             static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(1)).at(0) ) );
        
        //============================================ Join the element lists
        ret                                           = joinElementsFromDifferentGroups ( &first, &second, matrixTolerance, true );
        
        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret, matrixTolerance ) ) { return ( ret ); }
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
            const FloatingPoint< proshade_double > lhs1 ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ), rhs1 ( 3.0 );
            if ( lhs1.AlmostEquals ( rhs1 ) )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(1),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(2),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(3),
                                                                                                  static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ) );
                
                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, matrixTolerance, false );
            }
        }
        
        //============================================ Generate elements for all three C2 axes second
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            const FloatingPoint< proshade_double > lhs1 ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ), rhs1 ( 2.0 );
            if ( lhs1.AlmostEquals ( rhs1 ) )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(1),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(2),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(3),
                                                                                                  static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ) );
                
                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, matrixTolerance, false );
            }
        }
        
        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret, matrixTolerance ) ) { return ( ret ); }
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
            const FloatingPoint< proshade_double > lhs1 ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ), rhs1 ( 4.0 );
            if ( lhs1.AlmostEquals ( rhs1 ) )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(1),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(2),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(3),
                                                                                                  static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, matrixTolerance, false );
            }
        }

        //============================================ Generate elements for all four C3 axes first
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            const FloatingPoint< proshade_double > lhs1 ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ), rhs1 ( 3.0 );
            if ( lhs1.AlmostEquals ( rhs1 ) )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(1),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(2),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(3),
                                                                                                  static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, matrixTolerance, false );
            }
        }
        
        //============================================ Generate elements for all six C2 axes next
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            const FloatingPoint< proshade_double > lhs1 ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ), rhs1 ( 2.0 );
            if ( lhs1.AlmostEquals ( rhs1 ) )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(1),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(2),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(3),
                                                                                                  static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, matrixTolerance, false );
            }
        }

        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret, matrixTolerance ) ) { return ( ret ); }
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
            const FloatingPoint< proshade_double > lhs1 ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ), rhs1 ( 5.0 );
            if ( lhs1.AlmostEquals ( rhs1 ) )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(1),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(2),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(3),
                                                                                                  static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, matrixTolerance, false );
            }
        }
        
        //============================================ Generate elements for all ten C3 axes next
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            const FloatingPoint< proshade_double > lhs1 ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ), rhs1 ( 3.0 );
            if ( lhs1.AlmostEquals ( rhs1 ) )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(1),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(2),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(3),
                                                                                                  static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, matrixTolerance, false );
            }
        }
        
        //============================================ Generate elements for all fifteen C2 axes lastly
        for ( proshade_unsign grIt = 0; grIt < static_cast<proshade_unsign> ( axesList.size() ); grIt++ )
        {
            //======================================== If this is a C3 axis
            const FloatingPoint< proshade_double > lhs1 ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ), rhs1 ( 2.0 );
            if ( lhs1.AlmostEquals ( rhs1 ) )
            {
                //==================================== Generate the elements
                std::vector<std::vector< proshade_double > > els = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(1),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(2),
                                                                                                  settings->allDetectedCAxes.at(axesList.at(grIt)).at(3),
                                                                                                  static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ) );

                //==================================== Join the elements to any already found
                ret                                   = joinElementsFromDifferentGroups ( &els, &ret, matrixTolerance, false );
            }
        }

        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret, matrixTolerance ) ) { return ( ret ); }
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
            std::vector<std::vector< proshade_double > > els = computeGroupElementsForGroup ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(1),
                                                                                              settings->allDetectedCAxes.at(axesList.at(grIt)).at(2),
                                                                                              settings->allDetectedCAxes.at(axesList.at(grIt)).at(3),
                                                                                              static_cast< proshade_signed > ( settings->allDetectedCAxes.at(axesList.at(grIt)).at(0) ) );
            
            //======================================== Join the elements to any already found
            ret                                       = joinElementsFromDifferentGroups ( &els, &ret, matrixTolerance, true );
        }
        
        //============================================ Check the element to form a group
        if ( checkElementsFormGroup ( &ret, matrixTolerance ) ) { return ( ret ); }
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

}

/*! \brief This function copies the internal map into the supplied pointer, which it also allocates.
 
    This function is provided so that the user can provide a pointer and have it allocated and filled with the map values.
 
    \param[in] saveTo A pointer where the internal map should be deep copied into.
    \param[in] verbose How loud the run should be?
 */
void ProSHADE_internal_data::ProSHADE_data::deepCopyMap ( proshade_double*& saveTo, proshade_signed verbose )
{
    //================================================ Sanity check
    if ( saveTo != nullptr )
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
        std::vector< proshade_double > comMove        = this->getMapCOMProcessChange ( );
        ssHlp << std::endl << "Detected " << settings->recommendedSymmetryType << " symmetry with fold " << settings->recommendedSymmetryFold << " about point [" << comMove.at(0) << " , " << comMove.at(1) << " , " << comMove.at(2) << "] away from centre of mass .";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        
        if ( settings->detectedSymmetry.size() > 0 )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << "  Fold       X           Y          Z           Angle        Height         Average FSC";
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( settings->detectedSymmetry.size() ); symIt++ )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << std::showpos << std::fixed << std::setprecision(0) << "   " << settings->detectedSymmetry.at(symIt)[0] << std::setprecision(5) << "     " << settings->detectedSymmetry.at(symIt)[1] << "   " << settings->detectedSymmetry.at(symIt)[2] << "   " << settings->detectedSymmetry.at(symIt)[3] << "     " << settings->detectedSymmetry.at(symIt)[4] << "      " << settings->detectedSymmetry.at(symIt)[5] << "      " << settings->detectedSymmetry.at(symIt)[6];
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        
        std::stringstream hlpSS3;
        ssHlp.clear(); ssHlp.str ( "" );
        hlpSS3 << std::endl << "To facilitate manual checking for symmetries, the following is a list of all detected C symmetries:";
        ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, hlpSS3.str() );
        
        if ( settings->allDetectedCAxes.size() > 0 )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << "  Fold       X           Y          Z           Angle        Height         Average FSC";
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        for ( proshade_unsign symIt = 0; symIt < static_cast<proshade_unsign> ( settings->allDetectedCAxes.size() ); symIt++ )
        {
            ssHlp.clear(); ssHlp.str ( "" );
            ssHlp << std::showpos << std::fixed << std::setprecision(0) << "   " << settings->allDetectedCAxes.at(symIt)[0] << std::setprecision(5) << "     " << settings->allDetectedCAxes.at(symIt)[1] << "   " << settings->allDetectedCAxes.at(symIt)[2] << "   " << settings->allDetectedCAxes.at(symIt)[3] << "     " << settings->allDetectedCAxes.at(symIt)[4] << "      " << settings->allDetectedCAxes.at(symIt)[5] << "      " << settings->allDetectedCAxes.at(symIt)[6];
            ProSHADE_internal_messages::printProgressMessage ( settings->verbose, 0, ssHlp.str() );
        }
        
    }
    
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
                mapIt                                 = zIt  + static_cast< proshade_signed > ( this->zDimIndices ) * ( yIt  + static_cast< proshade_signed > ( this->yDimIndices ) * xIt  );
                
                //==================================== Use only positive density
                if ( this->internalMap[mapIt] <= 0.0 ) { continue; }
                
                //==================================== Compute Index COM
                this->xCom                           += this->internalMap[mapIt] * static_cast<proshade_double> ( xIt + this->xFrom );
                this->yCom                           += this->internalMap[mapIt] * static_cast<proshade_double> ( yIt + this->yFrom );
                this->zCom                           += this->internalMap[mapIt] * static_cast<proshade_double> ( zIt + this->zFrom );
                totNonZeroPoints                     += this->internalMap[mapIt];
            }
        }
    }
    
    this->xCom                                       /= totNonZeroPoints;
    this->yCom                                       /= totNonZeroPoints;
    this->zCom                                       /= totNonZeroPoints;
    
    //================================================ Convert to real world
    this->xCom         = static_cast< proshade_double > ( ( static_cast< proshade_single > ( this->xFrom ) * ( this->xDimSizeOriginal / static_cast< proshade_single > ( this->xDimIndicesOriginal ) ) ) +
                                                        ( ( static_cast< proshade_single > ( this->xCom ) - static_cast< proshade_single > ( this->xFrom ) ) *
                                                        ( static_cast< proshade_single > ( this->xDimSizeOriginal ) / static_cast< proshade_single > ( this->xDimIndicesOriginal ) ) ) );
    this->yCom         = static_cast< proshade_double > ( ( static_cast< proshade_single > ( this->yFrom ) * ( this->yDimSizeOriginal / static_cast< proshade_single > ( this->yDimIndicesOriginal ) ) ) +
                                                        ( ( static_cast< proshade_single > ( this->yCom ) - static_cast< proshade_single > ( this->yFrom ) ) *
                                                        ( static_cast< proshade_single > ( this->yDimSizeOriginal ) / static_cast< proshade_single > ( this->yDimIndicesOriginal ) ) ) );
    this->zCom         = static_cast< proshade_double > ( ( static_cast< proshade_single > ( this->zFrom ) * ( this->zDimSizeOriginal / static_cast< proshade_single > ( this->zDimIndicesOriginal ) ) ) +
                                                        ( ( static_cast< proshade_single > ( this->zCom ) - static_cast< proshade_single > ( this->zFrom ) ) *
                                                        ( static_cast< proshade_single > ( this->zDimSizeOriginal ) / static_cast< proshade_single > ( this->zDimIndicesOriginal ) ) ) );
    
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

/*! \brief This function removes phase from the map, effectively converting it to Patterson map.
 
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
    fftw_plan forward                                 = fftw_plan_dft_3d ( static_cast< int > ( this->xDimIndices ), static_cast< int > ( this->yDimIndices ), static_cast< int > ( this->zDimIndices ),
                                                                           pattersonMap, mapCoeffs, FFTW_FORWARD,  FFTW_ESTIMATE );
    fftw_plan inverse                                 = fftw_plan_dft_3d ( static_cast< int > ( this->xDimIndices ), static_cast< int > ( this->yDimIndices ), static_cast< int > ( this->zDimIndices ),
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
                mapIt                                 = zIt  + static_cast< proshade_signed > ( this->zDimIndices ) * ( yIt  + static_cast< proshade_signed > ( this->yDimIndices ) * xIt  );
                patIt                                 = patZ + static_cast< proshade_signed > ( this->zDimIndices ) * ( patY + static_cast< proshade_signed > ( this->yDimIndices ) * patX );
                
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
    return                                            ( &this->sphericalHarmonics[shell][seanindex ( static_cast< int > ( order ) - static_cast< int > ( band ),
                                                                                                     static_cast< int > ( band ),
                                                                                                     static_cast< int > ( this->spheres[shell]->getLocalBandwidth() ) )][0] );
    
}

/*! \brief This function allows access to the private internal imaginary spherical harmonics values.
 
    \param[out] X Pointer to the value of the internal private spherical harmonics imaginary value of the given index.
 */
proshade_double* ProSHADE_internal_data::ProSHADE_data::getImagSphHarmValue ( proshade_unsign band, proshade_unsign order, proshade_unsign shell )
{
    //================================================ Done
    return                                            ( &this->sphericalHarmonics[shell][seanindex ( static_cast< int > ( order ) - static_cast< int > ( band ),
                                                                                                     static_cast< int > ( band ),
                                                                                                     static_cast< int > ( this->spheres[shell]->getLocalBandwidth() ) )][1] );
    
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
proshade_single ProSHADE_internal_data::ProSHADE_data::getSpherePosValue ( proshade_unsign shell )
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

/*! \brief This function allows access to the translation caused by structure processing.

    \param[out] mapCOMProcessChange Vector of the distances in Angstroms that the structure has been moved internally.
*/
std::vector< proshade_double > ProSHADE_internal_data::ProSHADE_data::getMapCOMProcessChange ( void )
{
    //================================================ Initialise local variables
    std::vector< proshade_double > ret;
    
    //================================================ Save the values
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, this->mapCOMProcessChangeX );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, this->mapCOMProcessChangeY );
    ProSHADE_internal_misc::addToDoubleVector         ( &ret, this->mapCOMProcessChangeZ );
    
    //================================================ Return the requested value
    return                                            ( ret );
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
    return                                            ( static_cast<int> ( so3CoefLoc ( static_cast< int > ( order1 ), static_cast< int > ( order2 ), static_cast< int > ( band ), static_cast< int > ( this->getMaxBand() ) ) ) );
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
    ProSHADE_internal_maths::getEulerZXZFromSOFTPosition ( static_cast< int > ( this->getMaxBand() ), aI, bI, gI, &eA, &eB, &eG );
    
    //================================================ Prepare internal rotation matrix memory
    proshade_double* rMat                             = nullptr;
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
    \param[in] eulA The Euler alpha angle value, by which the moving structure is to be rotated by.
    \param[in] eulB The Euler beta angle value, by which the moving structure is to be rotated by.
    \param[in] eulG The Euler gamma angle value, by which the moving structure is to be rotated by.
    \param[in] rotCentre The rotation centre position.
    \param[in] ultimateTranslation The final translation as determined by the translation function.
*/
void ProSHADE_internal_data::ProSHADE_data::writeOutOverlayFiles ( ProSHADE_settings* settings, proshade_double eulA, proshade_double eulB, proshade_double eulG, std::vector< proshade_double >* rotCentre, std::vector< proshade_double >* ultimateTranslation )
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
        this->writePdb                                ( fNameHlp.str(), eulA, eulB, eulG, ultimateTranslation->at(0), ultimateTranslation->at(1), ultimateTranslation->at(2), settings->firstModelOnly );
    }
    
    //================================================ Write out the json file with the results
    ProSHADE_internal_io::writeRotationTranslationJSON ( rotCentre->at(0), rotCentre->at(1), rotCentre->at(2),
                                                         eulA, eulB, eulG,
                                                         ultimateTranslation->at(0), ultimateTranslation->at(1), ultimateTranslation->at(2),
                                                         settings->rotTrsJSONFile );
    
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
void ProSHADE_internal_data::ProSHADE_data::reportOverlayResults ( ProSHADE_settings* settings, std::vector < proshade_double >* rotationCentre, std::vector < proshade_double >* eulerAngles, std::vector < proshade_double >* finalTranslation )
{
    //================================================ Empty line
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, "" );
    
    //================================================ Write out rotation centre translation results
    std::stringstream rotCen; rotCen << std::setprecision (3) << std::showpos << "The rotation centre to origin translation vector is:  " << -rotationCentre->at(0) << "     " << -rotationCentre->at(1) << "     " << -rotationCentre->at(2);
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, rotCen.str() );
    
    //================================================ Write out rotation matrix about origin
    proshade_double* rotMat                           = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( eulerAngles->at(0), eulerAngles->at(1), eulerAngles->at(2), rotMat );
    
    std::stringstream rotMatSS;  
    rotMatSS << std::setprecision (3) << std::showpos << "The rotation matrix about origin is                 : " << rotMat[0] << "     " << rotMat[1] << "     " << rotMat[2] << std::endl;
    rotMatSS << std::setprecision (3) << std::showpos << "                                                    : " << rotMat[3] << "     " << rotMat[4] << "     " << rotMat[5] << std::endl;
    rotMatSS << std::setprecision (3) << std::showpos << "                                                    : " << rotMat[6] << "     " << rotMat[7] << "     " << rotMat[8];
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, rotMatSS.str() );
    
    delete[] rotMat;
    
    //================================================ Write out origin to overlay translation results
    std::stringstream finTrs; finTrs << std::setprecision (3) << std::showpos << "The rotation centre to overlay translation vector is: " << finalTranslation->at(0) << "     " << finalTranslation->at(1) << "     " << finalTranslation->at(2);
    ProSHADE_internal_messages::printProgressMessage  ( settings->verbose, 0, finTrs.str() );
    
    //================================================ Done
    return ;
    
}
