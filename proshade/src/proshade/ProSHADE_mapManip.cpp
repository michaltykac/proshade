/*! \file ProSHADE_mapManip.cpp
    \brief This source file contains the functions required for internal map manipulation for various purposes
 
    The functions in this source file are grouped here as they all change either the internal map or the gemmi library objects from which map will be computed. These
    functions include, for example, phase-removal, trilinear interpolation and changing map boundaries functionalities.
 
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
#include "ProSHADE_mapManip.hpp"

//==================================================== Define round for C++98
/*! \brief Calls the appropriate version of round function depending on compiler version.
 
    \param[in] x A decimal point number to be rounded.
    \param[out] X The rounded number.
 */
proshade_signed ProSHADE_internal_mapManip::myRound ( proshade_double x )
{
#if __cplusplus >= 201103L
    return                                            ( static_cast< proshade_signed > ( std::round ( x ) ) );
#else
    return                                            ( static_cast< proshade_signed > ( round ( x ) ) );
#endif
}

//==================================================== Define round for C++98
/*! \brief Calls the appropriate version of round function depending on compiler version.
 
    \param[in] x A decimal point number to be rounded.
    \param[out] X The rounded number.
 */
proshade_signed ProSHADE_internal_mapManip::myRound ( proshade_single x )
{
#if __cplusplus >= 201103L
    return                                            ( static_cast< proshade_signed > ( std::round ( x ) ) );
#else
    return                                            ( static_cast< proshade_signed > ( round ( x ) ) );
#endif
}

/*! \brief Function for finding the PDB file ranges.
 
    This function does a quick read-through the PDB file and reports the x, y and z to and from values. This is used to determine if these
    need to be changed for proper theoretical map computation operation.
 
    \param[in] pdbFile A gemmi::Structure object read in from the input file.
    \param[in] xFrom Address to a variable to save the x axis minimal atom position.
    \param[in] xTo Address to a variable to save the x axis maximum atom position.
    \param[in] yFrom Address to a variable to save the y axis minimal atom position.
    \param[in] yTo Address to a variable to save the y axis maximum atom position.
    \param[in] zFrom Address to a variable to save the z axis minimal atom position.
    \param[in] zTo Address to a variable to save the z axis maximum atom position.
    \param[in] firstModel Should only the first, or all models be used?
 
    \warning This function ignores hydrogen atoms, as they will not be used in map computation by Gemmi and their inclusion causes
    mismatches between the map and the co-ordinates as they now contain different contents.
 */
void ProSHADE_internal_mapManip::determinePDBRanges ( gemmi::Structure pdbFile, proshade_single* xFrom, proshade_single* xTo, proshade_single* yFrom, proshade_single* yTo, proshade_single* zFrom, proshade_single* zTo, bool firstModel )
{
    //================================================ Initialise structure crawl
    bool firstAtom                                    = true;
  
    //================================================ Use the first model, if it exists
    if ( pdbFile.models.size() > 0 )
    {
        //============================================ For each model
        for ( proshade_unsign sIt = 0; sIt < static_cast<proshade_unsign> ( pdbFile.models.size() ); sIt++ )
        {
            //======================================== Check if multiple models are allowed
            if ( firstModel && ( sIt != 0 ) ) { break; }
            
            //======================================== Get model
            gemmi::Model model                        = pdbFile.models.at(sIt);
        
            //======================================== For each chain
            for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( model.chains.size() ); mIt++ )
            {
                //==================================== Get chain
                gemmi::Chain chain                    = model.chains.at(mIt);
                
                //==================================== For each residue
                for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( chain.residues.size() ); rIt++ )
                {
                    //================================ Get residue
                    gemmi::Residue residue            = chain.residues.at(rIt);
                    
                    //================================ For each atom
                    for ( proshade_unsign aIt = 0; aIt < static_cast<proshade_unsign> ( residue.atoms.size() ); aIt++ )
                    {
                        //============================ Get atom
                        gemmi::Atom atom              = residue.atoms.at(aIt);
                        
                        //============================ Ignore hydrogens, map computations ignore them anyway and inclusion here causes map - co-ordinate mismatches.
                        if ( atom.is_hydrogen() ) { continue; }
                        
                        //============================ Find the coordinate ranges
                        if ( firstAtom )
                        {
                            *xTo                      = static_cast<proshade_single> ( atom.pos.x );
                            *xFrom                    = static_cast<proshade_single> ( atom.pos.x );
                            *yTo                      = static_cast<proshade_single> ( atom.pos.y );
                            *yFrom                    = static_cast<proshade_single> ( atom.pos.y );
                            *zTo                      = static_cast<proshade_single> ( atom.pos.z );
                            *zFrom                    = static_cast<proshade_single> ( atom.pos.z );
                            firstAtom                 = false;
                        }
                        else
                        {
                            if ( static_cast<proshade_single> ( atom.pos.x ) > *xTo   ) { *xTo   = static_cast<proshade_single> ( atom.pos.x ); }
                            if ( static_cast<proshade_single> ( atom.pos.x ) < *xFrom ) { *xFrom = static_cast<proshade_single> ( atom.pos.x ); }
                            if ( static_cast<proshade_single> ( atom.pos.y ) > *yTo   ) { *yTo   = static_cast<proshade_single> ( atom.pos.y ); }
                            if ( static_cast<proshade_single> ( atom.pos.y ) < *yFrom ) { *yFrom = static_cast<proshade_single> ( atom.pos.y ); }
                            if ( static_cast<proshade_single> ( atom.pos.z ) > *zTo   ) { *zTo   = static_cast<proshade_single> ( atom.pos.z ); }
                            if ( static_cast<proshade_single> ( atom.pos.z ) < *zFrom ) { *zFrom = static_cast<proshade_single> ( atom.pos.z ); }
                        }
                    }
                }
            }
        }
    }
    else
    {
        std::stringstream hlpSS;
        hlpSS << "Found 0 models in input file " << pdbFile.name << ".\n                    : This suggests that the input co-ordinate file is\n                    : corrupted or mis-formatted.";
        throw ProSHADE_exception ( "Found no model in co-ordinate file.", "EP00050", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the Centre of Mass for the co-ordinate file.

    This function takes the gemmi::Structure object read from a co-odinate file and procceds to compute the Centre of Mass of this object.

    \param[in] pdbFile A gemmi::Structure object read in from input file.
    \param[in] xCom A pointer to proshade_double variable where the PDB file COM along the X-axis will be saved.
    \param[in] yCom A pointer to proshade_double variable where the PDB file COM along the Y-axis will be saved.
    \param[in] zCom A pointer to proshade_double variable where the PDB file COM along the Z-axis will be saved.
    \param[in] firstModel Should only the first, or all models be used?
*/
void ProSHADE_internal_mapManip::findPDBCOMValues ( gemmi::Structure pdbFile, proshade_double *xCom, proshade_double *yCom, proshade_double *zCom, bool firstModel )
{
    //================================================ Initialise structure crawl
    proshade_double totAtoms                          = 0.0;
   *xCom                                              = 0.0;
   *yCom                                              = 0.0;
   *zCom                                              = 0.0;
    
    //================================================ Use the first model, if it exists
    if ( pdbFile.models.size() > 0 )
    {
        //============================================ For each model
        for ( proshade_unsign sIt = 0; sIt < static_cast<proshade_unsign> ( pdbFile.models.size() ); sIt++ )
        {
            //======================================== Get model
            gemmi::Model model                        = pdbFile.models.at(sIt);
            
            //======================================== Check if multiple models are allowed
            if ( firstModel && ( sIt != 0 ) ) { break; }
            
            //======================================== For each chain
            for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( model.chains.size() ); mIt++ )
            {
                //==================================== Get chain
                gemmi::Chain chain                    = model.chains.at(mIt);
                
                //==================================== For each residue
                for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( chain.residues.size() ); rIt++ )
                {
                    //================================ Get residue
                    gemmi::Residue residue            = chain.residues.at(rIt);
                    
                    //================================ For each atom
                    for ( proshade_unsign aIt = 0; aIt < static_cast<proshade_unsign> ( residue.atoms.size() ); aIt++ )
                    {
                        //============================ Get atom
                        gemmi::Atom atom              = residue.atoms.at(aIt);
                        
                        //============================ Save the COM sums
                       *xCom                         += atom.pos.x * atom.element.weight();
                       *yCom                         += atom.pos.y * atom.element.weight();
                       *zCom                         += atom.pos.z * atom.element.weight();
                        totAtoms                     += atom.element.weight();
                    }
                }
            }
        }
    }
    else
    {
        std::stringstream hlpSS;
        hlpSS << "Found 0 models in input file " << pdbFile.name << ".\n                    : This suggests that the input co-ordinate file is\n                    : corrupted or mis-formatted.";
        throw ProSHADE_exception ( "Found no model in co-ordinate file.", "EP00050", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Normalise sums to COM
   *xCom                                             /= totAtoms;
   *yCom                                             /= totAtoms;
   *zCom                                             /= totAtoms;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function finds the Centre of Mass for a map.

    This function takes the proshade map object and procceds to compute the Centre of Mass of this object in Angstroms in real space.

    \param[in] map A pointer to prshade density map.
    \param[in] xCom A pointer to proshade_double variable where the MAP file COM along the X-axis will be saved.
    \param[in] yCom A pointer to proshade_double variable where the MAP file COM along the Y-axis will be saved.
    \param[in] zCom A pointer to proshade_double variable where the MAP file COM along the Z-axis will be saved.
    \param[in] xAngs The map size in Angstroms along the X axis.
    \param[in] yAngs The map size in Angstroms along the Y axis.
    \param[in] zAngs The map size in Angstroms along the Z axis.
    \param[in] xFrom The initial index of the x dimension of the map.
    \param[in] xTo The terminal index of the x dimension of the map.
    \param[in] yFrom The initial index of the y dimension of the map.
    \param[in] yTo The terminal index of the y dimension of the map.
    \param[in] zFrom The initial index of the z dimension of the map.
    \param[in] zTo The terminal index of the z dimension of the map.
*/
void ProSHADE_internal_mapManip::findMAPCOMValues ( proshade_double* map, proshade_double *xCom, proshade_double *yCom, proshade_double *zCom, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_signed xFrom, proshade_signed xTo, proshade_signed yFrom, proshade_signed yTo, proshade_signed zFrom, proshade_signed zTo )
{
    //================================================ Initialise computation
    proshade_double totDensity                        = 0.0;
   *xCom                                              = 0.0;
   *yCom                                              = 0.0;
   *zCom                                              = 0.0;
    proshade_signed arrPos                            = 0;
    proshade_single xSampRate                         = xAngs / static_cast< proshade_single > ( xTo - xFrom );
    proshade_single ySampRate                         = yAngs / static_cast< proshade_single > ( yTo - yFrom );
    proshade_single zSampRate                         = zAngs / static_cast< proshade_single > ( zTo - zFrom );
    
    //================================================ For each map point
    for ( proshade_signed xIt = xFrom; xIt <= xTo; xIt++ )
    {
        for ( proshade_signed yIt = yFrom; yIt <= yTo; yIt++ )
        {
            for ( proshade_signed zIt = zFrom; zIt <= zTo; zIt++ )
            {
                arrPos                                = (zIt-zFrom) + ( zTo - zFrom + 1 ) * ( ( yIt - yFrom ) + ( yTo - yFrom + 1 ) * ( xIt - xFrom ) );
                const FloatingPoint< proshade_double > lhs ( map[arrPos] );
                if ( !lhs.AlmostEquals ( lhs ) )       { map[arrPos] = 0.0; continue; }
                
                if ( map[arrPos] > 0.0 )
                {
                    totDensity                       += map[arrPos];
                    *xCom                            += static_cast<proshade_double> ( static_cast< proshade_single > ( xIt ) * xSampRate ) * map[arrPos];
                    *yCom                            += static_cast<proshade_double> ( static_cast< proshade_single > ( yIt ) * ySampRate ) * map[arrPos];
                    *zCom                            += static_cast<proshade_double> ( static_cast< proshade_single > ( zIt ) * zSampRate ) * map[arrPos];
                }
            }
        }
    }
    
    //================================================ Normalise sums to COM
   *xCom                                             /= totDensity;
   *yCom                                             /= totDensity;
   *zCom                                             /= totDensity;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for rotating the PDB file co-ordinates by Euler angles.
 
    This function takes the three Euler angles and a pointer to a gemmi::Structure and it then proceeds to compute the rotation matrix from the Euler angles. This matrix is then applied to
    the co-ordinates in the gemmi::Structure in a way so that the rotation is done over the centre of the co-ordinates, but the co-ordinate positions stay unchanged.
 
    \param[in] pdbFile Pointer to a gemmi::Structure object which will have its co-ordinates rotated.
    \param[in] euA The Euler angle alpha by which the co-ordinates should be rotated.
    \param[in] euB The Euler angle beta by which the co-ordinates should be rotated.
    \param[in] euG The Euler angle gamma by which the co-ordinates should be rotated.
    \param[in] xCom The x-axis position around which the rotation should be applied.
    \param[in] yCom The y-axis position around which the rotation should be applied.
    \param[in] zCom The z-axis position around which the rotation should be applied.
    \param[in] firstModel Should only the first, or all models be used?
 */
void ProSHADE_internal_mapManip::rotatePDBCoordinates ( gemmi::Structure *pdbFile, proshade_double euA, proshade_double euB, proshade_double euG, proshade_double xCom,
proshade_double yCom, proshade_double zCom, bool firstModel )
{
    //================================================ Convert Euler angles to rotation matrix
    proshade_double *rotMat                           = new proshade_double[9];
    ProSHADE_internal_misc::checkMemoryAllocation     ( rotMat, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_maths::getRotationMatrixFromEulerZXZAngles ( euG, euB, euA, rotMat );
    
    //================================================ Initialise internal variables
    proshade_double xTmp, yTmp, zTmp;
    
    //================================================ Use the first model, if it exists
    if ( pdbFile->models.size() > 0 )
    {
        //============================================ For each model
        for ( proshade_unsign sIt = 0; sIt < static_cast<proshade_unsign> ( pdbFile->models.size() ); sIt++ )
        {
            //======================================== Get model
            gemmi::Model *model                       = &pdbFile->models.at(sIt);
            
            //======================================== Check if multiple models are allowed
            if ( firstModel && ( sIt != 0 ) ) { break; }
            
            //======================================== For each chain
            for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( model->chains.size() ); mIt++ )
            {
                //==================================== Get chain
                gemmi::Chain *chain                   = &model->chains.at(mIt);
                
                //==================================== For each residue
                for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( chain->residues.size() ); rIt++ )
                {
                    //================================ Get residue
                    gemmi::Residue *residue           = &chain->residues.at(rIt);
                    
                    //================================ For each atom
                    for ( proshade_unsign aIt = 0; aIt < static_cast<proshade_unsign> ( residue->atoms.size() ); aIt++ )
                    {
                        //============================ Get atom
                        gemmi::Atom *atom             = &residue->atoms.at(aIt);
                        
                        //============================ Move to mid-point
                        xTmp                          = static_cast< proshade_double > ( atom->pos.x - xCom );
                        yTmp                          = static_cast< proshade_double > ( atom->pos.y - yCom );
                        zTmp                          = static_cast< proshade_double > ( atom->pos.z - zCom );
                        
                        //============================ Rotate the atom position
                        atom->pos.x                   = ( xTmp * rotMat[0] ) + ( yTmp * rotMat[1] ) + ( zTmp * rotMat[2] );
                        atom->pos.y                   = ( xTmp * rotMat[3] ) + ( yTmp * rotMat[4] ) + ( zTmp * rotMat[5] );
                        atom->pos.z                   = ( xTmp * rotMat[6] ) + ( yTmp * rotMat[7] ) + ( zTmp * rotMat[8] );
                        
                        //============================ Move back
                        atom->pos.x                   = atom->pos.x + xCom;
                        atom->pos.y                   = atom->pos.y + yCom;
                        atom->pos.z                   = atom->pos.z + zCom;
                    }
                }
            }
        }
    }
    else
    {
        std::stringstream hlpSS;
        hlpSS << "Found 0 models in input file " << pdbFile->name << ".\n                    : This suggests that the input co-ordinate file is\n                    : corrupted or mis-formatted.";
        throw ProSHADE_exception ( "Found no model in co-ordinate file.", "EP00050", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Release memory
    delete[] rotMat;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for translating the PDB file co-ordinates by given distances in Angstroms.
 
    This function simply iterates through the given structure object and adds the required shift to all atoms in the first model along the three axes.
 
    \param[in] pdbFile Pointer to a gemmi::Structure object whose co-ordinates are to be translated.
    \param[in] transX The
    \param[in] transY The
    \param[in] transZ The
    \param[in] firstModel Should only the first, or all models be used?
 */
void ProSHADE_internal_mapManip::translatePDBCoordinates ( gemmi::Structure *pdbFile, proshade_double transX, proshade_double transY, proshade_double transZ, bool firstModel )
{
    //================================================ Use the first model, if it exists
    if ( pdbFile->models.size() > 0 )
    {
        //============================================ For each model
        for ( proshade_unsign sIt = 0; sIt < static_cast<proshade_unsign> ( pdbFile->models.size() ); sIt++ )
        {
            //======================================== Check if multiple models are allowed
            if ( firstModel && ( sIt != 0 ) ) { break; }
            
            //======================================== Get model
            gemmi::Model *model                       = &pdbFile->models.at(sIt);
            
            //======================================== For each chain
            for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( model->chains.size() ); mIt++ )
            {
                //==================================== Get chain
                gemmi::Chain *chain                   = &model->chains.at(mIt);
                
                //==================================== For each residue
                for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( chain->residues.size() ); rIt++ )
                {
                    //================================ Get residue
                    gemmi::Residue *residue           = &chain->residues.at(rIt);
                    
                    //================================ For each atom
                    for ( proshade_unsign aIt = 0; aIt < static_cast<proshade_unsign> ( residue->atoms.size() ); aIt++ )
                    {
                        //============================ Get atom
                        gemmi::Atom *atom             = &residue->atoms.at(aIt);
                        
                        //============================ Translate
                        atom->pos.x                  += transX;
                        atom->pos.y                  += transY;
                        atom->pos.z                  += transZ;
                    }
                }
            }
        }
    }
    else
    {
        std::stringstream hlpSS;
        hlpSS << "Found 0 models in input file " << pdbFile->name << ".\n                    : This suggests that the input co-ordinate file is\n                    : corrupted or mis-formatted.";
        throw ProSHADE_exception ( "Found no model in co-ordinate file.", "EP00050", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for changing the PDB B-factor values to a specific single value.
 
    This function does a quick read-through the PDB file and changes all the atom B-factor values to a single, pre-set value. This is
    important for PDB files which have all atoms with B-factor value 0, as these then produce bad maps, which fail further processing
    by ProSHADE.
 
    \param[in] pdbFile A pointer to gemmi::Structure object read in from the input file.
    \param[in] newBFactorValue The value with which all atom B-factor values should be replaced.
    \param[in] firstModel Should only the first, or all models be used?
 */
void ProSHADE_internal_mapManip::changePDBBFactors ( gemmi::Structure *pdbFile, proshade_double newBFactorValue, bool firstModel )
{
    //================================================ Use the first model, if it exists
    if ( pdbFile->models.size() > 0 )
    {
        //============================================ For each model
        for ( proshade_unsign sIt = 0; sIt < static_cast<proshade_unsign> ( pdbFile->models.size() ); sIt++ )
        {
            //======================================== Check if multiple models are allowed
            if ( firstModel && ( sIt != 0 ) ) { break; }
            
            //======================================== Get model
            gemmi::Model *model                       = &pdbFile->models.at(sIt);
            
            //======================================== For each chain
            for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( model->chains.size() ); mIt++ )
            {
                //==================================== Get chain
                gemmi::Chain *chain                   = &model->chains.at(mIt);
                
                //==================================== For each residue
                for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( chain->residues.size() ); rIt++ )
                {
                    //================================ Get residue
                    gemmi::Residue *residue           = &chain->residues.at(rIt);
                    
                    //================================ For each atom
                    for ( proshade_unsign aIt = 0; aIt < static_cast<proshade_unsign> ( residue->atoms.size() ); aIt++ )
                    {
                        //============================ Get atom
                        gemmi::Atom *atom             = &residue->atoms.at(aIt);
                        
                        //============================ Change the B-factors
                        atom->b_iso                   = static_cast< float > ( newBFactorValue );
                    }
                }
            }
        }
    }
    else
    {
        std::stringstream hlpSS;
        hlpSS << "Found 0 models in input file " << pdbFile->name << ".\n                    : This suggests that the input co-ordinate file is\n                    : corrupted or mis-formatted.";
        throw ProSHADE_exception ( "Found no model in co-ordinate file.", "EP00050", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function removed all waters from PDB input file.
 
    This function does a quick read-through the PDB file and deletes all water molecules from the PDB file. This is important as
    comparing two similar structures, one with hundreds of waters (with high occupancy in some cases) and one without will lead
    to the structures being different in ProSHADE's eyes.
 
    \param[in] pdbFile A pointer to gemmi::Structure object read in from the input file.
    \param[in] firstModel Should only the first, or all models be used?
 */
void ProSHADE_internal_mapManip::removeWaters ( gemmi::Structure *pdbFile, bool firstModel )
{
    //================================================ Use the first model, if it exists
    if ( pdbFile->models.size() > 0 )
    {
        //============================================ For each model
        for ( proshade_unsign sIt = 0; sIt < static_cast<proshade_unsign> ( pdbFile->models.size() ); sIt++ )
        {
            //======================================== Check if multiple models are allowed
            if ( firstModel && ( sIt != 0 ) ) { break; }
            
            //======================================== Get model
            gemmi::Model *model                       = &pdbFile->models.at(sIt);
            
            //======================================== For each chain
            for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( model->chains.size() ); mIt++ )
            {
                //==================================== Get chain
                gemmi::Chain *chain                   = &model->chains.at(mIt);
                
                //==================================== Initialise del vector
                std::vector< proshade_unsign > delVec;
                
                //==================================== For each residue
                for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( chain->residues.size() ); rIt++ )
                {
                    //================================ Get residue
                    gemmi::Residue *residue           = &chain->residues.at(rIt);
                    
                    //================================ If residue is water
                    if ( residue->is_water() )
                    {
                        ProSHADE_internal_misc::addToUnsignVector ( &delVec, rIt );
                    }
                }
                
                //==================================== Delete from end to avoid indexing issues
                std::sort                             ( delVec.begin(), delVec.end(), std::greater<int>() );
                for ( proshade_unsign vecIt = 0; vecIt < static_cast<proshade_unsign> ( delVec.size() ); vecIt++ )
                {
                    chain->residues.erase             ( chain->residues.begin() + static_cast< long int > ( delVec.at(vecIt) ) );
                }
            }
        }
    }
    else
    {
        std::stringstream hlpSS;
        hlpSS << "Found 0 models in input file " << pdbFile->name << ".\n                    : This suggests that the input co-ordinate file is\n                    : corrupted or mis-formatted.";
        throw ProSHADE_exception ( "Found no model in co-ordinate file.", "EP00050", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for moving co-ordinate atoms to better suit theoretical map computation.
 
    This function translates all atoms by a given x, y and z distances. This is required as theoretical map computation can only output map cells starting from
    0, 0, 0 and therefore to avoid density being in corners for PDB atoms not located in posivite axes, the atoms need to be moved. This effect should be
    reversed later.
 
    \param[in] pdbFile A pointer to the gemmi::Structure object as read in from the input file.
    \param[in] xMov How many angstroms should the atoms be moved along the x axis.
    \param[in] yMov How many angstroms should the atoms be moved along the y axis.
    \param[in] zMov How many angstroms should the atoms be moved along the z axis.
    \param[in] firstModel Should only the first, or all models be used?
 */
void ProSHADE_internal_mapManip::movePDBForMapCalc ( gemmi::Structure *pdbFile, proshade_single xMov, proshade_single yMov, proshade_single zMov, bool firstModel )
{
    //================================================ Use the first model, if it exists
    if ( pdbFile->models.size() > 0 )
    {
        //============================================ For each model
        for ( proshade_unsign sIt = 0; sIt < static_cast<proshade_unsign> ( pdbFile->models.size() ); sIt++ )
        {
            //======================================== Check if multiple models are allowed
            if ( firstModel && ( sIt != 0 ) ) { break; }
            
            //======================================== Get model
            gemmi::Model *model                       = &pdbFile->models.at(sIt);
            
            //======================================== For each chain
            for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( model->chains.size() ); mIt++ )
            {
                //==================================== Get chain
                gemmi::Chain *chain                   = &model->chains.at(mIt);
                
                //==================================== For each residue
                for ( proshade_unsign rIt = 0; rIt < static_cast<proshade_unsign> ( chain->residues.size() ); rIt++ )
                {
                    //================================ Get residue
                    gemmi::Residue *residue           = &chain->residues.at(rIt);
                    
                    //================================ For each atom
                    for ( proshade_unsign aIt = 0; aIt < static_cast<proshade_unsign> ( residue->atoms.size() ); aIt++ )
                    {
                        //============================ Get atom
                        gemmi::Atom *atom             = &residue->atoms.at(aIt);
                        
                        //============================ Move the atoms
                        atom->pos                     = gemmi::Position ( atom->pos.x + static_cast< proshade_double > ( xMov ), atom->pos.y + static_cast< proshade_double > ( yMov ), atom->pos.z + static_cast< proshade_double > ( zMov ) );
                    }
                }
            }

        }
    }
    else
    {
        std::stringstream hlpSS;
        hlpSS << "Found 0 models in input file " << pdbFile->name << ".\n                    : This suggests that the input co-ordinate file is\n                    : corrupted or mis-formatted.";
        throw ProSHADE_exception ( "Found no model in co-ordinate file.", "EP00050", __FILE__, __LINE__, __func__, hlpSS.str() );
    }
    
    //================================================ Done
    return ;
}

/*! \brief This function generates a theoretical map from co-ordinate input files.
 
    This function makes use of the Gemmi internal functionality to firstly create a properly sized cell object, then to check the input co-ordinate file for containing known elements as well as
    elements for which Gemmi knows the form factors. Then, this function proceeds to compute the F's using Cromer & Libermann method (from Gemmi) to finally compute the theoretical
    map using these. This map is then copied to the variable supplied in the second argument, presumably the ProSHADE internal map variable.
 
    \param[in] pdbFile A gemmi::Structure object read in from the input file.
    \param[in] map Pointer reference to a variable to save the map data.
    \param[in] requestedResolution Map resolution to which the map should be computed.
    \param[in] xCell The size of the map cell to be created.
    \param[in] yCell The size of the map cell to be created.
    \param[in] zCell The size of the map cell to be created.
    \param[in] xTo Pointer to variable where the map size along the x-axis in indices will be saved.
    \param[in] yTo Pointer to variable where the map size along the y-axis in indices will be saved.
    \param[in] zTo Pointer to variable where the map size along the z-axis in indices will be saved.
    \param[in] forceP1 Should the P1 spacegroup be forced?
    \param[in] firstModel Should only the first, or all models be used?
 
    \warning By default, this function will force the P1 spacegroup!
 */
void ProSHADE_internal_mapManip::generateMapFromPDB ( gemmi::Structure pdbFile, proshade_double*& map, proshade_single requestedResolution, proshade_single xCell, proshade_single yCell, proshade_single zCell, proshade_signed* xTo, proshade_signed* yTo, proshade_signed* zTo, bool forceP1, bool firstModel )
{
    //================================================ Set cell dimensions from the increased ranges (we need to add some space) and re-calculate cell properties
    if ( forceP1 ) { pdbFile.cell = gemmi::UnitCell(); }
    pdbFile.cell.a                                    = static_cast< proshade_double > ( xCell );
    pdbFile.cell.b                                    = static_cast< proshade_double > ( yCell );
    pdbFile.cell.c                                    = static_cast< proshade_double > ( zCell );
    pdbFile.cell.calculate_properties                 ( );

    //================================================ Get elements in Gemmi format
    std::string totElString;
    for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( pdbFile.models.size() ); mIt++ )
    {
        //============================================ Check if multiple models are allowed
        if ( firstModel && ( mIt != 0 ) )
        {
            std::stringstream hlpSS;
            hlpSS << "!!! ProSHADE WARNING !!! Found multiple models (" << pdbFile.models.size() << ") in input file " << pdbFile.name << ", while the settings state that only the first PDB file model should be used. If all models should be used, please supply ProSHADE with the \"-x\" option.";
            ProSHADE_internal_messages::printWarningMessage ( 0, hlpSS.str(), "WP00055" );
            break;
        }
        
        std::string hlpStr                            = pdbFile.models[mIt].present_elements ( ).to_string<char,std::char_traits<char>,std::allocator<char> >();
        totElString                                   = totElString + hlpStr;
    }
    std::bitset< static_cast< size_t > ( gemmi::El::END )> present_elems ( totElString );
    
    //================================================ Sanity checks
    if ( present_elems[static_cast<int> ( gemmi::El::X )] )
    {
        throw ProSHADE_exception ( "Found unknown element in input file.", "EP00051", __FILE__, __LINE__, __func__, "Gemmi library does not recognise some of the elements in\n                    : the co-ordinate file. Please check the file for not being\n                    : corrupted and containing standard elements." );
    }
    
    for ( proshade_unsign elIt = 0; elIt < static_cast<proshade_unsign> ( present_elems.size() ); elIt++ )
    {
      if ( present_elems[elIt] && !gemmi::IT92<double>::has ( static_cast<gemmi::El> ( elIt ) ) )
      {
          std::stringstream hlpSS;
          hlpSS << "Missing form factor for element " << element_name ( static_cast<gemmi::El> ( elIt ) );
          throw ProSHADE_exception ( hlpSS.str().c_str(), "EP00052", __FILE__, __LINE__, __func__, "Gemmi library does not have a form factor value for this\n                    : reported element. Please report this to the author." );
      }
    }
    
    //================================================ Compute the f's
    double wavelength                                 = 10.0;
    double energy                                     = gemmi::hc() / wavelength;
    
    //================================================ Create the density calculator object and fill it in
    gemmi::DensityCalculator<gemmi::IT92<double>, float> dencalc;
    
    dencalc.d_min                                     = static_cast< double > ( requestedResolution );
    for ( size_t elIt = 0; elIt < present_elems.size(); elIt++ ) { if ( present_elems[elIt] ) { dencalc.addends.set ( static_cast< gemmi::El > ( elIt ), static_cast< float > ( gemmi::cromer_liberman ( static_cast< int > ( elIt ), energy, nullptr ) ) ); } }
    dencalc.set_grid_cell_and_spacegroup              ( pdbFile );
    
    //================================================ Force P1 spacegroup
    if ( forceP1 ) { dencalc.grid.spacegroup          = &gemmi::get_spacegroup_p1(); }
    
    //================================================ Compute the theoretical map for each model
    dencalc.grid.data.clear                           ( );
    dencalc.grid.set_size_from_spacing                ( dencalc.d_min / ( 2.0 * dencalc.rate), true );
    for ( proshade_unsign mIt = 0; mIt < static_cast<proshade_unsign> ( pdbFile.models.size() ); mIt++ )
    {
        if ( firstModel && ( mIt != 0 ) ) { break; }
        dencalc.add_model_density_to_grid             ( pdbFile.models[mIt] );
        dencalc.grid.symmetrize                       ( [](float a, float b) { return a + b; } );
    }
    
    //================================================ Get the map
    const gemmi::Grid<float>& grid                    = dencalc.grid;
    
    //================================================ Save the map dimensions
   *xTo                                               = grid.nu;
   *yTo                                               = grid.nv;
   *zTo                                               = grid.nw;

    //================================================ Copy the gemmi::Grid to my map format
    map                                               = new proshade_double [(*xTo) * (*yTo) * (*zTo)];
    ProSHADE_internal_misc::checkMemoryAllocation     ( map, __FILE__, __LINE__, __func__ );

    proshade_signed arrPos                            = 0;
    for ( proshade_signed uIt = 0; uIt < (*xTo); uIt++ )
    {
        for ( proshade_signed vIt = 0; vIt < (*yTo); vIt++ )
        {
            for ( proshade_signed wIt = 0; wIt < (*zTo); wIt++ )
            {
                arrPos                                = wIt  + (*zTo) * ( vIt  + (*yTo) * uIt );
                map[arrPos]                           = static_cast< proshade_double > ( grid.get_value_q( static_cast< int > ( uIt ), static_cast< int > ( vIt ), static_cast< int > ( wIt ) ) );
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for moving map back to original PDB location by changing the indices.
 
    This function translates the map by changing the to and from index values so that the location of the map will be
    the same as the location of the original PDB file. This most likely cannot be done exactly as it allows only
    movement by increments of the sampling rate, but it is a decent start.
 
    \param[in] xMov Pointer to the NEGATIVE value by how many angstroms should the x axis be moved.
    \param[in] yMov Pointer to the NEGATIVE value by how many angstroms should the y axis be moved.
    \param[in] zMov Pointer to the NEGATIVE value by how many angstroms should the z axis be moved.
    \param[in] xAngs How many angstroms are there along the x dimension.
    \param[in] yAngs How many angstroms are there along the y dimension.
    \param[in] zAngs How many angstroms are there along the z dimension.
    \param[in] xFrom The initial index of the x dimension of the map.
    \param[in] xTo The terminal index of the x dimension of the map.
    \param[in] yFrom The initial index of the y dimension of the map.
    \param[in] yTo The terminal index of the y dimension of the map.
    \param[in] zFrom The initial index of the z dimension of the map.
    \param[in] zTo The terminal index of the z dimension of the map.
    \param[in] xOrigin The first value of the x axis index.
    \param[in] yOrigin The first value of the y axis index.
    \param[in] zOrigin The first value of the z axis index.
 */
void ProSHADE_internal_mapManip::moveMapByIndices ( proshade_single* xMov, proshade_single* yMov, proshade_single* zMov, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_signed* xFrom, proshade_signed* xTo, proshade_signed* yFrom, proshade_signed* yTo, proshade_signed* zFrom, proshade_signed* zTo, proshade_signed* xOrigin, proshade_signed* yOrigin, proshade_signed* zOrigin )
{
    //================================================ Compute movement in indices
    proshade_single xIndMove                          = std::floor ( -(*xMov) / ( xAngs / ( static_cast< proshade_single > ( *xTo ) - static_cast< proshade_single > ( *xFrom ) + 1.0f ) ) );
    proshade_single yIndMove                          = std::floor ( -(*yMov) / ( yAngs / ( static_cast< proshade_single > ( *yTo ) - static_cast< proshade_single > ( *yFrom ) + 1.0f ) ) );
    proshade_single zIndMove                          = std::floor ( -(*zMov) / ( zAngs / ( static_cast< proshade_single > ( *zTo ) - static_cast< proshade_single > ( *zFrom ) + 1.0f ) ) );
    
    //================================================ Set the movs to the remainder
    *xMov                                             =  -( *xMov ) - ( xIndMove * ( xAngs / ( static_cast< proshade_single > ( *xTo ) - static_cast< proshade_single > ( *xFrom ) + 1.0f ) ) );
    *yMov                                             =  -( *yMov ) - ( yIndMove * ( yAngs / ( static_cast< proshade_single > ( *yTo ) - static_cast< proshade_single > ( *yFrom ) + 1.0f ) ) );
    *zMov                                             =  -( *zMov ) - ( zIndMove * ( zAngs / ( static_cast< proshade_single > ( *zTo ) - static_cast< proshade_single > ( *zFrom ) + 1.0f ) ) );
    
    //================================================ Move indices by as much
    *xFrom                                           += static_cast< proshade_signed > ( xIndMove );
    *xTo                                             += static_cast< proshade_signed > ( xIndMove );
    *yFrom                                           += static_cast< proshade_signed > ( yIndMove );
    *yTo                                             += static_cast< proshade_signed > ( yIndMove );
    *zFrom                                           += static_cast< proshade_signed > ( zIndMove );
    *zTo                                             += static_cast< proshade_signed > ( zIndMove );
    
    //================================================ And set origin to reflect the changes
    *xOrigin                                          = *xFrom;
    *yOrigin                                          = *yFrom;
    *zOrigin                                          = *zFrom;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for moving map back to original PDB location by using Fourier transformation.
 
    This function translates the map by changing the phase information of the map Fourier transform and then
    doing the inverse Fourier. This allows for movement smaller than one index distance, but should not be
    used over long distances (could move out of boundaries and meet pariodicity problem such as map from
    different cell now being moved into this cell).
 
    \param[in] map A Reference Pointer to the map which should be shifted.
    \param[in] xMov The NEGATIVE value by how many angstroms should the x axis be moved.
    \param[in] yMov The NEGATIVE value by how many angstroms should the y axis be moved.
    \param[in] zMov The NEGATIVE value by how many angstroms should the z axis be moved.
    \param[in] xAngs How many angstroms are there along the x dimension.
    \param[in] yAngs How many angstroms are there along the y dimensionzAngs
    \param[in] zAngs How many angstroms are there along the z dimension.
    \param[in] xDim How many indices are there along the x dimension.
    \param[in] yDim How many indices are there along the y dimension.
    \param[in] zDim How many indices are there along the z dimension.
 */
void ProSHADE_internal_mapManip::moveMapByFourier ( proshade_double*& map, proshade_single xMov, proshade_single yMov, proshade_single zMov, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim )
{
    //================================================ Local variables initialisation
    proshade_unsign arrayPos                          = 0;
    
    //================================================ Create Fourier map variable
    fftw_complex *fCoeffs                             = new fftw_complex [xDim * yDim * zDim];
    fftw_complex *translatedMap                       = new fftw_complex [xDim * yDim * zDim];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( fCoeffs,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( translatedMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Create plans
    fftw_plan planForwardFourier                      = fftw_plan_dft_3d ( static_cast< int > ( xDim ), static_cast< int > ( yDim ), static_cast< int > ( zDim ), translatedMap, fCoeffs, FFTW_FORWARD,  FFTW_ESTIMATE );
    fftw_plan planBackwardFourier                     = fftw_plan_dft_3d ( static_cast< int > ( xDim ), static_cast< int > ( yDim ), static_cast< int > ( zDim ), fCoeffs, translatedMap, FFTW_BACKWARD, FFTW_ESTIMATE );
    
    //================================================ Copy map to complex format
    for ( proshade_unsign uIt = 0; uIt < static_cast< proshade_unsign > ( xDim ); uIt++ )
    {
        for ( proshade_unsign vIt = 0; vIt < static_cast< proshade_unsign > ( yDim ); vIt++ )
        {
            for ( proshade_unsign wIt = 0; wIt < static_cast< proshade_unsign > ( zDim ); wIt++ )
            {
                arrayPos                              = wIt + static_cast< proshade_unsign > ( zDim ) * ( vIt + static_cast< proshade_unsign > ( yDim ) * uIt );
                
                const FloatingPoint< proshade_double > lhs ( map[arrayPos] ), rhs ( map[arrayPos] );
                if ( lhs.AlmostEquals ( rhs ) )       { translatedMap[arrayPos][0] = map[arrayPos]; }
                else                                  { translatedMap[arrayPos][0] = 0.0; }
                translatedMap[arrayPos][1]            = 0.0;
            }
        }
    }
    
    //================================================ Compute Forward Fourier
    fftw_execute                                      ( planForwardFourier );
    
    //================================================ Shift the Fourier coefficients
    proshade_double *weight                           = nullptr;
    moveMapByFourierInReci                            ( fCoeffs, weight, xMov, yMov, zMov, xAngs, yAngs, zAngs, xDim, yDim, zDim );

    //================================================ Compute inverse Fourier
    fftw_execute                                      ( planBackwardFourier );
    
    //================================================ Copy back to map
    for ( proshade_unsign uIt = 0; uIt < static_cast< proshade_unsign > ( xDim ); uIt++ )
    {
        for ( proshade_unsign vIt = 0; vIt < static_cast< proshade_unsign > ( yDim ); vIt++ )
        {
            for ( proshade_unsign wIt = 0; wIt < static_cast< proshade_unsign > ( zDim ); wIt++ )
            {
                arrayPos                              = wIt + static_cast< proshade_unsign > ( zDim ) * ( vIt + static_cast< proshade_unsign > ( yDim ) * uIt );
                map[arrayPos]                         = translatedMap[arrayPos][0];
            }
        }
    }
    
    //================================================ Release memory
    fftw_destroy_plan                                 ( planForwardFourier );
    fftw_destroy_plan                                 ( planBackwardFourier );
    delete[] fCoeffs;
    delete[] translatedMap;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for moving map using Fourier coefficients in the reciprocal space only.
 
    This function translates the map by changing the phase information of the map Fourier transform. Contrary to the
    moveMapByFourier() function, it does not compute the inverse Fourier transform
 
    \param[in] coeffs A Reference Pointer to the Fourier decomposition coefficients of the map to be moved.
    \param[in] weights An array of weights to be applied to the shifted coefficients.
    \param[in] xMov The NEGATIVE value by how many angstroms should the x axis be moved.
    \param[in] yMov The NEGATIVE value by how many angstroms should the y axis be moved.
    \param[in] zMov The NEGATIVE value by how many angstroms should the z axis be moved.
    \param[in] xAngs How many angstroms are there along the x dimension.
    \param[in] yAngs How many angstroms are there along the y dimensionzAngs
    \param[in] zAngs How many angstroms are there along the z dimension.
    \param[in] xDim How many indices are there along the x dimension.
    \param[in] yDim How many indices are there along the y dimension.
    \param[in] zDim How many indices are there along the z dimension.
 */
void ProSHADE_internal_mapManip::moveMapByFourierInReci ( proshade_complex*& coeffs, proshade_double*& weights, proshade_single xMov, proshade_single yMov, proshade_single zMov, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim )
{
    //================================================ Local variables initialisation
    proshade_unsign arrayPos                          = 0;
    proshade_signed h, k, l;
    proshade_double real                              = 0.0;
    proshade_double imag                              = 0.0;
    proshade_double trCoeffReal, trCoeffImag;
    proshade_double normFactor                        = static_cast< proshade_double > ( xDim * yDim * zDim );
    proshade_double exponent                          = 0.0;
    proshade_double hlpArrReal;
    proshade_double hlpArrImag;
    
    //================================================ If no weights are given, set them to 1.0, otherwise use supplied
    proshade_double* wght                             = new proshade_double[xDim * yDim * zDim];
    ProSHADE_internal_misc::checkMemoryAllocation     ( wght, __FILE__, __LINE__, __func__ );
    if ( weights == nullptr ) { for ( size_t iter = 0; iter < static_cast< size_t > ( xDim * yDim * zDim ); iter++ ) { wght[iter] = 1.0; } }
    else                      { for ( size_t iter = 0; iter < static_cast< size_t > ( xDim * yDim * zDim ); iter++ ) { wght[iter] = weights[iter]; } }
    
    //================================================ Add the required shift
    for ( proshade_unsign uIt = 0; uIt < static_cast<proshade_unsign> ( xDim ); uIt++ )
    {
        for ( proshade_unsign vIt = 0; vIt < static_cast<proshade_unsign> ( yDim ); vIt++ )
        {
            for ( proshade_unsign wIt = 0; wIt < static_cast<proshade_unsign> ( zDim ); wIt++ )
            {
                //==================================== Var init
                arrayPos                              = wIt + static_cast< proshade_unsign > ( zDim ) * ( vIt + static_cast< proshade_unsign > ( yDim ) * uIt );
                real                                  = coeffs[arrayPos][0];
                imag                                  = coeffs[arrayPos][1];
                
                //==================================== Convert 0-max indices to HKL
                if ( uIt > static_cast< proshade_unsign > ( (xDim+1) / 2) ) { h = static_cast < proshade_signed > ( uIt ) - xDim; } else { h = static_cast < proshade_signed > ( uIt ); }
                if ( vIt > static_cast< proshade_unsign > ( (yDim+1) / 2) ) { k = static_cast < proshade_signed > ( vIt ) - yDim; } else { k = static_cast < proshade_signed > ( vIt ); }
                if ( wIt > static_cast< proshade_unsign > ( (zDim+1) / 2) ) { l = static_cast < proshade_signed > ( wIt ) - zDim; } else { l = static_cast < proshade_signed > ( wIt ); }
                
                //==================================== Get translation coefficient change
                exponent                              = ( ( ( static_cast <proshade_double> ( h ) / static_cast <proshade_double> ( xAngs ) ) * static_cast< proshade_double > ( -xMov ) ) +
                                                          ( ( static_cast <proshade_double> ( k ) / static_cast <proshade_double> ( yAngs ) ) * static_cast< proshade_double > ( -yMov ) ) +
                                                          ( ( static_cast <proshade_double> ( l ) / static_cast <proshade_double> ( zAngs ) ) * static_cast< proshade_double > ( -zMov ) ) ) * 2.0 * M_PI;
                        
                trCoeffReal                           = cos ( exponent );
                trCoeffImag                           = sin ( exponent );
                ProSHADE_internal_maths::complexMultiplication ( &real, &imag, &trCoeffReal, &trCoeffImag, &hlpArrReal, &hlpArrImag );
                
                //==================================== Save the translated coefficient value and apply weights
                coeffs[arrayPos][0]                   = ( hlpArrReal / normFactor ) * wght[arrayPos];
                coeffs[arrayPos][1]                   = ( hlpArrImag / normFactor ) * wght[arrayPos];
            }
        }
    }
    
    //================================================ Release weights
    delete[] wght;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for blurring/sharpening maps.
 
    This function takes the internal map representation information from the parameters as well as the internal map itself
    and proceeds to compute blue/sharpen the map by changing its overall B-factors by adding a specific value. If this value
    is negative, then the map is sharpened, while if the value is positive, the map is blurred.
 
    \param[in] map A Reference Pointer to the map which should be blurred/sharpened.
    \param[in] blurredMap A Reference Pointer to the variable which will store the modified map.
    \param[in] xDimS The number of indices along the x axis of the map.
    \param[in] yDimS The number of indices along the y axis of the map.
    \param[in] zDimS The number of indices along the z axis of the map.
    \param[in] xAngs The size of the x dimension of the map in angstroms.
    \param[in] yAngs The size of the y dimension of the map in angstroms.
    \param[in] zAngs The size of the z dimension of the map in angstroms.
    \param[in] blurringFactor The amount of B-factor change that should be applied to the map. (I.e. increasing the map overall B-factors by 20 will be done by supplying 20 as a value for this parameter).
 */
void ProSHADE_internal_mapManip::blurSharpenMap ( proshade_double*& map, proshade_double*& blurredMap, proshade_unsign xDimS, proshade_unsign yDimS, proshade_unsign zDimS, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_single blurringFactor )
{
    //================================================ Set local variables
    proshade_signed xDim                              = static_cast< proshade_signed > ( xDimS );
    proshade_signed yDim                              = static_cast< proshade_signed > ( yDimS );
    proshade_signed zDim                              = static_cast< proshade_signed > ( zDimS );
    proshade_double real, imag, S, mag, phase;
    proshade_signed h, k, l;
    proshade_unsign arrayPos                          = 0;
    proshade_double normFactor                        = static_cast<proshade_double> ( xDim * yDim * zDim );
    
    //================================================ Copy map for processing
    fftw_complex* mapCoeffs                           = new fftw_complex[xDim * yDim * zDim];
    fftw_complex* mapMask                             = new fftw_complex[xDim * yDim * zDim];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( mapCoeffs, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( mapMask,   __FILE__, __LINE__, __func__ );
    
    //================================================ Copy data to mask
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> (xDim * yDim * zDim); iter++ )
    {
        mapMask[iter][0]                              = map[iter];
        mapMask[iter][1]                              = 0.0;
    }
    
    //================================================ Prepare FFTW plans
    fftw_plan forward                                 = fftw_plan_dft_3d ( static_cast< int > ( xDim ), static_cast< int > ( yDim ), static_cast< int > ( zDim ), mapMask, mapCoeffs, FFTW_FORWARD,  FFTW_ESTIMATE );
    fftw_plan inverse                                 = fftw_plan_dft_3d ( static_cast< int > ( xDim ), static_cast< int > ( yDim ), static_cast< int > ( zDim ), mapCoeffs, mapMask, FFTW_BACKWARD, FFTW_ESTIMATE );
    
    //================================================ Run forward Fourier
    fftw_execute                                      ( forward );
    
    //================================================ Blur the coeffs
    for ( proshade_unsign uIt = 0; uIt < static_cast<proshade_unsign> ( xDim ); uIt++ )
    {
        for ( proshade_unsign vIt = 0; vIt < static_cast<proshade_unsign> ( yDim ); vIt++ )
        {
            for ( proshade_unsign wIt = 0; wIt < static_cast<proshade_unsign> ( zDim ); wIt++ )
            {
                //==================================== Var init
                arrayPos                              = wIt + static_cast< proshade_unsign > ( zDim ) * ( vIt + static_cast< proshade_unsign > ( yDim ) * uIt );
                real                                  = mapCoeffs[arrayPos][0];
                imag                                  = mapCoeffs[arrayPos][1];
                
                //==================================== Convert to HKL
                if ( uIt > static_cast< proshade_unsign > ( (xDim+1) / 2) ) { h = static_cast < proshade_signed > ( uIt ) - xDim; } else { h = static_cast < proshade_signed > ( uIt ); }
                if ( vIt > static_cast< proshade_unsign > ( (yDim+1) / 2) ) { k = static_cast < proshade_signed > ( vIt ) - yDim; } else { k = static_cast < proshade_signed > ( vIt ); }
                if ( wIt > static_cast< proshade_unsign > ( (zDim+1) / 2) ) { l = static_cast < proshade_signed > ( wIt ) - zDim; } else { l = static_cast < proshade_signed > ( wIt ); }
                
                //====================================Get magnitude and phase with mask parameters
                S                                     = ( pow( static_cast< proshade_double > ( h ) / static_cast< proshade_double > ( xAngs ), 2.0 ) +
                                                          pow( static_cast< proshade_double > ( k ) / static_cast< proshade_double > ( yAngs ), 2.0 ) +
                                                          pow( static_cast< proshade_double > ( l ) / static_cast< proshade_double > ( zAngs ), 2.0 ) );
                mag                                   = std::sqrt ( (real*real) + (imag*imag) ) * std::exp ( - ( ( static_cast< proshade_double > ( blurringFactor ) * S ) / 4.0 ) );
                phase                                 = std::atan2 ( imag, real );
                
                //==================================== Save the mask data
                mapCoeffs[arrayPos][0]                = ( mag * cos(phase) ) / normFactor;
                mapCoeffs[arrayPos][1]                = ( mag * sin(phase) ) / normFactor;
            }
        }
    }
    
    //================================================ Run inverse Fourier
    fftw_execute                                      ( inverse );
    
    //================================================ Save the results
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> (xDim * yDim * zDim); iter++ )
    {
        blurredMap[iter]                              = mapMask[iter][0];
    }
    
    //================================================ Release memory
    delete[] mapMask;
    delete[] mapCoeffs;
    
    //================================================ Delete FFTW plans
    fftw_destroy_plan                                 ( forward );
    fftw_destroy_plan                                 ( inverse );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for computing mask from blurred map.
 
    This function takes a blurred map and proceeds to compute the cut-off threshold for masking. It then applies
    this mask to the blurred map and applies the masking to the second input map. The assumption is that this
    second map is the map before blurring, as non-zero mask points will not be changed. Careful about this!!! It,
    however, does not output the mask.
 
    \param[in] blurMap A Reference Pointer to the map which has been blurred for mask computation.
    \param[in] outMap A Reference Pointer to the map which will be masked - this map should be the map before blurring.
    \param[in] xDim The number of indices along the x axis of the map.
    \param[in] yDim The number of indices along the y axis of the map.
    \param[in] zDim The number of indices along the z axis of the map.
    \param[in] noIQRs The number of inter-quartile ranges from the median which should be used to compute the threshold for masking.
 */
void ProSHADE_internal_mapManip::getMaskFromBlurr ( proshade_double*& blurMap, proshade_double*& outMap, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_single noIQRs )
{
    //================================================ Initialise vector for map data
    std::vector<proshade_double> mapVals ( xDim * yDim * zDim, 0.0 );
    
    //================================================ Save map values in vector
    for ( proshade_unsign iter = 0; iter < ( xDim * yDim * zDim ); iter++ )
    {
        mapVals.at(iter)                              = blurMap[iter];
    }
    
    //================================================ Find median and IQRs
    proshade_double* medAndIQR                        = new proshade_double[2];
    ProSHADE_internal_maths::vectorMedianAndIQR ( &mapVals, medAndIQR );
    
    //================================================ Find the threshold
    proshade_double maskThreshold                     = medAndIQR[0] + ( medAndIQR[1] * static_cast<proshade_double> ( noIQRs ) );
    
    //================================================ Apply threshold
    for ( proshade_unsign iter = 0; iter < ( xDim * yDim * zDim ); iter++ )
    {
        if ( blurMap[iter] < maskThreshold )
        {
            outMap[iter]                              = 0.0;
            blurMap[iter]                             = 0.0;
        }
    }
    
    //================================================ Release vector values
    mapVals.clear                                     ( );
    
    //================================================ Release memory
    delete[] medAndIQR;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for finding the map boundaries enclosing positive only values.
 
    This function takes a map and searches for minimum and maximum map indices in all three dimensions, which enclose
    all the non-zero map values.
 
    \param[in] map A pointer to the map for which the bounds are to be found.
    \param[in] xDim The number of indices along the x axis of the map.
    \param[in] yDim The number of indices along the y axis of the map.
    \param[in] zDim The number of indices along the z axis of the map.
    \param[in] ret A proshade_unsign pointer to array of 6 to which the results will be saved (0 = minX; 1 = maxX; 2 = minY; 3 = maxY; 4 - minZ; 5 = maxZ).
 */
void ProSHADE_internal_mapManip::getNonZeroBounds ( proshade_double* map, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim, proshade_signed*& ret )
{
    //================================================ Initialise local variables
    proshade_signed arrayPos                          = 0;
    
    //================================================ Initialise result variable
    ret[0]                                            = xDim;
    ret[1]                                            = 0;
    ret[2]                                            = yDim;
    ret[3]                                            = 0;
    ret[4]                                            = zDim;
    ret[5]                                            = 0;
    
    //================================================ Iterate through map and check bounds
    for ( proshade_signed xIt = 0; xIt < xDim; xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < yDim; yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < zDim; zIt++ )
            {
                //==================================== Var init
                arrayPos                              = zIt + zDim * ( yIt + yDim * xIt );
                
                //==================================== Check bounds
                if ( map[arrayPos] > 0.001 )
                {
                    if ( xIt < ret[0] )               { ret[0] = xIt; }
                    if ( xIt > ret[1] )               { ret[1] = xIt; }
                    if ( yIt < ret[2] )               { ret[2] = yIt; }
                    if ( yIt > ret[3] )               { ret[3] = yIt; }
                    if ( zIt < ret[4] )               { ret[4] = zIt; }
                    if ( zIt > ret[5] )               { ret[5] = zIt; }
                }
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function takes a set of bounds and adds a given number of Angstroms to them.
 
    This function adds an extra distance to a set of bounds as given by the number of angstroms to be added. If the resulting bounds
    are larger than the maximum map bounds, the maxima and minima are returned instead, thus partially ignoring the requested extra
    space. !!! If the resulting index ranges would be odd, extra one idice is added to make them even. !!!
 
    \param[in] xDim The number of indices along the x axis of the map.
    \param[in] yDim The number of indices along the y axis of the map.
    \param[in] zDim The number of indices along the z axis of the map.
    \param[in] xAngs The size of the x dimension of the map in angstroms.
    \param[in] yAngs The size of the y dimension of the map in angstroms.
    \param[in] zAngs The size of the z dimension of the map in angstroms.
    \param[in] bounds A proshade_unsign pointer reference to array of 6 which has the current bounds (0 = minX; 1 = maxX; 2 = minY; 3 = maxY; 4 - minZ; 5 = maxZ).
    \param[in] extraSpace The number of angstroms to be added to each dimension (both to start and to end).
 */
void ProSHADE_internal_mapManip::addExtraBoundSpace ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_signed*& bounds, proshade_single extraSpace )
{
    //================================================ Convert angstrom distance to indices
    proshade_signed xExtraInds                        = ProSHADE_internal_mapManip::myRound ( extraSpace / (  xAngs / static_cast<proshade_single> ( xDim ) ) );
    proshade_signed yExtraInds                        = ProSHADE_internal_mapManip::myRound ( extraSpace / (  yAngs / static_cast<proshade_single> ( yDim ) ) );
    proshade_signed zExtraInds                        = ProSHADE_internal_mapManip::myRound ( extraSpace / (  zAngs / static_cast<proshade_single> ( zDim ) ) );
    
    //=============================================== Changed bounds even if exceeding physical map - this will be deal with in the map creation part
    bounds[0]                                         = bounds[0] - xExtraInds;
    bounds[1]                                         = bounds[1] + xExtraInds;
    bounds[2]                                         = bounds[2] - yExtraInds;
    bounds[3]                                         = bounds[3] + yExtraInds;
    bounds[4]                                         = bounds[4] - zExtraInds;
    bounds[5]                                         = bounds[5] + zExtraInds;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function re-samples a map to conform to given resolution using tri-linear interpolation.
 
    This function takes a map and resolution value and it proceeds to create a new map, which has sampling resolution/2 and is
    large enough to contain the original map. It then proceeds to interpolate the new map values from the old map values, re-writing
    the old map once interpolation is done.
 
    \param[in] map A Reference Pointer to the map for which the bounds are to be found.
    \param[in] resolution The required resolution value.
    \param[in] xDimS The number of indices along the x axis of the map.
    \param[in] yDimS The number of indices along the y axis of the map.
    \param[in] zDimS The number of indices along the z axis of the map.
    \param[in] xAngs The size of the x dimension of the map in angstroms.
    \param[in] yAngs The size of the y dimension of the map in angstroms.
    \param[in] zAngs The size of the z dimension of the map in angstroms.
 \param[in] corrs Pointer reference to proshade_single array of 6 values with the following meaning: 0 = xAdd; 1 = yAdd; 2 = zAdd; 3 = newXAng; 4 = newYAng;  5 = newZAng
 */
void ProSHADE_internal_mapManip::reSampleMapToResolutionTrilinear ( proshade_double*& map, proshade_single resolution, proshade_unsign xDimS, proshade_unsign yDimS, proshade_unsign zDimS, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_single*& corrs )
{
    //================================================ Sanity check - the resolution needs to be set
    if ( resolution <= 0.0f )
    {
        throw ProSHADE_exception ( "Requested resolution not set for map re-sampling.", "EM00015", __FILE__, __LINE__, __func__, "There is no resolution value set, but map re-sampling to\n                    : this unset resolution value is required. This error\n                    : occurs when a task with no resolution requirement is\n                    : requested on a map data and the map resolution change is\n                    : set to \'on\'. Either supply a resolution value, or do not\n                    : re-sample the map." );
    }
    
    //================================================ Initialise local variables
    proshade_signed xDim                              = static_cast<proshade_signed> ( xDimS );
    proshade_signed yDim                              = static_cast<proshade_signed> ( yDimS );
    proshade_signed zDim                              = static_cast<proshade_signed> ( zDimS );
    proshade_single oldXSample                        = ( xAngs / static_cast<proshade_single> ( xDim ) );
    proshade_single oldYSample                        = ( yAngs / static_cast<proshade_single> ( yDim ) );
    proshade_single oldZSample                        = ( zAngs / static_cast<proshade_single> ( zDim ) );
    proshade_single newXSample                        = static_cast< proshade_single  > ( resolution / 2.0f );
    proshade_single newYSample                        = static_cast< proshade_single  > ( resolution / 2.0f );
    proshade_single newZSample                        = static_cast< proshade_single  > ( resolution / 2.0f );
    
    //================================================ Compute required grid size
    proshade_signed newXDim                           = static_cast<proshade_signed> ( std::ceil ( xAngs / newXSample ) );
    proshade_signed newYDim                           = static_cast<proshade_signed> ( std::ceil ( yAngs / newYSample ) );
    proshade_signed newZDim                           = static_cast<proshade_signed> ( std::ceil ( zAngs / newZSample ) );
    
    //================================================ Create a new map variable
    proshade_double* newMap                           = new proshade_double [newXDim * newYDim * newZDim];
    
    //================================================ For each new map point
    proshade_signed xBottom = 0, xTop, yBottom = 0, yTop, zBottom = 0, zTop, oldMapIndex, newMapIndex;
    std::vector<proshade_double> c000                 = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c001                 = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c010                 = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c011                 = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c100                 = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c101                 = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c110                 = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c111                 = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c00                  = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c01                  = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c10                  = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c11                  = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c0                   = std::vector<proshade_double> ( 4, 0.0 );
    std::vector<proshade_double> c1                   = std::vector<proshade_double> ( 4, 0.0 );
    proshade_double xRelative, yRelative, zRelative;
    
    for ( proshade_signed xIt = 0; xIt < newXDim; xIt++ )
    {
        for ( proshade_signed yIt = 0; yIt < newYDim; yIt++ )
        {
            for ( proshade_signed zIt = 0; zIt < newZDim; zIt++ )
            {
                //==================================== Get this point's index
                newMapIndex                           = zIt + newZDim * ( yIt + newYDim * xIt );
                
                //==================================== Find this points bottom and top positions in the old map (including periodicity)
                for ( proshade_signed ox = 0; ox < ( static_cast< proshade_signed > ( xDimS ) - 1 ); ox++ ) { if ( ( ( static_cast< proshade_single > ( xIt ) * newXSample ) >= ( static_cast< proshade_single > ( ox ) * oldXSample ) ) && ( ( static_cast< proshade_single > ( xIt ) * newXSample ) <= ( ( static_cast< proshade_single > ( ox ) + 1 ) * oldXSample ) ) ) { xBottom = ox; break; } }
                for ( proshade_signed oy = 0; oy < ( static_cast< proshade_signed > ( yDimS ) - 1 ); oy++ ) { if ( ( ( static_cast< proshade_single > ( yIt ) * newYSample ) >= ( static_cast< proshade_single > ( oy ) * oldYSample ) ) && ( ( static_cast< proshade_single > ( yIt ) * newYSample ) <= ( ( static_cast< proshade_single > ( oy ) + 1 ) * oldYSample ) ) ) { yBottom = oy; break; } }
                for ( proshade_signed oz = 0; oz < ( static_cast< proshade_signed > ( zDimS ) - 1 ); oz++ ) { if ( ( ( static_cast< proshade_single > ( zIt ) * newZSample ) >= ( static_cast< proshade_single > ( oz ) * oldZSample ) ) && ( ( static_cast< proshade_single > ( zIt ) * newZSample ) <= ( ( static_cast< proshade_single > ( oz ) + 1 ) * oldZSample ) ) ) { zBottom = oz; break; } }
                xTop                                  = xBottom + 1;
                yTop                                  = yBottom + 1;
                zTop                                  = zBottom + 1;

                //==================================== Find the surrounding point's values from the original map
                oldMapIndex                           = zBottom + static_cast< proshade_signed > ( zDimS ) * ( yBottom + static_cast< proshade_signed > ( yDimS ) * xBottom );
                c000.at(0)                            = static_cast<proshade_double> ( xBottom ) * static_cast<proshade_double> ( oldXSample );
                c000.at(1)                            = static_cast<proshade_double> ( yBottom ) * static_cast<proshade_double> ( oldYSample );
                c000.at(2)                            = static_cast<proshade_double> ( zBottom ) * static_cast<proshade_double> ( oldZSample );
                c000.at(3)                            = static_cast<proshade_double> ( map[oldMapIndex] );
                
                oldMapIndex                           = zTop    + static_cast< proshade_signed > ( zDimS ) * ( yBottom + static_cast< proshade_signed > ( yDimS ) * xBottom );
                c001.at(0)                            = static_cast<proshade_double> ( xBottom ) * static_cast<proshade_double> ( oldXSample );
                c001.at(1)                            = static_cast<proshade_double> ( yBottom ) * static_cast<proshade_double> ( oldYSample );
                c001.at(2)                            = static_cast<proshade_double> ( zTop    ) * static_cast<proshade_double> ( oldZSample );
                c001.at(3)                            = static_cast<proshade_double> ( map[oldMapIndex] );
                
                oldMapIndex                           = zBottom + static_cast< proshade_signed > ( zDimS ) * ( yTop    + static_cast< proshade_signed > ( yDimS ) * xBottom );
                c010.at(0)                            = static_cast<proshade_double> ( xBottom ) * static_cast<proshade_double> ( oldXSample );
                c010.at(1)                            = static_cast<proshade_double> ( yTop    ) * static_cast<proshade_double> ( oldYSample );
                c010.at(2)                            = static_cast<proshade_double> ( zBottom ) * static_cast<proshade_double> ( oldZSample );
                c010.at(3)                            = static_cast<proshade_double> ( map[oldMapIndex] );
                
                oldMapIndex                           = zTop    + static_cast< proshade_signed > ( zDimS ) * ( yTop    + static_cast< proshade_signed > ( yDimS ) * xBottom );
                c011.at(0)                            = static_cast<proshade_double> ( xBottom ) * static_cast<proshade_double> ( oldXSample );
                c011.at(1)                            = static_cast<proshade_double> ( yTop    ) * static_cast<proshade_double> ( oldYSample );
                c011.at(2)                            = static_cast<proshade_double> ( zTop    ) * static_cast<proshade_double> ( oldZSample );
                c011.at(3)                            = static_cast<proshade_double> ( map[oldMapIndex] );
                
                oldMapIndex                           = zBottom + static_cast< proshade_signed > ( zDimS ) * ( yBottom + static_cast< proshade_signed > ( yDimS ) * xTop    );
                c100.at(0)                            = static_cast<proshade_double> ( xTop    ) * static_cast<proshade_double> ( oldXSample );
                c100.at(1)                            = static_cast<proshade_double> ( yBottom ) * static_cast<proshade_double> ( oldYSample );
                c100.at(2)                            = static_cast<proshade_double> ( zBottom ) * static_cast<proshade_double> ( oldZSample );
                c100.at(3)                            = static_cast<proshade_double> ( map[oldMapIndex] );
                
                oldMapIndex                           = zTop    + static_cast< proshade_signed > ( zDimS ) * ( yBottom + static_cast< proshade_signed > ( yDimS ) * xTop    );
                c101.at(0)                            = static_cast<proshade_double> ( xTop    ) * static_cast<proshade_double> ( oldXSample );
                c101.at(1)                            = static_cast<proshade_double> ( yBottom ) * static_cast<proshade_double> ( oldYSample );
                c101.at(2)                            = static_cast<proshade_double> ( zTop    ) * static_cast<proshade_double> ( oldZSample );
                c101.at(3)                            = static_cast<proshade_double> ( map[oldMapIndex] );
                
                oldMapIndex                           = zBottom + static_cast< proshade_signed > ( zDimS ) * ( yTop    + static_cast< proshade_signed > ( yDimS ) * xTop    );
                c110.at(0)                            = static_cast<proshade_double> ( xTop    ) * static_cast<proshade_double> ( oldXSample );
                c110.at(1)                            = static_cast<proshade_double> ( yTop    ) * static_cast<proshade_double> ( oldYSample );
                c110.at(2)                            = static_cast<proshade_double> ( zBottom ) * static_cast<proshade_double> ( oldZSample );
                c110.at(3)                            = static_cast<proshade_double> ( map[oldMapIndex] );
                
                oldMapIndex                           = zTop    + static_cast< proshade_signed > ( zDimS ) * ( yTop    + static_cast< proshade_signed > ( yDimS ) * xTop    );
                c111.at(0)                            = static_cast<proshade_double> ( xTop    ) * static_cast<proshade_double> ( oldXSample );
                c111.at(1)                            = static_cast<proshade_double> ( yTop    ) * static_cast<proshade_double> ( oldYSample );
                c111.at(2)                            = static_cast<proshade_double> ( zTop    ) * static_cast<proshade_double> ( oldZSample );
                c111.at(3)                            = static_cast<proshade_double> ( map[oldMapIndex] );
                
                //==================================== Interpolate to the new grid along X
                xRelative                             = ( ( static_cast<proshade_double> ( xIt ) * static_cast<proshade_double> ( newXSample ) ) - ( static_cast<proshade_double> ( xBottom ) * static_cast<proshade_double> ( oldXSample ) ) ) / ( ( static_cast<proshade_double> ( xTop ) * static_cast<proshade_double> ( oldXSample ) ) - ( static_cast<proshade_double> ( xBottom ) * static_cast<proshade_double> ( oldXSample ) ) );
                
                //==================================== Interpolate for the less less point
                c00.at(0)                             = ( static_cast< proshade_double > ( newXSample ) * xRelative ) + c000.at(0);
                c00.at(1)                             = c000.at(1);
                c00.at(2)                             = c000.at(2);
                c00.at(3)                             = ( c000.at(3) * ( 1.0 - xRelative ) ) + ( c100.at(3) * xRelative );
                
                //==================================== Interpolate for the less more point
                c01.at(0)                             = ( static_cast< proshade_double > ( newXSample ) * xRelative ) + c001.at(0);
                c01.at(1)                             = c001.at(1);
                c01.at(2)                             = c001.at(2);
                c01.at(3)                             = ( c001.at(3) * ( 1.0 - xRelative ) ) + ( c101.at(3) * xRelative );
                
                //==================================== Interpolate for the more less point
                c10.at(0)                             = ( static_cast< proshade_double > ( newXSample ) * xRelative ) + c010.at(0);
                c10.at(1)                             = c010.at(1);
                c10.at(2)                             = c010.at(2);
                c10.at(3)                             = ( c010.at(3) * ( 1.0 - xRelative ) ) + ( c110.at(3) * xRelative );
                
                //==================================== Interpolate for the more more point
                c11.at(0)                             = ( static_cast< proshade_double > ( newXSample ) * xRelative ) + c011.at(0);
                c11.at(1)                             = c011.at(1);
                c11.at(2)                             = c011.at(2);
                c11.at(3)                             = ( c011.at(3) * ( 1.0 - xRelative ) ) + ( c111.at(3) * xRelative );
                
                //==================================== Interpolate to the new grid along Y
                yRelative                             =  ( ( static_cast<proshade_double> ( yIt ) * static_cast<proshade_double> ( newYSample ) ) - ( static_cast<proshade_double> ( yBottom ) * static_cast<proshade_double> ( oldYSample ) ) ) / ( ( static_cast<proshade_double> ( yTop ) * static_cast<proshade_double> ( oldYSample ) ) - ( static_cast<proshade_double> ( yBottom ) * static_cast<proshade_double> ( oldYSample ) ) );
                
                //==================================== Interpolate for the less point
                c0.at(0)                              = c00.at(0);
                c0.at(1)                              = ( static_cast< proshade_double > ( newYSample ) * yRelative ) + c00.at(1);
                c0.at(2)                              = c00.at(2);
                c0.at(3)                              = ( c00.at(3) * ( 1.0 - yRelative ) ) + ( c10.at(3) * yRelative );
                
                //==================================== Interpolate for the more point
                c1.at(0)                              = c01.at(0);
                c1.at(1)                              = ( static_cast< proshade_double > ( newYSample ) * yRelative ) + c01.at(1);
                c1.at(2)                              = c01.at(2);
                c1.at(3)                              = ( c01.at(3) * ( 1.0 - yRelative ) ) + ( c11.at(3) * yRelative );
                
                //==================================== Interpolate to the new grid along Z
                zRelative                             = ( ( static_cast<proshade_double> ( zIt ) * static_cast< proshade_double > ( newZSample ) ) - ( static_cast<proshade_double> ( zBottom ) * static_cast<proshade_double> ( oldZSample ) ) ) / static_cast< proshade_double > ( ( static_cast<proshade_double> ( zTop ) * static_cast<proshade_double> ( oldZSample ) ) - ( static_cast<proshade_double> ( zBottom ) * static_cast<proshade_double> ( oldZSample ) ) );
                newMap[newMapIndex]                   = ( c0.at(3) * ( 1.0 - zRelative ) ) + ( c1.at(3) * zRelative );
            }
        }
    }
    
    //================================================ Delete old map and allocate new memory
    delete[] map;
    map                                               = new proshade_double [newXDim * newYDim * newZDim];
    
    //================================================ Copy map
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( newXDim * newYDim * newZDim ); iter++ )
    {
        map[iter]                                     = newMap[iter];
    }
    
    //================================================ Release memory
    delete[] newMap;
    
    //================================================ Define change in indices and return it
    corrs[0]                                          = static_cast< proshade_single > ( newXDim - xDim );
    corrs[1]                                          = static_cast< proshade_single > ( newYDim - yDim );
    corrs[2]                                          = static_cast< proshade_single > ( newZDim - zDim );
    corrs[3]                                          = static_cast< proshade_single > ( newXDim ) * static_cast< proshade_single > ( newXSample );
    corrs[4]                                          = static_cast< proshade_single > ( newYDim ) * static_cast< proshade_single > ( newYSample );
    corrs[5]                                          = static_cast< proshade_single > ( newZDim ) * static_cast< proshade_single > ( newZSample );
    
    //======================================== Done
    return ;
    
}

/*! \brief This function re-samples a map to conform to given resolution using Fourier.
 
    This function re-samples the internal map to a given resolutution by removing or zero-padding the Fourier (reciprocal space) coefficients and computing the inverse
    Fourier transform. This is the default option for map re-sampling, should it be required by the user.
 
    \param[in] map A Reference Pointer to the map for which the bounds are to be found.
    \param[in] resolution The required resolution value.
    \param[in] xDimS The number of indices along the x axis of the map.
    \param[in] yDimS The number of indices along the y axis of the map.
    \param[in] zDimS The number of indices along the z axis of the map.
    \param[in] xAngs The size of the x dimension of the map in angstroms.
    \param[in] yAngs The size of the y dimension of the map in angstroms.
    \param[in] zAngs The size of the z dimension of the map in angstroms.
    \param[in] corrs Pointer reference to proshade_single array of 6 values with the following meaning: 0 = xAdd; 1 = yAdd; 2 = zAdd; 3 = newXAng; 4 = newYAng;  5 = newZAng
    \param[in] xFrom The initial index of the x dimension of the map.
    \param[in] xTo The terminal index of the x dimension of the map.
    \param[in] yFrom The initial index of the y dimension of the map.
    \param[in] yTo The terminal index of the y dimension of the map.
    \param[in] zFrom The initial index of the z dimension of the map.
    \param[in] zTo The terminal index of the z dimension of the map.
    \param[in] xOrigin The first value of the x axis index.
    \param[in] yOrigin The first value of the y axis index.
    \param[in] zOrigin The first value of the z axis index.
 */
void ProSHADE_internal_mapManip::reSampleMapToResolutionFourier ( proshade_double*& map, proshade_single resolution, proshade_unsign xDimS, proshade_unsign yDimS, proshade_unsign zDimS, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_single*& corrs )
{
    //================================================ Sanity check - the resolution needs to be set
    if ( resolution <= 0.0f )
    {
        throw ProSHADE_exception ( "Requested resolution not set for map re-sampling.", "EM00015", __FILE__, __LINE__, __func__, "There is no resolution value set, but map re-sampling to\n                    : this unset resolution value is required. This error\n                    : occurs when a task with no resolution requirement is\n                    : requested on a map data and the map resolution change is\n                    : set to \'on\'. Either supply a resolution value, or do not\n                    : re-sample the map." );
    }
    
    //================================================ Initialise variables
    proshade_unsign newXDim                           = static_cast<proshade_unsign> ( ProSHADE_internal_mapManip::myRound ( xAngs / ( resolution / 2.0f ) ) );
    proshade_unsign newYDim                           = static_cast<proshade_unsign> ( ProSHADE_internal_mapManip::myRound ( yAngs / ( resolution / 2.0f ) ) );
    proshade_unsign newZDim                           = static_cast<proshade_unsign> ( ProSHADE_internal_mapManip::myRound ( zAngs / ( resolution / 2.0f ) ) );
    
    if ( newXDim % 2 != 0 ) { newXDim += 1; }
    if ( newYDim % 2 != 0 ) { newYDim += 1; }
    if ( newZDim % 2 != 0 ) { newZDim += 1; }
 
    proshade_signed preXChange, preYChange, preZChange;
    if ( ( xDimS % 2 ) == 0 ) { preXChange = static_cast< proshade_signed  > ( std::ceil  ( ( static_cast<proshade_signed> ( xDimS ) - static_cast<proshade_signed> ( newXDim ) ) / 2 ) ); }
    else                      { preXChange = static_cast< proshade_signed  > ( std::floor ( ( static_cast<proshade_signed> ( xDimS ) - static_cast<proshade_signed> ( newXDim ) ) / 2 ) ); }
    if ( ( yDimS % 2 ) == 0 ) { preYChange = static_cast< proshade_signed  > ( std::ceil  ( ( static_cast<proshade_signed> ( yDimS ) - static_cast<proshade_signed> ( newYDim ) ) / 2 ) ); }
    else                      { preYChange = static_cast< proshade_signed  > ( std::floor ( ( static_cast<proshade_signed> ( yDimS ) - static_cast<proshade_signed> ( newYDim ) ) / 2 ) ); }
    if ( ( zDimS % 2 ) == 0 ) { preZChange = static_cast< proshade_signed  > ( std::ceil  ( ( static_cast<proshade_signed> ( zDimS ) - static_cast<proshade_signed> ( newZDim ) ) / 2 ) ); }
    else                      { preZChange = static_cast< proshade_signed  > ( std::floor ( ( static_cast<proshade_signed> ( zDimS ) - static_cast<proshade_signed> ( newZDim ) ) / 2 ) ); }
    
    proshade_signed postXChange                       = static_cast<proshade_signed> ( xDimS ) - ( preXChange + static_cast<proshade_signed> ( newXDim ) );
    proshade_signed postYChange                       = static_cast<proshade_signed> ( yDimS ) - ( preYChange + static_cast<proshade_signed> ( newYDim ) );
    proshade_signed postZChange                       = static_cast<proshade_signed> ( zDimS ) - ( preZChange + static_cast<proshade_signed> ( newZDim ) );
    
    proshade_unsign origSizeArr = 0, newSizeArr = 0;
    proshade_double normFactor                        = static_cast<proshade_double> ( xDimS * yDimS * zDimS );
    
    //================================================ Manage memory
    fftw_complex *origMap, *fCoeffs, *newFCoeffs, *newMap;
    fftw_plan planForwardFourier, planBackwardRescaledFourier;
    allocateResolutionFourierMemory                   ( origMap, fCoeffs, newFCoeffs, newMap, planForwardFourier, planBackwardRescaledFourier,
                                                        xDimS, yDimS, zDimS, newXDim, newYDim, newZDim );
    
    //================================================ Fill maps with data and zeroes
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( xDimS * yDimS * zDimS ); iter++ ) { origMap[iter][0] = map[iter]; origMap[iter][1] = 0.0; }
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( newXDim * newYDim * newZDim ); iter++ ) { newFCoeffs[iter][0] = 0.0; newFCoeffs[iter][1] = 0.0; }
    
    //================================================ Get the Fourier coeffs
    fftw_execute                                      ( planForwardFourier );
    
    //================================================ Change the order of Fourier coefficients
    changeFourierOrder                                ( fCoeffs, static_cast< proshade_signed > ( xDimS ), static_cast< proshade_signed > ( yDimS ), static_cast< proshade_signed > ( zDimS ), true );
    
    //================================================ Re-sample the coefficients by removing high frequencies or adding these with 0 values
    for ( proshade_unsign xIt = 0; xIt < newXDim; xIt++ )
    {
        for ( proshade_unsign yIt = 0; yIt < newYDim; yIt++ )
        {
            for ( proshade_unsign zIt = 0; zIt < newZDim; zIt++ )
            {
                //==================================== Find the array positions
                origSizeArr                           = ( ( zIt + static_cast< proshade_unsign > ( preZChange ) ) + zDimS   *
                                                        ( ( yIt + static_cast< proshade_unsign > ( preYChange ) ) + yDimS   *
                                                          ( xIt + static_cast< proshade_unsign > ( preXChange ) ) ) );
                newSizeArr                            = zIt                + newZDim * ( yIt                + newYDim * xIt                );
                
                //==================================== If original coefficient for this new coefficient position exists, copy
                if ( ( ( -1  < static_cast< proshade_signed > ( xIt )   + preXChange ) && ( -1  < static_cast<proshade_signed> ( yIt )   + preYChange ) && ( -1  < static_cast<proshade_signed> ( zIt )   + preZChange ) ) &&
                     ( ( xIt < newXDim + static_cast<proshade_unsign> ( postXChange ) ) && ( yIt < newYDim + static_cast<proshade_unsign> ( postYChange ) ) && ( zIt < newZDim + static_cast<proshade_unsign> ( postZChange ) ) ) )
                {
                    //================================ Copy the Fourier coeff
                    newFCoeffs[newSizeArr][0]         = fCoeffs[origSizeArr][0] / normFactor;
                    newFCoeffs[newSizeArr][1]         = fCoeffs[origSizeArr][1] / normFactor;
                }
            }
        }
    }
    
    //================================================ Change the order of the re-sampled Fourier coefficients
    changeFourierOrder                                ( newFCoeffs, static_cast< proshade_signed > ( newXDim ), static_cast< proshade_signed > ( newYDim ), static_cast< proshade_signed > ( newZDim ), false );

    //================================================ Get the new map from the re-sized Fourier coefficients
    fftw_execute                                      ( planBackwardRescaledFourier );

    //================================================ Delete the old map and create a new, re-sized one. Then copy the new map values into this new map memory.
    delete map;
    map                                               = new proshade_double [newXDim * newYDim * newZDim];
    ProSHADE_internal_misc::checkMemoryAllocation     ( map, __FILE__, __LINE__, __func__ );
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( newXDim * newYDim * newZDim ); iter++ ) { map[iter] = newMap[iter][0]; }
    
    //================================================ Release memory
    releaseResolutionFourierMemory                    ( origMap, fCoeffs, newFCoeffs, newMap, planForwardFourier, planBackwardRescaledFourier );
    
    //================================================ Define change in indices and return it
    corrs[0]                                          = static_cast< proshade_single > ( newXDim ) - static_cast< proshade_single > ( xDimS );
    corrs[1]                                          = static_cast< proshade_single > ( newYDim ) - static_cast< proshade_single > ( yDimS );
    corrs[2]                                          = static_cast< proshade_single > ( newZDim ) - static_cast< proshade_single > ( zDimS );
    corrs[3]                                          = static_cast< proshade_single > ( newXDim ) * static_cast< proshade_single > ( resolution / 2.0f );
    corrs[4]                                          = static_cast< proshade_single > ( newYDim ) * static_cast< proshade_single > ( resolution / 2.0f );
    corrs[5]                                          = static_cast< proshade_single > ( newZDim ) * static_cast< proshade_single > ( resolution / 2.0f );
    
    //======================================== Done
    return ;
    
}

/*! \brief This function allocates and checks the allocatio of the memory required by the Fourier resampling.
 
    This function allocates the memory required for the Fourier space re-sampling of density maps. It allocates the original map and original map coefficients arrays (these need to be of the fftw_complex type) as well as the
    re-sampled map and re-sampled map coefficients. It then also proceeds to check the memory allocation and creating the FFTW transform plans for both, the forward and re-sampled backward operations.
 
    \param[in] origMap A Reference pointer to an array where the original map data will be stored.
    \param[in] fCoeffs A Reference pointer to an array where the original map Fourier coefficients data will be stored.
    \param[in] newFCoeffs A Reference pointer to an array where the re-sampled map Fourier coefficients data will be stored.
    \param[in] newMap A Reference pointer to an array where the re-sampled map data will be stored.
    \param[in] planForwardFourier FFTW_plan which will compute the original map to Fourier coefficients transform.
    \param[in] planBackwardRescaledFourier FFTW_plan which will compute the re-sampled Fourier coefficients to re-sampled map transform.
    \param[in] xDimOld The number of indices along the x-axis of the original map.
    \param[in] yDimOld The number of indices along the y-axis of the original map.
    \param[in] zDimOld The number of indices along the z-axis of the original map.
    \param[in] xDimNew The number of indices along the x-axis of the re-sampled map.
    \param[in] yDimNew The number of indices along the y-axis of the re-sampled map.
    \param[in] zDimNew The number of indices along the z-axis of the re-sampled map.
 */
void ProSHADE_internal_mapManip::allocateResolutionFourierMemory ( fftw_complex*& origMap, fftw_complex*& fCoeffs, fftw_complex*& newFCoeffs, fftw_complex*& newMap, fftw_plan& planForwardFourier, fftw_plan& planBackwardRescaledFourier, proshade_unsign xDimOld, proshade_unsign yDimOld, proshade_unsign zDimOld, proshade_unsign xDimNew, proshade_unsign yDimNew, proshade_unsign zDimNew )
{
    //================================================ Initialise memory
    origMap                                           = new fftw_complex [xDimOld * yDimOld * zDimOld];
    fCoeffs                                           = new fftw_complex [xDimOld * yDimOld * zDimOld];
    newFCoeffs                                        = new fftw_complex [xDimNew * yDimNew * zDimNew];
    newMap                                            = new fftw_complex [xDimNew * yDimNew * zDimNew];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( origMap,          __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( fCoeffs,          __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( newFCoeffs,       __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation     ( newMap,           __FILE__, __LINE__, __func__ );
    
    //================================================ Create plans
    planForwardFourier                                = fftw_plan_dft_3d ( static_cast< int > ( xDimOld ), static_cast< int > ( yDimOld ), static_cast< int > ( zDimOld ), origMap,    fCoeffs, FFTW_FORWARD,  FFTW_ESTIMATE );
    planBackwardRescaledFourier                       = fftw_plan_dft_3d ( static_cast< int > ( xDimNew ), static_cast< int > ( yDimNew ), static_cast< int > ( zDimNew ), newFCoeffs, newMap,  FFTW_BACKWARD, FFTW_ESTIMATE );
    
    //================================================ Done
    return ;
    
}

/*! \brief This function releases the memory required by the Fourier resampling.
 
    This function simply deletes all the memory allocated by the allocateResolutionFourierMemory() function.
 
    \param[in] origMap A Reference pointer to an array where the original map data were stored.
    \param[in] fCoeffs A Reference pointer to an array where the original map Fourier coefficients data were stored.
    \param[in] newFCoeffs A Reference pointer to an array where the re-sampled map Fourier coefficients data were stored.
    \param[in] newMap A Reference pointer to an array where the re-sampled map data were stored.
    \param[in] planForwardFourier FFTW_plan which computed the original map to Fourier coefficients transform.
    \param[in] planBackwardRescaledFourier FFTW_plan which computed the re-sampled Fourier coefficients to re-sampled map transform.
 */
void ProSHADE_internal_mapManip::releaseResolutionFourierMemory ( fftw_complex*& origMap, fftw_complex*& fCoeffs, fftw_complex*& newFCoeffs, fftw_complex*& newMap, fftw_plan& planForwardFourier, fftw_plan& planBackwardRescaledFourier )
{
    //================================================ Delete the FFTW plans
    fftw_destroy_plan                                 ( planForwardFourier );
    fftw_destroy_plan                                 ( planBackwardRescaledFourier );
    
    //================================================ Delete the complex arrays
    delete[] origMap;
    delete[] fCoeffs;
    delete[] newFCoeffs;
    delete[] newMap;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function changes the order of Fourier coefficients in a 3D array between positive first (default) and negative first (mass centered at xMax/2, yMax/2, zMax/2 instead of 0,0,0)
 
    This function firstly determines the start and end of the positive and negative Fourier coefficients for all three axes (it assumes 3D array); this works for both odd and even axis sizes. To do this properly, in the case of
    odd axis sizes, the function needs to know whether we are changing the order from FFTW defaul positive first, or if we are changing the other way around - the last parameter serves the purpose of passing this information.
    Finally, the function then switches the order of the Fourier coefficients as requested using a temporary array that it allocates and deletes within itself.
 
    \param[in] fCoeffs A Reference pointer to an array where the Fourier coefficients for re-ordering are stored.
    \param[in] xDim The size of the x-axis dimension of the fCoeffs array.
    \param[in] yDim The size of the x-axis dimension of the fCoeffs array.
    \param[in] zDim The size of the x-axis dimension of the fCoeffs array.
    \param[in] negativeFirst Should the coefficients be stored negarive first (TRUE), or are we reversing already re-ordered array (FALSE)?
 */
void ProSHADE_internal_mapManip::changeFourierOrder ( fftw_complex*& fCoeffs, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim, bool negativeFirst )
{
    //================================================ Initialise local variables
    proshade_signed h = 0, k = 0, l = 0, origSizeArr = 0, newSizeArr = 0;
    proshade_signed xSeq1FreqStart, ySeq1FreqStart, zSeq1FreqStart, xSeq2FreqStart, ySeq2FreqStart, zSeq2FreqStart;
    
    //================================================ Find the positive and negative indices cot-offs
    if ( negativeFirst )
    {
        if ( ( xDim % 2 ) == 0 ) { xSeq1FreqStart = xDim / 2; xSeq2FreqStart = xDim / 2; } else { xSeq1FreqStart = (xDim / 2) + 1; xSeq2FreqStart = xDim / 2; }
        if ( ( yDim % 2 ) == 0 ) { ySeq1FreqStart = yDim / 2; ySeq2FreqStart = yDim / 2; } else { ySeq1FreqStart = (yDim / 2) + 1; ySeq2FreqStart = yDim / 2; }
        if ( ( zDim % 2 ) == 0 ) { zSeq1FreqStart = zDim / 2; zSeq2FreqStart = zDim / 2; } else { zSeq1FreqStart = (zDim / 2) + 1; zSeq2FreqStart = zDim / 2; }
    }
    else
    {
        if ( ( xDim % 2 ) == 0 ) { xSeq1FreqStart = xDim / 2; xSeq2FreqStart = xDim / 2; } else { xSeq1FreqStart = (xDim / 2); xSeq2FreqStart = xDim / 2 + 1; }
        if ( ( yDim % 2 ) == 0 ) { ySeq1FreqStart = yDim / 2; ySeq2FreqStart = yDim / 2; } else { ySeq1FreqStart = (yDim / 2); ySeq2FreqStart = yDim / 2 + 1; }
        if ( ( zDim % 2 ) == 0 ) { zSeq1FreqStart = zDim / 2; zSeq2FreqStart = zDim / 2; } else { zSeq1FreqStart = (zDim / 2); zSeq2FreqStart = zDim / 2 + 1; }
    }
        
    //================================================ Allocate helper array memory
    fftw_complex *hlpFCoeffs                          = new fftw_complex [xDim * yDim * zDim];
    ProSHADE_internal_misc::checkMemoryAllocation     ( hlpFCoeffs, __FILE__, __LINE__, __func__ );
    
    //================================================ Change the coefficients order
    for ( proshade_signed xIt = 0; xIt < xDim; xIt++ )
    {
        //============================================ Find x frequency
        if ( xIt < xSeq1FreqStart ) { h = xIt + xSeq2FreqStart; } else { h = xIt - xSeq1FreqStart; }
        for ( proshade_signed yIt = 0; yIt < yDim; yIt++ )
        {
            //======================================== Find y frequency
            if ( yIt < ySeq1FreqStart ) { k = yIt + ySeq2FreqStart; } else { k = yIt - ySeq1FreqStart; }
            
            for ( proshade_signed zIt = 0; zIt < zDim; zIt++ )
            {
                //==================================== Find z frequency
                if ( zIt < zSeq1FreqStart ) { l = zIt + zSeq2FreqStart; } else { l = zIt - zSeq1FreqStart; }
                
                //==================================== Find array positions
                newSizeArr                            = l   + zDim * ( k   + yDim * h   );
                origSizeArr                           = zIt + zDim * ( yIt + yDim * xIt );
                
                //==================================== Copy vals
                hlpFCoeffs[newSizeArr][0]             = fCoeffs[origSizeArr][0];
                hlpFCoeffs[newSizeArr][1]             = fCoeffs[origSizeArr][1];
            }
        }
    }
    
    //================================================ Copy the helper array to the input Fourier coefficients array
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( xDim * yDim * zDim ); iter++ ) { fCoeffs[iter][0] = hlpFCoeffs[iter][0]; fCoeffs[iter][1] = hlpFCoeffs[iter][1]; }
    
    //================================================ Release helper array memory
    delete[] hlpFCoeffs;
    
    //================================================ Done
    return ;
    
}

/*! \brief This function removes the phase from reciprocal (frequency) map.
 
    This function takes an already FFTW-ed map and its dimensions as the input and proceeds to remove the phase
    from the map. It writes over the map and does not release any memory - it is the role of the calling function
    to deal with both these features.
 
    \param[in] mapCoeffs A Reference Pointer to the frequency map, from which phase is to be removed.
    \param[in] xDim The number of indices along the x-axis of the input map.
    \param[in] yDim The number of indices along the y-axis of the input map.
    \param[in] zDim The number of indices along the z-axis of the input map.
 */
void ProSHADE_internal_mapManip::removeMapPhase ( fftw_complex*& mapCoeffs, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim )
{
    //================================================ Set local variables
    proshade_double real, imag, mag, phase;
    proshade_unsign arrayPos                          = 0;
    proshade_double normFactor                        = static_cast<proshade_double> ( xDim * yDim * zDim );
    
    //================================================ Iterate through the map
    for ( proshade_unsign uIt = 0; uIt < xDim; uIt++ )
    {
        for ( proshade_unsign vIt = 0; vIt < yDim; vIt++ )
        {
            for ( proshade_unsign wIt = 0; wIt < zDim; wIt++ )
            {
                //==================================== Var init
                arrayPos                              = wIt + zDim * ( vIt + yDim * uIt );
                real                                  = mapCoeffs[arrayPos][0];
                imag                                  = mapCoeffs[arrayPos][1];

                //==================================== Get magnitude and phase with mask parameters
                mag                                   = std::sqrt ( (real*real) + (imag*imag) );;
                phase                                 = 0.0; // This would be std::atan2 ( imag, real ); - but here we remove the phase.

                //==================================== Save the phaseless data
                mapCoeffs[arrayPos][0]                = ( mag * cos(phase) ) / normFactor;
                mapCoeffs[arrayPos][1]                = ( mag * sin(phase) ) / normFactor;
            }
        }
    }

    //================================================ Done
    return ;
    
}

/*! \brief Function for creating "fake" half-maps.
 
    This function takes the internal map and an empty map array and proceeds to find all neighbours within the kernel distance to each
    map points. It then computes the average of these neighbours and saves this average as the value of the "fake" half-map. These are
    then useul for map masking, as the correlation between the "fake" half-maps and the original maps give nice masks according to Rangana
    - see him about how well these work.
 
    \param[in] map A Reference Pointer to the map which should be blurred/sharpened.
    \param[in] blurredMap A Reference Pointer to the variable which will store the modified map.
    \param[in] xDimS The number of indices along the x axis of the map.
    \param[in] yDimS The number of indices along the y axis of the map.
    \param[in] zDimS The number of indices along the z axis of the map.
    \param[in] fakeMapKernel The amount of neighbours in any direction whose average is to be used to get the current point.
 */
void ProSHADE_internal_mapManip::getFakeHalfMap ( proshade_double*& map, proshade_double*& fakeHalfMap, proshade_unsign xDimS, proshade_unsign yDimS, proshade_unsign zDimS, proshade_signed fakeMapKernel )
{
    //================================================ Set local variables
    proshade_signed xDim                              = static_cast< proshade_signed > ( xDimS );
    proshade_signed yDim                              = static_cast< proshade_signed > ( yDimS );
    proshade_signed zDim                              = static_cast< proshade_signed > ( zDimS );
    proshade_signed currentPos, neighArrPos, neighXPos, neighYPos, neighZPos;
    proshade_double neighSum;
    proshade_double neighCount                        = pow ( ( ( fakeMapKernel * 2 ) + 1 ), 3.0 ) - 1.0;

    //================================================ Blur the coeffs
    for ( proshade_signed uIt = 0; uIt < xDim; uIt++ )
    {
        for ( proshade_signed vIt = 0; vIt < yDim; vIt++ )
        {
            for ( proshade_signed wIt = 0; wIt < zDim; wIt++ )
            {
                //==================================== Var init
                currentPos                            = wIt + zDim * ( vIt + yDim * uIt );
                neighSum                              = 0.0;
                
                //==================================== Average neighbours
                for ( proshade_signed xCh = -fakeMapKernel; xCh <= +fakeMapKernel; xCh++ )
                {
                    for ( proshade_signed yCh = -fakeMapKernel; yCh <= +fakeMapKernel; yCh++ )
                    {
                        for ( proshade_signed zCh = -fakeMapKernel; zCh <= +fakeMapKernel; zCh++ )
                        {
                            if ( ( xCh == 0 ) && ( yCh == 0 ) && ( zCh == 0 ) ) { continue; }
                            
                            //======================== Find the nieghbour peak indices (with periodicity)
                            neighXPos                 = uIt + xCh; if ( neighXPos >= xDim ) { neighXPos -= xDim; }; if ( neighXPos < 0 ) { neighXPos += xDim; }
                            neighYPos                 = vIt + yCh; if ( neighYPos >= yDim ) { neighYPos -= yDim; }; if ( neighYPos < 0 ) { neighYPos += yDim; }
                            neighZPos                 = wIt + zCh; if ( neighZPos >= zDim ) { neighZPos -= zDim; }; if ( neighZPos < 0 ) { neighZPos += zDim; }
                            neighArrPos               = neighZPos + zDim * ( neighYPos + yDim * neighXPos );
                            
                            //======================== Add to average
                            neighSum                 += map[neighArrPos];
                        }
                    }
                }
                
                //==================================== Save the average to "fake" half-map
                fakeHalfMap[currentPos]               = neighSum / neighCount;
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for creating the correlation mask.
 
    This function is currently not used and should probably be deleted. The proper implementation of this masking approach should be available in the EMDA software.
 
    \param[in] map A Reference Pointer to the map which should be blurred/sharpened.
    \param[in] blurredMap A Reference Pointer to the variable which stores the fake half-map.
    \param[in] correlationMask A Reference Pointer to empty map where the mask will be saved.
    \param[in] xDimS The number of indices along the x axis of the map.
    \param[in] yDimS The number of indices along the y axis of the map.
    \param[in] zDimS The number of indices along the z axis of the map.
    \param[in] corrMaskKernel The amount of neighbours in any direction whose correlation is to be used to get the current points correlation.
 */
void ProSHADE_internal_mapManip::getCorrelationMapMask ( proshade_double*& map, proshade_double*& fakeHalfMap, proshade_double*& correlationMask, proshade_unsign xDimS, proshade_unsign yDimS, proshade_unsign zDimS, proshade_signed corrMaskKernel )
{
    //================================================ Set local variables
    proshade_signed xDim = static_cast< proshade_signed > ( xDimS ), yDim = static_cast< proshade_signed > ( yDimS ), zDim = static_cast< proshade_signed > ( zDimS ), currentPos, neighArrPos, neighXPos, neighYPos, neighZPos, corrIter;
    proshade_unsign noCorrVals                        = static_cast<proshade_unsign> ( pow ( ( ( corrMaskKernel * 2 ) + 1 ), 3 ) );
    
    //================================================ Alocate memory
    proshade_double *origMap                          = new proshade_double [noCorrVals];
    proshade_double *fakeHM                           = new proshade_double [noCorrVals];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation ( origMap, __FILE__, __LINE__, __func__ );
    ProSHADE_internal_misc::checkMemoryAllocation ( fakeHM,  __FILE__, __LINE__, __func__ );
    
    //================================================ Blur the coeffs
    for ( proshade_signed uIt = 0; uIt < xDim; uIt++ )
    {
        for ( proshade_signed vIt = 0; vIt < yDim; vIt++ )
        {
            for ( proshade_signed wIt = 0; wIt < zDim; wIt++ )
            {
                //==================================== Var init
                currentPos                            = wIt + zDim * ( vIt + yDim * uIt );
                corrIter                              = 0;
                
                //==================================== Average neighbours
                for ( proshade_signed xCh = -corrMaskKernel; xCh <= +corrMaskKernel; xCh++ )
                {
                    for ( proshade_signed yCh = -corrMaskKernel; yCh <= +corrMaskKernel; yCh++ )
                    {
                        for ( proshade_signed zCh = -corrMaskKernel; zCh <= +corrMaskKernel; zCh++ )
                        {
                            //======================== Find the nieghbour peak indices (with periodicity)
                            neighXPos                 = uIt + xCh; if ( neighXPos >= xDim ) { neighXPos -= xDim; }; if ( neighXPos < 0 ) { neighXPos += xDim; }
                            neighYPos                 = vIt + yCh; if ( neighYPos >= yDim ) { neighYPos -= yDim; }; if ( neighYPos < 0 ) { neighYPos += yDim; }
                            neighZPos                 = wIt + zCh; if ( neighZPos >= zDim ) { neighZPos -= zDim; }; if ( neighZPos < 0 ) { neighZPos += zDim; }
                            neighArrPos               = neighZPos + zDim * ( neighYPos + yDim * neighXPos );
                            
                            //======================== Add to correct arrays
                            origMap[corrIter]         = map[neighArrPos];
                            fakeHM[corrIter]          = fakeHalfMap[neighArrPos];
                            corrIter                 += 1;
                        }
                    }
                }
                
                //==================================== Save the correlation comparison result
                correlationMask[currentPos]           = ProSHADE_internal_maths::pearsonCorrCoeff ( origMap, fakeHM, noCorrVals );
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function converts distance in Angstroms to distance in map indices.
 
    This function finds out how many indices are required to cover a space of size "dist" in Angstroms.
    If you need this rounded in any way (ceil, floor, ...), just apply appropriate function to the output
    of this function.
 
    \param[in] xDim The number of map indices along the x-axis.
    \param[in] yDim The number of map indices along the y-axis.
    \param[in] zDim The number of map indices along the z-axis.
    \param[in] xAngs The map size in Angstroms along the x-axis.
    \param[in] yAngs The map size in Angstroms along the y-axis.
    \param[in] zAngs The map size in Angstroms along the z-axis.
    \param[in] dist The distance in Angstroms to be converted
 */
proshade_single ProSHADE_internal_mapManip::getIndicesFromAngstroms ( proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_single dist )
{
    //================================================ Compute
    proshade_single ret                               = static_cast< proshade_single > ( ProSHADE_internal_mapManip::myRound ( std::max ( dist / ( xAngs / static_cast<proshade_single> ( xDim ) ),
                                                                                                                               std::max ( dist / ( yAngs / static_cast<proshade_single> ( yDim ) ),
                                                                                                                               dist / ( zAngs / static_cast<proshade_single> ( zDim ) ) ) ) ) );

    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief This function connects blobs in mask.
 
    \param[in] mask A pointer reference to mask map in which blobs should be connected.
    \param[in] xDim The number of map indices along the x-axis.
    \param[in] yDim The number of map indices along the y-axis.
    \param[in] zDim The number of map indices along the z-axis.
    \param[in] xAngs The map size in Angstroms along the x-axis.
    \param[in] yAngs The map size in Angstroms along the y-axis.
    \param[in] zAngs The map size in Angstroms along the z-axis.
    \param[in] maskThres The threshold which will be used for applying mask.
 */
void ProSHADE_internal_mapManip::connectMaskBlobs ( proshade_double*& mask, proshade_signed xDim, proshade_signed yDim, proshade_signed zDim, proshade_single xAngs, proshade_single yAngs, proshade_single zAngs, proshade_single maskThres )
{
    //================================================ Initialise variables
    proshade_double* hlpMap                           = new proshade_double[xDim * yDim * zDim];
    proshade_signed addSurroundingPoints              = static_cast< proshade_signed > ( std::max ( 3L, static_cast<proshade_signed> ( std::ceil ( getIndicesFromAngstroms( static_cast< proshade_unsign > ( xDim ), static_cast< proshade_unsign > ( yDim ), static_cast< proshade_unsign > ( zDim ), xAngs, yAngs, zAngs, static_cast< proshade_single > ( std::max( xAngs, std::max( yAngs, zAngs ) ) * 0.1f ) ) ) ) ) );
    proshade_signed currPos, neighXPos, neighYPos, neighZPos, neighArrPos;
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( hlpMap, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy the mask map
    for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( xDim * yDim * zDim ); iter++ ) { hlpMap[iter] = mask[iter]; }
    
    //================================================ Repeat as many times as needed
    for ( proshade_signed it = 0; it < addSurroundingPoints; it++ )
    {
        //============================================ For each x, y and z
        for ( proshade_signed xIt = 0; xIt < xDim; xIt++ )
        {
            for ( proshade_signed yIt = 0; yIt < yDim; yIt++ )
            {
                for ( proshade_signed zIt = 0; zIt < zDim; zIt++ )
                {
                    //================================ Current position
                    currPos                           = zIt + zDim * ( yIt + yDim * xIt );
                    
                    //================================ If zero, ignore
                    if ( hlpMap[currPos] < static_cast< proshade_double > ( maskThres ) ) { continue; }
                    
                    //================================ Check neighbours
                    for ( proshade_signed xCh = -1; xCh <= +1; xCh++ )
                    {
                        for ( proshade_signed yCh = -1; yCh <= +1; yCh++ )
                        {
                            for ( proshade_signed zCh = -1; zCh <= +1; zCh++ )
                            {
                                if ( ( xCh == 0 ) && ( yCh == 0 ) && ( zCh == 0 ) ) { continue; }
                                
                                //==================== Find the nieghbour indices (without periodicity)
                                neighXPos             = xIt + xCh; if ( neighXPos < 0 ) { continue; } if ( neighXPos >= xDim ) { continue; }
                                neighYPos             = yIt + yCh; if ( neighYPos < 0 ) { continue; } if ( neighYPos >= yDim ) { continue; }
                                neighZPos             = zIt + zCh; if ( neighZPos < 0 ) { continue; } if ( neighZPos >= zDim ) { continue; }
                                neighArrPos           = neighZPos + zDim * ( neighYPos + yDim * neighXPos );
                                
                                //==================== Add to mask if this point is below it (as it is a neighbour to a point which is part of the mask)
                                if ( hlpMap[neighArrPos] < static_cast< proshade_double > ( maskThres ) ) { mask[neighArrPos] = static_cast< proshade_double > ( maskThres ); }
                            }
                        }
                    }
                }
            }
        }
        
        //============================================ Now copy the updated mask to the temporary helper, unless last iteration
        for ( proshade_unsign iter = 0; iter < static_cast<proshade_unsign> ( xDim * yDim * zDim ); iter++ ) { hlpMap[iter] = mask[iter]; }
    }
    
    //================================================ Release memory
    delete[] hlpMap;
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for modifying boundaries to a mathematically more pleasant values.
 
    This function does two things. Firstly, it attempts to find even boundaries which have a small sum of the prime factors - i.e.
    this is good for further processing. Secondly, it attempts to make boundaries with sizes within some threshold the same. This
    is purely human thing :-).
 
    \param[in] bounds A proshade_signed pointer to array of 6 in which the current bounds are (0 = minX; 1 = maxX; 2 = minY; 3 = maxY; 4 - minZ; 5 = maxZ).
    \param[in] xDim The number of indices along the x axis of the map.
    \param[in] yDim The number of indices along the y axis of the map.
    \param[in] zDim The number of indices along the z axis of the map.
    \param[in] boundsDiffThres Number of indices by which boudaries in different dimensions can differ and still be made the same.
 */
void ProSHADE_internal_mapManip::beautifyBoundaries ( proshade_signed*& bounds, proshade_unsign xDim, proshade_unsign yDim, proshade_unsign zDim, proshade_signed boundsDiffThres )
{
    //================================================ If extra bounds space pushed the bounds over the physical map, freely add up to 10 indices over the current value
    while ( bounds[1] >= static_cast<proshade_signed> ( xDim ) ) { xDim += 10; }
    while ( bounds[3] >= static_cast<proshade_signed> ( yDim ) ) { yDim += 10; }
    while ( bounds[5] >= static_cast<proshade_signed> ( zDim ) ) { zDim += 10; }
    
    //================================================ Find if better bouds exist in terms of prime numbers
    proshade_signed addToX                            = betterClosePrimeFactors ( bounds[1] - bounds[0] + 1, static_cast< proshade_signed > ( xDim ) );
    proshade_signed addToY                            = betterClosePrimeFactors ( bounds[3] - bounds[2] + 1, static_cast< proshade_signed > ( yDim ) );
    proshade_signed addToZ                            = betterClosePrimeFactors ( bounds[5] - bounds[4] + 1, static_cast< proshade_signed > ( zDim ) );
    
    //================================================ Figure if similar sizes should not be forced to be the same
    proshade_signed XtoY                              = std::abs ( addToX - addToY );
    proshade_signed XtoZ                              = std::abs ( addToX - addToZ );
    proshade_signed YtoZ                              = std::abs ( addToY - addToZ );
    
    if ( ( ( XtoY < boundsDiffThres ) && ( XtoZ < boundsDiffThres ) ) ||
         ( ( XtoY < boundsDiffThres ) && ( YtoZ < boundsDiffThres ) ) ||
         ( ( XtoZ < boundsDiffThres ) && ( YtoZ < boundsDiffThres ) ) )
    {
        //============================================ Dealing with chain ( a <= b <= c type ) where at least two dista are smaller than threshold
        proshade_signed maxSize                       = std::max ( addToX, std::max ( addToY, addToZ ) );
        addToX                                        = maxSize;
        addToY                                        = maxSize;
        addToZ                                        = maxSize;
    }
    else
    {
        //============================================ Only two may be similar enough
        if ( XtoY <= boundsDiffThres )
        {
            proshade_signed maxSize                   = std::max ( addToX, addToY );
            addToX                                    = maxSize;
            addToY                                    = maxSize;
        }
        if ( XtoZ <= boundsDiffThres )
        {
            proshade_signed maxSize                   = std::max ( addToX, addToZ );
            addToX                                    = maxSize;
            addToZ                                    = maxSize;
        }
        if ( YtoZ <= boundsDiffThres )
        {
            proshade_signed maxSize                   = std::max ( addToY, addToZ );
            addToY                                    = maxSize;
            addToZ                                    = maxSize;
        }
    }
    
    //================================================ Add the extra space appropriately
    distributeSpaceToBoundaries                       ( bounds[0], bounds[1], ( bounds[1] - bounds[0] + 1 ), addToX );
    distributeSpaceToBoundaries                       ( bounds[2], bounds[3], ( bounds[3] - bounds[2] + 1 ), addToY );
    distributeSpaceToBoundaries                       ( bounds[4], bounds[5], ( bounds[5] - bounds[4] + 1 ), addToZ );
    
    //================================================ Done
    return ;
    
}

/*! \brief Function for finding close numbers with better prime factors.
 
    This function takes a range of integer numbers and proceeds to find if a close number in the range to the range start
    does not have better (smaller sum of) prime factors, taking into account the distance of this better number to the
    starting position. This is useful for finding boundaries with nice prime number sizes.
 
    \param[in] fromRange Which number the search will be from and possibilities compared to.
    \param[in] toRange The maximum number to which the posibilities should be searched.
    \param[out] X The value in the range which has better prime factors, or the first value if none other has.
 */
proshade_signed ProSHADE_internal_mapManip::betterClosePrimeFactors ( proshade_signed fromRange, proshade_signed toRange )
{
    //================================================ Initialise variables
    proshade_signed ret                               = fromRange;
    std::vector < proshade_signed > posibles, hlp;
    proshade_signed sum;
    
    //================================================ Find the sums of prime factors for each number in the whole range
    for ( proshade_signed iter = fromRange; iter < toRange; iter++ )
    {
        sum                                           = 0;
        hlp                                           = ProSHADE_internal_maths::primeFactorsDecomp ( iter+1 );
        for ( proshade_unsign i = 0; i < static_cast<proshade_unsign> ( hlp.size() ); i++ ) { sum += hlp.at(i); }
        hlp.clear                                     ( );
        ProSHADE_internal_misc::addToSignedVector     ( &posibles, sum );
    }
    
    //================================================ Find a better number
    for ( proshade_signed iter = fromRange; iter < toRange; iter++ )
    {
        //============================================ Ignore odd numbers
        if ( iter %2 != 0 ) { continue; }
        
        //============================================ Better number needs to have lower prime factor sum, but also needs not to be too far
        if ( posibles.at( static_cast< size_t > ( iter - fromRange ) ) < ( posibles.at( static_cast< size_t > ( ret - fromRange ) ) - ( iter - ret ) ) ) { ret = iter; }
    }
    
    //================================================ In the case of large prime number, just add one for even number
    if ( ( ret % 2 != 0 ) && ( ret < ( toRange - 1 ) ) ) { ret += 1; }
    
    //================================================ Done
    return                                            ( ret );
    
}

/*! \brief Function for adding space to boundaries within the map confines.
 
    This function takes the current boudaries and attempts to add the extra space specified by the parameters 3 and 4 equally
    to both the top and bottom boundaries, making sure the bottom boundary does not go belowe 0 and the top boundary does not
    exceed the final parameter.
 
    \param[in] minBound Reference to where the lower boundary value is held.
    \param[in] maxBound Reference to where the upper boundary value is held.
    \param[in] oldBoundRange The range between the original boundaries.
    \param[in] newBoundRange The range which should be between the new boundaries.
 */
void ProSHADE_internal_mapManip::distributeSpaceToBoundaries ( proshade_signed& minBound, proshade_signed& maxBound, proshade_signed oldBoundRange, proshade_signed newBoundRange )
{
    //================================================ Only bother when sometings is to be added
    if ( newBoundRange > oldBoundRange )
    {
        //============================================ How much are we adding?
        proshade_signed distributeThis                = newBoundRange - oldBoundRange;
        
        while ( distributeThis != 0 )
        {
            minBound                                 -= 1;
            distributeThis                           -= 1;
            
            if ( distributeThis != 0 )
            {
                maxBound                             += 1;
                distributeThis                       -= 1;
            }
        }
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function copies an old map to a new map with different boundaries.
 
    This function takes old and new structure details and the original structure map. It then proceed to iterate through the
    new structure map, checking if this particular point is inside the original structure's map and if it is, it will copy the
    value of this point from the original to the new map. Otherwise, it will set the new map's point value to 0.0. This is used
    when creating a new structure from an old structure given new bounds, which can be larger or smaller than the originals.
 
    \param[in] xFrom The starting x-axis index of the new map.
    \param[in] xTo The last x-axis index of the new map.
    \param[in] yFrom The starting y-axis index of the new map.
    \param[in] yTo The last y-axis index of the new map.
    \param[in] zFrom The starting z-axis index of the new map.
    \param[in] zTo The last z-axis index of the new map.
    \param[in] origXFrom The starting x-axis index of the original (copied) map.
    \param[in] origYFrom The starting y-axis index of the original (copied) map.
    \param[in] origZFrom The starting z-axis index of the original (copied) map.
    \param[in] yDimIndices The size of the y-axis dimension in the new map.
    \param[in] zDimIndices The size of the z-axis dimension in the new map.
    \param[in] origXDimIndices The size of the x-axis dimension in the old map.
    \param[in] origYDimIndices The size of the y-axis dimension in the old map.
    \param[in] origZDimIndices The size of the z-axis dimension in the old map.
    \param[in] newMap Pointer reference to the array where new map will be saved to.
    \param[in] origMap Pointer to the array where the original map can be accessed.
 */
void ProSHADE_internal_mapManip::copyMapByBounds ( proshade_signed xFrom, proshade_signed xTo, proshade_signed yFrom, proshade_signed yTo, proshade_signed zFrom, proshade_signed zTo, proshade_signed origXFrom, proshade_signed origYFrom, proshade_signed origZFrom, proshade_unsign yDimIndices, proshade_unsign zDimIndices, proshade_unsign origXDimIndices, proshade_unsign origYDimIndices, proshade_unsign origZDimIndices, proshade_double*& newMap, proshade_double* origMap )
{
    //================================================ Initialise variables
    proshade_signed newMapIndex, oldMapIndex, oldX, oldY, oldZ, newX, newY, newZ;
    
    //=============================================== For all values in the new map
    for ( proshade_signed xIt = xFrom; xIt <= xTo; xIt++ )
    {
        //============================================ Index position init
        newX                                          = ( xIt - xFrom );
        oldX                                          = ( newX + ( xFrom - origXFrom ) );
        
        for ( proshade_signed yIt = yFrom; yIt <= yTo; yIt++ )
        {
            //======================================== Index position init
            newY                                      = ( yIt - yFrom );
            oldY                                      = ( newY + ( yFrom - origYFrom ) );
            
            for ( proshade_signed zIt = zFrom; zIt <= zTo; zIt++ )
            {
                //==================================== Index position init
                newZ                                  = ( zIt - zFrom );
                oldZ                                  = ( newZ + ( zFrom - origZFrom ) );
                
                //==================================== Index arrays
                newMapIndex                           = newZ + static_cast< proshade_signed > ( zDimIndices ) * ( newY + static_cast< proshade_signed > ( yDimIndices ) * newX );
                oldMapIndex                           = oldZ + static_cast< proshade_signed > ( origZDimIndices ) * ( oldY + static_cast< proshade_signed > ( origYDimIndices ) * oldX );
                
                //============================ Check if we are adding new point
                if ( ( ( oldX < 0 ) || ( oldX >= static_cast< proshade_signed > ( origXDimIndices ) ) ) ||
                     ( ( oldY < 0 ) || ( oldY >= static_cast< proshade_signed > ( origYDimIndices ) ) ) ||
                     ( ( oldZ < 0 ) || ( oldZ >= static_cast< proshade_signed > ( origZDimIndices ) ) ) )
                {
                    //================================ Yes, this is a new point, no known value for it
                    newMap[newMapIndex]               = 0.0;
                    continue;
                }
                
                //==================================== If we got here, this should be within the old map, so just copy  value
                newMap[newMapIndex]                   = origMap[oldMapIndex];
            }
        }
    }
    
    //================================================ Done
    return ;
    
}
