######################################################
######################################################
#   \file getMaskedSymmetryPredictions_centre.py
#   \brief This file shows how ProSHADE can be used in conjunction with EMDA to improve the symmetry predictions.
#
#   This code allows a simple symmetry detection run on a long list of EMDB structures (assuming these are downloaded
#   and available locally). It requires and input file with the following format:
#
#   EMDB_CODE[space]USED_SYMMETRY[space]STRUCTURE_RESOLUTION[space]PDB_MODEL_CODE
#
#   e.g.
#
#   EMD-0001 C1 3.4 6gh5
#
#   The code will parse this file and use all the settings supplied in the "Global settings" section to iteratively
#   run ProSHADE symmetry detection on each structure. If there is a mask file available, it will use it, otherwise,
#   the code will use EMDA to compute mask from map and use this mask instead.
#
#   The results are then saved in two files, the condensed one containing the determined symmetry as well as the used
#   symmetry for each EMDB entry, while the full file contains all the detected C axes for each entry.
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
#   \version   0.7.6.2
#   \date      DEC 2021
######################################################
######################################################


######################################################
### Import modules
### ==============
###
### This is where Python modules are loaded.
###

### System modules
import os
import subprocess
import numpy
import proshade
import time
import mrcfile
import shutil
import time
import sys

### Import ProSHADE
import proshade

### Import mrcfile
import mrcfile

### Import EMDA
import emda.emda_methods as emda_methods


######################################################
### Global settings
### ===============
###
### This is where all the settings are given.
###

resolution                                            = 8.0
minimalAllowedResolution                              = 20.0
mapReSampling                                         = True
symmetryCentering                                     = True
comCentering                                          = False
verbosity                                             = -1
inputFileName                                         = "../symmetryPredictions/emdb_spa_210329.dat"
outputFileName                                        = "results_allKnownEMDB_SYM_resol-"
EMDBDataPath                                          = "/Users/mysak/BioCEV/proshade/xx_EMDBSymmetry"
unreleasedIDsList                                     = [ "EMD-10163", "EMD-10165", "EMD-10166", "EMD-10168", "EMD-10169", "EMD-10170", "EMD-10174", "EMD-21320", "EMD-4320", "EMD-4522", "EMD-4523", "EMD-4524", "EMD-4606", "EMD-4607", "EMD-4718", "EMD-5039", "EMD-6758", "EMD-8144", "EMD-8145" ]
tooLargeIDsList                                       = [ "EMD-0174", "EMD-11111", "EMD-20091", "EMD-21648", "EMD-0880", "EMD-11008", "EMD-0436", "EMD-11040", "EMD-0618" ]

######################################################
### Local settings
### ==============
###
### This is where all the internal settings are set,
### no user manipulation is required.
###

startFrom                                             = 0
resolutionFilename                                    = resolution
outResCondensed                                       = 0
outResAxes                                            = 0


######################################################
### Local functions
### ===============
###
### These are function written for simpler coding, but
### which should be local to this code.
###

"""
This function opens the output files, if this is done the first time, then
for writing, otherwise for appending.
"""
def openOutputFiles ( cntr ):
    ### Declare local variables
    outResCondensed                                   = ""
    outResAxes                                        = ""

    ### Open files for output
    if cntr == 0:
        outResCondensed                               = open ( outputFileName + str ( resolutionFilename ) + "_condensed.txt", "w")
        outResAxes                                    = open ( outputFileName + str ( resolutionFilename ) + ".txt", "w")
    else:
        outResCondensed                               = open ( outputFileName + str ( resolutionFilename ) + "_condensed.txt", "a")
        outResAxes                                    = open ( outputFileName + str ( resolutionFilename ) + ".txt", "a")
        
    return                                            ( outResCondensed, outResAxes )

"""
This function closes the output files. This is done after each symmetry prediction
to make sure results are written out after each step.
"""
def closeOutputFiles ( outResCondensed, outResAxes ):
    ## Close output files
    outResCondensed.close                             ( )
    outResAxes.close                                  ( )

"""
This function reads in and parses the input file into an array of arrays of values.
"""
def readInEMDBList ( filename ):
    ### Create local variables
    retList                                           = []

    ### Read in the EMDB list from the file
    emdbIDList                                        = open ( filename, 'r' )
    IDs                                               = emdbIDList.readlines ( )
    emdbIDList.close                                  ( )
    
    ### Remove C1 and None symmetries
    symIDs                                            = []
    for ln in IDs:
        id_hlp                                        = ln.split(' ')[1]
        if id_hlp != "C1" and id_hlp != "None":
            symIDs.append                             ( ln )

    ### Parse the input lines
    for idLine in symIDs:
    
        ### Get the ID, depositted symmetry and resolution (we do not need the PDB file ID)
        id                                            = str ( idLine.split(' ')[0] )
        declaredSym                                   = str ( idLine.split(' ')[1] )
        declaredRes                                   = str ( idLine.split(' ')[2] )
        cleanID                                       = str ( id.split ('EMD-')[1] )
        
        ### Add to return
        retList.append                                ( [ cleanID, id, declaredSym, declaredRes ] )

    ### Return the results
    return                                            ( retList )
    
"""
This function checks structure ID for being on the unreleased (and thus unusable) list.
"""
def skipUnreleasedIDs ( fullID ):
    ### Check for list belonging
    if ( fullID in unreleasedIDsList ):
        return                                        ( True )
    else:
        return                                        ( False )
        
"""
This function checks structure ID for being on the too large (and thus unusable) list.
"""
def skipTooLargeIDs ( fullID ):
    ### Check for list belonging
    if ( fullID in tooLargeIDsList ):
        return                                        ( True )
    else:
        return                                        ( False )

"""
This function returns the map file name from the id and the data path
"""
def findMapFile ( id, dataPath ):
    ### Find the map file
    mapPath                                           = os.path.join ( os.path.join( dataPath, "EMD-" + str( id ) ), "emd_" + str( id ) + ".map.gz" )
    if not os.path.isfile ( mapPath ):
        sys.exit                                      ( "Failed to find the map " + str( mapPath ) )
    else:
        return                                        ( mapPath )
        
"""
This function returns the mask file name from the id and the data path, if it exists
"""
def findMaskFile ( id, dataPath ):
    ### Find the mask file
    maskPath                                          = os.path.join ( os.path.join( dataPath, "EMD-" + str( id ) ), "emd_" + str( id ) + "_mask.map.gz" )
    if not os.path.isfile ( maskPath ):
        ### Report progress
        print                                         ( " ... " + "Mask not found." )
        return                                        ( "" )
    else:
        ### Report progress
        print                                         ( " ... " + "Mask found." )
        return                                        ( maskPath )

"""
This function returns the resolution to which the ProSHADE computation is to be done. It
takes into account the minimal required resolution, the required resolution and the map
resolution to make the decision. If this structure is not to be considered, it returns -1.
"""
def checkMapResolution ( mapFile, declRes ):
    ### Set local variables
    locRes                                            = resolution

    ### Read in map using mrc file
    mrc                                               = mrcfile.open ( mapFile, mode = "r+" )

    ### Find dims in Angstroms
    uc                                                = mrc.header.cella[ [ "x", "y", "z" ] ]
    cell                                              = uc.view( ( "f4", 3 ) )
    
    ### Find dims in indices
    dimSize                                           = mrc.data.shape

    ### Compute map resolution
    mapResol                                          = numpy.max ( ( cell / dimSize ) * 2.0 )
    
    ### If map resolution below decent resolution, report
    if mapResol > minimalAllowedResolution:
        return                                        ( -1.0 )
        
    ### Check map resolution against the minimal required resolution
    if ( locRes < float ( mapResol ) ):
        locRes                                        = float ( mapResol )

    ### Check against reported resolution
    if ( locRes < float ( declaredRes ) ):
        locRes                                        = float ( declaredRes )
        
    ### Done
    return                                            ( locRes, dimSize[0] * dimSize[1] * dimSize[2] )

"""
This function reads in the density map data from the file and proceeds to use
EMDA to compute the density mask, saving it into file "mapmask.mrc".
"""
def maskMapUsingEMDA ( mapFile, locRes ):
    ### Read in the data using mrcfile
    mrc                                               = mrcfile.open ( mapFile, mode = "r+" )
    
    ### Parse out data from the mrcfile object
    cell                                              = numpy.zeros ( 6, dtype="float" )
    uc                                                = numpy.array(mrc.header.cella)
    uc                                                = mrc.header.cella[ [ "x", "y", "z" ] ]
    cell[:3]                                          = uc.view( ( "f4", 3 ) )
    cell[3:]                                          = float ( 90.0 )
    origin                                            = [mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart]
    arr                                               = mrc.data
    
    ### Close the file to stop changes from being written onto the disc
    mrc.close                                         ( )
    
    ### Compute mask
    mapmask                                           = emda_methods.mask_from_map ( uc=cell, arr=arr, kern=4, resol=locRes, filter='butterworth', prob=0.99, itr=3, orig=origin )
    
"""
This function runs ProSHADE symmetry detection and returns the recommented symmetry type,
fold and a list of all detected C axes.
"""
def runProSHADESymmetry ( mapFile, maskFile, resol, chngSampl, cntrMap, symCenMap, verbosity ):
    ### Run proshade symmetry detection on the map
    pSet                                              = proshade.ProSHADE_settings ( )

    ### Set up the run
    pSet.task                                         = proshade.Symmetry
    pSet.setResolution                                ( resol )
    pSet.setMapResolutionChange                       ( chngSampl )
    pSet.setMapCentering                              ( cntrMap )
    pSet.verbose                                      = verbosity
    pSet.setAppliedMaskFilename                       ( maskFile )
    pSet.setSymmetryCentreSearch                      ( symCenMap )

    ### Print major settings
    print ( " ... Running: res = " + str( resol ) + " mask = " + str( maskFile ) )

    ### Read in the structure
    pStruct                                           = proshade.ProSHADE_data ( )
    pStruct.readInStructure                           ( mapFile, 0, pSet )

    ### Should we search for symmetry centre?
    if pSet.findSymCentre:
        ### Start centre detection - create the settings objects for the phaseless and phased centre detection runs
        rotCenSettingsPhased                          = proshade.ProSHADE_settings ( pSet )
        rotCenSettingsUnphased                        = proshade.ProSHADE_settings ( pSet )
    
        ### Enforce the necessary settings
        rotCenSettingsPhased.messageShift             = 1;
        rotCenSettingsPhased.moveToCOM                = False;
        rotCenSettingsUnphased.usePhase               = False;
        rotCenSettingsUnphased.requestedSymmetryType  = "onlyC";
        rotCenSettingsUnphased.moveToCOM              = False;
        rotCenSettingsUnphased.addExtraSpace          = pSet.addExtraSpace * 5.0;
            
        ###  Read in the structure and find all symmetries without using phase information
        symStr                                        = proshade.ProSHADE_data ( )
        symStr.readInStructure                        ( mapFile, 0, rotCenSettingsUnphased )
        symStr.processInternalMap                     ( rotCenSettingsUnphased )
        symStr.mapToSpheres                           ( rotCenSettingsUnphased )
        symStr.computeSphericalHarmonics              ( rotCenSettingsUnphased )
        symStr.computeRotationFunction                ( rotCenSettingsUnphased )
        symStr.detectSymmetryInStructure              ( rotCenSettingsUnphased )
    
        ### Find reliable symmetries in the Patterson map
        relSym                                        = proshade.findReliableUnphasedSymmetries ( rotCenSettingsUnphased, rotCenSettingsUnphased.verbose, rotCenSettingsUnphased.messageShift, rotCenSettingsUnphased.axisErrTolerance );
        
        ### Are there any reasonable symmetries?
        if len ( relSym ) != 0:
        
            ### If there is a D, then optimise it
            if len ( relSym ) == 2:
                proshade.optimiseDGroupAngleFromAxesHeights ( relSym, symStr, rotCenSettingsUnphased )
                
            ### Get all axes
            allCAxes                                  = symStr.getAllCSyms ( rotCenSettingsUnphased )
                
            ### Read in with phases
            del symStr
            symStr                                    = proshade.ProSHADE_data ( )
            symStr.readInStructure                    ( mapFile, 0, rotCenSettingsPhased )
            symStr.processInternalMap                 ( rotCenSettingsPhased )
            
            ### If single axis, determine point closest to COM
            if len( relSym ) == 1:
            
                ### Find the line and point on it closest to COM
                point1                                = proshade.findPointFromTranslations ( rotCenSettingsPhased, symStr, allCAxes, relSym[0] )
                axis1                                 = numpy.array( [ allCAxes[relSym[0]][1], allCAxes[relSym[0]][2], allCAxes[relSym[0]][3]] )
                COM                                   = proshade.findMAPCOMValues ( symStr )
                
                xBoxCentre                            = ( ( symStr.xTo - symStr.xFrom ) / 2 ) + symStr.xFrom
                yBoxCentre                            = ( ( symStr.yTo - symStr.yFrom ) / 2 ) + symStr.yFrom
                zBoxCentre                            = ( ( symStr.zTo - symStr.zFrom ) / 2 ) + symStr.zFrom
                
                COMFromBoxCen                         = numpy.zeros ( 3 )
                COMFromBoxCen[0]                      = xBoxCentre - ( COM[0] / ( symStr.xDimSize / symStr.xDimIndices ) )
                COMFromBoxCen[1]                      = yBoxCentre - ( COM[1] / ( symStr.yDimSize / symStr.yDimIndices ) )
                COMFromBoxCen[2]                      = zBoxCentre - ( COM[2] / ( symStr.zDimSize / symStr.zDimIndices ) )
                
                alpha1                                = numpy.dot ( point1 - COMFromBoxCen, axis1 ) / numpy.dot ( axis1, axis1 )
                
                cpVec                                 = numpy.zeros ( 3 )
                cpVec[0]                              = point1[0] + ( alpha1 * axis1[0] )
                cpVec[1]                              = point1[1] + ( alpha1 * axis1[1] )
                cpVec[2]                              = point1[2] + ( alpha1 * axis1[2] )
                
                pSet.setSymmetryCentrePosition        ( cpVec )
            
            ### If dihedral, find exact point
            if len( relSym ) == 2:
                
                ### Find the point
                point1                                = proshade.findPointFromTranslations ( rotCenSettingsPhased, symStr, allCAxes, relSym[0] )
                point2                                = proshade.findPointFromTranslations ( rotCenSettingsPhased, symStr, allCAxes, relSym[1] )
                
                axis1                                 = numpy.array( [ allCAxes[relSym[0]][1], allCAxes[relSym[0]][2], allCAxes[relSym[0]][3]] )
                axis2                                 = numpy.array( [ allCAxes[relSym[1]][1], allCAxes[relSym[1]][2], allCAxes[relSym[1]][3]] )
                
                tangentToAxes                         = numpy.cross ( axis1, axis2 )
                correctedSecondAxis                   = numpy.cross ( axis1, tangentToAxes )
                correctedFirstAxis                    = numpy.cross ( axis2, tangentToAxes )
                
                alpha1                                = numpy.dot ( point2 - point1, correctedFirstAxis  ) / numpy.dot ( axis1, correctedFirstAxis  )
                alpha2                                = numpy.dot ( point1 - point2, correctedSecondAxis ) / numpy.dot ( axis2, correctedSecondAxis )
                
                cpVec                                 = numpy.zeros ( 3 )
                cpVec[0]                              = ( ( point1[0] + ( alpha1 * axis1[0] ) ) + ( point2[0] + ( alpha2 * axis2[0] ) ) ) / 2.0
                cpVec[1]                              = ( ( point1[1] + ( alpha1 * axis1[1] ) ) + ( point2[1] + ( alpha2 * axis2[1] ) ) ) / 2.0
                cpVec[2]                              = ( ( point1[2] + ( alpha1 * axis1[2] ) ) + ( point2[2] + ( alpha2 * axis2[2] ) ) ) / 2.0
                
                pSet.setSymmetryCentrePosition        ( cpVec )
            
        ### If no symmetries, just be done
        else:
            print ( "No symmetry found in the Patterson map. Will try detecting symmetry over the centre of the box now..." )
            del symStr

    ### Do all the computations
    pStruct.processInternalMap                        ( pSet )
    pStruct.mapToSpheres                              ( pSet )
    pStruct.computeSphericalHarmonics                 ( pSet )
    pStruct.computeRotationFunction                   ( pSet )
    pStruct.detectSymmetryInStructure                 ( pSet )

    ### Retrieve results
    recSymmetryType                                   = pStruct.getRecommendedSymmetryType ( pSet )
    recSymmetryFold                                   = pStruct.getRecommendedSymmetryFold ( pSet )
    allCAxes                                          = pStruct.getAllCSyms ( pSet )
    
    ### Convert results for return
    retList                                           = []
    retList.append                                    ( recSymmetryType )
    retList.append                                    ( recSymmetryFold )
    retList.append                                    ( allCAxes )
    
    ### Release memory
    del pStruct
    del pSet
    
    ### Done
    return                                            ( retList )
    

######################################################
### Run the symmetry detection
### ==========================
###
### This code actually does the computations :-).
###

### Get the list of input structures
symIDs                                                = readInEMDBList ( inputFileName )

### For each entry on the list
counter                                               = 1
for entry in symIDs:

    ### Clean variables
    locRes                                            = resolution
    mapPath                                           = ""
    maskPath                                          = ""

    ### Rename for reading ease
    id                                                = entry[0]
    fullID                                            = entry[1]
    declaredSym                                       = entry[2]
    declaredRes                                       = entry[3]
    
    ### Ignore if so required
    if startFrom > counter:
        counter                                       = counter + 1
        continue

    ### Use only entries with decent resolution
    if float ( declaredRes ) > minimalAllowedResolution:
        counter                                       = counter + 1
        continue

    ### Ignore not yet released entries
    if skipUnreleasedIDs ( fullID ):
        counter                                       = counter + 1
        continue

    ### If map too large, ignore for now
    if skipTooLargeIDs ( fullID ):
        counter                                       = counter + 1
        continue

    ### Report progress
    print                                             ( str ( id ) + " ( " + str ( counter ) + " out of " + str ( len ( symIDs ) ) + " ):" )

    ### Figure map filename
    mapPath                                           = findMapFile  ( id, EMDBDataPath )

    ### Decide computation resolution
    ( compResolution, mapVolumeInds )                 = checkMapResolution ( mapPath, declaredRes )
    if compResolution == -1.0:
        counter                                       = counter + 1
        continue
        
    ### Report progress
    print                                             ( " ... " + "The computation resolution will be " + str( compResolution ) )

    ### Figure mask filename and if none, use EMDA
    maskPath                                          = findMaskFile ( id, EMDBDataPath )
    if ( maskPath == "" ):
        maskMapUsingEMDA                              ( mapPath, compResolution )
        maskPath                                      = "mapmask.mrc"

    ### Start timer
    startTime                                         = time.time ( )

    ### Run symmetry detection
    symRes                                            = runProSHADESymmetry ( mapPath, maskPath, compResolution, mapReSampling, comCentering, symmetryCentering, verbosity )

    ### Stop timer
    stopTime                                          = time.time ( )

    ### Report progress
    print                                             ( " ... Symmetry detection complete ( time taken: " + str( stopTime - startTime ) + " )." )

    ### Open files for output
    ( outResCondensed, outResAxes )                   = openOutputFiles ( counter )

    ### Write results
    if ( symRes[0] == "I" ) or ( symRes[0] == "O" ) or ( symRes[0] == "T" ):
        outResCondensed.write                         ( str( id ) + "\t" + str( declaredSym ) + "\t" + str( symRes[0]  ) + "\t" + str( stopTime - startTime ) + "\n" )
    else:
        outResCondensed.write                         ( str( id ) + "\t" + str( declaredSym ) + "\t" + str( symRes[0]  ) + str( symRes[1] ) + "\t" + str( stopTime - startTime ) + "\n" )
    outResAxes.write                                  ( str( id ) + " :\n==========\n" )
    for ax in range ( 0, len ( symRes[2] ) ):
        outResAxes.write                              ( str ( symRes[2][ax][0] ) + "\t" + str ( symRes[2][ax][1] ) + "\t" + str( symRes[2][ax][2] ) + "\t" + str( symRes[2][ax][3] ) + "\t" + str( symRes[2][ax][4] ) + "\t" + str( symRes[2][ax][5] ) + "\t" + str( symRes[2][ax][6] ) + "\t" + str( mapVolumeInds ) + "\t" + str( compResolution ) + "\n" )
    outResAxes.write                                  ( "\n" )
    
    ### Close output files
    closeOutputFiles                                  ( outResCondensed, outResAxes )
    
    ### Move counter
    counter                                           = counter + 1

    ### End of symmetry detection for this structure

### Done

