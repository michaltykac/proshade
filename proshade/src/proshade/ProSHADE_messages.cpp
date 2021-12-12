/*! \file ProSHADE_messages.cpp
    \brief This source file contains all user message functions.
 
    The functions defined in this source file are used by ProSHADE to report various things to the user either using stdout (for progress messages) or the stderr
    (for error messages).
 
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
#include "ProSHADE_messages.hpp"

/*! \brief Wellcome message printing.
 
    This function prints the first message that appears when ProSHADE is being run.
 
    \param[in] verbose Signed int value stating how loud the user requested the ProSHADE run to be.
 */
void ProSHADE_internal_messages::printWellcomeMessage ( proshade_signed verbose )
{
    if ( verbose >= 0 )
    {
        std::cout << "ProSHADE " << PROSHADE_VERSION << ":" << std::endl << "============================" << std::endl << std::endl << std::flush;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief Final message printing.
 
    This function prints the last message that appears when ProSHADE is being run, assuming the run ends with no problems.
 
    \param[in] verbose Signed int value stating how loud the user requested the ProSHADE run to be.
 */
void ProSHADE_internal_messages::printTerminateMessage ( proshade_signed verbose )
{
    if ( verbose >= 0 )
    {
        std::cout << std::endl << "======================" << std::endl << "ProSHADE run complete." << std::endl << "Time taken: " << std::clock() / CLOCKS_PER_SEC << " seconds." << std::endl << "======================" << std::endl << std::endl << std::flush;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief General stdout message printing.
 
    This function is used to print progress messages to the stdout. It takes the verbosity level, the message inportance level and the
    message as such and decides whether the message should be printed, doing so if required.
 
    \param[in] verbose Int value stating how loud the user requested the ProSHADE run to be.
    \param[in] messageLevel Int value stating how important the message is.
    \param[in] message String of the actual message to be displayed.
    \param[in] messageShift How much should  the message be shifted (i.e. if this a sub-process of larger task, set this to 1 instead of 0).
 */
void ProSHADE_internal_messages::printProgressMessage ( proshade_signed verbose, proshade_signed messageLevel, std::string message, proshade_signed messageShift )
{
    if ( verbose >= messageLevel )
    {
        if ( messageLevel > 0 )
        {
            std::cout << " ";
        }
        
        for ( proshade_signed iter = 0; iter < ( messageLevel + messageShift ); iter++ )
        {
            std::cout << "... ";
        }
        
        std::cout << message << std::endl << std::flush;
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief General stderr message printing (used for warnings).
 
    This function is used to print warnings to the user through the stderr (cerr) stream. It will not terminate the program
    run, as it is just a warning. Note that errors (exceptions) are handled elsewhere.
 
    \param[in] verbose Int value stating how loud the user requested the ProSHADE run to be.
    \param[in] message String of the actual message to be displayed.
    \param[in] message String of the warning code to be displayed.
 */
void ProSHADE_internal_messages::printWarningMessage ( proshade_signed verbose, std::string message, std::string warnCode )
{
    if ( verbose >= -2 )
    {
        std::cerr << std::endl << message << std::endl << std::flush;
        std::cerr << " ... CODE: " << warnCode << std::endl << std::endl << std::flush;
        
    }
    
    //================================================ Done
    return ;
    
}

/*! \brief This function prints the help screen in the case -h is called, or if command line arguments cannot be parsed.
 
 */
void ProSHADE_internal_messages::printHelp [[noreturn]] ( void )
{
    //================================================ Print the help screen
    std::cout << "ProSHADE " << PROSHADE_VERSION << ":" << std::endl << "==========================" << std::endl << std::endl << std::flush;
    std::cout << "                                                                                " << std::endl << std::flush;
    std::cout << "DESCRIPTION:                                                                    " << std::endl;
    std::cout << "   ProSHADE is a library of functionalities for computing distances between     " << std::endl;
    std::cout << " two three-dimensional macromolecular structures, finding symmetry in a par-    " << std::endl;
    std::cout << " ticular structure and general structure manipulation including detection of    " << std::endl;
    std::cout << " optimal rotation and translation for structure overlay. The undelying method   " << std::endl;
    std::cout << " is spherical harmonics decomposition of theoretical or experimental density    " << std::endl;
    std::cout << " maps.                                                                          " << std::endl;
    std::cout << "   ProSHADE is available as an executable, C++ library and Python module.       " << std::endl;
    std::cout << " The following dialogue describes the usage of the executable, for help with    " << std::endl;
    std::cout << " the library or module, see the examples and the documentation available in     " << std::endl;
    std::cout << " the installation (for CCP-EM in the checkout) directory.                       " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "AUTHORS:                                                                        " << std::endl;
    std::cout << "   Michal Tykac                                         <tykacm@ibt.cas.cz>     " << std::endl;
    std::cout << "   Garib N. Murshudov                                                           " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "MODES:                                                                          " << std::endl;
    std::cout << "    There are several different modes in which ProSHADE can be run, depen-      " << std::endl;
    std::cout << " ding on the required functionality. The selection between these modes is       " << std::endl;
    std::cout << " done using the following options.                                              " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -D or --distances                                                           " << std::endl;
    std::cout << "            The shape distances will be computed between the first supplied     " << std::endl;
    std::cout << "            structure and all other structures. Requires at least two           " << std::endl;
    std::cout << "            structures.                                                         " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -M or --mapManip                                                            " << std::endl;
    std::cout << "            The input maps will be re-boxed using the internal map masking and  " << std::endl;
    std::cout << "            boundary finding procedures. Requires at least one structure. In    " << std::endl;
    std::cout << "            case of co-ordinate input, these will be converted to map using     " << std::endl;
    std::cout << "            Gemmi library.                                                      " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -S or --symmetry                                                            " << std::endl;
    std::cout << "            Detect if any C, D, T or I symmetries are present in all supplied   " << std::endl;
    std::cout << "            structures.                                                         " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -O or --strOverlay                                                          " << std::endl;
    std::cout << "            Given two structures, find the optimal overlay using the            " << std::endl;
    std::cout << "            rotation and translation functions. The first structure is          " << std::endl;
    std::cout << "            always unchanged, while a rotated and translated version of the     " << std::endl;
    std::cout << "            second structure will be written to the \'--overlayFile\' option      " << std::endl;
    std::cout << "            path or its default value \'./movedStructure\' file.                  " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "ARGUMENTS:                                                                      " << std::endl;
    std::cout << "    The following options can be used to to supply information and values       " << std::endl;
    std::cout << " to be used when executing the functionality - i.e. they all require some       " << std::endl;
    std::cout << " argument to follow. Some of the settings may be mandatory for some, but not    " << std::endl;
    std::cout << " other modes. See the MODES section above.                                      " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -! or --verbose                                 [DEFAULT:            1]     " << std::endl;
    std::cout << "            The verbosity of the run. Accepted values are from 0 to 4 with      " << std::endl;
    std::cout << "            increasing amount of lines being printed.                           " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -f or --file                                    [DEFAULT:         NONE]     " << std::endl;
    std::cout << "            File name (including path) of the input coordinate or map file.     " << std::endl;
    std::cout << "            For multiple files, use the option multiple times.                  " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -u                                              [DEFAULT:         TRUE]     " << std::endl;
    std::cout << "            Switch the value of this boolean. If true, any input PDB files      " << std::endl;
    std::cout << "            will be forced to have P1 spacegroup, CRYST1 value otherwise.       " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -w                                              [DEFAULT:         TRUE]     " << std::endl;
    std::cout << "            Switch the value of this boolean. If true, all water molecules      " << std::endl;
    std::cout << "            in input PDB files will be removed.                                 " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -x                                              [DEFAULT:         TRUE]     " << std::endl;
    std::cout << "            Switch the value of this boolean. If true, only the first PDB       " << std::endl;
    std::cout << "            file model will be used, all models will be used otherwise.         " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -r or --resolution                              [DEFAULT:         NONE]     " << std::endl;
    std::cout << "            The resolution to which the calculations are to be done and to      " << std::endl;
    std::cout << "            which PDB files theoretical maps will be sampled to.                " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -b or --bandwidth                               [DEFAULT:         AUTO]     " << std::endl;
    std::cout << "            The bandwidth to which spherical harmonics decomposition shoud      " << std::endl;
    std::cout << "            be computed to. For automatic determination supply 0 or nothing.    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -s or --sphereDists                             [DEFAULT:         AUTO]     " << std::endl;
    std::cout << "            The distance in Angstroms between any two concentric spheres to     " << std::endl;
    std::cout << "            which the internal map representation will be mapped to. Use        " << std::endl;
    std::cout << "            0.0 for automatic determination.                                    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -e or --extraSpace                              [DEFAULT:         10.0]     " << std::endl;
    std::cout << "            The supplied number of Angstroms will be added to the structure     " << std::endl;
    std::cout << "            internal map representation in order to avoid clashes with          " << std::endl;
    std::cout << "            adjacent cells.                                                     " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -H or --coordExtraSpace                         [DEFAULT:         10.0]     " << std::endl;
    std::cout << "            The extra space added to any co-ordinates before computing their    " << std::endl;
    std::cout << "            theoretical density map to make sure no atom is at the boundary.    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -i or --integOrder                              [DEFAULT:         AUTO]     " << std::endl;
    std::cout << "            The order to which the Gauss-Legendre integration should be         " << std::endl;
    std::cout << "            computed to. For automatic determination use 0.                     " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -t or --taylorCap                               [DEFAULT:           10]     " << std::endl;
    std::cout << "            The cap on Taylor series calculation (used for the Legendre         " << std::endl;
    std::cout << "            polynomial computation). No automatic determination available       " << std::endl;
    std::cout << "            yet.                                                                " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -d or --pdbTempFact                             [DEFAULT:         -1.0]     " << std::endl;
    std::cout << "            Some PDB files have issues with B-factors (like all being 0.0)      " << std::endl;
    std::cout << "            and to allow simple dealing with this, if this value is >= 0.0,     " << std::endl;
    std::cout << "            then all B-factors of all PDB inputs will be set to this value.     " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -F or --keepNegDens                             [DEFAULT:         TRUE]     " << std::endl;
    std::cout << "            Some input files have negative density that causes differences      " << std::endl;
    std::cout << "            to be detected between structures that appear identical. By         " << std::endl;
    std::cout << "            default negative density is removed, this option keeps it in.       " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --maskFile                                      [DEFAULT: \"./maskFile\"]     " << std::endl;
    std::cout << "            The filename to which the mask will be saved to. The extension      " << std::endl;
    std::cout << "            will be added as well as the structure index (order of input).      " << std::endl;
    std::cout << "    -G or --applyMask                                 [DEFAULT:         \"\"]     " << std::endl;
    std::cout << "            This option allows supplying a map mask, which will be applied      " << std::endl;
    std::cout << "            before any other processing in the input map. This option is only   " << std::endl;
    std::cout << "            available for the map manipulation and symmetry detection tasks.    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -z or --fourierWeights                                    [DEFAULT: \"\"]     " << std::endl;
    std::cout << "            This option allows the user to supply a map file which has values   " << std::endl;
    std::cout << "            that will be applied to the Fourier transform of the input map. It  " << std::endl;
    std::cout << "            is expected to be in the same format as standard coefficients order " << std::endl;
    std::cout << "            of FFTW (version 3) - i.e. f_0, ..., f_N/2, f_-N/2+1, ..., f_-1.    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --maskBlurring                                  [DEFAULT:        350.0]     " << std::endl;
    std::cout << "            The B-factor (temperature factor) increase, which should be         " << std::endl;
    std::cout << "            applied to blurr the map for its subsequent masking.                " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --maskThreshold                                 [DEFAULT:          3.0]     " << std::endl;
    std::cout << "            The number of inter-quartile ranges from median which will be       " << std::endl;
    std::cout << "            used to determine the map masking threshold.                        " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --boundsSpace                                   [DEFAULT:          3.0]     " << std::endl;
    std::cout << "            The number of angstroms to be added to all re-boxing determined     " << std::endl;
    std::cout << "            bounds to make sure no important surface information is lost.       " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --boundsThreshold                               [DEFAULT:            0]     " << std::endl;
    std::cout << "            The number of indices which can be added to a dimension in order    " << std::endl;
    std::cout << "            to make two dimension sizes the same.                               " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -g or --reBoxedFilename                         [DEFAULT:    \"reBoxed\"]     " << std::endl;
    std::cout << "            The file name to which the re-boxed structure will be saved to.     " << std::endl;
    std::cout << "            The extension will be added as well as the structure index (order   " << std::endl;
    std::cout << "            of input).                                                          " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --EnLWeight                                     [DEFAULT:          1.0]     " << std::endl;
    std::cout << "            The exponential weight to be applied to the shell distance for      " << std::endl;
    std::cout << "            the energy levels descriptor.                                       " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --peakNeigh                                     [DEFAULT:            1]     " << std::endl;
    std::cout << "            Number of points in each dimension that need to be lower for        " << std::endl;
    std::cout << "            peak to be detected.                                                " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --peakThres                                     [DEFAULT:         AUTO]     " << std::endl;
    std::cout << "            Number of IQRs from median for small peaks threshold for remo-      " << std::endl;
    std::cout << "            ving small peaks.                                                   " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --missAxThres                                   [DEFAULT:          0.3]     " << std::endl;
    std::cout << "            The fraction of axes that can be missing for missing axes           " << std::endl;
    std::cout << "            search to be initiated.                                             " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --sameAxComp                                     [DEFAULT:        0.01]     " << std::endl;
    std::cout << "            The difference in dot product of two vectors for them to be         " << std::endl;
    std::cout << "            still considered to be the same.                                    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --axisComBeh or -q                               [DEFAULT:        TRUE]     " << std::endl;
    std::cout << "            Should the maximum difference in dot product of two vectors for     " << std::endl;
    std::cout << "            them to be still considered to be the same decrease with fold of    " << std::endl;
    std::cout << "            tested symmetry?                                                    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --bicubSearch or -A                              [DEFAULT:        TRUE]     " << std::endl;
    std::cout << "            Should the bi-cubic interpolation for sphere peaks be used to       " << std::endl;
    std::cout << "            improve the axis by searching between grid indices?                 " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --maxSymPrime or -B                              [DEFAULT:          30]     " << std::endl;
    std::cout << "            The automated symmetry search starts by looking for prime number    " << std::endl;
    std::cout << "            folds and then for multiples of any folds found. This sets the      " << std::endl;
    std::cout << "            maximum prime number to use in the search.                          " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --minPeakHeight or -o                           [DEFAULT:          0.3]     " << std::endl;
    std::cout << "            The minimum average peak height for symmetry axis to be still       " << std::endl;
    std::cout << "            considered as \"real\" for the symmetry detection.                    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --fscThres or -C                               [DEFAULT:          0.33]     " << std::endl;
    std::cout << "            The Fourier Shell Correlation value the axes need to achieve in     " << std::endl;
    std::cout << "            order to be considered \"real\" by the symmetry detection algorithm.  " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --peakMinThres or -E                           [DEFAULT:          0.75]     " << std::endl;
    std::cout << "            The average peak height value the axes need to achieve in order     " << std::endl;
    std::cout << "            to be considered \"possible\" by the symmetry detection algorithm.    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --reqSym                                        [DEFAULT:           \"\"]     " << std::endl;
    std::cout << "            This is where the user states any particular symmetry he is         " << std::endl;
    std::cout << "            interested in detecting. The way to specify any symmetry is to      " << std::endl;
    std::cout << "            first use the letter for the symmetry type (C = cyclic, D = di-     " << std::endl;
    std::cout << "            hedral, T = tetrahedral, O = octahedral or I = icosahedral) and     " << std::endl;
    std::cout << "            for C and D symmetries to follow them with the requested fold.      " << std::endl;
    std::cout << "            I.e. \"C4\" means cyclic symmetry with fold 4.                        " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --overlayFile                          [DEFAULT:      \"movedStructure\"]     " << std::endl;
    std::cout << "            Filename to where the translated and rotated moving structure       " << std::endl;
    std::cout << "            with optimal placement relative to the static structure will be     " << std::endl;
    std::cout << "            saved to. Extension will be added automatically                     " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --overlayJSONFile or -y      [DEFAULT: \"movedStructureOperations.json\"]     " << std::endl;
    std::cout << "            Filename to where the translations and rotation operations required " << std::endl;
    std::cout << "            for moving the \"moving\" structure to overlay the \"static\" structure " << std::endl;
    std::cout << "            will be written into.                                               " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "FLAGS:                                                                          " << std::endl;
    std::cout << "    The following options can be used to override the default values and        " << std::endl;
    std::cout << " specify the execution path.                                                    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --invertMap                                     [DEFAULT:        FALSE]     " << std::endl;
    std::cout << "            In the case of getting the wrong hand of the structure during its   " << std::endl;
    std::cout << "            creation, ProSHADE can switch these by inverting the map (i.e.      " << std::endl;
    std::cout << "            any x,y,z = -x,-y,-z )                                              " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --normalise                                     [DEFAULT:        FALSE]     " << std::endl;
    std::cout << "            Should the internal map and any written out maps be normalised      " << std::endl;
    std::cout << "            to mean 0.0 and standard deviation 1.0?                             " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --mask                                          [DEFAULT:        FALSE]     " << std::endl;
    std::cout << "            Should the internal map be masked by blurring by \"maskBlurring\"     " << std::endl;
    std::cout << "            and thresholding by \"maskThreshold\" IQRs? NOTE: This option will    " << std::endl;
    std::cout << "            not save the mask. For saving the mask, use \"saveMask\".             " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --saveMask                                      [DEFAULT:        FALSE]     " << std::endl;
    std::cout << "            This option will do map masking as well as save the used mask       " << std::endl;
    std::cout << "            to the filename supplied in \"maskFile\".                             " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    -R or --mapReboxing                             [DEFAULT:        FALSE]     " << std::endl;
    std::cout << "            Should the map be re-boxed? If yes, please note that the masking    " << std::endl;
    std::cout << "            will be turned on automatically.                                    " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --sameBoundaries                                [DEFAULT:        FALSE]     " << std::endl;
    std::cout << "            Use same boundaries for multiple maps? Useful for half-maps.        " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --center or -c                                  [DEFAULT:        FALSE]     " << std::endl;
    std::cout << "            Should the map be moved to centre of mass (COM) before process-     " << std::endl;
    std::cout << "            ing? This will not affect map overlay positioning.                  " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --changeMapResol or -j                          [DEFAULT:         TRUE]     " << std::endl;
    std::cout << "            Should the map sampling be changed to fit the required resolu-      " << std::endl;
    std::cout << "            tion, or should it be left alone? Uses Fourier space re-sampling.   " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --changeMapTriLin or -a                         [DEFAULT:        FALSE]     " << std::endl;
    std::cout << "            Should the map sampling be changed to fit the required resolu-      " << std::endl;
    std::cout << "            tion, or should it be left alone? Uses tri-linear interpolation re- " << std::endl;
    std::cout << "            sampling.                                                           " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --noPhase or -p                                 [DEFAULT:         TRUE]     " << std::endl;
    std::cout << "            This option forces removal of phase from the internal map rep-      " << std::endl;
    std::cout << "            resentation.                                                        " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --progressive or -k                             [DEFAULT:        FALSE]     " << std::endl;
    std::cout << "            In this case, only the optimal number of bands will be used for     " << std::endl;
    std::cout << "            each sphere, depending on its size. This is in cotrast to the       " << std::endl;
    std::cout << "            default all spheres same number of bands setting.                   " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --noEnL or -l                                   [DEFAULT:         TRUE]     " << std::endl;
    std::cout << "            Is the computation of the energy levels descriptor required?        " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --noTrS or -m                                   [DEFAULT:         TRUE]     " << std::endl;
    std::cout << "            Is the computation of the trace sigma descriptor required?          " << std::endl;
    std::cout << "                                                                                " << std::endl;
    std::cout << "    --noFRF or -n                                   [DEFAULT:         TRUE]     " << std::endl;
    std::cout << "            Is the computation of the full rotation function descriptor         " << std::endl;
    std::cout << "            required?                                                           " << std::endl;
    std::cout << "                                                    [DEFAUlT:        FALSE]     " << std::endl;
    std::cout << "    -I or --symCentre                                                           " << std::endl;
    std::cout << "            Should symmetry centre be sought using phaseless map symmetry       " << std::endl;
    std::cout << "            before symmetry detection task is done?                             " << std::endl;
    std::cout << std::endl << std::flush;
    
    //================================================ Done
    exit                                              ( EXIT_SUCCESS );
    
}
