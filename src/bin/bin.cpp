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
    \date      JUL 2020
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
 * \section download Obtaining ProSHADE
 *
 * The most recent stable version of ProSHADE is available from the \e master branch of the GitHub repository https://github.com/michaltykac/proshade, from where it can be cloned using the git application or downloaded manually using the interface. More advanced users may be interested
 * in obtaining the \e development or the \e experimental branches, which are available from the same link. The \e experimental branch is where I do all new development and it may or may not be currently compilable and working properly, while the \e development branch should always compile, but is more
 * likely to contain bugs and issues as it is the code before proper testing.
 *
 * \section index Index
 *
 * 1) \ref intro
 *
 * 2) \ref download
 *
 * 3) \ref index
 *
 * 4) \ref install
 *
 * 4.1) \ref stdSys
 *
 * 4.2) \ref installBehaviour
 *
 * 4.3) \ref otherDependencies
 *
 * 4.4) \ref installcode
 *
 * 4.5) \ref uninstall
 *
 * 5) \ref proshadeBinary
 *
 * 5.1) \ref symDetection
 *
 * 5.2) \ref distDetection
 *
 * 5.3) \ref reboxingUsage
 *
 * 5.4) \ref overlayExample
 *
 * 6) \ref libuse
 *
 * 6.1) \ref liblink
 *
 * 6.2) \ref libexamples
 *
 * 7) \ref pyusage
 *
 * 7.1) \ref pyinstall
 *
 * 7.2) \ref pyexamples
 *
 * \section install Installation
 *
 * The installation of the ProSHADE software should be done using the CMake system and the supplied CMakeLists.txt file. The minimual requiered version of CMake is 2.6, however, python modules and single source file
 * compilation will not be available unless CMake version 3.1 or higher is used. The CMakeLists.txt file assumes the standard system dependencies are installed in the system folders; for a full list of standard system dependencies,
 * please see the section \ref stdSys.
 *
 * Once all of the standard system dependencies are installed CMake can be run to create the make files. There are several options that can be used to modify the default behaviour of the installation; these typically drive the
 * installation locations and dependencies paths in the case of non-standard dependency location. Please see the section \ref installBehaviour for details as to how to use these options and what do they do.
 *
 * Please note that while the ProSHADE code is C++98 standard compatible, some of the dependencies do require at least partial support for the C++11 standard.
 *
 * \subsection stdSys Standard System Dependencies
 *
 * Generally, the following list of standard system libraries and utilities are required for successfull installation of ProSHADE on Linux systems. On MacOS systems, most of these should be installed by default except where specifically stated:
 * - \b gcc
 * - \b g++ (on Ubuntu and Debian) or \b gcc-c++ (on CentOS and SuSe)
 * - \b gfortran (on Ubuntu and Debian) or \b gcc-gfortran (on CentOS and SuSe )
 * - \b make
 * - \b cmake
 * - \b m4
 * - \b fftw3-dev (on Ubuntu, Debian and SuSe) or \b fftw3-devel (on CentOS)
 * - \b libblas-dev (on Ubuntu and Debian) or \b blas-devel (on CentOS and SuSe)
 * - \b liblapack-dev (on Ubuntu and Debian) or \b lapack-devel (on CentOS) or \b lapack-dev (on SuSe)
 * - \b python (on Ubuntu, Debian and SuSe) or \b python2 (on CentOS)
 * - \b python-pip (on Ubuntu, Debian and SuSe) or \b python2-pip (on CentOS)
 * - \b python-dev (on Ubuntu, Debian and SuSe) or \b python2-devel (on CentOS)
 * - \b python3
 * - \b python3-pip
 * - \b python3-dev (on Ubuntu, Debian and SuSe) or \b python3-devel
 * - \b swig
 * - \b git
 * - \b numpy (installed using pip or pip3 separately for python2.x and python3.x)
 *
 * \subsection installBehaviour CMake options
 *
 * \b -DINSTALL_LOCALLY=ON or \b OFF
 * - This option is used to decide whether all the installed ProSHADE components are installed in the local source directory (value \b ON ) or whether they are instead installed in the system folders (value \b OFF ). This option
 * applies to the binary, the C++ library, the python2 and python3 modules (which are installed in the appropriate site-packages folder if the option is \b OFF ) and the headers as well.
 *
 * \b -DINSTALL_BIN_DIR=/path
 * - This option is used to manually specify the folder to which the ProSHADE binary shold be installed into.
 *
 * \b -DINSTALL_LIB_DIR=/path
 * - This option is used to manually specify the folder to which the ProSHADE C++ library should be installed into.
 *
 * \b -DINSTALL_INC_DIR=/path
 * - This option is used to specify the folder to which the ProSHADE header files required by the library should be installed into.
 *
 * \b -DCUSTOM_FFTW3_LIB_PATH=/path
 * - This option is used to supply the path to the libfftw3.a/so/dylib in the case where ProSHADE CMake installation fails to detect the FFTW3 dependency. This is typically the case when FFTW3 is installed outside of the
 * standard FFTW3 installation locations.
 *
 * \b -CUSTOM_FFTW3_INC_PATH=/path
 * - This option is used to supply the path to the fftw3.h file in the case where ProSHADE CMake installation fails to detect the FFTW3 dependency. This is typically the case when FFTW3 is installed outside of the
 * standard FFTW3 installation locations.
 *
 * \b -DCUSTOM_LAPACK_LIB_PATH=/path
 * - This option is used to supply the path to the liblapack.a/so/dylib in the case where ProSHADE CMake installation fails to detect the LAPACK dependency. This is typically the case when the LAPACK is installed
 * outside of the standard LAPACK installation locations.
 *
 * \subsection otherDependencies Other dependencies
 *
 * ProSHADE also depends on the \e Gemmi, \e libccp4, \e MMDB2 and \e SOFT2.0 libraries. Since the installation of these libraries is non-trivial and does require some user input, these libraries are
 * supplied with the ProSHADE code and will be installed locally by the ProSHADE CMake installation. Please note that these dependencies do have their own licences (the CCP4 licence, the GPL licence, ...) and therefore
 * this may limit the ProSHADE usage for some users beyond the ProSHADE copyright and licence itself.
 *
 * \subsection installcode Install
 *
 * In order to install ProSHADE, first please check that all the \ref stdSys are installed, preferably using a package management system such as \e apt or \e yum. Next, please navigate to any folder to which you would like to write the install
 * files; some find it useful to create a \c build folder in the ProSHADE folder in order to keep the install files in the same location as the source codes. Then, issue the following set of commands, setting the \c \path\to\ProSHADE
 * to the correct path on your system and adding any required \ref installBehaviour to the first command. Please note that \c sudo may be required for the \c make \c install command if you are installing into the system folders.
 *
 * \c cmake \c \path\to\ProSHADE
 *
 * \c make
 *
 * \c make \c install
 *
 * \subsection uninstall Uninstall
 *
 * To remove the installed ProSHADE components, the command \c make \c remove needs to be issued to the makefile originally created by the CMake call. Please note that \c sudo may need to be used if the installation was
 * done into the system folders and your current user does not have admin rights.
 *
 * \section proshadeBinary Using the ProSHADE binary
 *
 * The ProSHADE tool was developed in a modular fashion and as the usage slightly changes depending on the functionality that is required. Nonetheless, care has been taken to
 * make sure that identical or closely related features are controlled by the same command line arguments in all cases. Moreover, the GNU command line options standard
 * have been adhered to (through the getOpts library) and therefore the users familiar with other command line tools should find the entering of command line arguments
 * simple. The following subsections relate to examples of using different functionalities; for a full list of command line options, please use the \c --help command line option of the
 * ProSHADE binary.
 *
 * \subsection symDetection Symmetry Detection
 *
 * In order to detect symmetry in either a coordinate input file, or in a map input file, the ProSHADE executable needs to be supplied with the option \p -S or \p --symmetry and it will also
 * require a single input file to be supplied using the \p -f option. These two options are the only mandatory options, although there are many optional values that the user can supply to supersede
 * the default values and therefore modify the operation fo the ProSHADE executable to fit their purpose.
 *
 * One particular option regarding the symmetry detection mode should be noted; the \p --sym (or \p -u) option allows the user to state which symmetry they believe to exist in the structure. The allowed
 * values for this command line argument are "Cx", "Dx", "T", "O" and "I", where the \e x should be an integer number specifying the fold of the requested symmetry. When this option is used, it removes the
 * default behaviour of returning the highest detected symmetry and instead the symmetry requested by the user is returned, if it can be found in the structure.
 *
 * Another noteworthy option is the \c --center or \c -c option, which  tells ProSHADE to center the internal map representation over the centre of co-ordinates before running any processing of the
 * map. This may be important as ProSHADE detects symmetries over the centre of the co-ordinates and therefore a non-centered map (map which does not have the centre of mass at the centre of co-ordinates) will
 * be found to have no symmetries even if these are present, just not over the co-ordinate centre.
 *
 * To demonstrate how the tool can be run and the standard output for the symmetry mode of operation, the current version of the ProSHADE executable was used to detect the
 * symmetry of a density map of the bacteriophage T4 portal protein with the PDB accession code 3JA7 (EMDB accession code 6324), which has the \a C12 symmetry. The visualisation of the structure is
 * shown in the following figure, while the output of the ProSHADE tool follows:
 *
 * \image html ProSHADE_3JA7.jpg width=500cm
 *
*\code{.sh}
 $: ./proshade -S -f ./emd_6324.map -r 12 -a -k
 ProSHADE 0.7.3 (JUL 2020):
 ==========================

  ... Starting to read the structure: ./emd_6324.map
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Map centering not requested.
  ... Adding extra 10 angstroms.
  ... Phase information retained in the data.
  ... Starting sphere mapping procedure.
  ... Preparing spherical harmonics environment.
  ... Starting spherical harmonics decomposition.
  ... Starting self-rotation function computation.
  ... Starting C symmetry detection.
  ... Starting D symmetry detection.
  ... Starting I symmetry detection.
  ... Starting O symmetry detection.
  ... Starting T symmetry detection.
 Detected C symmetry with fold 12 .
  ...   Fold       X           Y          Z           Angle        Height
  ...    12   -0.014426   +0.016779   +0.9994       +0.5236      +0.34321

 ======================
 ProSHADE run complete.
 Time taken: 3 seconds.
 ======================
*\endcode
 *
 * \subsection distDetection Shape similarity distances
 *
 * The distances computation mode is signalled to the ProSHADE executable by the command line argument \p -D or \p --distances. This mode requires two or more structures to be
 * supplied using the \p -f command line option. At least two structures are mandatory for the ProSHADE tool to proceed. Moreover, the resolution of the structures to which the comparison should be done
 * needs to be supplied using the \p -r option. This resolution does not need to be the real resolution to which the structure(s) were solved, but rather reflects the amount of details which should be taken into
 * accout when comparing shapes. Therefore, higher resolution comparison will focus more on details of the shapes, while lower resolution comparison will focus more on the overall shape ignoring the minor
 * details. Please note that the results are calculated only for the first structure against all the remaining structures, \b not for all against all distance matrix.
 *
 * There are a number of useful options for the shape distances computation, please consult the \p --help dialogue for their complete listing.
 *
 * To demonstrate the output of the ProSHADE software tool for computing distances between structure shapes, the distances between the BALBES protein domains 1BFO_A_dom_1 and
 * 1H8N_A_dom_1 (which have similar shape) and the 3IGU_A_dom_1 domain which has a different shape, as can be seen from the following figure - the first two domains are
 * both in cluster a), while the last domain is from the cluster b). The output of the ProSHADE software tool is then shown below:
 *
 * \image html ProSHADE_dists.png width=500cm
 *
 *\code{.sh}
  $: ./proshade -D -f ./1BFO_A_dom_1.pdb -f ./1H8N_A_dom_1.pdb -f ./3IGU_A_dom_1.pdb -r 6 -k
  ProSHADE 0.7.3 (JUL 2020):
  ==========================

   ... Starting to read the structure: ./1BFO_A_dom_1.pdb
   ... Map inversion (mirror image) not requested.
   ... Map normalisation not requested.
   ... Masking not requested.
   ... Map centering not requested.
   ... Adding extra 10 angstroms.
   ... Phase information retained in the data.
   ... Starting sphere mapping procedure.
   ... Preparing spherical harmonics environment.
   ... Starting spherical harmonics decomposition.
   ... Starting to read the structure: ./1H8N_A_dom_1.pdb
   ... Map inversion (mirror image) not requested.
   ... Map normalisation not requested.
   ... Masking not requested.
   ... Map centering not requested.
   ... Adding extra 10 angstroms.
   ... Phase information retained in the data.
   ... Starting sphere mapping procedure.
   ... Preparing spherical harmonics environment.
   ... Starting spherical harmonics decomposition.
   ... Starting energy levels distance computation.
   ... Starting trace sigma distance computation.
   ... Starting rotation function distance computation.
  Distances between ./1BFO_A_dom_1.pdb and ./1H8N_A_dom_1.pdb
  Energy levels distance    : 0.790525
  Trace sigma distance      : 0.954362
  Rotation function distance: 0.793044
   ... Starting to read the structure: ./3IGU_A_dom_1.pdb
   ... Map inversion (mirror image) not requested.
   ... Map normalisation not requested.
   ... Masking not requested.
   ... Map centering not requested.
   ... Adding extra 10 angstroms.
   ... Phase information retained in the data.
   ... Starting sphere mapping procedure.
   ... Preparing spherical harmonics environment.
   ... Starting spherical harmonics decomposition.
   ... Starting energy levels distance computation.
   ... Starting trace sigma distance computation.
   ... Starting rotation function distance computation.
  Distances between ./1BFO_A_dom_1.pdb and ./3IGU_A_dom_1.pdb
  Energy levels distance    : 0.603907
  Trace sigma distance      : 0.745459
  Rotation function distance: 0.572001

  ======================
  ProSHADE run complete.
  Time taken: 3 seconds.
  ======================
 *\endcode
 *
 * \subsection reboxingUsage Re-boxing structures
 *
 * Another useful feature of the ProSHADE tool is the re-boxing of macromolecular density maps. This mode is signalled to the ProSHADE tool by the command line option \p -M or \p --mapManip followed by the \p -R option to specify that the
 * required map manipulations include re-boxing. Furthermore, a single map structure file needs to be supplied after the \p -f flag. In this mode, ProSHADE will attempt to find a suitable map mask by blurring the map (increasing the overall B-factors).
 * Consequently, it will use the map boundaries to create a new, hopefully smaller, box to which the appropriate part of the map will be copied.
 *
 * This ProSHADE functionality can be combinaed with other map manipulations, which include the map invertion (signalled by the \p --invertMap option and useful for cases where map reconstruction software mistakes the hands of the structure),
 * the map normalisation (signalled by the \p --normalise option, which makes sure the map mean is 0 and standard deviation is 1), centering of centre of mass to the centre of co-ordinates (using the \p --center or \p -c option) or the phase
 * removal (creating Patterson maps using the \p --noPhase or \p -p options).
 *
 * The location and filename of where this new map should be saved can be specified using the \p --reBoxedFilename (or \p -g ) command line option followed by the filename.
 *
 * The following snippet shows the output of the ProSHADE tool when used to re-box the TubZ-Bt four-stranded filament structure (EMDB accession code 5762 and PDB accession code 3J4S), where the original volume can be decreased to 46.9% of
 * the original structure volume and thus any linear processing of such structure will be more than twice faster and the original. The original TubZ-Bt four-stranded filament structure box is shown in the following figure as semi-transparent grey, while the
 * new box is shown in non-transparent yellow.
 *
 * \image html ProSHADE_rebox.png width=500cm
 *
 *\code{.sh}
 $ ./proshade -MR -f ./emd_5762.map
 ProSHADE 0.7.3 (JUL 2020):
 ==========================

  ... Starting to read the structure: ./emd_5762.map
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Computing mask.
  ... Map centering not requested.
  ... Adding extra 10 angstroms.
  ... Phase information retained in the data.
  ... Finding new boundaries.
  ... Creating new structure according to the new  bounds.
  ... Saving the re-boxed map into reBoxed_0.map

 ======================
 ProSHADE run complete.
 Time taken: 7 seconds.
 ======================
 \endcode
 *
 * \subsection overlayExample Optimal rotation and translation
 *
 * In order to find the rotation and translation which optimally overlays (or fits) one structure into another, be them PDB files or maps (and any combination thereof), the ProSHADE tool can be used in the Overlay
 * mode. This is signalled to the ProSHADE tool binary by the command line option \p --strOverlay or the \p -O and this mode requires exactly two structure files to be supplied using the \p -f command line options. The order of the
 * two files does matter, as the second file will always be moved to match the first structure, which will remain static.
 *
 * Due to the requirement for the second stucture movement and rotation, it is worth noting that the structure may need to be re-sampled and/or moved to the same viewing position as the first structure. This is
 * done so that only the internal representation is modified, but never the input file. However, when the overlay structure is outputted (a non-default name can be specified by the \p --overlayFile command line option) the header of this
 * output file may differ from the second structure header. Furthermore, if there is no extra space around the structure, movement and rotation may move pieces of the structure through the box boundaries to the
 * other side of the box. To avoid this, please use the \p --extraSpace option to add some extra space around the structure.
 *
 * As an example of the Overlay mode, we will be matching a single PDB structure (1BFO_A_dom_1 from the BALBES database, original structure code 1BFO) shown in part a) of the following figure to another PDB structure, this time the
 * 1H8N_A_dom_1 structure from the BALBES database, shown in part b) of the following figure. Please note that ProSHADE can fit any allowed input (map or co-ordinates) to any allowed input, it is just this example which uses two PDB files.
 * Part c) of the figure then shows the match obtained by the internal map representation of the moving structure (1H8N_A_dom_1) after rotation and translation with the static structure (1BFO_A_dom_1). Finally, part d) then shows the original
 * static structure (1BFO_A_dom_1) in brown and the rotated and translated version of the moving structure (1H8N_A_dom_1) in blue. Please note that the optimal rotation matrix and translation vector are written into the output when verbosity
 * (\p --verbose) is increased to at least 3, but are better accessed programatically (see the following sections) if you are interested in using these further.
 *
 * \image html ProSHADE_overlay.jpg width=500cm
 *
 *\code{.sh}
 $ ./proshade -O -f ./1BFO_A_dom_1.pdb -f ./1H8N_A_dom_1.pdb -r 4 -k
 ProSHADE 0.7.3 (JUL 2020):
 ==========================

  ... Starting to read the structure: ./1BFO_A_dom_1.pdb
  ... Starting to read the structure: ./1H8N_A_dom_1.pdb
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Map centering not requested.
  ... Adding extra 10 angstroms.
  ... Centering map onto its COM.
  ... Phase information removed from the data.
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Map centering not requested.
  ... Adding extra 10 angstroms.
  ... Centering map onto its COM.
  ... Phase information removed from the data.
  ... Starting sphere mapping procedure.
  ... Preparing spherical harmonics environment.
  ... Starting sphere mapping procedure.
  ... Preparing spherical harmonics environment.
  ... Starting spherical harmonics decomposition.
  ... Starting spherical harmonics decomposition.
  ... Starting rotation function computation.
  ... Starting to read the structure: ./1BFO_A_dom_1.pdb
  ... Starting to read the structure: ./1H8N_A_dom_1.pdb
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Map centering not requested.
  ... Adding extra 10 angstroms.
  ... Phase information retained in the data.
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Map centering not requested.
  ... Adding extra 10 angstroms.
  ... Phase information retained in the data.
  ... Starting sphere mapping procedure.
  ... Preparing spherical harmonics environment.
  ... Starting spherical harmonics decomposition.
  ... Starting translation function computation.

 ======================
 ProSHADE run complete.
 Time taken: 8 seconds.
 ======================
 \endcode
 *
 * \section libuse Using the ProSHADE library
 *
 * ProSHADE allows more programmatic access to its functionality through a C++ dynamic library, which is compiled at the same time as the binary is made. This library can be linked to any C++ project to allow direct access to the ProSHADE objects, functions and results. This section discusses how the ProSHADE
 * library can be linked against and how the basic objects can be accessed.
 *
 * \subsection liblink Linking against the ProSHADE library
 *
 * The ProSHADE library can be linked as any other C++ library, that is by using the \p -lproshade option when calling the compiler (tested on \e clang and \e g++ ) and including the header file (\p ProSHADE.hpp ). However, as the \p ProSHADE.hpp header file includes header files from some of the dependencies, any
 * C++ project linking against the ProSHADE library will need to provide their paths to the compiler. Moreover, if the ProSHADE library was not installed in the system folders (which are by defaul in the compiler paths), any project linking against the ProSHADE library will also need to provide the path to the libproshade.a/so/dylib
 * library file and the RPATH to the same location. The following list states all the paths that may be required for a successfull compilation against the ProSHADE library:
 *
 * - \b -I/path/to/proshade/install/include This path is required for the cmaplib dependency header file to be located correctly. It may not needed if cmaplib is installed into system folders, i.e. when ProSHADE was installed with the CMake -DINSTALL_LOCALLY=FALSE option.
 * - \b -I/path/to/proshade/extern/soft-2.0/include This path is required for the SOFT2.0 dependency header file to be located correctly (it is confusingly called fftw_wrapper.h).
 * - \b -L/path/to/proshade/install/lib This is the path the where libproshade.a/so/dylib is installed. If ProSHADE was installed using the CMake -DINSTALL_LOCALLY=FALSE option, then this path may already be available to the compiler and it may not be needed.
 * - \b -Wl, \b -rpath, \b /path/to/proshade/install/lib or \b -rpath \b /path/to/proshade/install/lib This compiler option will be required if the proshade library was not installed into a system folder which is already included in the project's RPATH.
 *
 * Overall, a compilation of a C++ project linking against the ProSHADE library may look like the following code:
 *
 *\code{.sh}
 $ clang ./proshadeProject.cpp -I/path/to/proshade/install/include \
                               -I/path/to/proshade/extern/soft-2.0/include \
                               -L/path/to/proshade/install/lib \
                               -lproshade \
                               -rpath /path/to/proshade/install/lib \
                               -o ./proshadeProject
 \endcode
 *
 * or
 *
 *\code{.sh}
 $ g++ ./proshadeProject.cpp -I/path/to/proshade/install/include \
                             -I/path/to/proshade/extern/soft-2.0/include \
                             -L/path/to/proshade/install/lib \
                             -lproshade -Wl,-rpath,/path/to/proshade/install/lib \
                             -o ./proshadeProject
 \endcode
 *
 *\subsection libexamples Examples of ProSHADE library usage
 *
 * There are several examples of C++ code which makes use of the ProSHADE dynamic library to compute the standard ProSHADE functionalities and access the results programmatically (i.e. without the need for parsing any log files).
 *
 * \b Simple \b access
 *
 * The examples are avaialbe in the \b /path/to/proshade/examples/libproshade folder and are divided into two categories of four examples. The source files with names starting with \e simpleAccess_... provide a \e black \e box experience similar to using ProSHADE binary. The user firstly
 * creates a \p ProSHADE_settings object, which provides all the variables that can be set in order to drive which ProSHADE functionality is required and how it should be done. Next, the user needs to create the \p ProSHADE_run object, whose constructor takes the already created and
 * filled \p ProSHADE_setings object as its only argument. This constructor will then proceed to compute all required information according to the settings object and return when complete. While the computation is being done, the execution is with the ProSHADE library and any C++ project
 * using this mode will be waiting for the ProSHADE library to finish. Once the computation is complete, the execution will be returned to the calling C++ project and the results will be accessible through public functions of the \p ProSHADE_run object. The following code shows a very simple
 * example of how ProSHADE can be run in this mode, but for more specific examples the users should review the \e simpleAccess_... example files.
 *
 *\code{.cpp}
 #include "ProSHADE.hpp"

 int main ( int argc, char **argv )
 {
     //================================================ Create the settings object
     ProSHADE_Task task                                = Distances;
     ProSHADE_settings* settings                       = new ProSHADE_settings ( task );

     //================================================ Set the settings object up
     settings->setResolution                           ( 10.0 );
     settings->addStructure                            ( "./str1.pdb" );
     settings->addStructure                            ( "./str2.pdb" );
 
     //================================================ Run ProSHADE. This may take some time, depending on what computations are required.
     ProSHADE_run* runProshade                         = new ProSHADE_run ( settings );
 
     //================================================ Access the results
     std::vector< proshade_double > energyDistances    = runProshade->getEnergyLevelsVector     ( );
     std::vector< proshade_double > traceDistances     = runProshade->getTraceSigmaVector       ( );
     std::vector< proshade_double > rotFunDistances    = runProshade->getRotationFunctionVector ( );

     //================================================ Release the memory
     delete runProshade;
     delete settings;
 
     //================================================ Done
     return EXIT_SUCCESS;
 }
 \endcode
 *
 * \b Advanced \b access
 *
 * The second set of examples of usage of the ProSHADE library are the source files with names starting with \e advancedAccess_... . These files provide examples of how individual ProSHADE functions can be arranged to provide the results of the main ProSHADE functionalities. Using the ProSHADE tool in the manner shown in these example
 * codes gives the user more control over the execution and it also allows the user to modify the behaviour directly. On the other hand, using ProSHADE in this way required a bit more understanding than the simple \e black \e box approach and this documentation should be helpful for all who wish to use ProSHADE this way. Interested users
 * are advised to review all the \e advancedAccess_... source files as well as the following simple example code.
 *
 *\code{.cpp}
 #include "ProSHADE.hpp"

 int main ( int argc, char **argv )
 {
     //================================================ Create the settings object
     ProSHADE_Task task                                = Symmetry;
     ProSHADE_settings* settings                       = new ProSHADE_settings ( task );

     //================================================ Create the structure objects
     ProSHADE_internal_data::ProSHADE_data* simpleSym  = new ProSHADE_internal_data::ProSHADE_data ( settings );
 
     //================================================ Read in the structures
     simpleSym->readInStructure                        ( "./emd_6324.map", 0, settings );
 
     //================================================ Process internal map
     simpleSym->processInternalMap                     ( settings );

     //================================================ Map to spheres
     simpleSym->mapToSpheres                           ( settings );
 
     //================================================ Compute spherical harmonics decompostion
     simpleSym->computeSphericalHarmonics              ( settings );
 
     //================================================ Compute self-rotation function
     simpleSym->getRotationFunction                    ( settings );
 
     //================================================ Detect the recommended symmetry
     std::vector< proshade_double* > symAxes;
     simpleSym->detectSymmetryInStructure              ( settings, &symAxes );
     std::string symmetryType                          = simpleSym->getRecommendedSymmetryType ( settings );
     proshade_unsign symmetryFold                      = simpleSym->getRecommendedSymmetryFold ( settings );
 
     //================================================ Write out the symmetry detection results
     std::cout << "Detected symmetry: " << symmetryType << "-" << symmetryFold << " with axes:" << std::endl;
     for ( proshade_unsign axIt = 0; axIt < static_cast<proshade_unsign> ( symAxes.size() ); axIt++ )
     {
        std::cout << "Symmetry axis number " << axIt << std::endl;
        std::cout << " ... Fold             " << symAxes.at(axIt)[0] << std::endl;
        std::cout << " ... XYZ:             " << symAxes.at(axIt)[1] << " ; " << symAxes.at(axIt)[2] << " ; " << symAxes.at(axIt)[3] << std::endl;
        std::cout << " ... Angle (radians): " << symAxes.at(axIt)[4] << std::endl;
        std::cout << " ... Axis peak:       " << symAxes.at(axIt)[5] << std::endl;
     }
 
     //================================================ Release the memory
     delete simpleSym;
     delete settings;
 
     //================================================ Done
     return EXIT_SUCCESS;
 }
 \endcode
 *
 *\section pyusage Using the Python modules
 *
 * ProSHADE also provides python2 and python3 modules which allow the same programmatical access to the ProSHADE tool as the dynamic C++ library. These modules are produced using the SWIG tool and contain all the functionality of the dynamic library. Furthermore,
 * both the modules (the python2 and the python3 versions) support the Numpy arrays and do require it to be installed.
 *
 * \subsection pyinstall Python modules installation notes
 *
 * There are several caveats that the user should be aware of before using the python modules. This section will discuss these, but if there are any issues installing the modules, please contact the author.
 *
 * \b Automatic \b installation
 * - The python modules will be installed automatically as long as the CMake version on your system is version \b 3.1 or higher. Having this CMake version and missing any of the python specific dependencies (\e SWIG, \e python, \e python3,
 * \e Numpy, ...) will cause build errors. If you do not want python modules installed, you can modify the CMakeLists.txt script, but this is recommended only for experienced CMake users.
 *
 * \b Python \b versions
 * - The CMake installation scripts use the \e which \e python or \e which \e python2 commands to detect the specific python2 interpreter version and the \e which \e python3 command to detect the specific python3 version. The
 * python library (libpython) appropriate for python version selected in this manner will then be used and linked against by the ProSHADE modules, which will then work only for this python interpreter. Therefore, if you are using multiple python versions,
 * please make sure that the \e which command in the command line points to the required python version before installing ProSHADE.
 *
 * \b Module \b locations
 * - The ProSHADE python modules will be installed into the correct \b site-packages folder for the selected python version if the CMake option -DINSTALL_LOCALLY=FALSE is used. If not, then the modules will be installed locally in the \b install folder in
 * the ProSHADE folder cloned by git. If the latter is the case, then you will need to suppy your python interpreter with the path to the location of the module before you can import the module - this can be done by calling the \e sys.path.append function of the
 * \e sys module. Please note that if the module import fails or ProSHADE gives segmentation error upon simple creation of the first ProSHADE object, then it is likely that the module was built for different python interpreter version.
 *
 * \subsection pyexamples Python module examples
 *
 * Similarly to the ProSHADE dynamic library, the python code examples are available in the \b /path/to/proshade/examples/python2 (or \b python3 ) folders. The examples are basically identical between the python2 and python3 folders, so just review the examples for the version of
 * python that you intend to use.
 *
 * \b Simple \b access
 *
 * Similarly to the dynamic library case, there are three types of examples available for the python modules. The first set of examples (files named \e simpleAccess_... ) show the \e black \e box experience, which is similar to using ProSHADE binary.
 * The user needs to create the \b ProSHADE_settings object and can then supply it with all the settings values which will then drive the ProSHADE computations. The same settings are available in the python modules as in the ProSHADE library; the example code below shows only a
 * small selection of these (for full selection, please see the example files). Next, the user creates the \b ProSHADE_run object, constructor of which takes the settings object as its only argument and then proceeds to do all computations required by the settings in the settings object. The
 * computations are done on this one line and if they take time, the execution of the script will be halted until ProSHADE is done computing. Once completed, the results can be retrieved from the \b ProSHADE_run object using the public accessor functions; the example code below shows
 * how the symmetry functionality can be run and results retrieved.
 *
 * \code{.py}
 """ Import the system modules """
 import sys
 import numpy
 
 """ Import ProSHADE """
 sys.path.append                               ( "/path/to/proshade/install/python3" ) # only needed if ProSHADE was installed locally
 import proshade
 
 """ Create the ProSHADE_settings object """
 pSet                                          = proshade.ProSHADE_settings ()
 
 """ Set the settings for Symmetry detection """
 # Required values
 pSet.task                                     = proshade.Symmetry                     # The task ProSHADE is expected to perform
 pSet.verbose                                  = -1                                    # How verbose should the run be? Use -1 for absolute silence.
 pSet.setResolution                            ( 8.0 )                                 # The resolution to which computations are to be done.
 pSet.addStructure                             ( "./C2.pdb" )                          # The path to the structure to be processed.
 
 """ Run the Symmetry task """
 pRun                                          = proshade.ProSHADE_run ( pSet )        # Do the computations. This takes most of the time.

 """ Get the recommended symmetry """
 detectedSymType                               = proshade.getDetectedSymmetryType ( pRun ) # Retrieve the results. Takes minimum time as these are already computed
 detectedSymFold                               = proshade.getDetectedSymmetryFold ( pRun ) # Retrieve the results. Takes minimum time as these are already computed
 detectedSymAxes                               = proshade.getDetectedSymmetryAxes ( pRun ) # Retrieve the results. Takes minimum time as these are already computed

 """ Print results """
 print ( "Detected symmetry " + str( detectedSymType ) + "-" + str( detectedSymFold ) + " with axes: " )
 print ( "Fold      x         y         z       Angle     Height" )
 for iter in range ( 0, len( detectedSymAxes ) ):
      print ( "  %s    %+1.3f    %+1.3f    %+1.3f    %+1.3f    %+1.4f" % ( detectedSymAxes[iter][0], detectedSymAxes[iter][1], detectedSymAxes[iter][2], detectedSymAxes[iter][3], detectedSymAxes[iter][4], detectedSymAxes[iter][5] ) )
 
 """ Delete the C++ pointers """
 del pRun
 del pSet
 \endcode
 *
 * \b Advanced \b access
 *
 * If the user needs more control over the ProSHADE exectution, or simply wants any behaviour not simply available by variables in the settings object, then there are the \e advancedAccess_... examples. These showcase the ability to call internal ProSHADE functions and by ordering them correctly,
 * achieving the requested functionality. This usage of the python modules does required a better understanding the of the ProSHADE tool and the functions it implements. This documentation is a good starting point as to which functions are available and the ProSHADE_tasks.cpp source file shows how
 * the internal functions can be arranged to achieve the standard ProSHADE tasks. Please be aware that most of the functions do require that a pre-requisite function is run before it, but not all of these pre-requisites have their own warning or error messages. Therefore, if any code causes segmentation error
 * (which usually kills the python interpreter), it is likely that you forgot to call some pre-requisite function.
 *
 * The following code is an example of how this approach to the ProSHADE python module can be used to compute the shape-wise distances between two structures. After importing the required modules, the code creates the setings object and sets the basic settings (for a full list of settings, please see the
 * example files). It then proceeds to create the \b ProSHADE_data objects for each structure, reads in the structures from files on the hard-drive (PDB and MAP formats are supported, the mmCIF should work as well). Next, the code processes in data - this is where map centering, masking, normalisation,
 * axis inversion, etc. happens - onto an internal ProSHADE data representation. This representation can then be mapped onto a set of concentric spheres, which can in turn have their spherical harmonics decomposition computed. Once this is done, the shape distances can be computed using the three functions
 * shown.
 *
 * \code{.py}
 """ Import the system modules """
 import sys
 import numpy
 
 """ Import ProSHADE """
 sys.path.append                               ( "/path/to/proshade/install/python3" ) # only needed if ProSHADE was installed locally
 import proshade
 
 """ Create the ProSHADE_settings object """
 pSet                                          = proshade.ProSHADE_settings ()
 
 """ Set the settings for Distances detection """
 pSet.task                                     = proshade.Distances
 pSet.verbose                                  = 1
 pSet.setResolution                            ( 6.0 )
 
 """ Create the structure objects """
 pStruct1                                      = proshade.ProSHADE_data ( pSet )
 pStruct2                                      = proshade.ProSHADE_data ( pSet )

 """ Read in the structure """
 pStruct1.readInStructure                      ( "./C2.pdb", 0, pSet )
 pStruct2.readInStructure                      ( "./testMap.map", 1, pSet )

 """ Process maps into internal representations """
 pStruct1.processInternalMap                   ( pSet )
 pStruct2.processInternalMap                   ( pSet )

 """ Map internal representations to spheres """
 pStruct1.mapToSpheres                         ( pSet )
 pStruct2.mapToSpheres                         ( pSet )

 """ Compute spherical harmonics of the mapped data """
 pStruct1.computeSphericalHarmonics            ( pSet )
 pStruct2.computeSphericalHarmonics            ( pSet )
 
 """ Get the results """
 energyLevelsDescriptor                        = proshade.computeEnergyLevelsDescriptor    ( pStruct1, pStruct2, pSet )
 traceSigmaDescriptor                          = proshade.computeTraceSigmaDescriptor      ( pStruct1, pStruct2, pSet )
 fullRotationFunctionDescriptor                = proshade.computeRotationunctionDescriptor ( pStruct1, pStruct2, pSet )

 """ Print the results """
 print ( "The energy levels distance is          %+1.3f" % ( energyLevelsDescriptor ) )
 print ( "The trace sigma distance is            %+1.3f" % ( traceSigmaDescriptor ) )
 print ( "The rotation function distance is      %+1.3f" % ( fullRotationFunctionDescriptor ) )

 
 """ Delete the C++ pointers """
 del pStruct1
 del pStruct2
 del pSet
 \endcode
 *
 * One of the advantages of this mode of operating the ProSHADE python modules is that the execution only takes the time required to compute the specific computation the function requires and therefore if the user only needs some preliminary results, or can prepare the data for execution later, this
 * is all allowed by this mode.
 *
 * \b Direct \b access
 *
 * This is the most advanced mode of using the ProSHADE tool, as it allows direct access into the internal ProSHADE working. This in turn allows supplying non-standard values as well as retrieving any partial results for alternative processing by a different code; however, it also requires the deepest
 * understanding of how ProSHADE works, what data are available at which point in the execution and it may require some data format manipulations on the side of the executing code. The following tutorial as well as this documentation should be a good starting point as well as the \b directAccess.py file.
 *
 * In order to showcase this approach, we will import the required modules:
 *
 * \code{.py}
 """ Import the system modules """
 import sys
 import numpy
 
 """ Import ProSHADE """
 sys.path.append                               ( "/path/to/proshade/install/python3" ) # only needed if ProSHADE was installed locally
 import proshade
 \endcode
 *
 * \e Reading \e a  \e structure \e from \e file
 *
 * The first step of most ProSHADE workflows will be reading a structure (be it co-ordinates or map) from a file on the hard-drive. This can be done in the same manner as shown in the \b advanced \b access section of this tutorial, that is: Firstly we create the \b ProSHADE_settings object, which
 * needs to be filled with the initial data. It does not really matter which task you select at this stage, but it may affect some of the default values and therefore it is recommended to use the correct taks. Next, the \b ProSHADE_data object is created and finally the structure is read in. Please note that
 * on some systems using relative paths may not work and it may result in ProSHADE error stating that the file type cannot be recognised. If this is the case, please use the full path.
 *
 * \code{.py}
 """ Create the ProSHADE_settings object """
 pSet                                          = proshade.ProSHADE_settings ()
 
 """ Set the settings for Distances detection """
 pSet.task                                     = proshade.Distances
 pSet.verbose                                  = 1
 pSet.setResolution                            ( 6.0 )
 
 """ Create the structure objects """
 pStruct1                                      = proshade.ProSHADE_data ( pSet )
 
 """ Read in the structure """
 pStruct1.readInStructure                      ( "./C2.pdb", 0, pSet )
 
 """ Delete the structure as we will proceed with a different object created later """
 del pStruct1
 \endcode
 *
 * \e Creating \e ProSHADE_data \e object \e from \e map
 *
 * Alteratively, \b ProSHADE_data object can be created from an already existing map and some of the basic map information data. As an example,  we will create a 1D numpy.array, which will hold the density values of a map that we would like to supply to ProSHADE. Of course this array can
 * be the result of any other python module, the only requirement is that the data type is 1D numpy.array with the XYZ axis order.
 *
 * \code{.py}
  """ Set the density map description values """
  xDimIndices                                   = 100                          # Number of indices along the x-axis.
  yDimIndices                                   = 120                          # Number of indices along the y-axis.
  zDimIndices                                   = 60                           # Number of indices along the z-axis.
  xDimAngstroms                                 = xDimIndices * 1.3            # X-axis size in Angstroms.
  yDimAngstroms                                 = yDimIndices * 1.3            # Y-axis size in Angstroms.
  zDimAngstroms                                 = zDimIndices * 1.3            # Z-axis size in Angstroms.
  xFrom                                         = int ( -xDimIndices/2 )       # Starting index of the x-axis.
  yFrom                                         = int ( -yDimIndices/2 )       # Starting index of the y-axis.
  zFrom                                         = int ( -zDimIndices/2 )       # Starting index of the z-axis.
  xTo                                           = int ( (xDimIndices/2)-1 )    # Last index of the x-axis.
  yTo                                           = int ( (yDimIndices/2)-1 )    # Last index of the y-axis.
  zTo                                           = int ( (zDimIndices/2)-1 )    # Last index of the z-axis.
  ord                                           = 0                            # The order of the MAP file (the binary writting style, leave 0 unless you strongly prefer different type).
 
  """ Create an example map (this will be a ball in the middle of the map) """
  testMap = numpy.empty ( [ ( xDimIndices * yDimIndices * zDimIndices ) ] )
  for xIt in range( 0, xDimIndices ):
      for yIt in range( 0, yDimIndices ):
          for zIt in range( 0, zDimIndices ):
              """ Compute the 1D array index """
              ind                               = zIt + zDimIndices * ( yIt + yDimIndices * xIt )
 
              """ Save density value for this point """
              testMap[ind]                      = 1.0 / ( numpy.sqrt( numpy.power ( (xDimIndices/2) - xIt, 2.0 ) +
                                                                      numpy.power ( (yDimIndices/2) - yIt, 2.0 ) +
                                                                      numpy.power ( (zDimIndices/2) - zIt, 2.0 ) ) + 0.01 )
 \endcode
 *
 * with an example map created as an 1D numpy.array, it can now be supplied to a \b ProSHADE_data object, which will then be in the same state as if the data were read in from a file. This can be done with the following call:
 *
 * \code{.py}
 """ Create ProSHADE_data object from numpy.array """
 pStruct                                       = proshade.ProSHADE_data ( pSet,                # The settings object
                                                                          "python_map_test",   # The map title
                                                                          testMap,             # The numpy.array object
                                                                          xDimAngstroms,       # The size of x-axis in Angstroms
                                                                          yDimAngstroms,       # The size of y-axis in Angstroms
 *                                                                        zDimAngstroms,       # The size of z-axis in Angstroms
 *                                                                        xDimIndices,         # The number of indices along the x-axis
 *                                                                        yDimIndices,         # The number of indices along the y-axis
 *                                                                        zDimIndices,         # The number of indices along the z-axis
 *                                                                        xFrom,               # The starting x-axis index
 *                                                                        yFrom,               # The starting y-axis index
 *                                                                        zFrom,               # The starting z-axis index
 *                                                                        xTo,                 # The last x-axis index
 *                                                                        yTo,                 # The last y-axis index
 *                                                                        zTo,                 # The last z-axis index
 *                                                                        ord )                # The map file binary order
 \endcode
 *
 * There are several assumption that the \b ProSHADE_data constructor shown above makes and not all of these are currently checked with a warning or error message. Some of these are described in the \b directAccess.py file, but the most common things to consider are the following:
 * - The constructor assumes that the angles between all three axes are 90 degrees. If this is not the case, it is the users responsibility to transform and re-sample the map before submitting it to ProSHADE (support for non-orthogonal maps is on the TODO list).
 * - The map dimensions needs to be the same as the x/y/zDimIndices variables.
 * - The following equality should hold: x/y/zTo - x/y/zFrom + 1 = x/y/zDimIndices.
 * - The constructor assumes that axis grids are equal to indices and that the origins are equal to the starting indices of each axis.
 * - The axis order is XYZ
 *
 * If some of these assumptions do not hold, the \b ProSHADE_data object is likely to still be created, but it is the users responsibility to change the \b pStruct (ProSHADE_data) object internal variables to reflect the reality or face the consequences.
 *
 * \e 3D \e arrays
 *
 *  It is possible that instead of 1D arrays as shown above, some other python module would work with maps using 3D arrays. This poses an issue for ProSHADE, as it is pythonised using SWIG and Numpy.i typedefs, which apparently do not support passing of 3D arrays between the two languages. To bridge this issue,
 *  ProSHADE provides a function called \e convert3Dto1DArray(), which can be used to convert standard 3D numpy.arrays to ProSHADE compatible 1D numpy.arrays. However, as python seems to have problems with multiple for loops, this function may be time consuming compared to other much more complex
 *  functions. Nonetheless, the following code will create a different \b ProSHADE_data object from a 3D numpy.array using the conversion function. It will also delete this object, as the example here will continue with the \b pStruct object created before.
 *
 * \code{.py}
 testMap3D = numpy.empty ( ( xDimIndices, yDimIndices, zDimIndices ) )
 for xIt in range( 0, xDimIndices ):
     for yIt in range( 0, yDimIndices ):
         for zIt in range( 0, zDimIndices ):
             testMap3D[xIt][yIt][zIt] = 1.0 / ( numpy.sqrt( numpy.power ( (xDimIndices/2) - xIt, 2.0 ) +
                                                            numpy.power ( (yDimIndices/2) - yIt, 2.0 ) +
                                                            numpy.power ( (zDimIndices/2) - zIt, 2.0 ) ) + 0.01 )

 pStruct2                                      = proshade.ProSHADE_data ( pSet,
                                                                          "python_map_test",
                                                                          proshade.convert3Dto1DArray ( testMap3D ),
                                                                          xDimAngstroms,
                                                                          yDimAngstroms,
                                                                          zDimAngstroms,
                                                                          xDimIndices,
                                                                          yDimIndices,
                                                                          zDimIndices,
                                                                          xFrom,
                                                                          yFrom,
                                                                          zFrom,
                                                                          xTo,
                                                                          yTo,
                                                                          zTo,
                                                                          ord )

 del pStruct2
 \endcode
 *
 * \e Writing \e out \e maps
 *
 * Now, it is possible to write out the internal map representation at any point since the \b ProSHADE_data object has been filled in. To demonstrate this functionality, it is immediately possible to issue the following command, which will write an \b initialMap.map file
 * onto the hard-drive.
 *
 * \code{.py}
 """ Write out the initial map """
 pStruct.writeMap                              ( "initialMap.map" )
 \endcode
 *
 * \e Getting \e back \e the \e ProSHADE \e internal \e representation \e map
 *
 * Here, we will demonstrate how the user can access the ProSHADE internal representation map from the \b ProSHADE_data object. Please note that this is not limitted to the initial map, this will work for any \b ProSHADE_data object which has map data in any stage of
 * ProSHADE computations. The user has a choice between a 1D and 3D numpy array maps being returned by ProSHADE; the indexing of the 1D map is the same as above, that is [ z + pStruct.zDimIndices * ( y + pStruct.yDimIndices * x ) ]. Please note that there is the issue
 * with 3D maps and therefore getting a 3D map may be slower (approximately 0.5 seconds per average sized map) as compared to getting a 1D map. The following code shows how the maps can be retrieved back to python:
 *
 * \code{.py}
 """ Get back the ProSHADE internal map representation """
 initialMapArray1D                             = proshade.getMapPython1D ( pStruct )
 initialMapArray3D                             = proshade.getMapPython3D ( pStruct )
 \endcode
 *
 * \e Changing \e the \e already \e supplied \e map
 *
 * A situation may also occur, where the already loaded map needs to be processed, for example, if you want to read in a map file, but then you want to remove some part of the map using algorithm not available in ProSHADE or in the case where you want to apply some ProSHADE processing, then do some more
 * processing with algorithms not available in ProSHADE and then return to ProSHADE for some more calculations. To do this, ProSHADE allows you to modify the retrieved (or any other) map in any way; in the following example we will remove the last 3 indices along the y-axis dimension and divide (or multiply) all
 * density values by 2. Next, it is now the users responsibility to modify the ProSHADE object to reflect the new size and position of the map and finally, the new map can be pushed into ProSHADE using the \e setNewMapPythonXD()) functions (where X can be 1 or 3). The following example shows how this can be done for both, the 1D and the 3D maps.
 *
 * \code{.py}
 """ Create a new map arrays with the y-dimension decreased by 3 indices """
 newMapArr1D                                   = numpy.empty ( ( pStruct.xDimIndices * ( pStruct.yDimIndices - 3 ) * pStruct.zDimIndices ) )
 newMapArr3D                                   = numpy.empty ( [ pStruct.xDimIndices,  ( pStruct.yDimIndices - 3 ),  pStruct.zDimIndices ] )
 
 """ Fill in the new arrays using the previously obtained ProSHADE internal maps """
 for xIt in range( 0, pStruct.xDimIndices ):
     for yIt in range( 0, pStruct.yDimIndices ):
         for zIt in range( 0, pStruct.zDimIndices ):
             if yIt >= ( pStruct.yDimIndices - 3 ):
                 continue
             arrPos                            = zIt + pStruct.zDimIndices * ( yIt  + pStruct.yDimIndices * xIt );
             newMapPos                         = zIt + pStruct.zDimIndices * ( yIt  + (pStruct.yDimIndices-3) * xIt );
             newMapArr1D[newMapPos]            = initialMapArray1D[arrPos] * 2
             newMapArr3D[xIt][yIt][zIt]        = initialMapArray1D[arrPos] / 2

 """ Now change the ProSHADE_data structure to reflect the new map """
 pStruct.yDimSize                              = pStruct.yDimSize - ( ( pStruct.yDimSize / pStruct.yDimIndices ) * 3 )
 pStruct.yDimIndices                           = pStruct.yDimIndices - 3
 pStruct.yGridIndices                          = pStruct.yDimIndices
 pStruct.yTo                                   = pStruct.yTo - 3

 """ And now force ProSHADE to change the internal map """
 proshade.setNewMapPython1D                    ( pStruct, newMapArr1D )
 proshade.setNewMapPython3D                    ( pStruct, newMapArr3D )
 \endcode
 *
 * \e Initial \e map \e procesing
 *
 * Once the map is read into the \b ProSHADE_data object, it needs to be processed in order to make sure ProSHADE will be able to use it for any further computations. While processing, ProSHADE offers the following map modifications through the \b ProSHADE_setting object variables: map invertion (this
 * will invert the map indices along each dimension), map normalisation (making the map density values mean 0 and standard deviation 1), map masking (computing a map mask by blurring and then setting mask as all values above threshold), map centering (moving the map into the centre of mass), adding extra
 * space (in case the map density is close to map edge, what can lead to map artefacts and lower accuracy of further processing) and map phase removal (removing the phase of the map density, effectively producing Patterson maps). The user can chose any, all or none of these, but the processing function needs
 * to be called before any further processing is possible. The following example code showcases how some of the processing functionalities can be chosen and how the map can be processed.
 *
 * \code{.py}
 """ Set which processing should be done """
 pSet.setMapInversion                          ( False )  # Should all map positions x,y,z be swapped to -x,-y,-z? Use this only if your helices have the wrong hand ...
 pSet.setNormalisation                         ( False )  # Should internal map representation be normalised to mean 0 and standard deviation 1?
 pSet.setMasking                               ( False )  # Should maps be masked by blurring?
 pSet.setMapCentering                          ( False )  # Move structure COM to the centre of map box?
 pSet.setExtraSpace                            ( 10.0 )   # Extra space in Angs to be added when creating internap map representation.
 pSet.setPhaseUsage                            ( True )   # Use full maps, or Patterson-like maps?
 
 """ Process the internal map representation """
 pStruct.processInternalMap                    ( pSet )
 \endcode
 *
 * \e Computing \e standard \e ProSHADE \e tasks
 *
 * If the user now wants to use ProSHADE to compute some of the standard ProSHADE taks, \e i.e. Distances computation, Symmetry detection, Re-boxing or Map overlay, it is recommended that the user proceeds in the same fashion as shown in the the \e advancedAccess_... example files. Moreover, these are also
 * demonstrated in the \b directAccess.py file available in the examples folder. Therefore, none of these tasks will be shown here in a step-wise manner; instead, the rest of this tutorial will focus on how partial information and results can be obtained from ProSHADE, should the user want to use them in a way that is not
 * currently supported by ProSHADE.
 *
 * \e Computing \e the \e spherical \e harmonics \e decomposition
 *
 * ProSHADE can compute the spherical harmonics decomposition of the internal map. However, instead of using the spherical-Bessel functions, it firstly creates a set of concentric spheres centered on the centre of indices (xDimIndices/2, yDimIndices/2, zDimIndices/2) point and spaced 2 indices apart, then it maps
 * the density map data onto these spheres and then it computes the spherical harmonics decomposition on each of these spheres independently. There is quite a few settings that relate to the spherical harmonics decompostion computation, such as the bandwidth of the computation, the sphere placement and spacing,
 * the resolution on the spheres, etc.; these arre mostly inter-related and ProSHADE will set them up automatically, unless the user specifies otherwise. Since these are quite technical, the interested users are referred to the second chapter of my Ph.D. thesis, which specifies all the technical details:
 * https://www.repository.cam.ac.uk/handle/1810/284410 . To issue this computation, please use the functions shown in the following example code:
 *
 * \code{.py}
 """ Map the internal map representation onto a set of concentric spheres """
 pStruct.mapToSpheres                          ( pSet )
 
 """ Compute the spherical harmonics decomposition """
 pStruct.computeSphericalHarmonics             ( pSet )
 \endcode
 *
 * If the user is interested in the spherical harmonics values (and possibly does not need any further computations from ProSHADE), these can be accessed using the function showcased below. It is worth noting that the organisation of the spherical harmonics is as follows: Each concentric shell has a two dimensional array of values,
 * where the row is the spherical harmonics band, while the column is the spherical harmonics order. Please note that there are only 2 * band + 1 orders for each band and therefore the array is not rectangular. Instead, it is recommended that the users take advantage of the supplied \e sphericalHarmonicsIndex() function, which takes
 * the order, the band and the shell number as its arguments (in this order) and returns the index of the spherical harmonics value in the retrieved spherical harmonics array. Moreover, please note that the spherical harmonics are complex numbers.
 *
 * \code{.py}
 """ Obtain all the spherical harmonics values for all shells """
 sphericalHarmonics                            = proshade.getSphericalHarmonics ( pStruct )

 """ Retrieve s specific value for shell 3, band 4 and order -2 """
 Shell3Band4OrderMin2Value                     = sphericalHarmonics[3][ pStruct.sphericalHarmonicsIndex ( -2, 4, 3 ) ] # Order -2, band 4, shell 3.
 \endcode
 *
 * \e Computing \e the \e self-rotation \e function
 *
 * ProSHADE also allows computing the self-rotation function. More specifically, it firstly computes the so called E matrices, which are matrices of the integral over all the concentric spheres of the spherical harmonics coefficients of order1 and order2, or in mathematical (LaTeX) form: Integral _0 ^rMAX ( c^lm * c'^lm ). It then
 * proceeds to normalise these E matrices, resulting in the SO(3) decomposition (Wigner D based decomposition) coefficients. Finally, by computing the inverse SO(3) Fourier transform (SOFT) on these coefficients, ProSHADE obtains the self-rotation function. In order to isue this computaion, the following code can be used:
 *
 * \code{.py}
 """ Compute the self-rotation function """
 pStruct.getRotationFunction                   ( pSet )
 \endcode
 *
 * Once the self-rotation function is computed, ProSHADE allows the user to access all of its interim results as well as the rotation function map. Specifically the E matrices, which are ordered by the band, order1 and order2 (in this order) can be obtained as shown in the following code. The E matrices are 3D arrays, which suffer from the
 * different number of orders for different bands feature of spherical harmonics. Therefore, the band dimensions of the arrays are zero padded; furthermore, as the order indexing goes from -band to +band, but the array indexing starts from zero, the correction to the array indexes is necessary. Regarding the SO(3) coefficients, they have
 * the same technical structure as the E matrices; however, due to some development issues, accessing them is done in a different manner, using the \e so3CoeffsArrayIndex() function as shown in the example code below. Finally, the self-rotation function map can be accessed as a 1D or 3D numpy.array using the functions shown in the
 * example code below. Again, as passing 3D matrices is not possible using SWIG and Numpy.i typedefs, the 3D version of the function is much slower than the 1D version. Assuming the user would find a point in the map for which he would like to know the rotation matrix (the indices of the self-rotation function are related to the Euler
 * angles in a non-trivial manner), ProSHADE provides the \e getRotationMatrixFromRotFunIndices() function also shown in the following example code. This function returns the rotation matrix belonging to the given self-rotation function indices in a numpy.array format.
 *
 * \code{.py}
 """ Obtain the E matrices """
 eMat                                          = proshade.getEMatrix ( pStruct )
 Band4OrderOneMin2OrderTwo3EMatrixValue        = eMat[4][-2+4][3+4] # Band = 4, Order1 = -2 and Order2 = 3
 
 """ Obtain the SO(3) coefficients """
 so3Coeffs                                     = proshade.getSO3Coeffs( pStruct )
 so3CoeffsOrderOneMin1OrderTwo3Band5           = so3Coeffs[pStruct.so3CoeffsArrayIndex ( -1, 3, 5 )] # Order1 = -1; Order2 = 3, Band = 5
 
 """ Obtain the self-rotation function """
 selfRotationFunction1D                        = proshade.getRotationFunction1D ( pStruct )
 selfRotationFunction3D                        = proshade.getRotationFunction3D ( pStruct )
 
 """ Convert self-rotation function indices to rotation matrix """
 rotMat                                        = proshade.getRotationMatrixFromRotFunIndices ( pStruct, 10, 11, 7 )
 \endcode
 *
 * \e Computing \e the \e optimal \e rotation \e function
 *
 * A related ProSHADE functionality is the computation of an optimal rotation function for two input structures. In the standard ProSHADE tasks, this is done for two phase-less structure maps (the phase is removed to achive identical centering on the maps) in order to find the optimal rotation, which overlays the two maps, but the user is free to
 * call this function for any two \b ProSHADE_data objects which both have their spherical harmonics values computed. To do this, we will create two new \b ProSHADE_data objects, read in some structures, process them, map them onto spheres, compute their spherical harmonics values and then we call the
 * \e getOverlayRotationFunction(). This function works similarly to the \e getRotationFunction() used above, but it uses spherical harmonics coefficients from two different structures as opposed to the same structures.
 *
 * \code{.py}
 """ Modify the settings object for optimal rotation function computation """
 pSet.task                                     = proshade.OverlayMap
 pSet.verbose                                  = 1
 pSet.requestedResolution                      = 8.0;
 pSet.usePhase                                 = False;
 pSet.changeMapResolution                      = True;
 pSet.maskMap                                  = False;
 pSet.moveToCOM                                = False;
 pSet.normaliseMap                             = False;
 
 """ Create the two new ProSHADE_data objects """
 pStruct_static                                = proshade.ProSHADE_data ( pSet )
 pStruct_moving                                = proshade.ProSHADE_data ( pSet )
 
 """ Read in two structures """
 pStruct_static.readInStructure                ( "./test1_rotTrs.map", 0, pSet )
 pStruct_moving.readInStructure                ( "./test1.pdb", 1, pSet )

 """ Process, map and get spherical harmonics for both structures """
 pStruct_static.processInternalMap             ( pSet )
 pStruct_moving.processInternalMap             ( pSet )

 pStruct_static.mapToSpheres                   ( pSet )
 pStruct_moving.mapToSpheres                   ( pSet )

 pStruct_static.computeSphericalHarmonics      ( pSet )
 pStruct_moving.computeSphericalHarmonics      ( pSet )
 
 """ Compute the optimal rotation function """
 pStruct_moving.getOverlayRotationFunction     ( pSet, pStruct_static )
 \endcode
 *
 * Now, in order to access the optimal rotation function, the user can use the same functions as for accessing the self-rotation function above, \e i.e. the \e getRotationFunction1D() and the \e getRotationFunction3D() functions, both with the moving structure as their only argument. Moreover, the same function for converting
 * the rotation function indices to the appropriate rotation matrix can also be used as shown in the following example code.
 *
 * \code{.py}
 """ Obtain the self-rotation function """
 rotationFunction1D                            = proshade.getRotationFunction1D ( pStruct_moving )
 rotationFunction3D                            = proshade.getRotationFunction3D ( pStruct_moving )
 
 """ Convert self-rotation function indices to rotation matrix """
 rotMat                                        = proshade.getRotationMatrixFromRotFunIndices ( pStruct_moving, 10, 11, 7 )
 \endcode
 *
 * \e Finding \e the \e optimal \e rotation
 *
 * Once the optimal rotation map is computed, the user may be interested in the highest value in the map and the corresponding rotation matrix (or Euler angles), as these will represent the rotation which overlays most of the two structures (within the error of the map
 * sampling). To facilitate this taks, ProSHADE contains the \e getBestRotationMapPeaksEulerAngles() function, which finds the highest peak in the map and returns the associated Euler angles. The following example code demonstrates how to use this function as well
 * as how to obtain the the rotation matrix from the Euler angles using ProSHADE.
 *
 * \code{.py}
 """ Find the highest peak in the map, associated Euler angles and rotation matrix """
 optimalRotationAngles                         = pStruct_moving.getBestRotationMapPeaksEulerAngles ( pSet )
 optimalRotationMatrix                         = proshade.getRotationMatrixFromEulerZXZ ( optimalRotationAngles )
 \endcode
 *
 * \e Rotating \e the \e internal \e map \e representation
 *
 * Once the optimal rotation angles are obtained, it is the next logical step to rotate the structure by these angles to get the two structures in identical orientation. This can also be done with ProSHADE function \e rotateMap(), which works with the Euler angles as
 * reported by ProSHADE. The rotation is done using the spherical harmonics coefficients, which are multiplied by the Wigner D matrices for the required rotation and the resulting rotated coefficients are then inverted back and interpolated to a new map. This process
 * has two side effects: Firstly, the resulting maps tend to suffer from minor artefacts resulting from the sequence termination errors and the interpolation to and from spheres to cartesian co-ordinates. And secondly, the input maps need to have their spherical harmonics
 * coefficients computed. Therefore, this approach is not recommended for any maps that are to be deposited or fitted into, but they are sufficient for computation of most ProSHADE standard tasks as the shape is largely identical.
 *
 * In terms of this tutorial, since we have already computed the optimal rotation between two structures, we will continue to show how this result can be used to rotate a new structure. This will allow us to demonstrate the next functionality of ProSHADE in the later sections
 * of this tutorial in a more streamlined fashion. To cause \b ProSHADE_data map rotation, the function in the example code can be used.
 *
 * \code{.py}
 """ Delete the old structure objects so that new can be created """
 del pStruct_static
 del pStruct_moving
 
 """ Change the settings to use phased maps """
 pSet.usePhase                                 = True
 pSet.changeMapResolution                      = True
 
 """ Create the two new ProSHADE_data objects """
 pStruct_static                                = proshade.ProSHADE_data ( pSet )
 pStruct_moving                                = proshade.ProSHADE_data ( pSet )
 
 """ Read in two structures """
 pStruct_static.readInStructure                ( "./test1_rotTrs.map", 0, pSet )
 pStruct_moving.readInStructure                ( "./test1.pdb", 1, pSet )

 """ Process the static structure to allow further examples """
 pStruct_static.processInternalMap             ( pSet )

 """ Process the moving structure and compute the spherical harmonics to allow rotation - this is not needed for the static structure """
 pStruct_moving.processInternalMap             ( pSet )
 pStruct_moving.mapToSpheres                   ( pSet )
 pStruct_moving.computeSphericalHarmonics      ( pSet )
 
 """ Rotate the moving structure """
 pStruct_moving.rotateMap                      ( pSet, optimalRotationAngles[0], optimalRotationAngles[1], optimalRotationAngles[2] )
 \endcode
 *
 * \e Computing \e the \e translation \e function
 *
 * Similarly to the rotation function, the user may be interested in the optimal translation required to overlay two structures. ProSHADE can compute such an optimal translation using the translation function; however, in order to compute it, it requires that the two internal map representation have the same
 * dimensions in terms of map indices and map sampling; identical map sampling is achieved by setting the \b changeMapResolution setting to true. Still, the identical number of indices will not generally be the case, ProSHADE provides a padding function, which can add zeroes around the internal
 * representation map to make sure that it has given dimensions. Therefore, in order to compute the translation function, it is required that the two structures are modified by the \e zeroPaddToDims() function to both have the same dimensions; the higher of the two structures are chosen in order to avoid
 * loss of information.
 *
 * \code{.py}
 """ Add zeroes around he structure to achieve given number of indicel along each dimension """
 pStruct_static.zeroPaddToDims                 ( int ( numpy.max ( [ pStruct_static.getXDim(), pStruct_moving.getXDim() ] ) ),
                                                 int ( numpy.max ( [ pStruct_static.getYDim(), pStruct_moving.getYDim() ] ) ),
                                                 int ( numpy.max ( [ pStruct_static.getZDim(), pStruct_moving.getZDim() ] ) ) )
 pStruct_moving.zeroPaddToDims                 ( int ( numpy.max ( [ pStruct_static.getXDim(), pStruct_moving.getXDim() ] ) ),
                                                 int ( numpy.max ( [ pStruct_static.getYDim(), pStruct_moving.getYDim() ] ) ),
                                                 int ( numpy.max ( [ pStruct_static.getZDim(), pStruct_moving.getZDim() ] ) ) )
 \endcode
 *
 * Once the structures have the same dimensions, it is possible to compute the translation function. This function will compute the Fourier transforms of both maps, combine the Fourier coefficients and compute the inverse Fourier transform on the resulting combined coefficients map, thus obtaining the
 * translation map. Once computed, this map can be accessed from the ProSHADE python module as shown in the following example code, again keeping in mind that the 3D version takes considerably longer to obtain than the 1D version.
 *
 * \code{.py}
 """ Compute the translation function """
 pStruct_moving.computeTranslationMap          ( pStruct_static )
 
 """ Access the translation map as 1D or 3D numpy.array """
 translationMap1D                              = proshade.getTranslationFunction1D ( pStruct_moving )
 translationMap3D                              = proshade.getTranslationFunction3D ( pStruct_moving )
 \endcode
 *
 * Also, similarly to the rotation function, ProSHADE provides a useful function for detecting the highest peak in the translation map and computing the corresponding translation in Angstroms. This is then demonstrated in the following example code:
 *
 * \code{.py}
 """ Find the optimal translation vector """
 optimalTranslationVector                      = pStruct_moving.getBestTranslationMapPeaksAngstrom ( pStruct_static )
 \endcode
 *
 *
 * \e Translating \e the \e internal \e representation
 *
 * Once the optimal translation vector is computed, it makes sense that ProSHADE should also be able to apply it to the internal map representation. Therefore, the function \e translateMap() is provided to facilitate this task. The translation is done in two steps, firstly, ProSHADE
 * simply modifies the starting indices and axes origins of the map to minimise the movement of the map in the cell by moving the cell as a whole. Next, the remaining translation is then done in the frequency domain (by modifing the Fouries coefficients) of the internal representation.
 *
  * \code{.py}
 """ Translate the internal representation """
 pStruct_moving.translateMap                   ( pSet, optimalTranslationVector[0], optimalTranslationVector[1], optimalTranslationVector[2] )
 \endcode
 *
 * \e Writing \e out \e resulting \e structures
 *
 * Finally, it is worth noting that while the MAP formatted data can be written out of the \b ProSHADE_data object at any time (albeit their quality may be decreased if the rotation was applied as discussed in the rotating internal representation map section), ProSHADE can also write
 * out the co-ordinate data for input structures, which were read in from a co-ordinate file. Please note that ProSHADE cannot generate co-ordinate data from maps, the co-ordinate data need to pre-exist ProSHADE run. Nonetheless, in the case of, for example, finding the optimal rotation
 * and translation of one structure to overlay with another structure, the user may be interested in writing out the modified co-ordinates. To do this, ProSHADE contains the \e writePdb() function, which needs to be supplied with the file name, the required rotation and translation and it
 * will write out the PDB file with these modifications applied.
 *
 * Also, please note that it is the users responsibility to add any rotations or translations together and to supply this function with the correct cumulative values. The example code below shows how a rotated and translated PDB file can be outputted by ProSHADE.
 *
 * \code{.py}
 """ Translate the internal representation """
 pStruct_moving.writePdb                       ( "overlayed.pdb",
                                                 optimalRotationAngles[0],
                                                 optimalRotationAngles[1],
                                                 optimalRotationAngles[2],
                                                 optimalTranslationVector[0],
                                                 optimalTranslationVector[1],
                                                 optimalTranslationVector[2] )
 \endcode
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
