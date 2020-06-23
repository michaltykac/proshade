ProSHADE
========

Protein Shape Description and Symmetry Detection

# Introduction

ProSHADE is a C++ language library and an associated tool providing functionalities for working with structural biology molecular structures. The library implements functions for computing shape-wise structural distances between pairs of molecules, detecting symmetry over the centre of mass of a single structure, map re-sizing as well as matching density maps and PDB coordinate files into one another. The executable implemented in the bin.cpp file then allows easy access to these functionalities without the need for library linking, while the python modules provide easy access to the functionality from the python language. For help on how the executable should be used, refer to the -h option of it. For more details about the functionalities, see below

# Index

- [Introduction](#introduction)
- [Index](#index)
- [Installation](#installation)
    - [Standard System Dependencies](#standard-system-dependencies)
    - [CMake options](#cmake-options)
    - [Other dependencies](#other-dependencies)
    - [Install](#install)
    - [Uninstall](#Uninstall)
- [Using the ProSHADE binary](#using-the-proshade-binary)
    - [Symmetry Detection](#symmetry-detection)
    - [Shape similarity distances](#shape-similarity-distances)
    - [Re-boxing structures](#re-boxing-structures)
    - [Optimal rotation and translation](#optimal-rotation-and-translation)
- [Using the ProSHADE library](#using-the-proshade-library)

# Installation

The installation of the ProSHADE software should be done using the CMake system and the supplied CMakeLists.txt file. The minimual requiered version of CMake is 2.6, however, python modules and single source file compilation will not be available unless CMake version 3.1 or higher is used. The CMakeLists.txt file assumes the standard system dependencies are installed in the system folders; for a full list of standard system dependencies, please see the section [Standard System Dependencies](#standard-system-dependencies).
 
 Once all of the standard system dependencies are installed CMake can be run to create the make files. There are several options that can be used to modify the default behaviour of the installation; these typically drive the installation locations and dependencies paths in the case of non-standard dependency location. Please see the section [CMake options](#cmake-options) for details as to how to use these options and what do they do.
 
 Please note that while the ProSHADE code is C++ 98 standard compatible, some of the dependencies do require at least partial support for the C++ 11 standard.
 
 ## Standard System Dependencies
 
 Generally, the following list of standard system libraries and utilities are required for successfull installation of ProSHADE on Linux systems. On MacOS systems, most of these should be installed by default except where specifically stated:
  - **gcc**
  - **g++** (on Ubuntu and Debian) or **gcc**-c++ (on CentOS and SuSe)
  - **gfortran** (on Ubuntu and Debian) or **gcc-gfortran** (on CentOS and SuSe )
  - **make**
  - **cmake**
  - **m4**
  - **fftw3-dev** (on Ubuntu, Debian and SuSe) or **fftw3-devel** (on CentOS)
  - **libblas-dev** (on Ubuntu and Debian) or **blas-devel** (on CentOS and SuSe)
  - **liblapack-dev** (on Ubuntu and Debian) or **lapack-devel** (on CentOS) or **lapack-dev** (on SuSe)
  - **python** (on Ubuntu, Debian and SuSe) or **python2** (on CentOS)
  - **python-pip** (on Ubuntu, Debian and SuSe) or **python2-pip** (on CentOS)
  - **python-dev** (on Ubuntu, Debian and SuSe) or **python2-devel** (on CentOS)
  - **python3**
  - **python3-pip**
  - **python3-dev** (on Ubuntu, Debian and SuSe) or **python3-devel**
  - **swig**
  - **git**
  - **numpy** (installed using pip or pip3 separately for python2.x and python3.x)
 
 ## CMake options
 
  **-DINSTALL_LOCALLY=ON** or **OFF**
  - This option is used to decide whether all the installed ProSHADE components are installed in the local source directory (value \b ON ) or whether they are instead installed in the system folders (value \b OFF ). This option applies to the binary, the C++ library, the python2 and python3 modules (which are installed in the appropriate site-packages folder if the option is \b OFF ) and the headers as well.
 
  **-DINSTALL_BIN_DIR=/path**
  - This option is used to manually specify the folder to which the ProSHADE binary shold be installed into.
 
  **-DINSTALL_LIB_DIR=/path**
  - This option is used to manually specify the folder to which the ProSHADE C++ library should be installed into.
 
  **-DINSTALL_INC_DIR=/path**
  - This option is used to specify the folder to which the ProSHADE header files required by the library should be installed into.
 
  **-DCUSTOM_FFTW3_LIB_PATH=/path**
  - This option is used to supply the path to the libfftw3.a/so/dylib in the case where ProSHADE CMake installation fails to detect the FFTW3 dependency. This is typically the case when FFTW3 is installed outside of the standard FFTW3 installation locations.
 
  **-CUSTOM_FFTW3_INC_PATH=/path**
  - This option is used to supply the path to the fftw3.h file in the case where ProSHADE CMake installation fails to detect the FFTW3 dependency. This is typically the case when FFTW3 is installed outside of the standard FFTW3 installation locations.
 
  **-DCUSTOM_LAPACK_LIB_PATH=/path**
  - This option is used to supply the path to the liblapack.a/so/dylib in the case where ProSHADE CMake installation fails to detect the LAPACK dependency. This is typically the case when the LAPACK is installed outside of the standard LAPACK installation locations.
 
 ## Other dependencies
 
  ProSHADE also depends on the *libccp4*, *MMDB2*, *FFTW2*, *clipper* and *SOFT2.0* libraries. Since the installation of these libraries is non-trivial and does require some user input, these libraries are supplied with the ProSHADE code and will be installed locally by the ProSHADE CMake installation. Please note that these dependencies do have their own licences (the CCP4 licence, the GPL licence, ...) and therefore this may limit the ProSHADE usage for some users beyond the ProSHADE copyright and licence itself.
 
 ## Install
 
  In order to install ProSHADE, first please check that all the [Standard System Dependencies](#standard-system-dependencies) are installed, preferably using a package management system such as *apt* or *yum*. Next, please navigate to any folder to which you would like to write the install files; some find it useful to create a ```build``` folder in the ProSHADE folder in order to keep the install files in the same location as the source codes. Then, issue the following set of commands, setting the ```\path\to\ProSHADE``` to the correct path on your system and adding any required [CMake options](#cmake-options) to the first command. Please note that ```sudo``` may be required for the ```make install``` command if you are installing into the system folders.
 
 ```
  cmake \path\to\ProSHADE
  make
  make install
```
 
 ## Uninstall
 
  To remove the installed ProSHADE components, the command ```make remove``` needs to be issued to the makefile originally created by the CMake call. Please note that ```sudo``` may need to be used if the installation was done into the system folders and your current user does not have admin rights.
 
 # Using the ProSHADE binary

The ProSHADE tool was developed in a modular fashion and as the usage slightly changes depending on the functionality that is required. Nonetheless, care has been taken to make sure that identical or closely related features are controlled by the same command line arguments in all cases. Moreover, the GNU command line options standard have been adhered to (through the getOpts library) and therefore the users familiar with other command line tools should find the entering of command line arguments simple. The following subsections relate to examples of using different functionalities; for a full list of command line options, please use the **--help** command line option of the ProSHADE binary.

## Symmetry Detection

In order to detect symmetry in either a coordinate input file, or in a map input file, the ProSHADE executable needs to be supplied with the option **-S** or **--symmetry** and it will also require a single input file to be supplied using the **-f** option. These two options are the only mandatory options, although there are many optional values that the user can supply to supersede the default values and therefore modify the operation fo the ProSHADE executable to fit their purpose.

One particular option regarding the symmetry detection mode should be noted; the **--sym** (or **-u**) option allows the user to state which symmetry they believe to exist in the structure. The allowed values for this command line argument are "Cx", "Dx", "T", "O" and "I", where the *x* should be an integer number specifying the fold of the requested symmetry. When this option is used, it removes the default behaviour of returning the highest detected symmetry and instead the symmetry requested by the user is returned, if it can be found in the structure.

Another noteworthy option is the **--center** or **-c** option, which  tells ProSHADE to center the internal map representation over the centre of co-ordinates before running any processing of the map. This may be important as ProSHADE detects symmetries over the centre of the co-ordinates and therefore a non-centered map (map which does not have the centre of mass at the centre of co-ordinates) will be found to have no symmetries even if these are present, just not over the co-ordinate centre.

To demonstrate how the tool can be run and the standard output for the symmetry mode of operation, the current version of the ProSHADE executable was used to detect the
symmetry of a density map of the bacteriophage T4 portal protein with the PDB accession code 3JA7, which has the \a C12 symmetry. The visualisation of the structure is
shown in the following figure, while the output of the ProSHADE tool follows:

![T4 Portal Protein](https://github.com/michaltykac/proshade/blob/experimental/documentation/ProSHADE_3JA7.jpg)

```
 $: proshade -S -f ~/Downloads/3ja7.pdb -c --sym C12
 ProSHADE 0.7.3 (JUN 2020):
 ==========================

  ... Starting to read the structure: ./3ja7.pdb
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Centering map onto its COM.
  ... Adding extra 10 angstroms.
  ... Phase information retained in the data.
  ... Starting sphere mapping procedure.
  ... Preparing spherical harmonics environment.
  ... Starting spherical harmonics decomposition.
  ... Starting self-rotation function computation.
  ... Starting C symmetry detection.
 Detected C symmetry with fold 12 .

 ======================
 ProSHADE run complete.
 Time taken: 6 seconds.
 ======================
```

## Shape similarity distances

 The distances computation mode is signalled to the ProSHADE executable by the command line argument **-D** or **--distances**. This mode requires two or more structures to be supplied using the **-f** command line option. At least two structures are mandatory for the ProSHADE tool to proceed. Moreover, the resolution of the structures to which the comparison should be done needs to be supplied using the **-r** option. This resolution does not need to be the real resolution to which the structure(s) were solved, but rather reflects the amount of details which should be taken into accout when comparing shapes. Therefore, higher resolution comparison will focus more on details of the shapes, while lower resolution comparison will focus more on the overall shape ignoring the minor details. Please note that the results are calculated only for the first structure against all the remaining structures, **not** for all against all distance matrix.

 There are a number of useful options for the shape distances computation, please consult the **--help** dialogue for their complete listing.

 To demonstrate the output of the ProSHADE software tool for computing distances between structure shapes, the distances between the BALBES protein domains 1BFO_A_dom_1 and
 1H8N_A_dom_1 (which have similar shape) and the 3IGU_A_dom_1 domain which has a different shape, as can be seen from the following figure - the first two domains are
 both in cluster a), while the last domain is from the cluster b). The output of the ProSHADE software tool is then shown below:

![BALBES domain used for similarity detection](https://github.com/michaltykac/proshade/blob/experimental/documentation/ProSHADE_dists.png)

```
  $: proshade -D -f ./1BFO_A_dom_1.pdb -f ./1H8N_A_dom_1.pdb -f ./3IGU_A_dom_1.pdb -r 4
 ProSHADE 0.7.3 (JUN 2020):
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
 Energy levels distance    : 0.784872
 Trace sigma distance      : 0.805476
 Rotation function distance: 0.647008
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
 Energy levels distance    : 0.470989
 Trace sigma distance      : 0.53887
 Rotation function distance: 0.399456

 ======================
 ProSHADE run complete.
 Time taken: 0 seconds.
 ======================
```
 
 ## Re-boxing structures
 
 Another useful feature of the ProSHADE tool is the re-boxing of macromolecular density maps. This mode is signalled to the ProSHADE tool by the command line option **-M** or **--mapManip** followed by the **-R** option to specify that the required map manipulations include re-boxing. Furthermore, a single map structure file needs to be supplied after the **-f** flag. In this mode, ProSHADE will attempt to find a suitable map mask by blurring the map (increasing the overall B-factors). Consequently, it will use the map boundaries to create a new, hopefully smaller, box to which the appropriate part of the map will be copied.
 
 This ProSHADE functionality can be combinaed with other map manipulations, which include the map invertion (signalled by the **--invertMap** option and useful for cases where map reconstruction software mistakes the hands of the structure), the map normalisation (signalled by the **--normalise** option, which makes sure the map mean is 0 and standard deviation is 1), centering of centre of mass to the centre of co-ordinates (using the **--center** or **-c** option) or the phase removal (creating Patterson maps using the **--noPhase** or **-p** options).
 
 The location and filename of where this new map should be saved can be specified using the **--reBoxedFilename** (or **-g** ) command line option followed by the filename.
 
 The following snippet shows the output of the ProSHADE tool when used to re-box the TubZ-Bt four-stranded filament structure (EMD-5762), where the original volume can be decreased to 46.9% of the original structure volume and thus any linear processing of such structure will be more than twice faster and the original. The original TubZ-Bt four-stranded filament structure box is shown in the following figure as semi-transparent grey, while the new box is shown in non-transparent yellow.
 
![Re-boxing result for TubZ-Bt four-stranded filament](https://github.com/michaltykac/proshade/blob/experimental/documentation/ProSHADE_rebox.png)
 
```
 $ proshade --mapManip -R -f ./emd_5762.map
 ProSHADE 0.7.3 (JUN 2020):
 ==========================

  ... Starting to read the structure: /Users/mysak/Downloads/emd_5762.map
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
 Time taken: 48 seconds.
 ======================
```

## Optimal rotation and translation

 In order to find the rotation and translation which optimally overlays (or fits) one structure into another, be them PDB files or maps (and any combination thereof), the ProSHADE tool can be used in the Overlay mode. This is signalled to the ProSHADE tool binary by the command line option **--strOverlay** or the **-O** and this mode requires exactly two structure files to be supplied using the **-f** command line options. The order of the two files does matter, as the second file will always be moved to match the first structure, which will remain static.
 
 Due to the requirement for the second stucture movement and rotation, it is worth noting that the structure may need to be re-sampled and/or moved to the same viewing position as the first structure. This is done so that only the internal representation is modified, but never the input file. However, when the overlay structure is outputted (a non-default name can be specified by the **--overlayFile** command line option) the header of this output file may differ from the second structure header. Furthermore, if there is no extra space around the structure, movement and rotation may move pieces of the structure through the box boundaries to the other side of the box. To avoid this, please use the **--extraSpace** option to add some extra space around the structure.
 
 An example of the Overlay mode matching a single PDB structure (2A2Q_T_dom_2 from the BALBES database, original structure code 2A2Q) shown in part b) of the following figure to a density map computed with low resolution from this structure shown in part a) of the figure. Part c) then shows the match obtained by the internal map representations of both inputs, while part d) shows the final match of the moved and rotated PDB file to the original map input. The output and call of the ProSHADE tool is shown below. Please note that the optimal rotation matrix and translation vector are written into the output when verbosity (**--verbose**) is increased to at least 3, but are better accessed programatically (see the following sections) if you are interested in using these further.
 
 ![ProSHADE Overlay results for 2A2Q_T_dom_2.pdb](https://github.com/michaltykac/proshade/blob/experimental/documentation/ProSHADE_overlay.jpg)
 
```
 $ proshade -O -f ./2A2Q_T_dom_2_rotTrs.map -f ./2A2Q_T_dom_2.pdb --overlayFile moved.overlay -r 4 --extraSpace 25.0
 ProSHADE 0.7.3 (JUN 2020):
 ==========================

  ... Starting to read the structure: ./2A2Q_T_dom_2_rotTrs.map
  ... Starting to read the structure: ./2A2Q_T_dom_2.pdb
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Map centering not requested.
  ... Adding extra 25 angstroms.
  ... Centering map onto its COM.
  ... Phase information removed from the data.
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Map centering not requested.
  ... Adding extra 25 angstroms.
  ... Centering map onto its COM.
  ... Phase information removed from the data.
  ... Starting sphere mapping procedure.
  ... Preparing spherical harmonics environment.
  ... Starting sphere mapping procedure.
  ... Preparing spherical harmonics environment.
  ... Starting spherical harmonics decomposition.
  ... Starting spherical harmonics decomposition.
  ... Starting rotation function computation.
  ... Starting to read the structure: ./2A2Q_T_dom_2_rotTrs.map
  ... Starting to read the structure: ./2A2Q_T_dom_2.pdb
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Map centering not requested.
  ... Adding extra 25 angstroms.
  ... Phase information retained in the data.
  ... Map inversion (mirror image) not requested.
  ... Map normalisation not requested.
  ... Masking not requested.
  ... Map centering not requested.
  ... Adding extra 25 angstroms.
  ... Phase information retained in the data.
  ... Starting sphere mapping procedure.
  ... Preparing spherical harmonics environment.
  ... Starting spherical harmonics decomposition.
  ... Starting translation function computation.

 ======================
 ProSHADE run complete.
 Time taken: 17 seconds.
 ======================
```
# Using the ProSHADE library

Stay tuned for improved documentation comming soon!
