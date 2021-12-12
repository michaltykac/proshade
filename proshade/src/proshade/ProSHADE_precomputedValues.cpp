/*! \file ProSHADE_precomputedValues.cpp
    \brief This source file contains functions for the precomputed values classes.
 
    The pre-computed values are held in classes as C++ requires exit-time destructors to be provided for non-trivial types. This file then
    contains all the functions required to service these classes.
 
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
#include "ProSHADE_precomputedValues.hpp"

/*! \brief Constructor of the tetrahedronAxes class.
 */
ProSHADE_internal_precomputedVals::tetrahedronAxes::tetrahedronAxes ( )
{
    //================================================ Create the object with the pre-computed values
    this->tetrahedronAxesVals                         = std::vector < std::vector < proshade_double > > {
        { 3.0, 0.816496580927726,  0.000000000000000, -0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 1
        { 3.0,-0.816496580927726,  0.000000000000000, -0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 2
        { 3.0, 0.000000000000000,  0.816496580927726,  0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 3
        { 3.0, 0.000000000000000, -0.816496580927726,  0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 4
        
        { 2.0, 0.000000000000000,  0.000000000000000, -1.000000000000000, M_PI, 0.0 },                     // C2 axis 1
        { 2.0,-0.707106781186547,  0.707106781186547,  0.000000000000000, M_PI, 0.0 },                     // C2 axis 2
        { 2.0, 0.707106781186547,  0.707106781186547,  0.000000000000000, M_PI, 0.0 },                     // C2 axis 3
    };
    
    //================================================ Done
    
}

/*! \brief Destructor of the tetrahedronAxes class.
 */
ProSHADE_internal_precomputedVals::tetrahedronAxes::~tetrahedronAxes ( )
{
    //================================================ Delete the pre-computed values
    for ( size_t axIt = 0; axIt < this->tetrahedronAxesVals.size(); axIt++ )
    {
        this->tetrahedronAxesVals.at(axIt).clear      ( );
    }
    this->tetrahedronAxesVals.clear                   ( );
    
    //================================================ Done
    
}

/*! \brief Accessor for the tetrahedronAxesVals variable.
 
    This function is what gives access to the tetrahedronAxesVals variable and the meaning behind this class
 
    \param[in] axis The axis index of the axis for which value should be retrieved.
    \param[in] element The element index for the value in the axis to be retrieved.
    \param[out] val The value of the requested element of the requested axis.
 */
proshade_double ProSHADE_internal_precomputedVals::tetrahedronAxes::getValue ( proshade_unsign axis, proshade_unsign element )
{
    //================================================ Done
    return                                            ( this->tetrahedronAxesVals.at(axis).at(element) );
    
}

/*! \brief Accessor for the tetrahedronAxesVals variable number of axes.
 
    \param[out] size The number of axes in the tetrehadronAxesVals variable.
 */
proshade_unsign ProSHADE_internal_precomputedVals::tetrahedronAxes::getNoAxes ( )
{
    //================================================ Done
    return                                            ( static_cast< proshade_unsign > ( this->tetrahedronAxesVals.size() ) );
    
}

/*! \brief Constructor of the octahedronAxes class.
 */
ProSHADE_internal_precomputedVals::octahedronAxes::octahedronAxes ( )
{
    //================================================ Create the object with the pre-computed values
    this->octahedronAxesVals                          = std::vector < std::vector < proshade_double > > {
        { 4.0, 1.000000000000000,  0.000000000000000,  0.000000000000000, (2.0 * M_PI) / 4.0, 0.0 },       // C4 axis 1
        { 4.0, 0.000000000000000,  1.000000000000000,  0.000000000000000, (2.0 * M_PI) / 4.0, 0.0 },       // C4 axis 2
        { 4.0, 0.000000000000000,  0.000000000000000,  1.000000000000000, (2.0 * M_PI) / 4.0, 0.0 },       // C4 axis 3
    
        { 3.0, 0.577350269189626,  0.577350269189626,  0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 1
        { 3.0,-0.577350269189626,  0.577350269189626,  0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 2
        { 3.0, 0.577350269189626,  0.577350269189626, -0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 3
        { 3.0,-0.577350269189626,  0.577350269189626, -0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 4
        
        { 2.0, 0.707106781186547,  0.707106781186547,  0.000000000000000, M_PI, 0.0 },                     // C2 axis 1
        { 2.0, 0.000000000000000,  0.707106781186547,  0.707106781186547, M_PI, 0.0 },                     // C2 axis 2
        { 2.0, 0.707106781186547,  0.000000000000000,  0.707106781186547, M_PI, 0.0 },                     // C2 axis 3
        { 2.0,-0.707106781186547,  0.707106781186547,  0.000000000000000, M_PI, 0.0 },                     // C2 axis 4
        { 2.0,-0.707106781186547,  0.000000000000000,  0.707106781186547, M_PI, 0.0 },                     // C2 axis 5
        { 2.0, 0.000000000000000,  0.707106781186547, -0.707106781186547, M_PI, 0.0 }                      // C2 axis 6
    };
    
    //================================================ Done
    
}

/*! \brief Destructor of the octahedronAxes class.
 */
ProSHADE_internal_precomputedVals::octahedronAxes::~octahedronAxes ( )
{
    //================================================ Delete the pre-computed values
    for ( size_t axIt = 0; axIt < this->octahedronAxesVals.size(); axIt++ )
    {
        this->octahedronAxesVals.at(axIt).clear       ( );
    }
    this->octahedronAxesVals.clear                    ( );
    
    //================================================ Done
    
}

/*! \brief Accessor for the octahedronAxesVals variable.
 
    This function is what gives access to the octahedronAxesVals variable and the meaning behind this class
 
    \param[in] axis The axis index of the axis for which value should be retrieved.
    \param[in] element The element index for the value in the axis to be retrieved.
    \param[out] val The value of the requested element of the requested axis.
 */
proshade_double ProSHADE_internal_precomputedVals::octahedronAxes::getValue ( proshade_unsign axis, proshade_unsign element )
{
    //================================================ Done
    return                                            ( this->octahedronAxesVals.at(axis).at(element) );
    
}

/*! \brief Accessor for the octahedronAxesVals variable number of axes.
 
    \param[out] size The number of axes in the octahedronAxesVals variable.
 */
proshade_unsign ProSHADE_internal_precomputedVals::octahedronAxes::getNoAxes ( )
{
    //================================================ Done
    return                                            ( static_cast< proshade_unsign > ( this->octahedronAxesVals.size() ) );
    
}


/*! \brief Constructor of the icosahedronAxes class.
 */
ProSHADE_internal_precomputedVals::icosahedronAxes::icosahedronAxes ( )
{
    //================================================ Create the object with the pre-computed values
    this->icosahedronAxesVals                         = std::vector < std::vector < proshade_double > > {
        { 5.0, 0.850650808352040,  0.525731112119134,  0.000000000000000, (2.0 * M_PI) / 5.0, 0.0 },       // C5 axis 1
        { 5.0, 0.850650808352040, -0.525731112119134,  0.000000000000000, (2.0 * M_PI) / 5.0, 0.0 },       // C5 axis 2
        { 5.0, 0.525731112119134,  0.000000000000000,  0.850650808352040, (2.0 * M_PI) / 5.0, 0.0 },       // C5 axis 3
        { 5.0, 0.525731112119134,  0.000000000000000, -0.850650808352040, (2.0 * M_PI) / 5.0, 0.0 },       // C5 axis 4
        { 5.0, 0.000000000000000,  0.850650808352040,  0.525731112119134, (2.0 * M_PI) / 5.0, 0.0 },       // C5 axis 5
        { 5.0, 0.000000000000000,  0.850650808352040, -0.525731112119134, (2.0 * M_PI) / 5.0, 0.0 },       // C5 axis 6
        
        { 3.0, 0.577350269189626,  0.577350269189626,  0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 1
        { 3.0, 0.577350269189626,  0.577350269189626, -0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 2
        { 3.0, 0.577350269189626, -0.577350269189626,  0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 3
        { 3.0, 0.577350269189626, -0.577350269189626, -0.577350269189626, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 4
        { 3.0, 0.356822089773090,  0.934172358962716,  0.000000000000000, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 5
        { 3.0,-0.356822089773090,  0.934172358962716,  0.000000000000000, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 6
        { 3.0, 0.934172358962716,  0.000000000000000,  0.356822089773090, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 7
        { 3.0, 0.934172358962716,  0.000000000000000, -0.356822089773090, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 8
        { 3.0, 0.000000000000000,  0.356822089773090,  0.934172358962716, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 9
        { 3.0, 0.000000000000000,  0.356822089773090, -0.934172358962716, (2.0 * M_PI) / 3.0, 0.0 },       // C3 axis 10
        
        { 2.0, 0.500000000000000,  0.809016994374947,  0.309016994374947, M_PI, 0.0 },                     // C2 axis 1
        { 2.0, 0.309016994374947,  0.500000000000000,  0.809016994374947, M_PI, 0.0 },                     // C2 axis 2
        { 2.0, 0.809016994374947,  0.309016994374947,  0.500000000000000, M_PI, 0.0 },                     // C2 axis 3
        { 2.0, 0.809016994374947,  0.309016994374947, -0.500000000000000, M_PI, 0.0 },                     // C2 axis 4
        { 2.0, 0.309016994374947,  0.500000000000000, -0.809016994374947, M_PI, 0.0 },                     // C2 axis 5
        { 2.0, 0.500000000000000,  0.809016994374947, -0.309016994374947, M_PI, 0.0 },                     // C2 axis 6
        { 2.0, 0.809016994374947, -0.309016994374947,  0.500000000000000, M_PI, 0.0 },                     // C2 axis 7
        { 2.0, 0.309016994374947, -0.500000000000000,  0.809016994374947, M_PI, 0.0 },                     // C2 axis 8
        { 2.0, 0.500000000000000, -0.809016994374947,  0.309016994374947, M_PI, 0.0 },                     // C2 axis 9
        { 2.0, 0.500000000000000, -0.809016994374947, -0.309016994374947, M_PI, 0.0 },                     // C2 axis 10
        { 2.0, 0.309016994374947, -0.500000000000000, -0.809016994374947, M_PI, 0.0 },                     // C2 axis 11
        { 2.0, 0.809016994374947, -0.309016994374947, -0.500000000000000, M_PI, 0.0 },                     // C2 axis 12
        { 2.0, 0.500000000000000,  0.809016994374947, -0.309016994374947, M_PI, 0.0 },                     // C2 axis 13
        { 2.0, 0.000000000000000,  1.000000000000000,  0.000000000000000, M_PI, 0.0 },                     // C2 axis 14
        { 2.0, 0.500000000000000,  0.809016994374947,  0.309016994374947, M_PI, 0.0 }                      // C2 axis 15
    };
    
    //================================================ Done
    
}

/*! \brief Destructor of the icosahedronAxes class.
 */
ProSHADE_internal_precomputedVals::icosahedronAxes::~icosahedronAxes ( )
{
    //================================================ Delete the pre-computed values
    for ( size_t axIt = 0; axIt < this->icosahedronAxesVals.size(); axIt++ )
    {
        this->icosahedronAxesVals.at(axIt).clear      ( );
    }
    this->icosahedronAxesVals.clear                   ( );
    
    //================================================ Done
    
}

/*! \brief Accessor for the icosahedronAxesVals variable.
 
    This function is what gives access to the icosahedronAxesVals variable and the meaning behind this class
 
    \param[in] axis The axis index of the axis for which value should be retrieved.
    \param[in] element The element index for the value in the axis to be retrieved.
    \param[out] val The value of the requested element of the requested axis.
 */
proshade_double ProSHADE_internal_precomputedVals::icosahedronAxes::getValue ( proshade_unsign axis, proshade_unsign element )
{
    //================================================ Done
    return                                            ( this->icosahedronAxesVals.at(axis).at(element) );
    
}

/*! \brief Accessor for the octahedronAxesVals variable number of axes.
 
    \param[out] size The number of axes in the octahedronAxesVals variable.
 */
proshade_unsign ProSHADE_internal_precomputedVals::icosahedronAxes::getNoAxes ( )
{
    //================================================ Done
    return                                            ( static_cast< proshade_unsign > ( this->icosahedronAxesVals.size() ) );
    
}
