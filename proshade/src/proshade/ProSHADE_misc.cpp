/*! \file ProSHADE_misc.cpp
    \brief This source file contains all miscellaneous functions.
 
    The functions defined in here are used by ProSHADE to deal with minor tasks unrelated to anything
    in the main library of functions, for example C++ version independent vector addition code and such.
 
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
#include "ProSHADE_misc.hpp"

/*! \brief Adds the element to the vector.
 
    This function takes a pointer to a vector of strings and a single string element and adds this element to the end of
    the vector. The reason for this function is to make vector elongation C++ standard independent (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of strings which should be elongated.
    \param[in] elementToAdd String to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToStringVector ( std::vector < std::string >* vecToAddTo, std::string elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds the element to the vector.
 
    This function takes a pointer to a vector of proshade_single and a single proshade_single element and adds this element to the end of
    the vector. The reason for this function is to make vector elongation C++ standard independent (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of proshade_single's which should be elongated.
    \param[in] elementToAdd proshade_single to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToSingleVector ( std::vector < proshade_single >* vecToAddTo, proshade_single elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds the element to the vector.
 
    This function takes a pointer to a vector of proshade_double and a single proshade_double element and adds this element to the end of
    the vector. The reason for this function is to make vector elongation C++ standard independent (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of proshade_double's which should be elongated.
    \param[in] elementToAdd proshade_double to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToDoubleVector ( std::vector < proshade_double >* vecToAddTo, proshade_double elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds the element to the vector.
 
    This function takes a pointer to a vector of proshade_unsign and a single proshade_unsign element and adds this element to the end of
    the vector. The reason for this function is to make vector elongation C++ standard independent (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of proshade_unsign's which should be elongated.
    \param[in] elementToAdd proshade_unsign to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToUnsignVector ( std::vector < proshade_unsign >* vecToAddTo, proshade_unsign elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds the element to the vector.
 
    This function takes a pointer to a vector of proshade_signed and a single proshade_signed element and adds this element to the end of
    the vector. The reason for this function is to make vector elongation C++ standard independent (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of proshade_signed's which should be elongated.
    \param[in] elementToAdd proshade_signed to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToSignedVector ( std::vector < proshade_signed >* vecToAddTo, proshade_signed elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds the element to the vector.
 
    This function takes a pointer to a vector of complex<double> and a single complex<double> element and adds this element to the end of
    the vector. The reason for this function is to make vector elongation C++ standard independent (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of complex<double>'s which should be elongated.
    \param[in] elementToAdd complex<double> to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToDblPtrVector ( std::vector < proshade_double* >* vecToAddTo, proshade_double* elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds the element to the vector.
 
    This function takes a pointer to a vector of proshade_signed pointers and a single proshade_signed pointer and adds this
    element to the end of the vector. The reason for this function is to make vector elongation C++ standard independent
    (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of complex<double>'s which should be elongated.
    \param[in] elementToAdd complex<double> to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToSigPtrVector ( std::vector < proshade_signed* >* vecToAddTo, proshade_signed* elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds the element to the vector.
 
    This function takes a pointer to a vector of proshade_unsign pointers and a single proshade_unsign pointer and adds this
    element to the end of the vector. The reason for this function is to make vector elongation C++ standard independent
    (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of proshade_unsign pointer's which should be elongated.
    \param[in] elementToAdd proshade_unsign pointer to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToUnsPtrVector ( std::vector < proshade_unsign* >* vecToAddTo, proshade_unsign* elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds the element to the vector of vectors.
 
    This function takes a pointer to a vector of vectors of unsigns and a single vector of unsigns adds this element to the end of
    the vector of vectors. The reason for this function is to make vector elongation C++ standard independent (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of vectors of unsigns which should be elongated.
    \param[in] elementToAdd vector of unsigns to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToUnsignVectorVector ( std::vector < std::vector < proshade_unsign  > >* vecToAddTo, std::vector < proshade_unsign  > elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief Adds the element to the vector of vectors.
 
    This function takes a pointer to a vector of vectors of doubles and a single vector of doubles and adds this element to the end of
    the vector of vectors. The reason for this function is to make vector elongation C++ standard independent (push_back vs. emplace_back)
 
    \param[in] vecToAddTo Pointer to vector of vectors of doubles which should be elongated.
    \param[in] elementToAdd vector of doubles to be added to the back of the vector.
 */
void ProSHADE_internal_misc::addToDoubleVectorVector ( std::vector < std::vector < proshade_double  > >* vecToAddTo, std::vector < proshade_double  > elementToAdd )
{
    //================================================ Based on the compiler C++11 support, use the correct vector addition function.
#if __cplusplus >= 201103L
    vecToAddTo->emplace_back                          ( elementToAdd );
#else
    vecToAddTo->push_back                             ( elementToAdd );
#endif
    
    //================================================ Done
    return ;
    
}

/*! \brief This function compares two arrays of two based on the fifth number, sorting lowest first.
 
    \param[in] a The first array to compare.
    \param[in] b The second array to compare.
    \param[out] X Boolean whether the first is smaller than the second.
 */
bool ProSHADE_internal_misc::sortSymHlp ( const proshade_double* a, const proshade_double* b )
{
    //================================================ Compare
    return                                            ( a[5] < b[5] );
    
}

/*! \brief This function compares two arrays of two based on the fifth number, sorting highest first.
 
    \param[in] a The first array to compare.
    \param[in] b The second array to compare.
    \param[out] X Boolean whether the first is smaller than the second.
 */
bool ProSHADE_internal_misc::sortSymHlpInv ( const proshade_double* a, const proshade_double* b )
{
    //================================================ Compare
    return                                            ( a[5] > b[5] );
    
}

/*! \brief This function compares two vectors of arrays based on the sum of the fifth column, sorting highest first.
 
    \param[in] a The first vector to compare.
    \param[in] b The second vector to compare.
    \param[out] X Boolean whether the first is smaller than the second.
 */
bool ProSHADE_internal_misc::sortISymByPeak ( const std::vector< proshade_double* > a, const std::vector< proshade_double* > b )
{
    //================================================ If one is less or more than 31, it is wrong
    if ( a.size() != 31 ) { return ( false ); }
    if ( b.size() != 31 ) { return ( true  ); }
    
    //================================================ Find the sums
    proshade_double aSum = 0.0, bSum = 0.0;
    for ( size_t it = 0; it < 31; it++ )
    {
        aSum                                         += a.at(it)[5];
        bSum                                         += b.at(it)[5];
    }
    
    //================================================ Decide
    return                                            ( aSum > bSum );
    
} 

/*! \brief This function compares two arrays of the ProSHADE dihedral symmetry list based on combination of axes folds and heightst.
 
    \param[in] a The first array to compare.
    \param[in] b The second array to compare.
    \param[out] X Boolean whether the first is smaller than the second.
 */
bool ProSHADE_internal_misc::sortDSymHlpInv ( const proshade_double* a, const proshade_double* b )
{
    //================================================ Get weighted averages
    proshade_double aScore                            = ( ( a[0] * a[5] ) + ( a[6] * a[11] ) ) / ( a[0] + a[6] );
    proshade_double bScore                            = ( ( b[0] * b[5] ) + ( b[6] * b[11] ) ) / ( b[0] + b[6] );
    
    //================================================ Compare
    return                                            ( aScore > bScore );
    
}

/*! \brief This function compares two arrays of two based on the first number, sorting highest first.
 
    \param[in] a The first array to compare.
    \param[in] b The second array to compare.
    \param[out] X Boolean whether the first is smaller than the second.
 */
bool ProSHADE_internal_misc::sortSymInvFoldHlp ( const proshade_double* a, const proshade_double* b )
{
    //================================================ Compare
    return                                            ( a[0] > b[0] );
    
}

/*! \brief Does a deep copy of a double array to a vector of double arrays.
 
    This function deep copies a single symmetry axis ( proshade_double[7] ) into a vector of such arrays.
 
    \param[in] dblPtrVec Pointer to vector of proshade_double arrays to which the second argument is to be deep copied.
    \param[in] axis The proshade_double array to be copied to the first argument.
 */
void ProSHADE_internal_misc::deepCopyAxisToDblPtrVector ( std::vector< proshade_double* >* dblPtrVec, proshade_double* axis )
{
    //================================================ Allocate new memory
    proshade_double* symAx                            = new proshade_double[7];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( symAx, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy the 6 values
    symAx[0]                                          = axis[0];
    symAx[1]                                          = axis[1];
    symAx[2]                                          = axis[2];
    symAx[3]                                          = axis[3];
    symAx[4]                                          = axis[4];
    symAx[5]                                          = axis[5];
    symAx[6]                                          = axis[6];
    
    //================================================ Save
    ProSHADE_internal_misc::addToDblPtrVector         ( dblPtrVec, symAx );
    
    //================================================ Done
    return ;
    
}

/*! \brief Does a deep copy of a signed int array to a vector of signed int arrays.
 
    This function deep copies a single bounds array ( proshade_signed[6] ) into a vector of such arrays.
 
    \param[in] sigPtrVec Pointer to vector of proshade_signed arrays to which the second argument is to be deep copied.
    \param[in] xFrom Pointer to the index value from which the x dimension start.
    \param[in] xTo Pointer to the index value to which the x dimensions runs.
    \param[in] yFrom Pointer to the index value from which the y dimension start.
    \param[in] yTo Pointer to the index value to which the y dimensions runs.
    \param[in] zFrom Pointer to the index value from which the z dimension start.
    \param[in] zTo Pointer to the index value to which the z dimensions runs.
 */
void ProSHADE_internal_misc::deepCopyBoundsSigPtrVector ( std::vector < proshade_signed* >* sigPtrVec, proshade_signed* xFrom, proshade_signed* xTo, proshade_signed* yFrom, proshade_signed* yTo, proshade_signed* zFrom, proshade_signed* zTo  )
{
    //================================================ Allocate new memory
    proshade_signed* bounds                           = nullptr;
    bounds                                            = new proshade_signed[6];
    
    //================================================ Check memory allocation
    ProSHADE_internal_misc::checkMemoryAllocation     ( bounds, __FILE__, __LINE__, __func__ );
    
    //================================================ Copy the 6 values
    bounds[0]                                         = *xFrom;
    bounds[1]                                         = *xTo;
    bounds[2]                                         = *yFrom;
    bounds[3]                                         = *yTo;
    bounds[4]                                         = *zFrom;
    bounds[5]                                         = *zTo;
    
    //================================================ Save
    ProSHADE_internal_misc::addToSigPtrVector         ( sigPtrVec, bounds );
    
    //================================================ Done
    return ;
    
}
