/*! \file ProSHADE_misc.hpp
    \brief This header file declares all miscellaneous functions.
 
    The functions declared in here are used by ProSHADE to deal with minor tasks unrelated to anything
    in the main library of functions, such as memory allocation checking and C++ version independent vector
    extension.
 
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

//==================================================== Overinclusion protection
#ifndef PROSHADE_MISC
#define PROSHADE_MISC

//==================================================== ProSHADE_internal_messages Namespace
/*! \namespace ProSHADE_internal_misc
    \brief This namespace contains all the miscellaneous functions used throughout the code but unrelated to it.
 
    The ProSHADE_internal_misc namespace contains all the functions which are useful, but not really related to shape
    description or data processing.
 */
namespace ProSHADE_internal_misc
{
    void addToStringVector                            ( std::vector < std::string >* vecToAddTo, std::string elementToAdd );
    void addToSingleVector                            ( std::vector < proshade_single  >* vecToAddTo, proshade_single elementToAdd );
    void addToDoubleVector                            ( std::vector < proshade_double  >* vecToAddTo, proshade_double elementToAdd );
    void addToUnsignVector                            ( std::vector < proshade_unsign  >* vecToAddTo, proshade_unsign elementToAdd );
    void addToSignedVector                            ( std::vector < proshade_signed  >* vecToAddTo, proshade_signed elementToAdd );
    void addToDblPtrVector                            ( std::vector < proshade_double* >* vecToAddTo, proshade_double* elementToAdd );
    void addToSigPtrVector                            ( std::vector < proshade_signed* >* vecToAddTo, proshade_signed* elementToAdd );
    void addToUnsPtrVector                            ( std::vector < proshade_unsign* >* vecToAddTo, proshade_unsign* elementToAdd );
    void addToUnsignVectorVector                      ( std::vector < std::vector < proshade_unsign  > >* vecToAddTo, std::vector < proshade_unsign  > elementToAdd );
    void addToDoubleVectorVector                      ( std::vector < std::vector < proshade_double  > >* vecToAddTo, std::vector < proshade_double  > elementToAdd );
            
    bool sortSymHlp                                   ( const proshade_double* a, const proshade_double* b );
    bool sortSymHlpInv                                ( const proshade_double* a, const proshade_double* b );
    bool sortISymByPeak                               ( const std::vector< proshade_double* > a, const std::vector< proshade_double* > b );
    bool sortDSymHlpInv                               ( const proshade_double* a, const proshade_double* b );
    bool sortSymFoldHlp                               ( const proshade_double* a, const proshade_double* b );
    bool sortSymInvFoldHlp                            ( const proshade_double* a, const proshade_double* b );
        
    void deepCopyAxisToDblPtrVector                   ( std::vector < proshade_double* >* dblPtrVec, proshade_double* axis );
    void deepCopyBoundsSigPtrVector                   ( std::vector < proshade_signed* >* sigPtrVec, proshade_signed* xFrom, proshade_signed* xTo, proshade_signed* yFrom,
                                                        proshade_signed* yTo, proshade_signed* zFrom, proshade_signed* zTo );
    
/*! \brief Checks if memory was allocated properly.

    This function checks if the memory allocation has suceeded for a given pointer, printing error message if not.

    \param[in] checkVar Pointer to be checked.
 */
    
    template <class chVar> inline void checkMemoryAllocation ( chVar checkVar, std::string fileP, unsigned int lineP, std::string funcP, std::string infoP = "This error may occurs when ProSHADE requests memory to be\n                    : allocated to it and this operation fails. This could\n                    : happen when not enough memory is available, either due to\n                    : other processes using a lot of memory, or when the machine\n                    : does not have sufficient memory available. Re-run to see\n                    : if this problem persists." )
    {
        //============================================ Check against NULL
        if ( checkVar == nullptr )
        {
            throw ProSHADE_exception                  ( "Failed to allocate memory.", "E000007", fileP, lineP, funcP, infoP );
        }
        
        //============================================ Done
        return ;
        
    }
}

#endif
