/*! \file ProSHADE_misc.hpp
 \brief This header contains all miscellaneous functions.
 
 The functions declared in here are used by ProSHADE to deal with minor tasks unrelated to anything
 in the main library of functions.
 
 This file is part of the ProSHADE library for calculating
 shape descriptors and symmetry operators of protein structures.
 This is a prototype code, which is by no means complete or fully
 tested. Its use is at your own risk only. There is no quarantee
 that the results are correct.
 
 \author    Michal Tykac
 \author    Garib N. Murshudov
 \version   0.7.2
 \date      DEC 2019
 */

//============================================ ProSHADE
#include "ProSHADE_messages.hpp"

//============================================ Overinclusion protection
#ifndef __PROSHADE_MISC__
#define __PROSHADE_MISC__

//============================================ ProSHADE_internal_messages Namespace
/*! \namespace ProSHADE_internal_misc
 \brief This namespace contains all the miscellaneous functions used throughout the code but unrelated to it.
 
 The ProSHADE_internal_misc namespace contains all the functions which are useful, but not really related to shape
 description or data processing.
 */
namespace ProSHADE_internal_misc
{
    void addToStringVector                    ( std::vector < std::string >* vecToAddTo, std::string elementToAdd );
    void addToSingleVector                    ( std::vector < proshade_single  >* vecToAddTo, proshade_single elementToAdd );
    void addToDoubleVector                    ( std::vector < proshade_double  >* vecToAddTo, proshade_double elementToAdd );
    void addToUnsignVector                    ( std::vector < proshade_unsign  >* vecToAddTo, proshade_unsign elementToAdd );
    void addToSignedVector                    ( std::vector < proshade_signed  >* vecToAddTo, proshade_signed elementToAdd );
    void addToDblPtrVector                    ( std::vector < proshade_double* >* vecToAddTo, proshade_double* elementToAdd );
    void addToSigPtrVector                    ( std::vector < proshade_signed* >* vecToAddTo, proshade_signed* elementToAdd );
    void addToUnsignVectorVector              ( std::vector < std::vector < proshade_unsign  > >* vecToAddTo, std::vector < proshade_unsign  > elementToAdd );
    
    bool sortSymHlp                           ( const proshade_double* a, const proshade_double* b );
    bool sortSymHlpInv                        ( const proshade_double* a, const proshade_double* b );
    bool sortDSymHlpInv                       ( const proshade_double* a, const proshade_double* b );

    void deepCopyAxisToDblPtrVector           ( std::vector < proshade_double* >* dblPtrVec, proshade_double* axis );
    void deepCopyBoundsSigPtrVector           ( std::vector < proshade_signed* >* sigPtrVec, proshade_signed* xFrom, proshade_signed* xTo, proshade_signed* yFrom,
                                                proshade_signed* yTo, proshade_signed* zFrom, proshade_signed* zTo );
    
    /*! \brief Checks if memory was allocated properly.
     
     This function checks if the memory allocation has suceeded for a given pointer, printing error message if not.
     
     \param[in] checkVar Pointer to be checked.
     */
    
    template <class chVar> inline void checkMemoryAllocation ( chVar checkVar, std::string fileP, unsigned int lineP, std::string funcP, std::string infoP = "This error may occurs when ProSHADE requests memory to be\n                    : allocated to it and this operation fails. This could\n                    : happen when not enough memory is available, either due to\n                    : other processes using a lot of memory, or when the machine\n                    : does not have sufficient memory available. Re-run to see\n                    : if this problem persists." )
    {
        //======================================== Check against NULL
        if ( checkVar == NULL )
        {
            throw ProSHADE_exception ( "Failed to allocate memory.", "E000007", fileP, lineP, funcP, infoP );
        }
        
        //======================================== Done
        return ;
        
    }
}

#endif
