// rFunctions.hpp
// This file has C++ functions for use in writing R-compatible data output.
// WTS 2005-04: adapted from mhp-s-funcs.cpp by Michael Prager  -- December, 2002
#pragma once
#ifndef __HPP_RFUNCTIONS__
    #define __HPP_RFUNCTIONS__

namespace Rpr{
    //===========================================================================================
        // ADMB FUNCTION to write a matrix as part of an array structure
        // ----------Arguments----------
        // os       stream for output file.
        // xx      matrix of data to be written.
        // ----------Value--------------
        // dim    dims for writing the matrix to an array structure
    adstring writeDataToR(ostream& os, const dmatrix& xx);
    //===========================================================================================
        // ADMB FUNCTION to write a matrix as part of an R list
        // ----------Arguments----------
        // os       stream for output file.
        // xx       matrix of data to be written.
        // dimnames adstring with dimnames.
    void writeToR(ostream& os, const dmatrix& xx, adstring dimnames);
}//namespace Rpr

namespace R{
    //===========================================================================================
        // ADMB FUNCTION to write an ivector as part of an R list
        //
        // ----------Arguments----------
        // os       stream for output file.
        // xx       ivector to be written.
    void writeToR(ostream& os, const ivector& xx);
    //===========================================================================================
        // ADMB FUNCTION to write a vector as part of an R list
        //
        // ----------Arguments----------
        // os       stream for output file.
        // xx       vector of data to be written.
    void writeToR(ostream& os, const dvector& xx);
    //===========================================================================================
        // ADMB FUNCTION to write a vector as part of an R list
        //
        // ----------Arguments----------
        // os       stream for output file.
        // xx       vector of data to be written.
        // names    adstring to be written as names for R structure (comma-delimited, quoted if necessary)
    void writeToR(ostream& os, const dvector& xx, adstring names);
    //===========================================================================================
        // ADMB FUNCTION to write a vector as part of an R list
        //
        // ----------Arguments----------
        // os       stream for output file.
        // xx       vector of data to be written.
        // names    adstring_array to be written as names for R structure 
    void writeToR(ostream& os, const dvector& xx, adstring_array names);
    //===========================================================================================
        // ADMB FUNCTION to write a matrix as part of an R list
        // ----------Arguments----------
        // os       stream for output file.
        // xx      matrix of data to be written.
    void writeToR(ostream& os, const dmatrix& xx);
    //===========================================================================================
    // ADMB FUNCTION to write a matrix as part of an R list
    // ----------Arguments----------
    // os       stream for output file.
    // xx      matrix of data to be written.
    // colnames comma-delimited, single-quoted string of column names
    void writeToR(ostream& os, const dmatrix& xx, adstring colnames);
    //===========================================================================================
    // ADMB FUNCTION to write a matrix as part of an R list
    // ----------Arguments----------
    // os   stream for output file.
    // xx   matrix of data to be written.
    // n1   comma-delimited, single-quoted string of names for 1st (leftmost) index
    // n2   comma-delimited, single-quoted string of names for 2nd index
    void writeToR(ostream& os, const dmatrix& xx, adstring n1, adstring n2);
    //===========================================================================================
    // ADMB FUNCTION to write a matrix as part of an R list
    // ----------Arguments----------
    // os  stream for output file.
    // xx  matrix of data to be written.
    // n1  ivector to use as names for 1st (leftmost; row) index
    // n2  comma-delimited, single-quoted string of names for 2nd (column) index
    void writeToR(ostream& os, const dmatrix& xx, const ivector& n1, adstring n2);
    //===========================================================================================
    // ADMB FUNCTION to write a d3_array as part of an R list
    // ----------Arguments----------
    // os  stream for output file.
    // xx  data to be written.
    // n1  comma-delimited, single-quoted string of names for 1st (leftmost) index
    // n2  comma-delimited, single-quoted string of names for 2nd index
    // n3  comma-delimited, single-quoted string of names for 3rd index
    void writeToR(ostream& os, const d3_array& xx, adstring n1, adstring n2, adstring n3);
    //===========================================================================================
    // ADMB FUNCTION to write a d4_array as part of an R list
    // ----------Arguments----------
    // os  stream for output file.
    // xx  data to be written.
    // n1  comma-delimited, single-quoted string of names for 1st (leftmost) index
    // n2  comma-delimited, single-quoted string of names for 2nd index
    // n3  comma-delimited, single-quoted string of names for 3rd index
    // n4  comma-delimited, single-quoted string of names for 4th index
    void writeToR(ostream& os, const d4_array& xx, adstring n1, adstring n2, adstring n3, adstring n4);
    //===========================================================================================
    // ADMB FUNCTION to write a d5_array as part of an R list
    // ----------Arguments----------
    // os  stream for output file.
    // xx  data to be written.
    // n1  comma-delimited, single-quoted string of names for 1st (leftmost) index
    // n2  comma-delimited, single-quoted string of names for 2nd index
    // n3  comma-delimited, single-quoted string of names for 3rd index
    // n4  comma-delimited, single-quoted string of names for 4th index
    // n5  comma-delimited, single-quoted string of names for 5th index
    void writeToR(ostream& os, const d5_array& xx, adstring n1, adstring n2, adstring n3, adstring n4, adstring n5);
} //namespace R
#endif
