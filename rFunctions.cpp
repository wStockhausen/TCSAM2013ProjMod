// rFunctions.cpp
// This file has C++ functions for use in writing R-compatible data output.
// WTS 2005-04: adapted from mhp-s-funcs.cpp by Michael Prager  -- December, 2002
#include <admodel.h>
#include "ModelConstants.hpp"
#include "admbFunctions.hpp"
#include "rFunctions.hpp"
namespace Rpr{
    //===========================================================================================
    // ADMB FUNCTION to write a matrix as part of an array structure
    // ----------Arguments----------
    // os       stream for output file.
    // xx      matrix of data to be written.
    // ----------Value--------------
    // dim    dims for writing the matrix to an array structure
    adstring writeDataToR(ostream& os, const dmatrix& xx){
        ivector bds = wts::getBounds(xx);
        int ctr = 1;
        for (int j=bds(3);j<bds(4);j++) {
            for (int i=bds(1);i<=bds(2);i++)  {
                os<<xx(i,j)<<cc;
                if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
            }
        }
        for (int i=bds(1);i<bds(2);i++)  {
            os<<xx(i,bds(4))<<cc;
            if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
        }
        os<<xx(bds(2),bds(4));
        adstring dim = str(bds(2)-bds(1)+1)+cc+
                       str(bds(4)-bds(3)+1);
        return dim;
    }
    //===========================================================================================
    // ADMB FUNCTION to write a matrix as part of an R list
    // ----------Arguments----------
    // os       stream for output file.
    // xx       matrix of data to be written.
    // dimnames adstring with dimnames.
    void writeToR(ostream& os, const dmatrix& xx, adstring dimnames){
        os<<"structure(c(";
        adstring dims = Rpr::writeDataToR(os,xx);
        os<<"),"<<endl<<tb<<tb;    
        os<<"dimnames=list("<<dimnames<<"),";
        os<<"dim=c("<<dims<<"))";
    }
    
    //===========================================================================================
    // ADMB FUNCTION to write a d3_array as part of an array structure
    // ----------Arguments----------
    // os   stream for output file.
    // xx   d3_array to be written.
    // ----------Value--------------
    // dim  dims for writing the matrix to an array structure
    adstring writeDataToR(ostream& os, const d3_array& xx){
        ivector bds = wts::getBounds(xx);
        int ctr = 1;
        for (int k=bds(5);k<bds(6);k++) {
            for (int j=bds(3);j<=bds(4);j++) {
                for (int i=bds(1);i<=bds(2);i++)  {
                    os<<xx(i,j,k)<<cc;
                    if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
                }
            }
        }
        for (int j=bds(3);j<bds(4);j++) {
            for (int i=bds(1);i<=bds(2);i++)  {
                os<<xx(i,j,bds(6))<<cc;
                if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
            }
        }
        for (int i=bds(1);i<bds(2);i++)  {
            os<<xx(i,bds(4),bds(6))<<cc;
            if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
        }            
        os<<xx(bds(2),bds(4),bds(6));
        adstring dim = str(bds(2)-bds(1)+1)+cc+
                       str(bds(4)-bds(3)+1)+cc+
                       str(bds(6)-bds(5)+1);
        return dim;
    }
    //===========================================================================================
    // ADMB FUNCTION to write a d3_array as part of an R list
    // ----------Arguments----------
    // os       stream for output file.
    // xx       d3_array to be written.
    // dimnames adstring with dimnames.
    void writeToR(ostream& os, const d3_array& xx, adstring dimnames){
        os<<"structure(c(";
        adstring dims = Rpr::writeDataToR(os,xx);
        os<<"),"<<endl<<tb<<tb;    
        os<<"dimnames=list("<<dimnames<<"),";
        os<<"dim=c("<<dims<<"))";
    }
    
    //===========================================================================================
    // ADMB FUNCTION to write a d4_array as part of an array structure
    // ----------Arguments----------
    // os   stream for output file.
    // xx   d4_array to be written.
    // ----------Value--------------
    // dim    dims for writing the d4_array to an array structure
    adstring writeDataToR(ostream& os, const d4_array& xx){
        ivector bds = wts::getBounds(xx);
        int ctr = 1;
        for (int l=bds(7);l<bds(8);l++) {
            for (int k=bds(5);k<=bds(6);k++) {
                for (int j=bds(3);j<=bds(4);j++) {
                    for (int i=bds(1);i<=bds(2);i++)  {
                        os<<xx(i,j,k,l)<<cc;
                        if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
                    }
                }
            }
        }
        for (int k=bds(5);k<bds(6);k++) {
            for (int j=bds(3);j<=bds(4);j++) {
                for (int i=bds(1);i<=bds(2);i++)  {
                    os<<xx(i,j,k,bds(8))<<cc;
                    if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
                }
            }
        }
        for (int j=bds(3);j<bds(4);j++) {
            for (int i=bds(1);i<=bds(2);i++)  {
                os<<xx(i,j,bds(6),bds(8))<<cc;
                if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
            }
        }
        for (int i=bds(1);i<bds(2);i++)  {
            os<<xx(i,bds(4),bds(6),bds(8))<<cc;
            if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
        }            
        os<<xx(bds(2),bds(4),bds(6),bds(8));
        adstring dim = str(bds(2)-bds(1)+1)+cc+
                       str(bds(4)-bds(3)+1)+cc+
                       str(bds(6)-bds(5)+1)+cc+
                       str(bds(8)-bds(7)+1);
        return dim;
    }
    //===========================================================================================
    // ADMB FUNCTION to write a d4_array as part of an R list
    // ----------Arguments----------
    // os       stream for output file.
    // xx       d4_array to be written.
    // dimnames adstring with dimnames.
    void writeToR(ostream& os, const d4_array& xx, adstring dimnames){
        os<<"structure(c(";
        adstring dims = Rpr::writeDataToR(os,xx);
        os<<"),"<<endl<<tb<<tb;    
        os<<"dimnames=list("<<dimnames<<"),";
        os<<"dim=c("<<dims<<"))";
    }
    
    //===========================================================================================
    // ADMB FUNCTION to write a d5_array as part of an array structure
    // ----------Arguments----------
    // os   stream for output file.
    // xx   d5_array to be written.
    // ----------Value--------------
    // dim    dims for writing the d4_array to an array structure
    adstring writeDataToR(ostream& os, const d5_array& xx){
        ivector bds = wts::getBounds(xx);
        int ctr = 1;
        for (int m=bds(9);m<bds(10);m++) {
            for (int l=bds(7);l<=bds(8);l++) {
                for (int k=bds(5);k<=bds(6);k++) {
                    for (int j=bds(3);j<=bds(4);j++) {
                        for (int i=bds(1);i<=bds(2);i++)  {
                            os<<xx(i,j,k,l,m)<<cc;
                            if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
                        }
                    }
                }
            }
        }
        for (int l=bds(7);l<bds(8);l++) {
            for (int k=bds(5);k<=bds(6);k++) {
                for (int j=bds(3);j<=bds(4);j++) {
                    for (int i=bds(1);i<=bds(2);i++)  {
                        os<<xx(i,j,k,l,bds(10))<<cc;
                        if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
                    }
                }
            }
        }
        for (int k=bds(5);k<bds(6);k++) {
            for (int j=bds(3);j<=bds(4);j++) {
                for (int i=bds(1);i<=bds(2);i++)  {
                    os<<xx(i,j,k,bds(8),bds(10))<<cc;
                    if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
                }
            }
        }
        for (int j=bds(3);j<bds(4);j++) {
            for (int i=bds(1);i<=bds(2);i++)  {
                os<<xx(i,j,bds(6),bds(8),bds(10))<<cc;
                if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
            }
        }
        for (int i=bds(1);i<bds(2);i++)  {
            os<<xx(i,bds(4),bds(6),bds(8),bds(10))<<cc;
            if (++ctr>100){os<<endl<<tb<<tb; ctr=0;}
        }            
        os<<xx(bds(2),bds(4),bds(6),bds(8),bds(10));
        adstring dim = str(bds(2)-bds(1)+1)+cc+
                       str(bds(4)-bds(3)+1)+cc+
                       str(bds(6)-bds(5)+1)+cc+
                       str(bds(8)-bds(7)+1)+cc+
                       str(bds(10)-bds(9)+1);
        return dim;
    }
    //===========================================================================================
    // ADMB FUNCTION to write a d5_array as part of an R list
    // ----------Arguments----------
    // os       stream for output file.
    // xx       d5_array to be written.
    // dimnames adstring with dimnames.
    void writeToR(ostream& os, const d5_array& xx, adstring dimnames){
        os<<"structure(c(";
        adstring dims = Rpr::writeDataToR(os,xx);
        os<<"),"<<endl<<tb<<tb;    
        os<<"dimnames=list("<<dimnames<<"),";
        os<<"dim=c("<<dims<<"))";
    }
}//namespace Rpr

namespace R{
    //===========================================================================================
        // ADMB FUNCTION to write an ivector as part of an R list
        //
        // ----------Arguments----------
        // os       stream for output file.
        // xx       ivector to be written.
    void writeToR(ostream& os, const ivector& xx){
        int mn = xx.indexmin();
        int mx = xx.indexmax();
        os<<"structure(c(";
        for (int i=mn;i<mx;i++) os<<xx(i)<<cc;  os<<xx(mx)<<"),"<<endl<<tb<<tb;
        os<<"names="<<mn<<":"<<mx<<",dim=c("<<mx-mn+1<<"))";
    }
    //===========================================================================================
        // ADMB FUNCTION to write a vector as part of an R list
        //
        // ----------Arguments----------
        // os       stream for output file.
        // xx       vector of data to be written.
    void writeToR(ostream& os, const dvector& xx){
        int mn = xx.indexmin();
        int mx = xx.indexmax();
        os<<"structure(c(";
        for (int i=mn;i<mx;i++) os<<xx(i)<<cc;  os<<xx(mx)<<"),"<<endl<<tb<<tb;
        os<<"names="<<mn<<":"<<mx<<",dim=c("<<mx-mn+1<<"))";
    }
    //===========================================================================================
        // ADMB FUNCTION to write a vector as part of an R list
        //
        // ----------Arguments----------
        // os       stream for output file.
        // xx       vector of data to be written.
        // names    adstring to be written as names for R structure (comma-delimited, quoted if necessary)
    void writeToR(ostream& os, const dvector& xx, adstring names){
        int mn = xx.indexmin();
        int mx = xx.indexmax();
        os<<"structure(c(";
        for (int i=mn;i<mx;i++) os<<xx(i)<<cc;  os<<xx(mx)<<"),"<<endl<<tb<<tb;
        os<<"names=c("<<names<<"),dim=c("<<mx-mn+1<<"))";
    }
    //===========================================================================================
        // ADMB FUNCTION to write a vector as part of an R list
        //
        // ----------Arguments----------
        // os       stream for output file.
        // xx       vector of data to be written.
        // names    adstring_array to be written as names for R structure 
    void writeToR(ostream& os, const dvector& xx, adstring_array names){
        int mn = xx.indexmin();
        int mx = xx.indexmax();
        os<<"structure(c(";
        for (int i=mn;i<mx;i++) os<<xx(i)<<cc;  os<<xx(mx)<<"),"<<endl<<tb<<tb;
        os<<"names=c('"; for (int i=mn;i<mx;i++) os<<names(i)<<"','"; os<<names(mx)<<"'))";
    }
    //===========================================================================================
        // ADMB FUNCTION to write a matrix as part of an R list
        // ----------Arguments----------
        // os       stream for output file.
        // xx      matrix of data to be written.
    void writeToR(ostream& os, const dmatrix& xx){
        int mnI = xx.indexmin();
        int mxI = xx.indexmax();
        int mnJ = xx(mnI).indexmin();
        int mxJ = xx(mnJ).indexmax();
        adstring dimnames=str(mnI)+":"+str(mxI)+cc+str(mnJ)+":"+str(mxJ);
        Rpr::writeToR(os,xx,dimnames);
    }
    //===========================================================================================
    // ADMB FUNCTION to write a matrix as part of an R list
    // ----------Arguments----------
    // os       stream for output file.
    // xx      matrix of data to be written.
    // colnames comma-delimited string of column names
    void writeToR(ostream& os, const dmatrix& xx, adstring colnames){
        int mnI = xx.indexmin();
        int mxI = xx.indexmax();
        adstring dimnames=str(mnI)+":"+str(mxI)+cc+"c("+colnames+")";
        Rpr::writeToR(os,xx,dimnames);
    }
    //===========================================================================================
    // ADMB FUNCTION to write a matrix as part of an R list
    // ----------Arguments----------
    // os   stream for output file.
    // xx   matrix of data to be written.
    // n1   comma-delimited, single-quoted string of names for 1st (leftmost) index
    // n2   comma-delimited, single-quoted string of names for 2nd index
    void writeToR(ostream& os, const dmatrix& xx, adstring n1, adstring n2){
        adstring dimnames="c("+n1+")"+cc+"c("+n2+")";
        Rpr::writeToR(os,xx,dimnames);
    }
    //===========================================================================================
    // ADMB FUNCTION to write a matrix as part of an R list
    // ----------Arguments----------
    // os  stream for output file.
    // xx  matrix of data to be written.
    // n1  ivector to use as names for 1st (leftmost; row) index
    // n2  comma-delimited, single-quoted string of names for 2nd (column) index
    void writeToR(ostream& os, const dmatrix& xx, const ivector& n1, adstring n2){
        adstring rows = wts::to_qcsv(n1);
        R::writeToR(os,xx,rows,n2);
    }
    //===========================================================================================
    // ADMB FUNCTION to write a d3_array as part of an R list
    // ----------Arguments----------
    // os  stream for output file.
    // xx  data to be written.
    // n1  comma-delimited, single-quoted string of names for 1st (leftmost) index
    // n2  comma-delimited, single-quoted string of names for 2nd index
    // n3  comma-delimited, single-quoted string of names for 3rd index
    void writeToR(ostream& os, const d3_array& xx, adstring n1, adstring n2, adstring n3){
        adstring dimnames="c("+n1+")"+cc+"c("+n2+")"+cc+"c("+n3+")";
        Rpr::writeToR(os,xx,dimnames);
    }
    //===========================================================================================
    // ADMB FUNCTION to write a d4_array as part of an R list
    // ----------Arguments----------
    // os  stream for output file.
    // xx  data to be written.
    // n1  comma-delimited, single-quoted string of names for 1st (leftmost) index
    // n2  comma-delimited, single-quoted string of names for 2nd index
    // n3  comma-delimited, single-quoted string of names for 3rd index
    // n4  comma-delimited, single-quoted string of names for 4th index
    void writeToR(ostream& os, const d4_array& xx, adstring n1, adstring n2, adstring n3, adstring n4){
        adstring dimnames="c("+n1+")"+cc+"c("+n2+")"+cc+"c("+n3+")"+cc+"c("+n4+")";
        Rpr::writeToR(os,xx,dimnames);
    }
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
    void writeToR(ostream& os, const d5_array& xx, adstring n1, adstring n2, adstring n3, adstring n4, adstring n5){
        adstring dimnames="c("+n1+")"+cc+"c("+n2+")"+cc+"c("+n3+")"+cc+"c("+n4+")"+cc+"c("+n5+")";
        Rpr::writeToR(os,xx,dimnames);
    }
}//namespace R