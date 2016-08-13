//admbFunctions.cpp
#include <admodel.h>
#include "ModelConstants.hpp"
#include "admbFunctions.hpp"

int debugADMBFunctions= 0;
namespace wts{
    /****************************************************************
     * convert vectors to quoted csv string
     ***************************************************************/
    adstring to_qcsv(const ivector& v){
        int mn = v.indexmin();
        int mx = v.indexmax();
        adstring s = qt+str(v(mn))+qt;
        for (int i=mn;i<mx;i++) s = s+cc+qt+str(v(i+1))+qt;
        return s;
    }
    adstring to_qcsv(const dvector& v){
        int mn = v.indexmin();
        int mx = v.indexmax();
        adstring s = qt+str(v(mn))+qt;
        for (int i=mn;i<mx;i++) s = s+cc+qt+str(v(i+1))+qt;
        return s;
    }
    /****************************************************************
     * convert dvector to ivector.
     ***************************************************************/
    ivector to_ivector(const dvector& v){
        ivector iv(v.indexmin(),v.indexmax());
        for (int i=v.indexmin();i<=v.indexmax();i++) iv(i) = floor(v(i));
        return iv;
    }
    /****************************************************************
    * name      : copy                                              *
    * purpose   : create a deep copy of a variable                  *
    *NOTE: Does not work for ragged arrays                          *
    ****************************************************************/
    dvector  copy(dvector& v){
        ivector bnds = getBounds(v);
        dvector c(bnds(1),bnds(2));
        c = v;
        return c;
    }
    //-------------------------------------------------------------
    dmatrix  copy(dmatrix& v){
        ivector bnds = getBounds(v);
        dmatrix c(bnds(1),bnds(2),bnds(3),bnds(4));
        c = v;
        return c;
    }
    //-------------------------------------------------------------
    d3_array copy(d3_array& v){
        ivector bnds = getBounds(v);
        d3_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6));
        c = v;
        return c;
    }
    //-------------------------------------------------------------
    d4_array copy(d4_array& v){
        ivector bnds = getBounds(v);
        d4_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8));
        c = v;
        return c;
    }
    //-------------------------------------------------------------
    d5_array copy(d5_array& v){
        ivector bnds = getBounds(v);
        d5_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8),bnds(9),bnds(10));
        c = v;
        return c;
    }
    //-------------------------------------------------------------
    d6_array copy(d6_array& v){
        ivector bnds = getBounds(v);
        d6_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8),bnds(9),bnds(10),bnds(11),bnds(12));
        c = v;
        return c;
    }
    //-------------------------------------------------------------
    d7_array copy(d7_array& v){
        ivector bnds = getBounds(v);
        d7_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8),bnds(9),bnds(10),bnds(11),bnds(12),bnds(13),bnds(14));
        c = v;
        return c;
    }
    //-------------------------------------------------------------
    dvar_vector  copy(dvar_vector& v){
        RETURN_ARRAYS_INCREMENT();
        ivector bnds = getBounds(v);
        dvar_vector c(bnds(1),bnds(2));
        c = v;
        RETURN_ARRAYS_DECREMENT();
        return c;
    }
    //-------------------------------------------------------------
    dvar_matrix  copy(dvar_matrix& v){
        RETURN_ARRAYS_INCREMENT();
        ivector bnds = getBounds(v);
        dvar_matrix c(bnds(1),bnds(2),bnds(3),bnds(4));
        c = v;
        RETURN_ARRAYS_DECREMENT();
        return c;
    }
    //-------------------------------------------------------------
    dvar3_array copy(dvar3_array& v){
        RETURN_ARRAYS_INCREMENT();
        ivector bnds = getBounds(v);
        dvar3_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6));
        c = v;
        RETURN_ARRAYS_DECREMENT();
        return c;
    }
    //-------------------------------------------------------------
    dvar4_array copy(dvar4_array& v){
        RETURN_ARRAYS_INCREMENT();
        ivector bnds = getBounds(v);
        dvar4_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8));
        c = v;
        RETURN_ARRAYS_DECREMENT();
        return c;
    }
    //-------------------------------------------------------------
    dvar5_array copy(dvar5_array& v){
        RETURN_ARRAYS_INCREMENT();
        ivector bnds = getBounds(v);
        dvar5_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8),bnds(9),bnds(10));
        c = v;
        RETURN_ARRAYS_DECREMENT();
        return c;
    }
    //-------------------------------------------------------------
    dvar6_array copy(dvar6_array& v){
        RETURN_ARRAYS_INCREMENT();
        ivector bnds = getBounds(v);
        dvar6_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8),bnds(9),bnds(10),bnds(11),bnds(12));
        c = v;
        RETURN_ARRAYS_DECREMENT();
        return c;
    }
    //-------------------------------------------------------------
    dvar7_array copy(dvar7_array& v){
        RETURN_ARRAYS_INCREMENT();
        ivector bnds = getBounds(v);
        dvar7_array c(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8),bnds(9),bnds(10),bnds(11),bnds(12),bnds(13),bnds(14));
        c = v;
        RETURN_ARRAYS_DECREMENT();
        return c;
    }

    /****************************************************************
    * name      : adstring_matrix                                   *
    * purpose   : functionality for (possibly ragged) matrix of     *
    *               of adstrings                                    *
    ****************************************************************/
    bool adstring_matrix::debug = false;
    /****************************************************************
    * Constructor                                                   *
    ****************************************************************/
    adstring_matrix::adstring_matrix(){
        nAAs = 0;
        idxmn = 0;
        ppAAs = 0;
    }

    /****************************************************************
    * Constructor                                                   *
    ****************************************************************/
    adstring_matrix::adstring_matrix(int rwmn, int rwmx, int clmn, int clmx){
        allocate(rwmn,rwmx,clmn,clmx);
    }

    /****************************************************************
    * Constructor                                                   *
    ****************************************************************/
    adstring_matrix::adstring_matrix(int rwmn, int rwmx, ivector& clmns, int clmx){
        allocate(rwmn,rwmx,clmns,clmx);
    }

    /****************************************************************
    * Constructor                                                   *
    ****************************************************************/
    adstring_matrix::adstring_matrix(int rwmn, int rwmx, int clmn, ivector& clmxs){
        allocate(rwmn,rwmx,clmn,clmxs);
    }

    /****************************************************************
    * Constructor                                                   *
    ****************************************************************/
    adstring_matrix::adstring_matrix(int rwmn, int rwmx, ivector& clmns, ivector& clmxs){
        allocate(rwmn,rwmx,clmns,clmxs);
    }

    /****************************************************************
    * Constructor                                                   *
    ****************************************************************/
    void adstring_matrix::deallocate(void){
        if (ppAAs) {
            for (int r=0;r<nAAs;r++){
                delete ppAAs[r]; ppAAs[r]=0;
            }
            delete[] ppAAs; ppAAs = 0;
        }
        nAAs  = 0;
        idxmn = 0;
    }

    void adstring_matrix::allocate(int rwmn, int rwmx){
        deallocate();
        nAAs = rwmx-rwmn+1;
        if (nAAs<0) {
            nAAs = 0;
        } else {
            idxmn = rwmn;
            ppAAs = new adstring_array*[nAAs];
            for (int r=0;r<nAAs;r++) ppAAs[r] = 0;
        }
        if (debug){
            cout<<"adstring_matrix::allocate("<<rwmn<<","<<rwmx<<")"<<endl;
            cout<<"nAAs = "<<nAAs<<", ppAAs = "<<ppAAs<<endl;
        }
    }

    void adstring_matrix::allocate(int rwmn, int rwmx, int clmn, int clmx){
        deallocate();
        ivector clmns(rwmn,rwmx); clmns = clmn;
        ivector clmxs(rwmn,rwmx); clmxs = clmx;
        allocate(rwmn,rwmx,clmns,clmxs);
    }

    void adstring_matrix::allocate(int rwmn, int rwmx, ivector& clmns, int clmx){
        deallocate();
        ivector clmxs(rwmn,rwmx); clmxs = clmx;
        allocate(rwmn,rwmx,clmns,clmxs);
    }

    void adstring_matrix::allocate(int rwmn, int rwmx, int clmn, ivector& clmxs){
        deallocate();
        ivector clmns(rwmn,rwmx); clmns = clmn;
        allocate(rwmn,rwmx,clmns,clmxs);
    }

    void adstring_matrix::allocate(int rwmn, int rwmx, ivector& clmns, ivector& clmxs){
        allocate(rwmn,rwmx);
        if (rwmn!=clmns.indexmin()){
            cout<<"Error in adstring_matrix::allocate(int rwmn, int rwmx, ivector& clmns, ivector& clmxs)"<<endl;
            cout<<"min row index ("<<rwmn<<"must equal min min-columns index ("<<clmns.indexmin()<<")"<<endl;
            exit(0);
        }
        if (rwmn!=clmxs.indexmin()){
            cout<<"Error in adstring_matrix::allocate(int rwmn, int rwmx, ivector& clmns, ivector& clmxs)"<<endl;
            cout<<"min row index ("<<rwmn<<"must equal min max-columns index ("<<clmxs.indexmin()<<")"<<endl;
            exit(0);
        }
        if (rwmx!=clmns.indexmax()){
            cout<<"Error in adstring_matrix::allocate(int rwmn, int rwmx, ivector& clmns, ivector& clmxs)"<<endl;
            cout<<"max row index ("<<rwmn<<"must equal max min-columns index ("<<clmns.indexmax()<<")"<<endl;
            exit(0);
        }
        if (rwmx!=clmxs.indexmax()){
            cout<<"Error in adstring_matrix::allocate(int rwmn, int rwmx, ivector& clmns, ivector& clmxs)"<<endl;
            cout<<"max row index ("<<rwmx<<"must equal max max-columns index ("<<clmxs.indexmax()<<")"<<endl;
            exit(0);
        }
        this->clmns.allocate(rwmn,rwmx);
        this->clmxs.allocate(rwmn,rwmx);
        this->clmns = clmns;
        this->clmxs = clmxs;
        //note that the following means ppAAs[r] COULD BE 0 for some r even though adstring_arrays have been allocated
        //for other vlaues of r.  Thus need to check that ppAAs[r] is not a null ptr before accessing it.
        for (int r=0;r<nAAs;r++) if (clmns(r+rwmn)<=clmxs(r+rwmn)) ppAAs[r] = new adstring_array(clmns(r+rwmn),clmxs(r+rwmn));
        if (debug){
            cout<<"adstring_matrix::allocate(int rwmn, int rwmx, ivector& clmns, ivector& clmxs)"<<endl;
            cout<<"rwmn, rwmx = "<<rwmn<<", "<<rwmx<<endl;
            cout<<"clmns = "<<clmns<<",  clmxs = "<<clmxs<<endl;
            cout<<"ppAAs = "; for (int r=0;r<nAAs;r++) cout<<ppAAs[r]<<tb; cout<<endl;
        }
    }

    adstring& adstring_matrix::operator() (const int i, const int j){
        if (ppAAs){
            if ((idxmn<=i)&&(i-idxmn<nAAs)&&(ppAAs[i-idxmn])){
                if ((clmns(i)<=j)&&(j<=clmxs(i))){
                    return (*(ppAAs[i-idxmn]))(j);
                } else {
                    cout<<"Attempting illegal access in adstring_matrix::("<<i<<","<<j<<")"<<endl;
                    cout<<"Column index is invalid! Limits are "<<clmns(i)<<","<<clmxs(i)<<"."<<endl;
                }
            } else {
                cout<<"Attempting illegal access in adstring_matrix::("<<i<<","<<j<<")"<<endl;
                cout<<"Row index is invalid! Limits are "<<idxmn<<","<<idxmn+nAAs<<"."<<endl;
            }
        } else {
            cout<<"Attempting illegal access in adstring_matrix::("<<i<<","<<j<<")"<<endl;
            cout<<"Matrix is not allocated! nAAs = "<<nAAs<<", ppAAs = "<<ppAAs<<endl;
        }
        adstring* ptr = new adstring();
        return *ptr;
    }

    adstring_array& adstring_matrix::operator() (const int i){
        if (ppAAs){
            if ((idxmn<=i)&&(i-idxmn<nAAs)&&(ppAAs[i-idxmn])) {
                return *(ppAAs[i-idxmn]);
            } else {
                cout<<"Attempting illegal access in adstring_matrix::("<<i<<")"<<endl;
                cout<<"Row index is invalid! Limits are "<<idxmn<<","<<idxmn+nAAs<<"."<<endl;
            }
        } else {
            cout<<"Attempting illegal access in adstring_matrix::("<<i<<")"<<endl;
            cout<<"Matrix is not allocated!"<<endl;
        }
        adstring_array* ptr = new adstring_array();
        return *ptr;
    }

    adstring_array& adstring_matrix::operator[] (const int i){
        if (ppAAs){
            if ((idxmn<=i)&&(i-idxmn<nAAs)&&(ppAAs[i-idxmn])) {
                return *(ppAAs[i-idxmn]);
            } else {
                cout<<"Attempting illegal access in adstring_matrix::("<<i<<")"<<endl;
                cout<<"Row index is invalid! Limits are "<<idxmn<<","<<idxmn+nAAs<<"."<<endl;
            }
        } else {
            cout<<"Attempting illegal access in adstring_matrix::("<<i<<")"<<endl;
            cout<<"Matrix is not allocated!"<<endl;
        }
        adstring_array* ptr = new adstring_array();
        return *ptr;
    }

    void adstring_matrix::read(cifstream & is){
        if (ppAAs){
            for (int r=0;r<nAAs;r++){
                if (allocated(r)) is>>(*(ppAAs[r]));
            }
        }
    }

    void adstring_matrix::write(ostream & os){
        if (ppAAs){
            for (int r=0;r<nAAs;r++){
                for (int c=clmns(r+idxmn);c<=clmxs(r+idxmn);c++){
                    os<<(*(ppAAs[r]))(c)<<tb;
                }
                if (r<(nAAs-1)) os<<endl<<tb;
            }
        }
    }

    /****************************************************************
    * name      : CompareAdstrings                                  *
    * purpose   : Class with memeber function to determine if       *
    *               lhs < rhs for input adstrings lhs, rhs          *
    ****************************************************************/
    bool CompareAdstrings::debug = false;
    bool CompareAdstrings::operator() (const adstring& lhs, const adstring& rhs) const {
        if (debug) cout<<"CompareAdstrings::('"<<lhs<<"','"<<rhs<<"'): "<<&lhs<<"; "<<&rhs<<endl;
        int nl = lhs.size(); int nr = rhs.size();
        int nz = nr;//compare characters over min string size (assume rhs is shorter)
        if (nl<nr) nz = nl;
        int i = 1;
        bool cont = true;//flag to continue checking
        if (debug) cout<<"CompareAdstrings: "<<lhs<<" < "<<rhs<<"; "<<nl<<", "<<nr<<endl;
        while (cont&&(i<=nz)) {
            cont = (lhs(i)==rhs(i));
            if (debug) cout<<lhs(i)<<" == "<<rhs(i)<<" ?: "<<cont<<endl;
            i++;
        }
        if (debug) cout<<"i: "<<i<<", cont = "<<cont<<endl;
        if (!cont) {
            cont = lhs(i-1)<rhs(i-1);//equality failed above, so which is it?
        } else {
            //comparison was equal up to end of shorter string.
            //set shorter string < longer string
            if (nl<nr) {
                cont = true;//rhs is longer than lhs, so lhs<rhs
            } else {
                cont = false;
            }
        }
        if (debug) cout<<"CompareAdstrings: "<<lhs<<" < "<<rhs<<" ? "<<cont<<endl;
        return cont;
    }

    /****************************************************************
    * name      : IdentityMatrix                                    *
    * purpose   : returns an identity dmatrix with row/col          *
    *               indices running from mn to mx.                  *
    ****************************************************************/
    dmatrix IdentityMatrix(int mn,int mx){dmatrix m(mn,mx); m.initialize(); for (int i=mn;i<=mx;i++) m(i,i)=1.0; return m;}

    /****************************************************************
    * name      : testEquals                                        *
    * purpose   : test equality                                     *
    ****************************************************************/
    int testEquals(ivector& i1, ivector& i2){
        int mn = i1.indexmin();
        if (mn!=i2.indexmin()) return 0;
        int mx = i1.indexmax();
        if (mx!=i2.indexmax()) return 0;
        int rs = mn;
        while ((rs)&&(rs<=mn)) {
            if (i1(rs)==i2(rs)) rs++; else return 0;
        }
        return 1;
    }

    /*************************************************
    * name      : getBounds                          *
    * purpose   : return min, max indices for array  *
    *************************************************/
    ivector getBounds(const prevariable& o){
        ivector b(1,1);
        b(1) = 0;
        return(b);
    }
    /*************************************************
    * name      : getBounds                          *
    * purpose   : return min, max indices for array  *
    *************************************************/
    ivector getBounds(const ivector& o){
        ivector b(1,2);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const imatrix& o){
        ivector b(1,4);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,4)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const dvector& o){
        ivector b(1,2);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const dmatrix& o){
        ivector b(1,4);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,4)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const d3_array& o){
        ivector b(1,6);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,6)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const d4_array& o){
        ivector b(1,8);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,8)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const d5_array& o){
        ivector b(1,10);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,10)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const d6_array& o){
        ivector b(1,12);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,12)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const d7_array& o){
        ivector b(1,14);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,14)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const dvar_vector& o){
        if (debugADMBFunctions>0) cout<<"start getBounds(dvar_vector)"<<endl;
        ivector b(1,2);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        if (debugADMBFunctions>0) cout<<b<<endl;
        if (debugADMBFunctions>0) cout<<"end getBounds(dvar_vector)"<<endl;
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const dvar_matrix& o){
        if (debugADMBFunctions>0) cout<<"start getBounds(dvar_matrix)"<<endl;
        ivector b(1,4);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,4)) = getBounds(o(b(1)));
        if (debugADMBFunctions>0) cout<<b<<endl;
        if (debugADMBFunctions>0) cout<<"end getBounds(dvar_matrix)"<<endl;
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const dvar3_array& o){
        if (debugADMBFunctions>0) cout<<"start getBounds(dvar3_array)"<<endl;
        ivector b(1,6);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,6)) = getBounds(o(b(1)));
        if (debugADMBFunctions>0) cout<<b<<endl;
        if (debugADMBFunctions>0) cout<<"end getBounds(dvar3_array)"<<endl;
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const dvar4_array& o){
        ivector b(1,8);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,8)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const dvar5_array& o){
        ivector b(1,10);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,10)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const dvar6_array& o){
        ivector b(1,12);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,12)) = getBounds(o(b(1)));
        return(b);
    }
    //-------------------------------------------------------------
    ivector getBounds(const dvar7_array& o){
        ivector b(1,14);
        b(1) = o.indexmin();
        b(2) = o.indexmax();
        --(--b(3,14)) = getBounds(o(b(1)));
        return(b);
    }

    /*************************************************
    * name      : length                             *
    * purpose   : return number of elements in array *
    *************************************************/
    int length(const ivector& o){
        ivector b(1,2);
        b = getBounds(o);
        int i = b(2)-b(1)+1;
        return(i);
    }
    //-------------------------------------------------------------
    int length(const imatrix& o){
        ivector b(1,2);
        b = getBounds(o);
        int i = 1;
        for (int j=1;j<=2;j++) i *= (b(2*j)-b(2*j-1)+1);
        return(i);
    }
    //-------------------------------------------------------------
    int length(const dvector& o){
        ivector b(1,2);
        b = getBounds(o);
        int i = b(2)-b(1)+1;
        return(i);
    }
    //-------------------------------------------------------------
    int length(const dmatrix& o){
        ivector b(1,2);
        b = getBounds(o);
        int i = 1;
        for (int j=1;j<=2;j++) i *= (b(2*j)-b(2*j-1)+1);
        return(i);
    }
    //-------------------------------------------------------------
    int length(const dvar_vector& o){
        ivector b(1,2);
        b = getBounds(o);
        int i = b(2)-b(1)+1;
        return(i);
    }
    //-------------------------------------------------------------
    int length(const dvar_matrix& o){
        ivector b(1,2);
        b = getBounds(o);
        int i = 1;
        for (int j=1;j<=2;j++) i *= (b(2*j)-b(2*j-1)+1);
        return(i);
    }
    //-------------------------------------------------------------


    /********************************************
    * name      : mean                          *
    * purpose   : compute mean value of object  *
    ********************************************/
    double mean(dvector & v) {
        return sum(v)/double(v.indexmax()-v.indexmin()+1);
    }
    //-------------------------------------------------------------
    dvariable mean(dvar_vector & v) {
        return sum(v)/double(v.indexmax()-v.indexmin()+1);
    }
    //-------------------------------------------------------------


    /*************************************************
    * name      : value                              *
    * purpose   : return constant version of array   *
    *************************************************/
    d4_array value(const dvar4_array& o){
        if (debugADMBFunctions>0) cout<<"start value(dvar4_array)"<<endl;
        ivector b = getBounds(o);
        d4_array a(b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8));
        for (int r=b(1);r<=b(2);r++) a(r) = value(o(r));
        if (debugADMBFunctions>0) cout<<"   end value(dvar4_array)"<<endl;
        return(a);
    }
    //-------------------------------------------------------------
    d5_array value(const dvar5_array& o){
        if (debugADMBFunctions>0) cout<<"start value(dvar5_array)"<<endl;
        ivector b = getBounds(o);
        d5_array a(b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8),b(9),b(10));
        for (int r=b(1);r<=b(2);r++) a(r) = value(o(r));
        if (debugADMBFunctions>0) cout<<"   end value(dvar5_array)"<<endl;
        return(a);
    }
    //-------------------------------------------------------------
    d6_array value(const dvar6_array& o){
        if (debugADMBFunctions>0) cout<<"start value(dvar6_array)"<<endl;
        ivector b = getBounds(o);
        d6_array a(b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8),b(9),b(10),b(11),b(12));
        for (int r=b(1);r<=b(2);r++) a(r) = value(o(r));
        if (debugADMBFunctions>0) cout<<"   end value(dvar6_array)"<<endl;
        return(a);
    }
    //-------------------------------------------------------------
    d7_array value(const dvar7_array& o){
        if (debugADMBFunctions>0) cout<<"start value(dvar7_array)"<<endl;
        ivector b = getBounds(o);
        d7_array a(b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8),b(9),b(10),b(11),b(12),b(13),b(14));
        for (int r=b(1);r<=b(2);r++) a(r) = value(o(r));
        if (debugADMBFunctions>0) cout<<"   end value(dvar7_array)"<<endl;
        return(a);
    }
    //-------------------------------------------------------------

    /********************************************
    * name      : variance                      *
    * purpose   : compute variance of object    *
    ********************************************/
    double variance(dvector & v) {
        return norm2(v-mean(v))/double(length(v)-1);
    }
    //-------------------------------------------------------------
    dvariable variance(dvar_vector & v) {
        return norm2(v-mean(v))/double(length(v)-1);
    }
    //-------------------------------------------------------------


    /****************************************************************
    * name      : cdf_normal                                        *
    * purpose   : compute cdf of normal distribution                *
    *   parameters:                                                 *
    *       mu : location parameter (mean)                          *
    *       sd:  standard deviation                                 *
    *       lb : lower bound                                        *
    *       ub : upper bound                                        *
    ****************************************************************/
    double cdf_normal(const double& mu,const double& sd,const double& lb,const double& ub){
    //      if (debug20) cout<<endl<<"Compute cdf_normal"<<endl;
        double nf = 0.0;
        if (sd>0) {
            double zn = (lb-mu)/sd;
            double zx = (ub-mu)/sd;
            nf = cumd_norm(zx)-cumd_norm(zn);
        } else {
            if ((lb<=mu)&&(mu<=ub)) nf = 1.0;
        }
        return nf;
    }
    //-------------------------------------------------------------
    dvariable cdf_normal(const dvariable& mu,const double& sd,const double& lb,const double& ub){
        RETURN_ARRAYS_INCREMENT();
    //      if (debug20) cout<<endl<<"Compute cdf_normal"<<endl;
        dvariable nf = 0.0;
        if (sd>0) {
            dvariable zn = (lb-mu)/sd;
            dvariable zx = (ub-mu)/sd;
            nf = cumd_norm(zx)-cumd_norm(zn);
        } else {
            if ((lb<=mu)&&(mu<=ub)) nf = 1.0;
        }
        RETURN_ARRAYS_DECREMENT();
        return nf;
    }
    //-------------------------------------------------------------
    dvariable cdf_normal(const dvariable& mu,const dvariable& sd,const double& lb,const double& ub){
        RETURN_ARRAYS_INCREMENT();
    //        if (debug20) cout<<endl<<"Compute cdf_normal"<<endl;
        dvariable nf = 0.0;
        if (sd>0) {
            dvariable zn = (lb-mu)/sd;
            dvariable zx = (ub-mu)/sd;
            nf = cumd_norm(zx)-cumd_norm(zn);
        } else {
            if ((lb<=mu)&&(mu<=ub)) nf = 1.0;
        }
        RETURN_ARRAYS_DECREMENT();
        return nf;
    }

    /****************************************************************
    * name      : log_normal_density                                *
    * purpose   : compute log of normal distribution                *
    *   parameters:                                                 *
    *       z: normalized random deviate (z=(x-mu)/sdv)             *
    ****************************************************************/
    //-------------------------------------------------------------
    double    log_normal_density(const double& z){
        double d = -0.5*(log(2.0*PI)+square(z));
        return d;
    }
    //-------------------------------------------------------------
    dvariable log_normal_density(const prevariable& z){
        RETURN_ARRAYS_INCREMENT();
        dvariable d = -0.5*(log(2.0*PI)+square(z));
        RETURN_ARRAYS_DECREMENT();
        return d;
    }

    //-------------------------------------------------------------
    dvector log_normal_density(const dvector& z){
        RETURN_ARRAYS_INCREMENT();
        dvector d = -0.5*(log(2.0*PI)+elem_prod(z,z));
        RETURN_ARRAYS_DECREMENT();
        return d;
    }

    //-------------------------------------------------------------
    dvar_vector log_normal_density(const dvar_vector& z){
        RETURN_ARRAYS_INCREMENT();
        dvar_vector d = -0.5*(log(2.0*PI)+elem_prod(z,z));
        RETURN_ARRAYS_DECREMENT();
        return d;
    }

    /****************************************************************
    * name      : log_normal_density                                *
    * purpose   : compute log of normal distribution                *
    *   parameters:                                                 *
    *       x  : value at which to evaluate pdf                     *
    *       mu : mean                                               *
    *       sdv: standard deviation                                 *
    ****************************************************************/
    //-------------------------------------------------------------
    double    log_normal_density(const double& x,const double& mu,const double& sdv){
        double z = (x-mu)/sdv;
        return log_normal_density(z);
    }
    //-------------------------------------------------------------
    dvariable log_normal_density(const prevariable& x,const double& mu,const double& sdv){
        RETURN_ARRAYS_INCREMENT();
        dvariable z = (x-mu)/sdv;
        dvariable d = log_normal_density(z);
        RETURN_ARRAYS_DECREMENT();
        return d;
    }
    //-------------------------------------------------------------
    dvariable log_normal_density(const double& x,const prevariable& mu,const prevariable& sdv){
        RETURN_ARRAYS_INCREMENT();
        dvariable z = (x-mu)/sdv;
        dvariable d = log_normal_density(z);
        RETURN_ARRAYS_DECREMENT();
        return d;
    }
    //-------------------------------------------------------------
    dvariable log_normal_density(const prevariable& x,const prevariable& mu,const prevariable& sdv){
        RETURN_ARRAYS_INCREMENT();
        dvariable z = (x-mu)/sdv;
        dvariable d = log_normal_density(z);
        RETURN_ARRAYS_DECREMENT();
        return d;
    }
    //-------------------------------------------------------------
    dvector     log_normal_density(const dvector& x,const double& mu,const double& sdv){
        dvector z = (x-mu)/sdv;
        return log_normal_density(z);
    }
    //-------------------------------------------------------------
    dvar_vector log_normal_density(const dvar_vector& x,const double& mu,const double& sdv){
        RETURN_ARRAYS_INCREMENT();
        dvar_vector z = (x-mu)/sdv;
        dvar_vector d = log_normal_density(z);
        RETURN_ARRAYS_DECREMENT();
        return d;
    }

    /****************************************************************
    * name      : log_lognormal_density                             *
    * purpose   : compute log of lognormal distribution             *
    *   parameters:                                                 *
    *       x  : arithmetic-scale value                             *
    *       med: arithmetic-scale median                            *
    *       cv: arithmetic-scale cv (stdev(X)/mean(X))              *
    ****************************************************************/
    //-------------------------------------------------------------
    double    log_lognormal_density(const double& x,const double& med,const double& cv){
        double sdv = sqrt(log(1.0+square(cv)));
        double d = -0.5*log(2.0*PI*square(x*sdv))-0.5*square((log(x)-log(med))/sdv);
        return d;
    }
    //-------------------------------------------------------------
    dvariable log_lognormal_density(const prevariable& x,const double& med,const double& cv){
        RETURN_ARRAYS_INCREMENT();
        double sdv = sqrt(log(1.0+square(cv)));
        dvariable d = -0.5*log(2.0*PI*square(x*sdv))-0.5*square((log(x)-log(med))/sdv);
        RETURN_ARRAYS_DECREMENT();
        return d;
    }
    //-------------------------------------------------------------
    dvariable log_lognormal_density(const double& x,const prevariable& med,const prevariable& cv){
        RETURN_ARRAYS_INCREMENT();
        dvariable sdv = sqrt(log(1.0+square(cv)));
        dvariable d = -0.5*log(2.0*PI*square(x*sdv))-0.5*square((log(x)-log(med))/sdv);
        RETURN_ARRAYS_DECREMENT();
        return d;
    }
    //-------------------------------------------------------------
    dvariable log_lognormal_density(const prevariable& x,const prevariable& med,const prevariable& cv){
        RETURN_ARRAYS_INCREMENT();
        dvariable sdv = sqrt(log(1.0+square(cv)));
        dvariable d = -0.5*log(2.0*PI*square(x*sdv))-0.5*square((log(x)-log(med))/sdv);
        RETURN_ARRAYS_DECREMENT();
        return d;
    }
    //-------------------------------------------------------------
    dvector     log_lognormal_density(const dvector& x,const double& med,const double& cv){
        double  sdv = sqrt(log(1.0+square(cv)));
        dvector z = (log(x)-log(med))/sdv;
        dvector d = -0.5*log(2.0*PI)-log(x*sdv)-0.5*elem_prod(z,z);
        return d;
    }
    //-------------------------------------------------------------
    dvar_vector log_lognormal_density(const dvar_vector& x,const double& med,const double& cv){
        RETURN_ARRAYS_INCREMENT();
        double sdv = sqrt(log(1.0+square(cv)));
        dvar_vector z = (log(x)-log(med))/sdv;
        dvar_vector d = -0.5*log(2.0*PI)-log(x*sdv)-0.5*elem_prod(z,z);
        RETURN_ARRAYS_DECREMENT();
        return d;
    }

    /************************************************************************
    * name      : log_gamma_density                                         *
    * purpose   : compute log of gamma pdf                                  *
    * log_gamma_density(x,r,mu) = r*log(mu)-log_gamma(r)+(r-1)*log(x)-mu*x  *
    * gamma(x,r,mu) = (mu^r)/gamma(r) * x^(r-1) * exp(-mu*x)                *
    * This is SAME as Gelman et al., Bayesian Data Analysis                 *
    *   parameters:                                                         *
    *       x : value                                                       *
    *       r : ??                                                          *
    *       mu: rate (inverse scale) parameter                              *
    ************************************************************************/
    //-------------------------------------------------------------
    dvector log_gamma_density(const dvector& xv,const double& r,const double& mu){
        int mn = xv.indexmin();
        int mx = xv.indexmax();
        dvector d(mn,mx);
        for (int i=mn;i<=mx;i++) {
            d(i) = ::log_gamma_density(xv(i),r,mu);
        }
        return d;
    }
    //-------------------------------------------------------------
    dvar_vector log_gamma_density(const dvar_vector& xv,const double& r,const double& mu){
        RETURN_ARRAYS_INCREMENT();
        int mn = xv.indexmin();
        int mx = xv.indexmax();
        dvar_vector d(mn,mx);
        for (int i=mn;i<=mx;i++) {
            d(i) = log_gamma_density(xv(i),r,mu);
        }
        RETURN_ARRAYS_DECREMENT();
        return d;
    }
    //-------------------------------------------------------------
    dvar_vector log_gamma_density(const dvar_vector& xv,const dvariable& r,const dvariable& mu){
        cout<<"Starting log_gamma_density(dvar_vector&, dvariable&, dvariable&)"<<endl;
        RETURN_ARRAYS_INCREMENT();
        int mn = xv.indexmin();
        int mx = xv.indexmax();
        dvar_vector d(mn,mx);
        dvariable xp;
        for (int i=mn;i<=mx;i++) {
            cout<<xv(i)<<tb<<r<<tb<<mu<<tb;
            xp = xv(i)+(1.0e-10);
            d(i) = log_gamma_density(xp,r,mu);
            cout<<d(i)<<endl;
        }
        RETURN_ARRAYS_DECREMENT();
        cout<<"Finished log_gamma_density(dvar_vector&, dvariable&, dvariable&)"<<endl;
        return d;
    }
    /****************************************************************
    * name      : drawSampleLognormal                               *
    * purpose   : draw sample from Lognormal distribution           *
    *   parameters:                                                 *
    *       md : location parameter (arithmetic median)             *
    *       cv:  coefficient of variation     (sd/mean)             *
    ****************************************************************/
    double drawSampleLognormal(random_number_generator& rng, const double md, const double cv) {
        return md*mfexp(randn(rng)*sqrt(log(1+cv*cv)));
    }
    /****************************************************************
    * name      : drawSampleNegBinomial                             *
    * purpose   : draw sample from negative binomial distribution   *
    *   parameters:                                                 *
    *       mu : location parameter (??)                            *
    *       tau:  ??                                                *
    ****************************************************************/
    double drawSampleNegBinomial(random_number_generator& rng, const double mu, const double tau) {
        return randnegbinomial(mu,tau,rng);
    }
    /****************************************************************
    * name      : drawSampleNormal                                  *
    * purpose   : draw sample from Normal distribution              *
    *   parameters:                                                 *
    *       mu : location parameter (mean)                          *
    *       sd:  standard deviation                                 *
    ****************************************************************/
    double drawSampleNormal(random_number_generator& rng, const double mu, const double sd) {
        return mu+randn(rng)*sd;
    }
    /****************************************************************
    * name      : drawSamplePoisson                                 *
    * purpose   : draw sample from Normal distribution              *
    *   parameters:                                                 *
    *       lam: rate parameter                                     *
    ****************************************************************/
    double drawSamplePoisson(random_number_generator& rng, const double mu, const double sd) {
        return mu+randn(rng)*sd;
    }
    /****************************************************************
    * name      : drawSampleUniform                                 *
    * purpose   : draw sample from Uniform distribution             *
    *   parameters:                                                 *
    *       lb : lower bound                                        *
    *       ub:  upper bound                                        *
    ****************************************************************/
    double drawSampleUniform(random_number_generator& rng, const double lb, const double ub) {
        return lb+randn(rng)*(ub-lb);
    }

    /*************************************************
    * name      : getIndexVector                     *
    * purpose   : get indices for vector             *
    *************************************************/
    ivector getIndexVector(dvector& o) {
        ivector bnds = getBounds(o);
        ivector indx(bnds(1),bnds(2));
        indx.fill_seqadd(bnds(1),1);
        return indx;
    }
    //-------------------------------------------------------------

    ivector getIndexVector(dvar_vector& o) {
        ivector bnds = getBounds(o);
        ivector indx(bnds(1),bnds(2));
        indx.fill_seqadd(bnds(1),1);
        return indx;
    }
    //-------------------------------------------------------------


    /********************************************
    * standardized double functions             *
    ********************************************/
    double wts_none  (double x, dvector& consts){return  x;        }
    double wts_acos  (double x, dvector& consts){return  acos(x);  }
    double wts_asin  (double x, dvector& consts){return  asin(x);  }
    double wts_atan  (double x, dvector& consts){return  atan(x);  }
    double wts_cos   (double x, dvector& consts){return  cos(x);   }
    double wts_exp   (double x, dvector& consts){return  mfexp(x); }
    double wts_expneg(double x, dvector& consts){return -mfexp(x); }
    double wts_log   (double x, dvector& consts){return  log(x);   }
    double wts_logneg(double x, dvector& consts){return  log(-x);  }
    /***********************************************************
    wts_logistic
        transform:
            y = min+(max-min)/[1+exp(-x)]
        x range:
            (-Inf,+Inf)
        y range:
            (min, max)
        consts :
            {min, max}
        inverse transform:
            wts_logit(y,{min,max})
    ***********************************************************/
    double wts_logistic(double x, dvector& consts){
        return consts(1)+(consts(2)-consts(1))/(1+mfexp(-x));
    }

    /***********************************************************
    wts_logit
        transform:
            y = ln[(x-min)/(max-x)]
        x range:
            (min, max)
        y range:
            (-Inf,+Inf)
        consts :
            {min, max}
        inverse transform:
            wts_logistic(y,{min,max})
    ***********************************************************/
    double wts_logit   (double x, dvector& consts){
        return log((x-consts(1))/(consts(2)-x));
    }
    double wts_sin   (double x, dvector& consts){return sin(x);   }
    double wts_sqrt  (double x, dvector& consts){return sqrt(x);  }
    double wts_square(double x, dvector& consts){return square(x);}
    double wts_tan   (double x, dvector& consts){return tan(x);   }

    /********************************************
    * standardized prevariable functions        *
    ********************************************/
    dvariable wts_none  (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp= x;        RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_acos  (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp= acos(x);  RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_asin  (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp= asin(x);  RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_atan  (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp= atan(x);  RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_cos   (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp= cos(x);   RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_exp   (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp= mfexp(x); RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_expneg(const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp=-mfexp(x); RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_log   (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp= log(x);   RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_logneg(const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp= log(-x);  RETURN_ARRAYS_DECREMENT();return xp;}
    /***********************************************************
    wts_logistic
        transform:
            y = min+(max-min)/[1+exp(-x)]
        x range:
            (-Inf,+Inf)
        y range:
            (min, max)
        consts :
            {min, max}
        inverse transform:
            wts_logit(y,{min,max})
    ***********************************************************/
    dvariable wts_logistic(const prevariable& x, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable xp=consts(1)+(consts(2)-consts(1))/(1+mfexp(-x));
        RETURN_ARRAYS_DECREMENT();
        return xp;
    }
    /***********************************************************
    wts_logit
        transform:
            y = ln[(x-min)/(max-x)]
        x range:
            (min, max)
        y range:
            (-Inf,+Inf)
        consts :
            {min, max}
        inverse transform:
            wts_logistic(y,{min,max})
    ***********************************************************/
    dvariable wts_logit(const prevariable& x, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable xp=log((x-consts(1))/(consts(2)-x));
        RETURN_ARRAYS_DECREMENT();
        return xp;
    }
    dvariable wts_sin   (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp=sin(x);   RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_sqrt  (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp=sqrt(x);  RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_square(const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp=square(x);RETURN_ARRAYS_DECREMENT();return xp;}
    dvariable wts_tan   (const prevariable& x, dvector& consts){RETURN_ARRAYS_INCREMENT();dvariable xp=tan(x);   RETURN_ARRAYS_DECREMENT();return xp;}

    /********************************************
    * name      : logPDF_?????                  *
    * purpose   : compute ln(pdf(x)) for x      *
    *   params: vector of parameters            *
    *   consts: vector of constants             *
    ********************************************/
    /*----------------------------------------------------*
    *   name      : logPDF_constant                       *
    *   returns 0 as log of value for                     *
    *   improper constant pdf(x) = 1.                     *
    *----------------------------------------------------*/
    dvariable logPDF_constant(prevariable& x,dvar_vector params,dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable val=0.0;
        RETURN_ARRAYS_DECREMENT();
        return val;
    }
    // /*---------------------------------------------------*
    // *   params = alpha, beta                             *
    // *   consts = min, max                                *
    // *---------------------------------------------------*/
    // dvariable logPDF_beta(prevariable& x,dvar_vector params, dvector& consts){
    //     RETURN_ARRAYS_INCREMENT();
    //     dvariable alf  = params(1);
    //     dvariable bta  = params(2);
    //     double A = 0.0;
    //     double B = 1.0;
    //     if (consts.allocated()) {
    //         A = consts(1);
    //         B = consts(2);
    //     }
    //     dvariable y = (y-A)/(B-A);//scale to interval [0,1].
    //     dvariable logPDF = log(beta(alf,bta))+(alf-1)*log(x)+(bta-1)*log(1-x);
    //     RETURN_ARRAYS_DECREMENT();
    //     return logPDF;
    // }
    /*---------------------------------------------------*
    *   pdf(x) = 1/{pi*gamma*[1+((x-x0)/gamma)^2]}       *
    *   params = x0, gamma (location, scale)             *
    *   consts = none                                    *
    *---------------------------------------------------*/
    dvariable logPDF_cauchy(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable x0     = params(1);//location param
        dvariable gamma  = params(2);//scale param (>0)
        dvariable logPDF = -(log(PI*gamma)+log(1+square((x-x0)/gamma)));
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*-------------------------------------------*
    * x~Chisquare(nu) w/ nu dof.                 *
    *   params = nu (dof)                        *
    *   consts = <empty>                         *
    *-------------------------------------------*/
    dvariable logPDF_chisquare(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        /*---------------------------------------------------*
        * if X~Chisquare(nu) then X~Gamma(r=nu/2,mu=1/2)     *
        *       where mu is the rate (1/scale) parameter.    *
        *---------------------------------------------------*/
        dvariable r  = params(1)/2.0; //k in wikipedia article on gamma pdf
        dvariable mu = 0.5;           //1/theta in wikipedia article on gamma pdf
        dvariable logPDF = log_gamma_density(x,r,mu);//chi-square
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*-------------------------------------------*
    * x~Chisquare(nu) w/ nu dof.                 *
    *   params = nu (dof)                        *
    *   consts = <empty>                         *
    *-------------------------------------------*/
    dvar_vector logPDF_chisquare(dvar_vector& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
    //    cout<<"Starting logPDF_chisquare(dvar_vector&, dvar_vector, dvector&)"<<endl;
        /*---------------------------------------------------*
        * if X~Chisquare(nu) then X~Gamma(r=nu/2,mu=1/2)     *
        *       where mu is the rate (1/scale) parameter.    *
        *---------------------------------------------------*/
        dvariable r  = params(1)/2.0; //k in wikipedia article on gamma pdf
        dvariable mu = 0.5;           //1/theta in wikipedia article on gamma pdf
        dvar_vector logPDF(x.indexmin(),x.indexmax());
    //    cout<<r<<tb<<mu<<x<<endl;
        logPDF = log_gamma_density(x,r,mu);//chi-square
    //    cout<<logPDF<<tb<<logPDF.indexmin()<<tb<<logPDF.indexmax()<<endl;
    //    cout<<"Finished logPDF_chisquare(dvar_vector&, dvar_vector, dvector&)"<<endl;
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*---------------------------------------------------*
    * X~Chisquare(nu) w/ nu dof.                         *
    *   params = <empty>                                 *
    *   consts = nu (dof)                                *
    *---------------------------------------------------*/
    dvariable logPDF_chisqdevs(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable r  = consts(1)/2.0; //k in wikipedia article on gamma pdf
        dvariable mu = 0.5;           //1/theta in wikipedia article on gamma pdf
        dvariable logPDF = log_gamma_density(x,r,mu);//chi-square
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*---------------------------------------------------*
    * norm2(x/stdev)~Chisquare(nu)                       *
    *   params = stdev                                   *
    *   consts = <empty>                                 *
    *---------------------------------------------------*/
    dvariable logPDF_chisqdev(dvar_vector& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable X2 = norm2(x)/square(params(1));
        dvariable r  = (length(x)-1)/2.0; //k in wikipedia article on gamma pdf
        dvariable mu = 0.5;               //1/theta in wikipedia article on gamma pdf
        dvariable logPDF = log_gamma_density(X2,r,mu);//chi-square
    //    cout<<logPDF<<tb<<logPDF.indexmin()<<tb<<logPDF.indexmax()<<endl;
    //    cout<<"Finished logPDF_chisquare(dvar_vector&, dvar_vector, dvector&)"<<endl;
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*---------------------------------------------------*
    *   pdf(x) = (1/lambda)*exp(-x/lambda)               *
    *   params = lambda (scale)                          *
    *   consts = none                                    *
    *----------------------------------------------------*/
    dvariable logPDF_exponential(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable lambda = params(1);//scale param
        dvariable logPDF = -log(lambda)-lambda*x;
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    dvar_vector logPDF_exponential(dvar_vector& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable lambda = params(1);//scale param
        dvar_vector logPDF = -log(lambda)-lambda*x;
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*----------------------------------------------------------------------*
    * name      : logPDF_gamma                                              *
    * purpose   : compute log of gamma pdf                                  *
    * log_gamma_density(x,r,mu) = r*log(mu)-log_gamma(r)+(r-1)*log(x)-mu*x  *
    * gamma(x,r,mu) = (mu^r)/gamma(r) * x^(r-1) * exp(-mu*x)                *
    * This is SAME as Gelman et al., Bayesian Data Analysis                 *
    *       x : value                                                       *
    *       r : shape factor                                                *
    *       mu: rate (inverse scale) parameter                              *
    * inputs:                                                               *
    *   params = r,mu (shape, rate)                                         *
    *   consts = <empty>                                                    *
    *  OR                                                                   *
    *   params = <empty>                                                    *
    *   consts = r,mu (shape, rate)                                         *
    *-----------------------------------------------------------------------*/
    dvariable logPDF_gamma(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable logPDF;
        if (params.indexmax()==2) {
            dvariable r  = params(1); //shape parameter: k in wikipedia article on gamma pdf
            dvariable mu = params(2); //rate parameter : 1/theta in wikipedia article on gamma pdf
            logPDF = log_gamma_density(x,r,mu);
        } else {
            double r  = consts(1); //shape parameter: k in wikipedia article on gamma pdf
            double mu = consts(2); //rate parameter : 1/theta in wikipedia article on gamma pdf
            logPDF = log_gamma_density(x,r,mu);
        }
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    dvar_vector logPDF_gamma(dvar_vector& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvar_vector logPDF(x.indexmin(),x.indexmax());
        if (params.indexmax()==2) {
            dvariable r  = params(1); //shape parameter: k in wikipedia article on gamma pdf
            dvariable mu = params(2); //rate parameter : 1/theta in wikipedia article on gamma pdf
            logPDF = log_gamma_density(x,r,mu);
        } else {
            double r  = consts(1); //shape parameter: k in wikipedia article on gamma pdf
            double mu = consts(2); //rate parameter : 1/theta in wikipedia article on gamma pdf
            logPDF = log_gamma_density(x,r,mu);
        }
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*-----------------------------------------------*
    * if X~InvChisquare(nu) then Y=1/X~Chisquare(nu) *
    *   params = nu (dof)                            *
    *   consts = <empty>                             *
    *------------------------------------------------*/
    dvariable logPDF_invchisquare(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        //if X~InvChisquare(r) then 1/X~Chisquare(r)
        dvariable y = 1/x; //y is chisquare-distributed
        dvariable logPDF = logPDF_chisquare(y,params,consts);
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*---------------------------------------------*
    * if X~InvGamma(r,mu) then Y=1/X~Gamma(r,1/mu) *
    *   params = r,mu (shape, rate)                *
    *   consts = <empty>                           *
    *----------------------------------------------*/
    dvariable logPDF_invgamma(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        //if X~InvGamma(r,mu) then 1/X~Gamma(r,1/mu)
        dvariable y     = 1/x;         //gamma-distributed
        dvariable r     = params(1);   //shape parameter: same as gamma
        dvariable invmu = 1/params(2); //gamma rate parameter: recip of mu
        dvariable logPDF = log_gamma_density(y,r,invmu);
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*---------------------------------------------------*
    *   pdf(x) = sqrt(lambda/(2*pi*x^3))*                *
    *                    exp(-(lambda/(2*x))*(x/mu-1)^2) *
    *   params = mu, lambda (location, shape)            *
    *   consts = none                                    *
    *----------------------------------------------------*/
    dvariable logPDF_invgaussian(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable mu     = params(1);//location param
        dvariable lambda = params(2);//shape param
        dvariable logPDF = 0.5*(log(lambda)-log(2*PI)-3*log(x))-lambda/(2*x)*square(x/mu-1);
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*----------------------------------------------------------*
    *   pdf(x) =1/x*(2*pi*sigma^2)^-0.5 *                       *
    *                     exp(-0.5*((ln(x)-ln(mu))/sigma)^2)    *
    *   params = mu, cv (median, cv)                            *
    *   consts = none                                           *
    *-----------------------------------------------------------*/
    dvariable logPDF_lognormal(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable mu     = params(1);                 //median
        dvariable var    = log(1.0+square(params(2)));//log-scale variance
        dvariable logPDF = -log(x)-0.5*(log(2.0*PI)+log(var)+
                                            square(log(x)-log(mu))/var);
    //    cout<<"logPDF_lognormal: "<<x<<tb<<params<<tb<<sigma<<tb<<logPDF<<endl;
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*----------------------------------------------------------*
    *   pdf(x) = (2*pi*sigma^2)^-0.5 *                          *
    *                     exp(-0.5*((ln(x)-ln(mu))/sigma)^2)    *
    *   params = mu, cv (median, cv)                            *
    *   consts = none                                           *
    *-----------------------------------------------------------*/
    dvariable logPDF_logscale_normal(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable mu     = params(1);                 //median
        dvariable var    = log(1.0+square(params(2)));//log-scale variance
        dvariable logPDF = -0.5*(log(2.0*PI)+log(var)+
                                            square(log(x)-log(mu))/var);
    //    cout<<"logPDF_logscale_normal: "<<x<<tb<<params<<tb<<sigma<<tb<<logPDF<<endl;
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*------------------------------------------------------------*
    *   pdf(x) = (2*pi*sigma^2)^-0.5 * exp(-0.5*((x-mu)/sigma)^2) *
    *   params = mu, sigma (mean, stdev)                          *
    *   consts = none                                             *
    *-------------------------------------------------------------*/
    dvariable logPDF_normal(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable mu     = params(1);//location param
        dvariable sigma  = params(2);//scale param
        dvariable logPDF = -0.5*(log(2.0*PI)+2.0*log(sigma)+square((x-mu)/sigma));
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*-----------------------------------------------------*
    * if X~InvChisquare(nu,s^2) then                       *
    *    Y=1/X~Gamma(nu/2,nu/2*s^2)                        *
    * and                                                  *
    *   pdf_InvGamma(X;nu,s^2) =                           *
    *       (X^-2)*pdf_Gamma(1/X;alpha=nu/2,beta=nu/2*s^2) *
    * inputs:                                              *
    *   params = nu (dof), s                               *
    *   consts = <empty>                                   *
    * or                                                   *
    *   params = <empty>                                   *
    *   consts = nu (dof), s                               *
    *-----------------------------------------------------*/
    dvariable logPDF_scaled_invchisquare(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable logPDF;
        dvariable y = 1/x;
        if (params.indexmax()==2){
            dvariable a = params(1)/2;
            dvariable b = (params(1)/2)*square(params(2));
            logPDF = log_gamma_density(y,a,b)-2*log(x);
        } else {
            dvariable a = consts(1)/2;
            dvariable b = (consts(1)/2)*square(consts(2));
            logPDF = log_gamma_density(y,a,b)-2*log(x);
        }
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*-----------------------------------------------------*
    * x = log(1+cv^2)                                      *
    * if X~InvChisquare(nu,s^2) then                       *
    *    Y=1/X~Gamma(nu/2,nu/2*s^2)                        *
    * and                                                  *
    *   pdf_InvGamma(X;nu,s^2) =                           *
    *       (X^-2)*pdf_Gamma(1/X;scale=nu/2,rate=nu/2*s^2) *
    * inputs:                                              *
    *   params = nu (dof), CV(s) = sqrt(exp(s^2)-1)        *
    *   consts = <empty>                                   *
    * or                                                   *
    *   params = <empty>                                   *
    *   consts = nu (dof), CV(s) = sqrt(exp(s^2)-1)        *
    *-----------------------------------------------------*/
    dvariable logPDF_scaledCV_invchisquare(prevariable& cv,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable logPDF;
        dvariable x = log(1.0+square(cv));
        dvariable y = 1/x;
        if (params.indexmax()==2){
            dvariable a = params(1)/2;
            dvariable b = (params(1)/2)*log(1+square(params(2)));
            logPDF = log_gamma_density(y,a,b)+2*log(y);
        } else {
            dvariable a = consts(1)/2;
            dvariable b = (consts(1)/2)*log(1+square(consts(2)));
            logPDF = log_gamma_density(y,a,b)+2*log(y);
        }
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    /*------------------------------------------------------*
    *   pdf(x) = gamma((vu+1)/2)/(sqrt(nu*pi)*gamma(nu/2))* *
    *               [1+x^2/nu]^(nu+1)/2                     *
    *   params = nu (dof)                                   *
    *   consts = none                                       *
    *-------------------------------------------------------*/
    dvariable logPDF_t(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable nu     = params(1);//location param
        dvariable logPDF = gammln((nu+1.0)/2.0)-gammln(nu/2.0)-0.5*log(nu*PI)+
                            (nu+1.0)/2.0*log(1+square(x)/2.0);
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }

    dvariable logSquareWave(prevariable& x,double& min,double& max,double m){
        RETURN_ARRAYS_INCREMENT();
        dvariable logSW  = -m*(log(1+exp(-(x-min)))+log(1+exp(-(max-x))));
        RETURN_ARRAYS_DECREMENT();
        return logSW;
    }


    /*-------------------------------------------------------------*
    *   pdf(x) = f*(2*pi*sigma^2)^-0.5 * exp(-0.5*((x-mu)/sigma)^2)*
    *          on interval {min,max} [f is a coefficient so the    *
    *          integral from min to max = 1.                       *
    *   params = mu, sigma                                         *
    *   consts = min, max                                          *
    *NOTE: This gives the same value as logPDF_normal since f is a *
    *      constant and immaterial to derivative calcs.            *
    *-------------------------------------------------------------*/
    dvariable logPDF_truncated_normal(prevariable& x,dvar_vector params, dvector& consts){
        RETURN_ARRAYS_INCREMENT();
        dvariable mu     = params(1);//location param
        dvariable sigma  = params(2);//scale param
        double min = consts(1);
        double max = consts(2);
        dvariable logPDF = -0.5*(log(2.0*PI)+2.0*log(sigma)+square((x-mu)/sigma));//+logSquareWave(x,min,max);
        RETURN_ARRAYS_DECREMENT();
        return logPDF;
    }
    //-------------------------------------------------------------

    /********************************************
    * name      : samplePDF_?????               *
    * purpose   : draw sample from pdf(x))      *
    *   params: vector of parameters            *
    *   consts: vector of constants             *
    ********************************************/
    // /*-------------------------------------------------------------*
    // * name      : samplePDF_beta                                   *
    // * purpose   : draw sample from a beta distribution             *
    // *   params:                                                    *
    // *     1 alpha : shape parameter (>0)                           *
    // *     2 beta  : shape parameter (>0)                           *
    // *   consts:                                                    *
    // *     1 A : minimum (default=0)                                *
    // *     2 B : maximum (default=1)                                *
    // *-------------------------------------------------------------*/
    // double samplePDF_beta(random_number_generator& rng,dvector& params,dvector& consts){
    //     double p   = randu(rng);
    //     double alf = params(1);
    //     double bta = params(2);
    //     double val = betai(alf,bta,p);//TODO: this is ASSUMED to be the inverse beta cdf
    //     if (consts.allocated()) val = consts(1)+val*(consts(2)-consts(1));
    //     return val;
    // }
    /*-------------------------------------------------------------*
    * name      : samplePDF_cauchy                                 *
    * purpose   : draw sample from a cauchy distribution           *
    *   params:                                                    *
    *     1 x0     : location parameter (real)                     *
    *     2 gamma  : scale parameter (>0)                          *
    *-------------------------------------------------------------*/
    double samplePDF_cauchy(random_number_generator& rng,dvector& params,dvector& consts){
        double p  = randu(rng);
        double x0 = params(1);
        double gm = params(2);
        double val = x0+gm*tan(PI*(p-0.5));
        return val;
    }
    /*-------------------------------------------------------------*
    * name      : samplePDF_chisquare                              *
    * purpose   : draw sample from a chisquare distribution        *
    *   params:                                                    *
    *     1 k     : degrees of freedom (integer)                   *
    *-------------------------------------------------------------*/
    double samplePDF_chisquare(random_number_generator& rng,dvector& params,dvector& consts){
        int k = (int)params(1);
        dvector p(1,k);
        p.fill_randn(rng);
        double val = norm2(p);
        return val;
    }
    /*-------------------------------------------------------------*
    * name      : samplePDF_exponential                            *
    * purpose   : draw sample from an exponential distribution     *
    *   params:                                                    *
    *     1 lambda : rate parameter (real)                         *
    *-------------------------------------------------------------*/
    double samplePDF_exponential(random_number_generator& rng,dvector& params,dvector& consts){
        double p  = randu(rng);
        double lm = params(1);
        double val = -log(p)/lm;
        return val;
    }
    /*-------------------------------------------------------------*
    * name      : samplePDF_gamma                                  *
    * purpose   : draw sample from a gamma distribution            *
    *   params:                                                    *
    *     1 r     : shape parameter (>0)                           *
    *     2 mu    : rate (inverse scale) parameter (>0)            *
    *-------------------------------------------------------------*/
    double samplePDF_gamma(random_number_generator& rng,dvector& params,dvector& consts){
        double p  = randu(rng);
        double r  = params(1);
        double mu = params(2);
        double val = inv_cumd_gamma(p,r)/mu;//TODO: check inv_cumd_gamma
        return val;
    }
    /*---------------------------------------------------------------*
    * name      : samplePDF_invchisquare                             *
    * purpose   : draw sample from an inverse chisquare distribution *
    *   params:                                                      *
    *     1 k     : degrees of freedom (integer)                     *
    *---------------------------------------------------------------*/
    double samplePDF_invchisquare(random_number_generator& rng,dvector& params,dvector& consts){
        int k = (int)params(1);
        dvector p(1,k);
        p.fill_randn(rng);
        double val = 1.0/norm2(p);
        return val;
    }
    /*-------------------------------------------------------------*
    * name      : samplePDF_invgamma                               *
    * purpose   : draw sample from an inverse gamma distribution   *
    *   params:                                                    *
    *     1 r     : shape parameter (>0)                           *
    *     2 mu    : rate parameter (>0)                            *
    *-------------------------------------------------------------*/
    double samplePDF_invgamma(random_number_generator& rng,dvector& params,dvector& consts){
        dvector gparams(1,2);
        gparams(1) = params(1);    //r
        gparams(2) = 1.0/params(2);//1/mu
        double val = 1.0/samplePDF_gamma(rng,gparams,consts);//using Y~Gamma(r,1/mu) -> X=1/Y~InvGamma(r,mu)
        return val;
    }
    /*---------------------------------------------------------------*
    * name      : samplePDF_invgaussian                              *
    * purpose   : draw sample from an inverse gaussian distribution  *
    *   params:                                                      *
    *     1 mu     : location parameter (mean; >0)                   *
    *     2 lambda : shape parameter (>0)                            *
    *----------------------------------------------------------------*/
    double samplePDF_invgaussian(random_number_generator& rng,dvector& params,dvector& consts){
        double val = 0.0;
        double n = randn(rng);
        double z = randu(rng);
        double mu = params(1);
        double lm = params(2);
        double x = mu+(square(mu*n)-mu*sqrt((4*mu*lm*square(n))+square(mu*square(n))))/(2*lm);
        if (z<=mu/(mu+x)) val = x; else val = square(mu)/x;
        return val;
    }
    /*-------------------------------------------------------------*
    * name      : samplePDF_lognormal                              *
    * purpose   : draw sample from a lognormal distribution        *
    *   params:                                                    *
    *     1 med : location parameter (median)                      *
    *     2 cv  : standard deviation                               *
    *-------------------------------------------------------------*/
    double samplePDF_lognormal(random_number_generator& rng,dvector& params,dvector& consts){
        double med = params(1);
        double sdv = sqrt(log(1+square(params(2))));
        double val =  med*exp(randn(rng)*sdv);
        return val;
    }
    /*--------------------------------------------------------------*
    * name      : samplePDF_normal                                  *
    * purpose   : draw sample from a normal distribution            *
    *   params:                                                     *
    *     1 mu : location parameter (mean)                          *
    *     2 sd:  standard deviation                                 *
    *--------------------------------------------------------------*/
    double samplePDF_normal(random_number_generator& rng,dvector& params,dvector& consts){
        double mu = params(1);
        double sd = params(2);
        return mu+randn(rng)*sd;
    }
    /*-----------------------------------------------------*
    * name      : samplePDF_scaled_invchisquare            *
    * purpose   : draw sample from                         *
    *               scaled inverse chi-square distribution *
    * if X~InvChisquare(nu,sig^2) then                     *
    *    Y=1/X~Gamma(nu/2,nu/2*s^2)                        *
    * inputs:                                              *
    *   params = nu (dof), s                               *
    *   consts = <empty>                                   *
    * or                                                   *
    *   params = <empty>                                   *
    *   consts = nu (dof), s                               *
    *-----------------------------------------------------*/
    double samplePDF_scaled_invchisquare(random_number_generator& rng,dvector& params,dvector& consts){
        double p  = randu(rng);
        double r  = params(1)/2;                  //shape factor (=beta in Bayesian Data Analysis, p 575)
        double mu = params(1)/2*square(params(2));//rate (inverse scale) factor (=beta in Bayesian Data Analysis, p 575)
        double val = 1/(inv_cumd_gamma(p,r)/mu);
        return val;
    }
    /*-----------------------------------------------------*
    * name      : samplePDF_scaledCV_invchisquare          *
    * purpose   : draw sample from                         *
    *               scaled inverse chi-square distribution *
    *               as a cv                                *
    *------------------------------------------------------*
    * cv = sqrt(exp(x)-1)                                  *
    * if X~InvChisquare(nu,s^2) then                       *
    *    Y=1/X~Gamma(nu/2,nu/2*s^2)                        *
    * and                                                  *
    *   pdf_InvGamma(X;nu,s^2) =                           *
    *       (X^-2)*pdf_Gamma(1/X;scale=nu/2,rate=nu/2*s^2) *
    * inputs:                                              *
    *   params = nu (dof), CV = sqrt(exp(s^2)-1)           *
    *   consts = <empty>                                   *
    * or                                                   *
    *   params = <empty>                                   *
    *   consts = nu (dof), CV = sqrt(exp(s^2)-1)           *
    *-----------------------------------------------------*/
    double samplePDF_scaledCV_invchisquare(random_number_generator& rng,dvector& params,dvector& consts){
        double p  = randu(rng);
        double r  = params(1)/2;                         //shape factor (=beta in Bayesian Data Analysis, p 575)
        double mu = params(1)/2*log(1+square(params(2)));//rate (inverse scale) factor (=beta in Bayesian Data Analysis, p 575)
        double x = 1.0/(inv_cumd_gamma(p,r)/mu);
        double val = sqrt(exp(x)-1.0);                   //cv corresponding to x
        return val;
    }
    /*--------------------------------------------------------------*
    * name      : samplePDF_t                                       *
    * purpose   : draw sample from t distribution                   *
    *   params:                                                     *
    *     1 nu : degrees of freedom                                 *
    *--------------------------------------------------------------*/
    double samplePDF_t(random_number_generator& rng,dvector& params,dvector& consts){
        double nu = params(1);
        double z = randn(rng);
        double v = samplePDF_chisquare(rng,params,consts);
        double val = z/sqrt(v/nu);
        return val;
    }

    /*--------------------------------------------------------------*
    * name      : samplePDF_truncated_normal                        *
    * purpose   : draw sample from a truncated normal distribution  *
    *             on the open interval (min,max)                    *
    *   params:                                                     *
    *     1 mu : location parameter (mode of underlying normal)     *
    *     2 sd:  scale parameter    (stdv of underlying normal)     *
    *   consts:                                                     *
    *     1 min : min value                                         *
    *     2 max : max value                                         *
    *---------------------------------------------------------------*/
    double samplePDF_truncated_normal(random_number_generator& rng,dvector& params,dvector& consts){
        double val = params(1)+randn(rng)*params(2);
        while ((val<=consts(1))||(val>=consts(2))) {
            val = params(1)+randn(rng)*params(2);
        }
        return val;
    }
    //delta function    
    int deltafcn(int i, int j){return (i==j);}

    /******************************************************
     * Tests if a parameter is currently being estimated. *
    ******************************************************/
    bool isActive(param_init_number& p){return active(p);}
    bool isActive(param_init_vector& p){return active(p);}
    bool isActive(param_init_number_vector& p){
        bool res = false;
        for (int i=p.indexmin();i<=p.indexmax();i++) res = res||active(p(i));
        return res;
    }
    bool isActive(param_init_vector_vector& p){
        bool res = false;
        for (int i=p.indexmin();i<=p.indexmax();i++) res = res||active(p(i));
        return res;
    }
    
} // namespace wts