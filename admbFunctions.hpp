// admbFunctions.hpp
// file contains some general ADMB functions
#pragma once
#ifndef __HPP_ADMBFUNCTIONS__
    #define __HPP_ADMBFUNCTIONS__

    extern int debugADMBFunctions;
namespace wts{
    /****************************************************************
     * convert vectors to quoted csv string
     ***************************************************************/
    adstring to_qcsv(const ivector& v);
    adstring to_qcsv(const dvector& v);
    /****************************************************************
     * convert dvector to ivector.
     ***************************************************************/
    ivector to_ivector(const dvector& v);
    
    /****************************************************************
    * name      : copy                                              *
    * purpose   : create a deep copy of a variable                  *
    *NOTE: Does not work for ragged arrays                          *
    ****************************************************************/
    dvector  copy(dvector& v);
    dmatrix  copy(dmatrix& v);
    d3_array copy(d3_array& v);
    d4_array copy(d4_array& v);
    d5_array copy(d5_array& v);
    d6_array copy(d6_array& v);
    d7_array copy(d7_array& v);
    dvar_vector copy(dvar_vector& v);
    dvar_matrix copy(dvar_matrix& v);
    dvar3_array copy(dvar3_array& v);
    dvar4_array copy(dvar4_array& v);
    dvar5_array copy(dvar5_array& v);
    dvar6_array copy(dvar6_array& v);
    dvar7_array copy(dvar7_array& v);

    /****************************************************************
    * name      : adstring_matrix                                   *
    * purpose   : functionality for (possibly ragged) matrix of     *
    *               of adstrings                                    *
    ****************************************************************/
    class adstring_matrix {
        public:
            static bool debug;
        protected:
            int nAAs;      //number of "rows"
            int idxmn;     //min (row) index
            ivector clmns; //min column indices
            ivector clmxs; //max column indices
            adstring_array** ppAAs;
        protected:
            void allocate(int rwmn, int rwmx);
        public:
            adstring_matrix();
            adstring_matrix(int rwmn, int rwmx, int clmn, int clmx);
            adstring_matrix(int rwmn, int rwmx, int clmn, ivector& clmxs);
            adstring_matrix(int rwmn, int rwmx, ivector& clmns, int clmx);
            adstring_matrix(int rwmn, int rwmx, ivector& clmns, ivector& clmxs);
            ~adstring_matrix(){deallocate();}
            bool allocated(){return (ppAAs!=0);}
            bool allocated(int rw){if ((ppAAs)&&(idxmn<=rw)&&(rw-idxmn<nAAs)){return (ppAAs[rw-idxmn]!=0);} return false;}
            void allocate(int rwmn, int rwmx, int clmn, int clmx);
            void allocate(int rwmn, int rwmx, ivector& clmns, int clmx);
            void allocate(int rwmn, int rwmx, int clmn, ivector& clmxs);
            void allocate(int rwmn, int rwmx, ivector& clmns, ivector& clmxs);
            void deallocate(void);
            adstring& operator() (const int i, const int j);
            adstring_array& operator() (const int i);
            adstring_array& operator[] (const int i);
            int indexmin(){return idxmn;}
            int indexmax(){return idxmn+nAAs;}
            int size(){return nAAs;}
            void to_lower(){for (int i=0;i<nAAs;i++) (*ppAAs[i]).to_lower();}
            void read(cifstream & is);
            void write(ostream & os);
            friend cifstream & operator>> (cifstream & is, adstring_matrix & obj){obj.read(is); return is;}
            friend ostream & operator<< (ostream & os, adstring_matrix & obj){obj.write(os); return os;}
    };

    /****************************************************************
    * name      : CompareAdstrings                                  *
    * purpose   : Class with member function to determine if        *
    *               lhs < rhs for input adstrings lhs, rhs          *
    ****************************************************************/
    class CompareAdstrings {
        public:
            static bool debug;
            bool operator() (const adstring& lhs, const adstring& rhs) const;
    };

    /****************************************************************
    * name      : IdentityMatrix                                    *
    * purpose   : returns an identity dmatrix with row/col          *
    *               indices running from mn to mx.                  *
    ****************************************************************/
    dmatrix IdentityMatrix(int mn,int mx);

    /****************************************************************
    * name      : testEquals                                        *
    * purpose   : test equality                                     *
    ****************************************************************/
    int testEquals(ivector& i1, ivector& i2);

    /****************************************************************
    * name      : log_normal_density                                *
    * purpose   : compute log of normal distribution                *
    *   parameters:                                                 *
    *       z: normalized random deviate (z=(x-mu)/sdv              *
    ****************************************************************/
    double      log_normal_density(const double& z);
    dvariable   log_normal_density(const prevariable& z);
    dvector     log_normal_density(const dvector& z);
    dvar_vector log_normal_density(const dvar_vector& z);
    /****************************************************************
    * name      : log_normal_density                                *
    * purpose   : compute log of normal distribution                *
    *   parameters:                                                 *
    *       x  : arithmetic-scale value                             *
    *       mu : mean                                               *
    *       sdv: standard deviation                                 *
    ****************************************************************/
    double      log_normal_density(const double& x,const double& mu,const double& sdv);
    dvariable   log_normal_density(const prevariable& x,const double& mu,const double& sdv);
    dvariable   log_normal_density(const double& x,const prevariable& mu,const prevariable& sdv);
    dvariable   log_normal_density(const prevariable& x,const prevariable& mu,const prevariable& sdv);
    dvector     log_normal_density(const dvector& x,const double& mu,const double& sdv);
    dvar_vector log_normal_density(const dvar_vector& x,const double& mu,const double& sdv);

    /****************************************************************
    * name      : log_lognormal_density                             *
    * purpose   : compute log of lognormal distribution             *
    *   parameters:                                                 *
    *       x  : arithmetic-scale value                             *
    *       med: arithmetic-scale median                            *
    *       sdv: log-scale standard deviation                       *
    ****************************************************************/
    double      log_lognormal_density(const double& x,const double& med,const double& sdv);
    dvariable   log_lognormal_density(const prevariable& x,const double& med,const double& sdv);
    dvariable   log_lognormal_density(const double& x,const prevariable& med,const prevariable& sdv);
    dvariable   log_lognormal_density(const prevariable& x,const prevariable& med,const prevariable& sdv);
    dvector     log_lognormal_density(const dvector& x,const double& med,const double& sdv);
    dvar_vector log_lognormal_density(const dvar_vector& x,const double& med,const double& sdv);

/************************************************************************
* name      : log_gamma_density                                         *
* purpose   : compute log of gamma pdf                                  *
* log_gamma_density(x,r,mu) = r*log(mu)-log_gamma(r)+(r-1)*log(x)-mu*x  *
* gamma(x,r,mu) = (mu^r)/gamma(r) * x^(r-1) * exp(-mu*x)                *
* This is SAME as Gelman et al., Bayesian Data Analysis                 *
*   parameters:                                                         *
*       x : value                                                       *
*       r : shape factor                                                *
*       mu: rate (inverse scale) parameter                              *
************************************************************************/
    dvector     log_gamma_density(const dvector& x,const double& r,const double& mu);
    dvar_vector log_gamma_density(const dvar_vector& x,const double& r,const double& mu);
    dvar_vector log_gamma_density(const dvar_vector& x,const dvariable& r,const dvariable& mu);

    /****************************************************************
    * name      : cdf_normal                                        *
    * purpose   : compute cdf of normal distribution                *
    *   parameters:                                                 *
    *       mu : location parameter (mean)                          *
    *       sd:  standard deviation                                 *
    *       lb : lower bound                                        *
    *       ub : upper bound                                        *
    ****************************************************************/
    //-------------------------------------------------------------
    double cdf_normal(const double& mu,const double& sd,const double& lb,const double& ub);
    dvariable cdf_normal(const dvariable& mu,const double& sd,const double& lb,const double& ub);
    dvariable cdf_normal(const dvariable& mu,const dvariable& sd,const double& lb,const double& ub);

    /****************************************************************
    * name      : drawSampleLognormal                               *
    * purpose   : draw sample from Lognormal distribution           *
    *   parameters:                                                 *
    *       md : location parameter (arithmetic median)             *
    *       cv:  coefficient of variation     (sd/mean)             *
    ****************************************************************/
    double drawSampleLognormal(random_number_generator& rng, const double md, const double cv);
    /****************************************************************
    * name      : drawSampleNegBinomial                             *
    * purpose   : draw sample from negative binomial distribution   *
    *   parameters:                                                 *
    *       mu : location parameter (??)                            *
    *       tau:  ??                                                *
    ****************************************************************/
    double drawSampleNegBinomial(random_number_generator& rng, const double mu, const double tau);
    /****************************************************************
    * name      : drawSampleNormal                                  *
    * purpose   : draw sample from Normal distribution              *
    *   parameters:                                                 *
    *       mu : location parameter (mean)                          *
    *       sd:  standard deviation                                 *
    ****************************************************************/
    double drawSampleNormal(random_number_generator& rng, const double mu, const double sd);
    /****************************************************************
    * name      : drawSamplePoisson                                 *
    * purpose   : draw sample from Normal distribution              *
    *   parameters:                                                 *
    *       lam: rate parameter                                     *
    ****************************************************************/
    double drawSamplePoisson(random_number_generator& rng, const double mu, const double sd);
    /****************************************************************
    * name      : drawSampleUniform                                 *
    * purpose   : draw sample from Uniform distribution             *
    *   parameters:                                                 *
    *       lb : lower bound                                        *
    *       ub:  upper bound                                        *
    ****************************************************************/
    double drawSampleUniform(random_number_generator& rng, const double lb, const double ub);

    /*************************************************
    * name      : getBounds                          *
    * purpose   : return min, max indices for array  *
    *************************************************/
    ivector getBounds(const prevariable& o);
    ivector getBounds(const ivector& o);
    ivector getBounds(const imatrix& o);
    ivector getBounds(const dvector& o);
    ivector getBounds(const dmatrix& o);
    ivector getBounds(const d3_array& o);
    ivector getBounds(const d4_array& o);
    ivector getBounds(const d5_array& o);
    ivector getBounds(const d6_array& o);
    ivector getBounds(const d7_array& o);
    ivector getBounds(const dvar_vector& o);
    ivector getBounds(const dvar_matrix& o);
    ivector getBounds(const dvar3_array& o);
    ivector getBounds(const dvar4_array& o);
    ivector getBounds(const dvar5_array& o);
    ivector getBounds(const dvar6_array& o);
    ivector getBounds(const dvar7_array& o);

    /*************************************************
    * name      : getIndexVector                     *
    * purpose   : get vector of indices for vector   *
    *************************************************/
    ivector getIndexVector(dvector& o);
    ivector getIndexVector(dvar_vector& o);

    /*************************************************
    * name      : length                             *
    * purpose   : return number of elements in array *
    *************************************************/
    int length(const ivector& o);
    int length(const imatrix& o);
    int length(const dvector& o);
    int length(const dmatrix& o);
    int length(const dvar_vector& o);
    int length(const dvar_matrix& o);

    /********************************************
    * name      : mean                          *
    * purpose   : compute mean value of object  *
    ********************************************/
    double mean(dvector & v);
    dvariable mean(dvar_vector & v);

    /*************************************************
    * name      : value                              *
    * purpose   : return constant version of array   *
    *************************************************/
    //-------------------------------------------------------------
    d4_array value(const dvar4_array& o);
    d5_array value(const dvar5_array& o);
    d6_array value(const dvar6_array& o);
    d7_array value(const dvar7_array& o);

    /********************************************
    * name      : variance                      *
    * purpose   : compute variance of object    *
    ********************************************/
    double variance(dvector & v);
    dvariable variance(dvar_vector & v);

    /********************************************
    * standardized double functions             *
    ********************************************/
    double wts_none    (double x, dvector& consts);
    double wts_acos    (double x, dvector& consts);
    double wts_asin    (double x, dvector& consts);
    double wts_atan    (double x, dvector& consts);
    double wts_cos     (double x, dvector& consts);
    double wts_exp     (double x, dvector& consts);
    double wts_expneg  (double x, dvector& consts);
    double wts_log     (double x, dvector& consts);
    double wts_logneg  (double x, dvector& consts);
    double wts_logistic(double x, dvector& consts);
    double wts_logit   (double x, dvector& consts);
    double wts_sin     (double x, dvector& consts);
    double wts_sqrt    (double x, dvector& consts);
    double wts_square  (double x, dvector& consts);
    double wts_tan     (double x, dvector& consts);

    /********************************************
    * standardized prevariable functions        *
    ********************************************/
    dvariable wts_none    (const prevariable& x, dvector& consts);
    dvariable wts_acos    (const prevariable& x, dvector& consts);
    dvariable wts_asin    (const prevariable& x, dvector& consts);
    dvariable wts_atan    (const prevariable& x, dvector& consts);
    dvariable wts_cos     (const prevariable& x, dvector& consts);
    dvariable wts_exp     (const prevariable& x, dvector& consts);
    dvariable wts_expneg  (const prevariable& x, dvector& consts);
    dvariable wts_log     (const prevariable& x, dvector& consts);
    dvariable wts_logneg  (const prevariable& x, dvector& consts);
    dvariable wts_logistic(const prevariable& x, dvector& consts);
    dvariable wts_logit   (const prevariable& x, dvector& consts);
    dvariable wts_sin     (const prevariable& x, dvector& consts);
    dvariable wts_sqrt    (const prevariable& x, dvector& consts);
    dvariable wts_square  (const prevariable& x, dvector& consts);
    dvariable wts_tan     (const prevariable& x, dvector& consts);
    
    dvariable logSquareWave(prevariable& x,double& min,double& max,double m=100.0);


/********************************************
* name      : logPDF_?????                  *
* purpose   : compute ln(pdf(x)) for x      *
*   params: vector of parameters            *
*   consts: vector of constants             *
********************************************/
dvariable logPDF_constant    (prevariable& x,dvar_vector params,dvector& consts);
//dvariable logPDF_beta        (prevariable& x,dvar_vector params,dvector& consts);
dvariable logPDF_cauchy      (prevariable& x,dvar_vector params,dvector& consts);
/*-------------------------------------------*
* x~Chisquare(nu) w/ nu dof.                 *
*   params = nu (dof)                        *
*   consts = <empty>                         *
*-------------------------------------------*/
dvariable logPDF_chisquare(prevariable& x,dvar_vector params,dvector& consts);
dvar_vector logPDF_chisquare(dvar_vector& x,dvar_vector params,dvector& consts);
/*---------------------------------------------------*
* X~Chisquare(nu) w/ nu dof.                         *
*   params = <empty>                                 *
*   consts = nu (dof)                                *
*---------------------------------------------------*/
dvariable logPDF_chisqdevs(prevariable& x,dvar_vector params, dvector& consts);
/*---------------------------------------------------*
* norm2(x/stdev)~Chisquare(nu)                       *
*   params = stdev                                   *
*   consts = <empty>                                 *
*---------------------------------------------------*/
dvariable logPDF_chisqdevs(dvar_vector& x,dvar_vector params, dvector& consts);
/*---------------------------------------------------*
*   pdf(x) = (1/lambda)*exp(-x/lambda)               *
*   params = lambda (scale)                          *
*   consts = none                                    *
*----------------------------------------------------*/
dvariable   logPDF_exponential(prevariable& x,dvar_vector params,dvector& consts);
dvar_vector logPDF_exponential(dvar_vector& x,dvar_vector params,dvector& consts);
/*----------------------------------------------------------------------*
* name      : log_gamma_density                                         *
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
dvariable   logPDF_gamma(prevariable& x,dvar_vector params,dvector& consts);
dvar_vector logPDF_gamma(dvar_vector& x,dvar_vector params,dvector& consts);

dvariable logPDF_invchisquare       (prevariable& x,dvar_vector params,dvector& consts);
dvariable logPDF_invgamma           (prevariable& x,dvar_vector params,dvector& consts);
dvariable logPDF_invgaussian        (prevariable& x,dvar_vector params,dvector& consts);
dvariable logPDF_lognormal          (prevariable& x,dvar_vector params,dvector& consts);
dvariable logPDF_logscale_normal    (prevariable& x,dvar_vector params,dvector& consts);
dvariable logPDF_normal             (prevariable& x,dvar_vector params,dvector& consts);
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
dvariable logPDF_scaled_invchisquare(prevariable& x,dvar_vector params,dvector& consts);
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
dvariable logPDF_scaledCV_invchisquare(prevariable& cv,dvar_vector params, dvector& consts);
/*------------------------------------------------------*
*   pdf(x) = gamma((vu+1)/2)/(sqrt(nu*pi)*gamma(nu/2))* *
*               [1+x^2/nu]^(nu+1)/2                     *
*   params = nu (dof)                                   *
*   consts = none                                       *
*-------------------------------------------------------*/
dvariable logPDF_t(prevariable& x,dvar_vector params,dvector& consts);
/*-------------------------------------------------------------*
*   pdf(x) = f*(2*pi*sigma^2)^-0.5 * exp(-0.5*((x-mu)/sigma)^2)*
*          on interval {min,max} [f is a coefficient so the    *
*          integral from min to max = 1.                       *
*   params = mu, sigma                                         *
*   consts = min, max                                          *
*NOTE: This gives the same value as logPDF_normal since f is a *
*      constant and immaterial to derivative calcs.            *
*-------------------------------------------------------------*/
dvariable logPDF_truncated_normal(prevariable& x,dvar_vector params,dvector& consts);


    /********************************************
    * name      : samplePDF_?????               *
    * purpose   : draw sample from pdf(x))      *
    *   params: vector of parameters            *
    *   consts: vector of constants             *
    ********************************************/
//    double samplePDF_beta        (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_cauchy             (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_chisquare          (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_exponential        (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_gamma              (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_invchisquare       (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_invgamma           (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_invgaussian        (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_lognormal          (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_normal             (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_scaled_invchisquare(random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_scaledCV_invchisquare(random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_t                  (random_number_generator& rng,dvector& params,dvector& consts);
    double samplePDF_truncated_normal   (random_number_generator& rng,dvector& params,dvector& consts);

    //delta function    
    int deltafcn(int i, int j);

    /***********************************************************
     * Test if at least one of a set of parameters is          *
     * being estimated in the current phase.                   *
    ***********************************************************/
    bool isActive(param_init_number& p);
    bool isActive(param_init_vector& p);
    bool isActive(param_init_number_vector& p);
    bool isActive(param_init_vector_vector& p);
    
} //namespace wts
#endif
