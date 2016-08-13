//-------------------------------------------------------------------------------------
//Tanner Crab Projection Model
//  Developed by 
//      2012-09: Jack Turnock
//  Completely rewritten by
//      2013-07: William Stockhausen
//
// Input recruitments, n-at-size are in THOUSANDS of crab. Converted to MILLIONS of crab.
// Input spawning biomass in THOUSANDS t.
// Input weight-at-size in t/indiv. Converted to kg/indiv
//
//uses mature male biomass as spawning currency.
//
// Output recruitment and other number quantities in MILLIONS of crab.
// Output MMB and other biomass quantities        in THOUSANDS t.
//
//  20160813: 1. Put under version control at github (wStockhausen/TCSAM2013ProjMod)
//            2. Renamed tpl to TCSAM2013ProjMod.tpl.
//            3. moved associated header, source files to "include", "src" sub-folders.
//            4. Renamed echo file to TCSAM2013ProjMod.chk
//            5. Input recruitment units REMAIN 1000's of crabs (in accordance 
//                  with TCSAM2013 output file for projection model)
//            6. Initial test with 2015 input file worked fine.

//-------------------------------------------------------------
//-------------------------------------------------------------
GLOBALS_SECTION
    #include <admodel.h>
    #include "ModelConstants.hpp"
    #include "admbFunctions.hpp"
    #include "rFunctions.hpp"
    #include "StockRecruitFunctions.hpp"
    
    ofstream echo("TCSAM2013ProjMod.chk");
    ofstream mcmc;     //stream for mcmc output

    
    int on = 0;
    int flg = 0;
    int resp = 0;
    double pi = acos(-1.0);
    double smallVal = 0.00001;//small value
    
    int iSeed = 0;        //random number generator seed
    random_number_generator rng(111);
    
    int recLag = 4;    //assumed lag from spawning to recruitment
    
    const int nXs  = 2;//number of sexes            (FEMALE    = 1, MALE      = 2; defined in ModelConstants.cpp)
    const int nSCs = 2;//number of shell conditions (NEW_SHELL = 1, OLD_SHELL = 2; defined in ModelConstants.cpp)    
    const int nMSs = 2;//number of maturity states  (IMMATURE  = 1, MATURE    = 2; defined in ModelConstants.cpp)
    
    const int idxYPR = 9;//row index into dmatrix res = calcSBPRF(...) for yield-per-recruit (based on retained catch in directed Tanner crab fishery)
    
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
DATA_SECTION

 LOCAL_CALCS    
    //resample
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-rnSeed"))>-1) {
        if (on+1<argc) {
            iSeed=atoi(ad_comm::argv[on+1]);
        } else {
            cout<<"-------------------------------------------"<<endl;
            cout<<"Enter iSeed for random number generator: ";
            cin>>iSeed;
        }
        rng.reinitialize(iSeed);
        echo<<"#random number generator seed = "<<iSeed<<endl;
        echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
 END_CALCS   
 
    //data inputs
    init_int mMnYr  //assessment model start year
    init_int mMxYr  //assessment model end year
    int styr   //start year of projections
    !!styr = mMxYr;
    init_int nYrs   //number of years for projections
    init_int nSims  //number of simulations when using random recruitment   
    init_int nZs    //number of size bins in the model
    !!echo<<"mMnYr  = "<<mMnYr<<endl;
    !!echo<<"mMxYr  = "<<mMxYr<<endl;
    !!echo<<"styr  = "<<styr<<endl;
    !!echo<<"nYrs  = "<<nYrs<<endl;    
    !!echo<<"nSims = "<<nSims<<endl;
    !!echo<<"nZs   = "<<nZs<<endl;
    
    init_number tgtSBPR       //target fraction (XX) for spawning biomass per recruit, e.g. to get F50% this is 0.5
    init_number avgRecForBXX  //total (males+females) average recruitment for BXX/Bmsy calculation   (input in THOUSANDS of recruits)
    init_number inpCurrSpB;   //current spawning biomass from assessment model (THOUSANDS t)
    init_number bioCV;        //uncertainty in initial biomass for ABC calculation
    init_int inpHrvStrat      //harvest strategy to run      <-wts: make this a command line input?
    init_number ctrlAlpha     //control rule parameter: F/Fmsy = (B/Bmsy-alpha)/(1-alpha) for beta  <  B/Bmsy <=1  
    init_number ctrlBeta      //control rule parameter: F      = 0              for B/Bmsy<= beta 
    !!echo<<"tgtSBPR      = "<<tgtSBPR<<endl;
    !!echo<<"inpHrvStrat  = "<<inpHrvStrat<<endl;
    !!echo<<"avgRecForBXX = "<<avgRecForBXX<<endl;
    !!echo<<"inpCurrSpB   = "<<inpCurrSpB<<endl;
    !!echo<<"bioCV        = "<<bioCV<<endl;
    !!echo<<"ctrlAlpha    = "<<ctrlAlpha<<endl;
    !!echo<<"ctrlBeta     = "<<ctrlBeta<<endl;
    
    //convert recruit units from THOUSANDS to MILLIONS of indivs
    !!avgRecForBXX = avgRecForBXX/1000.0;
    
    //-----------------------------------------------------------------------------------------
    // stock-recruit relationship info
    //-----------------------------------------------------------------------------------------
    init_int srType         //flag for type of spawner-recruit curve: bev holt >0 or Ricker curve <0
    init_number recGamma    //autocorrelation parameter for recruiment variability
    init_number recDep      //SR depensation parameter (not used if < 0)
    init_number inpRecH     //steepness parameter
    init_number inpRecR0    //recruitment in THOUSANDS of individuals
    !!echo<<"srType      = "<<srType<<endl;      //flag for type of spawner-recruit curve: bev holt >0 or Ricker curve <0
    !!echo<<"recGamma    = "<<recGamma<<endl;    //autocorrelation parameter for recruiment variability
    !!echo<<"recDep      = "<<recDep<<endl;      //SR depensation parameter (not used if < 0)
    !!echo<<"inpRecH     = "<<inpRecH<<endl;     //steepness parameter
    !!echo<<"inpRecR0    = "<<inpRecR0<<endl;    //recruitment in 1000's of individuals
    
    //scale recruits from THOUSANDS to MILLIONS of indivs
    !!inpRecR0 = inpRecR0/1000.0;
    
    //-----------------------------------------------------------------------------------------
    //recruitment and spawning biomass time series
    //-----------------------------------------------------------------------------------------
    init_int mMnYrForRecAvg     //min assessment model year for averaging recruitment
    init_int mMxYrForRecAvg     //max assessment model year for averaging recruitment
    !!echo<<"mMnYrForRecAvg = "<<mMnYrForRecAvg<<endl;
    !!echo<<"mMxYrForRecAvg = "<<mMxYrForRecAvg<<endl;
    init_matrix asmtModRec(1,nXs,mMnYr,mMxYr) //recruitments from ass. model start year to endyr (1000's)
    init_vector asmtModSpB(mMnYr,mMxYr-1)     //male spawning biomass at matetime for ass. model strt year to endyr-1 for spawner recruit curve (1000's t)
    !!echo<<"asmtModRec(1,nXs,mMnYr,mMxYr) = "<<endl<<asmtModRec<<endl;
    !!echo<<"asmtModSpB(mMnYr,mMxYr-1)     = "<<endl<<asmtModSpB<<endl;
    
    //scale recruits from THOUSANDS to MILLIONS of indivs
    !!asmtModRec /= 1000.0;

    //-----------------------------------------------------------------------------------------
    // pop info in last year fo assessment model
    //-----------------------------------------------------------------------------------------        
    init_matrix natlength_inew_i(1,nXs,1,nZs); //1000's of indivs
    //!!echo<<natlength_inew_i<<endl;
    init_matrix natlength_iold_i(1,nXs,1,nZs); //1000's
    //!!echo<<natlength_iold_i<<endl;
    init_matrix natlength_mnew_i(1,nXs,1,nZs); //1000's
    //!!echo<<natlength_mnew_i<<endl;
    init_matrix natlength_mold_i(1,nXs,1,nZs); //1000's
    4darray inpNatZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);   //on July 1, year 1; 1's
 LOCAL_CALCS
    //convert to indiv's from  THOUSANDS to MILLIONS of indivs
    for (int x=1;x<=nXs;x++){
        inpNatZ_xsmz(x,NEW_SHELL,IMMATURE) = natlength_inew_i(x)/1000.0;
        inpNatZ_xsmz(x,OLD_SHELL,IMMATURE) = natlength_iold_i(x)/1000.0;
        inpNatZ_xsmz(x,NEW_SHELL,  MATURE) = natlength_mnew_i(x)/1000.0;
        inpNatZ_xsmz(x,OLD_SHELL,  MATURE) = natlength_mold_i(x)/1000.0;
    }
 END_CALCS

    //-----------------------------------------------------------------------------------------
    //fisheries info    
    //-----------------------------------------------------------------------------------------
    init_number midptFishery          //time of fishery relative to survey as fraction of year
    !!echo<<"midptFishery = "<<endl<<midptFishery<<endl;
    
    init_number inpFmTCF   //F for directed Tanner crab fishing mortality
    init_number inpFmSCF   //F for male and female snow fishing bycatch
    init_number inpFmRKF   //F for male and female red king fishing bycatch
    init_number inpFmGTF   //F for male and female trawl fishing bycatch
    !!echo<<"inpFmTCF = "<<endl<<inpFmTCF<<endl;
    !!echo<<"inpFmSCF = "<<endl<<inpFmSCF<<endl;
    !!echo<<"inpFmRKF = "<<endl<<inpFmRKF<<endl;
    !!echo<<"inpFmGTF = "<<endl<<inpFmGTF<<endl;
    
    init_matrix selTCF_TotMale(1,nSCs,1,nZs)     //average of last 4 years sel total male new old shell
    init_matrix selTCF_RetMale(1,nSCs,1,nZs)     //average of last 4 years sel retained curve male new old shell
    init_matrix selTCF_TotMaleEast(1,nSCs,1,nZs) //east total - set same as avg total
    init_matrix selTCF_RetMaleEast(1,nSCs,1,nZs) //east retained - set same as avg retained
    init_matrix selTCF_TotMaleWest(1,nSCs,1,nZs) //west total - set same as avg total
    init_matrix selTCF_RetMaleWest(1,nSCs,1,nZs) //west retained - set relative to total sel. - lower end moved 10mm, rest set between east retained and total retained with lower 10mm less than east retained
    init_vector selTCF_Female(1,nZs)                   //selectivity in directed pot fishery for females
    init_matrix selSCF(1,nXs,1,nZs)                    //selectivity in snow crab fishery females/males
    init_matrix selRKF(1,nXs,1,nZs)                    //selectivity in red king crab fishery females/males
    init_matrix selGTF(1,nXs,1,nZs)                    //selectivity in groundfish trawl fishery female/male
    !!echo<<"selTCF_TotMale(1,nSCs,1,nZs) = "<<endl<<selTCF_TotMale<<endl;     //average of last 4 years sel total male new old shell
    !!echo<<"selTCF_RetMale(1,nSCs,1,nZs) = "<<endl<<selTCF_RetMale<<endl;     //average of last 4 years sel retained curve male new old shell
    !!echo<<"selTCF_TotMaleEast(1,nSCs,1,nZs) = "<<endl<<selTCF_TotMaleEast<<endl; //east total - set same as avg total
    !!echo<<"selTCF_RetMaleEast(1,nSCs,1,nZs) = "<<endl<<selTCF_RetMaleEast<<endl; //east retained - set same as avg retained
    !!echo<<"selTCF_TotMaleWest(1,nSCs,1,nZs) = "<<endl<<selTCF_TotMaleWest<<endl; //west total - set same as avg total
    !!echo<<"selTCF_RetMaleWest(1,nSCs,1,nZs) = "<<endl<<selTCF_RetMaleWest<<endl; //west retained - set relative to total sel. - lower end moved 10mm, rest set between east retained and total retained with lower 10mm less than east retained
    !!echo<<"selTCF_Female(1,nZs) = "<<endl<<selTCF_Female<<endl;                  //selectivity in directed pot fishery for females
    !!echo<<"selSCF(1,nXs,1,nZs) = "<<endl<<selSCF<<endl;                          //selectivity in snow crab fishery females/males
    !!echo<<"selRKF(1,nXs,1,nZs) = "<<endl<<selRKF<<endl;                          //selectivity in red king crab fishery females/males
    !!echo<<"selGTF(1,nXs,1,nZs) = "<<endl<<selGTF<<endl;                          //selectivity in groundfish trawl fishery female/male
    
    //-----------------------------------------------------------------------------------------
    //biological "constants"  
    //-----------------------------------------------------------------------------------------
        //natural mortality  
    init_matrix M_f(1,nSCs,1,nMSs)
    init_matrix M_m(1,nSCs,1,nMSs)
    3darray M_xsm(1,nXs,1,nSCs,1,nMSs);//natural mortality
    !!M_xsm(FEMALE) = M_f;
    !!M_xsm(  MALE) = M_m;
    !!echo<<"M(females):"<<endl;
    !!echo<<M_xsm(FEMALE)<<endl;
    !!echo<<"M(males):"<<endl;
    !!echo<<M_xsm(MALE)<<endl;    
        //weight-at-size
    init_vector wtjf(1,nZs)          //immature female weights-at-size         (t/indiv)
    init_matrix wt(1,nXs,1,nZs)      //mature female, all male weights-at-size (t/indiv)
    !!echo<<"wtjf = "<<endl<<wtjf<<endl;
    !!echo<<"wt   = "<<endl<<wt<<endl;
    
    //convert weights/indiv from t/indiv to kg/indiv.
    //note that if w is in kg/indiv and n in MILLIONS indivs, then b = w*n in THOUSANDS t
    !!wtjf *= 1000.;
    !!wt   *= 1000.;
    
        //size transition matrix
    init_3darray tmZtoZ_xzz(1,nXs,1,nZs,1,nZs) //length to length transition matrix
    !!echo<<"tmZtoZ_xzz = "<<endl<<tmZtoZ_xzz<<endl;
    
    init_matrix prMatNS(1,nXs,1,nZs)  //maturity curve for new shell inidividuals
    init_matrix prMoltImm(1,nXs,1,nZs)//molting probability for immature individuals
    init_matrix prMoltMat(1,nXs,1,nZs)//molting probability for mature individuals
    init_number recPropAsMale         //proportion recruiting as male
    init_number recPropAsNewShell     //proportion of recruiting animals as new shell individuals
    init_vector recPropAtZ(1,nZs)     //numbers recruiting at size (assumed same by sex)
    !!echo<<"prMatNS(1,nXs,1,nZs)   = "<<endl<<prMatNS<<endl;
    !!echo<<"prMoltImm(1,nXs,1,nZs) = "<<endl<<prMoltImm<<endl;
    !!echo<<"prMoltMat(1,nXs,1,nZs) = "<<endl<<prMoltMat<<endl;
    !!echo<<"recPropAsMale     = "<<endl<<recPropAsMale<<endl;
    !!echo<<"recPropAsNewShell = "<<endl<<recPropAsNewShell<<endl;
    !!echo<<"recPropAtZ        = "<<endl<<recPropAtZ<<endl;
    
    init_vector propEast(1,nZs) //proportion of population east of 166   <-wts: should be generalized to sex/shell condition/maturity state
    !!echo<<"propEast = "<<propEast<<endl;
    
    vector propWest(1,nZs)      //proportion of population west of 166
 LOCAL_CALCS
    propWest = 1.0-propEast;
    echo<<"propWest = "<<propWest<<endl;
 END_CALCS   
    
    //these should probably be inputs (they're fixed below)
    number sigma_r                 //ln-scale std dev for recruitment 
    number avgRecBeta              //prior mean for pRecBeta 
    number sigRecBeta              //prior standard deviation for pRecBeta
    
    number medAsmtModSpB                //median male spawning biomass (THOUSANDS t)
    
//-------------------------------------------------------------------------------------
//--new quantities for population dynamics---------------------------------------------
//-------------------------------------------------------------------------------------   
    //spawning biomass metrics
    number currSpB;       //"current" spawning biomass for control rule evaluation                                       (THOUSANDS t)
    vector sbMT(1,nXs);   //spawning biomass at mate time                                                                (THOUSANDS t)
    vector sbpr0(1,nXs);  //spawning biomass per recruit for virgin stock on July 1                                      (kg/indiv or THOUSANDS t/MILLION indivs)
    vector sbpr0FT(1,nXs);//spawning biomass per recruit for virgin stock just prior to fishery (if it occurred)         (kg/indiv or THOUSANDS t/MILLION indivs)
    vector sbpr0MT(1,nXs);//spawning biomass per recruit for virgin stock at mating time (after fishery, if it occurred) (kg/indiv or THOUSANDS t/MILLION indivs)
    vector sbprf(1,nXs);  //spawning biomass per recruit at "f" on July 1                                                (kg/indiv or THOUSANDS t/MILLION indivs)
    vector sbprfFT(1,nXs);//spawning biomass per recruit at "f" just prior to fishery                                    (kg/indiv or THOUSANDS t/MILLION indivs)
    vector sbprfMT(1,nXs);//spawning biomass per recruit at "f" at mating time (after fishery)                           (kg/indiv or THOUSANDS t/MILLION indivs)
    number totB;          //total biomass on July 1                                                                      (kg/indiv or THOUSANDS t/MILLION indivs)
    
    //numbers at size
    3darray rec_xsz(1,nXs,1,nSCs,1,nZs);            //recruitment              (MILLION indivs)
    4darray nAtZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);   //on July 1                (MILLION indivs)
    4darray nAtZBF_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //just prior to fishery    (MILLION indivs)
    4darray nAtZAF_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //just after to fishery    (MILLION indivs)
    4darray nAtZMT_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //numbers at mating time (mateTime = time AFTER fishery occurs that mating occurs) (MILLION indivs)
    4darray nAtZ03_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //only those that will grow  (MILLION indivs)
    4darray nAtZ04_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //all after growth/no growth (MILLION indivs)
    4darray nAtZ05_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //after maturity transition  (MILLION indivs)
        
    4darray prvNatZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);//on July 1 from PREVIOUS iteration in while loop below (MILLION indivs)
    
    //fishing mortality and survival
    number fmTCF; //fully selected fishing mortality in tanner crab fishery
    number fmSCF; //fully selected fishing mortality in snow crab fishery
    number fmRKF; //fully selected fishing mortality in red king crab fishery
    number fmGTF; //fully selected fishing mortality in groundfish trawl fishery    
        //in directed Tanner crab fishery
    4darray tcTotF_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);     //directed Tanner crab fishery mortality on all crabs
    4darray tcTotF_East_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);//directed Tanner crab fishery mortality on all crab in East Region
    4darray tcTotF_West_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);//directed Tanner crab fishery mortality on all crab in West Region
    4darray tcRetF_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);     //directed Tanner crab fishery mortality on all retained crab        
    4darray tcRetF_East_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);//directed Tanner crab fishery mortality on retained crab in East Region
    4darray tcRetF_West_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);//directed Tanner crab fishery mortality on retained crab in West Region
        //in other fisheries
    4darray scF_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //snow crab bycatch (discard) fishing mortality
    4darray rkF_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //red king crab bycatch (discard) fishing mortality
    4darray gtF_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //groundfish trawl bycatch (discard) fishing mortality
        //total fishing mortality and survival
    4darray totF_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //total fishing mortality
    4darray totS_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //survival in fishery
       
    //arrays for catch calculations
    4darray tcRetCatZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);//retained catch-at-size in directed Tanner crab fishery     (MILLION indivs)
    4darray tcTotCatZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);//total catch-at-size in directed Tanner crab fishery        (MILLION indivs)
    4darray scCatZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);   //catch-at-size in snow crab crab fishery                    (MILLION indivs)
    4darray rkCatZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);   //catch-at-size in red king crab fishery                     (MILLION indivs)
    4darray gtCatZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);   //catch-at-size in groundfish trawl crab fishery             (MILLION indivs)
        
    //quantities for SBPR calculations  
    vector recPropBySex(1,nXs);                             //proportion recruiting by sex   
    vector recPropByShell(1,nSCs);                          //proportion recruiting by shell condition
    3darray wtp_xmz(1,nXs,1,nMSs,1,nZs);                    //weight at size by sex, maturity state        (kg/indiv or THOUSANDS t/MILLIONS indivs)
    5darray tmMolt_xmssz(1,nXs,1,nMSs,1,nSCs,1,nSCs,1,nZs); //transition probabilities among shell conditions based on probability crab WILL molt (and change shell condition)
    5darray tmMat_xsmmz(1,nXs,1,nSCs,1,nMSs,1,nMSs,1,nZs);  //transition probability from maturity state to maturity state
    
    number resAvgRec //average recruitment from inputs (MILLIONS indivs)
    
    //quantities for simulations
    number lnVarRec; //ln-scale recruitment variance   WTS: need to calc this!! what about sigma_r above? seems to be same thing
    
    //biological reference points
    vector resSBPR0;         //sbpr(FEMALE,MALE) at F=0  (kg/indiv or THOUSANDS t/MILLION recruits)
    number phi0              //equilibrium MMBPR at F=0  (kg/indiv or THOUSANDS t/MILLION recruits)
    
    matrix resSBPRFXX;       //sbpr info at FXX
    number FXX;              //retained fishing mortality that results in sbprF/sbpr0 = XX (XX=0.35, typically)
    number BXX;              //spawning biomass (MMB at mating time) at FXX using avgRecForBXX*XX*phi0
    number retYXX;           //retained yield at FXX using avgRecForBXX
    number totYXX;           //total    yield at FXX using avgRecForBXX
    
    vector resMSY            //msy info from SR calculation (note it is not allocated here)
    number Fmsy;             //Fmsy (may be FXX or from msy calc.)
    number Bmsy;             //Bmsy (may be BXX or from msy calc.)
    number retMSY;           //retained yield MSY  (may be totYXX or from msy calc.) (THOUSANDS t)
    number totMSY;           //total    yield MSY  (may be totYXX or from msy calc.) (THOUSANDS t)
    
    //harvest control rule outputs
    number Fofl; //Fofl based on control rule
    number OFL;  //OFL based on control rule   (THOUSANDS t)
    number ghl;  //SOA guideline harvest level
    
    //simulation output quantities
    matrix simsMMB(1,nSims,1-recLag,nYrs); //mature   male biomass (at mating time) (THOUSANDS t)
    matrix simsMFB(1,nSims,1,nYrs); //mature female biomass (at mating time)        (THOUSANDS t) 
    matrix simsRecDevs(1,nSims,1,nYrs);    //ln-scale recruit devs
    matrix simsRecruits(1,nSims,1,nYrs);   //simsRecruits                           (MILLIONS)
    matrix simsFofl(1,nSims,1,nYrs);       //retained Fofl
    matrix simsOFL(1,nSims,1,nYrs);        //retained OFL                           (THOUSANDS t)
    matrix simsGHL(1,nSims,1,nYrs);        //GHL                                    (THOUSANDS t)
    matrix simsRetF(1,nSims,1,nYrs);       //actual retained F
    matrix simsTotF(1,nSims,1,nYrs);       //actual total F
    matrix simsRetY(1,nSims,1,nYrs);       //retained yield                         (THOUSANDS t)
    matrix simsTotY(1,nSims,1,nYrs);       //total yield (dead  crab)               (THOUSANDS t)
    matrix simsTotB(1,nSims,1,nYrs);       //total biomass (July 1)                 (THOUSANDS t)
    
    
 LOCAL_CALCS    
    cout<<"finished DATA_SECTION"<<endl; 
//     cout<<"Enter 1 and hit return to continue >> ";
//     cin>>resp;
//     if (resp<1) exit(-1);
 END_CALCS

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
INITIALIZATION_SECTION
    pLnR0 2.
    pRecBeta -0.35
    pLnRecSigma 0.01

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
PARAMETER_SECTION
//parameters to be estimated are all ones that begin with init_ and have a positive
//need to have parameters in log space or logit (pRecBeta) to get things to converge

    init_number pLnR0(1)
    init_bounded_number pRecBeta(-20.,20.,2)
    init_bounded_number pLnRecSigma(-20.,20.,2)

    number recVar
    number sigmasq
    vector predRec(mMnYrForRecAvg,mMxYrForRecAvg)    
    
    number recR0
    number recH
    
    sdreport_number sdDummy;//need to declare this for mcmc calculations
    
    objective_function_value objFun
    
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
PRELIMINARY_CALCS_SECTION
    cout<<"Starting PRELIMINARY_CALCS_SECTION"<<endl;
    int nobs;
    
    fmSCF = inpFmSCF;
    fmRKF = inpFmRKF;
    fmGTF = inpFmGTF;         //F for male and female trawl fishing bycatch
    
    //calc average total recruits over specified time frame
    resAvgRec = sum(asmtModRec(FEMALE)(mMnYrForRecAvg,mMxYrForRecAvg)+asmtModRec(MALE)(mMnYrForRecAvg,mMxYrForRecAvg))/(mMxYrForRecAvg-mMnYrForRecAvg+1);
    cout<<"resAvgRec (1) (avg over "<<mMnYrForRecAvg<<":"<<mMxYrForRecAvg<<")= "<<resAvgRec<<endl;
    if (avgRecForBXX<0) avgRecForBXX = resAvgRec;//use average total for BXX calc.s, if not provided in data file
    
    //use average of simsRecruits where the effective sp biomass is > the median effective sp biomass
    resAvgRec = 0.0;
    int nrecs = 0;
    dvector w=sort(asmtModSpB(mMnYrForRecAvg-recLag,mMxYrForRecAvg-recLag));//wts:added the limited time frame here
    cout<<"sorted spawning biomass from assessment model = "<<endl<<tb<<w<<endl;
    //get the median asmtModSpB
    int nYrsAsmtMod = mMxYrForRecAvg-mMnYrForRecAvg+1;
    w.shift(1);//shift min index to 1
    if(2*(int(nYrsAsmtMod/2))==nYrsAsmtMod){ 
        //if nYrsAsmtMod even average two values
        cout<<" nYrsAsmtMod even "<<endl;
        medAsmtModSpB = (w(nYrsAsmtMod/2)+w((nYrsAsmtMod/2)+1))/2.0;
    } else {
        //if nYrsAsmtMod odd
        cout<<" nYrsAsmtMod odd "<<endl;
        medAsmtModSpB = w((nYrsAsmtMod+1)/2);
    }
    cout<<"medAsmtModSpB = "<<medAsmtModSpB<<endl;
    for(int i=(mMnYrForRecAvg+recLag);i<=mMxYrForRecAvg;i++){
        if(asmtModSpB(i-recLag) >= medAsmtModSpB){
            resAvgRec += asmtModRec(FEMALE,i)+asmtModRec(MALE,i);
            nrecs++;
        }
    }
    resAvgRec = resAvgRec/nrecs;
    cout<<"resAvgRec (2) = "<<resAvgRec<<endl;
    
    // prior mean and sd from Clark's paper pRecBeta = 0.3 (steepness of 0.66) and sigma of 0.6
    //red king crab spawner recruit curve has tau of about .23 which is steepness of about 0.52
    //which is pRecBeta of -0.4
    //avgRecBeta is prior mean for pRecBeta 
    //sigRecBeta is prior standard deviation
    //using average observed recruitment as prior for R0
    // no prior needed for variance of recruitment
    //  avgRecBeta = -0.6;    // this gives steepness of about 0.48
    //  avgRecBeta =  0.3;    //this is steepness of 0.66
    avgRecBeta = -0.4;     //this is steepness of 0.52
    sigRecBeta = 0.3;
    sigma_r    = 0.6;
    
    recR0 = mfexp(pLnR0);  //in MILLIONS indivs
    recH  = 0.8*(mfexp(pRecBeta)/(1+mfexp(pRecBeta)))+0.2;
    recVar = square(mfexp(pLnRecSigma));
    cout<<"Initial values"<<endl;
    cout<<"recR0  = "<<recR0<<endl;
    cout<<"recH   = "<<recH<<endl;
    cout<<"recVar = "<<recVar<<endl;
    
    /////////////////////////////////////////////////////////////
    cout<<endl<<"------------------------------------------------------"<<endl;
    cout<<"-------SPBR0 calculations (kg/indiv or THOUSANDS t/MILLION recruits)----------"<<endl;
    setBioQuants();        //initialize quantities for sbpr calculations
    resSBPR0 = calcSBPR0();//calculate sbpr (kg/recruit or THOUSANDS T/MILLIONS recruits) for virgin population
    cout<<"SBPR0 = "<<resSBPR0<<endl;
    phi0 = resSBPR0(MALE); //MMB/R at mating time for virgin stock (kg/recruit or THOUSANDS T/MILLIONS recruits)
    
    /////////////////////////////////////////////////////////////
    //Solve for F35%, SBPRF35
    cout<<endl<<"------------------------------------------------------"<<endl;
    cout<<"-------SPBRFXX% calculations (kg/indiv or THOUSANDS t/MILLION recruits)-----------"<<endl;
    {
        int optF_EW = 1;                        //use original approach for merging assumed E/W differences in directed fishing mortality
        FXX = findFxx(tgtSBPR,phi0,optF_EW);    //calculate F35%
        BXX = avgRecForBXX*(tgtSBPR*phi0);      //assumes adequate convegence so phiXX = tgtSBPR*phi0.
        
        int optCalcYPR = 1;                            //turn off ypr calculation
        resSBPRFXX = calcSBPRF(FXX,optCalcYPR,optF_EW);//calculate sbpr, ypr at F35%
        
        retYXX = avgRecForBXX*sum(resSBPRFXX( 9));//retained equilibrium yield at FXX.
        totYXX = avgRecForBXX*sum(resSBPRFXX(10));//   total equilibrium yield at FXX.
        
        cout<<tgtSBPR                <<tb<<tb<<"#target SBPR ratio"<<endl;
        cout<<resSBPRFXX(1,MALE)/phi0<<tb<<tb<<"#resulting ratio"<<endl;
        cout<<FXX                    <<tb<<tb<<"#FXX"<<endl;
        cout<<avgRecForBXX           <<tb<<tb<<"#avg. rec used for BXX (MILLIONS)"<<endl;
        cout<<inpCurrSpB             <<tb<<tb<<"#current spawning biomass (1000's t)"<<endl;
        cout<<phi0                   <<tb<<tb<<"#phi0   (KG/INDIV OR  (THOUSANDS T)/(MILLIONS recruits)"<<endl;
        cout<<BXX                    <<tb<<tb<<"#BXX    (THOUSANDS T)"<<endl;
        cout<<retYXX                 <<tb<<tb<<"#retYXX (THOUSANDS T)"<<endl;
        cout<<totYXX                 <<tb<<tb<<"#totYXX (THOUSANDS T)"<<endl;
        writeSBPRF(cout,resSBPRFXX);    
            
        ofstream post("calcs.SBPRF35.csv");
        post<<"avg. rec"                 <<cc<<avgRecForBXX<<endl;
        post<<"current B"                <<cc<<inpCurrSpB<<cc<<"1000's t"<<endl;
        post<<"target SBPR ratio"        <<cc<<tgtSBPR<<endl;
        post<<"resulting ratio"          <<cc<<resSBPRFXX(1,MALE)/phi0<<endl;
        post<<"FXX"                      <<cc<<FXX<<endl;
        post<<"SpB/R(F=0)"               <<cc<<phi0<<cc<<"kg/recruit"<<endl;
        post<<"BXX"                      <<cc<<BXX<<cc<<"1000's t"<<endl;
        post<<"retained catch biomass"   <<cc<<retYXX<<cc<<"1000's t"<<endl;
        post<<"total mortality (biomass)"<<cc<<totYXX<<cc<<"1000's t"<<endl;
        writeSBPRF(post,resSBPRFXX);        
        post.close();
    }
    /////////////////////////////////////////////////////////////    
    cout<<"Enter 1 and hit return to continue >> ";
    cin>>resp;
    if (resp<1) exit(-1);
    ///////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////
    //Calculate Tier 3 Fofl, OFL and ABC
    cout<<"#-------------------------------------------------------"<<endl;
    cout<<"#-------------------------------------------------------"<<endl;
    cout<<"#--Calculating Tier 3 Fofl, OFL and ABC based on avgRec, currSpB, bioCV, and initial n-at-size"<<endl;
    calcOFLandABC();
    /////////////////////////////////////////////////////////////    
    cout<<"Enter 1 and hit return to continue >> ";
    cin>>resp;
    if (resp<1) exit(-1);
    ///////////////////////////////////////////////////////////
    
    //calculate MSY quantities for input recR0 and H
    {
        // R0 for deterministic Fmsy calc
        recR0 = inpRecR0;                  //wts: inpRecR0 should NOT be on ln-scale here
        recH  = inpRecH;
        
        //calc msy
        cout<<endl<<endl<<"Calculating Fmsy, Bmsy, MSY----------------------------"<<endl;
        int optTCF_EW = 1;//Jack's way to combine E/W fishing mortalities
        resMSY = calcMSY(value(recR0), value(recH), phi0, optTCF_EW, srType);
        cout<<"# Fmsy    msy    Bmsy    eqR       phi      xx"<<endl;
        cout<<resMSY<<tb<<resMSY(5)/phi0<<endl;
        ofstream post("calcs.MSY.dat");
        post<<"# Fmsy    msy    Bmsy    eqR       phi      xx"<<endl;
        post<<resMSY<<tb<<resMSY(5)/phi0<<endl;
        post.close();
    }
    /////////////////////////////////////////////////////////////    
    cout<<"Enter 1 and hit return to continue >> ";
    cin>>resp;
    if (resp<1) exit(-1);
    ///////////////////////////////////////////////////////////
    
    //calc recH such that Fmsy = FXX%
    {
        cout<<endl<<endl<<"Calculating recH such that Fmsy = F35%----------------------------"<<endl;
        int optTCF_EW = 1;//Jack's way to combine E/W fishing mortalities
        
//         dmatrix resRecHXX = calcHs(FXX, value(recR0), phi0, optTCF_EW, srType);        
//         writeRecHSearch(cout, resRecHXX);
        
        double recH35XX = findH(FXX,value(recR0), phi0, optTCF_EW, srType);
        cout<<recH35XX<<tb<<"#recH yielding Fmsy = FXX%"<<endl;
        
        ofstream post("calcs.findH.dat");
        post<<recH35XX<<tb<<"#recH yielding Fmsy = FXX%"<<endl;
//        writeRecHSearch(post, resRecHXX);
        post.close();
    }
    /////////////////////////////////////////////////////////////    
    cout<<"Enter 1 and hit return to continue >> ";
    cin>>resp;
    if (resp<1) exit(-1);
    /////////////////////////////////////////////////////////////    
    
    //run simulations with input values for recR0 and recH
    {
        cout<<"Running simulations!"<<endl;
        doSimulations(nSims,nYrs,inpHrvStrat);
        //then what?
        cout<<"Finished running simulations!"<<endl;
    }    
    /////////////////////////////////////////////////////////////    
    cout<<"Enter 1 and hit return to continue >> ";
    cin>>resp;
    if (resp<1) exit(-1);
    /////////////////////////////////////////////////////////////    
    
    //if running mcmc evaluations, open mcmc file
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1) {
        openMCMCFile();
        cout<<"MCEVAL is on"<<endl;
    }
    
    cout<<"Finished PRELIMINARY_CALCS_SECTION"<<endl;
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
PROCEDURE_SECTION
    recR0  = mfexp(pLnR0);
    recH   = 0.8*(mfexp(pRecBeta)/(1.0+mfexp(pRecBeta)))+0.2;
    recVar = square(mfexp(pLnRecSigma));
    
    predRecruitment();
    
    evaluate_the_objective_function();
    
    if (mceval_phase()){
        doSimulations(1,nYrs,inpHrvStrat);//only run 1 simulation for each mceval
        writeMCMCtoR(mcmc);
    }

//-------------------------------------------------------------------------------------
FUNCTION predRecruitment
    int debug = 0;
    if (debug) cout<<"starting predRecruitment"<<endl;

    //four year lag time                         <-wts: why 4 year time lag and not 5?
    for(int i=mMnYrForRecAvg;i<=mMxYrForRecAvg;i++){
        predRec(i) = calcSRFunction(recR0,recH,phi0,asmtModSpB(i-recLag),srType);//predRec is total recruitment
    }
    
    if (debug) cout<<"finished predRecruitment"<<endl;
    
//-------------------------------------------------------------------------------------
FUNCTION evaluate_the_objective_function
    int debug = 0;
    if (debug) cout<<"starting evaluate_the_objective_function"<<endl;
    
    objFun.initialize();
    sigmasq.initialize();
    for(int i=mMnYrForRecAvg;i<=mMxYrForRecAvg;i++) {
        sigmasq += square(log((asmtModRec(FEMALE,i)+asmtModRec(MALE,i)))-log(predRec(i))+recVar/2.);//predRec is total recruitment
        // cout<<" sigmasq = "<<sigmasq<<endl;
    }
    objFun += sigmasq/(2.*recVar);
    //  cout<<" objFun 1 "<<objFun<<endl;
    //Post convergence estimate of residual variance
    sigmasq = sigmasq/double(mMxYrForRecAvg-mMnYrForRecAvg+1);
    
    objFun += double(mMxYrForRecAvg-mMnYrForRecAvg+1)*log(sqrt(recVar));
    // cout<<" objFun 2 "<<objFun<<endl;
    
    //variance is cv * prior mean, prior mean is observed recruitment, resAvgRec/2.
    
    objFun += log(sqrt(2.0*pi))+log(sigma_r*log(resAvgRec+1e-06)+1e-06)+ 0.5*square((pLnR0-log(resAvgRec+1e-06))/(sigma_r*log(resAvgRec+1e-06)));
    // cout<<" objFun 3 "<<objFun<<endl;
    
    objFun += log(sqrt(2.0*pi))+2.0*log(sigRecBeta)+square(pRecBeta-avgRecBeta)/(2.0*square(sigRecBeta)); 
    // cout<<" objFun 4 "<<objFun<<endl;
    
    if (debug) cout<<"finished evaluate_the_objective_function"<<endl;

//-------------------------------------------------------------------------------------
//  Calculates dmatrix of yields (not ypr) based on 
//  input directed fishing mortality rate.
//
//  Output:
//      dmatrix res from doFisheries(...,...)
//  Modified:
//      nAtZAF_xsmz
//
//-------------------------------------------------------------------------------------
FUNCTION dmatrix calcYield(double tcF)

    double scF = fmSCF;
    double rkF = fmRKF;
    double gtF = fmGTF;
    int optTCF_EW = 1;
    
    setFishingMortalityRates(tcF, scF, rkF, gtF, optTCF_EW);
    dmatrix res = doFisheries(nAtZBF_xsmz, nAtZAF_xsmz);
        
    return res;

//-------------------------------------------------------------------------------------
//  Find the directed Tanner crab fishing mortality that results in the given yield.
//     yldType =  1: totCat;
//     yldType =  2: totYld;
//     yldType =  3: retCat;
//     yldType =  4: retYld;
//     yldType =  5: tcTotCat;
//     yldType =  6: tcTotYld;
//     yldType =  7: tcRetCat;
//     yldType =  8: tcRetYld;
//     yldType =  9: scTotCat;
//     yldType = 10: scTotYld;
//     yldType = 11: rkTotCat;
//     yldType = 12: rkTotYld;
//     yldType = 13: gtTotCat;
//     yldType = 14: gtTotYld;
//
FUNCTION double findTCFforYield(double yld, double initTCF, int yldType)
    int debug = 1;
    if (debug) cout<<"Starting findTCFforYield("<<yld<<cc<<initTCF<<cc<<yldType<<")"<<endl;

    double dF  = 0.001;
    double tcF = initTCF;
    double tcFp;
    double tcFm;
    double yldp;
    double yld0;
    double yldm;
    double dydf;
    double delF;
    
    //loop to find tcF that gives desired yld
    double cvg     = 1.0;   //convergence indicator
    double convCri = 1.0e-6;//covergence criteria
    int itr = 1;            //iteration counter
    while ((convCri<cvg)&&(itr<=20)){
        yldm  = sum(calcYield(tcF-dF)(yldType));
        yld0  = sum(calcYield(tcF)(yldType));
        yldp  = sum(calcYield(tcF+dF)(yldType));
        dydf  = (yldp-yldm)/(2*dF); //  Newton-Raphson approximation for first derivitive
        delF  = (yld-yld0)/dydf;    //  update to tcF
        cvg   = fabs(delF);        // convergence indicator
        if (debug) cout <<"yld = "<<yld<<" yld0= "<<yld0<<" dydf= "<<dydf<<" cvg = "<<cvg<<" tcF = "<<tcF+delF<<endl;
        tcF += delF; //update tcF
        ++itr;
    }
            
    return tcF;

    if (debug) cout<<"finished findTCFforYield("<<yld<<cc<<initTCF<<cc<<yldType<<")"<<endl;

//-------------------------------------------------------------------------------------
//  Calculate sex-specific biomass corresponding to nAtZ_xsmz 4d array.
FUNCTION dvector calcBiomass(d4_array& n_xsmz)
    dvector b(1,nXs);
    b.initialize();
    for (int x=1;x<=nXs;x++){
        for (int sc=1;sc<=nSCs;sc++){
            for (int ms=1;ms<=nMSs;ms++) b(x) += wtp_xmz(x,ms)*n_xsmz(x,sc,ms);
        }
    }
    return b;

//-------------------------------------------------------------------------------------
FUNCTION calcOFLandABC
    //initialize "current" spawning stock biomass
    currSpB = inpCurrSpB;
    
    //For Tier 3
    double Bmsy = BXX;
    double Fmsy = FXX;//this is directed fishing
    
    //calculate true OFL
    //initialize numbers-at-size
    nAtZ_xsmz.initialize();
    for (int x=1;x<=nXs;x++){
        for (int s=1;s<=nSCs;s++){
            nAtZ_xsmz(x,s) = inpNatZ_xsmz(x,s);
        }
    }
    
    //calculate Fofl based on current nAtZ
    double cFofl = calcFofl(currSpB, Bmsy, Fmsy);//for directed fishing
    
    //Advance pop model to fisheries time
    double tcF = cFofl;
    double scF = fmSCF;
    double rkF = fmRKF;
    double gtF = fmGTF;
    int optTCF_EW = 1;
    advancePopToFisheries(nAtZ_xsmz, nAtZBF_xsmz);
    setFishingMortalityRates(tcF, scF, rkF, gtF, optTCF_EW);
    dmatrix res = doFisheries(nAtZBF_xsmz, nAtZAF_xsmz);
    
    //save results
    double cTotOFL = sum(res(2));//total yield OFL based on assessment model results
    
    //now determine distribution of OFLs based on uncertainty in initial numbers at size
    double bioVar = log(1+bioCV*bioCV);//ln-scale variance in biomass
    
    int nRepsOFL = 10000; //these should be model inputs
    double pstar = 0.49;
    
    dvector estSpB(1,nRepsOFL);
    dvector estFofl(1,nRepsOFL);
    dvector estTotOFL(1,nRepsOFL);
    for (int iRep=1;iRep<=nRepsOFL;iRep++){
        double rfac = mfexp(sqrt(bioVar)*randn(rng)-bioVar/2.0);//random factor for biomass estimation error w/ bias correction so input is mean
        //estimate of current spawning biomass
        estSpB(iRep) = rfac*inpCurrSpB;
        //initialize numbers-at-size
        nAtZ_xsmz.initialize();
        for (int x=1;x<=nXs;x++){
            for (int s=1;s<=nSCs;s++){
                nAtZ_xsmz(x,s) = rfac*inpNatZ_xsmz(x,s);
            }
        }
        //calculate Fofl based on estimated biomass
        estFofl(iRep) = calcFofl(estSpB(iRep), Bmsy, Fmsy);//for directed fishing
        //Advance pop model to fisheries time
        advancePopToFisheries(nAtZ_xsmz, nAtZBF_xsmz);
        //set fishing mortality rates
        setFishingMortalityRates(estFofl(iRep), scF, rkF, gtF, optTCF_EW);
        //do fisheries to calculate OFL (don't need to calculate population beyond this)
        dmatrix res = doFisheries(nAtZBF_xsmz, nAtZAF_xsmz);        
        //save results
        estTotOFL(iRep) = sum(res(2));
    }//iSims loop
    
    dvector sortEstOFL = sort(estTotOFL);
    double medOFL   = sortEstOFL(nRepsOFL/2);
    double pstarOFL = sortEstOFL((int) (pstar*nRepsOFL));
    double by       = medOFL-pstarOFL;
    double totABC   = (1.0-by)*medOFL;
    
    //calculate MMB at mate time corresponding to taking the OFL
    double nxtSpB;
    dmatrix ylds;
    {
        nAtZ_xsmz.initialize();
        for (int x=1;x<=nXs;x++){
            for (int s=1;s<=nSCs;s++){
                nAtZ_xsmz(x,s) = inpNatZ_xsmz(x,s);//initialize numbers-at-size to assessment model
            }
        }                
        advancePopToFisheries(nAtZ_xsmz, nAtZBF_xsmz);//Advance pop model to fisheries time
        //find tcF such that fishing yields median OFL as total catch
        int yldType = 2;//fit to total yield
        tcF = findTCFforYield(medOFL, Fmsy, yldType);
        ylds = calcYield(tcF); //calc pop model through fisheries at fishing mortality that results in median OFL as total catch
        cout<<"tcF = "<<tcF<<endl;
        cout<<"total yield = "<<sum(ylds(yldType))<<endl;
        //Advance pop model from after fisheries to end of year (don't need to add in recruits)
        dvector sbMT = advancePopAfterFisheries(0.0, nAtZAF_xsmz, nAtZ_xsmz);
        cout<<"total mature biomass at mate time = "<<sum(sbMT)<<endl;
        nxtSpB = sbMT(MALE);//MMB at mating time after fishing at median OFL
    }
    
    ofstream osOFL_ABC("resOFL_ABC.dat");
    writeOFL_ABC(osOFL_ABC,bioCV,Fmsy,Bmsy,currSpB,nxtSpB,cFofl,cTotOFL,medOFL,pstarOFL,totABC,ylds,estSpB,estFofl,estTotOFL);
    osOFL_ABC.close();    
    osOFL_ABC.open("resOFL_ABC.R");
    writeToR_OFL_ABC(osOFL_ABC,bioCV,Fmsy,Bmsy,currSpB,nxtSpB,cFofl,cTotOFL,medOFL,pstarOFL,totABC,ylds,estSpB,estFofl,sortEstOFL);
    osOFL_ABC.close();    

//-------------------------------------------------------------------------------------
FUNCTION void writeOFL_ABC(ostream& os, double bioCV, double Fmsy, double Bmsy, double currSpB, double nextSpB, double cFofl, double cOFL, 
                           double medOFL, double pstarOFL, double ABC, dmatrix& ylds, dvector& estSpB, dvector& estFofl, dvector& estOFL)
    os<<"bioCV   = "<<bioCV<<endl;
    os<<"Fmsy    = "<<Fmsy<<endl;
    os<<"Bmsy    = "<<Bmsy<<endl;
    os<<"currSpB = "<<currSpB<<endl;
    os<<"currB/Bmsy  = "<<currSpB/Bmsy<<endl;
    os<<"cFofl   = "<<cFofl<<endl;
    os<<"cOFL    = "<<cOFL<<endl;
    os<<"medOFL  = "<<medOFL<<endl;
    os<<"p* OFL  = "<<pstarOFL<<endl;
    os<<"ABC     = "<<ABC<<endl;
    os<<"nextB   = "<<nextSpB<<endl;
    os<<"nextB/Bmsy = "<<nextSpB/Bmsy<<endl;
    os<<"OFL yields = "<<endl<<ylds<<endl;
    os<<"estSpB  = "<<endl<<estSpB<<endl;
    os<<"estFofl = "<<endl<<estFofl<<endl;
    os<<"estOFL  = "<<endl<<estOFL<<endl;
    
//-------------------------------------------------------------------------------------
FUNCTION void writeToR_OFL_ABC(ostream& os, double bioCV, double Fmsy, double Bmsy, double currSpB, double nextSpB, double cFofl, double cOFL, 
                               double medOFL, double pstarOFL, double ABC, dmatrix& ylds, dvector& estSpB, dvector& estFofl, dvector& estOFL)
    os<<"res=list(";
    os<<"target.xx="<<tgtSBPR<<cc;
    os<<"phi0="<<phi0<<cc;
    os<<"FXX="<<FXX<<cc;
    os<<"avg.rec="<<avgRecForBXX<<cc;
    os<<"BXX="<<BXX<<cc;
    os<<"B.curr="<<inpCurrSpB<<cc;
    os<<"B.next = "<<nextSpB<<cc<<endl;
    os<<"retYXX="<<retYXX<<cc;
    os<<"totYXX="<<totYXX<<cc<<endl;
    os<<"bioCV  = "<<bioCV<<cc;
    os<<"Fmsy   = "<<Fmsy<<cc;
    os<<"Bmsy   = "<<Bmsy<<cc;
    os<<"currMMB    = "<<currSpB<<cc;
    os<<"nextMMB    = "<<nextSpB<<cc;
    os<<"cFofl   = "<<cFofl<<cc;
    os<<"cOFL    = "<<cOFL<<cc;
    os<<"medOFL = "<<medOFL<<cc;
    os<<"pstarOFL = "<<pstarOFL<<cc;
    os<<"ABC    = "<<ABC<<cc<<endl;
    os<<"yields="; R::writeToR(os,ylds); os<<cc<<endl;
    os<<"estSpB="; R::writeToR(os,estSpB); os<<cc<<endl;
    os<<"estFofl="; R::writeToR(os,estFofl); os<<cc<<endl;
    os<<"estOFL="; R::writeToR(os,estOFL); os<<endl;
    os<<")"<<endl;
    
//-------------------------------------------------------------------------------------
FUNCTION void doSimulations(int nSims, int nYrs, int harvStrat)
    int debug = 1;
    if (debug) cout<<"starting doSimulations("<<nSims<<cc<<nYrs<<cc<<harvStrat<<")"<<endl;
    simsMMB.initialize();        //mature   male biomass (at mating time)
    simsMFB.initialize();        //mature female biomass (at mating time)
    simsRecDevs.initialize();    //ln-scale recruit devs
    simsRecruits.initialize();   //simsRecruits
    simsFofl.initialize();       //retained Fofl
    simsOFL.initialize();        //retained OFL
    simsGHL.initialize();        //GHL
    simsRetF.initialize();       //actual retained F
    simsTotF.initialize();       //actual total F
    simsRetY.initialize();       //retained yield
    simsTotY.initialize();       //total yield (dead  crab)
    simsTotB.initialize();       //total biomass (July 1)
    
    //run projections 
    for (int sim=1;sim<=nSims;sim++){
        //run population dynamics model
        runPopDyMod(sim,nYrs,harvStrat);
        //assemble results
        //??
    }

    if (debug) cout<<"finished doSimulations("<<nSims<<cc<<nYrs<<cc<<harvStrat<<")"<<endl;
        
//-------------------------------------------------------------------------------------
FUNCTION void runPopDyMod(int iSim, int nYrs, int harvStrat)
    int debug = 1;
    if (debug) cout<<"starting runPopDyMod("<<iSim<<cc<<nYrs<<cc<<harvStrat<<")"<<endl;
    //initialize population dynamics model
    initPopDyMod(iSim, nYrs, harvStrat);
    for (int yr=1;yr<=nYrs;yr++){
        runPopDyModOneYear(iSim,yr,harvStrat);
    }

    if (debug) cout<<"finished runPopDyMod("<<iSim<<cc<<nYrs<<cc<<harvStrat<<")"<<endl;
    
//-------------------------------------------------------------------------------------
FUNCTION void initPopDyMod(int iSim, int nYrs, int harvStrat)
    int debug = 1;
    if (debug) cout<<"starting initPopDyMod("<<iSim<<cc<<nYrs<<cc<<harvStrat<<")"<<endl;
    
    //initialize "current" spawning stock biomass
    currSpB = inpCurrSpB;
    
    //initialize prior spawning stock biomass
    for (int i=1;i<=recLag;i++) simsMMB(iSim,1-i) = asmtModSpB(mMxYr-i);
    
    //initialize numbers-at-size
    for (int x=1;x<=nXs;x++){
        for (int s=1;s<=nSCs;s++){
            nAtZ_xsmz(x,s) = inpNatZ_xsmz(x,s);
        }
    }

    if (debug) cout<<"finished initPopDyMod("<<iSim<<cc<<nYrs<<cc<<harvStrat<<")"<<endl;
    
//-------------------------------------------------------------------------------------
FUNCTION void runPopDyModOneYear(int iSim, int yr, int harvStrat)
    int debug= 1;
    if (debug) cout<<"Starting runPopDyModOneYear("<<iSim<<cc<<yr<<cc<<harvStrat<<")"<<endl;
    
    //Advance pop model to fisheries time
    simsTotB(iSim,yr) = advancePopToFisheries(nAtZ_xsmz, nAtZBF_xsmz);

    //Determine Fofl, OFL, actual F based on harvest strategy/rebuilding plan.
    //Set fishing mortality rates and calculate catches.
    switch (harvStrat) {
        //Case 1:
        //  ??
        case 1: doRebuildingCase01(iSim, yr);
        break;
        
        //Case 2: 
        //  current ADFG harvest strategy and estimates of survey biomass
        case 2: doRebuildingCase02(iSim, yr);
        break;
        
        
        //Case 3:
        //  F=M_imm policy, use mature new shell male M_imm
        case 3: doRebuildingCase03(iSim, yr);
        break;
        
        //Tuer 4:
        //   using current SY overfishing definition
        case 4: doRebuildingCase04(iSim, yr);
        break;
    
        //Case 5:
        //  Tier 3 F% and spr%
        case 5: doRebuildingCase05(iSim, yr);
        break;    
        
        //Case 6:
        //  1. Get F from harvest control rule fmTCF is fmsy
        //      for 2003 harvest control rule set fmTCF = 0.3, bmsy=921.6 and ctrlAlpha = -0.35 (in data file)
        //  Note: to get exploitation rate on mature male biomass need to check for 58% cap and then
        //        calculate the catch - then have to convert to an F.
        //  Do the harvest strategy for the ending biomass from assessment model here then get F for next 
        //      year after get the biomass for next year.
        //  Set fmTCF=0.0 to get no fishing at all.
        case 6: doRebuildingCase06(iSim, yr);
        break;
           
        //Case 7: 
        //  1. Calculate both OFL control rule using tier 3 or fmsy and ADFG harvest strategy.
        //  2. Use whichever F is lower.
        //  3. Keep track of both F's for output to compare
        case 7:  doRebuildingCase07(iSim, yr);
        break;
        
        //Case 8: 
        //  average catch strategy   
        case 8: doRebuildingCase08(iSim, yr);
        break;
    
        //Case 10:
        // Use rebuilding rates until rebuilt in 2 years, then switch to other strategy
        case 10: doRebuildingCase10(iSim, yr);
        break;
    
        //Case 11: 
        //  switches strategies after rebuilt using tier 3 75% *F% 
        case 11: doRebuildingCase11(iSim, yr);
        break;
    
        //Case 12:
        //  if using tier 3 F% and spr% changing by year
        case 12: doRebuildingCase12(iSim, yr);
        break;
        
        //Case 13:
        //  This strategy uses 75% fofl in first year then switches to constant F after 
        case 13: doRebuildingCase13(iSim, yr);
        break;
        
        //Case 14:
        //  Groundfish bycatch only set to maximum for COBLZ catch strategy.
        //  Min bycatch is 4.5 mill crab to maximum 13.5 mill crab, 
        //    with bycatch estimated by 0.001133 times total abundance in numbers.
        case 14: doRebuildingCase14(iSim, yr);
        break;    
    }  //end of switch for harvest strategy
    
    //calculate recruitment
    double R = calcSRFunction(value(recR0),value(recH),phi0,simsMMB(yr-recLag,iSim),srType);
    if (yr==1) {
        simsRecDevs(iSim,1)  = sqrt(value(recVar))*randn(rng);//??: is this correct w/ recVar??
    } else {
        simsRecDevs(iSim,yr)  = recGamma*simsRecDevs(iSim,yr-1)+sqrt(1.0-square(recGamma))*sqrt(value(recVar))*randn(rng);//??: is this correct w/ recVar??
    }
    simsRecruits(iSim,yr) = R*mfexp(simsRecDevs(iSim,yr));
    
    //Advance pop model from after fisheries to end of year
    dvector sbMT = advancePopAfterFisheries(simsRecruits(iSim,yr), nAtZAF_xsmz, nAtZ_xsmz);
    simsMMB(iSim,yr) = sbMT(MALE);
    simsMFB(iSim,yr) = sbMT(FEMALE);
    
    
    if (debug) cout<<"finished runPopDyModOneYear("<<iSim<<cc<<yr<<cc<<harvStrat<<")"<<endl;
    
//-------------------------------------------------------------------------------------
//Fishing at directed F = 0.
FUNCTION void doRebuildingCase01(int iSim, int yr)
    int debug = 1;
    if (debug) cout<<"starting doRebuildingCase01("<<iSim<<cc<<yr<<")"<<endl;
    double tcF = 0.0;
    double scF = fmSCF;
    double rkF = fmRKF;
    double gtF = fmGTF;
    int optTCF_EW = 1;
    setFishingMortalityRates(tcF, scF, rkF, gtF, optTCF_EW);
    dmatrix res = doFisheries(nAtZBF_xsmz, nAtZAF_xsmz);
    //save results
    simsFofl(iSim,yr)  = Fofl;
    simsOFL(iSim,yr)   = OFL;
    simsGHL(iSim,yr)   = ghl;
    simsRetF(iSim,yr)  = tcF;
    simsTotF(iSim,yr)  = tcF+scF+rkF+gtF;
    simsTotY(iSim,yr)  = sum(res(2));
    simsRetY(iSim,yr)  = sum(res(4));
    if (debug) cout<<"finished doRebuildingCase01("<<iSim<<cc<<yr<<")"<<endl;

//-------------------------------------------------------------------------------------
//Fishing at Fofl = F35%
FUNCTION void doRebuildingCase02(int iSim, int yr)
    double tcF = FXX;
    double scF = fmSCF;
    double rkF = fmRKF;
    double gtF = fmGTF;
    int optTCF_EW = 1;
    setFishingMortalityRates(tcF, scF, rkF, gtF, optTCF_EW);
    dmatrix res = doFisheries(nAtZBF_xsmz, nAtZAF_xsmz);
    //save results
    simsFofl(iSim,yr)  = Fofl;
    simsOFL(iSim,yr)   = OFL;
    simsGHL(iSim,yr)   = ghl;
    simsRetF(iSim,yr)  = tcF;
    simsTotF(iSim,yr)  = tcF+scF+rkF+gtF;
    simsTotY(iSim,yr)  = sum(res(2));
    simsRetY(iSim,yr)  = sum(res(4));

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase03(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase04(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase05(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase06(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase07(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase08(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase09(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase10(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase11(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase12(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase13(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION void doRebuildingCase14(int iSim, int yr)

//-------------------------------------------------------------------------------------
FUNCTION double calcFofl(double B, double Bmsy, double Fmsy)
    double Fofl = 0.0;
    if (B/Bmsy<=ctrlBeta) return Fofl; //Fofl = 0
    //otherwise
    if (1.0<=(B/Bmsy)) {
        Fofl = Fmsy;
    } else {
        Fofl = Fmsy*((B/Bmsy)-ctrlAlpha)/(1-ctrlAlpha);
    }
    return Fofl;
    
//-------------------------------------------------------------------------------------
//From the current nAtZ_zsmx, this function calculates: 
//  nAtZBF_xsmz: the numbers-at-size just prior to (i.e., before) the fisheries
//
//Inputs: none
//--nAtZ_xsmz   = numbers at size on July 1 (d4_array) 
//--nAtZBF_xsmz = numbers at size on before fishery (d4_array). WILL BE MODIFIED ON COMPLETION
//
//Output: none
//--totB = total biomass on July 1
//
FUNCTION double advancePopToFisheries(d4_array& nAtZ_xsmz, d4_array& nAtZBF_xsmz)
    int dbg   = 0;
    int debug = 0;
    ofstream post;
    if (dbg) post.open("advancePopToFisheries.dat");
    if (debug) cout<<"Starting advancePopToFisheries"<<endl;
    
    //calculate total biomass on July 1 (NOT at mating time)
    totB = 0.0;
    for (int x=1;x<=nXs;x++) {
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                totB += nAtZ_xsmz(x,sc,ms)*wtp_xmz(x,ms);//note dot-product here
            }
        }
    }
    if (debug) cout<<"totB = "<<totB<<endl;
    
    //advance numbers at size by sex
    nAtZBF_xsmz.initialize();
    for (int x=1;x<=nXs;x++){
        //numbers just prior to fishery: nAtZ01_xsmz
        if (debug) cout<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                nAtZBF_xsmz(x,sc,ms) = mfexp(-midptFishery*M_xsm(x,sc,ms))*nAtZ_xsmz(x,sc,ms);
                if (debug) cout<<"nAtZBF("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZBF_xsmz(x,sc,ms)<<endl;
            }
        }
    }//sex loop
    
    if (dbg){
        post.close();
    }
    if (debug) cout<<"Finished advancePopToFisheries"<<endl;
    return totB;
           
//-------------------------------------------------------------------------------------
//Sets fishing mortality rates and survival in all fisheries
//  (tanner crab, snow crab, red king crab, groundfish trawl).
//
//Inputs:
//--tcF : fully-selected male fishing mortality rate in the directed tanner crab fishery
//--scF : fully-selected fishing mortality rate in the snow crab fishery
//--rkF : fully-selected fishing mortality rate in the red king crab fishery
//--gtF : fully-selected fishing mortality rate in the groundfish trawl fishery
//
//Output: none
//
//Modified:
//--tcTotF_xsmz
//--tcTotF_East_xsmz
//--tcTotF_West_xsmz
//--tcRetF_xsmz
//--tcRetF_East_xsmz
//--tcRetF_West_xsmz
//--scF_xsmz
//--rkF_xsmz
//--gtF_xsmz
//--totF_xsmz
//--totS_xsmz
//
FUNCTION void setFishingMortalityRates(double tcF, double scF, double rkF, double gtF, int optTCF_EW)
    int dbg   = 0;
    int debug = 0;
    if (debug) cout<<"Starting setFishingMortalityRates("<<tcF<<cc<<scF<<cc<<rkF<<cc<<gtF<<cc<<optTCF_EW<<")"<<endl;
    ofstream post;
    if (dbg) post.open("setFishingMortalityRates.dat");
    if (dbg) post<<"Starting setFishingMortalityRates("<<tcF<<cc<<scF<<cc<<rkF<<cc<<gtF<<cc<<optTCF_EW<<")"<<endl;
    
    //calculate fishing mortality rates        
    for(int x=1;x<=nXs;x++) {
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                //in directed Tanner crab fishery 
                if (x==FEMALE) {
                    //on all females (discards only)
                    tcTotF_xsmz(x,sc,ms)      = selTCF_Female*tcF;
                    tcTotF_East_xsmz(x,sc,ms) = elem_prod(propEast,selTCF_Female)*tcF;
                    tcTotF_West_xsmz(x,sc,ms) = elem_prod(propWest,selTCF_Female)*tcF;
                    tcRetF_xsmz(x,sc,ms)      = 0.0;
                    tcRetF_East_xsmz(x,sc,ms) = 0.0;
                    tcRetF_West_xsmz(x,sc,ms) = 0.0;
                } else if (x==MALE) {
                    //on all males (includes discards)
                    tcTotF_East_xsmz(x,sc,ms) = selTCF_TotMaleEast(sc)*tcF;
                    tcTotF_West_xsmz(x,sc,ms) = selTCF_TotMaleWest(sc)*tcF;
                    //on retained males
                    tcRetF_East_xsmz(x,sc,ms) = selTCF_RetMaleEast(sc)*tcF;
                    tcRetF_West_xsmz(x,sc,ms) = selTCF_RetMaleWest(sc)*tcF;
                    if (optTCF_EW==1) {
                        //original (Jack's) approach
                        tcTotF_xsmz(x,sc,ms) = -log(elem_prod(propEast,mfexp(-tcTotF_East_xsmz(x,sc,ms)))+
                                                    elem_prod(propWest,mfexp(-tcTotF_West_xsmz(x,sc,ms))));
                        tcRetF_xsmz(x,sc,ms) = -log(elem_prod(propEast,mfexp(-tcRetF_East_xsmz(x,sc,ms)))+
                                                    elem_prod(propWest,mfexp(-tcRetF_West_xsmz(x,sc,ms))));
                        for (int z=1;z<=nZs;z++){
                            if (tcTotF_xsmz(x,sc,ms,z)<1.0e-20) tcTotF_xsmz(x,sc,ms,z) = 0.0;
                            if (tcRetF_xsmz(x,sc,ms,z)<1.0e-20) tcRetF_xsmz(x,sc,ms,z) = 0.0;
                        }
                    } else if (optTCF_EW==2) {
                        tcTotF_xsmz(x,sc,ms) = elem_prod(propEast,tcTotF_East_xsmz(x,sc,ms))+
                                                elem_prod(propWest,tcTotF_West_xsmz(x,sc,ms));
                        tcRetF_xsmz(x,sc,ms) = elem_prod(propEast,tcRetF_East_xsmz(x,sc,ms))+
                                                elem_prod(propWest,tcRetF_West_xsmz(x,sc,ms));
                    }
                }
                //in other fisheries
                scF_xsmz(x,sc,ms) = selSCF(x)*scF;
                rkF_xsmz(x,sc,ms) = selRKF(x)*rkF;
                gtF_xsmz(x,sc,ms) = selGTF(x)*gtF;
                
                //total fishing mortality and survival
                totF_xsmz(x,sc,ms) = tcTotF_xsmz(x,sc,ms)+scF_xsmz(x,sc,ms)+rkF_xsmz(x,sc,ms)+gtF_xsmz(x,sc,ms);
                totS_xsmz(x,sc,ms) = mfexp(-totF_xsmz(x,sc,ms));
            }
        }
    }
    
    if (debug) {
        writeFishingMortalityRates(tcF,scF,rkF,gtF,post);
        post.close();
    }
    
    if (debug) cout<<"finished setFishingMortalityRates("<<tcF<<cc<<scF<<cc<<rkF<<cc<<gtF<<cc<<optTCF_EW<<")"<<endl;
    
//-------------------------------------------------------------------------------------
//From the current numbers-at-size before the fishery (nAtZBF_zsmx), this function calculates: 
//  nAtZAF_xsmz: the numbers-at-size just after the fisheries
//
//Inputs: 
//--nAtZBF_xsmz    = the numbers-at-size just before the fisheries.
//--nAtZAF_xsmz    = the numbers-at-size just after the fisheries.  WILL BE MODIFIED ON COMPLETION
//
//Output: 
//--dmatrix(1,14,1,nXs) with sex-specific catch, yield for the various fisheries
//     output( 1) = totCat;
//     output( 2) = totYld;
//     output( 3) = retCat;
//     output( 4) = retYld;
//     output( 5) = tcTotCat;
//     output( 6) = tcTotYld;
//     output( 7) = tcRetCat;
//     output( 8) = tcRetYld;
//     output( 9) = scTotCat;
//     output(10) = scTotYld;
//     output(11) = rkTotCat;
//     output(12) = rkTotYld;
//     output(13) = gtTotCat;
//     output(14) = gtTotYld;
//
FUNCTION dmatrix doFisheries(d4_array& nAtZBF_xsmz, d4_array& nAtZAF_xsmz)

    int dbg   = 0;
    int debug = 0;
    ofstream post;
    if (dbg) post.open("doFisheries.dat");
    if (debug) cout<<"Starting doFisheries"<<endl;    
    
    //quantities for sex-specific catch, yield
    dmatrix output(1,14,1,nXs);
    dvector totCat(1,nXs);  //catch
    dvector totYld(1,nXs);  //yield          (units??)
    dvector retCat(1,nXs);  //retained catch
    dvector retYld(1,nXs);  //retained yield (units??)
    dvector tcTotCat(1,nXs);  //catch in directed Tanner crab fishery
    dvector tcTotYld(1,nXs);  //yield in directed Tanner crab fishery (units??)
    dvector tcRetCat(1,nXs);  //retained catch in directed Tanner crab fishery
    dvector tcRetYld(1,nXs);  //retained yield in directed Tanner crab fishery (units??)
    dvector scTotCat(1,nXs);  //catch in snow crab fishery
    dvector scTotYld(1,nXs);  //yield in snow crab fishery (units??)
    dvector rkTotCat(1,nXs);  //catch in red king crab fishery
    dvector rkTotYld(1,nXs);  //yield in red king crab fishery (units??)
    dvector gtTotCat(1,nXs);  //catch in groundfish trawl crab fishery
    dvector gtTotYld(1,nXs);  //yield in groundfish trawl crab fishery (units??)
    
    nAtZAF_xsmz.initialize();
    for (int x=1;x<=nXs;x++) {
        //numbers after fishery: nAtZAF_xsmz
        if (debug) cout<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                nAtZAF_xsmz(x,sc,ms) = elem_prod(totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms));
                if (debug) cout<<"nAtZAF("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZAF_xsmz(x,sc,ms)<<endl;
            }
        }
        //catch-at-size, sex-specifc catch, sex-specific yield
        totCat(x) = 0.0;
        totYld(x) = 0.0;
        retCat(x) = 0.0;
        retYld(x) = 0.0;
        tcTotCat(x) = 0.0;
        tcTotYld(x) = 0.0;
        tcRetCat(x) = 0.0;
        tcRetYld(x) = 0.0;
        scTotCat(x) = 0.0;
        scTotYld(x) = 0.0;
        rkTotCat(x) = 0.0;
        rkTotYld(x) = 0.0;
        gtTotCat(x) = 0.0;
        gtTotYld(x) = 0.0;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                tcRetCatZ_xsmz(x,sc,ms) = elem_prod(elem_div(tcRetF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                    elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                tcTotCatZ_xsmz(x,sc,ms) = elem_prod(elem_div(tcTotF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                    elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                scCatZ_xsmz(x,sc,ms)    = elem_prod(elem_div(scF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                    elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                rkCatZ_xsmz(x,sc,ms)    = elem_prod(elem_div(rkF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                    elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                gtCatZ_xsmz(x,sc,ms)    = elem_prod(elem_div(gtF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                    elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                                                    
                totCat(x)   += sum(tcTotCatZ_xsmz(x,sc,ms)+scCatZ_xsmz(x,sc,ms)+rkCatZ_xsmz(x,sc,ms)+gtCatZ_xsmz(x,sc,ms));
                retCat(x)   += sum(tcRetCatZ_xsmz(x,sc,ms));
                tcTotCat(x) += sum(tcTotCatZ_xsmz(x,sc,ms));
                tcRetCat(x) += sum(tcRetCatZ_xsmz(x,sc,ms));
                scTotCat(x) += sum(scCatZ_xsmz(x,sc,ms));
                rkTotCat(x) += sum(rkCatZ_xsmz(x,sc,ms));
                gtTotCat(x) += sum(gtCatZ_xsmz(x,sc,ms));
                
                totYld(x)   += wtp_xmz(x,ms)*(tcTotCatZ_xsmz(x,sc,ms)+scCatZ_xsmz(x,sc,ms)+rkCatZ_xsmz(x,sc,ms)+gtCatZ_xsmz(x,sc,ms));
                retYld(x)   += wtp_xmz(x,ms)*tcRetCatZ_xsmz(x,sc,ms);
                tcTotYld(x) += wtp_xmz(x,ms)*tcTotCatZ_xsmz(x,sc,ms);
                tcRetYld(x) += wtp_xmz(x,ms)*tcRetCatZ_xsmz(x,sc,ms);
                scTotYld(x) += wtp_xmz(x,ms)*scCatZ_xsmz(x,sc,ms);
                rkTotYld(x) += wtp_xmz(x,ms)*rkCatZ_xsmz(x,sc,ms);
                gtTotYld(x) += wtp_xmz(x,ms)*gtCatZ_xsmz(x,sc,ms);//note dot-product here
            }
        }
        if (debug) {
            cout<<"totCat("<<x<<") = "<<totCat(x)<<endl;
            cout<<"totYld("<<x<<") = "<<totYld(x)<<endl;
            cout<<"retCat("<<x<<") = "<<retCat(x)<<endl;
            cout<<"retYld("<<x<<") = "<<retYld(x)<<endl;
            cout<<"tcTotCat("<<x<<") = "<<tcTotCat(x)<<endl;
            cout<<"tcTotYld("<<x<<") = "<<tcTotYld(x)<<endl;
            cout<<"tcRetCat("<<x<<") = "<<tcRetCat(x)<<endl;
            cout<<"tcRetYld("<<x<<") = "<<tcRetYld(x)<<endl;
            cout<<"scTotCat("<<x<<") = "<<scTotCat(x)<<endl;
            cout<<"scTotYld("<<x<<") = "<<scTotYld(x)<<endl;
            cout<<"rkTotCat("<<x<<") = "<<rkTotCat(x)<<endl;
            cout<<"rkTotYld("<<x<<") = "<<rkTotYld(x)<<endl;
            cout<<"gtTotCat("<<x<<") = "<<gtTotCat(x)<<endl;
            cout<<"gtTotYld("<<x<<") = "<<gtTotYld(x)<<endl;
        }
    }//sex loop
    output( 1) = totCat;
    output( 2) = totYld;
    output( 3) = retCat;
    output( 4) = retYld;
    output( 5) = tcTotCat;
    output( 6) = tcTotYld;
    output( 7) = tcRetCat;
    output( 8) = tcRetYld;
    output( 9) = scTotCat;
    output(10) = scTotYld;
    output(11) = rkTotCat;
    output(12) = rkTotYld;
    output(13) = gtTotCat;
    output(14) = gtTotYld;
    
    if (dbg) {
        post.close();
    }
    if (debug) cout<<"finished doFisheries"<<endl;
    
    return(output);
    
//-------------------------------------------------------------------------------------
//From the current numbers-at-size after the fishery (nAtZAF_zsmx), this function 
//  advances the population to the end of the current year/beginning of the next year.
//
//Inputs: none
//--recs        = total recruitment at end of year
//--nAtZAF_xsmz = numbers at size AFTER fishery (d4_array)
//--nAtZ_xsmz   = numbers at size at end of year/beginning of next year MODIFIED ON OUTPUT
//Output: 
//--sbMT = spawning biomass by sex at mating time
//
FUNCTION dvector advancePopAfterFisheries(double recs, d4_array& nAtZAF_xsmz, d4_array& nAtZ_xsmz)
    int dbg   = 1;
    int debug = 0;
    ofstream post;
    if (dbg) post.open("advancePopAfterFisheries.dat");
    if (debug) cout<<"Starting advancePopAfterFisheries"<<endl;
    
//     d4_array nAtZ03_xsmz = nAtZ_xsmz;//temporary storage
//     d4_array nAtZ04_xsmz = nAtZ_xsmz;//temporary storage
//     d4_array nAtZ05_xsmz = nAtZ_xsmz;//temporary storage
//     d4_array nAtZMT_xsmz = nAtZ_xsmz;//temporary storage
    
    double mateTime = 0.0;//time from end of fisheries to mating
    for (int x=1;x<=nXs;x++) {
        //number at mating time (after fishery, before growth)
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                nAtZMT_xsmz(x,sc,ms) = mfexp(-mateTime*M_xsm(x,sc,ms))*nAtZAF_xsmz(x,sc,ms);
                if (debug) cout<<"nAtZAF("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZAF_xsmz(x,sc,ms)<<endl;
                if (debug) cout<<"nAtZMT("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZMT_xsmz(x,sc,ms)<<endl;
            }
        }
        //numbers that will molt, just after mating, with change in shell condition: nAtZ03_xsmz
        if (debug) cout<<endl;
        for(int ms=1;ms<=nMSs;ms++) {
            for(int sc=1;sc<=nSCs;sc++) { //to sc
                nAtZ03_xsmz(x,sc,ms) = 0.0;
                for(int sp=1;sp<=nSCs;sp++) { //from sp
                    nAtZ03_xsmz(x,sc,ms) += elem_prod(tmMolt_xmssz(x,ms,sc,sp),nAtZMT_xsmz(x,sp,ms));
                }
                if (debug) cout<<"nAtZ03("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ03_xsmz(x,sc,ms)<<endl;
            }
        }
        
        //numbers after molting and growth: nAtZ04_xsmz
        if (debug) cout<<endl;
        for(int ms=1;ms<=nMSs;ms++) {
            for(int sc=1;sc<=nSCs;sc++) {
                //crabs that molted (became new shell) and grow
                nAtZ04_xsmz(x,NEW_SHELL,ms) = nAtZ03_xsmz(x,NEW_SHELL,ms)*tmZtoZ_xzz(x);       
                //crabs that didn't molt (and thus didn't grow)
                nAtZ04_xsmz(x,OLD_SHELL,ms) = nAtZ03_xsmz(x,OLD_SHELL,ms);
                if (debug) cout<<"nAtZ04("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ04_xsmz(x,sc,ms)<<endl;
            }
        }
        
        //apply maturity transitions from mp to ms: nAtZ05_xsmz
        if (debug) cout<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) { //to ms
                nAtZ05_xsmz(x,sc,ms) = 0.0;
                for(int mp=1;mp<=nMSs;mp++) { //from mp
                    nAtZ05_xsmz(x,sc,ms) += elem_prod(tmMat_xsmmz(x,sc,ms,mp),nAtZ04_xsmz(x,sc,mp));
                }
                if (debug) cout<<"nAtZ05("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ05_xsmz(x,sc,ms)<<endl;
            }
        }
        
        //apply end-of-year mortality: nAtZ_xsmz
        if (debug) cout<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                nAtZ_xsmz(x,sc,ms) = mfexp(-(1.0-(midptFishery+mateTime))*M_xsm(x,sc,ms))*nAtZ05_xsmz(x,sc,ms);
                if (debug) cout<<"nAtZ("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ_xsmz(x,sc,ms)<<endl;
            }
        }
        
        //add in recruitment <--wts: factor in actual number of simsRecruits
        for(int sc=1;sc<=nSCs;sc++) {
            rec_xsz(x,sc) = recs*recPropBySex(x)*recPropByShell(sc)*recPropAtZ;
            nAtZ_xsmz(x,sc,IMMATURE) += rec_xsz(x,sc);
        }
        
        if (debug) {
            cout<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    //post<<x<<sc<<ms<<nAtZ_xsmz(x,sc,ms)<<endl;
                    cout<<"nAtZ("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ_xsmz(x,sc,ms)<<endl;
                }
            }
        }
    } //sex loop
            
    // sp biomass for final: want units in 1000's t. 
    //  wtp in t/indiv, nAtZ in 1000's indiv => nAtZ*wtp yields 1000's t 
    if (debug) cout<<endl;
    for (int x=1;x<=nXs;x++) {
        sbMT(x) = 0.0;//initialize for sum
        for(int sc=1;sc<=nSCs;sc++) {
            sbMT(x) += nAtZMT_xsmz(x,sc,MATURE)*wtp_xmz(x,MATURE);//note dot-product here
        }
        if (debug) {
            cout<<"sbMT("<<x<<") = "<<sbMT(x)<<endl;
        }
    }
    
    if (dbg) {
        writePopQuantsAfterFisheries(post);
        post.close();
    }
    
    if (debug) {
        cout<<"Finished advancePopAfterFisheries"<<endl<<endl;
        cout<<"Press 1 and hit return to continue >> ";
        cin>>resp;
    }
    
    return sbMT;

//-------------------------------------------------------------------------------------
//finds h such that Fmsy = FXX%.
//
//Inputs:
//--FXX     = desired value of F such that Fmsy = FXX%
//--R0      = stock-recruit function parameter--equilibium recruitment for unfished stock (i.e., at F=0 in all fisheries)
//--h       = stock-recruit function parameter--steepness 
//--phi0    = unfished spawning stock biomass per recruit (MMB/R at F=0 in all fisheries)
//--optF_EW = flag for handling combination of directed fishing mortality E/W of lon=166W
//--srType  = flag indicating stock-recruit function
//
//Outputs:
//--double  = value of h that results in Fmsy = FXX%, or < 0 if none found
//
FUNCTION double findH(double FXX, double R0, double phi0, int optTCF_EW, int srType)
    int debug = 0;
    if (debug) cout<<"starting findH"<<endl;
    double h0     = 0.5;   //initial guess for H that yields Fmsy = F35%
    double dh     = 0.001; //increment in h for calculating dFmsy/dh
    double cvgCri = 1.0e-6;//convergence criteria
    int itr       = 0;     //iteration counter
    double dhp    = 10;    //increment for updating h0 (initial value is large to get into while loop)
    double Fm, Fp, dFdh;
    dvector resMSY(1,5);   //result vector from MSY calculation
    while ((cvgCri<fabs(dhp))&&(itr<200)) {
        itr++;
        if (debug) cout<<"-----xxH iteration = "<<itr<<".  "<<endl;
        resMSY = calcMSY(R0,h0-dh,phi0,optTCF_EW,srType);
        Fm     = resMSY(1);             //Fmsy at h0-dh
        resMSY = calcMSY(R0,h0+dh,phi0,optTCF_EW,srType);
        Fp     = resMSY(1);             //Fmsy at h0+dh
        dFdh   = (Fp-Fm)/(2.0*dh);      //approximate derivative dFmsy/dh at h0
        dhp    = (FXX-0.5*(Fm+Fp))/dFdh;//increment for updating h0
        if (debug) cout<<endl<<"-----xxH iteration results = "<<h0<<tb<<dhp<<tb<<FXX<<tb<<0.5*(Fm+Fp)<<tb<<0.5*(Fm+Fp)/FXX<<endl;
        h0 += dhp;//update h0
        if (h0<0.2) {h0 = 0.2; cout<<"found h<0.2!!"<<endl; break;}
        if (h0>1.0) {h0 = 1.0; cout<<"found h>1.0!!"<<endl; break;}
    }
    
    if (debug) cout<<"finished findH"<<endl;
    return(h0);


//-------------------------------------------------------------------------------------
//finds h such that Fmsy = FXX%.
//
//Inputs:
//--R0      = stock-recruit function parameter--equilibium recruitment for unfished stock (i.e., at F=0 in all fisheries)
//--h       = stock-recruit function parameter--steepness 
//--phi0    = unfished spawning stock biomass per recruit (MMB/R at F=0 in all fisheries)
//--optF_EW = flag for handling combination of directed fishing mortality E/W of lon=166W
//--srType  = flag indicating stock-recruit function
//
//Outputs:
//--double  = value of h that results in Fmsy = FXX%, or < 0 if none found
//
FUNCTION dmatrix calcHs(double FXX, double R0, double phi0, int optTCF_EW, int srType)
    int debug = 0;
    if (debug) cout<<"starting calcHs"<<endl;
    double nN = 21;
    dvector h(1,nN);
    dvector Fmsy(1,nN);
    double dh = (0.995-0.205)/(nN-1);
    dmatrix res(1,nN,1,6);
    for (int n=1;n<=nN;n++){
        res(n)   = 0.0;         
        h(n)     = 0.205+dh*(n-1);
        res(n,1) = h(n);
        res(n)(2,6).shift(1) = calcMSY(R0,h(n),phi0,optTCF_EW,srType);
        if (debug) cout<<n<<tb<<res(n)<<endl;
    }
    
    if (debug) cout<<"finished calcHs"<<endl;
    return(res);

//-------------------------------------------------------------------------------------
//calculates quantities related to Fmsy, msy, Bmsy, Rmsy (recruitment), SBPR(F=Fmsy).
//
//Inputs:
//--R0      = stock-recruit function parameter--equilibium recruitment for unfished stock (i.e., at F=0 in all fisheries)
//--h       = stock-recruit function parameter--steepness 
//--phi0    = unfished spawning stock biomass per recruit (MMB/R at F=0 in all fisheries)
//--optF_EW = flag for handling combination of directed fishing mortality E/W of lon=166W
//--srType  = flag indicating stock-recruit function
//
//Outputs:
//--dvector res(1,5);
//----res(1) = Fmsy 
//----res(2) = msy       : retained male catch biomass in directed fishery at Fmsy
//----res(3) = Bmsy      : male spawning stock biomass at Fmsy
//----res(4) = Rmsy      : equilibrium recruitment at Fmsy
//----res(5) = SBPR(Fmsy): male spawning stock biomass/recruit at Fmsy (
FUNCTION dvector calcMSY(double R0, double h, double phi0, int optTCF_EW, int srType)
    int debug = 0;
    if (debug) cout<<"start calcMSY"<<endl;
    double d  = 0.00001;//delta F for derivative calcs
    double F1 = 0.01;    //initial F
    double F2;
    double F3;
    double F1p;
    double yld1;
    double yld2;
    double yld3;
    double F1prime;
    double F2prime;
    double cvg     = 1.0;   //convergence indicator
    double convCri = 1.0e-6;//covergence criteria
    int itr = 1;            //iteration counter
    while ((convCri<cvg)&&(itr<=20)){
        F2    = F1+d;
        F3    = F1-d;
        yld1  = calcEqYield(F1,R0,h,phi0,optTCF_EW);
        yld2  = calcEqYield(F2,R0,h,phi0,optTCF_EW);
        yld3  = calcEqYield(F3,R0,h,phi0,optTCF_EW);
        F1prime= (yld2-yld3)/(2*d);         //  Newton-Raphson approximation for first derivitive
        F2prime=(yld3-(2*yld1)+yld2)/(d*d); // Newton-Raphson approximation for second derivitive
        F1p    = F1-F1prime/F2prime;        // incremented F for next time through the loop
        cvg    = fabs((F1prime/F2prime)/F1);// convergence indicator
        if (debug) cout <<"F1 = "<<F1<<" F1p= "<<F1p<<" F1prime= "<<F1prime<<" cvg = "<<cvg<<endl;
        if(F1p>0) {
            F1 = F1p;
        } else {
            F1=0.000001;
        }
        ++itr;
    }
    if (debug) cout<<"#--calcMSY: out of iteration loop. "<<"Fmsy = "<<F1<<"  cvg = "<<cvg<<endl;
    double Fmsy;
    double Bmsy;
    double msy;
    double phi;
    double ypr;
    double R_eq;
    if((F1p>0) && (F1>0.000001)) {
        Fmsy    = F1;
        int optCalcYPR = 1;//calculate ypr
        dmatrix res = calcSBPRF(Fmsy,optCalcYPR,optTCF_EW);
        phi     = res(1,MALE);                         //MMB/R at mating time
        ypr     = res(idxYPR,MALE);                   //yield-per-recruit for retained males
        R_eq    = calcEqRec(R0,h,phi0,phi/phi0,srType);//equilibrium recruitment at Fmsy
        Bmsy    = R_eq*phi;
        msy     = R_eq*ypr;
    } else {
        Fmsy = -1.0;
        Bmsy = -1.0;
        msy  = -1.0;
        R_eq = -1.0;
        phi  = -1.0;
    }
    dvector res(1,5);
    res(1) = Fmsy;
    res(2) = msy;
    res(3) = Bmsy;
    res(4) = R_eq;
    res(5) = phi;
    
    if (debug) cout<<"end calcMSY"<<endl;
    return(res);


//-------------------------------------------------------------------------------------
//Calculates equilibrium yield (based on retained male catch biomass) for given
//directed fishing mortality rate.
//
//Inputs:
//--tcF     = directed fishing mortality on fully-selected males
//--R0      = stock-recruit function parameter--equilibium recruitment for unfished stock (i.e., at F=0 in all fisheries)
//--h       = stock-recruit function parameter--steepness 
//--phi0    = unfished spawning stock biomass per recruit (MMB/R at F=0 in all fisheries)
//--optF_EW = flag for handling combination of directed fishing mortality E/W of lon=166W
//
//Output:
//--Fxx = directed fishing mortality value that reduces SBPR to xx of that for an unfished stock
//
FUNCTION double calcEqYield(double tcF, double R0, double h, double phi0, int optTCF_EW)
    int debug = 0;
    if (debug) cout<<"start calcEqYield"<<endl;
    // calculates yield for any value of F
    int optCalcYPR = 1;//calc cpr/ypr values
    dmatrix res = calcSBPRF(tcF,optCalcYPR,optTCF_EW);
    double phi  = res(1,MALE);                         //MMB/R at mating time
    double ypr  = res(idxYPR,MALE);                    //yield-per-recruit for retained males
    double R_eq = calcEqRec(R0,h,phi0,phi/phi0,srType);//equilibrium recruitment
    double yield = R_eq*ypr;                           //equilibrium yield for retained males 
    if (debug) cout<<"end calcEqYield"<<endl;
    return(yield);

//-------------------------------------------------------------------------------------
//calculates the directed fishing mortality that results in spawning stock biomass per recruit
//being reduced to xx times its value for an unfished stock
//
//Inputs:
//--xx = fraction of unfished spawning stock biomass per recruit
//--phi0 = unfished spawning stock biomass per recruit (MMB at F=0 in all fisheries)
//--optF_EW = flag for handling combination of directed fishing mortality E/W of lon=166W
//
//Output:
//--Fxx = directed fishing mortality value that reduces SBPR to xx of that for an unfished stock
//
FUNCTION double findFxx(double xx, double phi0, int optF_EW)
    int debug = 0;
    if (debug) cout<<"starting findFxx("<<xx<<","<<phi0<<")"<<endl;
        
    //iterate on directed Tanner crab fishing mortality (tcFxx) to find
    //tcFxx such that xx = sbprF/phi0
    double tcFxxp = 0.3;    //initial guess for tcFxx
    double xxp    = 0;      //sbprF/phi0 for tcFxxp
    double cvgCri = 1.0e-6; //convergence criteria
    double scl    = 0.0;    //scale factor for updating tcFxxp
    int optCalcYPR = 0;     //don't calc ypr in calcSBPRF
    int itr        = 0;     //iteration counter
    dmatrix res(1,1,1,nXs); //results of calcSBPRF(...)
    if (debug) cout<<endl;
    while ((cvgCri<fabs(1.0-scl))&&(itr<200)) {
        itr++;
        if (debug) cout<<"-----Fxx iteration = "<<itr<<".  "<<endl;
        res    = calcSBPRF(tcFxxp,optCalcYPR,optF_EW);
        xxp    = res(1,MALE)/phi0;//assumes sbpr currency is MMB at mating time
        scl    = xxp/xx;
        tcFxxp = scl * tcFxxp;
        if (debug) cout<<endl<<"-----Fxx iteration results = "<<res(1,MALE)<<" "<<xxp<<" "<<scl<<endl;
    }
    
    return tcFxxp;

//-------------------------------------------------------------------------------------
//Sets biological quantities required for SBPR and population dynamics calculations. 
//This function must be called once prior to calling calcSBPR0 or calcSBPRF.
//
FUNCTION setBioQuants
    int debug = 0;
    if (debug) cout<<endl<<endl<<"starting setBioQuants"<<endl;
        
    //proportion recruiting by sex
    recPropBySex(FEMALE) = 1.0-recPropAsMale;
    recPropBySex(  MALE) =     recPropAsMale;
        
    //proportion recruiting by shell condition
    recPropByShell(NEW_SHELL) =     recPropAsNewShell;
    recPropByShell(OLD_SHELL) = 1.0-recPropAsNewShell;
    
    //weight (t) at size by sex, maturity state
    wtp_xmz(FEMALE,IMMATURE) = wtjf;
    wtp_xmz(FEMALE,  MATURE) = wt(FEMALE);
    wtp_xmz(  MALE,IMMATURE) = wt(MALE);
    wtp_xmz(  MALE,  MATURE) = wt(MALE);
    
    //transition probabilities among shell conditions based on probability crab WILL molt (and change shell condition)
    for (int x=1;x<=nXs;x++) {
        tmMolt_xmssz(x,IMMATURE,NEW_SHELL,NEW_SHELL) =     prMoltImm(x);//molting occurs;        to new shell from new shell 
        tmMolt_xmssz(x,IMMATURE,OLD_SHELL,NEW_SHELL) = 1.0-prMoltImm(x);//molting doesn't occur; to old shell from new shell
        tmMolt_xmssz(x,IMMATURE,NEW_SHELL,OLD_SHELL) =     prMoltImm(x);//molting occurs;        to new shell from old shell (but there shouldn't be any old shell immatures)
        tmMolt_xmssz(x,IMMATURE,OLD_SHELL,OLD_SHELL) = 1.0-prMoltImm(x);//molting doesn't occur; to old shell from old shell (but there shouldn't be any old shell immatures)
        tmMolt_xmssz(x,  MATURE,NEW_SHELL,NEW_SHELL) =     prMoltMat(x);//molting occurs;        to new shell from new shell (for TC, this should be 0 identically) 
        tmMolt_xmssz(x,  MATURE,OLD_SHELL,NEW_SHELL) = 1.0-prMoltMat(x);//molting doesn't occur; to old shell from new shell (for TC, this should be 1 identically)
        tmMolt_xmssz(x,  MATURE,NEW_SHELL,OLD_SHELL) =              0.0;//molting occurs;        to new shell from old shell (can't happen:   matures don't molt)
        tmMolt_xmssz(x,  MATURE,OLD_SHELL,OLD_SHELL) =              1.0;//molting doesn't occur; to old shell from old shell (always happens: matures don't molt)
        if (debug) {
            cout<<endl;
            cout<<"prMoltImm("<<x<<") = "<<prMoltImm(x)<<endl;
            cout<<"prMoltMat("<<x<<") = "<<prMoltMat(x)<<endl;
            for (int ms=1;ms<=nMSs;ms++) {
                for (int sc=1;sc<=nSCs;sc++) {
                    for (int sp=1;sp<=nSCs;sp++) cout<<"tmMolt_xmssz("<<x<<","<<ms<<","<<sp<<","<<sc<<") = "<<tmMolt_xmssz(x,ms,sp,sc)<<endl;
                }
            }
        }
    }
    
    //transition probability from maturity state to maturity state
    for (int x=1;x<=nXs;x++) {
        tmMat_xsmmz(x,NEW_SHELL,IMMATURE,IMMATURE) = 1.0-prMatNS(x); //to immature from immature  (i.e., stays immature)
        tmMat_xsmmz(x,NEW_SHELL,  MATURE,IMMATURE) =     prMatNS(x); //to   mature from immature  (i.e., matures)
        tmMat_xsmmz(x,NEW_SHELL,IMMATURE,  MATURE) = 0.0;                 //to immature from mature (can't happen:   matures don't regress to immature)
        tmMat_xsmmz(x,NEW_SHELL,  MATURE,  MATURE) = 1.0;                 //to   mature from mature (always happens: matures don't regress to immature)
          
        tmMat_xsmmz(x,OLD_SHELL,IMMATURE,IMMATURE) = 0.0; //to immature from immature (but there shouldn't be any old shell immatures)
        tmMat_xsmmz(x,OLD_SHELL,  MATURE,IMMATURE) = 1.0; //to   mature from immature (but there shouldn't be any old shell immatures)
        tmMat_xsmmz(x,OLD_SHELL,IMMATURE,  MATURE) = 0.0; //to immature from mature (can't happen:   matures don't regress to immature)
        tmMat_xsmmz(x,OLD_SHELL,  MATURE,  MATURE) = 1.0; //to   mature from mature (always happens: matures don't regress to immature)
        if (debug) {
            cout<<endl;
            cout<<"prMatNS("<<x<<") = "<<prMatNS(x)<<endl;
            for (int sc=1;sc<=nSCs;sc++) {
                for (int ms=1;ms<=nMSs;ms++) {
                    for (int mp=1;mp<=nMSs;mp++) cout<<"tmMat_xsmmz("<<x<<","<<sc<<","<<mp<<","<<ms<<") = "<<tmMat_xsmmz(x,sc,mp,ms)<<endl;
                }
            }
        }
    }
    
    if (debug) {
        cout<<"Finished setBioQuants"<<endl<<endl;
        cout<<"Press 1 and hit return to continue >> ";
        cin>>resp;
    }

//-------------------------------------------------------------------------------------
//Calculates spawning biomass by sex per recruit for an unfished population.
// Multiply by TOTAL recruitment to get total spawning biomass by sex for an unfished population.
//
//Inputs: none
//
//Output:
// dvector sbpr0(1,nXs) based on sex-specific biomass at mating time.
//  spbr0(FEMALE) = female spawning biomass at mating time
//  spbr0(  MALE) =   male spawning biomass at mating time
//
FUNCTION dvector calcSBPR0(void)
    int dbg   = 0;
    int debug = 0;
    ofstream post;
    if (debug) cout<<endl<<endl<<"starting calcSBPR0"<<endl;    
    if (dbg) post.open("calcSBPR0.dat");
    
//     //numbers at size
//     d4_array nAtZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);   //on July 1
//     d4_array nAtZ01_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //just prior to fishery (if it was going to occur, which it doesn't)
//     d4_array nAtZ02_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //just after to fishery (if it was going to occur, which it doesn't)
//     d4_array nAtZ03_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //only those that will grow
//     d4_array nAtZ04_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //all after growth/no growth
//     d4_array nAtZ05_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //after maturity transition

//     d4_array nAtZFT_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //numbers just prior to fishery (same as nAtZ01)
//     d4_array nAtZMT_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //numbers at mating time (mateTime = time AFTER fishery occurs that mating occurs)
//         
//     d4_array prvNatZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);//on July 1 from PREVIOUS iteration in while loop below
//     
//     //spawning biomass metrics
//     dvector sbpr0(1,nXs);  //on July 1
//     dvector sbpr0FT(1,nXs);//just prior to fishery
//     dvector sbpr0MT(1,nXs);//at mating time (after fishery)
    
    //initialize numbers at size with recruitment
    if (debug) cout<<endl;
    for(int x=1;x<=nXs;x++) {
        for(int sc=1;sc<=nSCs;sc++) {
            nAtZ_xsmz(x,sc,IMMATURE)    = recPropBySex(x)*recPropByShell(sc)*recPropAtZ;
            nAtZ_xsmz(x,sc,  MATURE)    = 0.0;
            prvNatZ_xsmz(x,sc,IMMATURE) = 0.0;
            prvNatZ_xsmz(x,sc,  MATURE) = 0.0;
            if (debug) {
                cout<<"nAtZ("<<x<<","<<sc<<","<<IMMATURE<<")   = "<<nAtZ_xsmz(x,sc,IMMATURE)<<endl;
                cout<<"nAtZ("<<x<<","<<sc<<","<<  MATURE<<")   = "<<nAtZ_xsmz(x,sc,  MATURE)<<endl;
            }
        }
        if (dbg){
            nAtZ_xsmz(x,NEW_SHELL,IMMATURE) = 1.0;  //wtstest
            nAtZ_xsmz(x,NEW_SHELL,  MATURE) = 0.0;
            nAtZ_xsmz(x,OLD_SHELL,IMMATURE) = 0.0;
            nAtZ_xsmz(x,OLD_SHELL,  MATURE) = 0.0;
        }
    }
    // initial sp biomass
    for (int x=1;x<=nXs;x++) {
        sbpr0(x) = 0.0;//initialize for sum
        for(int sc=1;sc<=nSCs;sc++) {
            sbpr0(x) += nAtZ_xsmz(x,sc,MATURE)*wtp_xmz(x,MATURE);//note dot-product here
        }
        if (debug) cout<<"sbpr0("<<x<<") = "<<sbpr0(x)<<".  " ;
    }
    if (debug) cout<<endl;
        
    //run pop dy equations forward until stable equilibrium is achieved
    double cvgCri    = 1.0e-6; //criterion for convergence to stable equilibrium
    double cvg       = 100;    //initial convergence value to get into the while loop
    int itr          = 0;      //iteration counter
    while ((cvgCri<cvg)&&(itr<200)) {   //wtstest
        itr++;
        if (debug) cout<<endl<<"------------------"<<endl<<"iteration = "<<itr<<".  ";
        
        //advance numbers at size by sex
        for (int x=1;x<=nXs;x++){
            //numbers just prior to fishery
            if (debug) cout<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    nAtZBF_xsmz(x,sc,ms) = mfexp(-midptFishery*M_xsm(x,sc,ms))*nAtZ_xsmz(x,sc,ms);
                    if (debug) cout<<"nAtZBF("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZBF_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //numbers after fishery (same as before--no fishery)
            if (debug) cout<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    nAtZAF_xsmz(x,sc,ms) = nAtZBF_xsmz(x,sc,ms);
                    if (debug) cout<<"nAtZAF("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZAF_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //numbers that will molt, with change in shell condition
            if (debug) cout<<endl;
            for(int ms=1;ms<=nMSs;ms++) {
                for(int sc=1;sc<=nSCs;sc++) { //to sc
                    nAtZ03_xsmz(x,sc,ms) = 0.0;
                    for(int sp=1;sp<=nSCs;sp++) { //from sp
                        nAtZ03_xsmz(x,sc,ms) += elem_prod(tmMolt_xmssz(x,ms,sc,sp),nAtZAF_xsmz(x,sp,ms));
                    }
                    if (debug) cout<<"nAtZ03("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ03_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //numbers after molting and growth
            if (debug) cout<<endl;
            for(int ms=1;ms<=nMSs;ms++) {
                //crabs that molted (became new shell) and grow
                nAtZ04_xsmz(x,NEW_SHELL,ms) = nAtZ03_xsmz(x,NEW_SHELL,ms)*tmZtoZ_xzz(x);       
                //crabs that didn't molt (and thus didn't grow)
                nAtZ04_xsmz(x,OLD_SHELL,ms) = nAtZ03_xsmz(x,OLD_SHELL,ms);
                if (debug) for(int sc=1;sc<=nSCs;sc++) cout<<"nAtZ04("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ04_xsmz(x,sc,ms)<<endl;
            }
            
            //apply maturity transitions from mp to ms
            if (debug) cout<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) { //to ms
                    nAtZ05_xsmz(x,sc,ms) = 0.0;
                    for(int mp=1;mp<=nMSs;mp++) { //from mp
                        nAtZ05_xsmz(x,sc,ms) += elem_prod(tmMat_xsmmz(x,sc,ms,mp),nAtZ04_xsmz(x,sc,mp));
                    }
                    if (debug) cout<<"nAtZ05("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ05_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //apply end-of-year mortality
            if (debug) cout<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    nAtZ_xsmz(x,sc,ms) = mfexp(-(1-midptFishery)*M_xsm(x,sc,ms))*nAtZ05_xsmz(x,sc,ms);
                    if (debug) cout<<"nAtZ("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //add in recruitment
            for(int sc=1;sc<=nSCs;sc++) {
                nAtZ_xsmz(x,sc,IMMATURE) += recPropBySex(x)*recPropByShell(sc)*recPropAtZ;
            }
            
            if (debug){ 
                cout<<endl;
                for(int sc=1;sc<=nSCs;sc++) {
                    for(int ms=1;ms<=nMSs;ms++) {
                        //post<<x<<sc<<ms<<nAtZ_xsmz(x,sc,ms)<<endl;
                        cout<<"nAtZ("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ_xsmz(x,sc,ms)<<endl;
                    }
                }
            }
        } //sex loop
        
        //calc sp biomass
        for (int x=1;x<=nXs;x++) {
            sbpr0(x) = 0.0;//initialize for sum
            for(int sc=1;sc<=nSCs;sc++) {
                sbpr0(x) += nAtZ_xsmz(x,sc,MATURE)*wtp_xmz(x,MATURE);//note dot-product here
            }
            if (debug) cout<<"sbpr0("<<x<<") = "<<sbpr0(x)<<".  " ;
        }
        
        //check convergence and save current nAtZ
        cvg = 0.0;
        dvector num(1,nZs);
        dvector den(1,nZs);
        for (int x=1;x<=nXs;x++) {
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    num  = fabs(nAtZ_xsmz(x,sc,ms)-prvNatZ_xsmz(x,sc,ms));
                    den  = nAtZ_xsmz(x,sc,ms)+prvNatZ_xsmz(x,sc,ms)+1.0e-10;
                    cvg += sum(elem_div(num,den));
                    
                    prvNatZ_xsmz(x,sc,ms) = nAtZ_xsmz(x,sc,ms);
                }
            }
        }
        if (debug) {
            cout<<"cvg = "<<cvg<<endl;
            cout<<endl;
            cout<<"Press return to continue >> ";
            cin>>resp;
        }
    } //iteration loop
    
    if (debug){
        cout<<"Exited iteration loop after "<<itr<<" iterations."<<endl;
        cout<<"sbpr0 = "<<sbpr0<<endl;
        cout<<"Convergence value = "<<cvg<<endl;
    }
    
    //using final nAtZ's
    for (int x=1;x<=nXs;x++){
        //numbers at matetime
        if (debug) cout<<endl;
        double mateTime = 0.0;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                nAtZMT_xsmz(x,sc,ms) = mfexp(-mateTime*M_xsm(x,sc,ms))*nAtZAF_xsmz(x,sc,ms);
                if (debug) cout<<"nAtZMT("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZMT_xsmz(x,sc,ms)<<endl;
            }
        }
    }
            
    // sp biomass for final
    if (debug) cout<<endl;
    for (int x=1;x<=nXs;x++) {
        sbpr0(x)   = 0.0;//initialize for sum
        sbpr0FT(x) = 0.0;//initialize for sum
        sbpr0MT(x) = 0.0;//initialize for sum
        for(int sc=1;sc<=nSCs;sc++) {
            sbpr0(x)   += nAtZ_xsmz(x,sc,MATURE)  *wtp_xmz(x,MATURE);//note dot-product here
            sbpr0FT(x) += nAtZBF_xsmz(x,sc,MATURE)*wtp_xmz(x,MATURE);//note dot-product here
            sbpr0MT(x) += nAtZMT_xsmz(x,sc,MATURE)*wtp_xmz(x,MATURE);//note dot-product here
        }
        if (debug) {
            cout<<"sbpr0("<<x<<")   = "<<sbpr0(x)<<endl;
            cout<<"sbpr0FT("<<x<<") = "<<sbpr0FT(x)<<endl;
            cout<<"sbpr0MT("<<x<<") = "<<sbpr0MT(x)<<endl;
        }
            
    }
    
    if (dbg) {
        post<<endl;
        post<<"wtp_xmz(FEMALE,IMMATURE) = "<<wtp_xmz(FEMALE,IMMATURE)<<endl;
        post<<"wtp_xmz(FEMALE,MATURE)   = "<<wtp_xmz(FEMALE,  MATURE)<<endl;
        post<<"wtp_xmz(MALE,IMMATURE)   = "<<wtp_xmz(  MALE,IMMATURE)<<endl;
        post<<"wtp_xmz(MALE,MATURE)     = "<<wtp_xmz(  MALE,  MATURE)<<endl;
        
        for (int x=1;x<=nXs;x++) {
            post<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    post<<"M_xsm("<<x<<","<<sc<<","<<ms<<") = "<<M_xsm(x,sc,ms)<<tb<<mfexp(-midptFishery*M_xsm(x,sc,ms))<<endl;
                }
            }
        }
        for (int x=1;x<=nXs;x++) {
            post<<endl;
            post<<"prMoltImm("<<x<<") = "<<prMoltImm(x)<<endl;
            post<<"prMoltMat("<<x<<") = "<<prMoltMat(x)<<endl;
            for (int ms=1;ms<=nMSs;ms++) {
                for (int sc=1;sc<=nSCs;sc++) {
                    for (int sp=1;sp<=nSCs;sp++) post<<"tmMolt_xmssz("<<x<<","<<ms<<","<<sp<<","<<sc<<") = "<<tmMolt_xmssz(x,ms,sp,sc)<<endl;
                }
            }
        }
        for (int x=1;x<=nXs;x++) {
            post<<endl;
            post<<"prMatNS("<<x<<") = "<<prMatNS(x)<<endl;
            for (int sc=1;sc<=nSCs;sc++) {
                for (int ms=1;ms<=nMSs;ms++) {
                    for (int mp=1;mp<=nMSs;mp++) post<<"tmMat_xsmmz("<<x<<","<<sc<<","<<mp<<","<<ms<<") = "<<tmMat_xsmmz(x,sc,mp,ms)<<endl;
                }
            }
        }
        for (int x=1;x<=nXs;x++) {
            post<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    post<<"nAtZBF("<<x<<","<<sc<<","<<ms<<") = "<<nAtZBF_xsmz(x,sc,ms)<<endl;
                }
            }
        }
        writePopQuantsAfterFisheries(post);
        for (int x=1;x<=nXs;x++) {
            post<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                post<<"rec("<<x<<","<<sc<<") = "<<recPropBySex(x)*recPropByShell(sc)*recPropAtZ<<endl;
            }
        }
        for (int x=1;x<=nXs;x++) {
            post<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    post<<"nAtZ("<<x<<","<<sc<<","<<ms<<") = "<<nAtZ_xsmz(x,sc,ms)<<endl;
                }
            }
        }
        post.close();
    }
    if (debug) {
        cout<<"Finished calcSBPR0"<<endl<<endl;
        cout<<"Press 1 and hit return to continue >> ";
        cin>>resp;
    }
       
    return sbpr0MT; 

//-------------------------------------------------------------------------------------
//Calculates spawning biomass by sex per recruit for a fished population.
//Multiply by TOTAL recruitment to get spawning biomass by sex for the fished population.
//
//Inputs: 
//--tcF = directed fishing mortality on fully-selected males
//
//Output:
// if optCalcYPR = 0 (SBPRF only) 
//--dmatrix res(1,1,1,nXs) based on sex-specific biomass per recruit at mating time.
//----res(1,FEMALE) = female spawning biomass per recruit at mating time
//----res(1,  MALE) =   male spawning biomass per recruit at mating time
// if optCalcYPR = 1 (SBPRF + CPR/YPR calculations)
//--dmatrix res(1,15,1,nXs) based on sex-specific biomass per recruit, catch per recruit, and yield per recruit at mating time.
//----res( 1,FEMALE) = female spawning biomass per recruit at mating time
//----res( 1,  MALE) =   male spawning biomass per recruit at mating time
//----res( 2,FEMALE) = female retained catch per recruit  (all fisheries)
//----res( 2,  MALE) =   male retained catch per recruit  (all fisheries)
//----res( 3,FEMALE) = female total catch per recruit     (all fisheries)
//----res( 3,  MALE) =   male total catch per recruit     (all fisheries)
//----res( 4,FEMALE) = female retained catch per recruit in Tanner crab fishery
//----res( 4,  MALE) =   male retained catch per recruit in Tanner crab fishery
//----res( 5,FEMALE) = female total catch per recruit in Tanner crab fishery
//----res( 5,  MALE) =   male total catch per recruit in Tanner crab fishery
//----res( 6,FEMALE) = female total catch per recruit in snow crab fishery
//----res( 6,  MALE) =   male total catch per recruit in snow crab fishery
//----res( 7,FEMALE) = female total catch per recruit in red king crab fishery
//----res( 7,  MALE) =   male total catch per recruit in red king crab fishery
//----res( 8,FEMALE) = female total catch per recruit in groundfish trawl fishery
//----res( 8,  MALE) =   male total catch per recruit in groundfish trawl fishery
//----res( 9,FEMALE) = female retained yield per recruit (all fisheries)
//----res( 9,  MALE) =   male retained yield per recruit (all fisheries)
//----res(10,FEMALE) = female total yield per recruit    (all fisheries)
//----res(10,  MALE) =   male total yield per recruit    (all fisheries)
//----res(11,FEMALE) = female retained yield per recruit in Tanner crab fishery
//----res(11,  MALE) =   male retained yield per recruit in Tanner crab fishery
//----res(12,FEMALE) = female total yield per recruit in Tanner crab fishery
//----res(12,  MALE) =   male total yield per recruit in Tanner crab fishery
//----res(13,FEMALE) = female total yield per recruit in snow crab fishery
//----res(13,  MALE) =   male total yield per recruit in snow crab fishery
//----res(14,FEMALE) = female total yield per recruit in red king crab fishery
//----res(14,  MALE) =   male total yield per recruit in red king crab fishery
//----res(15,FEMALE) = female total yield per recruit in groundfish trawl fishery
//----res(15,  MALE) =   male total yield per recruit in groundfish trawl fishery
//
FUNCTION dmatrix calcSBPRF(double tcF, int optCalcYPR, int optTCF_EW)
    int dbg   = 0;
    int debug = 0;
    ofstream post;
    if (debug) cout<<endl<<endl<<"starting calcSBPRF("<<tcF<<")"<<endl;
    if (dbg) post.open("calcSBPRF.dat");
    
//     //numbers at size
//     d4_array nAtZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);   //on July 1
//     d4_array nAtZ01_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //just prior to fishery
//     d4_array nAtZ02_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //just after the fishery
//     d4_array nAtZ03_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //only those that will grow
//     d4_array nAtZ04_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //all after growth/no growth
//     d4_array nAtZ05_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //after maturity transition

//     d4_array nAtZFT_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //numbers just prior to fishery (same as nAtZ01)
//     d4_array nAtZMT_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs); //numbers at mating time (mateTime = time AFTER fishery occurs that mating occurs)
//         
//     d4_array prvNatZ_xsmz(1,nXs,1,nSCs,1,nMSs,1,nZs);//on July 1 from PREVIOUS iteration in while loop below
//     
//     //spawning biomass metrics
//     dvector sbprf(1,nXs);  //on July 1
//     dvector sbprfFT(1,nXs);//just prior to fishery
//     dvector sbprfMT(1,nXs);//at mating time (after fishery)
    
    //fishing mortality rates
    setFishingMortalityRates(tcF,fmSCF,fmRKF,fmGTF,optTCF_EW);
    
    //arrays for cpr, ypr calculations
    dvector totCPR;  //catch per recruit
    dvector retCPR;  //retained catch per recruit
    dvector totYPR;  //yield per recruit          (units??)
    dvector retYPR;  //retained yield per recruit (units??)
    dvector tcTotCPR;  //catch per recruit in directed Tanner crab fishery
    dvector tcRetCPR;  //retained catch per recruit in directed Tanner crab fishery
    dvector tcTotYPR;  //yield per recruit in directed Tanner crab fishery (units??)
    dvector tcRetYPR;  //retained yield per recruit in directed Tanner crab fishery (units??)
    dvector scTotCPR;  //catch per recruit in snow crab fishery
    dvector scTotYPR;  //yield per recruit in snow crab fishery (units??)
    dvector rkTotCPR;  //catch per recruit in red king crab fishery
    dvector rkTotYPR;  //yield per recruit in red king crab fishery (units??)
    dvector gtTotCPR;  //catch per recruit in groundfish trawl crab fishery
    dvector gtTotYPR;  //yield per recruit in groundfish trawl crab fishery (units??)
        
    //initialize numbers at size with recruitment
    if (debug) cout<<endl;
    for(int x=1;x<=nXs;x++) {
        for(int sc=1;sc<=nSCs;sc++) {
            nAtZ_xsmz(x,sc,IMMATURE)    = recPropBySex(x)*recPropByShell(sc)*recPropAtZ;
            nAtZ_xsmz(x,sc,  MATURE)    = 0.0;
            prvNatZ_xsmz(x,sc,IMMATURE) = 0.0;
            prvNatZ_xsmz(x,sc,  MATURE) = 0.0;
            if (debug) {
                cout<<"nAtZ("<<x<<","<<sc<<","<<IMMATURE<<") = "<<nAtZ_xsmz(x,sc,IMMATURE)<<endl;
                cout<<"nAtZ("<<x<<","<<sc<<","<<  MATURE<<") = "<<nAtZ_xsmz(x,sc,  MATURE)<<endl;
            }
        }
        if (dbg){
            nAtZ_xsmz(x,NEW_SHELL,IMMATURE) = 1.0;  //wtstest
            nAtZ_xsmz(x,NEW_SHELL,  MATURE) = 1.0;
            nAtZ_xsmz(x,OLD_SHELL,IMMATURE) = 1.0;
            nAtZ_xsmz(x,OLD_SHELL,  MATURE) = 1.0;
        }
    }
    // initial sp biomass
    for (int x=1;x<=nXs;x++) {
        sbprf(x) = 0.0;//initialize for sum
        for(int sc=1;sc<=nSCs;sc++) {
            sbprf(x) += nAtZ_xsmz(x,sc,MATURE)*wtp_xmz(x,MATURE);//note dot-product here
        }
        if (debug) cout<<"sbprf("<<x<<") = "<<sbprf(x)<<".  " ;
    }
    if (debug) cout<<endl;
        
    //run pop dy equations forward until stable equilibrium is achieved
    double cvgCri    = 1.0e-6; //criterion for convergence to stable equilibrium
    double cvg       = 100;    //initial convergence value to get into the while loop
    int itr          = 0;      //iteration counter
    while ((cvgCri<cvg)&&(itr<200)) {    //wtstest
        itr++;
        if (debug) cout<<endl<<"------------------"<<endl<<"iteration = "<<itr<<".  ";
        
        //advance numbers at size by sex
        for (int x=1;x<=nXs;x++){
            //numbers just prior to fishery: nAtZ01_xsmz
            if (debug) cout<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    nAtZBF_xsmz(x,sc,ms) = mfexp(-midptFishery*M_xsm(x,sc,ms))*nAtZ_xsmz(x,sc,ms);
                    if (debug) cout<<"nAtZBF("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZBF_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //numbers after fishery: nAtZ02_xsmz
            if (debug) cout<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    nAtZAF_xsmz(x,sc,ms) = elem_prod(totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms));
                    if (debug) cout<<"nAtZAF("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZAF_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //numbers that will molt, with change in shell condition: nAtZ03_xsmz
            if (debug) cout<<endl;
            for(int ms=1;ms<=nMSs;ms++) {
                for(int sc=1;sc<=nSCs;sc++) { //to sc
                    nAtZ03_xsmz(x,sc,ms) = 0.0;
                    for(int sp=1;sp<=nSCs;sp++) { //from sp
                        nAtZ03_xsmz(x,sc,ms) += elem_prod(tmMolt_xmssz(x,ms,sc,sp),nAtZAF_xsmz(x,sp,ms));
                    }
                    if (debug) cout<<"nAtZ03("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ03_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //numbers after molting and growth: nAtZ04_xsmz
            if (debug) cout<<endl;
            for(int ms=1;ms<=nMSs;ms++) {
                for(int sc=1;sc<=nSCs;sc++) {
                    //crabs that molted (became new shell) and grow
                    nAtZ04_xsmz(x,NEW_SHELL,ms) = nAtZ03_xsmz(x,NEW_SHELL,ms)*tmZtoZ_xzz(x);       
                    //crabs that didn't molt (and thus didn't grow)
                    nAtZ04_xsmz(x,OLD_SHELL,ms) = nAtZ03_xsmz(x,OLD_SHELL,ms);
                    if (debug) cout<<"nAtZ04("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ04_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //apply maturity transitions from mp to ms: nAtZ05_xsmz
            if (debug) cout<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) { //to ms
                    nAtZ05_xsmz(x,sc,ms) = 0.0;
                    for(int mp=1;mp<=nMSs;mp++) { //from mp
                        nAtZ05_xsmz(x,sc,ms) += elem_prod(tmMat_xsmmz(x,sc,ms,mp),nAtZ04_xsmz(x,sc,mp));
                    }
                    if (debug) cout<<"nAtZ05("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ05_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //apply end-of-year mortality: nAtZ_xsmz
            if (debug) cout<<endl;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    nAtZ_xsmz(x,sc,ms) = mfexp(-(1-midptFishery)*M_xsm(x,sc,ms))*nAtZ05_xsmz(x,sc,ms);
                    if (debug) cout<<"nAtZ("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ_xsmz(x,sc,ms)<<endl;
                }
            }
            
            //add in recruitment
            for(int sc=1;sc<=nSCs;sc++) {
                rec_xsz(x,sc) = recPropBySex(x)*recPropByShell(sc)*recPropAtZ;
                nAtZ_xsmz(x,sc,IMMATURE) += rec_xsz(x,sc);
            }
            
            if (debug) {
                cout<<endl;
                for(int sc=1;sc<=nSCs;sc++) {
                    for(int ms=1;ms<=nMSs;ms++) {
                        //post<<x<<sc<<ms<<nAtZ_xsmz(x,sc,ms)<<endl;
                        cout<<"nAtZ("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZ_xsmz(x,sc,ms)<<endl;
                    }
                }
            }
        } //sex loop
        
        //calc sp biomass
        for (int x=1;x<=nXs;x++) {
            sbprf(x) = 0.0;//initialize for sum
            for(int sc=1;sc<=nSCs;sc++) {
                sbprf(x) += nAtZ_xsmz(x,sc,MATURE)*wtp_xmz(x,MATURE);//note dot-product here
            }
            if (debug) cout<<"sbprf("<<x<<") = "<<sbprf(x)<<".  " ;
        }
        
        //check convergence and save current nAtZ
        cvg = 0.0;
        dvector num(1,nZs);
        dvector den(1,nZs);
        for (int x=1;x<=nXs;x++) {
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    num  = fabs(nAtZ_xsmz(x,sc,ms)-prvNatZ_xsmz(x,sc,ms));
                    den  = nAtZ_xsmz(x,sc,ms)+prvNatZ_xsmz(x,sc,ms)+1.0e-10;
                    cvg += sum(elem_div(num,den));
                    
                    prvNatZ_xsmz(x,sc,ms) = nAtZ_xsmz(x,sc,ms);
                }
            }
        }
        if (debug) {
            cout<<"cvg = "<<cvg<<endl;
            cout<<endl;
            cout<<"Press return to continue >> ";
            cin>>resp;
        }
    } //iteration loop
    
    if (debug){
        cout<<"Exited iteration loop after "<<itr<<" iterations."<<endl;
        cout<<"sbprf = "<<sbprf<<endl;
        cout<<"Convergence value = "<<cvg<<endl;
    }
    
    if (optCalcYPR) {
        //allocate arrays for totCPR and totYPR calculations
        totCPR.allocate(1,nXs);  //catch per recruit
        retCPR.allocate(1,nXs);  //retained catch per recruit
        totYPR.allocate(1,nXs);  //yield per recruit          (units??)
        retYPR.allocate(1,nXs);  //retained yield per recruit (units??)
        tcTotCPR.allocate(1,nXs);  //catch per recruit
        tcRetCPR.allocate(1,nXs);  //retained catch per recruit
        tcTotYPR.allocate(1,nXs);  //yield per recruit          (units??)
        tcRetYPR.allocate(1,nXs);  //retained yield per recruit (units??)
        scTotCPR.allocate(1,nXs);  //catch per recruit
        scTotYPR.allocate(1,nXs);  //yield per recruit          (units??)
        rkTotCPR.allocate(1,nXs);  //catch per recruit
        rkTotYPR.allocate(1,nXs);  //yield per recruit          (units??)
        gtTotCPR.allocate(1,nXs);  //catch per recruit
        gtTotYPR.allocate(1,nXs);  //yield per recruit          (units??)
    }
    
    //using final nAtZ's
    for (int x=1;x<=nXs;x++){
        //numbers at matetime
        if (debug) cout<<endl;
        double mateTime = 0.0;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                nAtZMT_xsmz(x,sc,ms) = mfexp(-mateTime*M_xsm(x,sc,ms))*nAtZAF_xsmz(x,sc,ms);
                if (debug) cout<<"nAtZMT("<<x<<","<<sc<<","<<ms<<")   = "<<nAtZMT_xsmz(x,sc,ms)<<endl;
            }
        }
        if (optCalcYPR) {
            totCPR(x) = 0.0;
            retCPR(x) = 0.0;
            totYPR(x) = 0.0;
            retYPR(x) = 0.0;
            tcTotCPR(x) = 0.0;
            tcRetCPR(x) = 0.0;
            tcTotYPR(x) = 0.0;
            tcRetYPR(x) = 0.0;
            scTotCPR(x) = 0.0;
            scTotYPR(x) = 0.0;
            rkTotCPR(x) = 0.0;
            rkTotYPR(x) = 0.0;
            gtTotCPR(x) = 0.0;
            gtTotYPR(x) = 0.0;
            for(int sc=1;sc<=nSCs;sc++) {
                for(int ms=1;ms<=nMSs;ms++) {
                    tcRetCatZ_xsmz(x,sc,ms) = elem_prod(elem_div(tcRetF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                        elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                    tcTotCatZ_xsmz(x,sc,ms) = elem_prod(elem_div(tcTotF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                        elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                    scCatZ_xsmz(x,sc,ms)    = elem_prod(elem_div(scF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                        elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                    rkCatZ_xsmz(x,sc,ms)    = elem_prod(elem_div(rkF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                        elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                    gtCatZ_xsmz(x,sc,ms)    = elem_prod(elem_div(gtF_xsmz(x,sc,ms),totF_xsmz(x,sc,ms)), 
                                                        elem_prod(1.0-totS_xsmz(x,sc,ms),nAtZBF_xsmz(x,sc,ms)));
                                                        
                    totCPR(x)   += sum(tcTotCatZ_xsmz(x,sc,ms)+scCatZ_xsmz(x,sc,ms)+rkCatZ_xsmz(x,sc,ms)+gtCatZ_xsmz(x,sc,ms));
                    retCPR(x)   += sum(tcRetCatZ_xsmz(x,sc,ms));
                    tcTotCPR(x) += sum(tcTotCatZ_xsmz(x,sc,ms));
                    tcRetCPR(x) += sum(tcRetCatZ_xsmz(x,sc,ms));
                    scTotCPR(x) += sum(scCatZ_xsmz(x,sc,ms));
                    rkTotCPR(x) += sum(rkCatZ_xsmz(x,sc,ms));
                    gtTotCPR(x) += sum(gtCatZ_xsmz(x,sc,ms));
                    
                    totYPR(x)   += wtp_xmz(x,ms)*(tcTotCatZ_xsmz(x,sc,ms)+scCatZ_xsmz(x,sc,ms)+rkCatZ_xsmz(x,sc,ms)+gtCatZ_xsmz(x,sc,ms));
                    retYPR(x)   += wtp_xmz(x,ms)*tcRetCatZ_xsmz(x,sc,ms);
                    tcTotYPR(x) += wtp_xmz(x,ms)*tcTotCatZ_xsmz(x,sc,ms);
                    tcRetYPR(x) += wtp_xmz(x,ms)*tcRetCatZ_xsmz(x,sc,ms);
                    scTotYPR(x) += wtp_xmz(x,ms)*scCatZ_xsmz(x,sc,ms);
                    rkTotYPR(x) += wtp_xmz(x,ms)*rkCatZ_xsmz(x,sc,ms);
                    gtTotYPR(x) += wtp_xmz(x,ms)*gtCatZ_xsmz(x,sc,ms);//note dot-product here
                }
            }
            if (debug) {
                cout<<"totCPR("<<x<<") = "<<totCPR(x)<<endl;
                cout<<"retCPR("<<x<<") = "<<retCPR(x)<<endl;
                cout<<"totYPR("<<x<<") = "<<totYPR(x)<<endl;
                cout<<"retYPR("<<x<<") = "<<retYPR(x)<<endl;
                cout<<"tcTotCPR("<<x<<") = "<<tcTotCPR(x)<<endl;
                cout<<"tcRetCPR("<<x<<") = "<<tcRetCPR(x)<<endl;
                cout<<"tcTotYPR("<<x<<") = "<<tcTotYPR(x)<<endl;
                cout<<"tcRetYPR("<<x<<") = "<<tcRetYPR(x)<<endl;
                cout<<"scTotCPR("<<x<<") = "<<scTotCPR(x)<<endl;
                cout<<"scTotYPR("<<x<<") = "<<scTotYPR(x)<<endl;
                cout<<"rkTotCPR("<<x<<") = "<<rkTotCPR(x)<<endl;
                cout<<"rkTotYPR("<<x<<") = "<<rkTotYPR(x)<<endl;
                cout<<"gtTotCPR("<<x<<") = "<<gtTotCPR(x)<<endl;
                cout<<"gtTotYPR("<<x<<") = "<<gtTotYPR(x)<<endl;
            }
        }        
    }
            
    // sp biomass for final
    if (debug) cout<<endl;
    for (int x=1;x<=nXs;x++) {
        sbprf(x)   = 0.0;//initialize for sum
        sbprfFT(x) = 0.0;//initialize for sum
        sbprfMT(x) = 0.0;//initialize for sum
        for(int sc=1;sc<=nSCs;sc++) {
            sbprf(x)   += nAtZ_xsmz(x,sc,MATURE)  *wtp_xmz(x,MATURE);//note dot-product here
            sbprfFT(x) += nAtZBF_xsmz(x,sc,MATURE)*wtp_xmz(x,MATURE);//note dot-product here
            sbprfMT(x) += nAtZMT_xsmz(x,sc,MATURE)*wtp_xmz(x,MATURE);//note dot-product here
        }
        if (debug) {
            cout<<"sbprf("<<x<<")   = "<<sbprf(x)<<endl;
            cout<<"sbprfFT("<<x<<") = "<<sbprfFT(x)<<endl;
            cout<<"sbprfMT("<<x<<") = "<<sbprfMT(x)<<endl;
        }
    }
    
    if (dbg) {
        writeBioQuants(post);
        writeFishingMortalityRates(tcF,fmSCF,fmRKF,fmGTF,post);
        writeFisheryQuants(post);
        writePopQuantsAfterFisheries(post);
        post.close();
    }
    
    dmatrix res;
    if (optCalcYPR){
        res.allocate(1,15,1,2);
        int k=1;
        res(k++) = sbprfMT;
        res(k++) = retCPR;
        res(k++) = totCPR;
        res(k++) = tcRetCPR;
        res(k++) = tcTotCPR;
        res(k++) = scTotCPR;
        res(k++) = rkTotCPR;
        res(k++) = gtTotCPR;
        res(k++) = retYPR;
        res(k++) = totYPR;
        res(k++) = tcRetYPR;
        res(k++) = tcTotYPR;
        res(k++) = scTotYPR;
        res(k++) = rkTotYPR;
        res(k++) = gtTotYPR;
    } else {
        res.allocate(1,1,1,2);
        res(1) = sbprfMT;
    }
    if (debug) {
        cout<<"Finished calcSBPRF"<<endl<<endl;
        cout<<"Press 1 and hit return to continue >> ";
        cin>>resp;
    }
        
    return res;

//-------------------------------------------------------------------------------------
FUNCTION void writeRecHSearch(ostream& post, dmatrix& resRecH)
        post<<"# recH   Fmsy    msy    Bmsy    eqR       phi      xx"<<endl;
        for (int i=resRecH.indexmin();i<=resRecH.indexmax();i++){
            post<<resRecH(i)<<tb<<resRecH(i,6)/phi0<<endl;
        }
        
//-------------------------------------------------------------------------------------
FUNCTION void writeSBPRF(ostream& post, dmatrix& mySBPRF)
    post<<"sbprf = "<<mySBPRF(1)<<endl;
    if (mySBPRF.indexmax()>1) {
        post<<endl;
        int k=2;
        post<<mySBPRF(k++)<<tb<<tb<<"retCPR   "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"totCPR   "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#tcRetCPR "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#tcTotCPR "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#scCPR    "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#rkCPR    "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#gtCPR    "<<endl;
        post<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#retYPR   "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#totYPR   "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#tcRetYPR "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#tcTotYPR "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#scYPR    "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#rkYPR    "<<endl;
        post<<mySBPRF(k++)<<tb<<tb<<"#gtYPR    "<<endl;
    }
    
//--------------------------------------------------------------------------
FUNCTION void writeBioQuants(ostream& post)
    post<<endl;
    post<<"wtp_xmz(FEMALE,IMMATURE) = "<<wtp_xmz(FEMALE,IMMATURE)<<endl;
    post<<"wtp_xmz(FEMALE,MATURE)   = "<<wtp_xmz(FEMALE,  MATURE)<<endl;
    post<<"wtp_xmz(MALE,IMMATURE)   = "<<wtp_xmz(  MALE,IMMATURE)<<endl;
    post<<"wtp_xmz(MALE,MATURE)     = "<<wtp_xmz(  MALE,  MATURE)<<endl;
    
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"M_xsm("<<x<<","<<sc<<","<<ms<<") = "<<M_xsm(x,sc,ms)<<tb<<mfexp(-midptFishery*M_xsm(x,sc,ms))<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        post<<"prMoltImm("<<x<<") = "<<prMoltImm(x)<<endl;
        post<<"prMoltMat("<<x<<") = "<<prMoltMat(x)<<endl;
        for (int ms=1;ms<=nMSs;ms++) {
            for (int sc=1;sc<=nSCs;sc++) {
                for (int sp=1;sp<=nSCs;sp++) post<<"tmMolt_xmssz("<<x<<","<<ms<<","<<sp<<","<<sc<<") = "<<tmMolt_xmssz(x,ms,sp,sc)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        post<<"prMatNS("<<x<<") = "<<prMatNS(x)<<endl;
        for (int sc=1;sc<=nSCs;sc++) {
            for (int ms=1;ms<=nMSs;ms++) {
                for (int mp=1;mp<=nMSs;mp++) post<<"tmMat_xsmmz("<<x<<","<<sc<<","<<mp<<","<<ms<<") = "<<tmMat_xsmmz(x,sc,mp,ms)<<endl;
            }
        }
    }
        
//--------------------------------------------------------------------------
FUNCTION void writePopQuantsBeforeFisheries(ostream& post)

//--------------------------------------------------------------------------
FUNCTION void writeFishingMortalityRates(double tcF, double scF, double rkF, double gtF,ostream& post)
    //fishing mortality rates    
    post<<endl<<"tcF = "<<tcF<<endl;    
    post<<endl<<"scF = "<<scF<<endl;    
    post<<endl<<"rkF = "<<rkF<<endl;    
    post<<endl<<"gtF = "<<gtF<<endl;    
    for(int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                post<<"tcTotF_East_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<tcTotF_East_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for(int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                post<<"tcTotF_West_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<tcTotF_West_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for(int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                post<<"tcTotF_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<tcTotF_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for(int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                post<<"tcRetF_East_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<tcRetF_East_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for(int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                post<<"tcRetF_West_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<tcRetF_West_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for(int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                post<<"tcRetF_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<tcRetF_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for(int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                post<<"scF_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<scF_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for(int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                post<<"rkF_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<rkF_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for(int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++)  {
            for(int ms=1;ms<=nMSs;ms++)  {
                post<<"gtF_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<gtF_xsmz(x,sc,ms)<<endl;
            }
        }
    }
        
//--------------------------------------------------------------------------
FUNCTION void writeFisheryQuants(ostream& post)
    //catches
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"tcRetCatZ_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<tcRetCatZ_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"tcTotCatZ_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<tcTotCatZ_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"scCatZ_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<scCatZ_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"rkCatZ_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<rkCatZ_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"gtCatZ_xsmz("<<x<<","<<sc<<","<<ms<<") = "<<gtCatZ_xsmz(x,sc,ms)<<endl;
            }
        }
    }

//--------------------------------------------------------------------------
FUNCTION void writePopQuantsAfterFisheries(ostream& post)
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"nAtZAF("<<x<<","<<sc<<","<<ms<<") = "<<nAtZAF_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"nAtZMT("<<x<<","<<sc<<","<<ms<<") = "<<nAtZMT_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"nAtZ03("<<x<<","<<sc<<","<<ms<<") = "<<nAtZ03_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"nAtZ04("<<x<<","<<sc<<","<<ms<<") = "<<nAtZ04_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"nAtZ05("<<x<<","<<sc<<","<<ms<<") = "<<nAtZ05_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            post<<"rec("<<x<<","<<sc<<") = "<<recPropBySex(x)*recPropByShell(sc)*recPropAtZ<<endl;
        }
    }
    for (int x=1;x<=nXs;x++) {
        post<<endl;
        for(int sc=1;sc<=nSCs;sc++) {
            for(int ms=1;ms<=nMSs;ms++) {
                post<<"nAtZ("<<x<<","<<sc<<","<<ms<<") = "<<nAtZ_xsmz(x,sc,ms)<<endl;
            }
        }
    }
    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION void writeParameterBounds(ostream& os, param_init_number& p, int toR)
    if (toR) {
        os<<p.name<<"=list("<<"type='param_init_number'"<<cc
                            <<"initial_value="<<p.initial_value<<cc
                            <<"phase="<<p.phase_start
                            <<")";
    } else {
        os<<p.name<<"phase="<<p.phase_start<<cc
                  <<"initial_value="<<p.initial_value<<endl;
    }
    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION void writeParameterBounds(ostream& os, param_init_bounded_number& p, int toR)
    if (toR) {
        os<<p.name<<"=list("<<"type='param_init_bounded_number'"<<cc
                            <<"initial_value="<<p.initial_value<<cc
                            <<"phase="<<p.phase_start<<cc
                            <<"bounds=c("<<p.minb<<cc<<p.maxb<<")"
                            <<")";
    } else {
        os<<p.name<<"phase="<<p.phase_start<<cc
                  <<"initial_value="<<p.initial_value<<cc
                  <<"bounds="<<p.minb<<cc<<p.maxb<<endl;
    }
    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION void writeParameterBounds(ostream& os, param_init_vector& p, int toR)
    if (toR) {
        os<<p.name<<"=list("<<"type='param_init_vector'"<<cc
                            <<"dims=c("<<p.indexmin()<<cc<<p.indexmax()<<")"<<cc
                            <<"initial_value="<<p.initial_value<<cc
                            <<"phase="<<p.phase_start
                            <<")";
    } else {
        os<<p.name<<"phase="<<p.phase_start<<cc
                  <<"initial_value="<<p.initial_value<<endl;
    }
    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION void writeParameterBounds(ostream& os, param_init_bounded_vector& p, int toR)
    if (toR) {
        os<<p.name<<"=list("<<"type='param_init_bounded_vector'"<<cc
                            <<"dims=c("<<p.indexmin()<<cc<<p.indexmax()<<")"<<cc
                            <<"initial_value="<<p.initial_value<<cc
                            <<"phase="<<p.phase_start<<cc
                            <<"bounds=c("<<p.minb<<cc<<p.maxb<<")"
                            <<")";
    } else {
        os<<p.name<<"phase="<<p.phase_start<<cc
                  <<"initial_value="<<p.initial_value<<cc
                  <<"bounds="<<p.minb<<cc<<p.maxb<<endl;
    }
    
// // ----------------------------------------------------------------------
// // ----------------------------------------------------------------------
// FUNCTION void writeParameterBounds(ostream& os, param_init_dev_vector& p, int toR)                        //wts: new
//     os<<p.name<<"=list("<<"type='param_init_dev_vector'"<<cc
//                         <<"dims=c("<<p.indexmin()<<cc<<p.indexmax()<<")"<<cc
//                         <<"initial_value="<<p.initial_value<<cc
//                         <<"phase="<<p.phase_start
//                         <<")";
//     
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
FUNCTION void writeParameterBounds(ostream& os, param_init_bounded_dev_vector& p, int toR)
    if (toR) {
        os<<p.name<<"=list("<<"type=param_init_bounded_dev_vector"<<cc
                            <<"dims=c("<<p.indexmin()<<cc<<p.indexmax()<<")"<<cc
                            <<"initial_value="<<p.initial_value<<cc
                            <<"phase="<<p.phase_start<<cc
                            <<"bounds=c("<<p.minb<<cc<<p.maxb<<")"
                            <<")";
    } else {
        os<<p.name<<"phase="<<p.phase_start<<cc
                  <<"initial_value="<<p.initial_value<<cc
                  <<"bounds="<<p.minb<<cc<<p.maxb<<endl;
    }
    

//-----------------------
FUNCTION openMCMCFile
    mcmc.open("TCProjMod_MCMC.r", ofstream::out|ofstream::trunc);
    mcmc<<"mcmc=list("<<endl;

//-----------------------
FUNCTION closeMCMCFile
    mcmc<<"dummy=NULL)"<<endl;//close list
    mcmc.close();//close stream
    
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
FUNCTION void writeMCMCtoR(ostream& os)
    os<<"list("; 
        writeParameterValuesToR(os); os<<","<<endl;
        writeSimToR(os,1);
    os<<"),";

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
FUNCTION void writeSimsToR(ostream& os)
    os<<"sims=list("; 
        for (int s=1;s<=nSims;s++){
            os<<"list(";
                writeSimToR(os,s); 
            os<<")"<<cc<<endl;
        }
    os<<"dummy=0)"<<endl;

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
FUNCTION void writeSimToR(ostream& os, int sim)
    os<<"sim=list(";
        os<<"sim="<<sim;//testing!!
    os<<")";

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
FUNCTION void writeParameterValuesToR(ostream& os)
    os<<"objFun="<<objFun<<cc;
    os<<"params=list(";
    os<<"pLnR0="<<pLnR0<<cc;
    os<<"pRecBeta="<<pRecBeta<<cc;
    os<<"pLnReSigma="<<pLnRecSigma<<cc;
    os<<"recR0="<<recR0<<cc;
    os<<"recH="<<recH<<cc;
    os<<"recVar="<<recVar;
    os<<")";
    
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
REPORT_SECTION
    report<<objFun     <<tb<<"#objFun"<<endl;
    report<<pLnR0      <<tb<<"#pLnR0"<<endl;
    report<<pRecBeta   <<tb<<"#pRecBeta"<<endl;
    report<<pLnRecSigma<<tb<<"#pLnReSigma"<<endl;
    report<<recR0      <<tb<<"#recR0"<<endl;
    report<<recH       <<tb<<"#recH"<<endl;
    report<<recVar     <<tb<<"#recVar"<<endl;

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
RUNTIME_SECTION
    //one number for each phase, if more phases then uses the last number
    maximum_function_evaluations 300,1000,1000,1000,1000,1000,1000,1000,3000
    convergence_criteria 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
TOP_OF_MAIN_SECTION
    arrmblsize = 1000000000;
    gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in the AUTODIF manual.
    gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);

//-------------------------------------------------------------
//-------------------------------------------------------------
FINAL_SECTION
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1) {
        closeMCMCFile();
    }
    echo<<"Model run completed"<<endl;
    cout<<"Model run completed"<<endl;
    echo.close();
