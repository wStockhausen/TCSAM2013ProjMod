/* 
 * File:   ModelConstants.hpp
 * Author: william.stockhausen
 *
 * Created on March 12, 2013, 7:12 AM
 */

#pragma once
#ifndef MODELCONSTANTS_HPP
    #define	MODELCONSTANTS_HPP

    extern const adstring cc;//comma
    extern const adstring qt;//quote
    extern const adstring tb;//tab

    extern const int OFF;
    extern const int ON;
    extern const adstring STR_OFF;
    extern const adstring STR_ON;
        
    extern const int INT_FALSE;
    extern const int INT_TRUE;
    extern const adstring STR_FALSE;
    extern const adstring STR_TRUE;
    
    extern const int IMMATURE;
    extern const int MATURE;
    extern const int ALL_MATURITY;
    extern const adstring STR_IMMATURE;
    extern const adstring STR_MATURE;
    extern const adstring STR_ALL_MATURITY;
    
    extern const int FEMALE;
    extern const int MALE;
    extern const int ALL_SEXES;
    extern const adstring STR_FEMALE;
    extern const adstring STR_MALE;
    extern const adstring STR_ALL_SEXES;
        
    extern const int NEW_SHELL;
    extern const int OLD_SHELL;
    extern const int ALL_SHELL;
    extern const adstring STR_NEW_SHELL;
    extern const adstring STR_OLD_SHELL;
    extern const adstring STR_ALL_SHELL;
    
    extern const adstring STR_CONSTANT;
    extern const adstring STR_BEVHOLT;
    extern const adstring STR_RICKER;
    extern const int SRTYPE_CONSTANT;
    extern const int SRTYPE_BEVHOLT;
    extern const int SRTYPE_RICKER;
//----------------------------------------------------------------------
//          consts
//----------------------------------------------------------------------
    class ModelConsts {
    public:
        static int getOnOffType(adstring s);
        static adstring getOnOffType(int i);
        
        static int getBooleanType(adstring s);
        static adstring getBooleanType(int i);
        
        static int getMaturityType(adstring s);
        static adstring getMaturityType(int i);
        
        static int getSexType(adstring s);
        static adstring getSexType(int i);
        
        static int getShellType(adstring s);
        static adstring getShellType(int i);
        
        //Stock-recruit function types
        static int getSRType(adstring s);
        static adstring getSRType(int i);
        
        static const adstring STR_VAR;
        static const adstring STR_STDV;
        static const adstring STR_CV;
        const static int SCLTYPE_VAR  = 0;//variances are given
        const static int SCLTYPE_STDV = 1;//std devs are given
        const static int SCLTYPE_CV   = 2;//cv's are given
        static int getScaleType(adstring sclType);
        static adstring getScaleType(int sclFlg);
        
        static double convertToStdDev(double sclVal, double mnVal, int sclFlg);
        static dvector convertToStdDev(dvector sclVal, dvector mnVal, int sclFlg);
};


#endif	/* MODELCONSTANTS_HPP */

