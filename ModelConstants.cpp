#include <admodel.h>
#include "ModelConstants.hpp"

const adstring cc = ", ";
const adstring qt = "'";
const adstring tb = "  ";

const int OFF = 0;
const int ON  = 1;
const adstring STR_OFF = "OFF";
const adstring STR_ON  = "ON";

const int INT_FALSE = 0;
const int INT_TRUE  = 1;
const adstring STR_FALSE = "FALSE";
const adstring STR_TRUE  = "TRUE";

const int NEW_SHELL = 1;
const int OLD_SHELL = 2;
const int ALL_SHELL = 3;
const adstring STR_NEW_SHELL = "NEW_SHELL";
const adstring STR_OLD_SHELL = "OLD_SHELL";
const adstring STR_ALL_SHELL = "ALL_SHELL";

const int FEMALE    = 1;
const int MALE      = 2;
const int ALL_SEXES = 3;
const adstring STR_FEMALE    = "FEMALE";
const adstring STR_MALE      = "MALE";
const adstring STR_ALL_SEXES = "ALL_SEXES";

const int IMMATURE     = 1;
const int MATURE       = 2;
const int ALL_MATURITY = 3;
const adstring STR_IMMATURE     = "IMMATURE";
const adstring STR_MATURE       = "MATURE";
const adstring STR_ALL_MATURITY = "ALL_MATURITY";

const int SRTYPE_CONSTANT = 0;
const int SRTYPE_BEVHOLT  = 1;
const int SRTYPE_RICKER   = 2;
const adstring STR_CONSTANT = "CONSTANT";
const adstring STR_BEVHOLT  = "BEVERTONHOLT";
const adstring STR_RICKER   = "RICKER";
//----------------------------------------------------------------------
//          ModelConsts
//----------------------------------------------------------------------
int ModelConsts::getOnOffType(adstring s){
    s.to_upper();
    if (s==STR_OFF) return OFF;
    if (s==STR_ON ) return ON;
    return -1;
}
adstring ModelConsts::getOnOffType(int i){
    if (i==OFF) return STR_OFF;
    if (i==ON ) return STR_ON;
    return adstring("");
}

int ModelConsts::getBooleanType(adstring s){
    s.to_upper();
    if (s==STR_FALSE) return INT_FALSE; else
    if (s==STR_TRUE)  return INT_TRUE;  else
    cout<<"Unrecognized BooleanType '"<<s<<"'"<<endl;
    return -1;
}
adstring ModelConsts::getBooleanType(int i){
    if (i==INT_FALSE) return STR_FALSE; else
    if (i==INT_TRUE)  return STR_TRUE;  else
    cout<<"Unrecognized BooleanType '"<<i<<"'"<<endl;
    return adstring("");
}

int ModelConsts::getMaturityType(adstring s){
    s.to_upper();
    if (s==STR_IMMATURE)     return IMMATURE;     else
    if (s==STR_MATURE)       return MATURE;       else
    if (s==STR_ALL_MATURITY) return ALL_MATURITY; else
    cout<<"Unrecognized MaturityType '"<<s<<"'"<<endl;
    return 0;
}
adstring ModelConsts::getMaturityType(int i){
    if (i==IMMATURE)     return STR_IMMATURE;     else
    if (i==MATURE)       return STR_MATURE;       else
    if (i==ALL_MATURITY) return STR_ALL_MATURITY; else
    cout<<"Unrecognized MaturityType '"<<i<<"'"<<endl;
    return adstring("");
}

int ModelConsts::getSexType(adstring s){
    s.to_upper();
    if (s==STR_FEMALE)    return FEMALE;    else
    if (s==STR_MALE)      return MALE;      else
    if (s==STR_ALL_SEXES) return ALL_SEXES; else
    cout<<"Unrecognized SexType '"<<s<<"'"<<endl;
    return 0;
}
adstring ModelConsts::getSexType(int i){
    if (i==FEMALE)    return STR_FEMALE;    else
    if (i==MALE)      return STR_MALE;      else
    if (i==ALL_SEXES) return STR_ALL_SEXES; else
    cout<<"Unrecognized SexType '"<<i<<"'"<<endl;
    return adstring("");
}

int ModelConsts::getShellType(adstring s){
    s.to_upper();
    if (s==STR_NEW_SHELL) return NEW_SHELL; else
    if (s==STR_OLD_SHELL) return OLD_SHELL; else
    if (s==STR_ALL_SHELL) return ALL_SHELL; else
    cout<<"Unrecognized ShellType '"<<s<<"'"<<endl;
    return 0;
}
adstring ModelConsts::getShellType(int i){
    if (i==NEW_SHELL) return STR_NEW_SHELL; else
    if (i==OLD_SHELL) return STR_OLD_SHELL; else
    if (i==ALL_SHELL) return STR_ALL_SHELL; else
    cout<<"Unrecognized ShellType '"<<i<<"'"<<endl;
    return adstring("");
}

int ModelConsts::getSRType(adstring s){
    s.to_lower();
    if (s==STR_CONSTANT) return SRTYPE_CONSTANT;
    if (s==STR_BEVHOLT)  return SRTYPE_BEVHOLT;
    if (s==STR_RICKER)   return SRTYPE_RICKER;
    return 0;
}
adstring ModelConsts::getSRType(int i){
    if (i==SRTYPE_CONSTANT) return STR_CONSTANT;
    if (i==SRTYPE_BEVHOLT)  return STR_BEVHOLT;
    if (i==SRTYPE_RICKER)   return STR_RICKER;
    return adstring("");
}

const adstring ModelConsts::STR_VAR  = "variance";
const adstring ModelConsts::STR_STDV = "std_dev";
const adstring ModelConsts::STR_CV   = "cv";
int ModelConsts::getScaleType(adstring sclType){
    if (sclType==STR_VAR)  return SCLTYPE_VAR;
    if (sclType==STR_STDV) return SCLTYPE_STDV;
    if (sclType==STR_CV)   return SCLTYPE_CV;
    return 0;
}
adstring ModelConsts::getScaleType(int sclFlg){
    if (sclFlg==SCLTYPE_VAR)  return STR_VAR;
    if (sclFlg==SCLTYPE_STDV) return STR_STDV;
    if (sclFlg==SCLTYPE_CV)   return STR_CV;
    return adstring("");
}

/***************************************************************
*   Converts from variance type indicated by sclFlg to         *
*   standard deviation.                                        *
***************************************************************/
double ModelConsts::convertToStdDev(double sclVal, double mnVal, int sclFlg){
    double sdv = 0.0;
    if (sclFlg==SCLTYPE_CV) {
        sdv = sclVal*mnVal;
    } else if (sclFlg==SCLTYPE_STDV) {
        sdv = sclVal;
    } else if (sclFlg==SCLTYPE_VAR) {
        sdv = sqrt(sclVal);
    }
    return sdv;
}

/***************************************************************
*   Converts from variance type indicated by sclFlg to         *
*   standard deviation.                                        *
***************************************************************/
dvector ModelConsts::convertToStdDev(dvector sclVal, dvector mnVal, int sclFlg){
    dvector sdv = 0.0*sclVal;
    if (sclFlg==SCLTYPE_CV) {
        sdv = elem_prod(sclVal,mnVal);
    } else if (sclFlg==SCLTYPE_STDV) {
        sdv = sclVal;
    } else if (sclFlg==SCLTYPE_VAR) {
        sdv = sqrt(sclVal);
    }
    return sdv;
}

