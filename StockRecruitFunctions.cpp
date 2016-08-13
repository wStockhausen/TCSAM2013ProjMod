#include <admodel.h>
#include "ModelConstants.hpp"
#include "StockRecruitFunctions.hpp"

//-------------------------------------------------------------------------------------
//calculates recruitment level based on a Beverton-Holt stock recruit function
//
//Inputs:
//--R0   = recruitment level for unfished stock
//--h    = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--spB  = spawning stock biomass
//
//Output:
//--recruitment level (value of the stock-recruit function)
double calcBevertonHolt(double R0, double h, double phi0, double spB){
    double alpha = 0.8*R0*h/(h-0.2);
    double beta  = 0.2*R0*phi0*(1.0-h)/(h-0.2);
    double rec = alpha*spB/(beta+spB);
    return rec;
}

//-------------------------------------------------------------------------------------
//calculates recruitment level based on a Beverton-Holt stock recruit function
//
//Inputs:
//--R0   = recruitment level for unfished stock
//--h    = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--spB  = spawning stock biomass
//
//Output:
//--recruitment level (value of the stock-recruit function)
dvariable calcBevertonHolt(dvariable& R0, dvariable& h, double phi0, double spB){
    RETURN_ARRAYS_INCREMENT();
    dvariable alpha = 0.8*R0*h/(h-0.2);
    dvariable beta  = 0.2*R0*phi0*(1.0-h)/(h-0.2);
    dvariable rec = alpha*spB/(beta+spB);
    RETURN_ARRAYS_DECREMENT();
    return rec;
}

//-------------------------------------------------------------------------------------
//calculates recruitment level based on a Beverton-Holt stock recruit function
//
//Inputs:
//--R0   = recruitment level for unfished stock
//--h    = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--spB  = spawning stock biomass
//
//Output:
//--recruitment level (value of the stock-recruit function)
double calcRicker(double R0, double h, double phi0, double spB){
    double alpha = pow(5*h,1.25)/phi0;
    double beta  = 0.2*R0*phi0*(1.0-h)/(h-2.0);
    double rec = alpha*spB*mfexp(-beta*spB);
    return rec;
}

//-------------------------------------------------------------------------------------
//calculates recruitment level based on a Beverton-Holt stock recruit function
//
//Inputs:
//--R0   = recruitment level for unfished stock
//--h    = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--spB  = spawning stock biomass
//
//Output:
//--recruitment level (value of the stock-recruit function)
dvariable calcRicker(dvariable& R0, dvariable& h, double phi0, double spB){
    RETURN_ARRAYS_INCREMENT();
    dvariable alpha = pow(5*h,1.25)/phi0;
    dvariable beta  = 0.2*R0*phi0*(1.0-h)/(h-2.0);
    dvariable rec = alpha*spB*mfexp(-beta*spB);
    RETURN_ARRAYS_DECREMENT();
    return rec;
}

//-------------------------------------------------------------------------------------
//Calculates recruitment level for spawning stock biomass = xx*spB0.
//
//Inputs:
//--R0   = recruitment level for unfished stock
//--h    = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--spB  = spawning stock biomass
//--srType = flag indicating stock recruit function (SRTYPE_RICKER, SRTYPE_BEVHOLT, or SRTYPE_CONSTANT)
//
//Output:
//--R at S based on the stock recruit function
double calcSRFunction(double R0, double h, double phi0, double spB, int srType){
    double rec;
    if (srType==SRTYPE_RICKER) {
        rec = calcBevertonHolt(R0,h,phi0,spB);
    } else if (srType==SRTYPE_BEVHOLT) {
        rec = calcRicker(R0,h,phi0,spB);
    } else if (srType==SRTYPE_CONSTANT) {
        rec = R0;                                     // assumption of a constant recruitment
    }
    return rec;
}

//-------------------------------------------------------------------------------------
//Calculates recruitment level for spawning stock biomass = xx*spB0.
//
//Inputs:
//--R0 (dvariable)   = recruitment level for unfished stock
//--h  (dvariable)   = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--spB  = spawning stock biomass
//--srType = flag indicating stock recruit function (SRTYPE_RICKER, SRTYPE_BEVHOLT, or SRTYPE_CONSTANT)
//
//Output:
//--R at S based on the stock recruit function
dvariable calcSRFunction(dvariable& R0, dvariable& h, double phi0, double spB, int srType){
    RETURN_ARRAYS_INCREMENT();
    dvariable rec;
    if (srType==SRTYPE_RICKER) {
        rec = calcBevertonHolt(R0,h,phi0,spB);
    } else if (srType==SRTYPE_BEVHOLT) {
        rec = calcRicker(R0,h,phi0,spB);
    } else if (srType==SRTYPE_CONSTANT) {
        rec = R0;                                     // assumption of a constant recruitment
    }
    RETURN_ARRAYS_DECREMENT();
    return rec;
}
    
//-------------------------------------------------------------------------------------
//Calculates equilibrium recruitment level for stock fished such that 
//spawning stock biomass/recruit equals xx*phi0.
//
//Inputs:
//--R0   = recruitment level for unfished stock
//--h    = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--xx   = spawning stock per recruit at F relative to the unfished size (or spawning stock biomass/recruit relative to unfished size)
//--srType = flag indicating stock recruit function (SRTYPE_RICKER, SRTYPE_BEVHOLT, or SRTYPE_CONSTANT)
//
//Output:
//--equilibrium R (NOT bias-adjusted)
double calcEqRec(double R0, double h, double phi0, double xx, int srType){
    double R_eq;
    if (srType==SRTYPE_RICKER) {
        double alpha = pow(5*h,1.25)/phi0;
        double beta  = 0.2*R0*phi0*(1.0-h)/(h-2.0);
        R_eq = -(log(1/(alpha*xx*phi0)))/(xx*phi0*beta);// Ricker formulation of equilibrium recruitment at each F
    } else if (srType==SRTYPE_BEVHOLT) {
        double alpha = 0.8*R0*h/(h-0.2);
        double beta  = 0.2*R0*phi0*(1.0-h)/(h-0.2);
        R_eq = (alpha*(xx*phi0) - beta)/(xx*phi0);     // BH formulation of equilibrium recruitment at F corresponding to xx
    } else if (srType==SRTYPE_CONSTANT) {
        R_eq = R0;                                     // assumption of a constant recruitment
    }
    return R_eq;
}

