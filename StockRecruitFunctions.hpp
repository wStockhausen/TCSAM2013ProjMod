/* 
 * File:   StockRecruitFunctions.hpp
 * Author: william.stockhausen
 *
 * Created on August 2, 2013, 7:00 PM
 */

#ifndef STOCKRECRUITFUNCTIONS_HPP
#define	STOCKRECRUITFUNCTIONS_HPP

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
double calcBevertonHolt(double R0, double h, double phi0, double spB);

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
dvariable calcBevertonHolt(dvariable& R0, dvariable& h, double phi0, double spB);

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
double calcRicker(double R0, double h, double phi0, double spB);

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
dvariable calcRicker(dvariable& R0, dvariable& h, double phi0, double spB);

//-------------------------------------------------------------------------------------
//Calculates recruitment level for spawning stock biomass spB.
//
//Inputs:
//--R0   = recruitment level for unfished stock
//--h    = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--spB  = spawning stock biomass
//--srType = flag indicating stock recruit function (SRTYPE_RICKER, SRTYPE_BEVHOLT, or SRTYPE_CONSTANT)
//
//Output:
//--R at spB based on the stock recruit function
double calcSRFunction(double R0, double h, double phi0, double spB, int srType);

//-------------------------------------------------------------------------------------
//Calculates recruitment level for spawning stock biomass spB.
//
//Inputs:
//--R0   = recruitment level for unfished stock
//--h    = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--spB  = spawning stock biomass
//--srType = flag indicating stock recruit function (SRTYPE_RICKER, SRTYPE_BEVHOLT, or SRTYPE_CONSTANT)
//
//Output:
//--R at spB based on the stock recruit function
dvariable calcSRFunction(dvariable& R0, dvariable& h, double phi0, double spB, int srType);
    
//-------------------------------------------------------------------------------------
//Calculates equilibrium recruitment level for stock fished such that 
//spawning stock biomass/recruit equals phi.
//
//Inputs:
//--R0   = recruitment level for unfished stock
//--h    = steepness for stock-recruit function
//--phi0 = spawning stock biomass/recruit for unfished stock
//--xx   = spawning stock size relative to the unfished size (or spawning stock biomass/recruit relative to unfished size) in range [0,1]
//--srType = flag indicating stock recruit function (SRTYPE_RICKER, SRTYPE_BEVHOLT, or SRTYPE_CONSTANT)
//
//Output:
//--equilibrium R (NOT bias-adjusted)
double calcEqRec(double R0, double h, double phi0, double xx, int srType);

#endif	/* STOCKRECRUITFUNCTIONS_HPP */

