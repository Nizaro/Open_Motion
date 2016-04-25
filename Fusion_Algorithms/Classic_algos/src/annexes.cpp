/*
 * annexes.cpp
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */



#include "annexes.h"

namespace om{


///////////////////////////////////////////////////////
/////          Tools display console              /////
///////////////////////////////////////////////////////





/* Convert double to string with specified number of places after the decimal
   and left padding. */
std::string prd(const double x, const int decDigits, const int width) {
    stringstream ss;
    ss << fixed << right;
    ss.fill(' ');        // fill space around displayed #
    ss.width(width);     // set  width around displayed #
    ss.precision(decDigits); // set # places after decimal
    ss << x;
    return ss.str();
}


/* Convert double to string with specified number of places after the decimal. */
std::string prd(const double x, const int decDigits) {
    stringstream ss;
    ss << fixed;
    ss.precision(decDigits); // set # places after decimal
    ss << x;
    return ss.str();
}

/*! Center-aligns string within a field of width w. Pads with blank spaces
    to enforce alignment. */
std::string center(const string s, const int w) {
    stringstream ss, spaces;
    int padding = w - s.size();                 // count excess room to pad
    for(int i=0; i<padding/2; ++i)
        spaces << " ";
    ss << spaces.str() << s << spaces.str();    // format with padding
    if(padding>0 && padding%2!=0)               // if odd #, add 1 space
        ss << " ";
    return ss.str();
}


///////////////////////////////////////////////////////
/////                 Others                      /////
///////////////////////////////////////////////////////


double signeOf(double x){

	return abs(x)<EPSILON ? 1.0 : x < 0.0 ? -1.0 : 1.0 ;

}




}/* namespace isf */
