#ifndef TIME_UTILS_H_NXEDSLDY
#define TIME_UTILS_H_NXEDSLDY

/*
 * =====================================================================================
 *
 *       Filename:  time_utils.hpp
 *
 *    Description:  File containing useful timing functions.
 *
 *        Version:  1.0
 *        Created:  17/02/16 11:02:23
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include "time.h"

struct timespec diff(struct timespec start, struct timespec end) ;
void startClock() ;
void stopClock() ;
float getElapsedTime() ;

#endif /* end of include guard: TIME_UTILS_H_NXEDSLDY */
