/*
 * =====================================================================================
 *
 *       Filename:  time_utils.cpp
 *
 *    Description:  Source code for time util functions.
 *
 *        Version:  1.0
 *        Created:  15/04/16 20:38:30
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include "time.h"
#include "time_utils.hpp"

static struct timespec startTime, endTime;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  diff
 *    Arguments:  stuct timespec start - The start time.
 *                stuct timespec end - The end time.
 *      Returns:  A timespec struct containing the time difference.
 *  Description:  Gets the difference in the time elapsed since start and end. Returns a
 *                time struct.
 * =====================================================================================
 */

struct timespec diff(struct timespec start, struct timespec end) {
		struct timespec temp;
			if ((end.tv_nsec-start.tv_nsec)<0) {
				temp.tv_sec = end.tv_sec-start.tv_sec-1;
				temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
			} else {
				temp.tv_sec = end.tv_sec-start.tv_sec;
				temp.tv_nsec = end.tv_nsec-start.tv_nsec;
			}
		return temp;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  startClock
 *  Description:  Starts the clock and saves time to startTime.
 * =====================================================================================
 */

void startClock() {
	clock_gettime(CLOCK_REALTIME, &startTime) ;
}		/* -----  end of function startClock  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  stopClock
 *  Description:  Stops the clock and saves time to endTime.
 * =====================================================================================
 */

void stopClock() {
	clock_gettime(CLOCK_REALTIME, &endTime) ;
}		/* -----  end of function stopClock  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getElapsedTime
 *      Returns:  The difference between the start time and end time.
 *  Description:  Returns the time difference between start time (startClock) and end time
 *                (stopClock).
 * =====================================================================================
 */

float getElapsedTime() {
	return (float)diff(startTime , endTime).tv_sec + (float) diff(startTime,endTime).tv_nsec/1E9 ;	
}		/* -----  end of function getElapsedTime  ----- */

