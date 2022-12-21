/*
 * debugging.h
 *
 *  Created on: 11 Dec. 2018
 *      Author: mc7
 */

#ifndef DEBUGGING_H_
#define DEBUGGING_H_

#if defined DEBUGGING
#define DEBUG(A) { if (_debugging) { A; } }
#else
#define DEBUG(A) { }
#endif


#endif /* DEBUGGING_H_ */
