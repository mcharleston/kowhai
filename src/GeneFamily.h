/*
 * GeneFamily.h
 *
 *  Created on: 21 Dec 2022
 *      Author: mac
 */

#ifndef GENEFAMILY_H_
#define GENEFAMILY_H_

#include "Tree.h"

namespace cospec {

class GeneFamily {
private:
	Tree* H;	// associated host tree
public:
	GeneFamily();
	virtual ~GeneFamily();
};

} /* namespace segdup */

#endif /* GENEFAMILY_H_ */
