/*
 * GeneFamily.h
 *
 *  Created on: 21 Dec 2022
 *      Author: mac
 */

#ifndef GENEFAMILY_H_
#define GENEFAMILY_H_

#include <map>
#include <string>

#include "Tree.h"

namespace kowhai {

class GeneFamily {
	/**
	 * Strictly speaking this isn't a gene family because you can create multiple gene families to go in here,
	 * each represented by a tree whose root may be associated anywhere in H, but it's either that or call it
	 * a GeneSet, which doesn't quite cut it.
	 * Also note that this is just one aspect of the cophylogeny problem: this class could also be used to keep
	 * track of multiple non-necessarily-independent parasites on a group of host species, say.
	 */
private:
	Tree* H;	// associated host tree
	std::map<std::string,Tree*> GF;
public:
	GeneFamily();
	virtual ~GeneFamily();


};

} /* namespace segdup */

#endif /* GENEFAMILY_H_ */
