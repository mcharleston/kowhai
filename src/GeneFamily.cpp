/*
 * GeneFamily.cpp
 *
 *  Created on: 21 Dec 2022
 *      Author: mac
 */

#include "../utility/appexception.h"
#include "GeneFamily.h"

namespace kowhai {

GeneFamily::GeneFamily() {
	// TODO Auto-generated constructor stub

}

GeneFamily::~GeneFamily() {
	// TODO Auto-generated destructor stub
}

void GeneFamily::addDependentTree(Tree *P, Node *h) {
	if (GF.count(P->getLabel())) {
		throw new app_exception("GeneFamily::addDependentTree(): cannot add this tree named "
			+ P->getLabel() + " as one of this name is already in the collection.");
	}
	GF[P->getLabel()] = P;
	P->getRoot()->setHost(h);
}

void GeneFamily::coevolve() {
	/**
	 * General Process:
	 * Begin with a complete Host tree H
	 * From the root of H proceed forwards in time toward the leaves.
	 * At each host node h_i, check all occupants / dependents for codivergence
	 *
	 */
}
} /* namespace kowhai */
