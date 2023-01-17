/*
 * Cophylogeny.h
 *
 *  Created on: 23 Dec 2022
 *      Author: mac
 */

#ifndef SRC_COPHYLOGENY_H_
#define SRC_COPHYLOGENY_H_

#include <map>
#include <set>

#include "Tree.h"

namespace kowhai {

class Cophylogeny {
private:
	Tree* H;	// the underlying host tree
	std::set<Tree*> PTrees;
//	std::map<int, std::set<Node*> > NodeByIndex;	// I need to come up with a better name for this
		// NodeByIndex is the map of eventIndex -> <all Nodes at that timepoint>
	std::map<double, std::set<Node*> > pNodeAtTime;
	int currentTimeIndex;
public:
	Cophylogeny() : H(nullptr), currentTimeIndex(0) {}
	virtual ~Cophylogeny() {}

	void coevolve();
	void cleverCoevolve();
	void correctCoevolve();

	Node* createParasiteRoot(Node *h, bool _onVertex);
	Node* createParasiteRoot(Node *h, double beforeTime);

//	void duplicate(Node* p, double t);

	inline Tree* getHostTree() { return H; }
	inline std::set<Tree*>& getParasiteTrees() { return PTrees; }

	void outputForSegdup(std::ostream &os);

	inline void setHostTree(Tree *T) { H = T; }
	void storeAssociationInfo();
};

std::ostream& operator<<(std::ostream& os, Cophylogeny &C);

} /* namespace kowhai */

#endif /* SRC_COPHYLOGENY_H_ */
