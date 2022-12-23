/*
 * Cophylogeny.cpp
 *
 *  Created on: 23 Dec 2022
 *      Author: mac
 */

#include <stdio.h>
#include <vector>

#include "../utility/debugging.h"
#include "../utility/myrandom.h"
#include "Cophylogeny.h"

using namespace std;

extern bool _debugging;

namespace kowhai {

void Cophylogeny::coevolve() {
	/**
	 * Use the same process as for growing a Yule tree but with additional actions:
	 * 		at each bifurcation of a node in H, call codiverge to see if the parasites also codiverge
	 * 		on each edge of H, check for next events:
	 * 			any of the associated parasites duplicating, going extinct or host-switching
	 * 			another bifurcation on that host edge
	 */

	bool _debugging(true);
	double currentTime(0.0);

	// Host tree nodes are sorted in time order from t=0 to t=t_max (the leaves)
	Node *h;
	int eventIndex(0);
	vector<Node*> orderedHostNodes = H->getOrderedNodes();
	int i(0);
	vector<string> remove;
	while (currentTime < H->getAge()) {
		// collect all parasites on the host nodes extant at this time point
		h = orderedHostNodes[i];
		DEBUG(cout << "Processing host node " + h->getLabel() << endl);
		remove.clear();
		if (h->isLeaf()) {
			DEBUG(cout << "Reached the leaves: time to stop." << endl);
			break;
		}
		for (auto a : h->getParasites()) {
			Node* p = a.second;
			if (p->onHostVertex()) {
				if (p->doesCodiverge()) {
					DEBUG(cout << p->getLabel() << " is codiverging on this host" << endl);
					p->codivergeWith(h);
				} else {
					// select one nascent host lineage at uniform random:
					Node* nuHost = h->getFirstChild();
					if (fran() < 0.5) {
						nuHost = nuHost->getSibling();
					}
					DEBUG(cout << p->getLabel() << " misses the boat and goes down new host lineage " << nuHost->getLabel() << endl);
					// put p on it & note it for removal from this host node (h):
					p->setHost(nuHost);
					p->onHostVertex() = true;
					nuHost->addParasite(p);
					remove.push_back(p->getLabel());
				}
			}
		}
		for (string& str : remove) {
			h->getParasites().erase(str);
		}
		++i;
	}
//	std::vector<Node*> L;
//	L.clear();
//	Node *root = new Node();
//	root->setLabel("r" + root->getLabel().erase(0,1));
//	root->setTree(H);
//	H->setRoot(root);
//	L.push_back(root);
//	int numLeaves = L.size();
//	double divTime = 0;
//	double height = 0;
//	if (targetNumLeaves < 2) {
//		return;
//	}
////	int vertexNumber(1);
//	while (numLeaves < targetNumLeaves) {
//		// Choose a leaf at random:
//		int idx = iran(L.size());
//		Node* x = L[idx];
//		// Bifurcate it with no branch lengths assigned:
//		x->bifurcate();
//		Node* y = x->getFirstChild();
//		Node *z = y->getSibling();
//		L[idx] = y; // replacing x in the array list
//		L.push_back(z);
//		numLeaves++;
//		// Add divergence time to all leaf branch lengths:
//		divTime = -log(fran()) / (birthRate * (double) numLeaves);
//		height += divTime;
//		for (Node* l : L) {
//			l->addToBranchLength(divTime);
//		}
//	}
//	// get a new divergence time for the last period:
//	divTime = -log(fran()) / (birthRate * (double) numLeaves);
//	// Now account for sampling: we assume we pick the tree any time during
//	// the period when there are this number of leaves:
//	double lastDivTime = divTime * fran();
//	height += lastDivTime;
//	for (Node* l : L) {
//		l->addToBranchLength(divTime);
//	}
//	gatherVertices();
}

void Cophylogeny::createParasiteRoot(Node* h, bool _onVertex) {
	/**
	 * Not sure now how to associate parasites with edges or vertices.
	 *
	 * The idea of a CophyMap (see segdup) handles it by keeping track of the event in the map information,
	 * but here I think I need to either
	 * 		keep a flag to say "the edge above this vertex" or "the vertex itself"
	 * 	or
	 * 		create a base "Location" class from which both Edge (which doesnt' currently exist as a class) and
	 * 		Node derive.  I think not, baby puppy.
	 */
	bool _debugging(true);
	DEBUG(cout << "Creating Parasite Root node on host ";
		if (_onVertex) {
			cout << "vertex ";
		} else {
			cout << "edge above ";
		}
		cout << h->getLabel() << endl);
	Node* p = new Node();
	Tree* P = new Tree(p);
	p->setHost(h);
	h->addParasite(p);
	p->onHostVertex() = _onVertex;
}

} /* namespace kowhai */
