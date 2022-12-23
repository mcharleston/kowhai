/*
 * Cophylogeny.cpp
 *
 *  Created on: 23 Dec 2022
 *      Author: mac
 */

#include <map>
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

	bool _debugging(false);
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
		DEBUG(cout << "Processing host node " + h->getLabel() << " at time " << h->getTime() << endl);
		remove.clear();
		if (h->isLeaf()) {
			DEBUG(cout << "Reached the leaves: time to stop." << endl);
			break;
		}
		DEBUG(cout << "This host node has " << h->getParasites().size() << " parasites." << endl);
		for (auto a : h->getParasites()) {
			Node* p = a.second;
			if (p->onHostVertex()) {
				if (p->doesCodiverge()) {
					DEBUG(cout << p->getLabel() << " codiverges with this host" << endl);
					p->codivergeWith(h);
					p->getFirstChild()->onHostVertex() = true;
					p->getFirstChild()->getSibling()->onHostVertex() = true;
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
	// Store the association information in the parasite tree(s):
	storeAssociationInfo();
}

Node* Cophylogeny::createParasiteRoot(Node* h, bool _onVertex) {
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
	PTrees.insert(P);
	cout << "Tree P with address " << P << " has been added to the cophylogeny" << endl;
	return p;
}

ostream& operator<<(ostream& os, Cophylogeny& C) {
	os << *(C.getHostTree());
	for (Tree* P : C.getParasiteTrees()) {
		P->setShowInfo(true);
		os << *P;
	}
	return os;
}

void Cophylogeny::storeAssociationInfo() {
	bool _debugging(false);
	for (Tree* P : PTrees) {
		P->gatherVertices();
		DEBUG(cout << "|V(P)| = " << P->getVertices().size() << endl);
		map<Node*, string>* info = P->getInfo();
		if (info == nullptr) {
			info = new map<Node*, string>();
		}
		for (auto as : P->getVertices()) {
			Node* p = as.second;
			DEBUG(cout << "parasite node " << p->getLabel() << endl);
			string str;
			if (p->isLeaf()) {
				str = as.first + ":" + p->getHost()->getLabel();
			} else {
				str = as.first + ":[";
				if (p->onHostVertex()) {
					str += eventSymbol[codivergence];
				} else {
					str += eventSymbol[duplication];
				}
				str += "]" + p->getHost()->getLabel();
			}
			DEBUG(cout << "info on " << p->getLabel() << ": " << str << endl);
			(*info)[p] = str;
		}
		P->setInfo(info);
	}
}

} /* namespace kowhai */
