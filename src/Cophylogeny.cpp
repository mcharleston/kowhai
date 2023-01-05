/*
 * Cophylogeny.cpp
 *
 *  Created on: 23 Dec 2022
 *      Author: mac
 */

#include <cmath>
#include <map>
#include <stdio.h>
#include <vector>

#include "../utility/appexception.h"
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
	vector<Node*> remove;
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
		// DUPLICATION, HOST SWITCHING and EXTINCTION:

		// CODIVERGENCE and MISSING THE BOAT:
		for (Node* p : h->getParasites()) {
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
					remove.push_back(p);
				}
			}
		}
		for (Node* p : remove) {
			h->getParasites().erase(p);
		}
		++i;
	}
	// Store the association information in the parasite tree(s):
	storeAssociationInfo();
}

void Cophylogeny::cleverCoevolve() {
	bool _debugging(true);
//	DEBUG(
//			for (auto pr : H->getVertices()) {
//				Node* h = pr.second;
//				cout << "\ttime of vertex " << h->getLabel() << " is " << h->getTime() << endl;
//			}
//			cout << (*H)
//	);
	double t = pNodeAtTime.begin()->first;	// the first event in time since map is *ordered*
	DEBUG(cout << "First event at time t = " << t << endl);
	double t_final = H->getAge();	// this is a time from 0 to the present
	DEBUG(cout << "Last event at time t = " << t_final << endl);
	double eventRate;
	vector<Node*> parasToRemove;
	set<Node*> availableHosts;
	for (std::map<double, std::set<Node*> >::iterator moment = pNodeAtTime.begin();
			moment != pNodeAtTime.end(); ++moment) {// currentPs : pNodeAtTime) {
		DEBUG(cout << "Time point: " << t << endl);
		DEBUG(cout << "|pNodeAtTime| = " << pNodeAtTime.size() << endl);
		set<Node*> occupants = moment->second;
//		set<Node*> occupants = currentPs.second;
		DEBUG(
				cout << "Occupants at this point: {";
				for (Node* p : occupants) {
					cout << " " << p->getLabel();
				}
				cout << " }" << endl;
		);
		Node* n = *(occupants.begin());
		Node* h = n->getHost();
		if (n->onHostVertex()) {
			DEBUG(cout << "This is on the host vertex " << h->getLabel() << endl);
			// do codivergence and lineage sorting
			for (Node* p : occupants) {
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
					parasToRemove.push_back(p);
				}
			}
			for (Node* q : parasToRemove) {
				h->getParasites().erase(q);	// may not need to do this.
			}
		} else {
			DEBUG(cout << "This is on the edge above " << h->getLabel() << endl);
			// check all active parasite nodes for other events
			eventRate = 0.0;
			availableHosts.clear();
			double nextHostEventTime(0.0);
			for (Node* p : occupants) {
				eventRate += p->getBirthRate();
				eventRate += p->getDeathRate();
				if (p->getHost()->hasParent()) {
					eventRate += p->getHostSwitchRate();
				}
				nextHostEventTime = min(nextHostEventTime, p->getHost()->getTime());
				availableHosts.insert(p->getHost());	// all the extant hosts
			}
			DEBUG(cout << "Time of next host node = " << nextHostEventTime << endl);
			DEBUG(cout << "Total event rate = " << eventRate << endl);
			double t_next = -log10(fran()) / eventRate;
			DEBUG(cout << "Next event is at time " << (t+t_next) << endl);
			if (t + t_next > nextHostEventTime) {
				DEBUG(cout << "Setting to the next time point, being a host node" << endl);
				t = nextHostEventTime;
				continue;	// no more events before the next host node
			}
			t += t_next;	// advance the current time
			DEBUG(cout << "Selecting the eventful lineage: ");
			// now select which node is going to do something:
			Node *q = getRandomElement<Node*>(occupants);
			DEBUG(cout << q->getLabel() << endl);
			// which event? :
			double pB = q->getBirthRate();
			double pS = (q->hasParent()) ? q->getHostSwitchRate() : 0.0;
			double pX = q->getDeathRate();
			double pTotal = pB + pS + pX;
			double ran = dran(pTotal);
//			pNodeAtTime[t].erase(q);	// this p-node won't be available for any more events at this time // XXX probably unnecessary to remove this
			if (ran < pB) {
				// DUPLICATION
				/*
				 *  	bifurcate p
				 *  	take this p off h and replace it with its two children
				 *  	add these two child ps to pNodeAtTime with the new time
				 */
				DEBUG(cout << "DUPLICATION event:" << endl);
				h = q->getHost();
				Tree* T = q->getTree();
				DEBUG(
						cout << "Tree before:" << endl << *T << endl;
				);
				q->bifurcate();
				DEBUG(
						cout << "Tree AFTER:" << endl << *T << endl;
				);
				q->setEvent(duplication);
				q->setTime(t);
				q->onHostVertex() = false;
				h->getParasites().erase(q);
				for (Node * c = q->getFirstChild(); c != nullptr; c = c->getSibling()) {
					c->setHost(h);
					h->getParasites().insert(c);
					c->setTime(t);
					c->onHostVertex() = false;
					pNodeAtTime[t].insert(c);	// now the edges descendant from q are available for events at the next timestep
					DEBUG(cout << "New parasite node " << c->getLabel() << " on host " << h->getLabel() << " at time " << t << endl);
					DEBUG(cout << "\t which has parent " << q->getLabel() << endl);
				}
				DEBUG(cout << "First child of parasite node " << q->getLabel() << " is " << q->getFirstChild()->getLabel() << endl);
			} else if (ran < pB + pS) {
				// HOST SWITCH
				DEBUG(cout << "HOST SWITCH event:" << endl);
				h = q->getHost();
				availableHosts.erase(h);
				if (availableHosts.size() == 0) {
					throw new app_exception("Cannot do a host switch, as there are no available host lineages!");
				}
				q->bifurcate();
				q->setEvent(duplication);
				q->setTime(t);
				q->onHostVertex() = false;
				h->getParasites().erase(q);
				Node *c = q->getFirstChild();
				c->setHost(h);
				h->getParasites().insert(c);
				c->setTime(t);
				c->onHostVertex() = false;
				pNodeAtTime[t].insert(c);
				c = c->getSibling();
				Node *nuhost = getRandomElement<Node*>(availableHosts);
				c->setHost(h);
				h->getParasites().insert(c);
				c->setTime(t);
				c->onHostVertex() = false;
				pNodeAtTime[t].insert(c);
			} else {
				// DEATH
				DEBUG(cout << "EXTINCTION (LOSS) event:" << endl);
				q->setTime(t);
				q->onHostVertex() = false;
				q->setEvent(loss);
				h->getParasites().erase(q);
			}
//			DEBUG(cout << *this);
//			DEBUG(q->getTree()->gatherVertices());
//			DEBUG(q->getTree()->gatherLeaves());
			DEBUG(cout << *(q->getTree()) << endl);
			DEBUG(cout << "|pNodeAtTime| = " << pNodeAtTime.size() << endl);
		}
		if (t >= t_final) {
			break;
		}
	}
}

Node* Cophylogeny::createParasiteRoot(Node* h, bool _onVertex) {
	/**
	 * Not sure now how to associate parasites with edges or vertices.
	 *
	 * The idea of a CophyMap (see segdup) handles it by keeping track of the event in the map information,
	 * but here I think I need to either
	 * 		keep a flag to say "the edge above this vertex" or "the vertex itself"
	 * 	or
	 * 		create a base "Location" class from which both Edge (which doesn't currently exist as a class) and
	 * 		Node derive.  I think not, baby puppy.
	 */
	bool _debugging(false);
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
	if (_onVertex) {
		pNodeAtTime[h->getTime()].insert(p);
	}
	PTrees.insert(P);
	DEBUG(cout << "Tree P with address " << P << " has been added to the cophylogeny" << endl);
	return p;
}
Node* Cophylogeny::createParasiteRoot(Node *h, double beforeBy) {
	bool _debugging(true);
// put a parasite node on this host node but back by <beforeTime>. Check boundaries.
	if (h->hasParent()) { // check "beforeBy" isn't longer than the branch length
		if (beforeBy >= h->getBranchLength()) {
			throw new app_exception("Node::createParasiteRoot(Node *h, double beforeBy): you are trying to put a parasite node too far back.");
		}
	}
	Node* p = createParasiteRoot(h, false);
	double t(h->getTime() - beforeBy);
	p->setTime(t);
	DEBUG(cout << "Setting parasite root at time t = " << t << " on host tree" << endl);
//	Node *q = new Node();
//	p->setFirstChild(q);
//	q->setHost(h);
//	q->onHostVertex() = false;
	pNodeAtTime[t].insert(p);
//	pNodeAtTime[h->getTime()].insert(q);	// XXX needs much testing!
	return p;
}

//void Cophylogeny::duplicate(Node* p, double t) {
//	// check that p has exactly one child node
//	if (p->isLeaf()) {
//		throw new app_exception("Attempting to duplicate but the parasite node needs a single child, but it's a leaf.");
//	}
//	Node* q = p->getFirstChild();
//	if (q->getSibling() != nullptr) {
//		throw new app_exception("Attempting to duplicate but the parasite node needs a single child, but has at least two.");
//	}
//	q->onHostVertex() = false;
//	q->setTime(t);
//	Node* h = p->getHost();
//	h->getParasites().erase(p);
//	h->getParasites().insert(q);	// XXX still debating whether I really need a map of the occupants; wouldn't a set of Node* do?
//}

void Cophylogeny::outputForSegdup(ostream& os) {
	os << "-S \"";
	H->writeNewick(os);
	os << "\"";
	for (auto P : PTrees) {
		os << " -G \"";
		P->writeNewick(os);
		os << "\" \"";
		for (map<string, Node*>::iterator iter = P->getLeaves().begin(); iter != P->getLeaves().end(); ) {
			Node* p = iter->second;
			os << p->getLabel() << ':' << p->getHost()->getLabel();
			++iter;
			if (iter != P->getLeaves().end()) {
				os << " ";
			}
		}
		os << "\"";
	}
	os << endl;
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
