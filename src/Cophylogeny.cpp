/*
 * Cophylogeny.cpp
 *
 *  Created on: 23 Dec 2022
 *      Author: mac
 */

#include <cmath>
#include <map>
#include <iomanip>
#include <stdio.h>
#include <vector>

#include "../utility/appexception.h"
#include "../utility/approx.h"
#include "../utility/debugging.h"
#include "../utility/myrandom.h"
#include "Cophylogeny.h"

using namespace std;

extern bool _debugging;

namespace kowhai {

extern double birthRate;
extern double deathRate;
extern double hostSwitchRate;
extern double jointDuplicationProb;
extern bool _hostDictatesRate;

void Cophylogeny::coevolve()
{
	// new try
	bool _debugging(true);
	set<Node*> active;
	double t_0(H->getAge());	// max possible time of anything in the Host tree
	// Find first "active" Parasite nodes: those that can do anything
	Tree* P = *(PTrees.begin());
//	active.insert(P->getRoot());
	double firstEventTime(t_0);
	/*
	 * ==================[t_0]
	 *    +--------------* Tree P1
	 *    \-----+----+---*
	 *          |    \---*
	 *          \----+---*
	 *               \---*
	 *    +--------------* Tree P2
	 *    \--------------*
	 *    ^t_min
	 */
	map<double, set<Node*>> hNodeTimes;
	for (Tree* P : PTrees) {
		Node* r = P->getRoot();
		double time = r->getTime();
		DEBUG(cout << "Time of parasite root " << r->getLabel() << " = " << time << endl);
		if (time < firstEventTime) {
			firstEventTime = time;
			active.clear();
			active.insert(r);
		} else if (time == firstEventTime) {
			active.insert(r);
		}
		hNodeTimes[time].insert(r);
	}
	DEBUG(cout << "COEVOLVE" << endl);
	DEBUG(
		cout << "Initial active set of parasites: { ";
		for (Node *p : active) {
			cout << p->getLabel() << '@' << p->getTime() << ' ';
		}
		cout << "}" << endl;
	);
	Node* h;
	H->getRoot()->storeNodeTimes(hNodeTimes);
	DEBUG(
		for (auto pr : hNodeTimes) {
			cout << pr.first << endl;
		}
	);
	hNodeTimes[firstEventTime].insert(active.begin(), active.end());
	double eventRate;
	set<Node*> availableHosts;
	unsigned int numCodivergences(0);
	unsigned int numDuplications(0);
	unsigned int numExtinctions(0);
	unsigned int numHostSwitches(0);
	unsigned int numJointDuplicationEvents(0);
	map<int, int> duplicationSizes;
	unsigned int numLineageSortingEvents(0);
	double t(firstEventTime);
	approx roughlyEqual(0.000001);	// functor to return true iff input numbers are within this tolerance of each other.
	for (map<double, set<Node*>>::iterator moment = hNodeTimes.begin(); moment != hNodeTimes.end(); ++moment) {
		auto nextMoment = moment;
		DEBUG(cout << "Current host time = " << moment->first << endl);
		++nextMoment;	// need the *next* moment in time
		double nextHostTimePoint = nextMoment->first;
		while (t < nextHostTimePoint) {
			DEBUG(cout << "\tcurrent time t=" << t << "; nextHostTime=" << nextHostTimePoint << endl);
			// get total rate:
			eventRate = 0.0;
			availableHosts.clear();
			if (_hostDictatesRate) {
				Node* p = *(active.begin());
				eventRate += birthRate;
				eventRate += deathRate;
				if (p->getHost()->hasParent()) {
					eventRate += hostSwitchRate;
				}
				for (Node* p : active) {
					availableHosts.insert(p->getHost());	// all the extant hosts
				}
			} else {
				for (Node* p : active) {
					eventRate += p->getBirthRate();
					eventRate += p->getDeathRate();
					if (p->getHost()->hasParent()) {
						eventRate += p->getHostSwitchRate();
					}
					availableHosts.insert(p->getHost());	// all the extant hosts
				}
			}
			DEBUG(
				cout << "\tActive set of parasites: { ";
				for (Node *p : active) {
					cout << p->getLabel() << ' ';
				}
				cout << "}" << endl;
			);
			DEBUG(cout << "\tTotal event rate (duplication + host switch + extinction) = " << eventRate << endl);
			double t_next = -log10(fran()) / eventRate;	// t_next is the amount of time from t to the next parasite event
			DEBUG(cout << "\tNext event would be at time " << (t+t_next) << endl);
			std::set<Node*> toActivate;
			std::set<Node*> toDeactivate;
			if (t + t_next > nextHostTimePoint) {
				t = nextHostTimePoint;
				DEBUG(cout << "\t\tThis is a HOST NODE-driven event at time " << t << endl);
				// CODIVERGENCE and LINEAGE SORTING
				for (Node* p: active) {
					h = p->getHost();
					if (h->isLeaf()) {
						continue;
					}
					DEBUG(cout << "\t\tp = " << p->getLabel() << "; h = " << h->getLabel() << endl);
//					DEBUG(cout << "\t\ttime of p = " << t << "; time of host = " << h->getTime() << endl);
					DEBUG(cout << "\t\tnext node(s): { ");
					DEBUG(
						for (Node* x : moment->second) {
							cout << x->getLabel() << ' ';
						}
						cout << '}' << endl;
					)
					if (moment->second.count(p)) {	// else could codiverge parasites early!
						hNodeTimes[t].erase(p);
						if (p->doesCodiverge()) {
							DEBUG(cout << "\t\tCODIVERGENCE of " << p->getLabel() << " with its host " << h->getLabel() << endl);
							p->codivergeWith(h);
							p->setEvent(codivergence);
							toDeactivate.insert(p);
							for (Node* c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
								toActivate.insert(c);
								hNodeTimes[c->getHost()->getTime()].insert(c);
	//							c->setEvent(codivergence);
							}
							++numCodivergences;
						} else {
							Node* nuHost = h->getFirstChild();
							if (fran() < 0.5) {
								nuHost = nuHost->getSibling();
							}
							p->setEvent(duplication);
							DEBUG(cout << "\t\tLINEAGE SORTING of " << p->getLabel() << " onto nascent host " << nuHost->getLabel() << endl);
							p->setHost(nuHost);
							p->onHostVertex() = true;
							p->setTime(nuHost->getTime());
							nuHost->addParasite(p);
							hNodeTimes[nuHost->getTime()].insert(p);
							++numLineageSortingEvents;
						}
					} else {
						DEBUG(cout << "\t\tThis p=" << p->getLabel() << " cannot codiverge or lineage sort as its host doesn't have the right time." << endl);
					}
				}
				if (roughlyEqual(t, H->getAge())) {
					break;
				}
			} else {
				// DUPLICATION, EXTINCTION, and HOST SWITCHING
				DEBUG(cout << "\t\tThis is an event that is INDEPENDENT of host events" << endl);
				Node *p = getRandomElement<Node*>(active);
				double pB = p->getBirthRate();
				double pS = (p->hasParent()) ? p->getHostSwitchRate() : 0.0;
				double pX = p->getDeathRate();
				double pTotal = pB + pS + pX;
				double ran = dran(pTotal);
				t += t_next;	// advance the current time
				if (ran < pB) {	// DUPLICATION
					h = p->getHost();
					DEBUG(cout << "\t\tDUPLICATION of " << p->getLabel() << " on its host " << h->getLabel() << endl);
					p->bifurcate();
					p->setTime(t);
					p->setEvent(duplication);
					toDeactivate.insert(p);	//active.erase(p);
					for (Node* c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
						toActivate.insert(c);	//active.insert(c);
	//					c->setEvent(duplication);
						c->setHost(h);
						c->setTime(t);
						hNodeTimes[h->getTime()].insert(c);
					}
					hNodeTimes[h->getTime()].erase(p);
					++numDuplications;
					// Now do "hitch-hiking" type duplication events, currently selected at uniform random:
					bool _jointEvent(false);
					int segmentalDupSize(1);
					for (Node* q : active) {
						if (p == q || (q->getHost() != h)) {
							continue;
						}
						if (dran() < jointDuplicationProb) {
							_jointEvent = true;
							q->bifurcate();
							q->setTime(t);
							q->setEvent(duplication);
							toDeactivate.insert(q);	//active.erase(p);
							for (Node* c = q->getFirstChild(); c != nullptr; c = c->getSibling()) {
								toActivate.insert(c);	//active.insert(c);
			//					c->setEvent(duplication);
								c->setHost(h);
								c->setTime(t);
								hNodeTimes[h->getTime()].insert(c);
							}
							hNodeTimes[h->getTime()].erase(q);
							++numDuplications;
							++segmentalDupSize;
						}
					}
					numJointDuplicationEvents += (_jointEvent) ? 1 : 0 ;
					duplicationSizes[segmentalDupSize] = duplicationSizes[segmentalDupSize] + 1;
				} else if (ran < pB + pS) {	// HOST SWITCH
					h = p->getHost();
					availableHosts.erase(h);
					if (availableHosts.size() == 0) {
						throw new app_exception("Cannot do a host switch, as there are no available host lineages!");
					}
					p->bifurcate();
					p->setTime(t);
					toDeactivate.insert(p);	//active.erase(p);
					Node* c  = p->getFirstChild();
	//				c->setEvent(duplication);
					toActivate.insert(c);	//active.insert(c);
					c->setHost(h);
					c->setTime(h->getTime());
					hNodeTimes[h->getTime()].insert(c);
					DEBUG(cout << "\t\t(duplication+) HOST SWITCH of " << p->getLabel() << " on its host " << h->getLabel()
							<< ", with " << c->getLabel() << " staying on host " << h->getLabel() << " and ");
					c = c->getSibling();
	//				c->setEvent(hostswitch);
					toActivate.insert(c);	//active.insert(c);
					Node *nuHost = getRandomElement<Node*>(availableHosts);
					c->setHost(nuHost);
					c->setTime(nuHost->getTime());
					hNodeTimes[nuHost->getTime()].insert(c);
					hNodeTimes[h->getTime()].erase(p);
					++numHostSwitches;
					DEBUG(cout << c->getLabel() << " jumping to host " << nuHost->getLabel() << endl);
				} else { 	// DEATH
					p->setTime(t);
	//				p->setEvent(death);
					toDeactivate.insert(p);	//active.erase(p);
					hNodeTimes[h->getTime()].erase(p);
					++numExtinctions;
					DEBUG(cout << "\t\tDEATH of " << p->getLabel() << " on its host " << h->getLabel() << endl);
				}
			}
			DEBUG(cout << *this);
			for (Node* kill : toDeactivate) {
				active.erase(kill);
			}
			for (Node* act : toActivate) {
				active.insert(act);
			}
		}
	}
	cout << "numCodivergences = " << numCodivergences << endl;
	cout << "numDuplication Events = " << numDuplications << endl;
	cout << "numExtinctions = " << numExtinctions << endl;
	cout << "numHostSwitches = " << numHostSwitches << endl;
	cout << "numLineageSortingEvents = " << numLineageSortingEvents << endl;
	cout << "numJointDuplicationEvents = " << numJointDuplicationEvents << endl;
	cout << "dupSize count" << endl;
	for (auto pr : duplicationSizes) {
		cout << std::setw(7) << pr.first << ' ' << pr.second << endl;
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
	occupantsAtTime[h->getTime()].insert(p);
//	if (_onVertex) {
//		occupantsAtTime[h->getTime()].insert(p);
//	}
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
	occupantsAtTime[t].insert(p);
//	occupantsAtTime[h->getTime()].insert(q);	// XXX needs much testing!
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
	os << "Host tree:" << endl << *(C.getHostTree());// << C.getHostTree()->details();
	C.storeAssociationInfo();
	for (Tree* P : C.getParasiteTrees()) {
		P->setShowInfo(true);
		os << P->getLabel() << ":" << endl << *P;// << P->details();
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
				str += eventSymbol[p->getEvent()];
//				if (p->onHostVertex()) {
//					str += eventSymbol[codivergence];
//				} else {
//					str += eventSymbol[duplication];
//				}
				str += "]" + p->getHost()->getLabel();
			}
			DEBUG(cout << "info on " << p->getLabel() << ": " << str << endl);
			(*info)[p] = str;
		}
		P->setInfo(info);
	}
}

} /* namespace kowhai */
