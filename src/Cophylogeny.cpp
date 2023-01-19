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
#include "../utility/approx.h"
#include "../utility/debugging.h"
#include "../utility/myrandom.h"
#include "Cophylogeny.h"

using namespace std;

extern bool _debugging;

namespace kowhai {

//void Cophylogeny::coevolve() {
//	/**
//	 * Use the same process as for growing a Yule tree but with additional actions:
//	 * 		at each bifurcation of a node in H, call codiverge to see if the parasites also codiverge
//	 * 		on each edge of H, check for next events:
//	 * 			any of the associated parasites duplicating, going extinct or host-switching
//	 * 			another bifurcation on that host edge
//	 */
//
//	bool _debugging(false);
//	double currentTime(0.0);
//
//	// Host tree nodes are sorted in time order from t=0 to t=t_max (the leaves)
//	Node *h;
////	int eventIndex(0);
//	vector<Node*> orderedHostNodes = H->getOrderedNodes();
//	int i(0);
//	vector<Node*> remove;
//	while (currentTime < H->getAge()) {
//		// collect all parasites on the host nodes extant at this time point
//		h = orderedHostNodes[i];
//		DEBUG(cout << "Processing host node " + h->getLabel() << " at time " << h->getTime() << endl);
//		remove.clear();
//		if (h->isLeaf()) {
//			DEBUG(cout << "Reached the leaves: time to stop." << endl);
//			break;
//		}
//		DEBUG(cout << "This host node has " << h->getParasites().size() << " parasites." << endl);
//		// DUPLICATION, HOST SWITCHING and EXTINCTION:
//
//		// CODIVERGENCE and MISSING THE BOAT:
//		for (Node* p : h->getParasites()) {
//			if (p->onHostVertex()) {
//				if (p->doesCodiverge()) {
//					DEBUG(cout << p->getLabel() << " codiverges with this host" << endl);
//					p->codivergeWith(h);
//					p->getFirstChild()->onHostVertex() = true;
//					p->getFirstChild()->getSibling()->onHostVertex() = true;
//				} else {
//					// select one nascent host lineage at uniform random:
//					Node* nuHost = h->getFirstChild();
//					if (fran() < 0.5) {
//						nuHost = nuHost->getSibling();
//					}
//					DEBUG(cout << p->getLabel() << " misses the boat and goes down new host lineage " << nuHost->getLabel() << endl);
//					// put p on it & note it for removal from this host node (h):
//					p->setHost(nuHost);
//					p->onHostVertex() = true;
//					nuHost->addParasite(p);
//					remove.push_back(p);
//				}
//			}
//		}
//		for (Node* p : remove) {
//			h->getParasites().erase(p);
//		}
//		++i;
//	}
//	// Store the association information in the parasite tree(s):
//	storeAssociationInfo();
//}
//
//void Cophylogeny::cleverCoevolve() {
//	bool _debugging(true);
////	DEBUG(
////			for (auto pr : H->getVertices()) {
////				Node* h = pr.second;
////				cout << "\ttime of vertex " << h->getLabel() << " is " << h->getTime() << endl;
////			}
////			cout << (*H)
////	);
//	double t = occupantsAtTime.begin()->first;	// the first event in time since map is *ordered*
//	DEBUG(cout << "First event at time t = " << t << endl);
//	double t_final = H->getAge();	// this is a time from 0 to the present
//	DEBUG(cout << "Last event at time t = " << t_final << endl);
//	double eventRate;
//	vector<Node*> parasToRemove;
//	set<Node*> availableHosts;
//	DEBUG(cout << "HOST TREE:" << endl << *H);
//	for (std::map<double, std::set<Node*> >::iterator moment = occupantsAtTime.begin();
//			moment != occupantsAtTime.end(); ++moment) {// currentPs : pNodeAtTime) {
//		double t_final = H->getAge();	// this is a time from 0 to the present
//		DEBUG(cout << "Time point: " << moment->first << endl);
//		DEBUG(
//				cout << "Nodes by time:" << endl;
//				for (auto pr : occupantsAtTime) {
//					cout << "\t" << pr.first << ": { ";
//					for (Node* p : pr.second) {
//						cout << p->getLabel() << " ";
//					}
//					cout <<"}" << endl;
//				}
//		);
////		DEBUG(cout << "|pNodeAtTime| = " << pNodeAtTime.size() << endl);
//		set<Node*> occupants = moment->second;
////		set<Node*> occupants = currentPs.second;
//		DEBUG(
//				cout << "Occupants at this point: {";
//				for (Node* p : occupants) {
//					cout << " " << p->getLabel();
//				}
//				cout << " }" << endl;
//		);
//		Node* n = *(occupants.begin());
//		Node* h = n->getHost();
//		set<Node*> nuNodes;
//		if (n->onHostVertex()) {
//			if (h->isLeaf()) {
//				continue;
//			}
//			// CODIVERGENCE or LINEAGE SORTING
//			DEBUG(cout << "This is on the host vertex " << h->getLabel() << endl);
//			double nextHostEventTime(t_final + 1.0);
//			// do codivergence and lineage sorting
//			for (Node* p : occupants) {
//				if (p->doesCodiverge()) {
//					DEBUG(cout << p->getLabel() << " codiverges with this host" << endl);
//					p->codivergeWith(h);
//					for (Node *c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
//						c->onHostVertex() = true;
//						nuNodes.insert(c);
//						c->setTime(c->getHost()->getTime());
//						nextHostEventTime = min(nextHostEventTime, c->getTime());
//						DEBUG(cout << "\t" << c->getLabel() <<":" << c->getHost()->getLabel() << "; t=" << c->getTime() << endl);
//					}
//				} else {
//					// select one nascent host lineage at uniform random:
//					Node* nuHost = h->getFirstChild();
//					if (fran() < 0.5) {
//						nuHost = nuHost->getSibling();
//					}
//					DEBUG(cout << p->getLabel() << " misses the boat and goes down new host lineage " << nuHost->getLabel() << endl);
//					// put p on it & note it for removal from this host node (h):
//					p->setHost(nuHost);
//					p->onHostVertex() = true;
//					nuHost->addParasite(p);
//					nuNodes.insert(p);
//					nextHostEventTime = min(nextHostEventTime, nuHost->getTime());
//				}
//			}
//			for (Node* nu : nuNodes) {
//				occupantsAtTime[nu->getHost()->getTime()].insert(nu);
//			}
//			t = nextHostEventTime;
////			for (Node* q : parasToRemove) {
////				h->getParasites().erase(q);	// may not need to do this.
////			}
//		} else {
//			DEBUG(cout << "This is on the edge above " << h->getLabel() << endl);
//			// check all active parasite nodes for other events
//			eventRate = 0.0;
//			availableHosts.clear();
//			double nextHostEventTime(0.0); // TODO This is going to always put the next time to 0!!!
//			for (Node* p : occupants) {
//				eventRate += p->getBirthRate();
//				eventRate += p->getDeathRate();
//				if (p->getHost()->hasParent()) {
//					eventRate += p->getHostSwitchRate();
//				}
//				nextHostEventTime = min(nextHostEventTime, p->getHost()->getTime());
//				availableHosts.insert(p->getHost());	// all the extant hosts
//			}
//			DEBUG(cout << "Time of next host node = " << nextHostEventTime << endl);
//			DEBUG(cout << "Total event rate = " << eventRate << endl);
//			double t_next = -log10(fran()) / eventRate;
//			DEBUG(cout << "Next event is at time " << (t+t_next) << endl);
//			if (t + t_next > nextHostEventTime) {
//				DEBUG(cout << "Setting to the next time point, being a host node" << endl);
//				throw new app_exception("There is more to do here: because the next event time is after the next host node, I need to set up for codivergence events next.");
//				t = nextHostEventTime;
//				continue;	// no more events before the next host node
//			}
//			t += t_next;	// advance the current time
//			DEBUG(cout << "Selecting the eventful lineage: ");
//			// now select which node is going to do something:
//			Node *q = getRandomElement<Node*>(occupants);
//			DEBUG(cout << q->getLabel() << endl);
//			// which event? :
//			double pB = q->getBirthRate();
//			double pS = (q->hasParent()) ? q->getHostSwitchRate() : 0.0;
//			double pX = q->getDeathRate();
//			double pTotal = pB + pS + pX;
//			double ran = dran(pTotal);
////			pNodeAtTime[t].erase(q);	// this p-node won't be available for any more events at this time // XXX probably unnecessary to remove this
//			if (ran < pB) {
//				// DUPLICATION
//				/*
//				 *  	bifurcate p
//				 *  	take this p off h and replace it with its two children
//				 *  	add these two child ps to pNodeAtTime with the new time
//				 */
//				DEBUG(cout << "DUPLICATION event:" << endl);
//				h = q->getHost();
//				Tree* T = q->getTree();
//				DEBUG(
//						cout << "Tree before:" << endl << *T << endl;
//				);
//				q->bifurcate();
//				DEBUG(
//						cout << "Tree AFTER:" << endl << *T << endl;
//				);
//				q->setEvent(duplication);
//				q->setTime(t);
//				q->onHostVertex() = false;
//				h->getParasites().erase(q);
//				DEBUG(cout << "Checking children of newly duplicated " << q->getLabel() << endl);
//				for (Node * c = q->getFirstChild(); c != nullptr; c = c->getSibling()) {
//					c->setHost(h);
//					h->getParasites().insert(c);
//					c->setTime(t);
//					c->onHostVertex() = false;
//					occupantsAtTime[t].insert(c);	// now the edges descendant from q are available for events at the next timestep
//					DEBUG(cout << "New parasite node " << c->getLabel() << " on host " << h->getLabel() << " at time " << t << endl);
//					DEBUG(cout << "\t which has parent " << q->getLabel() << endl);
//				}
//				DEBUG(cout << "First child of parasite node " << q->getLabel() << " is " << q->getFirstChild()->getLabel() << endl);
//			} else if (ran < pB + pS) {
//				// HOST SWITCH
//				DEBUG(cout << "HOST SWITCH event:" << endl);
//				h = q->getHost();
//				availableHosts.erase(h);
//				if (availableHosts.size() == 0) {
//					throw new app_exception("Cannot do a host switch, as there are no available host lineages!");
//				}
//				q->bifurcate();
//				q->setEvent(duplication);
//				q->setTime(t);
//				q->onHostVertex() = false;
//				h->getParasites().erase(q);
//				Node *c = q->getFirstChild();
//				c->setHost(h);
//				h->getParasites().insert(c);
//				c->setTime(t);
//				c->onHostVertex() = false;
//				occupantsAtTime[t].insert(c);
//				c = c->getSibling();
//				Node *nuhost = getRandomElement<Node*>(availableHosts);
//				c->setHost(nuhost);
//				nuhost->getParasites().insert(c);
//				c->setTime(t);
//				c->onHostVertex() = false;
//				occupantsAtTime[t].insert(c);
//			} else {
//				// DEATH
//				DEBUG(cout << "EXTINCTION (LOSS) event:" << endl);
//				q->setTime(t);
//				q->onHostVertex() = false;
//				q->setEvent(death);
//				h->getParasites().erase(q);
//			}
////			DEBUG(cout << *this);
////			DEBUG(q->getTree()->gatherVertices());
////			DEBUG(q->getTree()->gatherLeaves());
//			DEBUG(cout << "|pNodeAtTime| = " << occupantsAtTime.size() << endl);
//		}
//		DEBUG(cout << *this);
////		if (t >= t_final) {
////			DEBUG(cout << "This is at or beyond the age of the host tree so coevolution is stopping." << endl);
////			break;
////		}
//	}
//}
//
//void Cophylogeny::correctCoevolve() {
//	bool _debugging(true);
//	DEBUG(
//			for (auto pr : H->getVertices()) {
//				Node* h = pr.second;
//				cout << "\ttime of vertex " << h->getLabel() << " is " << h->getTime() << endl;
//			}
//			cout << (*H)
//	);
//	double t = occupantsAtTime.begin()->first;	// the first event in time since map is *ordered*
//	H->initialiseOccupants(occupantsAtTime);
//	DEBUG(cout << "First event at time t = " << t << endl);
//	double t_final = H->getAge();	// this is a time from 0 to the present
//	DEBUG(cout << "Last event at time t = " << t_final << endl);
//	double eventRate;
//	vector<Node*> parasToRemove;
//	set<Node*> availableHosts;
//	map<double, set<Node*>> hNodeTimes;
//	H->getRoot()->storeNodeTimes(hNodeTimes);
//	auto currentHosts = hNodeTimes.begin();
//	auto lastHostNodes = hNodeTimes.end();
//	double t_max = lastHostNodes->first;	// this is the time of the last host node.
////	Node* nextHost = hNodes.begin().second.begin();
//	DEBUG(cout << "HOST TREE:" << endl << *H);
//	set<Node*> nuNodes;
//	cout << "|pNodes| = " << occupantsAtTime.size() << endl;
//	for (std::map<double, std::set<Node*> >::iterator moment = occupantsAtTime.begin();
//			moment != occupantsAtTime.end(); ++moment) {
//		cout << "t=" << moment->first << "\t; pNodes={";
//		for (Node* p : moment->second) {
//			cout << " " << p->getLabel();
//		}
//		cout << " }" << endl;
//	}
//	for (std::map<double, std::set<Node*> >::iterator moment = occupantsAtTime.begin();
//			moment != occupantsAtTime.end(); ++moment) {
//
//		auto nextMoment = moment;
//		++nextMoment;	// need the *next* moment in time
//		double nextHostTimePoint = nextMoment->first;
//		set<Node*>& occupants = moment->second;
//		while (t < nextHostTimePoint) {
//			DEBUG(
//					cout << "Occupants at this point: {";
//					for (Node* p : occupants) {
//						cout << " " << p->getLabel();
//					}
//					cout << " }" << endl;
//			);
//			Node* n = *(occupants.begin());
//			Node* h = n->getHost();
//			DEBUG(cout << "This is on the edge above " << h->getLabel() << endl);
//			// check all active parasite nodes for other events
//			eventRate = 0.0;
//			availableHosts.clear();
//			for (Node* p : occupants) {
//				eventRate += p->getBirthRate();
//				eventRate += p->getDeathRate();
//				if (p->getHost()->hasParent()) {
//					eventRate += p->getHostSwitchRate();
//				}
//				availableHosts.insert(p->getHost());	// all the extant hosts
//			}
//			DEBUG(cout << "Time of next host node = " << nextHostTimePoint << endl);
//			DEBUG(cout << "Total event rate = " << eventRate << endl);
//			double t_next = -log10(fran()) / eventRate;	// t_next is the amount of time from t to the next parasite event
//			DEBUG(cout << "Next event is at time " << (t+t_next) << endl);
//			if (t + t_next > nextHostTimePoint) {
//				// do the codivergence / lineage sorting here
//				// CODIVERGENCE or LINEAGE SORTING
//				DEBUG(cout << "This is on the host vertex " << h->getLabel() << endl);
//				double nextHostEventTime(t_final + 1.0);
//				// do codivergence and lineage sorting
//				for (Node* p : occupants) {
//					if (p->doesCodiverge()) {
//						DEBUG(cout << p->getLabel() << " codiverges with this host" << endl);
//						p->codivergeWith(h);
//						for (Node *c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
//							c->onHostVertex() = true;
//							nuNodes.insert(c);
//							c->setTime(c->getHost()->getTime());
//							occupantsAtTime[c->getTime()].insert(c);
//							nextHostEventTime = min(nextHostEventTime, c->getTime());
//							DEBUG(cout << "\t" << c->getLabel() <<":" << c->getHost()->getLabel() << "; t=" << c->getTime() << endl);
//						}
//					} else {
//						// select one nascent host lineage at uniform random:
//						Node* nuHost = h->getFirstChild();
//						if (fran() < 0.5) {
//							nuHost = nuHost->getSibling();
//						}
//						DEBUG(cout << p->getLabel() << " misses the boat and goes down new host lineage " << nuHost->getLabel() << endl);
//						// put p on it & note it for removal from this host node (h):
//						p->setHost(nuHost);
//						p->onHostVertex() = true;
//						nuHost->addParasite(p);
//						occupantsAtTime[nuHost->getTime()].insert(p);
//						nuNodes.insert(p);
//						nextHostEventTime = min(nextHostEventTime, nuHost->getTime());
//					}
//				}
//				for (Node* nu : nuNodes) {
//					occupantsAtTime[nu->getHost()->getTime()].insert(nu);
//				}
//				t = nextHostEventTime;
//				DEBUG(cout << "Setting to the next time point, being a host node" << endl);
//				t = nextHostTimePoint;
//				if (currentHosts == hNodeTimes.end()) {
//					break;
//				}
//				++currentHosts;
//				continue;	// no more *host* events before the next host node
//			} else {
//				// choose an extant parasite lineage and do the events on that
//				DEBUG(cout << "Selecting the eventful lineage: ");
//				// now select which node is going to do something:
//				Node *q = getRandomElement<Node*>(occupants);
//				DEBUG(cout << q->getLabel() << endl);
//				// which event? :
//				double pB = q->getBirthRate();
//				double pS = (q->hasParent()) ? q->getHostSwitchRate() : 0.0;
//				double pX = q->getDeathRate();
//				double pTotal = pB + pS + pX;
//				double ran = dran(pTotal);
//				t += t_next;	// advance the current time
//	//			pNodeAtTime[t].erase(q);	// this p-node won't be available for any more events at this time // XXX probably unnecessary to remove this
//				if (ran < pB) {	// DUPLICATION
//					/*
//					 *  	bifurcate p
//					 *  	take this p off h and replace it with its two children
//					 *  	add these two child ps to occupantsAtTime with the new time
//					 */
//					DEBUG(cout << "DUPLICATION event:" << endl);
//					h = q->getHost();
//					Tree* T = q->getTree();
//					DEBUG(
//							cout << "Tree before:" << endl << *T << endl;
//					);
//					q->bifurcate();
//					DEBUG(
//							cout << "Tree AFTER:" << endl << *T << endl;
//					);
//					q->setEvent(duplication);
//					q->setTime(t);
//					q->onHostVertex() = false;
//					h->getParasites().erase(q);
//					DEBUG(cout << "Checking children of newly duplicated " << q->getLabel() << endl);
//					for (Node * c = q->getFirstChild(); c != nullptr; c = c->getSibling()) {
//						c->setHost(h);
//						h->getParasites().insert(c);
//						c->setTime(t);
//						c->onHostVertex() = false;
//						occupantsAtTime[t].insert(c);	// now the edges descendant from q are available for events at the next timestep
//						DEBUG(cout << "New parasite node " << c->getLabel() << " on host " << h->getLabel() << " at time " << t << endl);
//						DEBUG(cout << "\t which has parent " << q->getLabel() << endl);
//					}
//					DEBUG(cout << "First child of parasite node " << q->getLabel() << " is " << q->getFirstChild()->getLabel() << endl);
//				} else if (ran < pB + pS) {	// HOST SWITCH
//					DEBUG(cout << "HOST SWITCH event:" << endl);
//					h = q->getHost();
//					availableHosts.erase(h);
//					if (availableHosts.size() == 0) {
//						throw new app_exception("Cannot do a host switch, as there are no available host lineages!");
//					}
//					q->bifurcate();
//					q->setEvent(duplication);
//					q->setTime(t);
//					q->onHostVertex() = false;
//					h->getParasites().erase(q);
//					Node *c = q->getFirstChild();
//					c->setHost(h);
//					h->getParasites().insert(c);
//					c->setTime(t);
//					c->onHostVertex() = false;
//					occupantsAtTime[t].insert(c);
//					c = c->getSibling();
//					Node *nuhost = getRandomElement<Node*>(availableHosts);
//					c->setHost(nuhost);
//					nuhost->getParasites().insert(c);
//					c->setTime(t);
//					c->onHostVertex() = false;
//					occupantsAtTime[t].insert(c);
//				} else {	// DEATH
//					DEBUG(cout << "EXTINCTION (LOSS) event:" << endl);
//					q->setTime(t);
//					q->onHostVertex() = false;
//					q->setEvent(death);
//					h->getParasites().erase(q);
//				}
//	//			DEBUG(cout << *this);
//	//			DEBUG(q->getTree()->gatherVertices());
//	//			DEBUG(q->getTree()->gatherLeaves());
//				DEBUG(cout << "|occupantsAtTime| = " << occupantsAtTime.size() << endl);
//			}
//			if (t >= t_max) {
//				break;
//			}
//		}
//		// Get the time t_next of the next event for the parasites/genes
//		// While t_next is prior the next host/species event, do whatever it is:
//		//		duplication, extinction, host switch.
//		// If t_next is later than the next host/species event,
//		// 	handle the codivergence and lineage sorting, or extinctions, for that node.
//		DEBUG(cout << *this);
////		if (t >= t_final) {
////			DEBUG(cout << "This is at or beyond the age of the host tree so coevolution is stopping." << endl);
////			break;
////		}
//	}
//}

void Cophylogeny::coevolve()
{
	// new try
	bool _debugging(true);
	bool _hostDictatesRate(true);
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
	double hitchDuplicationProb(0.5);	// XXX magic number for parametrising later
	set<Node*> availableHosts;
	unsigned int numCodivergences(0);
	unsigned int numDuplications(0);
	unsigned int numExtinctions(0);
	unsigned int numHostSwitches(0);
	unsigned int numJointDuplicationEvents(0);
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
				eventRate += p->getBirthRate();
				eventRate += p->getDeathRate();
				if (p->getHost()->hasParent()) {
					eventRate += p->getHostSwitchRate();
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
					for (Node* q : active) {
						if (p == q || (q->getHost() != h)) {
							continue;
						}
						if (dran() < hitchDuplicationProb) {
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
						}
					}
					numJointDuplicationEvents += (_jointEvent) ? 1 : 0 ;
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
