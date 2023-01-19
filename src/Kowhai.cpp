/*
 * Kowhai.cpp
 *
 *  Created on: 21 Dec 2022
 *      Author: mac
 */

#include <cstring>
#include <fstream>
#include <string>

#include "Tree.h"
#include "Cophylogeny.h"

using namespace std;
using namespace kowhai;

bool _debugging(false);
bool _for_segdup(false);
string hline("=================================================================\n");
string kowhaiHelp("Kowhai Help");

void testTreeConstructionAndOutput() {
	Tree T;
	T.growYule(10);
	T.displayAs() = phylogram;
	cout << T;
	T.displayAs() = newick;
	cout << T << endl;
	T.displayBranchLengths() = true;
	cout << T << endl;
	T.displayInternalLabels() = true;
	cout << T << endl;
	Tree P;
	P.growYule(5);
	cout << P;
	Tree H("(A,(B,C))");
	cout << H;
}

void testCoevolveBirthModel() {
	Node* h = new Node();
	Tree T(h);
	T.growYule(10);
	cout << T << endl;
	Cophylogeny C;
	C.setHostTree(&T);
	Node *p = C.createParasiteRoot(T.getRoot(), true);
	p->getTree()->setCodivergenceProbability(0.5);
	Node *q = C.createParasiteRoot(T.getRoot(), true);
	q->getTree()->setCodivergenceProbability(0.8);
	C.coevolve();
	cout << C;
}

void testCoevolve() {
	cout << hline << "TESTING cleverCoevolve" << endl << hline;
	Node* h = new Node();
	Tree H(h);
	H.growYule(10);
	H.scaleTo(1.0);
	Cophylogeny C;
	C.setHostTree(&H);
	Node *p = C.createParasiteRoot(H.getRoot(), true);
	Tree *P = p->getTree();
	P->setCodivergenceProbability(0.6);
	P->setBirthRate(0.6);
	P->setDeathRate(0.1);
	P->setHostSwitchRate(0.5);
	P->calculateHeights();
	P->setShowInfo(true);
	P->setLabel("P");
//	Node* q = C.createParasiteRoot(H.getRoot(), true);
//	Tree *Q = q->getTree();
//	Q->setCodivergenceProbability(1.0);
//	Q->setBirthRate(1.0);
//	Q->setDeathRate(0.0);
//	Q->setHostSwitchRate(0.0);
//	Q->calculateHeights();
//	Q->setShowInfo(true);
//	Q->setLabel("Q");
	C.coevolve();
	C.storeAssociationInfo();
	cout << C;
	cout << "Output for segdup:" << endl;
	C.outputForSegdup(cout);
	cout << hline << "TESTING COMPLETE" << endl << hline;
}

void testHostSwitching() {
	cout << hline << "TESTING Host Switching" << endl << hline;
	Node* h = new Node();
	Tree H(h);
	H.growYule(10);
	Cophylogeny C;
	C.setHostTree(&H);
	Node *p = C.createParasiteRoot(H.getRoot(), true);
	Tree *P = p->getTree();
	P->setCodivergenceProbability(0.8);
	P->setBirthRate(20.0);
	P->setDeathRate(0.0);
	P->setHostSwitchRate(0.0);
	P->calculateHeights();
	P->setShowInfo(true);
	P->setLabel("P");
	C.coevolve();
	C.storeAssociationInfo();
	cout << C;
	cout << "Output for segdup:" << endl;
	C.outputForSegdup(cout);
	cout << hline << "TESTING COMPLETE" << endl << hline;
}

int main(int argn, char** argc) {
	if (argn < 2) {
		testCoevolve();
//		testCleverCoevolve();
//		cout << kowhaiHelp << endl;
		return 0;
	}
	bool _sim(false);
	int i(1);
	int numHosts(10);
	int numParasiteTrees(2);
	int numSamples(10);
	double codivProb(0.75);
	if (!strcmp(argc[i], "-sim")) {
		_sim = true;
		++i;
		while (i < argn) {
			if (!strcmp(argc[i], "-nH")) {
				++i;
				numHosts = atoi(argc[i]);
				cout << "Setting number of host species to simulate as " << numHosts << endl;
			} else if (!strcmp(argc[i], "-nP")) {
				++i;
				numParasiteTrees = atoi(argc[i]);
				cout << "Setting number of parasite trees to simulate as " << numParasiteTrees << endl;
			} else if (!strcmp(argc[i], "-s")) {
				++i;
				numSamples = atoi(argc[i]);
				cout << "Setting number of samples as " << numSamples << endl;
			} else if (!strcmp(argc[i], "-pC")) {
				++i;
				codivProb = atof(argc[i]);
				cout << "Setting codivergence probability as " << codivProb << endl;
			} else if (!strcmp(argc[i], "--for-segdup")) {
				_for_segdup = true;
			}
			++i;
		}
	} else {
		testCoevolveBirthModel();
	}
	ofstream fseg;
	if (_for_segdup) {
		fseg.open("for-segdup-from-kowhai.txt", std::ofstream::out);
	}
	if (_sim) {
		for (int s(0); s < numSamples; ++s) {
			Node::resetNodeCounter();
			Node* h = new Node();
			Tree H(h);
			H.growYule(numHosts);
			H.scaleTo(1.0);
			Cophylogeny C;
			C.setHostTree(&H);
			for (int a(0); a < numParasiteTrees; ++a) {
				Node* p = C.createParasiteRoot(H.getRoot(), true);
				p->getTree()->setCodivergenceProbability(codivProb);
			}
			C.coevolve();
			if (_for_segdup) {
				C.outputForSegdup(cout);
			} else {
				cout << C;
				C.outputForSegdup(fseg);
			}
		}

	}
	if (_for_segdup) {
		fseg.close();
	}
	return 0;
}

