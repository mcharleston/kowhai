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

int main(int argn, char** argc) {
	if (argn < 2) {
		cout << kowhaiHelp << endl;
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
			Node* h = new Node();
			Tree H(h);
			H.growYule(numHosts);
			Cophylogeny C;
			C.setHostTree(&H);
			for (int a(0); a < numParasiteTrees; ++a) {
				Node* p = C.createParasiteRoot(H.getRoot(), true);
				p->getTree()->setCodivergenceProbability(codivProb);
			}
			C.coevolve();
			if (_for_segdup) {
				fseg << "-S \"";
				H.writeNewick(fseg);
				fseg << "\"";
				for (auto P : C.getParasiteTrees()) {
					fseg << " -G \"";
					P->writeNewick(fseg);
					for (auto pr : P->getLeaves()) {
						Node* p = pr.second;
						fseg << pr.first << ':' << p->getHost()->getLabel() << " ";
					}
					fseg << "\"";
				}
				fseg << endl;
			}
			cout << C;
		}

	}
	if (_for_segdup) {
		fseg.close();
	}
	return 0;
}

