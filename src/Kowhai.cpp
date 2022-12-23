/*
 * Kowhai.cpp
 *
 *  Created on: 21 Dec 2022
 *      Author: mac
 */

#include <string>

#include "Tree.h"
#include "Cophylogeny.h"

using namespace std;
using namespace kowhai;

bool _debugging(false);

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

int main(int argn, char** argc) {
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
	return 0;
}

