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
	Cophylogeny C;
	C.setHostTree(new Tree("(A,(B,C))"));
	return 0;
}

