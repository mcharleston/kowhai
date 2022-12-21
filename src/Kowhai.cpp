/*
 * Kowhai.cpp
 *
 *  Created on: 21 Dec 2022
 *      Author: mac
 */

#include <string>

#include "Tree.h"

using namespace std;
using namespace kowhai;

string kowhaiHelp("Kowhai Help");

int main(int argn, char** argc) {
	Tree T;
	T.growYule(6);
	T.displayAs() = phylogram;
	cout << T;
	T.displayAs() = newick;
	cout << T << endl;
	return 0;
}

