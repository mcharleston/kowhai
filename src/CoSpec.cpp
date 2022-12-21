/*
 * CoSpec.cpp
 *
 *  Created on: 21 Dec 2022
 *      Author: mac
 */

#include <string>

#include "Tree.h"

using namespace std;
using namespace cospec;

string cospecHelp("CoSpec Help");

int main(int argn, char** argc) {
	Tree T;
	T.growYule(5);
	cout << T;
	return 0;
}

