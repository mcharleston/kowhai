/*
 * Kowhai.cpp
 *
 *  Created on: 21 Dec 2022
 *      Author: mac
 */

#include <cstring>
#include <fstream>
#include <iomanip>
#include <string>

#include "../utility/debugging.h"
#include "../utility/parser.h"

#include "Tree.h"
#include "Cophylogeny.h"

using namespace std;
using namespace kowhai;
using namespace parsing;

bool _debugging(false);
bool _for_multrec(false);
bool _for_segdup(false);
bool _verbose(false);

ofstream summaryfile;


string hline("=================================================================\n");
//string kowhaiHelp("Kowhai Help");

const bool _defHostDictatesRate(true);
const double defBirthRate(1.0);
const double defDeathRate(0.0);
const double defHostSwitchRate(0.0);
const double defJointDuplicationProb(0.8);
const double defProbCodivergence(1.0);
const double defProbSampling(1.0);
const int defNumHostLeaves(10);
const int defNumParasiteTrees(1);
const int defNumReplicates(1);

namespace kowhai {
bool _hostDictatesRate(_defHostDictatesRate);
double birthRate(defBirthRate);
double deathRate(defDeathRate);
double hostSwitchRate(defHostSwitchRate);
double codivProb(defProbCodivergence);
double jointDuplicationProb(defJointDuplicationProb);
int numHosts(defNumHostLeaves);
int numParasiteTrees(defNumParasiteTrees);
int numSamples(defNumReplicates);
}

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
//	delete p->getTree();
//	delete q->getTree();
}

void testCoevolve() {
	cout << hline << "TESTING Coevolve" << endl << hline;
	Node* h = new Node();
	Tree H(h);
	H.growYule(10);
	H.scaleTo(1.0);
	Cophylogeny C;
	C.setHostTree(&H);
	Node *p = C.createParasiteRoot(H.getRoot(), 0.2);

	// This failed,but I'd like to be able to do something like grow a tree P
	// by duplicating a while then putting all its leaves on the root of H:
//	h = p->getHost();
//	p->bifurcate();
//	double t(-0.1);
//	p->setTime(t);
//	p->setEvent(duplication);
//	for (Node* c = p->getFirstChild(); c != nullptr; c = c->getSibling()) {
//		c->setHost(h);
//		c->setTime(t);
//	}

	Tree *P = p->getTree();
	P->setCodivergenceProbability(0.7);
	P->setBirthRate(0.6);
	P->setDeathRate(0.0);
	P->setHostSwitchRate(0.0);
	P->calculateHeights();
	P->setShowInfo(true);
	P->setLabel("P");

	Node* q = C.createParasiteRoot(H.getRoot(), 0.2);
	Tree *Q = q->getTree();
	Q->setCodivergenceProbability(1.0);
	Q->setBirthRate(1.0);
	Q->setDeathRate(0.0);
	Q->setHostSwitchRate(0.0);
	Q->calculateHeights();
	Q->setShowInfo(true);
	Q->setLabel("Q");

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

string kowhaiHelp("Kowhai Help:\n"
		"\t>./kowhai [options]\n"
		"\t-h or --help\n\t\tto print this help message\n"
		"\t--sim [options] ;\n"
		"\t\tCurrently Kowhai only creates simple simulated cophylogenies, but later it will allow analysis!\n"
		"\t\t-nH <int>\n\t\t\tto set the number of LEAVES/tips in a host/species tree, "
		"generated under a Yule model (default value " + to_string(defNumHostLeaves) + ").\n"
		"\t\t-nP <int>\n\t\t\tto set the number of parasite/gene TREES in each replicate (default value " + to_string(defNumParasiteTrees) + ");\n"
		"\t\t-nR <int>\n\t\t\tto set the number of simulation replicates to do (default value 1);\n"
		"\t\t-pC <float>\n\t\t\tto set the probability of codivergence at each host node (default value " + to_string(defProbCodivergence) + ");\n"
		"\t\t-pJ <float>\n\t\t\tto set the probability of joint duplication (default value " + to_string(defJointDuplicationProb) + ");\n"
		"\t\t-rB <float>\n\t\t\tto set the birth / duplication rate in the dependent phylogenies (default value " + to_string(defBirthRate) + ");\n"
		"\t\t-rHS <float>\n\t\t\tto set the host switch rate in the dependent phylogenies (default value " + to_string(defHostSwitchRate) + ");\n"
		"\t\t-rX <float>\n\t\t\tto set the death rate in the dependent phylogenies (default value " + to_string(defDeathRate) + ");\n"
		"\t\t--for-multrec\n\t\t\tto provide output suitable for segdup input (default: OFF);\n"
		"\t\t--for-segdup\n\t\t\tto provide output suitable for segdup input (default: OFF);\n"
		"\t\t\tNote that the output file \"for-segdup-from-kowhai.txt\" is always produced anyway.\n"
		"\t\t--host-sets-rate\n\t\t\tto set the rates of the dependent phylogenies as determined by the HOST lineage (default: OFF)\n"
		"\t\t--verbose\n\t\t\tset verbosity to output miore stuff.\n"
	);

int main(int argn, char** argc) {
	bool _sim(false);
	int i(1);
	if (argn < 2) {
		cout << kowhaiHelp;
	} else if (!strcmp(argc[i], "-h") || !strcmp(argc[i], "--help")) {
		cout << kowhaiHelp;
//	} else {
//		testCoevolveBirthModel();
	} else if (!strcmp(argc[i], "--sim") || (argn < 2)) {
		_sim = true;
		++i;
		while (i < argn) {
			if (!strcmp(argc[i], "-nH")) {	//matchesIgnoreCase(argc[i], {"-nH"})) {
				++i;
				numHosts = atoi(argc[i]);
				DEBUG(cout << "Setting number of host SPECIES to simulate as " << numHosts << endl);
			} else if (!strcmp(argc[i], "-nP")) {
				++i;
				numParasiteTrees = atoi(argc[i]);
				DEBUG(cout << "Setting number of parasite TREES to simulate as " << numParasiteTrees << endl);
			} else if (!strcmp(argc[i], "-nR")) {
				++i;
				numSamples = atoi(argc[i]);
				DEBUG(cout << "Setting number of samples as " << numSamples << endl);
			} else if (!strcmp(argc[i], "-pC")) {
				++i;
				codivProb = atof(argc[i]);
				DEBUG(cout << "Setting codivergence probability as " << codivProb << endl);
			} else if (!strcmp(argc[i], "-pJ")) {
				++i;
				jointDuplicationProb = atof(argc[i]);
				DEBUG(cout << "Setting joint duplication probability switch rate as " << jointDuplicationProb << endl);
			} else if (!strcmp(argc[i], "-rB")) {
				++i;
				birthRate = atof(argc[i]);
				DEBUG(cout << "Setting birth/duplication rate as " << birthRate << endl);
			} else if (!strcmp(argc[i], "-rHS")) {
				++i;
				hostSwitchRate = atof(argc[i]);
				DEBUG(cout << "Setting host switch rate as " << hostSwitchRate << endl);
			} else if (!strcmp(argc[i], "-rX")) {
				++i;
				deathRate = atof(argc[i]);
				DEBUG(cout << "Setting death rate as " << deathRate << endl);
			} else if (!strcmp(argc[i], "-v")) {
				++i;
				_verbose = true;
			} else if (!strcmp(argc[i], "--for-multrec")) {
				_for_multrec = true;
				DEBUG(cout << "Output simulated cophylogeny for multrec" << endl);
			} else if (!strcmp(argc[i], "--for-segdup")) {
				_for_segdup = true;
				DEBUG(cout << "Output simulated cophylogeny for segdup" << endl);
			} else if (!strcmp(argc[i], "--host-sets-rate")) {
				_hostDictatesRate = true;
				DEBUG(cout << "Duplication rate is set by the host tree, not dependent lineages." << endl);
			} else if (!strcmp(argc[i], "--verbose")) {
				_verbose = true;
				DEBUG(cout << "Setting verbose output." << endl);
			} else if (!strcmp(argc[i], ";")) {
				break;
			} else {
				cout << "Sorry, cannot understand this argument: \"" << argc[i] << "\": please check your input?" << endl;
				cout << kowhaiHelp;
				return (0);
			}
			++i;
		}
	} // add more main choices here, like analysis, later.
	ofstream fseg;
	if (_for_segdup) {
		fseg.open("for-segdup-from-kowhai.txt", std::ofstream::out);
	}
	ofstream fmultrec;
	if (_for_multrec) {
		fmultrec.open("for-multrec-from-kowhai.txt", std::ofstream::out);
	}
	summaryfile.open("summary.csv", std::ios_base::app);
	summaryfile << "nCospec,nIndivididualDups,nAllDupEvents,nJointDups,maxDupSize,meanDupSize,nXtinc,nHostSwitch,nLineageSort\n";
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
				Tree* P = p->getTree();
				P->setLabel("G" + to_string(a+1));
				P->setCodivergenceProbability(codivProb);
			}
			C.coevolve();
			if (_verbose) {
				cout << C;
			}
			if (_for_segdup) {
				C.outputForSegdup(cout);
				C.outputForSegdup(fseg);
			}
			if (_for_multrec) {
				C.outputForMultRec(fmultrec);
				C.outputForMultRec(cout);
			}
		}

	}
	if (_for_multrec) {
		fmultrec.close();
	}
	if (_for_segdup) {
		fseg.close();
	}
	return 0;
}

