/*
 * NEXUSParser.cpp
 *
 *  Created on: 20 Jul 2016
 *      Author: mac
 */


#include <algorithm>
#include <exception>
#include <stdio.h>
#include <sstream>
#include "../utility/debugging.h"
#include "../utility/appexception.h"
#include "../utility/parser.h"
#include "KowhaiParser.h"
#include "Node.h"
#include "Tree.h"
#include "project.h"

using namespace std;
using namespace kowhai;

extern bool _debugging;

namespace parsing {

KowhaiParser::KowhaiParser(const std::string& fileName, Project* pr)
		: Parser(fileName), proj(pr) {
}

void KowhaiParser::parse() {
	/*
	 * Grammar of NEXUS file:
	 * <nexus> ::= "#NEXUS" { <nexusblock> }
	 * <nexusblock> ::= "begin" <blockname> { <setting> } ( "end" | "endblock" ) ";"
	 * <blockname> ::= <string>
	 * <setting> ::= <id> [ "=" ] <val> [","]
	 */
	TL.reset();
	eat('#');
	eat("nexus");
	while (TL.hasNext()) {
		parseNEXUSBlock();
	}
}

void KowhaiParser::parseBranchLength(Node* v) {
	if (matches(':')) {
		advance();
		double d = getDouble();
		v->setBranchLength(d);
		advance();
	}
}

void KowhaiParser::parseKowhaiBlock() {
	/**
	 * EBNF format (case insensitive throughout):
	 * KowhaiBlock := "Kowhai" { <hosttree> | <associatetree> | <eventcosts> }
	 * ..?
	 */
	bool _debugging = false;
	eat("Kowhai");
	ignore(';');
	while (hasNext() && !matches({ "end", "endblock" })) {
//		advance();
		if (matches("hosttree")) {
			Tree* T = new Tree();
			parseNewickSubtree(T->getRoot(), 'v');	// XXX should move this to Phylogeny class.
			T->compressTraverseWrite(cout);
			ignore(';');
		} else if (matches("speciestree")) {
			advance();
		}

		if (hasNext()) {
			advance();
		}
	}
	DEBUG(cout << "Completed Kowhai block" << endl);
	advance();
	ignore(';');
}

void KowhaiParser::parseNEXUSBlock() {
	skipComments();
	eat("begin");
	if (matches("Kowhai")) {
		parseKowhaiBlock();
	} else {
		skipBlock();
	}
}

void KowhaiParser::parseNewickTree(Tree* T) {
	/**
	 * Newick format grammar, adapted from http://evolution.genetics.washington.edu/phylip/newick_doc.html
	 *
	 * <tree>            ::= <descendant_list> [ <rootlabel> ] [ ':' <branchlength> ] ';'
	 * <descendant_list> ::= '(' <subtree> { ',' <subtree> } ')'
	 * <subtree>         ::= (
	 *                          <descendant_list> [ <internal_label> ] [ ':' <branchlength> ]
	 *                       |
	 *                          <leaf_label> [ ':' <branchlength> ]
	 *                       )
	 * <rootlabel>       ::= <label>
	 * <internal_label>  ::= <label>
	 * <leaf_label>      ::= <label>
	 * <label>           ::= <string>
	 * <branchlength>    ::= <number>
	 *
	 * Further, the Newick format can be preceded by an optional tree name:
	 * <treename>        ::= "tree" [ '=' ] <string>
	 */
	if (matches("tree")) {
		advance();
		string treeName = getString();
		advance();
		ignore('=');
	}
	Node* root = new Node();
	parseNewickSubtree(root, T->getPrefixChar());
	T->setRoot(root);
	T->gatherVertices();
	T->calculateHeights();
}

void KowhaiParser::parseNewickSubtree(Node* v, char prefix) {
	/**
	 * See NEXUSParser::parseNewickFormatTree for Newick format
	 */
try {
	bool _debugging(false);
//	static int numInternals(0);
	char internalLabel[16];
	DEBUG(cout << current() << endl);
	DEBUG(cout << "parsing subtree" << endl);
	if (matches('(')) {
		// internal node
		advance();
		Node* child = new Node();
		parseNewickSubtree(child, prefix);
		if (!child->isLeaf()) {
			if (isString()) {
				child->setLabel(getString());
				advance();
			}
		}
		parseBranchLength(child);
		v->setFirstChild(child);
//		if (v->getLabel() == "*") {
//			++numInternals;
//			sprintf(internalLabel, "%c%d", prefix, numInternals);
//			v->setLabel(internalLabel);
//		}
		DEBUG(cout << " set " << v->getLabel() << " as parent of " << child->getLabel() << endl);
		while (matches(',')) {
			// siblings!
			advance();
			Node* sib = new Node();
			parseNewickSubtree(sib, prefix);
			if (!child->isLeaf()) {
				if (isString()) {
					child->setLabel(getString());
					advance();
				} else {
//					++numInternals;
//					sprintf(internalLabel, "%c%d", prefix, numInternals);
//					child->setLabel(internalLabel);
				}
			}
			parseBranchLength(sib);
			sib->setParent(v);
			DEBUG(cout << " set " << v->getLabel() << " as parent of " << sib->getLabel() << endl);
			child->addSibling(sib);
			DEBUG(cout << " set " << sib->getLabel() << " as sibling of " << child->getLabel() << endl);
			child = sib;
		}
		eat(')');
		ignore(';');
	} else if (isString()) {
//		Node* leaf = new Node();
		v->setLabel(getString());
		advance();
//		parseBranchLength(v);
//		leaf->setParent(v);
//		v->addChild(leaf);
	}
} catch (exception& e) {
	cout << "Something has gone wrong with parsing the tree. Please check your input." << endl;
	cout << e.what() << endl;
}
}

void KowhaiParser::skipBlock() {
	while (!matches({"end", "endblock"})) {
		advance();
	}
	advance();
	eat(';');
}

} /* namespace parsing */
