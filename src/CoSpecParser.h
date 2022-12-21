/*
 * NEXUSParser.h
 *
 *  Created on: 20 Jul 2016
 *      Author: mac
 */

#ifndef COSPECPARSER_H_
#define COSPECPARSER_H_

#include "Node.h"
#include "../utility/parser.h"
#include "Tree.h"

namespace parsing {

const std::vector<std::string> NEXUSSuffixes = { ".nx", ".nxs", ".nex", ".nexus" };

class CoSpecParser : public Parser {
private:
	cospec::Project* proj;

public:
	virtual ~CoSpecParser() { }
	CoSpecParser(const std::string& fileName, cospec::Project* pr);
	CoSpecParser(TokenList& tl) : Parser(tl), proj(nullptr) {}

	void parse();
	void parseBranchLength(cospec::Node* v);
	void parseCoSpecBlock();
	void parseNEXUSBlock();
	void parseNewickTree(cospec::Tree* T);
	void parseNewickSubtree(cospec::Node* v, char prefix);
	void parseTaxaBlock();

	void skipBlock();
};

} /* namespace parsing */

#endif /* SEGDUPPARSER_H_ */
