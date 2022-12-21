/*
 * NEXUSParser.h
 *
 *  Created on: 20 Jul 2016
 *      Author: mac
 */

#ifndef KOWHAIPARSER_H_
#define KOWHAIPARSER_H_

#include "Node.h"
#include "../utility/parser.h"
#include "Tree.h"

namespace parsing {

const std::vector<std::string> NEXUSSuffixes = { ".nx", ".nxs", ".nex", ".nexus" };

class KowhaiParser : public Parser {
private:
	Kowhai::Project* proj;

public:
	virtual ~KowhaiParser() { }
	KowhaiParser(const std::string& fileName, Kowhai::Project* pr);
	KowhaiParser(TokenList& tl) : Parser(tl), proj(nullptr) {}

	void parse();
	void parseBranchLength(Kowhai::Node* v);
	void parseKowhaiBlock();
	void parseNEXUSBlock();
	void parseNewickTree(Kowhai::Tree* T);
	void parseNewickSubtree(Kowhai::Node* v, char prefix);
	void parseTaxaBlock();

	void skipBlock();
};

} /* namespace parsing */

#endif /* SEGDUPPARSER_H_ */
