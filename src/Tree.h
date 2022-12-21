/*
 * Tree.h
 *
 *  Created on: 1 Aug 2022
 *      Author: mac
 */

#ifndef TREE_H_
#define TREE_H_

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>

#include "Node.h"

namespace cospec {

class Tree {
private:
	std::string label;
	Node* root;
	std::map<std::pair<Node*, Node*>, short> distUp;
	std::map<std::string, Node*> V;	// vertices by name
	int labelSpace;	// how much does each vertex label need?
	int numVertices;
	bool _showInfo;
	char prefix;
	std::map<Node*, std::string>* info; // TODO this needs to be turned into a pointer, referring back to maybe an original in a cophymap or null.
//	std::map<Node*, eventType> event;

public:
	Tree() : label("untitled"), root(nullptr), labelSpace(0), numVertices(-1), _showInfo(false), prefix('v'), info(nullptr) {}
	Tree(Node* r) : label("untitled"), root(r), labelSpace(0), numVertices(-1), _showInfo(false), prefix('v'), info(nullptr) {}
	Tree(char prefix, std::string);
	Tree(std::string str);
	~Tree() { delete root; }

	Node* at(const std::string& str) { return V.at(str); }

	void calcAncestry();
	void calculateHeights();
	inline void calculateHeights(Node* v) { root->calcHeight(); }
	void compressTraverseWrite(std::ostream& os);
	void compressTraverseWrite(std::ostream& os, Node* v);
	void constructFromNewickString(std::string str);
//	void compressTraverseWriteOld(std::ostream& os, Node* v, bool _showAssociations = false);

	void gatherVertices();
	int getDistUp(Node* lower, Node* upper);
	inline unsigned int getHeight() { return root->getHeight(); }
	std::map<Node*, std::string>* getInfo() { return info; }
	std::string getLabel() { return label; }
	Node* getLCAofChildren(Node *p);
	int getMaxLabelWidth(Node* v);
	inline int getNumEdges() const { return V.size() - 1; }
	inline int getNumVertices() const { return V.size(); }
	char getPrefixChar() const { return prefix; }
	Node* getRoot() { return root; }
	inline bool getShowInfo() const { return _showInfo; }
	std::map<std::string, Node*>& getVertices() { if (V.size() == 0) { gatherVertices(); } return V; }

	bool isAncestralTo(Node* x, Node* y);

	Node* LCA(std::set<Node*> V);
	Node* LCA(Node* u, Node* v);
	Node* LCA(const std::string& ustr, const std::string& vstr);

	Node* operator[](const std::string& str) { return V[str]; }
	Tree& operator=(const std::string& str);

	void putInternalVertices(std::set<Node*>& IV);

	void setInfo(std::map<Node*, std::string>* inf) { info = inf; }
	void setLabel(const std::string& str) { label = str; }
	void setNodeLabel(const std::string& str, const std::string& newlabel);
	void setNodeLabel(Node* n, const std::string& newlabel);
	void setRoot(Node* r) { root = r; }
	void setShowInfo(bool b) { _showInfo = b; }

};

std::ostream& operator<<(std::ostream& os, Tree& T);

} /* namespace segdup */

#endif /* TREE_H_ */
