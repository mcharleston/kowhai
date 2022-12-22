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

namespace kowhai {

enum TreeDisplayFormat {
	phylogram,
	newick
};

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
	double birthRate;
	std::map<Node*, std::string>* info;
	TreeDisplayFormat displayFormat;
	bool _displayBranchLengths;
	bool _displayInternalLabels;
public:
	Tree() : label("untitled"), root(nullptr), labelSpace(0), numVertices(-1), _showInfo(false), prefix('v'), info(nullptr) {
		birthRate = 1.0;
		displayFormat = phylogram;
		_displayBranchLengths = false;
		_displayInternalLabels = false;
	}
	Tree(Node* r) : label("untitled"), root(r), labelSpace(0), numVertices(-1), _showInfo(false), prefix('v'), info(nullptr) {
		birthRate = 1.0;
		displayFormat = phylogram;
		_displayBranchLengths = false;
		_displayInternalLabels = false;
	}
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

	TreeDisplayFormat& displayAs() { return displayFormat; }
	const TreeDisplayFormat& displayAs() const { return displayFormat; }
	inline bool& displayBranchLengths() { return _displayBranchLengths; }
	inline const bool& displayBranchLengths() const { return _displayBranchLengths; }
	inline bool& displayInternalLabels() { return _displayInternalLabels; }
	const inline bool& displayInternalLabels() const { return _displayInternalLabels; }

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
	void growYule(int numLeaves);

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

	void writeNewick(std::ostream& os);

};

std::ostream& operator<<(std::ostream& os, Tree& T);

} /* namespace kowhai */

#endif /* TREE_H_ */
