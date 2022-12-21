/*
 * Node.h
 *
 *  Created on: 1 Aug 2022
 *      Author: mac
 */

#ifndef NODE_H_
#define NODE_H_

#include <iostream>
#include <map>
#include <set>
#include <string>


namespace kowhai {

class Tree;

const short duplicationEvent(0);
const short codivergenceEvent(1);

enum eventType {
	codivergence,
	duplication,
	loss,
	noevent
};
const char eventSymbol[4]{'<','=','x',' '};

class Node {

	friend class CophyMap;
	friend class NodeMap;
	friend class CophyMultiMap;	// so many friends... probably should cut down.

private:
	std::string label;
	Node* parent;
	Node* firstChild;
	Node* sibling;
	Node* associate;
	Tree* T;
//	eventType event;
	int depth;
	int height;
	double branchLength;
	mutable bool _visited;
public:
	Node() : label("*"),
				parent(nullptr), firstChild(nullptr), sibling(nullptr), associate(nullptr),
				T(nullptr), depth(-1), height(-1), branchLength(0.0), _visited(false) {}
	Node(std::string str) : label(str),
				parent(nullptr), firstChild(nullptr), sibling(nullptr), associate(nullptr),
				T(nullptr), depth(-1), height(-1), branchLength(0.0), _visited(false) {}
	virtual ~Node();

	void addChild(Node* c);
	void addSibling(Node *s);
	void addSubtreeVertices(std::map<std::string, Node*>& V);

	void bifurcate(std::string a, std::string b);
	void bifurcate(Node* u, Node* v);

	void calcDepth();
	void calcHeight();

	std::string describeEvent();

	inline double getBranchLength() const { return branchLength; }
	int getDepth();
//	inline eventType& getEvent() { return event; }
	int getHeight();
	inline const std::string& getLabel() const { return label; }
	inline Node* getFirstChild() { return firstChild; }
	inline std::string& getLabel() { return label; }
	inline Node* getParent() const { return parent; }
	inline Node* getSibling() { return sibling; }
	const Tree* getTree() const { return T; }
	Tree* getTree() { return T; }

	bool hasParent() const { return parent!=nullptr; }

	bool isAncestralTo(Node* other);
	void inferEvents();

	inline bool isFirstChild() { if (parent==nullptr) return false; return (parent->firstChild==this); }
	inline bool isLeaf() const { return firstChild==nullptr; }

	Node* next();

	bool operator<=(Node& o) { return this->isAncestralTo(&o); }

	void putChildren(std::set<Node*>& children);

	inline void setBranchLength(double d) { branchLength = d; }
	void setFirstChild(Node* c);
	inline void setLabel(std::string str) { label = str; }
	inline void setParent(Node* p) { parent = p; }
	inline void setSibling(Node* sib) { sibling = sib; sib->parent = parent; }
	inline void setTree(Tree *tr) { T = tr; }
};


std::ostream& operator<<(std::ostream& os, Node& n) ;

} /* namespace segdup */

#endif /* NODE_H_ */
