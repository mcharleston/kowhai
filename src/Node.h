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

class Cophylogeny;
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
	Node* host;	// e.g. if this is a parasite or gene, the host species;
	bool _onHostVertex;
	std::map<std::string, Node*> parasites;	// e.g. if this is a host, the parasites or genes on it
	Tree* T;
	Cophylogeny* CoP;
	int depth;
	int height;
	int timeIndex;
	double time;	// beginning at time t=0 at the root
	double branchLength;
	mutable bool _visited;
public:
	Node();
	Node(std::string str) : label(str),
				parent(nullptr), firstChild(nullptr), sibling(nullptr), host(nullptr), _onHostVertex(false),
				T(nullptr), CoP(nullptr), depth(-1), height(-1), timeIndex(-1), time(0.0), branchLength(0.0), _visited(false) {}
	virtual ~Node();

	void addParasite(Node* p) { parasites[p->getLabel()] = p; }
	void addChild(Node* c);
	inline void addToBranchLength(double d) { branchLength += d; }
	void addSibling(Node *s);
	void addSubtreeVertices(std::map<std::string, Node*>& V);

	void bifurcate();
	void bifurcate(std::string a, std::string b);
	void bifurcate(Node* u, Node* v);

	void calcDepth();
	void calcHeight();
	inline void clearParasites() { parasites.clear(); }
	void codivergeWith(Node* h);

	void diverge();
	bool doesCodiverge();

	double getBirthRate() const;
	double getBranchLength() const { return branchLength; }
	double getDeathRate() const;
	int getDepth();
	int getHeight();
	inline const std::string& getLabel() const { return label; }
	inline Node* getFirstChild() { return firstChild; }
	inline Node* getHost() { return host; }
	inline const Node* getHost() const { return host; }
	inline std::string& getLabel() { return label; }
	inline Node* getParasite(std::string str) { return parasites[str]; }
	inline std::map<std::string, Node*>& getParasites() { return parasites; }
	inline const std::map<std::string, Node*>& getParasites() const { return parasites; }
	inline Node* getParent() const { return parent; }
	inline Node* getSibling() { return sibling; }
	inline double getTime() const { return time; }
	inline int getTimeIndex() { return timeIndex; }
	const Tree* getTree() const { return T; }
	Tree* getTree() { return T; }

	inline bool hasHost() const { return host!=nullptr;}
	inline bool hasParent() const { return parent!=nullptr; }

	bool isAncestralTo(Node* other);
	void inferEvents();

	inline bool isFirstChild() { if (parent==nullptr) return false; return (parent->firstChild==this); }
	inline bool isLeaf() const { return firstChild==nullptr; }

	Node* next();

	inline bool& onHostVertex() { return _onHostVertex; }
	inline const bool& onHostVertex() const { return _onHostVertex; }
	bool operator<=(Node& o) { return this->isAncestralTo(&o); }

	void putChildren(std::set<Node*>& children);

	inline void setBranchLength(double d) { branchLength = d; }
	void setFirstChild(Node* c);
	inline void setHost(Node* h) { host = h; }
	inline void setLabel(std::string str) { label = str; }
	inline void setParent(Node* p) { parent = p; }
	inline void setSibling(Node* sib) { sibling = sib; sib->parent = parent; }
	inline void setTime(double d) { time = d; }
	inline void setTimeIndex(int idx) { timeIndex = idx; }
	inline void setTree(Tree *tr) { T = tr; }

	void writeNewick(std::ostream& os);
};


std::ostream& operator<<(std::ostream& os, Node& n) ;

} /* namespace kowhai */

#endif /* NODE_H_ */
