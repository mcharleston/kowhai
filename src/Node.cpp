/*
 * Node.cpp
 *
 *  Created on: 1 Aug 2022
 *      Author: mac
 */

#include <map>
#include <set>
#include <stdio.h>

#include "Node.h"
#include "Tree.h"
#include "../utility/appexception.h"
#include "../utility/debugging.h"
#include "../utility/myrandom.h"

using namespace std;

extern bool _debugging;

namespace kowhai {

int Node::nodeCounter = 0;

Node::Node() : _extant(false), parent(nullptr), firstChild(nullptr), sibling(nullptr), host(nullptr), _onHostVertex(false),
		T(nullptr), CoP(nullptr), depth(-1), height(-1), timeIndex(-1), time(0.0), branchLength(0.0), _visited(false) {
	++nodeCounter;
	string str = "v" + to_string(nodeCounter);
	setLabel(str);
}

Node::Node(const Node& n) {
	throw new app_exception("This copy constructor should not be called!");
}


Node::~Node() {
	// Destroy the children first!
	if (!isLeaf()) {
		for (Node* c(firstChild); c != nullptr; c = c->sibling) {
			delete c;
		}
	}
}

void Node::addChild(Node* c) {
	if (isLeaf()) {
		firstChild = c;
	} else {
		Node* n = firstChild;
		while (n->sibling != nullptr) {
			n = n->sibling;
		}
		n->sibling = c;
		c->parent = this;
	}
}

void Node::addSibling(Node *s) {
	if (sibling != nullptr) {
		throw new app_exception("Adding a sibling kills previous sibling: can't condone this!");
	}
	sibling = s;
}

void Node::addSubtreeVertices(map<string,Node*>& V) {
	V[label] = this;
	for (Node* c(firstChild); c != nullptr; c = c->sibling) {
		c->addSubtreeVertices(V);
	}
}

void Node::bifurcate() {
	Node* u = new Node();
	Node* v = new Node();
	bifurcate(u, v);
}

void Node::bifurcate(string a, string b) {
	if (firstChild != nullptr) {
		throw new app_exception("Cannot bifurcate this node as it already has children!");
	}
	Node* u = new Node(a);
	Node* v = new Node(b);
	bifurcate(u, v);
}

void Node::bifurcate(Node* u, Node* v) {
	bool _debugging(true);
	if (firstChild != nullptr) {
		throw new app_exception("Node::bifurcate(Node *, Node*): attempting to bifurcate with existing child nodes!");
	}
	firstChild = u;
	u->sibling = v;
	u->parent = this;
	u->T = T;
	v->parent = this;
	v->T = T;
	T->getVertices()[u->getLabel()] = u;
	T->getVertices()[v->getLabel()] = v;
	height = 1;
	Node *par = parent;
	int h = height;
	while (par != nullptr) {
		par->height = max(par->height, h+1);
		h = par->height;
		par = par->parent;
	}
}

void Node::calcDepth() {
	/**
	 * The depth of a node is the length of the (shortest) path from it to the root.
	 * In a tree this path is of course unique.
	 */
	if (parent == nullptr) {
		depth = 0;
		return;
	}
	Node* par = parent;
	depth = par->getDepth() + 1;
}

void Node::calcHeight() {
	/**
	 * The height of a node is the maximum length of the unique path from the node to each of its descendant leaves.
	 */
	height = 0;
	if (isLeaf()) {
		height = 0;
	} else {
		int h(0);
		for (Node* c = firstChild; c != nullptr; c = c->sibling) {
			h = std::max<int>(h, c->getHeight()+1);
		}
		height = h;
	}
}

void Node::codivergeWith(Node* h) {
	bool _debugging(true);
	if (h->isLeaf()) {
		throw new app_exception("Node::codivergeWith(h): host node h=" + h->getLabel() + " is a leaf" );
	}
	bifurcate();
	DEBUG(cout << "codivergeWith: new children are " << firstChild->getLabel() << " and " << firstChild->getSibling()->getLabel() << endl);
	Node* x = h->firstChild;
	for (Node *c = firstChild; c != nullptr; c= c->sibling) {
		c->setHost(x);
		c->onHostVertex() = true;
		x->addParasite(c);
		DEBUG(cout << "Adding parasite " << c->getLabel() << " to host " << x->getLabel() << endl);
		DEBUG(cout << c->getLabel() << ":" << x->getLabel() << endl);
		x = x->sibling;
	}
}

bool Node::doesCodiverge() {
	bool _debugging(true);
	double prob = T->getCodivergenceProbability();
	double r(fran());
	DEBUG(cout << "Codivergence probability: " << prob << "; r = " << r << endl);
	return (r < prob);
}

void Node::diverge() {
	/**
	 * As opposed to Node::bifurcate, which is simply a tree construction process, Node::diverge allows dependent
	 * (parasite/gene) nodes to also bifurcate if they are extant at the same time.
	 * Because it's a terrible idea to attempt to compare floats or doubles for equality in this context,"at the
	 * same time" means "at the same integer index of all events in the whole cophylogeny thing".
	 *
	 * Procedure:
	 * 		bifurcate
	 * 		look at underlying Cophylogeny and find all other nodes at this time index
	 * 		for all such nodes, codiverge or miss the boat.
	 * 			(at this point there's no plan to include widespread parasites because I don't know how to model it)
	 * 			These new nodes get the same time index as their parent; this will be updated when *they* get an event.
	 */
	bifurcate();
	for (Node *p : parasites) {
		if (p->getTimeIndex() == timeIndex) {
			if (p->doesCodiverge()) {
				p->bifurcate();
			} else {

			}
		}
	}
}

double Node::getBirthRate() const {
	return T->getBirthRate();
}

double Node::getDeathRate() const {
	return T->getDeathRate();
}

int Node::getDepth() {
	if (depth < 0) {
		calcDepth();
	}
	return depth;
}

int Node::getHeight() {
//	cout << "Node::getHeight()..." << endl;
	if (height < 0) {
		calcHeight();
	}
	return height;
}

double Node::getHostSwitchRate() const {
	return T->getHostSwitchRate();
}

double Node::getMaxDescendantTime() const {
	if (isLeaf()) {
		return time;
	}
	double t = time;
	for (Node* c = firstChild; c != nullptr; c = c->sibling) {
		t = std::max(t, c->getMaxDescendantTime());
	}
	return t;
}

void Node::initialiseOccupants(std::map<double, std::set<Node*>>& occ) {
	occ[time];
	for (Node* c = firstChild; c != nullptr; c = c->sibling) {
		c->initialiseOccupants(occ);
	}
}

bool Node::isAncestralTo(Node* other) {
	return (T->isAncestralTo(this, other));
}

Node* Node::next() {
	Node* nextNode = nullptr;
	if (firstChild != nullptr) {
		return firstChild;
	}
	if (sibling != nullptr) {
		return sibling;
	}
	nextNode = this;
	do {
		if (nextNode->parent == nullptr) {
			return nullptr;
		}
		nextNode = nextNode->parent;
	} while (nextNode->sibling == nullptr);
	return nextNode->sibling;
}

void Node::putChildren(set<Node*>& children) {
	for (Node* c = firstChild; c != nullptr; c = c->sibling) {
		children.insert(c);
	}
}

void Node::scaleBy(double d) {
	time *= d;
	for (Node* c = firstChild; c != nullptr; c = c->sibling) {
		c->scaleBy(d);
	}
}

void Node::setFirstChild(Node* c) {
	firstChild = c;
	c->parent = this;
}

void Node::storeNodeTimes(std::map<double, std::set<Node*>>& M) {
	std::set<Node*>& S = M[time];
	S.insert(this);
	for (Node* c = firstChild; c != nullptr; c = c->getSibling()) {
		c->storeNodeTimes(M);
	}
}

void Node::writeNewick(ostream& os) {
	if (isLeaf()) {
		os << label;
		if (T->displayBranchLengths()) {
			os << ':' << to_string(branchLength);
		}
	} else {
		os << '(';
		Node* c = firstChild;
		c->writeNewick(os);
		c = c->sibling;
		while (c != nullptr) {
			os << ',';
			c->writeNewick(os);
			c = c->sibling;
		}
		os << ')';
		if (T->displayInternalLabels()) {
			os << label;
		}
		if (T->displayBranchLengths()) {
			os << ':' << to_string(branchLength);
		}
	}
}

std::ostream& operator<<(std::ostream& os, Node& n) {
	if (n.isLeaf()) {
		os << n.getLabel();
	} else {
		os << '\"' << n.getLabel() << "\":";
		os << '(';
		for (Node* c = n.getFirstChild(); c != nullptr; c = c->getSibling()) {
			os << *c;
			if (c->getSibling() != nullptr) {
				os << ',';
			}
		}
		os << ')';
	}
	return os;
}

} /* namespace kowhai */
