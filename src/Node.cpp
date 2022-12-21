/*
 * Node.cpp
 *
 *  Created on: 1 Aug 2022
 *      Author: mac
 */

#include <set>
#include <stdio.h>

#include "Node.h"
#include "Tree.h"
#include "../utility/appexception.h"

using namespace std;

namespace cospec {

Node::~Node() {
	// Destroy the children first!
	if (!isLeaf()) {
		for (Node* c(firstChild); c != nullptr; c = c->sibling) {
//			cout << "Deleting child node " << c->getLabel() << endl;
			delete c;
		}
	}
//	cout << "Deleting this node " << label << endl;
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

void Node::bifurcate(string a, string b) {
	if (firstChild != nullptr) {
		throw new app_exception("Cannot bifurcate this node as it already has children!");
	}
	Node* u = new Node(a);
	Node* v = new Node(b);
	firstChild = u;
	u->sibling = v;
	u->parent = this;
	u->T = T;
	v->parent = this;
	v->T = T;
}

void Node::bifurcate(Node* u, Node* v) {
	firstChild = u;
	u->sibling = v;
	u->parent = this;
	u->T = T;
	v->parent = this;
	v->T = T;
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

void Node::setFirstChild(Node* c) {
	firstChild = c;
	c->parent = this;
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

} /* namespace segdup */
