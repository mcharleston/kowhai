/*
 * project.h
 *
 *  Created on: 25 Jul 2016
 *      Author: mac
 */

#ifndef SRC_PROJECT_H_
#define SRC_PROJECT_H_

#include <set>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "../utility/parser.h"

namespace cospec {

class Project {
private:
	std::string name;	/// name of this project; used as base name for file outputs, e.g. name-landscape.csv, name-graph.dot.
	parsing::Parser* parser;
	int sampleSize;

	// Data:
//	std::vector<Phylogeny> trees;
	size_t ntax;		/// number of taxa (= tips = leaves)
	std::map<std::string, int> leafLabels;


	// For alignments:
	std::set<std::string> taxon;

//	std::ostream os; // main output stream

public:
	Project();
	virtual ~Project();

	Project(std::string filename);

	size_t getNumTaxa() const { return ntax; }
	void initialise();
	void read(const std::string& fileName);
	void run();
	void setName(const std::string& base) { name = base; }
	void setNumTaxa(int n) { ntax = n; }
	void setSampleSize(int n) { sampleSize = n; }
};

} /* namespace flatbush */

#endif /* SRC_PROJECT_H_ */
