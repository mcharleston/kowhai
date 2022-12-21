/*
 * project.cpp
 *
 *  Created on: 25 Jul 2016
 *      Author: mac
 */

#include <algorithm>
#include <chrono>
#include <ctime>
#include <exception>
#include <iomanip>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <valarray>
#include <vector>
#include "project.h"
#include "../utility/appexception.h"
#include "../utility/debugging.h"
#include "../utility/parser.h"
#include "Tree.h"
#include "KowhaiParser.h"

using namespace std;

extern bool _debugging;
extern bool _silent;
extern bool _verbose;

namespace kowhai {

bool matchesIgnoreCase(string str, const std::vector<std::string>& patterns) {
	// Convert initializer_list to vector and call the (string, vector) overload.
	bool _debugging = false;
	DEBUG(cout << "matchesIgnoreCase(" << str << ", patterns)" << endl);
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	vector<string> vec(patterns.begin(), patterns.end());
//	return matchesIgnoreCase(str, vec);
	for (string s : patterns) {
		transform(s.begin(), s.end(), s.begin(), ::tolower);
		DEBUG(cout << "\tcompared string s=" << s << endl);
		if (str.compare(s) == 0) {
			return true;
		}
	}
	return false;
}

Project::Project() : name("unnamed project"), parser(nullptr), sampleSize(1000),
		ntax(0)
{
	initialise();
}

Project::Project(std::string filename) : parser(nullptr), sampleSize(1000)
{
	// load the file with this filename, tokenize it, then parse & create all the defined objects
	read(filename);
}

Project::~Project() {

}

void Project::initialise() {
	char buffer[20];
	struct tm *timeinfo;
	time_t rawtime;
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(buffer, 20, "%Y-%m-%d-%H-%M-%S", timeinfo);
	string timeStr = buffer;
	name = "flatbush-project-" + timeStr;
}

void Project::run() {
	bool _debugging(true);
	DEBUG(cout << "Running now" << endl);
	// Do the Thing!
	if (!_silent) {
		cout << "Done." << endl;
	}
}

void Project::read(const std::string& filename) {
	bool _debugging(true);
	DEBUG(cout << "input filename = " << filename << endl);
	size_t pos = filename.find_last_of(".");
	string suffix = filename.substr(pos);
	if (matchesIgnoreCase(suffix, parsing::NEXUSSuffixes)) {
		// parse as NEXUS
		DEBUG(cout << "parsing " << filename << " as NEXUSformat file" << endl);
		parser = new parsing::KowhaiParser(filename, this);
	} else {
		stringstream ss;
		ss << "Cannot recognise the file name \'" << filename << "\': giving up.";
		throw new app_exception(ss.str());
	}

	DEBUG(cout << "setting project" << endl);
	parser->setProject(this);
	DEBUG(cout << "project set; now beginning parsing" << endl);
	parser->parse();
	DEBUG(cout << "Parsing complete." << endl);
}

} /* namespace flatbush */
