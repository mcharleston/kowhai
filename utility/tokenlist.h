/*
 * tokenlist.h
 *
 *  Created on: 2 Jun 2016
 *      Author: mc7
 */

#ifndef TOKENLIST_H_
#define TOKENLIST_H_

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "token.h"

namespace parsing {

class TokenList {
private:
	unsigned int lineNumber;
	char ignoreLineChar = '%';
	char beginCommentChar = '[';
	char endCommentChar = ']';
	std::string idchars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
public:
	std::vector<parsing::Token> toks;
	std::vector<parsing::Token>::iterator current;

	virtual ~TokenList();
	TokenList() : lineNumber(0) {}
	TokenList(int argn, char** argc);

	void advance();

	inline std::vector<parsing::Token>::const_iterator begin() const { return toks.begin(); }
	inline std::vector<parsing::Token>::iterator begin() { return toks.begin(); }
	inline std::vector<parsing::Token>::const_iterator end() const { return toks.end(); }
	inline std::vector<parsing::Token>::iterator end() { return toks.end(); }

	inline bool isBool() const { return (current->getType() == parsing::boolT); }
	inline bool isChar() const { return (current->getType() == parsing::charT); }
	inline bool isComment() const { return (current->getType() == parsing::commentT); }
	inline bool isInt() const { return (current->getType() == parsing::intT); }
	inline bool isDouble() const { return (current->getType() == parsing::doubleT); }
	inline bool isString() const { return (current->getType() == parsing::stringT); }
	inline bool isIntRange() const { return (current->getType() == parsing::intRangeT); }

	inline char getBeginCommentChar() const { return beginCommentChar; }
	inline char getEndCommentChar() const { return endCommentChar; }

	inline bool hasNext() const { return (current+1 != toks.end()); }
	inline bool hasNextBool() const { return (hasNext() && (current+1)->type == parsing::boolT); }
	inline bool hasNextChar() const { return (hasNext() && (current+1)->type == parsing::charT); }
	inline bool hasNextInt() const { return (hasNext() && (current+1)->type == parsing::intT); }
	inline bool hasNextDouble() const { return (hasNext() && (current+1)->type == parsing::doubleT); }

	TokenList& operator=(const TokenList& tl);
	inline int numTokens() const { return toks.size(); }

	inline void reset() { current = toks.begin(); }


	inline void setBeginCommentChar(char ch) { beginCommentChar = ch; }
	inline void setEndCommentChar(char ch) { endCommentChar = ch; }

	void tokenize(const std::string& fname);
	void tokenize(std::istream& fin);
	void tokenize(std::vector<std::string>& args);
};

} // namespace parsing

#endif /* TOKENLIST_H_ */
