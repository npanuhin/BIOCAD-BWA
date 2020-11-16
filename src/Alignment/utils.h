#ifndef UTILS_H
#define UTILS_H

std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ");

std::vector<char> readSeq(std::ifstream& input);

std::vector<char> string2vector(std::string& s);

std::string vector2string(std::vector<char>& a);

#endif
