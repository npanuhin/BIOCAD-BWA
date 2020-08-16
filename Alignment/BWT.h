#ifndef BWT_H
#define BWT_H

// BWT for vector<char>
std::vector<char> bwt(std::vector<char> &s);

// Inverse BWT for vector<char>
std::vector<char> ibwt(std::vector<char> &s);

// Compression for vector<char>
std::vector<char> compress(std::vector<char> &s);

// For testing
int bwtTest();

#endif