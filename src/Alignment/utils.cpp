#pragma GCC optimize("O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,tune=native")
#pragma GCC target("avx2")

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
typedef long long LL;typedef long long ll;typedef unsigned long long uLL;typedef unsigned long long ull;typedef unsigned int uint;typedef unsigned char uchar;typedef long double ld;typedef long double Ld;

string& trim(string& str, const string& chars = "\t\n\v\f\r ") {
	str.erase(0, str.find_first_not_of(chars));
	str.erase(str.find_last_not_of(chars) + 1);
    return str;
}

vector<char> readSeq(ifstream& input) {
	string line;
	vector<char> text;
	while (getline(input, line)) {
		trim(line);
		for (int i = 0; i < line.size(); ++i) {
			text.push_back(line[i]);
		}
	}
	return text;
}

vector<char> string2vector(string& s) {
	vector<char> res(s.begin(), s.end());
	return res;
}

string vector2string(vector<char>& a) {
	string res(a.begin(), a.end());
	return res;
}