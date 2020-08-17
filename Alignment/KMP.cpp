#pragma GCC optimize("O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,tune=native")
#pragma GCC target("avx2")

#include <iostream>
#include <fstream>
#include <vector>

#include "utils.h"

using namespace std;
typedef long long LL;typedef long long ll;typedef unsigned long long uLL;typedef unsigned long long ull;typedef unsigned int uint;typedef unsigned char uchar;typedef long double ld;typedef long double Ld;

vector<int> prefixFunction(vector<char>& s) {
	int n = (int)s.size();
	vector<int> pi(n);
	for (int i = 1; i < n; ++i) {
		int j = pi[i - 1];
		while (j > 0 && s[i] != s[j])
			j = pi[j - 1];
		if (s[i] == s[j])  ++j;
		pi[i] = j;
	}
	return pi;
}

vector<int> findAllSubstrings(vector<char>& query, vector<char>& text) {
	int n1 = (int)query.size(), n1_2 = n1 * 2, n2 = (int)text.size(), n3 = n1 + n2 + 1;

	vector<char> prefix_fun_text(query.begin(), query.end());
	prefix_fun_text.reserve(n3);
	prefix_fun_text.push_back('#');
	prefix_fun_text.insert(prefix_fun_text.end(), text.begin(), text.end());
	vector<int> prefix_fun = prefixFunction(prefix_fun_text);

	vector<int> res;
	for (int i = n1; i < n3; ++i) {
		if (prefix_fun[i] == n1) res.push_back(i - n1_2);
	}
	return res;
}

int kmpTest() {
	ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);

	cout << "Starting KMP test..." << endl;

	ifstream input("../src/BWT/BWT_test_large.txt");
	vector<char> text = readSeq(input);
	input.close();

	vector<char> search = {'T', 'T', 'A'};

	vector<int> substrings = findAllSubstrings(search, text);

	cout << "Test End" << endl << endl;
}