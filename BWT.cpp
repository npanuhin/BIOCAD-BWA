#pragma GCC optimize("O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,tune=native")
#pragma GCC target("avx2")

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;
typedef long long LL;typedef long long ll;typedef unsigned long long uLL;typedef unsigned long long ull;typedef unsigned int uint;typedef unsigned char uchar;typedef long double ld;typedef long double Ld;

// template <class T>
// vector<T> cut(vector<T> &a, int start, int end) {vector<T> r;while(max(start, 0)<min((int)a.size(),end)){r.emplace_back(a[start++]);}return r;}
// #define endl '\n'

string& trim(string& str, const string& chars = "\t\n\v\f\r ") {
	str.erase(0, str.find_first_not_of(chars));
	str.erase(str.find_last_not_of(chars) + 1);
    return str;
}

// BWT for vector<char>
vector<char> bwt(vector<char> &s) {
	int n = s.size() + 1;

	s.push_back('$');

	vector<int> starts(n);
	for (int i = 0; i < n; ++i) starts[i] = i;

	sort(starts.begin(), starts.end(), [&s, n](int x, int y) {
		int res;
		for (int i = 0; i < n; ++i) {
			res = s[(x + i) % n] - s[(y + i) % n];
			if (res != 0) return res < 0;
		}
		return false;
	});

	vector<char> res(n);
	for (int i = 0; i < n; ++i) res[i] = s[(starts[i] + n - 1) % n];

	return res;
}

// The same for string
vector<char> bwt(string &s) {
	vector<char> tmp(s.begin(), s.end());
	return bwt(tmp);
}

// Inverse BWT for vector<char>
vector<char> ibwt(vector<char> &s) {
	int n = s.size(),
		k = distance(s.begin(), find(s.begin(), s.end(), '$'));

	vector<pair<char, int>> permutation(n);
	for (int i = 0; i < n; ++i) {
		permutation[i].first = s[i];
		permutation[i].second = i;
	}
	sort(permutation.begin(), permutation.end());

	vector<char> res(n - 1);
	for (int i = 0; i < n - 1; ++i) {
		res[i] = permutation[k].first;
		k = permutation[k].second;
	}

	return res;
}

// The same for string
vector<char> ibwt(string &s) {
	vector<char> tmp(s.begin(), s.end());
	return ibwt(tmp);
}

// Compression for vector<char>
vector<char> compress(vector<char> &s) {
	int n = s.size(), count = 1;

	vector<char> res {s[0]};

	for (int pos = 1; pos < n; ++pos, ++count) {
		if (s[pos - 1] != s[pos]) {
			if (count > 1) {
				string count_str = to_string(count);
				for (int i = 0; i < count_str.size(); ++i) {
					res.push_back(count_str[i]);
				}
			}
			count = 0;
			res.push_back(s[pos]);
		}
	}
	if (count > 1) {
		string count_str = to_string(count);
		for (int i = 0; i < count_str.size(); ++i) {
			res.push_back(count_str[i]);
		}
	}
	return res;
}

// The same for string
vector<char> compress(string &s) {
	vector<char> tmp(s.begin(), s.end());
	return compress(tmp);
}


int main() {

	ifstream input("src/BWT_test_large.txt");
	string line;
	vector<char> text_vector;

	cout << "Starting..." << endl;

	while (getline(input, line)) {
		trim(line);
		for (int i = 0; i < line.size(); ++i) {
			text_vector.push_back(line[i]);
		}
	}

	// text_vector = cut(text_vector, 0, (int)1e8);

	string text(text_vector.begin(), text_vector.end());

	cout << "Encrypting..." << endl;

	auto encrypted_vector = bwt(text_vector);

	cout << "Encrypted" << endl;

	string encrypted_text(encrypted_vector.begin(), encrypted_vector.end());

	cout << "Compressing..." << endl;
	
	auto compressed_vector = compress(encrypted_vector);

	cout << "Compressed" << endl;

	string compressed_text(compressed_vector.begin(), compressed_vector.end());

	cout << "Decrypting..." << endl;
	
	auto decrypted_vector = ibwt(encrypted_vector);

	cout << "Decrypted" << endl;

	string decrypted_text(decrypted_vector.begin(), decrypted_vector.end());


	// cout << text.size() << " -> " << encrypted_vector.size() << " -> " << compressed_vector.size() << endl;

	// cout << text << endl << endl;
	// cout << encrypted_text << endl << endl;
	// cout << compressed_text << endl << endl;
	// cout << decrypted_text << endl << endl;


	auto compressed_raw_vector = compress(text_vector);

	string compressed_raw_text(compressed_raw_vector.begin(), compressed_raw_vector.end());

	cout << text.size() << " -> " << encrypted_vector.size() << " -> " << compressed_vector.size() << " / " << compressed_raw_text.size() << endl;
}
