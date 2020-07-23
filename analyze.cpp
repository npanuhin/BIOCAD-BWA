#pragma GCC optimize("O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,tune=native")
#pragma GCC target("avx2")

// #include <bits/stdc++.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <cmath>
#include <set>
#include <map>
#include <iomanip>
#include <cstring>
#include <iterator>
#include <bitset>
#include <deque>
#include <queue>
#include <limits.h>

using namespace std;
typedef long long LL;typedef long long ll;typedef unsigned long long uLL;typedef unsigned long long ull;typedef unsigned int uint;typedef unsigned char uchar;typedef long double ld;typedef long double Ld;

template <class T>
T _read(){T x;cin>>x;return x;}
template <class T>
long long sum(vector<T> &a){long long x=0;for(int _i=0;_i<a.size();++_i)x+=a[_i];return x;}
template<class InputIt, class T>
typename iterator_traits<InputIt>::difference_type
rcount(InputIt first, InputIt last, const T& value){typename iterator_traits<InputIt>::difference_type ret=0;for(;first!=last;++first){if(*first!=value)ret++;}return ret;}

#define all(_x) (_x).begin(), (_x).end()
#define rall(_x) (_x).end(), (_x).begin()
#define sortv(_x) sort((_x).begin(), (_x).end())
#define get(_type, _x) _type _x; cin >> (_x);
#define print(_x) cout << (_x) << endl;
#define read(_type) _read <_type> ()
#define maxv(_v) (max_element(_v.begin(), _v.end()) - _v.begin())
#define minv(_v) (min_element(_v.begin(), _v.end()) - _v.begin())

bool contains(string &s, char x) {for(size_t i=0;i<s.size();i++){if(s[i]==x){return true;}}return false;}
template <class T>
bool contains(vector<T> &a, T x) {for (size_t i=0;i<a.size();i++){if(a[i]==x){return true;}}return false;}
string cut(string &s, int start, int end) {string r;while(start<end){r+=s[start++];}return r;}
template <class T>
vector<T> cut(vector<T> &a, int start, int end) {vector<T> r;while(start<end){r.emplace_back(a[start++]);}return r;}

int gcd(int a, int b) {while(b){a%=b;swap(a,b);}return a;}
long long gcd(long long a, long long b) {while(b){a%=b;swap(a,b);}return a;}

int gcd(int a, int b, int &x, int &y) {if(a == 0){x = 0;y = 1;return b;}int x1, y1;int d = gcd (b%a, a, x1, y1);x = y1 - (b / a) * x1;y = x1;return d;}
long long gcd(long long a, long long b, long long &x, long long &y) {if(a == 0){x = 0;y = 1;return b;}long long x1, y1;long long d = gcd (b%a, a, x1, y1);x = y1 - (b / a) * x1;y = x1;return d;}

int lcm(int a, int b) {return a/gcd(a, b)*b;}
long long lcm(long long a, long long b) {return a/gcd(a, b)*b;}

long long binpow(long long a, unsigned int n) {long long r=1;while(n){if(n&1)r*=a;a*=a;n>>=1;}return r;}

bool isPrime (long long num) {
	if(num<=1||num%2==0){return false;}int upperLimit=static_cast<int>(sqrt((double)static_cast<double>(num)) + 1);
	for (int divisor=3;divisor<=upperLimit;divisor+=2){if(num%divisor==0){return false;}}return true;
}

// #define endl '\n'

int main() {
	ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);

	// string genome1_path = "samples/small/source.fasta", genome2_path = "samples/small/deletion.fasta";
	string genome1_path = "samples/large/large_genome1.fasta", genome2_path = "samples/large/large_genome2.fasta";

	// READING

	ifstream genome1_file(genome1_path);
	ifstream genome2_file(genome2_path);

	vector<char> genome1;
	vector<char> genome2;

	string line;
	int line_num;

	line_num = 0;
	while (getline(genome1_file, line)) {
		if (line_num > 0) {
			for (size_t i = 0; i < line.size(); ++i) {
				genome1.push_back(line[i]);
			}
		}
		++line_num;
	}

	line_num = 0;
	while (getline(genome2_file, line)) {
		if (line_num > 0) {
			for (size_t i = 0; i < line.size(); ++i) {
				genome2.push_back(line[i]);
			}
		}
		++line_num;
	}

	genome1_file.close();
	genome2_file.close();

	int n1 = genome1.size(), n2 = genome2.size(), length = max(n1, n2);

	cout << "Reading compete" << endl << "Genome1: " << n1 << endl << "Genome2: " << n2 << endl << "length: " << length << endl;

	//------------------------------------------------

	bool tmp_matrix_array[(int)1e4] = {false};
	bool matrix_array[(int)1e4][(int)1e4] = {false};

	// vector<vector<bool>> matrix(n1, tmp_matrix_vector);

	// for (int i = 0; i < n1; ++i) {
	// 	for (int j = 0; j < n2; ++j) {
	// 		if (genome1[i] == genome2[i]) {
	// 			matrix[i][j] = true;
	// 		}
	// 	}
	// }

	// ofstream fout("cppstudio.txt");
	// fout << "Работа с файлами в С++";
	// fout.close();
}
