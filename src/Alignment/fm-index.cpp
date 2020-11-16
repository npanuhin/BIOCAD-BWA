#pragma GCC optimize("O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,tune=native")
#pragma GCC target("avx2")

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

#include "utils.h"
#include "BWT.h"

using namespace std;
typedef long long LL;typedef long long ll;typedef unsigned long long uLL;typedef unsigned long long ull;typedef unsigned int uint;typedef unsigned char uchar;typedef long double ld;typedef long double Ld;

int fmIndexTest() {
	ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);

	cout << "Starting FM-index test..." << endl;

	ifstream input("../src/BWT/BWT_test_large.txt");
	vector<char> text = readSeq(input);
	input.close();

	cout << "Test End" << endl << endl;;
}
