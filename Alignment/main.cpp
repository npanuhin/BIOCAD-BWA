#pragma GCC optimize("O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,tune=native")
#pragma GCC target("avx2")

#include <iostream>
#include <fstream>
#include <vector>

#include "BWT.h"
#include "KNP.h"
#include "fm-index.h"
#include "k-mers.h"

using namespace std;
typedef long long LL;typedef long long ll;typedef unsigned long long uLL;typedef unsigned long long ull;typedef unsigned int uint;typedef unsigned char uchar;typedef long double ld;typedef long double Ld;

int main() {
	// bwtTest();
	// knpTest();
	// fmIndexTest();
	k_MersTest();
}