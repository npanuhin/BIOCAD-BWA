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

#include "utils.h"

using namespace std;

const int K = 3,
	  	  HSSP = 1,
	  	  SCORE_THRESOLD = 1;  // Warning! Others will not work!!!

int getKey(vector<char>::iterator start, vector<char>::iterator end) {
	int res = 0;
	for (vector<char>::iterator i = start; i < end; ++i) {
		res *= 4;
		res += (int)*i;
	}
	return res;
}

vector<char> reverseKey(int key) {
	vector<char> res;
	int tmp_key = key;
	for (int i = 0; i < K; ++i) {
		res.push_back((char)(tmp_key % 4));
		tmp_key /= 4;
	}
	reverse(res.begin(), res.end());
	return res;
}

void convertACGT(vector<char>& text) {
	for (int i = 0; i < (int)text.size(); ++i) {
		if (text[i] == 'A') {
			text[i] = 0;
		} else if (text[i] == 'C') {
			text[i] = 1;
		} else if (text[i] == 'G') {
			text[i] = 2;
		} else if (text[i] == 'T') {
			text[i] = 3;
		}
	}
}

int k_MersTest() {
	ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);

	cout << "Starting k-mers test..." << endl;

	ifstream input_text2("samples/genome2.txt");
	vector<char> text = readSeq(input_text2);
	input_text2.close();
	convertACGT(text);

	vector<int> degrees_of_4 = {1};
	for (int i = 1; i <= K; ++i) degrees_of_4.push_back(degrees_of_4[i - 1] * 4);

	vector<pair<int, vector<int>>> pos;
	int key = getKey(text.begin(), text.begin() + K - 1);
	key += 3 * degrees_of_4[K - 1];

	bool found;
	int found_ind;
	vector<int> tmp_v_int;

	for (int l = 0; l < (int)text.size() - K; ++l) {

		key = (key % degrees_of_4[K - 1]) * 4 + text[l + K];

		found = false;
		for (found_ind = 0; found_ind < (int)pos.size(); ++found_ind) {
			if (pos[found_ind].first == key) {
				found = true;
				break;
			}
		}
		if (!found) {
			found_ind = pos.size();
			pos.emplace_back(key, tmp_v_int);
		}
		pos[found_ind].second.push_back(l);
	}

	sort(pos.begin(), pos.end());

	cout << "pos size: " << pos.size() << endl;

	for (int i = 0; i < (int)pos.size(); ++i) {
		cout << pos[i].first << ' ';
		vector<char> reversed_key = reverseKey(pos[i].first);

		for (int j = 0; j < (int)reversed_key.size(); ++j) {
			if (reversed_key[j] == 0) {
				reversed_key[j] = 'A';
			} else if (reversed_key[j] == 1) {
				reversed_key[j] = 'C';
			} else if (reversed_key[j] == 2) {
				reversed_key[j] = 'G';
			} else if (reversed_key[j] == 3) {
				reversed_key[j] = 'T';
			}
			cout << reversed_key[j];
		}
		cout << ' ' << pos[i].second.size() << endl;
	}

	ifstream input_text1("samples/genome1.txt");
	vector<char> text2 = readSeq(input_text1);
	input_text1.close();
	convertACGT(text2);

	int left, right, middle;

	// vector<vector<int>> hssp_matches(text2.size() - K);
	// vector<pair<int, int>> seeds;

	for (int l = 0; l < (int)text2.size() - K; ++l) {
		vector<int> keys = {getKey(text2.begin() + l, text2.begin() + l + K)};

		for (int key_ind = 0; key_ind < (int)keys.size(); ++key_ind) {
			key = keys[key_ind];

			left = 0; right = pos.size();
			while (right - left > 1) {
				middle = (left + right) / 2;

				if (pos[middle].first == key) {
					for (int i = 0; i < (int)pos[middle].second.size(); ++i) {
						// cout << l << ' ' << pos[middle].second[i] << endl;
						// seeds.emplace_back(l, pos[middle].second[i]);
						int x = l, y = pos[middle].second[i], cur_score = 0;
						bool reached_score = false;

						while (!reached_score || cur_score >= SCORE_THRESOLD) {
							// TODO
						}
					}
					
					// hssp_matches[l].insert(hssp_matches[l].end(), pos[middle].second.begin(), pos[middle].second.end());
					break;

				} else if (key < pos[middle].first) {
					right = middle;
				} else {
					left = middle + 1;
				}
			}
		}
	}

	// for (int i = 0; i < hssp_matches.size(); ++i) {
	// 	cout << hssp_matches[i].size() << endl;
	// }

	cout << "Test End" << endl << endl;
}
