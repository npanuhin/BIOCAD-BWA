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
	  	  SCORE_THRESOLD = 3;

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

	ifstream input_database("samples/database.txt");
	vector<char> database = readSeq(input_database);
	int database_size = database.size();
	input_database.close();
	convertACGT(database);

	vector<int> degrees_of_4 = {1};
	for (int i = 1; i <= K; ++i) degrees_of_4.push_back(degrees_of_4[i - 1] * 4);

	vector<pair<int, vector<int>>> pos;
	int key = getKey(database.begin(), database.begin() + K - 1);
	key += 3 * degrees_of_4[K - 1];

	bool found;
	int found_ind;
	vector<int> tmp_v_int;

	for (int l = 0; l < database_size - K; ++l) {

		key = (key % degrees_of_4[K - 1]) * 4 + database[l + K - 1];

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

		// for (int f = 0; f < pos[i].second.size(); ++f) {
		// 	cout << pos[i].second[f] << ' ';
		// }
		// cout << endl;
	}

	ifstream input_query("samples/query.txt");
	vector<char> query = readSeq(input_query);
	int query_size = query.size();
	input_query.close();
	convertACGT(query);

	vector<int> matrix_line_one(query_size + 100), matrix_line_two(query_size + 100);
	matrix_line_one[0] = 0; matrix_line_two[0] = 0;
	vector<int>::iterator matrix_line_left = matrix_line_one.begin(), matrix_line_right = matrix_line_two.begin(), matrix_line_tmp;
	bool matrix_column_reached_score, up_set;

	int left, right, middle, y, x, up, bottom, highest_score, cur_value;
	// long long count = 0;

	cout << "Seeding..." << endl;

	ofstream result_file("result.txt");

	for (int l = 0; l < query_size - K; ++l) {

		// if (l % 10 == 0) {
		// 	cout << "l = " << l << endl;
		// }
		cout << "l = " << l << endl;

		vector<int> keys = {getKey(query.begin() + l, query.begin() + l + K)};

		for (int key_ind = 0; key_ind < (int)keys.size(); ++key_ind) {
			key = keys[key_ind];

			left = 0; right = pos.size();
			while (right - left > 1) {
				middle = (left + right) / 2;

				if (pos[middle].first == key) {
					for (int i = 0; i < (int)pos[middle].second.size(); ++i) {
						// cout << l << ' ' << pos[middle].second[i] << endl;
						y = l;
						x = pos[middle].second[i];
						up = y; bottom = y;
						highest_score = 0;

						*(matrix_line_left + bottom + 1) = 0;
						*(matrix_line_left + up) = 0;
						*(matrix_line_right + up) = 0;

						// cout << "x = " << x << " y = " << y << " | ";
						// ++count;

						result_file << x << ' ' << y << " -> ";

						matrix_column_reached_score = true;

						for (x = pos[middle].second[i]; matrix_column_reached_score && x < database_size; ++x) {
							matrix_column_reached_score = false;

							// cout << endl << "New x = " << x << endl;

							up_set = false;
							cur_value = -1;

							for (y = up; (cur_value != 0 || y != bottom) && y < query_size; ++y) {
								if (bottom == y) {
									*(matrix_line_left + bottom + 1) = 0;
									++bottom;
								}

								cur_value = max(max(
									*(matrix_line_left + y + 1) - 1, // (y + 1)
									*(matrix_line_right + y) - 1 // (y + 1) - 1
								),
									0
								);

								if (query[y] == database[x]) {
									cur_value = max(
										cur_value,
										*(matrix_line_left + y) + 1 // (y + 1) - 1
									);
								} else {
									cur_value = max(
										cur_value,
										*(matrix_line_left + y) - 1 // (y + 1) - 1
									);
								}

								if (cur_value <= highest_score - SCORE_THRESOLD) {
									cur_value = 0;
								} else {
									highest_score = max(highest_score, cur_value);
								}

								if (cur_value != 0) {
									// cout << '(' << y << ',' << cur_value << ')';
									matrix_column_reached_score = true;

									if (!up_set) {
										up_set = true;
										up = max(up, y);
									}
								}

								// cout << "y = " << y << endl << *(matrix_line_left + y) << ' ' << *(matrix_line_right + y) << endl << *(matrix_line_left + y + 1) << ' ' << cur_value << endl << endl;

								*(matrix_line_right + y + 1) = cur_value; // (y + 1)

								// if (cur_value == 0 && y + 1 == bottom) {
								// 	// cout << "break" << ' ';
								// 	// cout << '(' << *(matrix_line_right + y) << ')';
								// 	break;
								// }
							}

							matrix_line_tmp = matrix_line_left;
							matrix_line_left = matrix_line_right;
							matrix_line_right = matrix_line_tmp;

							// cout << x << '|' << matrix_column_reached_score << ' ';

							// if (!matrix_column_reached_score) break;

							// break;
						}

						result_file << x << ' ' << y << " (" << highest_score << ')' << endl << endl;

						// cout << endl << "------------------------------------------" << endl << endl;
					}
					break;

				} else if (key < pos[middle].first) {
					right = middle;
				} else {
					left = middle + 1;
				}
			}
		}
	}

	// cout << count << endl;

	// for (int i = 0; i < hssp_matches.size(); ++i) {
	// 	cout << hssp_matches[i].size() << endl;
	// }

	result_file.close();

	cout << "Test End" << endl << endl;
}


// Ищем геном1 в геноме2