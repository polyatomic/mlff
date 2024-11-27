#pragma once

#include <vector>
#include <string>

using std::vector;
using std::string;

#define MAX_THREADS 16

void sketch_matrices(int N, vector<string>& mat_names, int *ids_selected_total, int ranks[], int n_test_rows, int nthreads, int *selected_variables[], int max_cols, bool select_from_leading_matrix);
