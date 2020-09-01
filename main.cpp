#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include "toric_code.h"

int main() {
	int n = 5;
	double p = 1.5e-3;
	ToricCodeBareAncilla model(n, p);
	model.build_circuit();
	model.build_decoder_graph_one_round();
	return 0;
}
