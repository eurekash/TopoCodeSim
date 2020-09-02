#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include "toric_code.h"

int main(int argc, char **argv) {
	int n;
	int T;
	double p;

	sscanf(argv[1], "%d", &n);
	sscanf(argv[2], "%lf", &p);
	sscanf(argv[3], "%d", &T);
	
	ToricCodeCluster model(n, p);
	model.build_circuit();
	model.build_decoder_graph();
	int count = 0;

	for (int t = 0; t < T; t++) {
		if (model.simulate())  count++;
	}
	printf("%.8f\n", (double) count / T);

	return 0;
}
