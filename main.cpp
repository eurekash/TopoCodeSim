#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
//#include "toric_code.h"
#include "toric_ancilla_block.h"

int main(int argc, char **argv) {
	int m,n,R;
	int T;
	double p;

	sscanf(argv[1], "%d", &m);
	sscanf(argv[2], "%d", &n);
	sscanf(argv[3], "%d", &R);
	sscanf(argv[4], "%lf", &p);
	sscanf(argv[5], "%d", &T);
	
	//ToricCodeCluster model(n, p);
	ToricAncBlk model(m, n, R, p);
	model.build_circuit();
	model.build_decoder_graph();
	int count = 0;

	for (int t = 0; t < T; t++) {
		if (model.simulate())  count++;
	}
	printf("%.8f\n", (double) count / T);

	return 0;
}
