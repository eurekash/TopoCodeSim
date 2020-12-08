#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include "toric_code.h"
//#include "toric_ancilla_block.h"

int main(int argc, char **argv) {
	//int , n, offset;
	int k1, k2;
	double p;
	int T;

	sscanf(argv[1], "%d", &k1);
	sscanf(argv[2], "%d", &k2);
	sscanf(argv[3], "%lf", &p);
	sscanf(argv[4], "%d", &T);
	
	int m = 3*k1, n = m*k2, nround = 2*n;
	ToricCodeBlock model(m, n, nround, k1, p);
	model.build_circuit();
	model.build_decoder_graph();

	int count = 0;

	for (int t = 0; t < T; t++) {
		if (model.simulate())  count++;
	}
	printf("%.8f\n", (double) count / T);

	return 0;
}
