#include "toric_code.h"
#include "flags.h"
#include <cassert>
#include <algorithm>


inline double B(const double &p1, const double &p2) {
//	return p1 * (1-p2) + p2 * (1-p1);
	return p1 + p2;
}

inline int I(int n, int i, int j) {
	return n * i + j;
}

inline int L(int n, int i, int j) {
	return j > 0? I(n, i, j-1): I(n, i, n-1);
}

inline int R(int n, int i, int j) {
	return j < n-1? I(n, i, j+1): I(n, i, 0);
}

inline int U(int n, int i, int j) {
	return i > 0? I(n, i-1, j): I(n, n-1, j);
}

inline int D(int n, int i, int j) {
	return i < n-1? I(n, i+1, j): I(n, 0, j);
}

ToricCode :: ToricCode(int n, double p) {
	this->n = n;
	this->p = p;
	n2 = n * n;
	extractor = NULL;
	edges_x = new EdgeTable[n2];
	edges_z = new EdgeTable[n2];
}

void ToricCode :: build_circuit() {}

ToricCodeBareAncilla :: ToricCodeBareAncilla(int n, double p)
	: ToricCode(n, p)
{
}

void ToricCode :: build_decoder_graph_one_round() {
	std::vector<int> X0, X1, Z0, Z1;
	double p_err;
	printf("build_decoder_graph_one_round()\n");

	extractor->init_enumerator();

	while (extractor->enumerate(X0, X1, Z0, Z1, p_err)) {
		int u, v;
		std::pair<char, int> w;
		//X excitation
		if (X0.size() + X1.size() > 0) {
			if (X0.size() >= 2)   X1 = X0;
			if (X1.size() >= 2)  {
				u = X1[0],  v = X1[1],  w = std::make_pair(0, v);
			} else if (X0.size() == 1 && X1.size() == 1) {  //time-like string error
				u = X0[0],  v = X1[0],  w = std::make_pair(1, v);
			}
			edges_x[u][w] = B(edges_x[u][w], p_err);
		}
		//Z excitation
		if (Z0.size() + Z1.size() > 0) {
			if (Z0.size() >= 2)   Z1 = Z0;
			if (Z1.size() >= 2)  {
				u = Z1[0],  v = Z1[1],  w = std::make_pair(0, v);
			} else if (Z0.size() == 1 && Z1.size() == 1) {  //time-like string error
				u = Z0[0],  v = Z1[0],  w = std::make_pair(1, v);
			}
			edges_z[u][w] = B(edges_z[u][w], p_err);
		}
	}

#ifdef DEBUG
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (auto &e: edges_x[I(n,i,j)]) {
				printf("(0,%d,%d) <-> (%d,%d,%d), p = %.6f\n", i, j, e.first.first, e.first.second / n, e.first.second % n, e.second);
			}
		}
	}

#endif
}

void ToricCodeBareAncilla :: build_circuit() {
	extractor = new Extractor(2*n2, n2, n2);

	for (int i = 0; i < n2; i++) {
		extractor->set_syndrome_bit_x(i);
	}

	for (int i = 0; i < n2; i++) {
		extractor->set_syndrome_bit_z(i);
	}

	//Init
	for (int i = 0; i < n2; i++) {
		extractor->add_phase_flip(1, i, 2*p/3);
		extractor->add_bit_flip(2, i, 2*p/3);
	}
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(0, i, p);
	}

	//North
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, I(n,i,j)+n2);
			extractor->add_2qubit_depo(1, I(n,i,j), 0, I(n,i,j)+n2, p);
			extractor->add_CNOT(0, U(n,i,j), 2, I(n,i,j));
			extractor->add_2qubit_depo(0, U(n,i,j), 2, I(n,i,j), p);
		}
	}

	//West
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, L(n,i,j));
			extractor->add_2qubit_depo(1, I(n,i,j), 0, L(n,i,j), p);
			extractor->add_CNOT(0, I(n,i,j)+n2, 2, I(n,i,j));
			extractor->add_2qubit_depo(0, I(n,i,j)+n2, 2, I(n,i,j), p);
		}
	}

	//East
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, I(n,i,j));
			extractor->add_2qubit_depo(1, I(n,i,j), 0, I(n,i,j), p);
			extractor->add_CNOT(0, R(n,i,j)+n2, 2, I(n,i,j));
			extractor->add_2qubit_depo(0, R(n,i,j)+n2, 2, I(n,i,j), p);
		}
	}

	//South
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, D(n,i,j)+n2);
			extractor->add_2qubit_depo(1, I(n,i,j), 0, D(n,i,j)+n2, p);
			extractor->add_CNOT(0, I(n,i,j), 2, I(n,i,j));
			extractor->add_2qubit_depo(0, I(n,i,j), 2, I(n,i,j), p);
		}
	}

	//measure
	for (int i = 0; i < n2; i++) {
		extractor->add_phase_flip(1, i, 2*p/3);
		extractor->add_bit_flip(2, i, 2*p/3);
	}
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(0, i, p);
	}
}


