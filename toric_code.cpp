#include "toric_code.h"
#include <cassert>
#include <algorithm>

#define DEBUG

inline double B(const double &p1, const double &p2) {
	return p1 * (1-p2) + p2 * (1-p1);
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

ToricCodeBareAncilla :: ToricCodeBareAncilla(int n, double p)
	: ToricCode(n, p)
{
}

void ToricCode :: build_decoder_one_round() {
	std::vector<int> X0, X1, Z0, Z1;
	double p_err;

	extractor->init_enumerator();

	while (extractor->enumerate(X0, X1, Z0, Z1, p_err)) {
		int u, v;
		std::pair<char, int> w;
		//X excitation
		assert(X0.size() + X1.size() <= 4);
		assert(X0.size() < 4);
		if (X0.size() == 2 && X1.size() == 2)  {  //time-like 4-point excitation
			u = X0[0],  v = X0[1],  w = std::make_pair(0, v);
#ifdef DEBUG
			assert(u < v);
#endif
			edges_x[u][w] = B(edges_x[u][w], p_err);
		} else if (X1.size() == 4) {   //space-like 4-point excitation
		
#ifdef DEBUG
			assert(X1[0] < X1[1]);
			assert(X1[1] < X1[2]);
			assert(X1[2] < X1[3]);
#endif
			u = X1[0],  v = X1[1],  w = std::make_pair(0, v);
			edges_x[u][w] = B(edges_x[u][w], p_err);

			u = X1[2],  v = X1[3],  w = std::make_pair(0, v);
			edges_x[u][w] = B(edges_x[u][w], p_err);
		} else {
			if (X0.size() == 2)  X1 = X0;
			if (X1.size() == 2) {    //space-like string error
				u = X1[0],  v = X1[1],  w = std::make_pair(0, v);
				edges_x[u][w] = B(edges_x[u][w], p_err);
			} else if (X0.size() == 1 && X1.size() == 1) {  //time-like string error
				u = X0[0],  v = X1[0],  w = std::make_pair(1, v);
				edges_x[u][w] = B(edges_x[u][w], p_err);
			}
		}
		//Z excitation
		assert(Z0.size() + Z1.size() <= 4);
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
		extractor->add_phase_flip(1, i, p);
		extractor->add_bit_flip(2, i, p);
	}
}


