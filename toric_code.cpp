#include "toric_code.h"
#include "flags.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cstdlib>

inline double H(const double &p) {
	return log( (1-p) / p );
}

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
	edges[0] = edges[1] = NULL;
	decoder[0] = decoder[1] = NULL;
}

void ToricCode :: build_circuit() {}

ToricCodeBareAncilla :: ToricCodeBareAncilla(int n, double p)
	: ToricCode(n, p)
{
}

ToricCodeBareAncillaNSWE :: ToricCodeBareAncillaNSWE(int n, double p)
	: ToricCode(n, p)
{
}

ToricCodeCluster :: ToricCodeCluster(int n, double p)
	: ToricCode(n, p)
{
}

void ToricCode :: build_decoder_graph_one_round() {
	std::vector<int> E0[2], E1[2];
	double p_err;
#ifdef DEBUG
	printf("build_decoder_graph_one_round()\n");
#endif

	extractor->init_enumerator();
	edges[0] = new EdgeTable[n2];
	edges[1] = new EdgeTable[n2];

	while (extractor->enumerate(E0[0], E1[0], E0[1], E1[1], p_err)) {
		int u, v;
		std::pair<char, int> w;
		for (int t = 0; t < 2; t++) {
			assert((E0[t].size() + E1[t].size()) % 2 ==0);
			if (E0[t].size() + E1[t].size() == 0)  continue;
			if (E0[t].size() == 2 && E1[t].size() == 2) {
				u = E0[t][0],  v = E0[t][1],  w = std::make_pair(1, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
				u = E1[t][0],  v = E1[t][1],  w = std::make_pair(1, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
			} else if (E0[t].size() == 4) {
				u = E0[t][0],  v = E0[t][1],  w = std::make_pair(0, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
				u = E0[t][2],  v = E0[t][3],  w = std::make_pair(0, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
			} else if (E1[t].size() == 4) {
				u = E1[t][0],  v = E1[t][1],  w = std::make_pair(0, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
				u = E1[t][2],  v = E1[t][3],  w = std::make_pair(0, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
			} else if (E0[t].size() == 2) {
				u = E0[t][0],  v = E0[t][1],  w = std::make_pair(0, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
			} else if (E1[t].size() == 2) {
				u = E1[t][0],  v = E1[t][1],  w = std::make_pair(0, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
			} else if (E0[t].size() == 1 && E1[t].size() == 1) {
				u = E0[t][0],  v = E1[t][0],  w = std::make_pair(1, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
			}
		}
	}

	/*
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (auto &e: edges[0][I(n,i,j)]) {
				printf("(0,%d,%d) <-> (%d,%d,%d), p = %.6f\n", i, j, e.first.first, e.first.second / n, e.first.second % n, e.second);
			}
		}
	}

	*/
}

void ToricCode :: build_decoder_graph() 
{
	build_decoder_graph_one_round();
	for (int o = 0; o < 2; o++) {
		decoder[o] = new Decoder( std::vector<bool> ( (n+1)*n*n, false ) );
		for (int t = 0; t <= n; t++) {
			for (int u = 0; u < n2; u++) {
				for (auto &p: edges[o][u]) {
					int tp = p.first.first;
					int v  = p.first.second;
					double w = H(p.second);
					if (t+tp <= n) {
						decoder[o]->add_edge(I(n2,t,u), I(n2,t+tp,v), w);
					}
				}
			}
		}
	}
}

bool ToricCode :: simulate()
{
	std::vector<int>  data[2];
	std::vector<int>  anc[2];
	std::vector<int>  anc_last_round[2];

	for (int o = 0; o < 2; o++) {
		decoder[o]->reset();
		data[o].resize(2*n2, false);
		anc[o].resize(n2, false);
		anc_last_round[o].resize(n2, false);
	}

	int counter[2] = {0,0};

	for (int t = 0; t <= n; t++) {
		extractor->execute(data[0], data[1], anc[0], anc[1], t < n);
		for (int o = 0; o < 2; o++) {
			for (int i = 0; i < n2; i++) {
				if (anc_last_round[o][i] ^ anc[o][i]) {
					decoder[o]->excite(I(n2, t, i));
					++counter[o];
				}
			}
			anc_last_round[o] = anc[o];
		}
	}

	//check logical errors
	int cut_hori = 1,  cut_vert = n-1;
	int offset_hori = 0,  offset_vert = n2;

	for (int o = 0; o < 2; o++) {
		bool err_hori = false, err_vert = false;
		std::vector < std::pair<int, int> > correct = decoder[o]->decode();

		for (int i = 0; i < n; i++) {
			if (data[o^1][ I(n, i, 0) + offset_hori ])   err_hori = !err_hori;
			if (data[o^1][ I(n, 0, i) + offset_vert ])   err_vert = !err_vert;
		}

		for (auto &p: correct) {
			int u = p.first % n2,  v = p.second % n2,
				x1 = u / n,  y1 = u % n,  x2 = v / n,  y2 = v % n;

			if ((y1 == 0 && y2 == cut_hori) || (y1 == cut_hori && y2 == 0))   
				err_hori = !err_hori;
			if ((x1 == 0 && x2 == cut_vert) || (x1 == cut_vert && x2 == 0))   
				err_vert = !err_vert;
		}

		if (err_hori || err_vert)  return true;

		std::swap(offset_hori, offset_vert);
		std::swap(cut_hori, cut_vert);
	}

	return false;
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
	
	/*
	for (int i = 0; i < n2; i++) {
		extractor->add_1qubit_depo(1, i, p);
		extractor->add_1qubit_depo(2, i, p);
	}
	*/
	/*
	for (int i = 0; i < n2; i++) {
		extractor->add_phase_flip(1, i, p/10);
		extractor->add_bit_flip(2, i, p/10);
	}
	*/
	/*
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(0, i, p);
	}*/

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
	/*
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(0, i, p);
	}*/
}

void ToricCodeBareAncillaNSWE :: build_circuit() {
	extractor = new Extractor(2*n2, n2, n2);

	for (int i = 0; i < n2; i++) {
		extractor->set_syndrome_bit_x(i);
	}

	for (int i = 0; i < n2; i++) {
		extractor->set_syndrome_bit_z(i);
	}

	//Init
	/*
	for (int i = 0; i < n2; i++) {
		extractor->add_phase_flip(1, i, p/10);
		extractor->add_bit_flip(2, i, p/10);
	}
	*/
	/*
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(0, i, p);
	}*/

	//X North
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, I(n,i,j)+n2);
			extractor->add_2qubit_depo(1, I(n,i,j), 0, I(n,i,j)+n2, p);
		}
	}
	//X South
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, D(n,i,j)+n2);
			extractor->add_2qubit_depo(1, I(n,i,j), 0, D(n,i,j)+n2, p);
		}
	}
	//X West
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, L(n,i,j));
			extractor->add_2qubit_depo(1, I(n,i,j), 0, L(n,i,j), p);
		}
	}
	//X East
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, I(n,i,j));
			extractor->add_2qubit_depo(1, I(n,i,j), 0, I(n,i,j), p);
		}
	}

	//Z north
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(0, U(n,i,j), 2, I(n,i,j));
			extractor->add_2qubit_depo(0, U(n,i,j), 2, I(n,i,j), p);
		}
	}
	//Z south
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(0, I(n,i,j), 2, I(n,i,j));
			extractor->add_2qubit_depo(0, I(n,i,j), 2, I(n,i,j), p);
		}
	}

	//Z west
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(0, I(n,i,j)+n2, 2, I(n,i,j));
			extractor->add_2qubit_depo(0, I(n,i,j)+n2, 2, I(n,i,j), p);
		}
	}

	//Z east
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(0, R(n,i,j)+n2, 2, I(n,i,j));
			extractor->add_2qubit_depo(0, R(n,i,j)+n2, 2, I(n,i,j), p);
		}
	}

	//measure
	for (int i = 0; i < n2; i++) {
		extractor->add_phase_flip(1, i, p);
		extractor->add_bit_flip(2, i, p);
	}
	/*
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(0, i, p);
	}*/
}

void ToricCodeCluster :: build_circuit() 
{
	extractor = new Extractor(2*n2, 2*n2, 2*n2);

	for (int i = 0; i < n2; i++) {
		extractor->set_syndrome_bit_x(i);
		extractor->set_syndrome_bit_z(i);
	}



	//Cluster states for X stabilizers
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 1, I(n,i,j)+n2);
			extractor->add_2qubit_depo(1, I(n,i,j), 1, I(n,i,j)+n2, p/10);
		}
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 1, L(n,i,j)+n2);
			extractor->add_2qubit_depo(1, I(n,i,j), 1, L(n,i,j)+n2, p/10);
		}
	}
	//Cluster states for Z stabilizers
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(2, R(n,i,j)+n2, 2, I(n,i,j));
			extractor->add_2qubit_depo(2, R(n,i,j)+n2, 2, I(n,i,j), p/10);
		}
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(2, I(n,i,j)+n2, 2, I(n,i,j));
			extractor->add_2qubit_depo(2, I(n,i,j)+n2, 2, I(n,i,j), p/10);
		}
	}

	/*
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(1, i, p);
		extractor->add_1qubit_depo(2, i, p);
	}
	*/

	//X stabilizer transversal gates
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j)+n2, 0, I(n,i,j));
			extractor->add_2qubit_depo(1, I(n,i,j)+n2, 0, I(n,i,j), p);
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, I(n,i,j)+n2);
			extractor->add_2qubit_depo(1, I(n,i,j), 0, I(n,i,j)+n2, p);
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 0, D(n,i,j)+n2);
			extractor->add_2qubit_depo(1, I(n,i,j), 0, D(n,i,j)+n2, p);
		}
	}
	//Z stabilizer transversal gates
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(0, I(n,i,j)+n2, 2, I(n,i,j)+n2);
			extractor->add_2qubit_depo(0, I(n,i,j)+n2, 2, I(n,i,j)+n2, p);
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(0, U(n,i,j), 2, I(n,i,j));
			extractor->add_2qubit_depo(0, U(n,i,j), 2, I(n,i,j), p);
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(0, I(n,i,j), 2, I(n,i,j));
			extractor->add_2qubit_depo(0, I(n,i,j), 2, I(n,i,j), p);
		}
	}

	//measurement errors
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_phase_flip(1, i, p/10);
		extractor->add_bit_flip(2, i, p/10);
	}

	//X stabilizers ideal decoding
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 1, L(n,i,j)+n2);
		}
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(1, I(n,i,j), 1, I(n,i,j)+n2);
		}
	}
	//Z stabilizers ideal decoding
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(2, I(n,i,j)+n2, 2, I(n,i,j));
		}
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(2, R(n,i,j)+n2, 2, I(n,i,j));
		}
	}
}

/*
void ToricCode3qubitCat :: build_circuit() {
	assert(n % 2 == 0);
	extractor = new Extractor(2*n2, n2*3/2, n2*3/2);


}
*/

