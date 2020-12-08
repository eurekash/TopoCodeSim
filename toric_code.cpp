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

inline int offset(int n, int u, int offset) {
	int i = u / n, j = u % n;
	i = (i + offset) % n;
	j = (j + offset) % n;
	return I(n, i, j);
}

ToricCode :: ToricCode(int n, double p, int D, int offset_per_round) {
	this->n = n;
	this->p = p;
	this->D = D;
	this->offset_per_round = offset_per_round;
	n2 = n * n;
	extractor = NULL;
	edges[0] = edges[1] = NULL;
	decoder[0] = decoder[1] = NULL;
}

void ToricCode :: build_circuit() {}

/*
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
*/

ToricCodeBlock :: ToricCodeBlock(int m, int n, int D, int offset, double p)
	: ToricCode(n, p, D, offset)
{
	this->m = m;
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
				u = E0[t][0],  v = E0[t][1],  w = std::make_pair(0, v);
				edges[t][u][w] = B(edges[t][u][w], p_err);
				u = E1[t][0],  v = E1[t][1],  w = std::make_pair(0, v);
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
}

void ToricCode :: build_decoder_graph() 
{
	build_decoder_graph_one_round();
	for (int o = 0; o < 2; o++) {
		decoder[o] = new Decoder( std::vector<bool> ( (D+1)*n2, false ) );
		decoder[o]->init_shortest_path(n2);
		for (T = 0; ; T++) {
			bool success = decoder[o]->shortest_path(T*n2, D);
			int offset_t = T * offset_per_round;
			for (int u = 0; u < n2; u++) {
				int u1 = offset(n, u, offset_t);
				for (auto &p: edges[o][u]) {
					int tp = p.first.first;
					int v  = p.first.second;
					double w = H(p.second);
					int v1 = offset(n, v, offset_t);
					if (tp == 0 || !success) {
						decoder[o]->add_edge(I(n2, T, u1), I(n2, T+tp, v1), w);
					}
				}
			}
			if (success)  break;
		}
	}
}

bool ToricCode :: simulate()
{
	std::vector<int>  data[2], data_offset[2];
	std::vector<int>  anc[2], anc_offset[2];
	std::vector<int>  anc_last_round[2];

	for (int o = 0; o < 2; o++) {
		decoder[o]->reset();
		data[o].resize(2*n2, 0);
		data_offset[o].resize(2*n2, 0);
		anc[o].resize(n2, 0);
		anc_offset[o].resize(n2, 0);
		anc_last_round[o].resize(n2, 0);
	}

	for (int t = 0; t <= T; t++) {
		//Offset the torus
		int offset_t = (n-t*offset_per_round%n)%n;
		for (int o = 0; o < 2; o++) {
			for (int i = 0; i < n2; i++) {
				int j = offset(n, i, offset_t);
				data_offset[o][j] = data[o][i];
				data_offset[o][j+n2] = data[o][i+n2];
			}
		}
		extractor->execute(data_offset[0], data_offset[1], anc_offset[0], anc_offset[1], t < T);
		//shift back
		for (int o = 0; o < 2; o++) {
			for (int i = 0; i < n2; i++) {
				int j = offset(n, i, offset_t);
				data[o][i] = data_offset[o][j];
				data[o][i+n2] = data_offset[o][j+n2];
				anc[o][i] = anc_offset[o][j];
			}
		}
		for (int o = 0; o < 2; o++) {
			for (int i = 0; i < n2; i++) {
				if (anc_last_round[o][i] ^ anc[o][i]) {
					decoder[o]->excite(I(n2, t, i));
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

void ToricCodeBlock :: build_circuit() {
	int n_anc_qubits = 2*m*(m+1)*(n/m)*(n/m);
	extractor = new Extractor(
			2*n2,  //2n^2 data qubits
			n2 + n_anc_qubits,  //x ancilla qubits
			n2 + n_anc_qubits   //z ancilla qubits
		);


	for (int i = 0; i < n2; i++) {
		extractor->set_syndrome_bit_x(i);
		extractor->set_syndrome_bit_z(i);
	}

	
	//init ancilla error
	
	for (int i = 0; i < n_anc_qubits; i++) {
		extractor->add_1qubit_depo(1, i+n2, p);
		extractor->add_1qubit_depo(2, i+n2, p);
	}

	//Z stabilizer measurements

	int B = m*(m+1);

	for (int offsetx = 0, offset_anc = n2; offsetx < n; offsetx += m) {
		//offset_anc: offset of ancilla qubit id
		//first n2 qubits are syndrome bits
		for (int offsety = 0; offsety < n; offsety += m) {
			//horizontal
			for (int x = 0; x <= m; x++) {
				for (int y = 0; y < m; y++) {
					//data qubit: (x+offsetx, y+offsety)
					//ancilla qubit: I(m,x,y)
					int dat_id = U(n, (x+offsetx)%n, y+offsety);
					int anc_id = I(m,x,y) + offset_anc;
					//printf("[Z HORIZONTAL] data: %d(%d,%d), ancilla: %d(%d,%d)\n", dat_id, (x+offsetx)%n, (y+offsety)%n, anc_id, x, y);
					extractor->add_CNOT(0, dat_id, 2, anc_id);
					extractor->add_2qubit_depo(0, dat_id, 2, anc_id, p);
					//syndrome bit evaluation
					if (x > 0) {
						extractor->add_CNOT(2, anc_id, 2, I(n, x+offsetx-1, y+offsety));
					}
					if (x < m) {
						extractor->add_CNOT(2, anc_id, 2, I(n, x+offsetx, y+offsety));
					}
				}
			}
			offset_anc += B;
			//vertical
			for (int x = 0; x < m; x++) {
				for (int y = 0; y <= m; y++)  {
					int dat_id = I(n, x+offsetx, (y+offsety)%n) + n2;
					int anc_id = I(m+1, x, y) + offset_anc;
					//printf("[Z VERTICAL] data: %d(%d,%d), ancilla: %d(%d,%d)\n", dat_id, (x+offsetx)%n, (y+offsety)%n, anc_id, x, y);
					extractor->add_CNOT(0, dat_id, 2, anc_id);
					extractor->add_2qubit_depo(0, dat_id, 2, anc_id, p);
					//syndrome bit evaluation
					if (y > 0) {
						extractor->add_CNOT(2, anc_id, 2, I(n, x+offsetx, y+offsety-1));
					}
					if (y < m) {
						extractor->add_CNOT(2, anc_id, 2, I(n, x+offsetx, y+offsety));
					}
				}
			}
			offset_anc += B;
		}
	}

	//X stabilizer measurements
	for (int offsetx = 0, offset_anc = n2; offsetx < n; offsetx += m) {
		for (int offsety = 0; offsety < n; offsety += m) {
			//vertical
			for (int x = 0; x <= m; x++) {
				for (int y = 0; y < m; y++) {
					int dat_id = I(n, (x+offsetx)%n, y+offsety) + n2;
					int anc_id = I(m, x, y) + offset_anc;
					//printf("[X VERTICAL] data: %d(%d,%d), ancilla: %d(%d,%d)\n", dat_id, (x+offsetx)%n, (y+offsety)%n, anc_id, x, y);
					extractor->add_CNOT(1, anc_id, 0, dat_id);
					extractor->add_2qubit_depo(1, anc_id, 0, dat_id, p);
					if (x > 0) {
						extractor->add_CNOT(1, I(n, x+offsetx-1, y+offsety), 1, anc_id);
					}
					if (x < m) {
						extractor->add_CNOT(1, I(n, x+offsetx, y+offsety), 1, anc_id);
					}
				}
			}
			offset_anc += B;
			//horizontal
			for (int x = 0; x < m; x++) {
				for (int y = 0; y <= m; y++) {
					int dat_id = L(n, x+offsetx, (y+offsety)%n);
					int anc_id = I(m+1, x, y) + offset_anc;
					//printf("[X HORIZONTAL] data: %d(%d,%d), ancilla: %d(%d,%d)\n", dat_id, (x+offsetx)%n, (y+offsety)%n, anc_id, x, y);
					extractor->add_CNOT(1, anc_id, 0, dat_id);
					extractor->add_2qubit_depo(1, anc_id, 0, dat_id, p);
					if (y > 0) {
						extractor->add_CNOT(1, I(n, x+offsetx, y+offsety-1), 1, anc_id);
					}
					if (y < m) {
						extractor->add_CNOT(1, I(n, x+offsetx, y+offsety), 1, anc_id);
					}
				}
			}
			offset_anc += B;
		}
	}

	//measurement error
	
	for (int i = 0; i < n_anc_qubits; i++) {
		extractor->add_1qubit_depo(1, i+n2, p);
		extractor->add_1qubit_depo(2, i+n2, p);
	}
}


/*
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
		extractor->add_1qubit_depo(1, i, p);
		extractor->add_1qubit_depo(2, i, p);
	}
	for (int i = 0; i < n2; i++) {
		extractor->add_phase_flip(1, i, p/10);
		extractor->add_bit_flip(2, i, p/10);
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
		extractor->add_phase_flip(1, i, p);
		extractor->add_bit_flip(2, i, p);
	}
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(0, i, p);
	}
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
	for (int i = 0; i < n2; i++) {
		extractor->add_phase_flip(1, i, p/10);
		extractor->add_bit_flip(2, i, p/10);
	}
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(0, i, p);
	}

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
	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(0, i, p);
	}
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

	for (int i = 0; i < 2*n2; i++) {
		extractor->add_1qubit_depo(1, i, p);
		extractor->add_1qubit_depo(2, i, p);
	}

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
*/
/*
void ToricCode3qubitCat :: build_circuit() {
	assert(n % 2 == 0);
	extractor = new Extractor(2*n2, n2*3/2, n2*3/2);


}
*/
