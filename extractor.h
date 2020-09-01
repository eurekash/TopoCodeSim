#include "circuit.h"

class Extractor {

public:
	Extractor(int data_block_size, int anc_x_block_size, int anc_z_block_size);

	void set_syndrome_bit_x(int id_anc_qubit);
	void set_syndrome_bit_z(int id_anc_qubit);
	void add_CNOT(int type_ctrl, int id_ctrl, int type_target, int id_target);
	void add_2qubit_depo(int type1, int id1, int type2, int id2, double p);
	void add_1qubit_depo(int type, int id, double p);
	void add_bit_flip(int type, int id, double p);
	void add_phase_flip(int type, int id, double p);

	void execute();

	void init_enumerator();
	bool enumerate(std::vector<int> &X0,
				   std::vector<int> &X1,
				   std::vector<int> &Z0,
				   std::vector<int> &Z1,
				   double &p);
	
private:

	void back_propagation();

	int data_block_size;
	int anc_x_block_size;
	int anc_z_block_size;

	Wire **data_start;
	Wire **data_end;
	Wire **anc_x_start;
	Wire **anc_x_end;
	Wire **anc_z_start;
	Wire **anc_z_end;

	std::vector<Gate*> gates;

	std::vector<Noise*> errors;
	int err_id, err_choice;

	std::vector<int>  syndrome_bit_x;
	std::vector<int>  syndrome_bit_z;

};
