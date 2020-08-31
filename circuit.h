#include <set>
#include <map>
#include <utility>
#include <vector>
#include "random.h"

typedef std::pair<char,int>  pii;

struct Wire {
	int x, z;
	std::vector<pii> excite_x;
	std::vector<pii> excite_z;
	Wire();
};

struct Gate {
	virtual void forward();
	virtual void backward();
	virtual void reset_output();
};

struct CNOT: Gate {
	CNOT(Wire *&ctrl, Wire *&target);
	void forward()  override;
	void backward() override;
	void reset_output()   override;
	Wire *input_ctrl;
	Wire *input_target;
	Wire *output_ctrl;
	Wire *output_target;
};

/*
struct Hadamard: Gate {
	Hadamard(Wire *&target);
	void forward()  override;
	void backward() override;
	void reset_output()   override;
	Wire *input;
	Wire *output;
};
*/

struct Noise {
	Noise();
	std::vector<double> prob;
	int sample();
	int num_errors();
	double probability(int);
	virtual std::vector<pii> excitations(int);
	virtual void apply();
};

struct TwoQubitDepo: Noise {
	Wire *q1, *q2;
	TwoQubitDepo(Wire *q1, Wire *q2, double p);
	std::vector<pii> excitations(int t) override;
	void apply() override;
};

struct OneQubitDepo: Noise {
	Wire *q;
	OneQubitDepo(Wire *q, double p);
	std::vector<pii> excitations(int t) override;
	void apply() override;
};

struct PhaseFlip: Noise {
	Wire *q;
	PhaseFlip(Wire *q, double p);
	std::vector<pii> excitations(int t)  override;
	void apply() override;
};

struct BitFlip: Noise {
	Wire *q;
	BitFlip(Wire *q, double p);
	std::vector<pii> excitations(int t)  override;
	void apply() override;
};

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

	void back_propagation();
	void execute();

	void init_enumerator();
	bool enumerate(std::vector<int> &X0,
				   std::vector<int> &X1,
				   std::vector<int> &Z0,
				   std::vector<int> &Z1,
				   double &p);
	
private:


	int data_block_size;
	int anc_x_block_size;
	int anc_z_block_size;

	Wire **data_start;
	Wire **data_end;
	Wire **anc_x_start;
	Wire **anc_x_end;
	Wire **anc_z_start;
	Wire **anc_z_end;

	std::vector<Gate> gates;

	std::vector<Noise> errors;
	std::vector<Noise>::iterator error_pointer;
	int error_id;

	std::vector<int>  syndrome_bit_x;
	std::vector<int>  syndrome_bit_z;

};

