#include <set>
#include <utility>
#include <vector>
#include "random.h"

namespace Excitation {
	typedef std::set < std::pair<char, int> >  Set;

	Set XOR(const Set &A, const Set &B);
}

using Excitation :: Set;
using Excitation :: XOR;

struct Wire {
	int x, z;
	Set excite_x;
	Set excite_z;
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
	virtual void apply();
};

struct TwoQubitDepo: Noise {
	Wire *q1, *q2;
	void apply() override;
};

struct OneQubitDepo: Noise {
	Wire *q;
	void apply() override;
};

struct PhaseFlip: Noise {
	Wire *q;
	void apply() override;
};

struct BitFlip: Noise {
	Wire *q;
	void apply() override;
};

class Extractor {

public:
	Extractor(int data_block_size, int anc_x_block_size, int anc_z_block_size);

	void set_syndrome_bit_x(int id_anc_qubit);
	void set_syndrome_bit_z(int id_anc_qubit);
	void add_CNOT(int type_ctrl, int id_ctrl, int type_target, int id_target);
	void execute();
	
	
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

	std::vector<Gate> gates;

	std::vector<int>  syndrome_bit_x;
	std::vector<int>  syndrome_bit_z;

};

