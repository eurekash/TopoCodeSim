#include "circuit.h"
#include <cassert>
#include <algorithm>
#include <iterator>

template <typename T>
T XOR(const T &A, const T &B) {
	T result;
	std::set_symmetric_difference(
			A.begin(), A.end(),
			B.begin(), B.end(),
			std::back_inserter(result));
	return result;
}

Wire :: Wire() {
	x = z = 0;
	excite_x.clear();
	excite_z.clear();
}

void Gate :: forward() {}
void Gate :: backward() {}
void Gate :: reset_output() {}

Noise :: Noise(): prob(0)  {}

int Noise :: num_errors() {
	return prob.size();
}

std::vector<pii> Noise :: excitations(int)  {
	return std::vector<pii>();
}

int Noise :: sample() {
	double p = Random::get_instance() -> random_number();
	double sum = 0;
	int nchoices = prob.size();
	for (int i = 0; i < nchoices; i++) {
		sum += prob[i];
		if (p <= sum)  return i;
	}
	return nchoices;
}

void Noise :: apply() {}

double Noise :: probability(int t) {
	assert(t >= 0 && t < prob.size());
	return prob[t];
}

TwoQubitDepo :: TwoQubitDepo(Wire *q1, Wire *q2, double p) {
	this->q1 = q1;
	this->q2 = q2;
	for (int i = 0; i < 15; i++) {
		prob.emplace_back(p / 15);
	}
}


std::vector<pii> TwoQubitDepo :: excitations(int t) 
{
	assert(t < 15);
	std::vector<pii> result;
	switch (t / 4) {
		case 0:
			result = q1->excite_x;
			break;
		case 1:
			result = XOR(q1->excite_x, q1->excite_z);
			break;
		case 2:
			result = q1->excite_z;
			break;
		default:
			break;
	}

	switch (t % 4) {
		case 0:
			result = XOR(result, q2->excite_x);
			break;
		case 1:
			result = XOR(result, XOR(q2->excite_x, q2->excite_z));
			break;
		case 2:
			result = XOR(result, q2->excite_z);
			break;
		default:
			break;
	}

	return result;
}

void TwoQubitDepo :: apply() {
	int t = sample();
	switch (t / 4) {
		case 0: 
			q1->x ^= 1;
			break;
		case 1:
			q1->x ^= 1, q1->z ^= 1;
			break;
		case 2:
			q1->z ^= 1; 
			break;
		default:
			break;
	}

	switch (t % 4) {
		case 0: 
			q2->x ^= 1;
			break;
		case 1:
			q2->x ^= 1, q2->z ^= 1;
			break;
		case 2:
			q2->z ^= 1;
			break;
		default:
			break;
	}
}


OneQubitDepo :: OneQubitDepo(Wire *q, double p) {
	this->q = q;
	for (int i = 0; i < 3; i++) {
		prob.emplace_back(p / 3);
	}
}

void OneQubitDepo :: apply() {
	int t = sample();
	switch (t) {
		case 0:
			q->x ^= 1;
			break;
		case 1:
			q->x ^= 1;
			q->z ^= 1;
			break;
		case 2:
			q->z ^= 1;
			break;
		default:
			break;
	}
}

std::vector<pii> OneQubitDepo :: excitations(int t) {
	assert(t < 3);
	std::vector<pii> result;
	switch (t) {
		case 0:
			result = q->excite_x;
			break;
		case 1:
			result = XOR(q->excite_x, q->excite_z);
			break;
		case 2:
			result = q->excite_z;
			break;
		default:
			break;
	}
	return result;
}

PhaseFlip :: PhaseFlip(Wire *q, double p) 
{
	this->q = q;
	prob.emplace_back(p);
}

void PhaseFlip :: apply() {
	if (sample() == 0) {
		q->z ^= 1;
	}
}

std::vector<pii> PhaseFlip :: excitations(int t) {
	assert(t == 0);
	return q->excite_z;
}

BitFlip :: BitFlip(Wire *q, double p) {
	this->q = q;
	prob.emplace_back(p);
}

std::vector<pii> BitFlip :: excitations(int t) {
	assert(t == 0);
	return q->excite_x;
}


void BitFlip :: apply() {
	if (sample() == 0) {
		q->x ^= 1;
	}
}

CNOT :: CNOT(Wire *&ctrl, Wire *&target) {
	input_ctrl = ctrl;
	input_target = target;
	ctrl = output_ctrl = new Wire();
	target = output_target = new Wire();
}

void CNOT :: forward() {
	output_ctrl->x ^= input_ctrl->x;
	output_ctrl->z ^= input_ctrl->z ^ input_target->z;
	output_target->x ^= input_ctrl->x ^ input_target->x;
	output_target->z ^= input_target->z;
}

void CNOT :: reset_output() {
	output_ctrl->x = 0;
	output_ctrl->z = 0;
	output_target->x = 0;
	output_target->z = 0;
}

void CNOT :: backward() {
	input_ctrl->excite_x = XOR(output_ctrl->excite_x, output_target->excite_x);
	input_ctrl->excite_z = output_ctrl->excite_z;
	input_target->excite_x = output_target->excite_x;
	input_target->excite_z = XOR(output_ctrl->excite_z, output_target->excite_z);
}

Extractor :: Extractor(int data_block_size, int anc_x_block_size, int anc_z_block_size):
gates(0), syndrome_bit_x(0), syndrome_bit_z(0), errors(0)
{
	this->data_block_size = data_block_size;
	this->anc_x_block_size = anc_x_block_size;
	this->anc_z_block_size = anc_z_block_size;

	data_start  = new Wire*[data_block_size];
	data_end    = new Wire*[data_block_size];
	for (int i = 0; i < data_block_size; i++) {
		data_start[i] = data_end[i] = new Wire();
	}

	anc_x_start = new Wire*[anc_x_block_size];
	anc_x_end   = new Wire*[anc_x_block_size];
	for (int i = 0; i < anc_x_block_size; i++) {
		anc_x_start[i] = anc_x_end[i] = new Wire();
	}

	anc_z_start = new Wire*[anc_z_block_size];
	anc_z_end   = new Wire*[anc_z_block_size];
	for (int i = 0; i < anc_z_block_size; i++) {
		anc_z_start[i] = anc_z_end[i] = new Wire();
	}

}

void Extractor :: add_CNOT(int type_ctrl, int id_ctrl, int type_target, int id_target) {
	Wire *&ctrl = (type_ctrl == 0)? data_end[id_ctrl] : 
		          (type_ctrl == 1)? anc_x_end[id_ctrl]: 
				  					anc_z_end[id_ctrl]; 
	Wire *&target = (type_target == 0)? data_end[id_target]: 
					(type_target == 1)? anc_x_end[id_target]: 
										anc_z_end[id_ctrl];

	gates.emplace_back(CNOT(ctrl, target));
}

void Extractor :: add_2qubit_depo(int type1, int id1, int type2, int id2, double p) {
	Wire *q1 = (type1 == 0)? data_end[id1]:
			   (type1 == 1)? anc_x_end[id1]:
			   				 anc_z_end[id1];

	Wire *q2 = (type2 == 0)? data_end[id2]:
			   (type2 == 1)? anc_x_end[id2]:
			   				 anc_z_end[id2];

	errors.emplace_back(TwoQubitDepo(q1, q2, p));
}

void Extractor :: add_1qubit_depo(int type, int id, double p) {
	Wire *q = (type == 0)? data_end[id]:
			  (type == 1)? anc_x_end[id]:
			   			   anc_z_end[id];
	errors.emplace_back(OneQubitDepo(q, p));
}

void Extractor :: add_phase_flip(int type, int id, double p) {
	Wire *q = (type == 0)? data_end[id]:
			  (type == 1)? anc_x_end[id]:
			   			   anc_z_end[id];

	errors.emplace_back(PhaseFlip(q, p));
}


void Extractor :: add_bit_flip(int type, int id, double p) {
	Wire *q = (type == 0)? data_end[id]:
			  (type == 1)? anc_x_end[id]:
			   			   anc_z_end[id];

	errors.emplace_back(BitFlip(q, p));
}


void Extractor :: set_syndrome_bit_x(int id_anc_qubit) {
	syndrome_bit_x.emplace_back(id_anc_qubit);
}

void Extractor :: set_syndrome_bit_z(int id_anc_qubit) {
	syndrome_bit_z.emplace_back(id_anc_qubit);
}

void Extractor :: back_propagation() {
	//Step 1
	for (int i = 0; i < syndrome_bit_x.size(); i++) {
		anc_x_end[syndrome_bit_x[i]]->excite_z.emplace_back(std::make_pair(0x2, i));  //bitmask = 10
	}
	for (int i = 0; i < syndrome_bit_z.size(); i++) {
		anc_z_end[syndrome_bit_z[i]]->excite_x.emplace_back(std::make_pair(0x3, i));  //bitmask = 11
	}

	for (std::vector<Gate>::reverse_iterator g = gates.rbegin(); g != gates.rend(); ++g) {
		g->backward();
	}

	//Step 0
	//
	for (int i = 0; i < data_block_size; i++) {
		data_end[i]->excite_x = data_start[i]->excite_x;
		data_end[i]->excite_z = data_start[i]->excite_z;
	}

	for (int i = 0; i < anc_x_block_size; i++) {
		anc_x_end[i]->excite_x.clear();
		anc_x_end[i]->excite_z.clear();
	}

	for (int i = 0; i < anc_z_block_size; i++) {
		anc_z_end[i]->excite_x.clear();
		anc_z_end[i]->excite_z.clear();
	}

	for (int i = 0; i < syndrome_bit_x.size(); i++) {
		anc_x_end[syndrome_bit_x[i]]->excite_z.emplace_back(std::make_pair(0x0, i));  //bitmask = 00
	}
	for (int i = 0; i < syndrome_bit_z.size(); i++) {
		anc_z_end[syndrome_bit_z[i]]->excite_x.emplace_back(std::make_pair(0x1, i));  //bitmask = 01
	}

	for (std::vector<Gate>::reverse_iterator g = gates.rbegin(); g != gates.rend(); ++g) {
		g->backward();
	}
}

void Extractor :: execute() {
	//clear all the Wires 
	for (auto &g: gates)  g.reset_output();

}

void Extractor :: init_enumerator() {
	error_pointer = errors.begin();
	error_id = 0;
}

bool Extractor::enumerate(
		std::vector<int> &X0,
		std::vector<int> &X1,
		std::vector<int> &Z0,
		std::vector<int> &Z1,
		double &p )
{
	if (error_pointer == errors.end())  return false;

	p = error_pointer->probability(error_id);
	std::vector<pii> S = error_pointer->excitations(error_id);
	

	X0.clear(), X1.clear();
	Z0.clear(), Z1.clear();

	for (auto &e: S) {
		switch (e.first) {
			case 0x0: 
				X0.emplace_back(e.second);
				break;
			case 0x1:
				Z0.emplace_back(e.second);
				break;
			case 0x2:
				X1.emplace_back(e.second);
				break;
			case 0x3:
				Z1.emplace_back(e.second);
				break;
		}
	}

	X1 = XOR(X0, X1);
	Z1 = XOR(Z0, Z1);

	if (++error_id == error_pointer->num_errors()) {
		++error_pointer;
		error_id = 0;
	}
	return true;
}

