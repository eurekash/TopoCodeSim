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


