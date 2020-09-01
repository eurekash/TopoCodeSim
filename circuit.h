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


