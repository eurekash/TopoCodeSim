#include "extractor.h"
#include "decoder.h"
#include <map>
#include <utility>

typedef std::map< std::pair<char, int>, double > EdgeTable;

class ToricCode {

public:
	ToricCode(int n, double p, int nrounds, int offset_per_round = 0);
	virtual void build_circuit();
	void build_decoder_graph();
	bool simulate();
	
protected:
	void build_decoder_graph_one_round();

	int n, n2;
	int T;
	int offset_per_round;
	double p;
	Extractor *extractor;

	EdgeTable *edges[2];
	Decoder   *decoder[2];
};

/*
class ToricCodeBareAncilla: public ToricCode {
public:
	ToricCodeBareAncilla(int n, double p);
	void build_circuit() override;

};

class ToricCodeBareAncillaNSWE: public ToricCode {
public:
	ToricCodeBareAncillaNSWE(int n, double p);
	void build_circuit() override;

};

class ToricCodeCluster: public ToricCode {
public:
	ToricCodeCluster(int n, double p);
	void build_circuit() override;
};
*/

class ToricCodeBlock: public ToricCode {
public:
	ToricCodeBlock(int k1, int k2, double p);
	void build_circuit() override;
private:
	int k1, k2, m, T;
};
