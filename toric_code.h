#include "circuit.h"
#include <map>
#include <utility>

typedef std::map< std::pair<char, int>, double > EdgeTable;
class ToricCode {

public:
	ToricCode(int n, double p);
	virtual void build_circuit();
	void build_decoder_one_round();

protected:
	int n, n2;
	double p;
	Extractor *extractor;
	EdgeTable *edges_x;
	EdgeTable *edges_z;
};

class ToricCodeBareAncilla: public ToricCode {
public:
	ToricCodeBareAncilla(int n, double p);
	void build_circuit() override;

};



