#include "extractor.h"
#include "decoder.h"
#include <map>
#include <utility>

typedef std::map< std::pair<char, int>, double > EdgeTable;


class ToricAncBlk {

public:
	ToricAncBlk(int m, int n, int T, double p);
	void build_circuit();
	void build_decoder_graph();
	bool simulate();
	
protected:
	void build_decoder_graph_one_round();

	int m;
	int n;
	int T;
	double p;
	Extractor *extractor_x;
	Extractor *extractor_z;

	int num_hori;
	int num_vert;
	int num_xstabs;
	int num_zstabs;
	EdgeTable *edges;
	Decoder   *decoder;
};


