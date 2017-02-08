//#define OMP
#define DISPLAY
#include <chrono>
//#include "bestiau.hpp"
#include "coverage.hpp"
#include "coverage.hpp"
#include "include/gaga/gaga.hpp"
#include "include/grgen/int.hpp"
#include "include/grgen/grn.hpp"
#include "ship.hpp"
#include "signalxp.hpp"

int main(int argc, char** argv) {
	using GRN_t = GRN<Classic>;
	// using XP = DoubleFreqXP;
	// using XP = LowPassXP;
	// using XP = FlappyGRNXP;
	// using XP = CoverageXP;
	using XP = ShipEscape::shipXP;

	if (argc > 1) {
		std::ifstream t(argv[1]);
		std::stringstream buffer;
		buffer << t.rdbuf();
		GRN_t g(buffer.str());
		GAGA::Individual<GRN_t> ind(g);
		XP::evaluate(ind, true);
		std::cerr << "Fitnesses : " << std::endl;
		for (auto& f : ind.fitnesses) {
			std::cerr << " - " << f.first << " : " << f.second << std::endl;
		}
	}
	return 0;
}
