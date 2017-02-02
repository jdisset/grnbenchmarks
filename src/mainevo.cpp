#define OMP
#include <chrono>
#include "coverage.hpp"
#include "include/gaga/gaga.hpp"
#include "include/grgen/classic.hpp"
#include "include/grgen/grn.hpp"
#include "signalxp.hpp"
#include "bestiau.hpp"

int main(int argc, char** argv) {
	using GRN_t = GRN<Classic>;
	using XP = FlappyGRNXP;
	// using XP = LowPassXP;

	const int popSize = 500;

	// speciation
	for (int c = 0; c < 20; ++c) {
		const int minSpecieSize = 15;
		GAGA::GA<GRN_t> ga(argc, argv);
		ga.setVerbosity(1);
		ga.setPopSaveInterval(0);
		ga.setGenSaveInterval(10);
		ga.setSaveIndStats(true);

		ga.setPopSize(popSize);
		ga.setMutationProba(0.75);
		ga.setCrossoverProba(0.25);
		ga.setNbElites(1);
		ga.setMinSpecieSize(minSpecieSize);
		ga.setSpeciationThreshold(0.3);
		//ga.enableSpeciation();
		ga.setSpeciationThresholdIncrement(0.01);
		ga.setMaxSpeciationThreshold(0.8);
		ga.setMinSpeciationThreshold(0.01);
		ga.setIndDistanceFunction(
		    [](const auto& a, const auto& b) { return GRN_t::getDistance(a.dna, b.dna); });
		ga.setIndDistanceFunction([](const auto& a, const auto& b) {
			double aS = a.dna.getNbProteins();
			double bS = b.dna.getNbProteins();
			return abs(aS - bS) / (aS + bS);
		});
		ga.setEvaluator(
		    [](auto& ind) {
			    XP::evaluate(ind, false);
			    ind.stats["grnSize"] =
			        ind.dna.getFirstOutputIndex() - ind.dna.getFirstRegulIndex();
			  }
		    );

		std::ostringstream fName;
		fName << "../evos/Flappy-P" << ga.getPopSize() << "-M" << ga.getMutationProba()
		      << "-C" << ga.getCrossoverProba();
		if (ga.speciationEnabled()) {
			fName << "-ST" << ga.getSpeciationThreshold() << "-STI"
			      << ga.getSpeciationThresholdIncrement() << "-SMS" << ga.getMinSpecieSize();
		}
		fName << "/";
		ga.setSaveFolder(fName.str());

		std::vector<GRN_t> seeds;
		for (size_t i = 0; i < popSize / (minSpecieSize * 2); ++i) {
			seeds.push_back(XP::randomInit<GRN_t>(1));
		}
		ga.initPopulation([&]() {
			std::uniform_int_distribution<size_t> d(0, seeds.size() - 1);
			auto offspring = seeds[d(ga.globalRand)];
			offspring.mutate();
			return offspring;
		});

		auto t0 = std::chrono::high_resolution_clock::now();
		ga.step(400);
		ga.saveBests(1);
		auto t1 = std::chrono::high_resolution_clock::now();
		double t = std::chrono::duration<double>(t1 - t0).count();
		std::cerr << "TIME = " << t << "s" << std::endl;
	}

	return 0;
}
