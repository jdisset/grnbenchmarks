//#define OMP
#include <string>
#include <chrono>
#include "include/gaga/gaga.hpp"
#include "include/grgen/classic.hpp"
#include "include/grgen/real.hpp"
#include "include/grgen/grn.hpp"
#include "signalxp.hpp"
#include "coverage.hpp"
#include "bestiau.hpp"
// #include "ship.hpp"

// using XP = DoubleFreqXP;
// using XP = LowPassXP;
using XP = FlappyGRNXP;
// using XP = CoverageXP;
// using XP = ShipEscape::shipXP;

template <typename GRN_t>
void runGA(int argc, char** argv, string exp_name, int ff, int ii, int nn){
  const int minSpecieSize = 15;
	const int popSize = 500;
  GAGA::GA<GRN_t> ga(argc, argv);
  ga.setVerbosity(0);
  ga.setPopSaveInterval(0);
  ga.setGenSaveInterval(0);
  ga.setSaveIndStats(true);

  ga.setPopSize(popSize);
  ga.setMutationProba(0.75);
  ga.setCrossoverProba(0.25);
  ga.setNbElites(1);
  ga.setMinSpecieSize(minSpecieSize);
  ga.setSpeciationThreshold(0.3);
  ga.enableSpeciation();
  ga.setSpeciationThresholdIncrement(0.01);
  ga.setMaxSpeciationThreshold(0.8);
  ga.setMinSpeciationThreshold(0.01);

  ga.setIndDistanceFunction(
      [](const auto& a, const auto& b) { return GRN_t::getDistance(a.dna, b.dna); });
  ga.setEvaluator(
      [](auto& i) {
        XP::evaluate(i, false);
        i.stats["grnSize"] = i.dna.getNbProteins();
      }, exp_name);

  std::vector<GRN_t> seeds;
  for (size_t i = 0; i < popSize /(minSpecieSize * 2); ++i) {
    seeds.push_back(XP::randomInit<GRN_t>(1, ff, ii, nn));
  }

  std::ostringstream fName;
  fName << "evos/";
  ga.setSaveFolder(fName.str());

  ga.initPopulation([&]() {
    std::uniform_int_distribution<size_t> d(0, seeds.size() - 1);
    auto offspring = seeds[d(ga.globalRand)];
    offspring.mutate();
    return offspring;
  });

  ga.step(300);
  ga.saveBests(1);
}

int main(int argc, char** argv) {

  runGA<GRN<RealC>>(argc, argv, "s1_a1_f0_n0", 0, 1, 0);

	return 0;
}
