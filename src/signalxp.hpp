#ifndef DBLFREQ_HPP
#define DBLFREQ_HPP
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

struct DoubleFreqXP {
	template <typename G> static G randomInit(size_t nbReguls = 1) {
		G g;
		g.randomParams();
		g.addRandomProtein(G::ProteinType_t::input, "i");
		g.addRandomProtein(G::ProteinType_t::output, "o");
		g.randomReguls(nbReguls);
		return g;
	}

	template <typename I> static void evaluate(I &ind, bool dbg = false) {
		ind.fitnesses["fitness"] = dblFreq(ind.dna, 125.0, 1000, dbg) +
		                           dblFreq(ind.dna, 500.0, 1000, dbg) + zero(ind.dna, 1000);
	}

	template <typename Ctrl>
	static double zero(Ctrl &ctrl, int nStepMax, double scaleCoef = 3.0) {
		double fitness = 0.0;
		ctrl.reset();
		ctrl.setInputConcentration("i", 0.0);
		ctrl.step(25);
		for (int nStep = 0; nStep < nStepMax; nStep++) {
			ctrl.setInputConcentration("i", 0.0);
			ctrl.step();
			double dt = 0.0;
			fitness +=
			    std::abs(dt - std::min(1.0, ctrl.getOutputConcentration("o") * scaleCoef)) *
			    (1.0 + dt);
		}
		return -fitness;
	}

	template <typename Ctrl>
	static double dblFreq(Ctrl &ctrl, double halfFreq, int nStepMax, bool dbg,
	                      double scaleCoef = 3.0) {
		double fitness = 0.0;
		ctrl.reset();
		ctrl.setInputConcentration("i", 0.0);
		ctrl.step(25);
		std::array<double, 10> ots{};
		int nbEvents = 0;
		int lastEvent = 0;

		if (dbg) std::cerr << "T;Signal;Desired;Obtained" << std::endl;
		for (int nStep = 0; nStep < nStepMax; nStep++) {
			double sig = std::sin(nStep * M_PI / halfFreq - M_PI * 0.5) * 0.5 + 0.5;
			ctrl.setInputConcentration("i", sig);
			ctrl.step();
			double dt = std::sin(2.0 * nStep * M_PI / halfFreq - M_PI * 0.5) * 0.5 + 0.5;

			for (unsigned int nst = 0; nst < ots.size() - 1; ++nst) ots[nst] = ots[nst + 1];
			double obtained = ctrl.getOutputConcentration("o") * scaleCoef;
			ots[ots.size() - 1] = std::min(1.0, obtained);
			if (nStep - lastEvent > static_cast<int>(ots.size()) &&
			    ((ots[0] < 0.5 && ots[ots.size() - 1] >= 0.5) ||
			     (ots[0] >= 0.5 && ots[ots.size() - 1] < 0.5))) {
				++nbEvents;
				lastEvent = nStep;
			}
			fitness += std::abs(dt - std::min(1.0, ots[ots.size() - 1])) * (1.0 + dt);
			if (dbg)
				std::cerr << nStep << ";" << sig << ";" << dt << ";" << obtained << std::endl;
		}
		double nbEventDesired =
		    static_cast<int>(2.0 * static_cast<double>(nStepMax) / halfFreq);
		if (nbEvents == 0 || nbEvents > nbEventDesired * 2.0) {
			return -fitness;
		} else {
			return -fitness / (2.0 -
			                   static_cast<double>(std::abs(nbEvents - nbEventDesired)) /
			                       static_cast<double>(nbEventDesired));
		}
	}

	static double getY(double t, double p) {
		return 0.5 * std::sin(2.0 * M_PI * t / p - 0.5 * M_PI) + 0.5;
	}
};

struct LowPassXP {
	static const int initSteps = 25;
	static const int CUT_OFF = 50;

	template <typename G> static G randomInit(size_t nbReguls = 0) {
		G g;
		g.randomParams();
		g.addRandomProtein(G::ProteinType_t::input, "i");
		g.addRandomProtein(G::ProteinType_t::output, "o");
		g.randomReguls(nbReguls);
		return g;
	}

	template <typename I> static void evaluate(I &ind, bool dbg = false) {
		const double dt = 0.001;
		ind.fitnesses["fitness"] =
		    lowPass(ind.dna, {{{7.0, 0.7}, {250.0, 0.2}, {1250.0, 0.1}}}, 1.0, dt, dbg) +
		    lowPass(ind.dna, {{{17.0, 0.4}, {350.0, 0.2}, {1100.0, 0.2}, {2000.0, 0.2}}}, 1.0,
		            dt, dbg) +
		    zero(ind.dna, 0.5, dt);
	}

	template <typename G> static double zero(G &g, double T, double dt) {
		double fitness = 0.0;
		g.reset();
		g.setInputConcentration("i", 0.0);
		g.step(initSteps);
		for (double t = 0; t < T; t += dt) {
			g.step();
			fitness += std::pow(g.getOutputConcentration("o"), 2.0);
		}
		return -fitness;
	}

	template <typename G>
	static double lowPass(G &g, const std::vector<std::pair<double, double>> &frequencies,
	                      double T, double dt, bool dbg) {
		if (dbg) std::cerr << "T;Signal;Desired;Obtained" << std::endl;
		double error = 0.0;
		g.reset();
		g.setInputConcentration("i", 0.0);
		g.step(initSteps);
		for (double t = 0; t < T; t += dt) {
			double sig = getWaveValue(t, frequencies);
			g.setInputConcentration("i", sig);
			g.step();
			double desired = 0.0;
			if (frequencies.size() > 0 && frequencies[0].first < static_cast<double>(CUT_OFF))
				desired = getWaveValue(t, {{frequencies[0]}});
			double obtained = g.getOutputConcentration("o");
			error += std::pow(obtained - desired, 2.0);
			if (dbg)
				std::cerr << t << ";" << sig << ";" << desired << ";" << obtained << std::endl;
		}
		return -error;
	}

	static double getWaveValue(double t,
	                           const std::vector<std::pair<double, double>> &frequencies) {
		// frequencies = pair of Freq + Amplitudes
		double res = 0;
		for (auto &f : frequencies) res += f.second * cos(2.0 * M_PI * f.first * t);
		return (res + 1.0) * 0.5;
	}
};
#endif
