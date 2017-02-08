#ifndef REALC_HPP
#define REALC_HPP
#include <iostream>
#include <array>
#include <vector>
#include <utility>
#include <unordered_map>
#include "common.h"
#include "protein.hpp"

using namespace std;

struct RealC {
	// we use 3 coordinates proteins (id, enh, inh)
	using Protein_t = Protein<3, double, 0, 1>;

	// we need 2 parameters (beta, alpha)
	static constexpr unsigned int nbParams = 2;
	// and we produce 2 dimensional signatures (enhnance, inhibit)
	static constexpr unsigned int nbSignatureParams = 2;

	static const array<pair<double, double>, nbParams> paramsLimits() {
		return {{{0.5, 2.0}, {0.5, 2.0}}};
	}

	// helpers for proteins coords
	static inline double& getId(Protein_t& p) { return p.coords[0]; }
	static inline double& getEnh(Protein_t& p) { return p.coords[1]; }
	static inline double& getInh(Protein_t& p) { return p.coords[2]; }

	// aliases for ProteinType
	static constexpr ProteinType pinput = ProteinType::input;
	static constexpr ProteinType pregul = ProteinType::regul;
	static constexpr ProteinType poutput = ProteinType::output;

	double maxEnhance = 0.0, maxInhibit = 0.0;

  RealC() {}

	template <typename GRN> void updateSignatures(GRN& grn) {
		grn.signatures.clear();
		grn.signatures.resize(grn.actualProteins.size());
		for (size_t i = 0; i < grn.actualProteins.size(); ++i) {
			grn.signatures[i].resize(grn.actualProteins.size());
			for (size_t j = 0; j < grn.actualProteins.size(); ++j) {
				auto& p0 = grn.actualProteins[i];
				auto& p1 = grn.actualProteins[j];
				grn.signatures[i][j] = {
				    {static_cast<double>(abs(getEnh(p0) - getId(p1))),
				     static_cast<double>(abs(getInh(p0) - getId(p1)))}};
				if (grn.signatures[i][j][0] > maxEnhance) maxEnhance = grn.signatures[i][j][0];
				if (grn.signatures[i][j][1] > maxInhibit) maxInhibit = grn.signatures[i][j][1];
			}
		}
		// std::cerr << "maxEnh = " << maxEnhance << ", maxInh = " << maxInhibit << std::endl;
		for (size_t i = 0; i < grn.actualProteins.size(); ++i) {
			for (size_t j = 0; j < grn.actualProteins.size(); ++j) {
        double inside_enh = 0;
        double inside_inh = 0;
        switch(grn.impl) {
        case 0:
				  inside_enh = -grn.params[0]*(grn.signatures[i][j][0]);
				  inside_inh = -grn.params[0]*(grn.signatures[i][j][1]);
          break;
        case 1:
				  inside_enh = grn.params[0]*(1.0-grn.signatures[i][j][0])-maxEnhance;
				  inside_inh = grn.params[0]*(1.0-grn.signatures[i][j][1])-maxInhibit;
          break;
        case 2:
				  inside_enh = grn.params[0]*(1.0-grn.signatures[i][j][0]-maxEnhance);
				  inside_inh = grn.params[0]*(1.0-grn.signatures[i][j][1]-maxInhibit);
          break;
        }
        switch(grn.func) {
        case 0:
          grn.signatures[i][j] = {{exp(inside_enh), exp(inside_inh)}};
          break;
        case 1:
          grn.signatures[i][j] = {{tanh(inside_enh)+1, tanh(inside_inh)+1}};
          break;
        case 2:
          grn.signatures[i][j] = {{2/(1+exp(-inside_enh)), 2/(1+exp(-inside_inh))}};
          break;
        }
			}
		}
	}

	template <typename GRN> void step(GRN& grn, unsigned int nbSteps) {
		for (auto s = 0u; s < nbSteps; ++s) {
			std::vector<double> nextProteins;  // only reguls & outputs concentrations
			nextProteins.reserve(grn.getNbProteins() - grn.getProteinSize(ProteinType::input));
			const auto firstOutputId = grn.getFirstOutputIndex();
			for (size_t j = grn.getFirstRegulIndex(); j < grn.getNbProteins(); ++j) {
				double enh = 0.0, inh = 0.0;
				for (size_t k = 0; k < firstOutputId; ++k) {
					enh += grn.actualProteins[k].c * grn.signatures[k][j][0];
					inh += grn.actualProteins[k].c * grn.signatures[k][j][1];
				}
				nextProteins.push_back(
				    max(0.0, grn.actualProteins[j].c +
				                 (grn.params[1] / static_cast<double>(grn.getNbProteins())) *
				                     (enh - inh)));
			}
			// Normalizing regul & output proteins concentrations
      if (grn.norm == 0) {
        double sumConcentration = 0.0;
        for (auto i : nextProteins) {
          sumConcentration += i;
        }
        if (sumConcentration > 0) {
          for (auto& i : nextProteins) {
            i /= sumConcentration;
          }
        }
      } else { // capping
        for (auto& i : nextProteins) {
          i = min(1.0, i);
        }
      }
			auto firstRegulIndex = grn.getFirstRegulIndex();
			for (size_t i = firstRegulIndex; i < grn.getNbProteins(); ++i) {
				grn.actualProteins[i].c = nextProteins[i - firstRegulIndex];
			}
		}
	}
};
#endif
