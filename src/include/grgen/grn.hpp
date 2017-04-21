#ifndef GENERICGRN_HPP
#define GENERICGRN_HPP
#include <assert.h>
#include <array>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include "common.h"

using std::array;
using std::vector;
using std::map;
using std::string;
using std::pair;
using std::ostringstream;

template <typename Implem> class GRN {
	friend Implem;
	struct GAConfiguration {
		// crossover
		static constexpr double ALIGN_TRESHOLD = 0.15;
		static constexpr double APPEND_NON_ALIGNED = 1.0;
		static constexpr unsigned int MAX_REGULS = 50;

		// mutation
		double MODIF_RATE = 0.25;
		double ADD_RATE = 0.5;
		double DEL_RATE = 0.25;
	};

	friend Implem;

 public:
	using Protein = typename Implem::Protein_t;
	using ProteinType_t = ProteinType;
	using json = nlohmann::json;
	using InfluenceVec = array<double, Implem::nbSignatureParams>;
	template <typename A, typename B> using umap = std::unordered_map<A, B>;
	GAConfiguration config;

 protected:
	array<double, Implem::nbParams> params{};  // alpha, beta, ...
	array<std::unordered_map<string, size_t>, 3> proteinsRefs;
	vector<Protein> actualProteins;
	vector<vector<InfluenceVec>>
	    signatures;  // stores the influence of one protein onto the others (p0
	// influences
	int currentStep = 0;
	Implem implem;
  int func;
  int impl;
  int norm;

 public:
	GRN() { updateSignatures(); }

	GRN(int ff, int ii, int nn) : func(ff), impl(ii), norm(nn) {
    updateSignatures();
  }

	GRN(const GRN& grn)
	    : params(grn.params),
	      proteinsRefs(grn.proteinsRefs),
	      actualProteins(grn.actualProteins),
	      currentStep(grn.currentStep),
        impl(grn.impl),
        func(grn.func),
        norm(grn.norm) {
		updateSignatures();
	}

	/**************************************
	 *          UPDATES
	 *************************************/
	void updateSignatures() {
		// first we order all proteins (inputs | reguls | outputs)
		vector<size_t> orderedID;
		for (auto& t : proteinsRefs) {
			for (auto& p : t) {
				orderedID.push_back(p.second);
				p.second = orderedID.size() - 1;
			}
		}
		assert(orderedID.size() == actualProteins.size());
		auto buffer = actualProteins;
		for (size_t i = 0; i < actualProteins.size(); ++i) {
			actualProteins[i] = buffer[orderedID[i]];
		}
		implem.updateSignatures(*this);
	}

	vector<vector<InfluenceVec>> getSignatures() { return signatures; }

	void step(unsigned int nbSteps = 1) { implem.step(*this, nbSteps); }

	/**************************************
	 *               GET
	 *************************************/
	inline double getProteinConcentration(const string& name, const ProteinType t) const {
		try {
			return actualProteins[proteinsRefs[to_underlying(t)].at(name)].c;
		} catch (...) {
			std::cerr << "Exception raised in getProteinConcentration for name = " << name
			          << ", proteintype = " << to_underlying(t) << std::endl;
			std::cerr << "size = " << actualProteins.size() << std::endl;
			exit(0);
		}
	}
	double getOutputConcentration(const string& name) const {
		return getProteinConcentration(name, ProteinType::output);
	}

	size_t getFirstRegulIndex() const { return getProteinSize(ProteinType::input); }
	size_t getFirstOutputIndex() const {
		return getProteinSize(ProteinType::input) + getProteinSize(ProteinType::regul);
	}

	array<double, Implem::nbParams> getParams() const { return params; }

	inline size_t getProteinSize(ProteinType t) const {
		return proteinsRefs[to_underlying(t)].size();
	}

	size_t getNbProteins() const { return actualProteins.size(); }

	int getCurrentStep() const { return currentStep; }

	Protein& getProtein(ProteinType t, const string& name) {
		size_t id = proteinsRefs[to_underlying(t)].at(name);
		if (id >= actualProteins.size()) throw std::runtime_error("unknown protein");
		return actualProteins[id];
	}

	Protein getProtein_const(ProteinType t, const string& name) const {
		size_t id = proteinsRefs[to_underlying(t)].at(name);
		if (id >= actualProteins.size()) throw std::runtime_error("unknown protein");
		return actualProteins[id];
	}
	vector<Protein> getActualProteinsCopy() const { return actualProteins; }
	vector<Protein>& getActualProteins() { return actualProteins; }
	const vector<Protein>& getActualProteinsConst() const { return actualProteins; }

	/**************************************
	 *               SET
	 *************************************/
	inline void reset() {
		for (auto& p : actualProteins) p.reset();
	}

	void setParam(size_t i, double val) {
		if (i < params.size()) params[i] = val;
	}

	void setProteinConcentration(const string& name, ProteinType t, double c) {
		getProtein(t, name).c = c;
	}

	void setInputConcentration(const string& name, double c) {
		getProtein(ProteinType::input, name).c = c;
	}

	vector<string> getProteinNames(ProteinType t) const {
		vector<string> res;
		for (auto& p : proteinsRefs[(size_t)t]) {
			res.push_back(p.first);
		}
		return res;
	}

	void randomParams() {
		array<pair<double, double>, Implem::nbParams> limits = Implem::paramsLimits();
		for (size_t i = 0; i < Implem::nbParams; ++i) {
			std::uniform_real_distribution<double> distrib(limits[i].first, limits[i].second);
			params[i] = distrib(grnRand);
		}
		updateSignatures();
	}

	/**************************************
	 *          ADDING PROTEINS
	 *************************************/
	void addProtein(const ProteinType t, const string& name, const Protein& p) {
		actualProteins.push_back(p);
		proteinsRefs[to_underlying(t)].insert(make_pair(name, actualProteins.size() - 1u));
		updateSignatures();
	}

	void addRandomProtein(const ProteinType t, const string& name) {
		addProtein(t, name, Protein());
	}

	void addProteins(map<string, Protein>& prots, const ProteinType t) {
		for (auto& p : prots) {
			addProtein(t, p.first, p.second);
		}
	}

	void deleteProtein(size_t id) {
		actualProteins.erase(actualProteins.begin() + static_cast<long>(id));
		// we need to decrement every protein ref after it
		for (auto& t : proteinsRefs) {
			for (auto it = t.begin(); it != t.end();) {
				if ((*it).second == id) {
					it = t.erase(it);
				} else {
					if ((*it).second > id) {
						(*it).second--;
					}
					++it;
				}
			}
		}
	}

	void randomReguls(size_t n) {
		while (proteinsRefs[to_underlying(ProteinType::regul)].size() > 0)
			deleteProtein((proteinsRefs[to_underlying(ProteinType::regul)].begin())->second);

		ostringstream name;
		for (size_t i = 0; i < n; ++i) {
			name.str("");
			name.clear();
			name << "r" << i;
			addProtein(ProteinType::regul, name.str(), Protein());
		}
		updateSignatures();
	}

	void updateRegulNames() {
		int id = 0;
		std::unordered_map<string, size_t> newReguls;
		for (auto& i : proteinsRefs[to_underlying(ProteinType::regul)]) {
			ostringstream name;
			name << "r" << id++;
			newReguls[name.str()] = i.second;
		}
		proteinsRefs[to_underlying(ProteinType::regul)] = newReguls;
	};

	/**************************************
	 *       MUTATION & CROSSOVER
	 *************************************/
	void mutate() {
		std::uniform_real_distribution<double> dReal(0.0, 1.0);
		double dTot = config.MODIF_RATE + config.ADD_RATE + config.DEL_RATE;
		double diceRoll = dReal(grnRand);
		if (diceRoll < config.MODIF_RATE / dTot) {
			// modification (of either a param or a protein)
			double v = 3.0 / static_cast<double>(actualProteins.size() + params.size());
			for (auto& p : actualProteins) {
				if (dReal(grnRand) < v) p.mutate();
			}
			for (size_t paramId = 0; paramId < params.size(); ++paramId) {
				auto limits = Implem::paramsLimits();
				std::uniform_real_distribution<double> distrib(limits[paramId].first,
				                                               limits[paramId].second);
				if (dReal(grnRand) < v) params[paramId] = distrib(grnRand);
			}
		} else if (diceRoll < (config.MODIF_RATE + config.ADD_RATE) / dTot) {
			// we add a new regulatory protein
			ostringstream name;
			name << "r" << getProteinSize(ProteinType::regul);
			addProtein(ProteinType::regul, name.str(), Protein());
		} else {
			// we delete one regulatory protein
			if (getProteinSize(ProteinType::regul) > 0) {
				std::uniform_int_distribution<int> dRegul(getFirstRegulIndex(),
				                                          getFirstOutputIndex() - 1);
				deleteProtein(dRegul(grnRand));
				updateRegulNames();
			}
		}
		updateSignatures();
	}

	GRN crossover(const GRN& other) { return GRN::crossover(*this, other); }

	static double getDistance(const GRN& g0, const GRN& g1) {
		assert(g0.proteinsRefs.size() == g1.proteinsRefs.size());
		assert(g0.params.size() == g1.params.size());
		assert(g0.proteinsRefs[to_underlying(ProteinType::input)].size() ==
		       g1.proteinsRefs[to_underlying(ProteinType::input)].size());
		assert(g0.proteinsRefs[to_underlying(ProteinType::output)].size() ==
		       g1.proteinsRefs[to_underlying(ProteinType::output)].size());

		double pDist = 0, iDist = 0, oDist = 0, rDist = 0;
		// abs(beta1 - beta0)/(bMax - bMin) , same for delta
		for (size_t i = 0; i < g0.params.size(); ++i) {
			pDist += abs(g0.params[i] - g1.params[i]) /
			         (Implem::paramsLimits()[i].second - Implem::paramsLimits()[i].first);
		}

		// inputs
		for (auto& i : g0.proteinsRefs[to_underlying(ProteinType::input)]) {
			iDist += g0.getProtein_const(ProteinType::input, i.first)
			             .getDistanceWith(g1.getProtein_const(ProteinType::input, i.first));
		}
		// outputs
		for (auto& i : g0.proteinsRefs[to_underlying(ProteinType::output)]) {
			oDist += g0.getProtein_const(ProteinType::output, i.first)
			             .getDistanceWith(g1.getProtein_const(ProteinType::output, i.first));
		}

		// faster version:
		std::vector<size_t> r0(g0.getFirstOutputIndex() - g0.getFirstRegulIndex());
		std::iota(std::begin(r0), std::end(r0), g0.getFirstRegulIndex());
		std::vector<size_t> r1(g1.getFirstOutputIndex() - g1.getFirstRegulIndex());
		std::iota(std::begin(r1), std::end(r1), g1.getFirstRegulIndex());
		while (r0.size() > 0 && r1.size() > 0) {
			pair<size_t, size_t> closest;
			double minDist = std::numeric_limits<double>::infinity();
			for (const auto& i : r0) {
				for (const auto& j : r1) {
					double dist = g0.getActualProteinsConst()[i].getDistanceWith(
					    g1.getActualProteinsConst()[j]);
					if (dist < minDist) {
						closest = {i, j};
						minDist = dist;
					}
				}
			}
			assert(minDist < std::numeric_limits<double>::infinity());
			rDist += minDist;
			r0.erase(std::remove(r0.begin(), r0.end(), closest.first), r0.end());
			r1.erase(std::remove(r1.begin(), r1.end(), closest.second), r1.end());
		}

		rDist +=
		    Protein::getMaxDistance() *
		    static_cast<double>(r0.size() + r1.size());  // we add the non aligned distances
		double normalizationCoef =
		    g0.params.size() + g0.proteinsRefs[to_underlying(ProteinType::input)].size() +
		    g0.proteinsRefs[to_underlying(ProteinType::output)].size() +
		    std::max(g0.proteinsRefs[to_underlying(ProteinType::regul)].size(),
		             g1.proteinsRefs[to_underlying(ProteinType::regul)].size());
		return (pDist + iDist + rDist + oDist) / normalizationCoef;
	}

	static GRN crossover(const GRN& g0, const GRN& g1) {
		assert(g0.proteinsRefs.size() == g1.proteinsRefs.size());
		assert(g0.params.size() == g1.params.size());
		assert(g0.proteinsRefs[to_underlying(ProteinType::input)].size() ==
		       g1.proteinsRefs[to_underlying(ProteinType::input)].size());
		assert(g0.proteinsRefs[to_underlying(ProteinType::output)].size() ==
		       g1.proteinsRefs[to_underlying(ProteinType::output)].size());
		GRN offspring;
    offspring.func = g0.func;
    offspring.impl = g0.impl;
    offspring.norm = g0.norm;
		assert(offspring.params.size() == g0.params.size());
		std::uniform_int_distribution<int> d5050(0, 1);
		std::uniform_real_distribution<double> dReal(0.0, 1.0);
		// 50/50 for params, inputs and outputs:

		// params:
		for (size_t i = 0; i < g0.params.size(); ++i) {
			offspring.params[i] = d5050(grnRand) ? g0.params[i] : g1.params[i];
		}

		// inputs
		for (auto& i : g0.proteinsRefs[to_underlying(ProteinType::input)]) {
			if (d5050(grnRand) == 1) {
				offspring.addProtein(ProteinType::input, i.first,
				                     g0.getProtein_const(ProteinType::input, i.first));
			} else {
				offspring.addProtein(ProteinType::input, i.first,
				                     g1.getProtein_const(ProteinType::input, i.first));
			}
		}

		// outputs
		for (auto& i : g0.proteinsRefs[to_underlying(ProteinType::output)]) {
			if (d5050(grnRand) == 1) {
				offspring.addProtein(ProteinType::output, i.first,
				                     g0.getProtein_const(ProteinType::output, i.first));
			} else {
				offspring.addProtein(ProteinType::output, i.first,
				                     g1.getProtein_const(ProteinType::output, i.first));
			}
		}

		// find closest pairs
		auto r0 = g0.getProteinNames(ProteinType::regul);
		auto r1 = g1.getProteinNames(ProteinType::regul);
		vector<pair<Protein, Protein>> aligned;  // first = g0's proteins, second = g1's
		double minDist = 0;
		while (minDist < GAConfiguration::ALIGN_TRESHOLD && r0.size() > 0 && r1.size() > 0 &&
		       aligned.size() < GAConfiguration::MAX_REGULS) {
			pair<string, string> closest;
			minDist = std::numeric_limits<double>::infinity();
			for (const auto& i : r0) {
				for (const auto& j : r1) {
					double dist = g0.getProtein_const(ProteinType::regul, i)
					                  .getDistanceWith(g1.getProtein_const(ProteinType::regul, j));
					if (dist < minDist) {
						closest = {i, j};
						minDist = dist;
					}
				}
			}
			if (minDist < GAConfiguration::ALIGN_TRESHOLD) {
				aligned.push_back({g0.getProtein_const(ProteinType::regul, closest.first),
				                   g1.getProtein_const(ProteinType::regul, closest.second)});
				r0.erase(std::remove(r0.begin(), r0.end(), closest.first), r0.end());
				r1.erase(std::remove(r1.begin(), r1.end(), closest.second), r1.end());
			}
		}
		// ProteinType::regul : 50/50 with aligned
		int id = offspring.getProteinSize(ProteinType::regul);
		assert(id == 0);
		for (auto& i : aligned) {
			ostringstream name;
			name << "r" << id++;
			if (d5050(grnRand))
				offspring.addProtein(ProteinType::regul, name.str(), i.first);
			else
				offspring.addProtein(ProteinType::regul, name.str(), i.second);
		}
		// append the non aligned ones either from r0 or r1
		// if (d5050(grnRand)) {
		auto diceRoll = d5050(grnRand);
		const auto& r = diceRoll ? r0 : r1;
		const auto& g = diceRoll ? g0 : g1;
		for (auto& i : r) {
			if (offspring.getProteinSize(ProteinType::regul) < GAConfiguration::MAX_REGULS) {
				ostringstream name;
				name << "r" << id++;
				offspring.addProtein(ProteinType::regul, name.str(),
				                     g.getProtein_const(ProteinType::regul, i));
			}
		}
		//}
		offspring.updateSignatures();
		return offspring;
	}

	/**************************************
	 *              JSON
	 *************************************/
	GRN(const string& js, int ff, int ii, int nn) {
		auto o = json::parse(js);
		assert(o.count("params"));
		json par = o.at("params");
		assert(par.size() == Implem::nbParams);
		size_t i = 0;
		for (auto& p : par) params[i++] = p.get<double>();
		assert(o.count("proteins"));
		for (size_t t = to_underlying(ProteinType::input);
		     t <= to_underlying(ProteinType::output); ++t) {
			json prots = o.at("proteins").at(typeToString((ProteinType)t));
			for (json::iterator it = prots.begin(); it != prots.end(); ++it) {
				addProtein((ProteinType)t, it.key(), Protein(it.value()));
			}
		}
    func = ff; impl = ii; norm = nn;
		updateSignatures();
	}

	string serialize() const {
		json protObj;
		for (size_t t = 0; t < proteinsRefs.size(); ++t) {
			json pr;
			for (auto p = proteinsRefs[t].begin(); p != proteinsRefs[t].end(); ++p) {
				pr[p->first] = actualProteins[p->second].toJSON();
			}
			protObj[typeToString((ProteinType)t)] = pr;
		}
		json o;
		o["proteins"] = protObj;
		o["params"] = params;
		return o.dump(2);
	}

	std::string typeToString(ProteinType t) const {
		switch (t) {
			case ProteinType::input:
				return "input";
			case ProteinType::regul:
				return "regul";
			case ProteinType::output:
				return "output";
		}
		return "unknown_type";
	}
};
#endif
