#ifndef COVERAGE_HPP
#define COVERAGE_HPP
#include <stdlib.h>
#include <chrono>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <thread>
#include <array>
#include <utility>

#define PURPLE "\033[35m"
#define PURPLEBOLD "\033[1;35m"
#define BLUE "\033[34m"
#define BLUEBOLD "\033[1;34m"
#define GREY "\033[30m"
#define GREYBOLD "\033[1;30m"
#define YELLOW "\033[33m"
#define YELLOWBOLD "\033[1;33m"
#define RED "\033[31m"
#define REDBOLD "\033[1;31m"
#define CYAN "\033[36m"
#define CYANBOLD "\033[1;36m"
#define GREEN "\033[32m"
#define GREENBOLD "\033[1;32m"
#define NORMAL "\033[0m"

using namespace std::chrono_literals;
enum CellType { Visited, NotVisited, Obstacle };
struct Vec2D {
	int x = 0, y = 0;
	Vec2D(int vx, int vy) : x(vx), y(vy) {}
	Vec2D operator+(const Vec2D& v) const { return Vec2D(v.x + x, v.y + y); }
	Vec2D operator*(double s) const { return Vec2D(x * s, y * s); }
};

template <typename M>
Vec2D getNewPos(const M& map, const Vec2D& position, const Vec2D& direction) {
	auto np = position + direction;
	if (np.x < 0)
		np.x += map.size();
	else if (np.x >= static_cast<int>(map.size()))
		np.x -= map.size();
	if (np.y < 0)
		np.y += map[0].size();
	else if (np.y >= static_cast<int>(map[0].size()))
		np.y -= map[0].size();
	return np;
}

template <typename M> struct Agent {
	M& map;
	Vec2D position{0, 0};
	Vec2D direction{1, 0};
	Agent(M& m) : map(m) { map[position.x][position.y] = CellType::Visited; }
	void go() {
		auto np = getNewPos(map, position, direction);
		if (map[np.x][np.y] != CellType::Obstacle) {
			map[np.x][np.y] = CellType::Visited;
			position = np;
		}
	}
};

template <typename M> std::tuple<int, int, int> countMap(const M& m) {
	int nbVisited = 0;
	int nbNotVisited = 0;
	int nbObstacles = 0;
	for (const auto& i : m) {
		for (const auto& j : i) {
			switch (j) {
				case CellType::Visited:
					++nbVisited;
					break;
				case CellType::NotVisited:
					++nbNotVisited;
					break;
				default:
					++nbObstacles;
					break;
			}
		}
	}
	return (std::tuple<int, int, int>){nbVisited, nbNotVisited, nbObstacles};
}

template <typename A> void printMap(const A& robot) {
	std::ostringstream o;
	int X = robot.map.size();
	int Y = robot.map[0].size();

	int nbVisited = 1;  // robot.position
	int nbNotVisited = 0;
	int nbObstacles = 0;

	for (int y = 0; y < Y; ++y) {
		for (int x = 0; x < X; ++x) {
			if (robot.position.x == x && robot.position.y == y)
				o << YELLOWBOLD << "♖ " << NORMAL;
			else {
				const auto& c = robot.map[x][y];
				if (c == CellType::Visited) {
					o << GREEN << "░ " << NORMAL;
					++nbVisited;
				} else if (c == CellType::NotVisited) {
					o << PURPLEBOLD << "░ " << NORMAL;
					++nbNotVisited;
				} else {
					o << RED << "█ " << NORMAL;
					++nbObstacles;
				}
			}
		}
		o << std::endl;
	}
	std::cerr << o.str() << std::endl
	          << std::endl
	          << nbVisited << " visited; " << nbNotVisited << " not visited; "
	          << nbObstacles << " obstacles" << std::endl;
}

template <size_t N, class T> std::array<T, N> make_array(const T& v) {
	std::array<T, N> ret;
	ret.fill(v);
	return ret;
}

template <size_t X, size_t Y> auto generateRandomMap(size_t nbObstacles, int seed) {
	auto m = make_array<X>(make_array<Y>(CellType::NotVisited));
	std::mt19937 gen(seed);
	std::uniform_int_distribution<> dX(0, X);
	std::uniform_int_distribution<> dY(0, Y);
	for (size_t i = 0; i < nbObstacles; ++i) m[dX(gen)][dY(gen)] = CellType::Obstacle;
	return m;
}

struct CoverageXP {
  template <typename G> static G randomInit(size_t nbReguls = 0, int ff = 0, int ii = 0, int nn = 0) {
		G g(ff, ii, nn);
		g.randomParams();
		// 4 inputs :  sumed reward = nbNotVisited in each direction in the next 3 cells, /3
		g.addRandomProtein(G::ProteinType_t::input, "nN");
		g.addRandomProtein(G::ProteinType_t::input, "nW");
		g.addRandomProtein(G::ProteinType_t::input, "nS");
		g.addRandomProtein(G::ProteinType_t::input, "nE");
		// 4 inputs : nb of obstacles in the 3 adjacent cells, each direction /3
		g.addRandomProtein(G::ProteinType_t::input, "oN");
		g.addRandomProtein(G::ProteinType_t::input, "oW");
		g.addRandomProtein(G::ProteinType_t::input, "oS");
		g.addRandomProtein(G::ProteinType_t::input, "oE");
		// 1 input : reward = nbVisited / nbCells
		g.addRandomProtein(G::ProteinType_t::input, "r");
		// 4 outputs : move N,W,S,E
		g.addRandomProtein(G::ProteinType_t::output, "N");
		g.addRandomProtein(G::ProteinType_t::output, "W");
		g.addRandomProtein(G::ProteinType_t::output, "S");
		g.addRandomProtein(G::ProteinType_t::output, "E");
		g.randomReguls(nbReguls);
		return g;
	}

	template <typename M>
	static std::pair<double, double> look(const M& m, const Vec2D& pos, const Vec2D& dir) {
		const double VIEWDIST = 3;
		double nbNotVisited = 0, nbObstacles = 0;
		for (double j = 0; j < VIEWDIST; ++j) {
			auto np = getNewPos(m, pos, dir * j);
			if (m[np.x][np.y] == CellType::NotVisited)
				++nbNotVisited;
			else if (m[np.x][np.y] == CellType::Obstacle)
				++nbObstacles;
		}
		nbNotVisited /= VIEWDIST;
		nbObstacles /= VIEWDIST;
		return {nbNotVisited, nbObstacles};
	}

	template <typename I> static void evaluate(I& ind, bool dbg = false) {
		const size_t X = 10;
		const size_t Y = 10;
		const int nbSteps = X * Y * 2;
		const int NMAP = 3;

		double fit = 0;
		auto& g = ind.dna;
		for (int s = 0; s < NMAP; ++s) {
			auto m = generateRandomMap<X, Y>(20, s * 100);
			Agent<decltype(m)> robot(m);
			for (int i = 0; i < nbSteps; ++i) {
				{
					auto l = look(robot.map, robot.position, Vec2D(0, 1));
					g.setInputConcentration("nN", l.first);
					g.setInputConcentration("oN", l.second);
				}
				{
					auto l = look(robot.map, robot.position, Vec2D(0, -1));
					g.setInputConcentration("nS", l.first);
					g.setInputConcentration("oS", l.second);
				}
				{
					auto l = look(robot.map, robot.position, Vec2D(1, 0));
					g.setInputConcentration("nE", l.first);
					g.setInputConcentration("oE", l.second);
				}
				{
					auto l = look(robot.map, robot.position, Vec2D(-1, 0));
					g.setInputConcentration("nW", l.first);
					g.setInputConcentration("oW", l.second);
				}
				auto c = countMap(robot.map);
				g.setInputConcentration(
				    "r", static_cast<double>(std::get<0>(c)) / static_cast<double>(X * Y));
				g.step();
				char maxProt = 'N';
				double maxConcentration = g.getOutputConcentration("N");
				{
					char prot = 'S';
					double con = g.getOutputConcentration(std::string(1, prot));
					if (con > maxConcentration) {
						maxConcentration = con;
						maxProt = prot;
					}
				}
				{
					char prot = 'W';
					double con = g.getOutputConcentration(std::string(1, prot));
					if (con > maxConcentration) {
						maxConcentration = con;
						maxProt = prot;
					}
				}
				{
					char prot = 'E';
					double con = g.getOutputConcentration(std::string(1, prot));
					if (con > maxConcentration) {
						maxConcentration = con;
						maxProt = prot;
					}
				}
				switch (maxProt) {
					case 'N':
						robot.direction = Vec2D(0, 1);
						robot.go();
						break;
					case 'S':
						robot.direction = Vec2D(0, -1);
						robot.go();
						break;
					case 'W':
						robot.direction = Vec2D(-1, 0);
						robot.go();
						break;
					case 'E':
						robot.direction = Vec2D(1, 0);
						robot.go();
						break;
				}
				if (dbg) {
					system("clear");
					printMap(robot);
					std::this_thread::sleep_for(0.05s);
				}
			}
			auto c = countMap(robot.map);
			fit += std::get<0>(c) / static_cast<double>(X * Y * NMAP);
		}
		ind.fitnesses["exploration"] = fit;
	}
};
#endif
