#ifndef BESTIAU_HPP
#define BESTIAU_HPP
#include <chrono>
#include <deque>
#include <memory>
#include <random>
#include <thread>
#include <vector>

#ifdef DISPLAY
#include <curses.h>

template <typename W> void display(const W &w, bool needClear = false) {
	start_color();
	init_color(COLOR_GREEN, 200, 800, 0);
	init_color(COLOR_RED, 900, 200, 0);
	init_color(COLOR_BLUE, 0, 450, 800);
	double squareV = pow(1.0 + 10.0 * w.bestiau.vit, 2);
	init_color(COLOR_YELLOW, 600 + 100 * squareV, 600 - 100 * squareV,
	           100 + (squareV * 0.2));
	init_pair(1, COLOR_WHITE, COLOR_BLUE);
	init_pair(2, COLOR_GREEN, COLOR_GREEN);
	init_pair(3, COLOR_RED, COLOR_RED);
	init_pair(4, COLOR_YELLOW, COLOR_YELLOW);
	bkgd(COLOR_PAIR(1));

	int row;
	int col;
	getmaxyx(stdscr, row, col);
	double whRatio = w.W / w.H;
	double crRatio = (double)col / (double)row;
	double coef = crRatio < whRatio ? (double)col / w.W : (double)row / w.H;

	if (needClear) {
		clear();
	}

	// bordures
	std::string bordures(w.W * coef, ' ');
	attron(COLOR_PAIR(2));
	for (size_t i = 0; i < w.epaisseurBordures * coef; ++i) {
		mvprintw(i, 0, bordures.c_str());
	}
	for (size_t i = 0; i < w.epaisseurBordures * coef; ++i) {
		mvprintw(w.H * coef - w.epaisseurBordures * coef + i, 0, bordures.c_str());
	}
	attroff(COLOR_PAIR(2));
	std::string clearObstacles(w.largeurObstacles * coef + 5, ' ');
	attron(COLOR_PAIR(1));
	for (size_t r = w.epaisseurBordures * coef + 1;
	     r < (w.H - w.epaisseurBordures) * coef - 1; ++r) {
		mvprintw(r, 0, clearObstacles.c_str());
	}
	std::string strObstacles(w.largeurObstacles * coef, ' ');
	for (auto &o : w.obstacles) {
		if (o.x < w.W - w.largeurObstacles) {
			for (size_t r = w.epaisseurBordures * coef + 1;
			     r < (w.H - w.epaisseurBordures) * coef - 1; ++r) {
				attron(COLOR_PAIR(1));
				mvprintw(r, o.x * coef, clearObstacles.c_str());
				attron(COLOR_PAIR(3));
				if (r < coef * (w.H - (o.y + w.hauteurPassage)) || r > coef * (w.H - o.y))
					mvprintw(r, o.x * coef, strObstacles.c_str());
			}
		}
	}
	attroff(COLOR_PAIR(3));
	// bestiau
	double deltaVit = 0.05;
	attron(COLOR_PAIR(4));
	mvprintw(coef * (w.H - w.bestiau.y), coef * w.bestiau.x, "  ");
	attron(COLOR_PAIR(1));
	mvprintw(0.1 * coef, 2.6 * coef, "%f", w.dist);
	attroff(COLOR_PAIR(1));
	refresh();
}
#endif

typedef std::chrono::duration<int, std::milli> ms;

enum Color { green, red, blue };

struct World {
	struct Obstacle {
		Obstacle(double Y) : y(Y) {}
		Obstacle(double X, double Y) : x(X), y(Y) {}
		Obstacle(Color c, int X, int Y) : color(c), x(X), y(Y) {}
		Color color = blue;
		double x = 0;
		double y = 0;
	};

	struct Bestiau {
		double x = 0.1;
		double y = 0.5;
		double prevy = y;
		double vit = 0;
		double forces = 0;
		double mass = 1.0;
		bool vivant = true;
		void updatePosition(double dt) {
			prevy = y;
			vit += dt * forces / mass;
			y += vit * dt;
		}
	};

	const double H = 1;
	const double W = 5;
	const double epaisseurBordures = 0.06;
	const double largeurObstacles = 0.18;
	const double hauteurPassage = 0.25;
	double dist = 0;
	Bestiau bestiau;
	const double g = -1.5;
	const double force = 15;
	double gravity = 2;
	std::default_random_engine gen{0};

	double dt = 1.0 / 25.0;
	double vitesseDefilement = 0.5;

	double ecartCible = 4;
	double prochainObstacle = 0;

	std::deque<Obstacle> obstacles;
	bool up = false;

	void testCollision() {
		if (bestiau.y <= epaisseurBordures || bestiau.y >= H - epaisseurBordures) {
			bestiau.vivant = false;
		} else if (obstacles.size()) {
			const Obstacle &o = obstacles.front();
			if (o.x <= bestiau.x && o.x + largeurObstacles >= bestiau.x &&
			    (o.y > bestiau.y || o.y + hauteurPassage < bestiau.y)) {
				bestiau.vivant = false;
			}
		}
	}
	void bestiauUp() { up = true; }
	void update() {
		if (up) {
			bestiau.forces = force;
			up = false;
		} else {
			bestiau.forces = -gravity;
		}
		bestiau.updatePosition(dt);
		testCollision();
		double dep = vitesseDefilement * dt;
		dist += dep;
		prochainObstacle -= dep;
		if (prochainObstacle <= 0) {
			std::uniform_real_distribution<double> udist(0, H - hauteurPassage);
			obstacles.push_back(Obstacle(W + 0.1, udist(gen)));
			std::normal_distribution<double> ndist(ecartCible, ecartCible / 4.0);
			prochainObstacle = ndist(gen);
		}
		vitesseDefilement += 0.05 * dt;
		for (auto &o : obstacles) o.x -= dep;
		while (obstacles.size() && obstacles.front().x < (0.0 - largeurObstacles * 2.0))
			obstacles.pop_front();
	}
};

struct FlappyGRNXP {
	static const constexpr int viewDist = 1;
	template <typename G> static G randomInit(size_t nbReguls = 1) {
		G g;
		g.randomParams();
		// input
		// bestiau's height
		// next obstacle's height
		// next obstacle's dist
		g.addRandomProtein(G::ProteinType_t::input, "h");
		g.addRandomProtein(G::ProteinType_t::input, "0h");
		g.addRandomProtein(G::ProteinType_t::input, "0d");
		if (viewDist > 1) {
			g.addRandomProtein(G::ProteinType_t::input, "1h");
			g.addRandomProtein(G::ProteinType_t::input, "1d");
		}
		// two outputs, changing ratio = press up
		g.addRandomProtein(G::ProteinType_t::output, "0");
		g.addRandomProtein(G::ProteinType_t::output, "1");
		g.randomReguls(nbReguls);
		return g;
	}

	template <typename I> static void evaluate(I &ind, bool dbg = false) {
#ifdef DISPLAY
		initscr();
		timeout(0);
		noecho();
		curs_set(FALSE);
#endif
		const int NRUN = 3;
		auto &g = ind.dna;
		double d = 0;
		for (int r = 0; r < NRUN; ++r) {
			World world;
			world.gen.seed(r * 100);
			while (world.bestiau.vivant) {
#ifdef DISPLAY
				std::this_thread::sleep_for(std::chrono::milliseconds(25));
				display(world, world.dist == 0);
#endif
				world.update();
				// update inputs
				g.setInputConcentration("h", world.bestiau.y);
				if (world.obstacles.size() > 0) {
					const auto &o = world.obstacles.front();
					g.setInputConcentration("0h", o.y + world.hauteurPassage * 0.5);
					g.setInputConcentration("0d", o.x / world.W);
					if (viewDist > 1) {
						if (world.obstacles.size() > 1) {
							const auto &o2 = world.obstacles[1];
							g.setInputConcentration("1h", o2.y + world.hauteurPassage * 0.5);
							g.setInputConcentration("1d", o2.x / world.W);
						} else {
							g.setInputConcentration("1h", 0);
							g.setInputConcentration("1d", 1);
						}
					}
				} else {
					g.setInputConcentration("0h", 0);
					g.setInputConcentration("0d", 1);
				}
				bool pbehind = g.getOutputConcentration("0") <= g.getOutputConcentration("1");
				g.step();
				bool behind = g.getOutputConcentration("0") <= g.getOutputConcentration("1");
				if (!behind && pbehind) world.bestiauUp();
			}
			d += world.dist;
		}
		ind.fitnesses["distance"] = d / static_cast<double>(NRUN);
#ifdef DISPLAY
		endwin();
#endif
	}
};

#endif
