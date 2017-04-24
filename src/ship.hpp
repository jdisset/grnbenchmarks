#ifndef SHIP_HPP
#define SHIP_HPP
#include <cmath>
#include <iostream>
#include <random>
#include <unordered_map>

#ifdef SHIP_DISPLAY
#include <QtCore/qmath.h>
#include <QtGui/QGuiApplication>
#include <QtGui/QMatrix4x4>
#include <QtGui/QOpenGLShaderProgram>
#include <QtGui/QScreen>
#include "viewer/shipwindow.hpp"
#endif

using namespace std;

#define CONTROLE 0.1
#define FRICTION 0.1

namespace ShipEscape {
struct V {
	V(){};
	V(double X, double Y) : x(X), y(Y) {}
	double x = 0;
	double y = 0;
	void rotate(double angle) {
		double cs = cos(angle);
		double sn = sin(angle);
		double px = x * cs - y * sn;
		double py = x * sn + y * cs;
		x = px;
		y = py;
	}
	void normalize() {
		double l = sqrt(x * x + y * y);
		if (l > 0) {
			x /= l;
			y /= l;
		}
	}
	double sqLength() const { return x * x + y * y; }
	V operator-(const V &v) const { return V(x - v.x, y - v.y); }
	V operator+(const V &v) const { return V(x + v.x, y + v.y); }
	V &operator+=(const V &v) {
		x += v.x;
		y += v.y;
		return *this;
	}
	double dot(const V &v) const { return x * v.x + y * v.y; }
};
V operator*(const V &v, const double &d) { return V(v.x * d, v.y * d); }
V operator/(const V &v, const double &d) { return V(v.x / d, v.y / d); }

struct Ship {
	double rotationSpeed = 1.0;
	V dimensions = V(1.3, 2.0);
	V position;
	V velocity;
	V forces;
	V orientation = V(0, 1);
	double thrustPower = 2500.0;
	bool thrusting = false;  // only useful for viewer
	Ship(){};
	void rotate(double dir, double dt) {
		dir = max(dir, -1.0);
		dir = min(dir, 1.0);
		orientation.rotate(dir * rotationSpeed * dt);
		orientation.normalize();
	}
	double getAngle() { return atan2(orientation.y, orientation.x) - atan2(1.0, 0.0); }
	void updatePosition(double dt) {
		double vit = sqrt(velocity.sqLength());
		V nV = velocity / vit;
		double r = max(0.0, nV.dot(orientation));
		if (vit > 0) {
			velocity = orientation * vit * r * CONTROLE + (velocity) * (1.0 - CONTROLE * r);
		}
		velocity = velocity + (forces - velocity * FRICTION) * dt;
		position += velocity * dt;
	}
	void thrust(double dt) {
		double actualQ = thrustPower * dt;
		forces = forces + orientation * actualQ;
		thrusting = true;
	}
};

struct Circle {
	Circle(const V &c, const double &r) : center(c), radius(r) {}
	V center;
	double radius = 0.2;
};

struct World {
	typedef V Vv;
	vector<Ship> ships;
	double gridSize = 10.0;
	unordered_map<int, vector<Circle>> obstacles;
	double obstacleDensity = 0.005;  // per surface unit
	double W = 80.0;
	double MAXH = 1.3 * W;
	double dt = 1.0 / 60.0;
	double currentTime = 0;
	double maxObstaclesRadius = 2.0;
	double maxCountdown = 6.0;
	double countdown = maxCountdown;
	double nextReset = 1;
	double prevReset = -1000;
	double step = 20;
	double coef = 1.0;
	double coefIncrement = 0.3;
	bool collided = false;
	const size_t NBSHIPS = 1;
	int seedOffset = 0;

	World() {
		ships.resize(NBSHIPS);
		ships.at(0).position = V(W / 2.0, 0);
	};
	int getSeed(int n) { return n * n + seedOffset; }

	double normalizedDistRay(V direction, double maxDist, const Ship &ship) {
		// direction must be normalized !!
		int gridCell = getGridPosition(ship.position.y);
		Circle *closest = nullptr;
		double closestDist = 1e30;
		// basis change
		// direction is Y, Xdir is X
		V Xdir(-direction.y, direction.x);

		int gridVisibility = (MAXH + 1.0) / gridSize;
		for (int shift = -gridVisibility; shift <= gridVisibility; ++shift) {
			if (obstacles.count(gridCell + shift)) {
				for (auto &o : obstacles[gridCell + shift]) {
					// newPos is o.center in direction basis
					V SE = o.center - ship.position;
					V newPos(SE.dot(Xdir), SE.dot(direction));
					if (newPos.y > 0 && (newPos.x - o.radius) * (newPos.x + o.radius) < 0) {
						// collision
						double dist = newPos.y - sqrt((o.radius * o.radius - newPos.x * newPos.x));
						if (dist > 0 && dist < closestDist) {
							closestDist = dist;
							closest = &o;
						}
					}
				}
			}
		}
		// lateral Walls
		if (direction.x != 0) {
			if (direction.x < 0) {
				// wall 0
				double Yintersect =
				    ship.position.y + direction.y * (-ship.position.x) / direction.x;
				double dist = (V(0, Yintersect) - ship.position).sqLength();
				if (dist < closestDist * closestDist) closestDist = sqrt(dist);
			} else {
				// wall 1
				double Yintersect =
				    ship.position.y + direction.y * (W - ship.position.x) / direction.x;
				double dist = (V(W, Yintersect) - ship.position).sqLength();
				if (dist < closestDist * closestDist) closestDist = sqrt(dist);
			}
		}
		// doors
		double apertureSize = W * countdown / maxCountdown;
		double closedSize = (W - apertureSize) * 0.5;
		if (direction.y != 0) {
			if (direction.y < 0) {
				// prevReset
				double Xintersect =
				    ship.position.x + direction.x * (prevReset - ship.position.y) / direction.y;
				double dist = (V(Xintersect, prevReset) - ship.position).sqLength();
				if (dist < closestDist * closestDist) closestDist = sqrt(dist);
			} else {
				// nextReset
				double Xintersect =
				    ship.position.x + direction.x * (nextReset - ship.position.y) / direction.y;
				if (Xintersect < closedSize || Xintersect > W - closedSize) {
					double dist = (V(Xintersect, nextReset) - ship.position).sqLength();
					if (dist < closestDist * closestDist) closestDist = sqrt(dist);
				}
			}
		}
		return min(closestDist, maxDist);
	}

	void updateObstacles() {
		uniform_real_distribution<double> dist(0.0, 1.0);
		int gridVisibility = (MAXH + 1.0) / gridSize;
		int nbObstacles = static_cast<int>(obstacleDensity * gridSize * W);
		// we need to generate all visible obstacles;
		int currentGridCell = static_cast<int>(floor(ships.at(0).position.y / gridSize));
		for (int visibleCell = currentGridCell - gridVisibility;
		     visibleCell <= currentGridCell + gridVisibility; ++visibleCell) {
			if (visibleCell >= 0 && !obstacles.count(visibleCell)) {
				// a potentially visible grid cell is empty, we need to fill it;
				std::mt19937 generator(getSeed(visibleCell));
				V bottomLeftCorner = V(0, 0) + V(0, 1) * visibleCell * gridSize;
				for (int i = 0; i < nbObstacles; ++i) {
					obstacles[visibleCell].push_back(Circle(
					    V(dist(generator) * W, dist(generator) * gridSize) + bottomLeftCorner,
					    max(0.3, dist(generator)) * maxObstaclesRadius));
				}
			}
		}
	}
	void update() {
		currentTime += dt;
		countdown -= dt;
		V prevPos = ships.at(0).position;
		for (auto &s : ships) {
			s.updatePosition(dt);
		}
		if (ships.at(0).position.y >= nextReset + ships.at(0).dimensions.y * 0.5) {
			prevReset = nextReset;
			nextReset += step * coef;
			coef += coefIncrement;
			countdown = maxCountdown;
		}
		updateObstacles();
		int gridPosition = getGridPosition(ships.at(0).position.y);
		if (ships.at(0).position.x < 0 || ships.at(0).position.x > W) collided = true;
		double apertureSize = W * countdown / maxCountdown;
		double closedSize = (W - apertureSize) * 0.5;
		if (ships.at(0).position.y < prevReset) collided = true;
		if (ships.at(0).position.y < nextReset + 0.5 &&
		    ships.at(0).position.y > nextReset - 0.5) {
			if (ships.at(0).position.x < closedSize || ships.at(0).position.x > W - closedSize)
				collided = true;
		}
		if (obstacles.count(gridPosition)) {
			for (auto &o : obstacles.at(gridPosition)) {
				if ((ships.at(0).position - o.center).sqLength() < pow(o.radius + 0.7, 2)) {
					collided = true;
				}
			}
		}
		ships.at(0).forces = V(0, 0);
	}

	int getGridPosition(double y) { return static_cast<int>(floor(y / gridSize)); }
};

struct shipXP {
	static const constexpr int NBLASERS = 11;
	template <typename G> static G randomInit(size_t nbReguls=1, int ff=0, int ii=0, int nn=0) {
		G g(ff, ii, nn);
		g.randomParams();
		g.addRandomProtein(G::ProteinType_t::input, "c");  // cos angle
		g.addRandomProtein(G::ProteinType_t::input, "s");  // sin angle
		for (int i = 0; i < NBLASERS; ++i)
			g.addRandomProtein(G::ProteinType_t::input, std::to_string(i));

		// turn left
		g.addRandomProtein(G::ProteinType_t::output, "l0");
		g.addRandomProtein(G::ProteinType_t::output, "l1");
		// turn right
		g.addRandomProtein(G::ProteinType_t::output, "r0");
		g.addRandomProtein(G::ProteinType_t::output, "r1");
		// thrust
		g.addRandomProtein(G::ProteinType_t::output, "t0");
		g.addRandomProtein(G::ProteinType_t::output, "t1");
		g.randomReguls(nbReguls);
		return g;
	}

	template <typename I> static void evaluate(I &ind, bool dbg = false) {
		const int NRUN = 2;
		const double TURNSPEED = 8.0;
		const double TETA = M_PI * 1.2;
		auto &g = ind.dna;
		double d = 0;
		for (int r = 0; r < NRUN; ++r) {
#ifdef SHIP_DISPLAY
			QSurfaceFormat f;
			f.setSamples(8);
			int argc = 0;
			QGuiApplication app(argc, nullptr);
#endif
			World world;
			world.seedOffset = r * 1000;
			const double maxDist = world.MAXH;
			auto &s = world.ships.at(0);
			bool finished = false;
			auto stepFunc = [&]() {
				auto dir = s.orientation;
				g.setInputConcentration("c", dir.x * 0.5 + 0.5);
				g.setInputConcentration("s", dir.y * 0.5 + 0.5);
				dir.rotate((-TETA / 2.0) - 0.5 * (TETA / static_cast<double>(NBLASERS)));
				for (int i = 0; i < NBLASERS; ++i) {
					dir.rotate(TETA / static_cast<double>(NBLASERS));
					dir.normalize();
					double dist = world.normalizedDistRay(dir, maxDist, s);
					g.setInputConcentration(std::to_string(i), dist);
				}
				g.step();
				bool tleft = g.getOutputConcentration("l0") > g.getOutputConcentration("l1");
				bool tright = g.getOutputConcentration("r0") > g.getOutputConcentration("r1");
				bool thrust = g.getOutputConcentration("t0") > g.getOutputConcentration("t1");
				if (tleft && !tright)
					s.rotate(1.0, world.dt * TURNSPEED);
				else if (!tleft && tright)
					s.rotate(-1.0, world.dt * TURNSPEED);
				if (thrust) s.thrust(world.dt);
				world.update();
				finished = world.collided || world.countdown <= 0;
				if (finished) {
#ifdef SHIP_DISPLAY
					exit(0);
					app.exit();
#endif
					return true;
				} else
					return false;
			};
#ifdef SHIP_DISPLAY
			ShipWindow<World> window(world, stepFunc);
			window.setFormat(f);
			window.resize(900, 900);
			window.show();
			window.setAnimating(true);
			app.exec();
#else
			while (!finished) {
				stepFunc();
			}
#endif
			d += s.position.y;
		}
		ind.fitnesses["distance"] = d / static_cast<double>(NRUN);

#ifdef SHIP_DISPLAY
		std::cerr << "Fitness = " << ind.fitnesses["distance"] << std::endl;
#endif
	}
};
}
#endif
