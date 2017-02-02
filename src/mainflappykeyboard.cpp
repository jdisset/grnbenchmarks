#define DISPLAY
#include <iostream>
#include <sstream>
#include "bestiau.hpp"

int main(int, char **) {
	const int upK = 122;
	const int dwK = 115;
	bool upPressed = false;
	World world;

	setlocale(LC_ALL, "en_US.UTF-8");
	initscr();
	timeout(0);
	noecho();
	curs_set(FALSE);

	bool prevUp = false;
	bool up = false;
	while (world.bestiau.vivant) {
		auto start = std::chrono::steady_clock::now();
		auto c = getch();
		up = c == upK;
		if (!prevUp && c == upK) world.bestiauUp();
		prevUp = up;
		flushinp();
		world.update();
		display(world);
		auto end = std::chrono::steady_clock::now();
		auto diff = end - start;
		auto target = ms((int)(world.dt * 1000));
		std::this_thread::sleep_for(
		    target - std::chrono::duration_cast<std::chrono::milliseconds>(diff));
	}
	endwin();
	std::cout << " SCORE = " << world.dist;
	return 0;
}
