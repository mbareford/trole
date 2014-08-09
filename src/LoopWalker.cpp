/**
 * @file LoopWalker.cpp
 *
 * This file defines the methods of the CLoopWalker class.
 */

#include "StdInclude.h"

#include "Exception.h"

#include "GSL.h"
#include "Profile.h"
#include "Threshold.h"
#include "RelaxedLoop.h"
#include "Loop.h"
#include "LoopWalker.h"
#include "ParameterSpace.h"


const double PI = CGSL::PI;


/**
 * All attributes are zeroed.
 */
CLoopWalker::CLoopWalker(void): m_initialised(false) {
	reset();
}

/**
 * Set attributes to zero.
 */
void CLoopWalker::reset(void) {
	m_step_size = 0.0;
	memset(&m_start, 0, sizeof(m_start));
	memset(&m_last_step, 0, sizeof(m_last_step));
	m_step_cnt = 0.0;

	m_started = false;
	m_correlated = false;
	m_init_std_dev = 0.0;
	m_std_dev = 0.0;
}

/**
 * Initialises the loop pointer and step size attributes.
 *
 * Other attributes are zeroed.
 *
 * If this function completes successfully the initialisation flag is
 * set to true.
 *
 * @param _step_size the size of each step in the walk.
 * @throw CException should any of the following occur:
 *        1. zero or negative step size.
 */
void CLoopWalker::init(double _step_size) {
	if (_step_size <= 0.0) {
		throw CException("Error! Zero or negative step size.");
	}

	m_initialised = false;
	reset();

	m_step_size = _step_size;

	m_initialised = true;
}

/**
 * Sets up the walk so that changes in y are correlated with changes in x.
 *
 * @param _init_std_dev the standard deviation of the normal distribution from which the starting x coordinate is chosen (0 is the mean)
 * @param _std_dev the standard deviation of the normal distribution from which all over coordinates are chosen
 *
 * @throw CException should any of the following occur:
 *        1. loop Walker has not been initialised;
 */
void CLoopWalker::correlate(double _init_std_dev, double _std_dev) {
	if (!m_initialised) {
		throw CException("Error! Loop Walker has not been initialised.");
	}

	m_correlated = true;
	m_init_std_dev = _init_std_dev;
	m_std_dev = _std_dev;
}

/**
 * Set a starting position for the walk.
 *
 * @param _x the x coordinate for starting position.
 * @param _y the y coordinate for starting position.
 *
 * @throw CException should any of the following occur:
 *        1. loop Walker has not been initialised.
 */
void CLoopWalker::start(double _x, double _y) {
	if (!m_initialised) {
		throw CException("Error! Loop Walker has not been initialised.");
	}

	m_started = false;
	m_step_cnt = 0;

	m_start.x = _x;
	m_start.y = _y;

	m_last_step.p1.x = m_start.x;
	m_last_step.p1.y = m_start.y;

	m_last_step.p2.x = m_last_step.p1.x;
	m_last_step.p2.y = m_last_step.p1.y;

	m_started = true;
}

/**
 * Restart a walk at the last starting position.
 *
 * @throw CException should any of the following occur:
 *        1. loop Walker has not been initialised.
 */
void CLoopWalker::restart(void) {
	if (!m_initialised) {
		throw CException("Error! Loop Walker has not been initialised.");
	}

	m_started = false;
	m_step_cnt = 0;

	m_last_step.p1.x = m_start.x;
	m_last_step.p1.y = m_start.y;

	m_last_step.p2.x = m_last_step.p1.x;
	m_last_step.p2.y = m_last_step.p1.y;

	m_started = true;
}

/**
 * Choose an appropriate starting position according to whether or not
 * walk is random or correlated.
 *
 * @param _min_x the minimum x coordinate for starting position.
 * @param _max_x the maximum x coordinate for starting position.
 * @param _min_y the minimum y coordinate for starting position.
 * @param _max_y the maximum y coordinate for starting position.
 *
 * @throw CException should any of the following occur:
 *        1. loop Walker has not been initialised;
 *        2. bounds for starting position are invalid.
 */
void CLoopWalker::start(double _min_x, double _max_x, double _min_y, double _max_y) {
	if (!m_initialised) {
		throw CException("Error! Loop Walker has not been initialised.");
	}

	if (_min_x >= _max_x || _min_y >= _max_y) {
		throw CException("Error! Bounds for starting position are invalid.");
	}


	m_started = false;
	m_step_cnt = 0;

	// ensure starting position is within bounds
	if (m_correlated) {
		double x, y;
		do {
			// mean x coordinate is 0; the higher the variance the more likely a loop emerges with twist
			x = CGSL::random_normal(0.0, m_init_std_dev);
			// ensure y coordinate is correlated with x coordinate
			y = CGSL::random_normal(x, m_std_dev);
		} while (x < _min_x || x > _max_x || y < _min_y || y > _max_y);

		m_start.x = x;
		m_start.y = y;
	}
	else {
		m_start.x = CGSL::random(_min_x, _max_x);
		m_start.y = CGSL::random(_min_y, _max_y);
	}

	m_last_step.p1.x = m_start.x;
	m_last_step.p1.y = m_start.y;
	m_last_step.p2.x = m_last_step.p1.x;
	m_last_step.p2.y = m_last_step.p1.y;

	m_started = true;
}


/**
 * Take one step.
 *
 * @throw CException should any of the following occur:
 *        1. Loop Walker has not started a walk.
 */
void CLoopWalker::step(void) {
	if (!m_started) {
		throw CException("Error! Loop Walker has not started a walk.");
	}

	// last position becomes the previous position
	m_last_step.p1.x = m_last_step.p2.x;
	m_last_step.p1.y = m_last_step.p2.y;

	if (m_correlated) {
		correlated_step();
	}
	else {
		random_step();
	}

	m_step_cnt++;
}

/**
 * Take one step in a randomly chosen direction.
 */
void CLoopWalker::random_step(void) {
	double dir = CGSL::random(0.0, 360.0);

	// determine the end point of the next step
	if (dir >= 0.0 && dir < 90.0) {
		// first quadrant
		m_last_step.p2.x += m_step_size*cos((PI/180.0)*dir);
		m_last_step.p2.y += m_step_size*sin((PI/180.0)*dir);
	}
	else if (dir >= 90.0 && dir < 180.0) {
		// second quadrant
		m_last_step.p2.x -= m_step_size*cos((PI/180.0)*(180.0-dir));
		m_last_step.p2.y += m_step_size*sin((PI/180.0)*(180.0-dir));
	}
	else if (dir >= 180.0 && dir < 270.0) {
		// third quadrant
		m_last_step.p2.x -= m_step_size*sin((PI/180.0)*(270.0-dir));
		m_last_step.p2.y -= m_step_size*cos((PI/180.0)*(270.0-dir));
	}
	else {
		// fourth quadrant
		m_last_step.p2.x += m_step_size*cos((PI/180.0)*(360.0-dir));
		m_last_step.p2.y -= m_step_size*sin((PI/180.0)*(360.0-dir));
	}
}

/**
 * Take one step such that the change in y is correlated with the change
 * in x. The change in x has a mean of m_step_size/sqrt(2).
 */
void CLoopWalker::correlated_step(void) {
	// choose dx from a normal distribution
	double dx = CGSL::random_normal(m_step_size/sqrt(2), m_std_dev);

	// randomly decide if the sign of the change should be negated
	if (CGSL::random_binary()) {
		dx *= -1.0;
	}	

	// choose dy from a normal distribution of mean dx
	double dy = CGSL::random_normal(dx, m_std_dev);

	m_last_step.p2.x += dx;
	m_last_step.p2.y += dy;
}

/**
 * Stop the walk so that one of the start methods must be called
 * before a new walk can begin.
 *
 * @throw CException should any of the following occur:
 *        1. Loop Walker has not started a walk.
 */
void CLoopWalker::stop(void) {
	if (!m_started) {
		throw CException("Error! Loop Walker has not started a walk.");
	}

	m_started = false;
	m_step_cnt = 0;
	memset(&m_start, 0, sizeof(m_start));
	memset(&m_last_step, 0, sizeof(m_last_step));
}


/**
 * Return the loop's current position.
 *
 * @throw CException should any of the following occur:
 *        1. Loop Walker has not started a walk.
 */
pt CLoopWalker::get_pos(void) const {
	if (!m_started) {
		throw CException("Error! Loop Walker has not started a walk.");
	}

	pt pos = {m_last_step.p2.x, m_last_step.p2.y};
	return pos;
}

/**
 * Return the coordinates of the last step taken by the loop.
 *
 * @throw CException should any of the following occur:
 *        1. Loop Walker has not started a walk.
 */
seg CLoopWalker::get_last_step(void) const {
	if (!m_started) {
		throw CException("Error! Loop Walker has not started a walk.");
	}

	return m_last_step;
}

/**
 * Return the number of steps taken during the loop's walk.
 *
 * @throw CException should any of the following occur:
 *        1. Loop Walker has not started a walk.
 */
unsigned long CLoopWalker::get_step_cnt() const {
	if (!m_started) {
		throw CException("Error! Loop Walker has not started a walk.");
	}

	return m_step_cnt;
}

