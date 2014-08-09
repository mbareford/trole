/**
 * @file LoopWalker.h
 *
 * This file declares the CLoopWalker class.
 */

#ifndef LOOPWALKER_H_
#define LOOPWALKER_H_

/**
 * The CLoopWalker class controls how a loop performs a walk through a
 * 2D parameter space. The walk always has a fixed step size; however,
 * the direction may be random or correlated. If the latter, the change
 * in one parameter is correlated with the change in the second parameter
 * during each step.
 */
class CLoopWalker {

public:
	CLoopWalker(void);

	void init(double _step_size);
	void correlate(double _init_std_dev, double _std_dev);

	void start(double _x, double _y);
	void start(double _min_x, double _max_x, double _min_y, double _max_y);
	void restart(void);
	void step(void);
	void stop(void);

	bool started(void) const { return m_started; }
	pt get_pos(void) const;
	seg get_last_step(void) const;
	unsigned long get_step_cnt() const;

private:
	/// if true the object has been initialised, i.e., the init method has been called successfully
	bool m_initialised;
	/// if true a walk has started, i.e., the start method has been called successfully
	bool m_started;

	/// the size of each step in the walk
	double m_step_size;

	/// the starting point
	pt m_start;

	/// the last step taken
	seg m_last_step;
	/// the number of steps taken since the last start
	unsigned long m_step_cnt;

	/// if true, the walk is not random, but changes in y are correlated with changes in x
	bool m_correlated;
	/// the standard deviation of the normal distribution from which the starting x coordinate is chosen (0 is the mean)
	double m_init_std_dev;
	/// the standard deviation of the normal distribution from which all over coordinates are chosen
	double m_std_dev;

	void reset(void);
	void random_step(void);
	void correlated_step(void);
};

#endif /* LOOPWALKER_H_ */
