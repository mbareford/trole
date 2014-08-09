/**
 * @file Profile.cpp
 *
 * This file defines the methods of the CProfile class.
 */

#include "StdInclude.h"

#include "Exception.h"
#include "Profile.h"

/**
 * The default constructor ensures that all the attributes are
 * properly initialised.
 */
CProfile::CProfile(void): m_initialised(false) {
	reset();
}

/**
 * Copy Constructor
 */
CProfile::CProfile(const CProfile& _pf): m_initialised(false) {
	unsigned long i, max_i = _pf.get_plot_cnt()-1;

	if (_pf.m_initialised) {
		init(_pf.get_plot_pos(0), _pf.get_plot_pos(max_i), max_i);

		for (i = 0; i <= max_i; i++) {
			plot(_pf.get_plot_val(i));
		}
	}
}

/**
 * The default destructor ensures that coordinate arrays are deleted.
 */
CProfile::~CProfile(void) {
	m_x.clear();
	m_y.clear();
}

/**
 * The assignment operator.
 */
CProfile& CProfile::operator=(const CProfile& _pf) {
	if (this != &_pf) { // protect against self-assignment
		// get rid of old data
		m_x.clear();
		m_y.clear();

		// assign the new data
		if (_pf.get_plot_cnt() > 0) {
			unsigned long i, max_i = _pf.get_plot_cnt()-1;

			init(_pf.get_plot_pos(0), _pf.get_plot_pos(max_i), max_i);

			for (i = 0; i <= max_i; i++) {
				plot(_pf.get_plot_val(i));
			}
		}
	}

    return *this;
}

/**
 * Set all the attributes to zero.
 */
void CProfile::reset(void) {
	m_min_x = m_max_x = 0.0;
	m_pt_count = m_pt_i = 0;

	m_x.clear();
	m_y.clear();
}

/**
 * Initialises the profile - the bounds of the x coordinates are set, as is the
 * profile resolution.
 *
 * @param _min_x the minimum x coordinate of profile.
 * @param _max_x the maximum x coordinate of profile.
 * @param _rs the profile resolution i.e., the number of points that make up the profile.
 * @throw CException should any of the following occur:
 *        1. _min_x >= _max_x;
 *        2. zero _rs.
 */
void CProfile::init(double _min_x, double _max_x, unsigned long _rs) {
	if (_min_x >= _max_x) {
		throw CException("Error! _min_x >= _max_x.");
	}

	if (0.0 == _rs) {
		throw CException("Error! zero _rs.");
	}

	// clear and zero all attributes
	m_initialised = false;
	reset();

	m_min_x = _min_x;
	m_max_x = _max_x;

	// allocate the memory for the coordinate arrays
	/////////////////////////////////////////////////////////////////////////////
	double inc = (m_max_x-m_min_x)/((double) _rs);

	m_pt_count = _rs+1;
	m_x.resize(m_pt_count);

	// initialise the x and y coordinate arrays
	unsigned long i;
	m_x[0] = m_min_x;
	for (i = 1; i < m_pt_count-1; i++) {
		double x = m_x[0] + i*inc;
		m_x[i] = x;
	}
	m_x[m_pt_count-1] = m_max_x;

	m_y.resize(m_pt_count);
	for (i = 0; i < m_pt_count; i++) {
		m_y[i] = 0.0;
	}
	/////////////////////////////////////////////////////////////////////////////

	m_initialised = true;
}

/**
 * Associates the specified plot value (or y coordinate) with the current plot position
 * (or x coordinate) - the plot position is subsequently incremented.
 *
 * @param _y the plot value.
 * @return the plot value.
 * @throw CException should object not be initialised.
 */
double CProfile::plot(double _y) {
	if (!m_initialised) {
		throw CException("Error! Profile uninitialised.");
	}

	if (m_pt_i < m_pt_count) {
		m_y[m_pt_i] = _y;
		m_pt_i++;
	}

	return _y;
}

/**
 * This method allows a profile object to be interpreted as something more
 * than a simple list of plot values.
 *
 * The x vector implicitly defines a list of bins that have the following bounds;
 * m_x[0] <= _x < m_x[1], m_x[1] <= _x < m_x[2],...,m_x[m_pt_count-2] <= _x < m_x[m_pt_count-1].
 * Notice that each bin corresponds to an index of the x vector.
 * m_x[0] corresponds to the m_x[0] <= _x < m_x[1] bin.
 * m_x[1] corresponds to the m_x[1] <= _x < m_x[2] bin.
 * m_x[m_pt_count-2] corresponds to the m_x[m_pt_count-2] <= _x < m_x[m_pt_count-1] bin.
 * m_y[0], m_y[1], m_y[m_pt_count-2] are the counts associated with these bins.
 * m_y[m_pt_count-1] is not used.
 *
 * So, this method locates the bin that _x falls into and then increments the count of that bin.
 *
 * @param _x the value to be associated with a bin.
 * @return the new count of the bin associated with _x.
 * @throw CException should object not be initialised.
 */
unsigned long CProfile::plot_count(double _x) {
	if (!m_initialised) {
		throw CException("Error! Profile uninitialised.");
	}

	unsigned long res = 0;

	if (_x >= m_x[0] && _x < m_x[m_pt_count-1]) {

		unsigned long i;

		for (i = 0; i < m_pt_count-1; i++) {
			if (_x >= m_x[i] && _x < m_x[i+1]) {
				break;
			}
		}

		m_y[i]++;
		res = (unsigned long) m_y[i];
	}

	return res;
}

/**
 * Associates a plot value (specified as a numerator and a denominator) with the current plot position -
 * the plot(double _y) method is only called should _y_den be non-zero otherwise a zero is returned
 * (plot position is still incremented).
 *
 *
 * @param _y_num the numerator of the plot value.
 * @param _y_den the denominator of the plot value.
 * @return the plot value.
 * @throw CException should object not be initialised.
 */
double CProfile::plot(double _y_num, double _y_den) {
	if (!m_initialised) {
		throw CException("Error! Profile uninitialised.");
	}

	double val = 0.0;

	if (0.0 == _y_den) {
		m_pt_i++;
	}
	else {
		val = plot(_y_num/_y_den);
	}

	return val;
}

/**
 * Returns the current plot position (or x coordinate) - the last
 * plot position is returned if the profile has been fully plotted.
 *
 * @throw CException should object not be initialised.
 */
double CProfile::get_plot_pos(void) const {
	if (!m_initialised) {
		throw CException("Error! Profile uninitialised.");
	}

	return done() ? m_x[m_pt_count-1] : m_x[m_pt_i];
}

/**
 * Returns the plot value associated with the specified plot position (or x coordinate).
 *
 * @param _x the plot position.
 * @throw CException should any of the following occur:
 *        1. object not initialised;
 *        2. _x outside the profile position boundaries.
 */
double CProfile::get_plot_val(double _x) const {
	if (!m_initialised) {
		throw CException("Error! Profile uninitialised.");
	}

	if (_x < m_min_x || _x > m_max_x) {
		throw CException("Error! Plot position out of bounds.");
	}

	unsigned long i = 0;
	double val = 0.0;
	bool found = false;

	for (i = 0; i < m_pt_count && !found; ++i) {
		if (_x == m_x[i]) {
			val = m_y[i];
			found = true;
		}
		else {
			if (i > 0 && _x < m_x[i]) {
				double frac;
				frac = (_x - m_x[i-1])/(m_x[i]-m_x[i-1]);
				val = m_y[i-1] + (frac*(m_y[i]-m_y[i-1]));
				found = true;
			}
		}
	}

	return val;
}

/**
 * Returns the plot position associated with the specified index.
 *
 * @param _i the plot index.
 * @throw CException should any of the following occur:
 *        1. object not initialised;
 *        2. _i outside the profile index boundaries.
 */
double CProfile::get_plot_pos(unsigned long _i) const {
	if (!m_initialised) {
		throw CException("Error! Profile uninitialised.");
	}

	if (_i > m_pt_count) {
		throw CException("Error! Plot index out of bounds.");
	}

	return m_x[_i];
}

/**
 * Returns the plot value associated with the specified index.
 *
 * @param _i the plot index.
 * @throw CException should any of the following occur:
 *        1. object not initialised;
 *        2. _i outside the profile index boundaries.
 */
double CProfile::get_plot_val(unsigned long _i) const {
	if (!m_initialised) {
		throw CException("Error! Profile uninitialised.");
	}

	if (_i > m_pt_count) {
		throw CException("Error! Plot index out of bounds.");
	}

	return m_y[_i];
}

/**
 * Returns the minimum plot value with a specified range
 * of x coordinates.
 *
 * @param _min_x the minimum x coordinate.
 * @param _max_x the maximum x coordinate.
 *
 * @throw CException should any of the following occur:
 *        1. object not initialised.
 *        2. parameters out of range.
 */
double CProfile::get_min_plot_val(double _min_x, double _max_x) const {
	if (!m_initialised) {
		throw CException("Error! Profile uninitialised.");
	}

	if (_min_x < m_min_x || _max_x > m_max_x) {
		throw CException("Error! Parameters out of range.");
	}

	bool min_y_initialised = false;
	double min_y = 0.0;
	unsigned long i = 0;

	for (i = 0; i < m_pt_count; i++) {
		if (m_x[i] >= _min_x && m_x[i] <= _max_x) {
			if (!min_y_initialised || min_y > m_y[i]) {
				min_y_initialised = true;
				min_y = m_y[i];
			}
		}
	}

	return min_y;
}

/**
 * Returns the maximum plot value with a specified range
 * of x coordinates.
 *
 * @param _min_x the minimum x coordinate.
 * @param _max_x the maximum x coordinate.
 *
 * @throw CException should any of the following occur:
 *        1. object not initialised.
 *        2. parameters out of range.
 */
double CProfile::get_max_plot_val(double _min_x, double _max_x) const {
	if (!m_initialised) {
		throw CException("Error! Profile uninitialised.");
	}

	if (_min_x < m_min_x || _max_x > m_max_x) {
		throw CException("Error! Parameters out of range.");
	}

	bool max_y_initialised = false;
	double max_y = 0.0;
	unsigned long i = 0;

	for (i = 0; i < m_pt_count; i++) {
		if (m_x[i] >= _min_x && m_x[i] <= _max_x) {
			if (!max_y_initialised || max_y < m_y[i]) {
				max_y_initialised = true;
				max_y = m_y[i];
			}
		}
	}

	return max_y;
}

/**
 * Returns the maximum difference between two profiles up to a specified
 * plot position.
 *
 * @param _ppf the pointer to the other profile
 * @param _weighted if true the values are weighted by position
 * @param _max_pos the maximum position for comparing profile values
 * @param _max_diff_pos the profile position of the maximum difference
 * @throw CException should any of the following occur:
 *        1. null profile object;
 *        2. profile(s) not initialised;
 *        3. profiles have different number of plot points;
 *        4. profiles have different plot position ranges;
 *        5. _max_pos outside range of profiles;
 *        6. profiles have different plot positions.
 */
double CProfile::get_max_diff(CProfile* _ppf, bool _weighted, double _max_pos, double& _max_diff_pos) {
	if (NULL == _ppf) {
		throw CException("Error! Null profile object.");
	}

	if (!m_initialised || !_ppf->is_initialised()) {
		throw CException("Error! Profile(s) uninitialised.");
	}

	if (m_pt_count != _ppf->get_plot_cnt()) {
		throw CException("Error! Profiles have different number of plot points.");
	}

	if (m_min_x != _ppf->get_plot_pos(0) || m_max_x != _ppf->get_plot_pos(_ppf->get_plot_cnt()-1)) {
		throw CException("Error! Profiles have different plot position ranges.");
	}

	if (_max_pos < m_min_x || _max_pos > m_max_x) {
		throw CException("Error! _max_pos outside range of profiles.");
	}

	unsigned long i = 0;
	double max_diff = 0.0;
	double pos = 0.0, diff = 0.0;

	for (i = 0; i < m_pt_count; i++) {
		if (m_x[i] != _ppf->get_plot_pos(i)) {
			throw CException("Error! Profiles have different plot positions.");
		}

		pos = m_x[i];
		if (pos > _max_pos) {
			break;
		}

		diff = fabs(m_y[i]-_ppf->get_plot_val(i));
		if (_weighted) {
			diff *= pos;
		}

		if (diff > max_diff) {
			max_diff = diff;
			_max_diff_pos = pos;
		}
	} // end of <for (i = 0; i < m_pt_count; i++)> loop

	return max_diff;
}
