/**
 * @file ParameterSpace.cpp
 *
 * This file defines the methods of the CParameterSpace class.
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

long CParameterSpace::m_count = 0;

/**
 * Default constructor merely ensures all attributes
 * are initialised.
 */
CParameterSpace::CParameterSpace(void) {
	m_thres.m_vertex.clear();
	m_alpha_thres.m_vertex.clear();

	m_thres_min_x = 0.0;
	m_thres_max_x = 0.0;
	m_thres_min_y = 0.0;
	m_thres_max_y = 0.0;
	m_thres_set = false;
	m_alpha_thres_set = false;

	m_count = 0;
}

/**
 * Initialise the threshold vector. The attribute m_alpha_thres
 * is also initialised if _ptfn is not NULL.
 *
 * @param _atfn the file name of the alpha threshold file
 * @param _ptfn the file name of the parameter threshold file
 * @throw CException should any of the following occur:
 *        1. cannot open alpha threshold file;
 *        2. cannot open parameter threshold file;
 *        3. cannot read threshold file;
 *        4. threshold has no vertices;
 *        5. cannot read alpha threshold file;
 *        6. parameter threshold does not have the same number of
 *           vertices as the alpha threshold.
 */
void CParameterSpace::set_threshold(const char* _atfn, const char* _ptfn/* = NULL*/) {
	ifstream atfs, ptfs;

	if (NULL == _ptfn) {
		ptfs.open(_atfn);
		if (!ptfs.is_open()) {
			throw CException("Error! Cannot open alpha threshold file.");
		}
	}
	else {
		ptfs.open(_ptfn);
		if (!ptfs.is_open()) {
			throw CException("Error! Cannot open parameter threshold file.");
		}
	}

	try {
		m_thres.m_vertex.clear();
		m_thres_min_x = 0.0;
		m_thres_max_x = 0.0;
		m_thres_min_y = 0.0;
		m_thres_max_y = 0.0;
		m_thres_set = false;
		m_alpha_thres_set = false;

		pt p = {};
		bool first_iter = true;

		while (!ptfs.eof()) {
			ptfs >> p.x >> p.y >> p.type;
			if (ptfs.eof()) {
				break;
			}
			else {
				if (!ptfs.good()) {
					throw CException("Error! Cannot read threshold file.");
				}
			}

			m_thres.m_vertex.push_back(p);

			if (first_iter) {
				m_thres_min_x = m_thres_max_x = p.x;
				m_thres_min_y = m_thres_max_y = p.y;
			}
			else {
				if (p.x < m_thres_min_x) {
					m_thres_min_x = p.x;
				}
				else if (p.x > m_thres_max_x) {
					m_thres_max_x = p.x;
				}

				if (p.y < m_thres_min_y) {
					m_thres_min_y = p.y;
				}
				else if (p.y > m_thres_max_y) {
					m_thres_max_y = p.y;
				}
			}

			first_iter = false;
		} // end of <while (!ptfs.eof())> loop

		if (0 == m_thres.m_vertex.size()) {
			throw CException("Error! Threshold has no vertices.");
		}

		ptfs.close();
		m_thres_set = true;

		if (NULL != _ptfn) {
			// read in alpha threshold so that any crossings of the parameter threshold
			// can be translated to alpha space
			m_alpha_thres.m_vertex.clear();

			atfs.open(_atfn);
			if (!atfs.is_open()) {
				throw CException("Error! Cannot open alpha threshold file.");
			}

			while (!atfs.eof()) {
				atfs >> p.x >> p.y >> p.type;
				if (atfs.eof()) {
					break;
				}
				else {
					if (!atfs.good()) {
						throw CException("Error! Cannot read alpha threshold file.");
					}
				}

				m_alpha_thres.m_vertex.push_back(p);
			}

			if (m_alpha_thres.m_vertex.size() != m_thres.m_vertex.size()) {
				throw CException("Error! Parameter threshold does not have the same number of vertices as the alpha threshold.");
			}

			atfs.close();

			m_alpha_thres_set = true;
		} // end of <if (NULL != _ptfn)> clause
	}
	catch (CException& x) {
		if (ptfs.is_open()) {
			ptfs.close();
		}

		if (atfs.is_open()) {
			atfs.close();
		}

		m_thres.m_vertex.clear();
		m_alpha_thres.m_vertex.clear();

		m_thres_set = false;
		m_alpha_thres_set = false;

		throw x;
	}
}

/**
 * A loop's position in parameter space is evolved until it crosses the
 * threshold. The loop undergoes a random walk through the parameter space:
 * each step of length _step_size has random direction.
 *
 * @param _ploop pointer to the loop walker that will traverse stability region
 * @param _crs_pt the point at which the loop crossed the parameter threshold
 * @param _alpha_crs_pt the point at which the loop crossed the alpha threshold
 * @throw CException should any of the following occur:
 *        1. the threshold has not been set;
 *        2. the starting point is not within the threshold;
 *        3. null loop walker pointer;
 *        4. unable to determine crossing point.
 */
void CParameterSpace::traverse_stability_region(CLoopWalker* _plw, pt& _crs_pt, pt& _alpha_crs_pt) {
	if (!m_thres_set) {
		throw CException("Error! No threshold.");
	}

	if (NULL == _plw) {
		throw CException("Error! Null loop walker pointer.");
	}

	if (!_plw->started()) {
		// start walk at a random point within threshold
		do {
			_plw->start(m_thres_min_x, m_thres_max_x, m_thres_min_y, m_thres_max_y);
		} while (!m_thres.contains(_plw->get_pos()));
	}
	// take steps in parameter space until instability threshold is crossed
	/////////////////////////////////////////////////////////////////////////////////////////////
	do {
		// take steps in parameter space until some threshold is crossed
		do {
			_plw->step();
		} while (m_thres.contains(_plw->get_pos()));

		unsigned long crs_i = m_thres.crossed_by(_plw->get_last_step(), _crs_pt);
		if (crs_i >= (m_thres.m_vertex.size()-1)) {
			throw CException("Error! Cannot find crossing point.");
		}

		if (m_alpha_thres_set) {
			m_alpha_thres.translate_pt(m_thres, crs_i, _crs_pt, _alpha_crs_pt);
		}
		else {
			// the parameter space is alpha space
			_alpha_crs_pt.x = _crs_pt.x;
			_alpha_crs_pt.y = _crs_pt.y;
			_alpha_crs_pt.type = _crs_pt.type;
			_alpha_crs_pt.ext1 = _crs_pt.ext1;
		}

		if (CThreshold::THRESHOLD_TYPE_INSTABILITY == _alpha_crs_pt.type) {
			// instability threshold has been crossed, so break out of do while loop
			break;
		}
		else {
			// bz reversal threshold has been crossed
			// return loop to previous start position
			m_count++;
			_plw->restart();
		}

	} while (1);
	/////////////////////////////////////////////////////////////////////////////////////////////
}


/**
 * Returns true if the specified point is within the stability region.
 *
 * @param _x the x coordinate
 * @param _y the y coordinate
 * @param _alpha if true, coordinates are given in alpha space
 * @throw CException should any of the following occur:
 *        1. the threshold has not been set;
 */
bool CParameterSpace::within_threshold(double _x, double _y, bool _alpha = false) {
	if (!m_thres_set) {
		throw CException("Error! No threshold.");
	}

	pt pos = {_x, _y};
	return _alpha && m_alpha_thres_set ? m_alpha_thres.contains(pos) : m_thres.contains(pos);
}

