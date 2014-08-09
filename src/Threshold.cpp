/**
 * @file Threshold.cpp
 *
 * This file defines the methods of the CThreshold class.
 */


#include "StdInclude.h"

#include "Exception.h"
#include "Threshold.h"


const double CThreshold::PARALLELISM_TOLERANCE = 0.00000001; // avoid division overflow

const unsigned short CThreshold::THRESHOLD_TYPE_UNKNOWN     = 0x00;
const unsigned short CThreshold::THRESHOLD_TYPE_INSTABILITY = 0x01;
const unsigned short CThreshold::THRESHOLD_TYPE_BZREVERSAL  = 0x02;


/**
 * This function initialises a point (p1) with the values of another
 * point (p2).
 *
 * @param _p1 the first point
 * @param _p2 the second point
 */
void CThreshold::init_pt(pt& _p1, pt& _p2) {
	do {
		_p1.x = _p2.x;
		_p1.y = _p2.y;
	} while (0);
}

/**
 * This function returns true if a two points are equal, otherwise false.
 *
 * @param _p1 the first point
 * @param _p2 the second point
 */
bool CThreshold::pts_equal(pt& _p1, pt& _p2) {
	bool equal = false;

	do {
		if (_p1.x == _p2.x && _p1.y == _p2.y) {
			equal = true;
		}

	} while (0);

	return equal;
}

/**
 * This function initialises a vector from two points.
 *
 * @param _p1 the first point
 * @param _p2 the second point
 * @param _v the vector
 */
void CThreshold::init_vec(pt& _p1, pt& _p2, vec& _v) {
	do {
		_v.x = _p2.x - _p1.x;
		_v.y = _p2.y - _p1.y;
	} while (0);
}

/**
 * This function multiplies two vectors together (using dot product
 * method) and returns the result.
 *
 * @param _p1: the first point
 * @param _p2: the second point
 */
double CThreshold::dot_vec(vec& _v1, vec& _v2) {
	double result = 0.0;

	do {
		result = (_v1.x * _v2.x) + (_v1.y * _v2.y);
	} while (0);

	return result;
}

/**
 * This function multiplies two vectors together (using cross product
 * method) and returns the result.
 *
 * @param _p1 the first point
 * @param _p2 the second point
 */
double CThreshold::cross_vec(vec& _v1, vec& _v2) {
	double result = 0.0;

	do {
		result = (_v1.x * _v2.y) - (_v1.y * _v2.x);
	} while (0);

	return result;
}

/**
 * Calculate the length between two points.
 *
 * @param _p1 the first point
 * @param _p2 the second point
 */
double CThreshold::get_len(const pt& _p1, const pt& _p2) {
	double dx = fabs(_p1.x - _p2.x);
	double dy = fabs(_p1.y - _p2.y);

	return sqrt(dx*dx + dy*dy);
}


/*
 * The implementions of the following methods are based on algorithms available
 * from the softsurfer website.
 *
 * point_within_segment()
 * point_within_threshold()
 * segment_intersection()
 * segment_intersect_threshold()
 *
 * The appropriate copyright notice follows.
 ********************************************************************************
 ** Copyright 2001, softSurfer (www.softsurfer.com)
 ** This code may be freely used and modified for any purpose
 ** providing that this copyright notice is included with it.
 ** SoftSurfer makes no warranty for this code, and cannot be held
 ** liable for any real or imagined damage resulting from its use.
 ** Users of this code must verify correctness for their application.
 ********************************************************************************
 */

/**
 * Determine if a point is inside a segment.
 *
 * @param _sg the segment
 * @param _pt the point
 * @return true if the point and segment do intersect, otherwise false
 */
bool CThreshold::pt_within_seg(seg& _sg, pt& _pt) {
	bool within = false;

	do {
		if (_sg.p1.x != _sg.p2.x) {    // S is not vertical
			if (_sg.p1.x <= _pt.x && _pt.x <= _sg.p2.x) {
				within = true;
			}
			else if (_sg.p1.x >= _pt.x && _pt.x >= _sg.p2.x) {
		    	within = true;
		    }
		}
		else {    // S is vertical, so test y coordinate
			if (_sg.p1.y <= _pt.y && _pt.y <= _sg.p2.y) {
				within = true;
			}
			else if (_sg.p1.y >= _pt.y && _pt.y >= _sg.p2.y) {
				within = true;
			}
		}
	} while (0);


    return within;
}

/**
 * This function returns true if a specified point falls within a
 * specified polygon, otherwise false.
 *
 * @param _pt the point
 */
bool CThreshold::contains(pt _pt) {
	double vt = 0.0;
	unsigned int cn = 0;  // the crossing number counter
	unsigned long i;

	do {
		// loop through all edges of the threshold
		for (i = 0; i < (m_vertex.size()-1); ++i) {
			// select downward and upward crossings
			if ((m_vertex[i].y <= _pt.y && m_vertex[i+1].y > _pt.y)
		        ||
		        (m_vertex[i].y > _pt.y && m_vertex[i+1].y <= _pt.y)) {

				// compute the actual edge-ray intersect x-coordinate
		        vt = (_pt.y - m_vertex[i].y) / (m_vertex[i+1].y - m_vertex[i].y);
		        if (_pt.x < (m_vertex[i].x + (vt*(m_vertex[i+1].x - m_vertex[i].x))))  {
		        	++cn;
		        }
		    }
		}

	} while (false); // end of main body loop

	return (cn & 1);
}

/**
 * Determine if one segment intersects another.
 *
 * @param _s1 the threshold edge
 * @param _s2 the segment
 * @param _i1 the possible intersection point
 * @param _i2 the possible intersection point
 * @return 0 if the segments do not intersect
 *         1 if one segment intersects another at precisely one point
 *         2 if one segment overlaps another, in which case there are two intersection points.
 * @throw CException should any of the following occur:
 *        1. zero v.y.
  */
unsigned short CThreshold::seg_intersection(seg& _s1, seg& _s2, pt& _i1, pt& _i2 ) {
	unsigned short intersect = 0;
	vec u, v, w, w2;
    double d, du, dv;
    double t, t0, t1;
    double si, ti;

    do {
    	init_vec(_s1.p1, _s1.p2, u);
		init_vec(_s2.p1, _s2.p2, v);
		init_vec(_s2.p1, _s1.p1, w);

    	d = cross_vec(u, v);

        if (fabs(d) < PARALLELISM_TOLERANCE) {
            // segments are parallel
        	if (cross_vec(u, w) != 0 || cross_vec(v, w) != 0) {
                // not colinear
        		break;
    	    }

        	// they are colinear or degenerate
    	    // check if they are degenerate points
    	    du = dot_vec(u, u);
    	    dv = dot_vec(v, v);

    	    if (0 == du && 0 == dv) {
    	    	// both segments are points
    	    	if (pts_equal(_s1.p1, _s2.p1)) {
    	    		// the same point
    	    		intersect = 1;
    	    		init_pt(_i1, _s1.p1);
    	        }

    	        break;
    	    }

    	    if (0 == du) {
    	    	// s1 is a single point
    	    	if (pt_within_seg(_s2, _s1.p1)) {
    	    		// is within s2
    	    		intersect = 1;
    	    		init_pt(_i1, _s1.p1);
    	    	}

    	        break;
    	    }

    	    if (0 == dv) {
    	    	// s2 a single point
    	    	if (pt_within_seg(_s1, _s2.p1)) {
    	    		// is within s1
    	    		intersect = 1;
    	    		init_pt(_i1, _s2.p1);
    	        }

    	        break;
    	    }

    	    // colinear segments - get overlap (or not)
    	    init_vec(_s2.p1, _s1.p2, w2);

    	    if (0 != v.x) {
    	    	t0 = w.x/v.x;
    	        t1 = w2.x/v.x;
    	    }
    	    else {
    	    	if (0.0 == v.y) {
    	    		throw CException("Error! Zero v.y.");
    	    	}

    	        t0 = w.y/v.y;
    	        t1 = w2.y/v.y;
    	    }

            // must have t0 smaller than t1
    	    if (t0 > t1) {
    	        // swap if not
    	    	t  = t0;
    	        t0 = t1;
    	        t1 = t;
    	    }

    	    if (t0 > 1 || t1 < 0) {
    	        break; // no overlap
    	    }

    	    t0 = t0 < 0 ? 0 : t0; // clip to min 0
    	    t1 = t1 > 1 ? 1 : t1; // clip to max 1

    	    if (t0 == t1) {
    	    	// intersection is a pt
    	    	intersect = 1;
    	    	_i1.x = _s2.p1.x + (t0*v.x);
    	    	_i1.y = _s2.p1.y + (t0*v.y);
    	    	break;
    	    }

    	    // they overlap in a valid subseg
    	    intersect = 2;
    	    _i1.x = _s2.p1.x + (t0*v.x);
    	    _i1.y = _s2.p1.y + (t0*v.y);
    	    _i2.x = _s2.p1.x + (t1*v.x);
    	    _i2.y = _s2.p1.y + (t1*v.y);
    	    break;
        } // end of <if (fabs(d) < PARALLELISM_TOLERANCE)> clause

        // the segs are skew and may intersect in a pt
        // get the intersect parameter for s1
        si = cross_vec(v, w)/d;
        if (si < 0 || si > 1) {
        	// no intersect with s1
        	break;
        }

        // get the intersect parameter for S2
        ti = cross_vec(u, w)/d;
        if (ti < 0 || ti > 1) {
        	// no intersect with s2
        	break;
        }

        intersect = 1;
        _i1.x = _s1.p1.x + (si*u.x);
        _i1.y = _s1.p1.y + (si*u.y);

    } while (0);

    return intersect;
}

/**
 * This method returns the m_vector index of the first point of the threshold edge that is
 * intersected by the specified segment, otherwise the method returns a value one higher than
 * the maximum m_vector index (i.e., m_vector.size())
 *
 * If the segment does intersect a threshold edge the intersection
 * point is returned by the _crs_pt parameter.
 *
 * NB: the threshold MUST be convex and have vertices oriented
 *     counterclockwise; this code does not verify these conditions.
 *
 * @param _seg the segment
 * @param _crs_pt the intersection point
 * @return the index of the first point of the edge that is intersected
 *         if no intersection is found the index of the last point of
 *         the last edge is returned
 * @throw CException should any of the following occur:
 *        1. segment is really a single point.
 */
unsigned long CThreshold::crossed_by(seg _seg, pt& _crs_pt) {
    unsigned long i = 0; // for loop index

    seg edge; // threshold edge
    pt i1, i2; // intersection points

    _crs_pt.x = 0.0;
    _crs_pt.y = 0.0;

	if (pts_equal(_seg.p1, _seg.p2)) {
		throw CException("Error! Segment is really a single point.");
	}

	// process each threshold edge
	for (i = 0; i < m_vertex.size()-1; i++) {
		edge.p1.x = m_vertex[i].x;
		edge.p1.y = m_vertex[i].y;
		edge.p2.x = m_vertex[i+1].x;
		edge.p2.y = m_vertex[i+1].y;

		if (seg_intersection(edge, _seg, i1, i2)) {
			// there is a valid intersection point
			_crs_pt.x = i1.x;
			_crs_pt.y = i1.y;

			if ((THRESHOLD_TYPE_BZREVERSAL == m_vertex[i].type && (THRESHOLD_TYPE_BZREVERSAL & m_vertex[i+1].type))
				||
				((THRESHOLD_TYPE_BZREVERSAL & m_vertex[i].type) && THRESHOLD_TYPE_BZREVERSAL == m_vertex[i+1].type)
				||
				((THRESHOLD_TYPE_BZREVERSAL & m_vertex[i].type) && (THRESHOLD_TYPE_BZREVERSAL & m_vertex[i+1].type))
				||
				(THRESHOLD_TYPE_BZREVERSAL == m_vertex[i].type || THRESHOLD_TYPE_BZREVERSAL == m_vertex[i+1].type)) {

				_crs_pt.type = THRESHOLD_TYPE_BZREVERSAL;
				_crs_pt.ext1 = 0.0;
			}
			else {
				_crs_pt.type = THRESHOLD_TYPE_INSTABILITY;

				// if ever there is a requirement for the relaxation radius to be
				// decided by where the instability threshold is crossed then the
				// commented-out code below may come in useful

				// calculate the relaxation radius for the intersection point
				//double edge_len = get_len(m_vertex[i], m_vertex[i+1]);
				//double ipt_to_crspt_len = get_len(m_vertex[i], _crs_pt);

				//_crs_pt.ext1 = m_vertex[i].ext1 + (ipt_to_crspt_len/edge_len)*(m_vertex[i+1].ext1 - m_vertex[i].ext1);
			}

			break;
		} // end of <if (seg_intersection(edge, _seg, i1, i2))> clause

	} // end of <for (i = 0; i < (m_vertex.size()-1) && !intersect; i++)> loop

	return i;
}

/**
 * Translate a point on one threshold to a point on this threshold. The specified threshold,
 * _src_thres, must have the same number of points as this threshold.
 *
 * @param _src_thres the source threshold
 * @param _src_pt_i the index of the vertex that occurs immediately before the point to be translated
 * @param _src_pt the point to be translated from the coordinates of the source threshold to the
 *                coordinates of this threshold
 * @param _trs_pt the translated point
 * @throw CException should any of the following occur:
 *        1. this threshold contains no vertices.
 *        2. source threshold does not have the same number of vertices as this threshold.
 */
void CThreshold::translate_pt(const CThreshold& _src_thres, unsigned long _src_pt_i, const pt& _src_pt, pt& _trs_pt) {
	if (0 == m_vertex.size()) {
		throw CException("Error! Threshold has no vertices.");
	}

	if (_src_thres.m_vertex.size() != m_vertex.size()) {
		throw CException("Error! Threshold does not have the same number of vertices as the source threshold.");
	}

	unsigned long i1 = _src_pt_i;
	unsigned long i2 = (i1 + 1 == _src_thres.m_vertex.size() ? 1 : i1 + 1);

	// source segment
	seg src_seg = {{_src_thres.m_vertex[i1].x, _src_thres.m_vertex[i1].y},
				   {_src_thres.m_vertex[i2].x, _src_thres.m_vertex[i2].y}};

	// translated segment
	seg trs_seg = {{m_vertex[i1].x, m_vertex[i1].y},
			       {m_vertex[i2].x, m_vertex[i2].y}};

	// do the translation
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double src_seg_len = get_len(src_seg.p1, src_seg.p2);
	double trs_seg_len = get_len(trs_seg.p1, trs_seg.p2);

	double trs_ds = 0.0;
	double trs_ds_x = 0.0, trs_ds_y = 0.0;

	double trs_x = fabs(trs_seg.p1.x - trs_seg.p2.x);
	double trs_y = fabs(trs_seg.p1.y - trs_seg.p2.y);
	double theta = 0.0;

	if (0.0 != src_seg_len && 0.0 != trs_seg_len) {
		trs_ds = (get_len(src_seg.p1, _src_pt)/src_seg_len)*trs_seg_len;
	}

	if (0.0 == trs_ds) {
		trs_ds_x = 0.0;
		trs_ds_y = 0.0;
	}
	else if (trs_seg.p2.x > trs_seg.p1.x) {
		if (trs_seg.p2.y > trs_seg.p1.y) {
			theta = atan2(trs_y, trs_x);
			trs_ds_x = trs_ds*cos(theta);
			trs_ds_y = trs_ds*sin(theta);
		}
		else if (trs_seg.p2.y < trs_seg.p1.y) {
			theta = atan2(trs_x, trs_y);
			trs_ds_x = trs_ds*sin(theta);
			trs_ds_y = trs_ds*cos(theta)*(-1.0);
		}
		else {
			// segment is from left to right
			trs_ds_x = trs_ds;
			trs_ds_y = 0.0;
		}
	}
	else if (trs_seg.p2.x < trs_seg.p1.x) {
		if (trs_seg.p2.y > trs_seg.p1.y) {
			theta = atan2(trs_y, trs_x);
			trs_ds_x = trs_ds*cos(theta)*(-1.0);
			trs_ds_y = trs_ds*sin(theta);
		}
		else if (trs_seg.p2.y < trs_seg.p1.y) {
			theta = atan2(trs_x, trs_y);
			trs_ds_x = trs_ds*sin(theta)*(-1.0);
			trs_ds_y = trs_ds*cos(theta)*(-1.0);
		}
		else {
			// segment is from right to left
			trs_ds_x = trs_ds*(-1.0);
			trs_ds_y = 0.0;
		}
	}
	else {
		// trs_seg.p2.x == trs_seg.p1.x
		if (trs_seg.p2.y > trs_seg.p1.y) {
			// segment is going straight up
			trs_ds_x = 0.0;
			trs_ds_y = trs_ds;
		}
		else if (trs_seg.p2.y < trs_seg.p1.y) {
			// segment is going straight down
			trs_ds_x = 0.0;
			trs_ds_y = trs_ds*(-1.0);
		}
		else {
			// segment is just a point
			trs_ds_x = 0.0;
			trs_ds_y = 0.0;
		}
	}

	_trs_pt.x = trs_seg.p1.x + trs_ds_x;
	_trs_pt.y = trs_seg.p1.y + trs_ds_y;
	_trs_pt.type = _src_pt.type;
	_trs_pt.ext1 = _src_pt.ext1;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
