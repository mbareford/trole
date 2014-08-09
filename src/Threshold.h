/**
 * @file Threshold.h
 *
 * This file declares the CThreshold class.
 */

#ifndef THRESHOLD_H_
#define THRESHOLD_H_

typedef struct {
	double x; // a1
	double y; // a2
	unsigned short type; // THRESHOLD_TYPE_INSTABILITY or THRESHOLD_TYPE_BZREVERSAL
	double ext1; // relaxation radius
} pt, vec;

typedef struct {
	pt p1;
	pt p2;
} seg;

/**
 * The CThreshold class specifies a 2 dimensional closed boundary.
 * Within the threshold a system is stable, without it is not.
 * The threshold is approximated by a polygon - it comprises a set of vertices that
 * must be in counter-clockwise order.
 */
class CThreshold {

public:
	CThreshold(void) {}
	~CThreshold(void) {}

	/// the vertices of the threshold
	vector<pt> m_vertex;

	bool contains(pt _pt);
	unsigned long crossed_by(seg _seg, pt& _crs_pt);

	void translate_pt(const CThreshold& _src_thres, unsigned long _src_pt_i, const pt& _src_pt, pt& _trs_pt);

	// constants for determining the threshold type of individual points in m_vertex
	static const unsigned short THRESHOLD_TYPE_UNKNOWN;
	static const unsigned short THRESHOLD_TYPE_INSTABILITY;
	static const unsigned short THRESHOLD_TYPE_BZREVERSAL;

private:
	/**
	 * this constant determines how close to parallel two vectors need to be before they
	 * are judged to be parallel
	 */
	static const double PARALLELISM_TOLERANCE;

	void init_pt(pt& _p1, pt& _p2);
	bool pts_equal(pt& _p1, pt& _p2);

	void init_vec(pt& _p1, pt& _p2, vec& _v);
	double dot_vec(vec& _v1, vec& _v2);
	double cross_vec(vec& _v1, vec& _v2);
	double vec_product(unsigned int type, pt& _p1, pt& _p2);

	double get_len(const pt& _p1, const pt& _p2);

	bool pt_within_seg(seg& _sg, pt& _pt);
	unsigned short seg_intersection(seg& _s1, seg& _s2, pt& _i1, pt& _i2);
};

#endif /* THRESHOLD_H_ */
