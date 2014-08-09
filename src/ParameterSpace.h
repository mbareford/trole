/**
 * @file ParameterSpace.h
 *
 * This file declares the CParameterSpace class.
 */

#ifndef PARAMETERSPACE_H_
#define PARAMETERSPACE_H_

/**
 * The CParameterSpace class defines the 2 dimensional space through which a
 * coronal loop evolves.
 */
class CParameterSpace {

public:
	CParameterSpace(void);
	~CParameterSpace(void) {}

	void set_threshold(const char* _atfn, const char* _ptfn = NULL);

	void traverse_stability_region(CLoopWalker* _plw, pt& _crs_pt, pt& _alpha_crs_pt);

	bool within_threshold(double _x, double _y, bool _alpha/*=false*/);

private:
	/// the instability threshold in 2 dimensional parameter space
	CThreshold m_thres;
	/**
	 * The instability threshold in 2 dimensional alpha space. If this attribute
	 * is used there must exist a one-to-one correspondence between the points
	 * of m_thres and m_alpha_thres.
	 */
	CThreshold m_alpha_thres;

	/// the minimum x coordinate of all the points in threshold
	double m_thres_min_x;
	/// the maximum x coordinate of all the points in threshold
	double m_thres_max_x;
	/// the minimum y coordinate of all the points in threshold
	double m_thres_min_y;
	/// the maximum y coordinate of all the points in threshold
	double m_thres_max_y;

	/**
	 * Is true if the calling of the method set_threshold has resulted in
	 * m_thres_set being set.
	 */
	bool m_thres_set;
	/**
	 * Is true if the calling of the method set_threshold has resulted in
	 * m_alpha_thres being set.
	 */
	bool m_alpha_thres_set;

	/// keep a count of how many times the loop crosses the bz reversal threshold
	static long m_count;
};

#endif /* PARAMETERSPACE_H_ */
