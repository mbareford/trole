/**
 * @file Profile.h
 *
 * This file declares the CProfile class.
 */

#ifndef PROFILE_H_
#define PROFILE_H_

/**
 * This class permits 2 dimensional profiles to be created - it is
 * possible to specify the resolution (or granularity) of these profiles.
 *
 * The minimum and maximum coordinates of the points within the profile
 * are also recorded.
 */
class CProfile {

public:
	CProfile(void);
	CProfile(const CProfile& _pf);
	~CProfile(void);

	CProfile& operator=(const CProfile& _pf);

	void reset(void);

	void init(double _min_x, double _max_x, unsigned long _res);
	bool is_initialised(void) { return m_initialised; }

	double plot(double _y);
	double plot(double _y_num, double _y_den);
	unsigned long plot_count(double _x);

	double get_plot_pos(void) const;
	double get_plot_val(double _x) const;

	double get_plot_pos(unsigned long _i) const;
	double get_plot_val(unsigned long _i) const;
	unsigned long get_plot_cnt(void) const { return m_pt_count; }

	bool done(void) const { return m_initialised && m_pt_i == m_pt_count; }
	bool done(double _x) const { return m_initialised && (m_pt_i == m_pt_count || _x < m_x[m_pt_i]); }

	double get_min_plot_val(double _min_x, double _max_x) const;
	double get_max_plot_val(double _min_x, double _max_x) const;

	double get_max_diff(CProfile* _ppf, bool _weighted, double _max_pos, double& _max_diff_pos);

private:
	/// initialisation flag - if true the init method has been successfully called for this object
	bool m_initialised;

	/// the minimum x coordinate
	double m_min_x;
	/// the maximum x coordinate
	double m_max_x;

	/// the total number of points within the profile
	unsigned long m_pt_count;
	/// the index of the current point
	unsigned long m_pt_i;

	/// pointer to array of x coordinates
	vector<double> m_x;
	/// pointer to array of y coordinates
	vector<double> m_y;

}; // end of CProfile class declaration


#endif /* PROFILE_H_ */
