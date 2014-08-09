/**
 * @file Loop.h
 *
 * This file declares the CLoop class.
 */

#ifndef LOOP_H_
#define LOOP_H_

/**
 * The CLoop class specifies a loop featuring 4 concentric regions: core, outer layer,
 * neutralisation layer and a potential envelope. Each of these regions has a unique alpha
 * value - the ratio of current density to magnetic field. The potential envelope has an
 * alpha value of zero.
 *
 * The functionality of this class is based on the following academic papers.
 *
 * [1] Browning, P. K., & Van der Linden, R. A. M., 2002, A&A, 400, 355
 * [2] Browning, P. K., et al. 2008, A&A, 485, 837
 *
 */
class CLoop {

public:
	CLoop(void);
	/// empty default destructor
	~CLoop(void) {}

	void reset(void);

	void init_dimensions(double _r1, double _r2, double _r3, double _r4, double _l);
	void init_field(double _a1, double _a2, bool _b1norm, bool _approx_nlalpha, bool _relaxable);
	void init_field(double _a1, double _a2, bool _approx_nlalpha, bool _relaxable);

	double calc_bz(double _r) const;
	double calc_bth(double _r) const;
	double calc_p(double _r) const;
	double calc_twist(double _r) const;
	double calc_baty_twist(double _r) const;
	double calc_shear(double _r) const;
	double calc_rsntk(double _r) const;
	double calc_psi(double _r) const;
	double calc_k(double _r) const;
	double calc_w(double _r) const;

	void profile(unsigned long _rs, bool pf_rx_scenarios = false);
	CRelaxedLoop relax(double _r) const;
	CRelaxedLoop relax(void) const;

	double get_arx(void) const;
	double get_wr(double _wr_limit) const;
	double get_twist(double _r) const;
	double get_twist_av(double _rmin, double _rmax) const;
	double get_baty_twist_av(double _rmin, double _rmax) const;
	double get_velli_twist_av(double _rmin, double _rmax) const;
	double get_shear(double _r) const;
	double get_shear_rms(double _rmin, double _rmax) const;
	double get_shear_mabs(double _rmin, double _rmax) const;

	void out_props(ofstream& _ofs) const;
	void out_field_pfs(ofstream& _ofs) const;
	void out_rx_pfs(ofstream& _ofs) const;

private:
	/// true if the init_dimensions method has been successfully called for this object since the last reset
	bool dimensions_initialised;
	/// true if the init_field_coeffs method has been successfully called for this object since the last reset
	bool coeffs_initialised;
	/// true if all field coefficients have been normalised
	bool coeffs_normalised;
	/// true if the init_field method has been successfully called for this object since the last reset
	bool field_initialised;
	/// true if the profile method has been successfully called for this object since the last reset
	bool profiled;
	/// true if the different relaxation scenarios have been profiled for this object since the last reset
	bool profiled_rx_scenarios;
	/// true if the int_twist and int_shear methods have been successfully called for this object since the last reset
	bool profiles_integrated;
	/// true if the relax method has been successfully called for this object since the last reset
	bool relaxed;
	/// false if there exists a current neutralisation layer, i.e., r3 > r2
	bool netcurrent;
	/// true if it is appropriate for the loop to be relaxed
	bool relaxable;

	/// the core radius
	double r1;
	/// the outer layer radius
	double r2;
	/// the neutralisation layer radius
	double r3;
	/// the potential envelope radius
	double r4;

	/// the length of the loop in units of m_r2
	double l;

	/// the alpha value within the core [0,m_r1]
	double a1;
	/// the alpha value within the outer layer [m_r1,m_r2]
	double a2;
	/// the alpha value within the neutralisation layer [m_r2,m_r3]
	double a3;

	/// the magnetic field coefficient in the core
	double b1;
	/// the magnetic field coefficients in the outer layer
	double b2, c2;
	/// the magnetic field coefficients in the neutralisation layer
	double b3, c3;
	/// the magnetic field coefficients in the potential envelope
	double b4, c4;

	/// the magnetic axial field profile
	CProfile pf_bz;
	/// the magnetic azimuthal field profile
	CProfile pf_bth;
	/// the magnetic twist profile ((l*bth)/(r*bz))
	CProfile pf_twist;
	/// the magnetic shear profile (bth/bz)
	CProfile pf_shear;
	/// the resonant k profile -(bth/(r*bz))
    CProfile pf_rsntk;

	/// the integrals (calculated analytically) of l*bth over the ranges [0,r1], [r1,r2], [r2,r3] and [r3,r4]
	double lbth_0r1, lbth_r1r2, lbth_r2r3, lbth_r3r4;
	/// the integrals (calculated analytically) of r*bz over the ranges [0,r1], [r1,r2], [r2,r3] and [r3,r4]
	double rbz_0r1, rbz_r1r2, rbz_r2r3, rbz_r3r4;
	/// the intergrals (calculated numerically and weighted by area) of the twist over the ranges [0,r1], [r1,r2], [r2,r3] and [r3,r4]
	double twist_0r1, twist_r1r2, twist_r2r3, twist_r3r4;
	/// the intergrals (calculated numerically) of the twist over the ranges [0,r1], [r1,r2], [r2,r3] and [r3,r4] calculated according to Baty 2001
	double baty_twist_0r1, baty_twist_r1r2, baty_twist_r2r3, baty_twist_r3r4;
	/// the intergrals (calculated numerically) of the shear squared over the ranges [0,r1], [r1,r2], [r2,r3] and [r3,r4]
	double shear_sq_0r1, shear_sq_r1r2, shear_sq_r2r3, shear_sq_r3r4;
	/// the intergrals (calculated numerically) of the absolute shear over the ranges [0,r1], [r1,r2], [r2,r3] and [r3,r4]
	double shear_abs_0r1, shear_abs_r1r2, shear_abs_r2r3, shear_abs_r3r4;


	/**
	 * The rx profiles contain data for five relaxation regimes.
	 *
	 * 1. Threshold (k/psi^2) 0-r4 is conserved against Relaxed (k/psi^2) 0-r4.
	 *    dW = Threshold (W) 0-r4 - Relaxed (W) 0-r4.
	 *    Field between loop and envelope of the relaxed state is continuous.
	 *
	 *
     * 2. Threshold (k/psi^2) 0-r4 is conserved against Relaxed (k/psi^2) 0-r4.
	 *    dW = Threshold (W) 0-r4 - Relaxed (W) 0-r4.
	 *    The field coefficients of the threshold envelope are retained during relaxation.
     *
     * 2.1. Same as 2. except Threshold (k/psi^2) 0-rx1 is conserved against Relaxed (k/psi^2) 0-rx1,
     *      where rx1 [r2:r4]
	 *      dW = Threshold (W) 0-rx1 - Relaxed (W) 0-rx1
	 *
	 *
     * 3. Threshold (k/psi^2) 0-r4 is conserved against Relaxed (k/psi^2) 0-r4.
	 *    dW = Threshold (W) 0-r4 - Relaxed (W) 0-r4.
	 *    Relaxed state features a neutralisation layer that ensures relaxed envelope has no azimuthal field.
	 *    The field coefficients of the threshold envelope are retained during relaxation.
	 *
	 * 3.1. Same as 3. except Threshold (k/psi^2) 0-rx1 is conserved against Relaxed (k/psi^2) 0-rx1,
     *      where rx1 [r2:r4]
	 *      dW = Threshold (W) 0-rx1 - Relaxed (W) 0-rx1
	 *
	 */
	/// The profiles that plot the relaxation alpha over a range of relaxation radii
	vector<CProfile> vpf_rx_a;
	/// The profiles that plot the relaxation b1 over a range of relaxation radii
	vector<CProfile> vpf_rx_b1;
	/// The profiles that plot the pressure difference (at a specific radial coordinate) over a range of relaxation radii
	vector<CProfile> vpf_rx_dp;
	/// The profiles that plot the average pressure difference (between relaxed loop and envelope) over a range of relaxation radii
	vector<CProfile> vpf_rx_davp;
	/// The profiles that plot the energy release for a range of relaxation radii
	vector<CProfile> vpf_rx_wr;
	/// The profiles that plot the q value for a range of relaxation radii
	vector<CProfile> vpf_rx_q;


	void init_field_coeffs(bool _b1norm, bool _approx_nlalpha);

	double calc_nlalpha(void);
	double approx_nlalpha(void);

	double calc_consrv_err(double _arx) const;
	CRelaxedLoop calc_rx(double _r,  double _rnl, double _psi, double _psi_r, double _consrv, double _consrv_r) const;

	void int_twist(void* _param);
	void int_baty_twist(void* _param);
	void int_velli_twist(void);
	void int_shear_sq(void* _param);
	void int_shear_abs(void* _param);

	static double calc_bz(double _r, void* _param);
	static double calc_bth(double _r, void* _param);
	static double calc_r3bth(double _a3, void* _param);
	static double calc_twist(double _r, void* _param);
	static double calc_baty_twist(double _r, void* _param);
	static double calc_shear_sq(double _r, void* _param);
	static double calc_shear_abs(double _r, void* _param);
	static double calc_consrv_err(double _arx, void* _param);

	double J0(double _param) const;
	double J1(double _param) const;
	double Y0(double _param) const;
	double Y1(double _param) const;
	double F0(double _param) const;
	double F1(double _param) const;
	double G0(double _param) const;
	double G1(double _param) const;

friend ofstream& operator<<(ofstream& _ofs, const CLoop& _lp);

}; // end of CLoop class declaration

#endif /* LOOP_H_ */



