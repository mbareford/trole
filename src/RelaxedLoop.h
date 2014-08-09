/**
 * @file RelaxedLoop.h
 *
 * This file declares the CRelaxedLoop class.
 */

#ifndef RELAXEDLOOP_H_
#define RELAXEDLOOP_H_

/**
 * The CRelaxedLoop class specifies a loop featuring at 2 concentric regions: core and
 * a neutralisation layer. Each of these regions has a unique alpha value - the ratio of current
 * density to magnetic field. The loop may also be placed within a potential envelope where the
 * alpha value is zero.
 *
 */
class CRelaxedLoop {

public:
	CRelaxedLoop(void);
	CRelaxedLoop(const CRelaxedLoop& _rx);
	/// empty default destructor
	~CRelaxedLoop(void) {}

	CRelaxedLoop& operator=(const CRelaxedLoop& _rx);

	void reset(bool _full = false);

	void init_dimensions(double _r1, double _r2, double _l);
	void init_dimensions(double _r1, double _r2, double _r3, double _l);
	void init_field(double _psi, double _psi_r, double _a1, bool _helical_test = true);
	void fix_envelope(double _b, double _c = 0.0);
	void neutralise_envelope(void);

	double calc_bz(double _r) const;
	double calc_bth(double _r) const;
	double calc_p(double _r) const;
	double calc_p_av(double _rmin, double _rmax) const;
	double calc_psi(double _r) const;
	double calc_k(double _r) const;
	double calc_w(double _r) const;
	double calc_q(double _r) const;

	void profile(unsigned long _rs);

	double get_a1(void) const;
	double get_b1(void) const;
	double get_r1(void) const;
	double get_r2(void) const;

	void out_props(ofstream& _ofs) const;
	void out_pfs(ofstream& _ofs) const;

	static void dump_pfs(double _r, const CRelaxedLoop& _rx);

	/**
	 * Critical helicity is achieved for a constant alpha region when the value of the
	 * Bessel arguments is greater than or equal to KCRIT_BESSEL_ARG.
	 */
	static const double KCRIT_BESSEL_ARG;

private:
	/// true if the init_dimensions method has been successfully called for this object since the last reset
	bool dimensions_initialised;
	/// true if the field coefficients for the potential envelope have been fixed
	bool envelope_fixed;
	/// true if the C coefficients for the potential envelope has been set to zero
	bool envelope_neutralised;
	/// true if the init_field_coeffs method has been successfully called for this object since the last reset
	bool coeffs_initialised;
	/// true if all field coefficients have been normalised such that axial flux is conserved
	bool coeffs_normalised;
	/// true if the init_field method has been successfully called for this object since the last reset
	bool field_initialised;
	/// true if the profile method has been successfully called for this object since the last reset
	bool profiled;

	/// true if the relaxed loop is placed within a potential envelope
	bool has_pe;
	/// true if the relaxed loop has a current neutralisation layer
	bool has_nl;
	/// true if the core region has a critical helicity i.e., the relaxed state features helical modes
	bool is_helical;

	/// the core radius
	double r1;
	/// the neutralisation layer radius (or potential envelope radius if neutralisation layer does not exist)
	double r2;
	/// the potential envelope radius (if neutralisation layer exists)
	double r3;

	/// the length of the loop in units of its unrelaxed radius
	double l;

	/// the alpha value within the core [0,m_r1]
	double a1;
	/// the alpha value within the neutralisation layer (if it exists) [m_r2,m_r3]
	double a2;

	/// the magnetic field coefficient in the core
	double b1;
	/// the magnetic field coefficients in the neutralisation layer (or envelope radius if neutralisation layer does not exist)
	double b2, c2;
	/// the magnetic field coefficients in the envelope (if neutralisation layer exists)
	double b3, c3;

	/// the fixed magnetic field coefficients for the potential envelope
	double bpe, cpe;

	/// the magnetic axial field profile
	CProfile pf_bz;
	/// the magnetic azimuthal field profile
	CProfile pf_bth;
	/// the magnetic twist profile ((l*bth)/(r*bz))
	CProfile pf_twist;
	/// the magnetic shear profile (bth/bz)
	CProfile pf_shear;

	// the dump_id is always incremented after every dump
	static unsigned int dump_id;


	void init_field_coeffs(double _psi, double _psi_r);

	double calc_nlalpha(void);

	bool validate_radius(double _r) const;

	static double calc_bz(double _r, void* _param);
	static double calc_bth(double _r, void* _param);
	static double calc_r2bth(double _an, void* _param);

	double J0(double _param) const;
	double J1(double _param) const;
	double Y0(double _param) const;
	double Y1(double _param) const;
	double F0(double _param) const;
	double F1(double _param) const;

friend ofstream& operator<<(ofstream& _ofs, const CRelaxedLoop& _lp);

}; // end of CRelaxedLoop class declaration

#endif /* RELAXEDLOOP_H_ */



