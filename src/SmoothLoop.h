/**
 * @file SmoothLoop.h
 *
 * This file declares the CSmoothLoop class.
 */

#ifndef SMOOTHLOOP_H_
#define SMOOTHLOOP_H_

/**
 * The CSmoothLoop class specifies a loop that has a smooth alpha profile that falls to
 * zero at the loop boundary (r1): i.e., the loop has zero net current.
 *
 * The functionality of this class is based on Hood, A. W., et al. 2009, A&A, 506, 913.
 *
 */
class CSmoothLoop {

public:
	CSmoothLoop(void);
	/// empty default destructor
	~CSmoothLoop(void) {}

	void reset(void);

	void init_dimensions(double _r1, double _r2, double _l);
	void init_field(double _lam, double _eta, bool _b1norm, bool _relaxable);
	void init_field(double _lam, double _eta, bool _relaxable);

	double calc_bz(double _r) const;
	double calc_bth(double _r) const;
	double calc_p(double _r) const;
	double calc_alpha(double _r) const;
	double calc_twist(double _r) const;
	double calc_shear(double _r) const;
	double calc_rsntk(double _r) const;
	double calc_psi(double _r) const;
	double calc_k(double _r, double _psi) const;
	double calc_w(double _r) const;

	void profile(unsigned long _rs, bool pf_rx_scenarios = false);
	CRelaxedLoop relax(double _r);
	CRelaxedLoop relax(void);

	
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
	/// true if the relax method has been successfully called for this object since the last reset
	bool relaxed;
	/// always false
	bool netcurrent;
	/// true if it is appropriate for the loop to be relaxed
	bool relaxable;

	/// the loop radius
	double r1;

    /// the envelope radius
	double r2;

	/// the length of the loop in units of r1
	double l;

	/// the magnetic field parameters
	double lam, eta;

    /// the magnetic field coefficients
	double b1, b2;

	/// the magnetic axial field profile
	CProfile pf_bz;
	/// the magnetic azimuthal field profile
	CProfile pf_bth;
	/// the magnetic twist profile ((l*bth)/(r*bz))
	CProfile pf_twist;
	/// the magnetic shear profile (bth/bz)
	CProfile pf_shear;
	/// the alpha profile (j.B/B^2)
	CProfile pf_alpha;
	/// the resonant k profile -(bth/(r*bz))
    CProfile pf_rsntk;


	/**
	 * The rx profiles contain data for a particular relaxation regime.
	 *
     * Threshold (k/psi^2) 0-rx1 is conserved against Relaxed (k/psi^2) 0-rx1, where rx1=[r1:r2]
	 * dW = Threshold (W) 0-rx1 - Relaxed (W) 0-rx1.
	 * The field coefficients of the threshold envelope are retained during relaxation.
	 *
	 */
	/// The profile that plots the relaxation alpha over a range of relaxation radii
	CProfile pf_rx_a;
	/// The profile that plots the relaxation b1 over a range of relaxation radii
	CProfile pf_rx_b1;
	/// The profile that plots the pressure difference (at a specific radial coordinate) over a range of relaxation radii
	CProfile pf_rx_dp;
	/// The profile that plots the average pressure difference (between relaxed loop and envelope) over a range of relaxation radii
	CProfile pf_rx_davp;
	/// The profile that plots the energy of the relaxed state for a range of relaxation radii
	CProfile pf_rx_w;
	/// The profile that plots the energy release for a range of relaxation radii
	CProfile pf_rx_wr;	
	/// The profile that plots the q value for a range of relaxation radii
	CProfile pf_rx_q;


	void init_field_coeffs(bool _b1norm);

	double calc_consrv_err(double _arx) const;
	CRelaxedLoop calc_rx(double _r, double _psi, double _psi_r, double _consrv, double _consrv_r) const;


	static double calc_bz(double _r, void* _param);
	static double calc_bth(double _r, void* _param);
	static double calc_twist(double _r, void* _param);
	static double calc_psi(double _r, void* _param);
	static double calc_k(double _r, void* _param);
	static double calc_consrv_err(double _arx, void* _param);


friend ofstream& operator<<(ofstream& _ofs, const CSmoothLoop& _lp);

}; // end of CLoop class declaration

#endif /* LOOP_H_ */



