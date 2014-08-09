/**
 * @file SmoothLoop.cpp
 *
 * This file defines the methods of the CSmoothLoop class.
 */

#include "StdInclude.h"

#include "Exception.h"
#include "GSL.h"
#include "NetBthCalculator.h"
#include "Profile.h"
#include "RelaxedLoop.h"
#include "RelaxationCalculator.h"
#include "SmoothLoop.h"


const double PI = CGSL::PI;


/**
 * All attributes are zeroed.
 */
CSmoothLoop::CSmoothLoop(void): dimensions_initialised(false), coeffs_initialised(false), coeffs_normalised(false), field_initialised(false),
                    profiled(false), profiled_rx_scenarios(false), relaxed(false),
                    netcurrent(false), relaxable(false) {
	reset();
}

/**
 * Set attributes to zero.
 */
void CSmoothLoop::reset(void) {
	if (!dimensions_initialised) {
		r1 = r2 = 0.0;
		l = 0.0;
	}

	if (!field_initialised) {
		lam = eta = 0.0;
		b1 = b2 = 0.0;
	}

	if (!profiled) {
		pf_bz.reset();
		pf_bth.reset();
		pf_twist.reset();
		pf_shear.reset();
		pf_alpha.reset();
		pf_rsntk.reset();
	}	
}

/**
 * Initialises the loop's dimensions - the radial boundaries and length.
 *
 * Other attributes are zeroed.
 *
 * If this function completes successfully the dimensions initialisation flag is
 * set to true.
 *
 * @param _r1 the loop radius.
 * @param _r2 the envelope radius.
 * @param _l the length of the loop in units of _r1.
 * @throw CException should any of the following occur:
 *        1. any of the radial boundaries has a zero value;
 *        2. radial boundaries in incorrect numerical order;
 *        3. length is not positive.
 */
void CSmoothLoop::init_dimensions(double _r1, double _r2, double _l) {
	if (0.0 == _r1 || 0.0 == _r2) {
		throw CException("Error! Zero radial boundary.");
	}

	if (_r1 > _r2) {
		throw CException("Error! Radial boundary out of order (must have 0 < r1 <= r2).");
	}

	if (_l <= 0.0) {
		throw CException("Error! Length is not positive.");
	}

	dimensions_initialised = false;
	coeffs_initialised = false;
	coeffs_normalised = false;
	field_initialised = false;
	profiled = false;
	relaxed = false;

	relaxable = false;

	reset();

	r1 = _r1;
	r2 = _r2;
	
	l = _l;

	dimensions_initialised = true;
}

/**
 * Calculates the magnetic field coefficient b1.
 *
 * This function is private, it is called by init_field.
 *
 * @param _b1norm if true b1 is normalised to 1.0, otherwise psi_0r2 is set to 1.0.
 */
void CSmoothLoop::init_field_coeffs(bool _b1norm) {
	coeffs_initialised = false;
	coeffs_normalised = false;

	// assume b1 is normalised
	b1 = 1.0;

	// so far all field coefficients are in units of b1
	coeffs_initialised = true;

	// if the total axial flux is normalised rescale all magnetic field coefficients
	if (!_b1norm) {
		void* param = this;

		b1 = 1.0/CGSL::integrate(CSmoothLoop::calc_psi, 0.0, r2, param);
	}

    double e1 = 2.0*eta;
    double e2 = e1 + 1.0;
    double e3 = e1/e2;
    b2 = sqrt(b1*b1 - e3*pow(r1,e2));

	coeffs_normalised = true;
}

/**
 * Initialises the loop object's lambda and eta values.
 *
 * If this function completes successfully the field initialisation flag is
 * set to true.
 *
 * @param _lam the lambda parameter used to define a smooth alpha profile.
 * @param _eta the eta parameter used to define a smooth alpha profile.
 * @param _b1norm if true b1 is normalised to 1.0, otherwise psi_0r2 is set to 1.0.
 * @param _relaxable is true if it is appropriate for this loop to be relaxed.
 * @throw CException should any of the following occur:
 *        1. loop dimensions uninitialised.
 */
void CSmoothLoop::init_field(double _lam, double _eta, bool _b1norm, bool _relaxable) {
	if (!dimensions_initialised) {
		throw CException("Error! Loop dimensions uninitialised.");
	}

	coeffs_initialised = false;
	coeffs_normalised = false;
	field_initialised = false;
	profiled = false;
	relaxed = false;

	relaxable = _relaxable;

	reset();

	lam = _lam;
	eta = _eta;

	init_field_coeffs(_b1norm);

	field_initialised = true;
}

/**
 * This function permits the field to be reinitialised. The normalisation
 * is kept the same.
 *
 * @param _lam the lambda parameter used to define a smooth alpha profile.
 * @param _eta the eta parameter used to define a smooth alpha profile.
 * @param _relaxable is true if it is appropriate for this loop to be relaxed.
 * @throw CException should any of the following occur:
 *        1. loop dimensions uninitialised.
 */
void CSmoothLoop::init_field(double _lam, double _eta, bool _relaxable) {
	if (!field_initialised) {
		throw CException("Error! Loop field uninitialised.");
	}

	init_field(_lam, _eta, (1.0 == b1), _relaxable);
}


/**
 * Calculate axial component of the magnetic field for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @return the axial component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 */
double CSmoothLoop::calc_bz(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r2) {
		throw CException("Error! Radial value out of bounds.");
	}

	double bz = 0.0;

	if (0.0 == r) {
		// on axis or potential loop or (within core and potential core)
		bz = b1;
	}
	else {
	  double l1 = lam*lam;
	  double e1 = 2.0*eta;
	  double e2 = e1 + 1.0;
	  double e3 = l1/e2;
		
	  if (r <= r1) {
	      // within core		
		  bz = sqrt(e3*pow(r1-r*r,e2) - l1*r*r*pow(r1-r*r,e1) + b1*b1 - e3*pow(r1,e2));
	  }
	  else {
		  // within potential envelope
		  bz = sqrt(b1*b1 - e3*pow(r1,e2));
	  }
	  
	}

	return bz;
}

/**
 * Calculate axial component of the magnetic field for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @param _param points to a loop object
 * @return the axial component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null loop object.
 */
double CSmoothLoop::calc_bz(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CSmoothLoop* loop = (CSmoothLoop*) (_param);
	if (NULL == loop) {
		throw CException("Error! NULL loop object.");
	}

	return loop->calc_bz(_r);
}

/**
 * Calculate azimuthal component of the magnetic field for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @return the azimuthal component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 */
double CSmoothLoop::calc_bth(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r2) {
		throw CException("Error! Radial value out of bounds.");
	}

	double bth = 0.0;

	if (0.0 == r) {
		// on axis
		bth = 0.0;
	}
	else if (r <= r1) {
		// within loop
		bth = lam*r*pow(r1-r*r, eta);
	}
	else {
		// within potential envelope
		bth = 0.0;
	}

	return bth;
}

/**
 * Calculate azimuthal component of the magnetic field for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @param _param points to a loop object
 * @return the azimuthal component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null loop object.
 */
double CSmoothLoop::calc_bth(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CSmoothLoop* loop = (CSmoothLoop*) (_param);
	if (NULL == loop) {
		throw CException("Error! NULL loop object.");
	}

	return loop->calc_bth(_r);
}

/**
 * Get the pressure of the magnetic field at the specified radial coordinate.
 *
 * @param _r the radial coordinate
 * @return the magnetic pressure at _r
 */
double CSmoothLoop::calc_p(double _r) const {
	double bz = calc_bz(_r);
	double bth = calc_bth(_r);

	return (bz*bz + bth*bth)/2.0;
}

/**
 * This is the function used by to calculate the twist at a given radius.
 *
 * @param _r the radial coordinate
 * @return the twist at _r
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 * @throw CSingularityException should any of the following occur:
 *        1. zero bz.
 */
double CSmoothLoop::calc_twist(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r2) {
		throw CException("Error! CLoop::calc_twist(), radial value out of bounds.");
	}

	if (0.0 == _r) {
		return l*calc_alpha(0.0)/2.0;
	}

	double tw = 0.0;

	if (r > 0.0 ) {
		double bz = calc_bz(r);
		double bth = calc_bth(r);

		if (0.0 == bz) {
			throw CSingularityException(r, "Error! CSmoothLoop::calc_twist(), zero bz.");
		}
		else {
			tw = (l*bth) / (r*bz);
		}
	}

	return tw;
}

/**
 * This is the function used by to calculate the shear at a given radius.
 *
 * @param _r the radial coordinate
 * @return the shear at _r
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 * @throw CSingularityException should any of the following occur:
 *        1. zero bz.
 */
double CSmoothLoop::calc_shear(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r2) {
		throw CException("Error! CSmoothLoop::calc_shear(), radial value out of bounds.");
	}

	double sh = 0.0;
	double bz = calc_bz(r);
	double bth = calc_bth(r);

	if (0.0 == bz) {
		throw CSingularityException(r, "Error! CSmoothLoop::calc_shear(), zero bz.");
	}
	else {
		sh = bth/bz;
	}

	return sh;
}

/**
 * Return the value of alpha (j.b/b^2) at radial coordinate _r.
 *
 * @param _r the radial coordinate that is the upper bound for the energy calculation
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. radial value out of bounds.
 */
double CSmoothLoop::calc_alpha(double _r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (_r < 0.0 || _r > r2) {
		throw CException("Error! Radial value out of bounds.");
	}

	if (0.0 == _r) {
		return 2.0*lam/b1;
	}

	if (0.0 == lam) {
		// potential loop
		return 0.0;
	}

	double alpha = 0.0;

	do {
		double r = _r < r1 ? _r : r1;
		
		// loop
   	    alpha = 2.0*lam*pow(r1-r*r,eta-1.0)*(r1 - r*r*(1.0+eta))/calc_bz(r);

		if (_r <= r1) {
			break;
		}

		r = _r < r2 ? _r : r2;

		// potential envelope
     	alpha = 0.0;
	
	} while (0);

	return alpha;
}

/**
 * This is the function used by to calculate the resonant k value at a given radius.
 *
 * @param _r the radial coordinate
 * @return the resonant k value at _r
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 * @throw CSingularityException should any of the following occur:
 *        1. zero bz.
 */
double CSmoothLoop::calc_rsntk(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r2) {
	    throw CException("Error! CSmoothLoop::calc_rsntk(), radial value out of bounds.");
	}

	if (0.0 == _r) {
		return -calc_alpha(0.0)/2.0;
	}

	double rsntk = 0.0;

    if (r > 0.0) {
	    double bz = calc_bz(r);
	    double bth = calc_bth(r);

		if (0.0 == bz) {
	        throw CSingularityException(r, "Error! CSmoothLoop::calc_rsntk(), zero bz.");
		}
		else {
	        rsntk = -bth/(r*bz);
		}
	}

	return rsntk;
}

/**
 * Calculate the axial magnetic flux for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @return the azimuthal component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 */
double CSmoothLoop::calc_psi(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r2) {
		throw CException("Error! CSmoothLoop::calc_psi(), radial value out of bounds.");
	}	

	return 2.0*PI*r*calc_bz(r);
}

/**
 * Calculate the axial magnetic flux for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @param _param points to a loop object
 * @return the azimuthal component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null loop object.
 */
double CSmoothLoop::calc_psi(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CSmoothLoop* loop = (CSmoothLoop*) (_param);
	if (NULL == loop) {
		throw CException("Error! NULL loop object.");
	}	

	return loop->calc_psi(_r);
}

/**
 * Calculate the magnetic helicity for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @param _psi the axial flux over 0-_r

 * @return the azimuthal component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 */
double CSmoothLoop::calc_k(double _r, double _psi) const {
	double r = _r;
	if (r < 0.0 || r > r2) {
		throw CException("Error! CSmoothLoop::calc_k(), radial value out of bounds.");
	}
    
	return 2.0*l*calc_bth(r)*_psi;
}

/**
 * Calculate the magnetic helicity for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @param _param points to a loop object
 * @return the azimuthal component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null loop object.
 */
double CSmoothLoop::calc_k(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CSmoothLoop* loop = (CSmoothLoop*) (_param);
	if (NULL == loop) {
		throw CException("Error! NULL loop object.");
	}

	return loop->calc_k(_r, CGSL::integrate(CSmoothLoop::calc_psi, 0.0, _r, _param));
}

/**
 * Return the energy between 0 and _r.
 *
 * @param _r the radial coordinate that is the upper bound for the energy calculation
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. radial value out of bounds.
 */
double CSmoothLoop::calc_w(double _r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (_r < 0.0 || _r > r2) {
		throw CException("Error! Radial value out of bounds.");
	}

	if (0.0 == _r) {
		return 0.0;
	}

	if (0.0 == lam) {
		return (l*PI*b1*b1*_r*_r)/2.0;
	}

	double w = 0.0;

	do {
		if (0.0 == lam) {
			// potential loop
			w = (b1*b1*_r*_r)/2.0;
			break;
		}

		double r = _r < r1 ? _r : r1;
		
		double l1 = lam*lam;
	    double e1 = 2.0*eta + 1.0;
	    double e2 = e1 + 1.0;  

		// loop
   	    w = (0.5*l1)*((pow(r1,e1)/e1 - (b1*b1)/l1)*(r1-r*r) - pow(r1-r*r,e2)/(e1*e2) - (pow(r1,e2)/e2 - (b1*b1/l1)*r1));		

		if (_r <= r1) {
			break;
		}

		r = _r < r2 ? _r : r2;

		// potential envelope
     	w += 0.5*(b1*b1 - (l1/e1)*pow(r1,e1))*(r*r - r1*r1);
	
	} while (0);

	return l*PI*w;
}

/**
 * When a loop undergoes a relaxation the quantity k/(psi^2) is conserved.
 * This function is used to find the value of the relaxation alpha that
 * conserves this quantity. It does this by returning the conservation error:
 * the difference in k/(psi^2) between the threshold and relaxed states.
 * If this quantity is conserved the difference should be zero.
 *
 * This function is used in conjunction with GSL root solver gsl_root_fsolver_brent.
 *
 * @param _a1 the prospective relaxed alpha value.
 * @param _param points to a CRelaxationCalculator object
 * @return zero if _rxa really does represent the relaxed state.
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null CRelaxationCalculator object.
 */
double CSmoothLoop::calc_consrv_err(double _a1, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CRelaxationCalculator* calc = (CRelaxationCalculator*) (_param);
	if (NULL == calc) {
		throw CException("Error! NULL CRelaxationCalculator object.");
	}

	return calc->consrv_err(_a1);
}

/**
 * Calculate the relaxed state i.e., the minimum magnetic energy of the loop
 * that conserves k/(psi*psi) and set the rx attribute to this state.
 *
 * @param _r relaxation radius
 * @param _psi the axial flux that is conserved by the relaxed loop
 * @param _psi_r the radial range of the axial flux that is conserved
 * @param _consrv the value k/(psi*psi) that is conserved by the relaxed loop
 * @param _consrv_r the radial range of the threshold state over which _consrv is calculated
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. invalid relaxation radius;
 * @throw CLoopException should any of the following occur:
 *        1. error finding a1 (relaxation).
 */
CRelaxedLoop CSmoothLoop::calc_rx(double _r, double _psi, double _psi_r, double _consrv, double _consrv_r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (_r <= 0.0 || _r > r2) {
		throw CException("Error! Invalid relaxation radius.");
	}

	CRelaxedLoop rx;

	// initialise the dimensions of the relaxed loop
	rx.init_dimensions(_r, r2, l);
	
	// initialise the relaxed state's potential envelope so that it has the
	// same axial field as that of the envelope for the threshold state
	rx.fix_envelope(b2);	

	// calculate relaxed state according to conservation parameters
	double rxa = 0.0;

	CRelaxationCalculator calc(_psi, _psi_r, _consrv, _consrv_r, rx);

	void* param = &calc;

	try {
		rxa = CGSL::find_root(CSmoothLoop::calc_consrv_err, 50.0, 0.01, 0, param);
	}
	catch (CGSLException& x) {
		throw CSmoothLoopException(lam, eta, "error finding a1 (relaxation)", x);
	}

	rx.init_field(_psi, _psi_r, rxa);

	return rx;
}


/**
 * Calculates four radial profiles - magnetic axial field, magnetic azimuthal field,
 * magnetic twist ((l*bth)/(r*bz)) and shear (bth/bz). Also creates the rx, dp and wr
 * profiles for different relaxation states.
 *
 * If this function completes successfully the profiled flag is set to true.
 *
 * @param _rs the profile resolution i.e., the number of points that make up the profile.
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. divide by zero error.
 */
void CSmoothLoop::profile(unsigned long _rs, bool pf_rx_scenarios/*=false*/) {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	profiled = false;
	profiled_rx_scenarios = false;
	
	reset();

	// construct field profiles
	//////////////////////////////////////////////////////////////////////////
	double r, bz, bth, tw, sh, alpha, rsntk;
	
	pf_bz.init(0.0, r2, _rs);
	pf_bth.init(0.0, r2, _rs);
	pf_twist.init(0.0, r2, _rs);
	pf_shear.init(0.0, r2, _rs);
	pf_alpha.init(0.0, r2, _rs);
	pf_rsntk.init(0.0, r2, _rs);

	do {
		r = pf_bz.get_plot_pos();

		// calculate axial and azimuthal magnetic components
		bz = calc_bz(r);
		bth = calc_bth(r);
		// calculate twist and shear
		tw = calc_twist(r);
		sh = calc_shear(r);
        // calculate j.B/B^2
        alpha = calc_alpha(r);
        // calculate the resonant k value for the given radius
        rsntk = calc_rsntk(r);

		pf_bz.plot(bz);
		pf_bth.plot(bth);
		pf_twist.plot(tw);
		pf_shear.plot(sh);
		pf_alpha.plot(alpha);
        pf_rsntk.plot(rsntk);

	} while (!pf_bz.done());
	//////////////////////////////////////////////////////////////////////////

	profiled = true;

	/////////////////////////////////////////////////////////////////////////////////////////////
	if (pf_rx_scenarios) {
		// construct arx, dp and wr profiles for the relaxation regimes
		/////////////////////////////////////////////////////////////////////////////////////////////
		
		// initialise the profiles for the relaxation scenarios
		double rmin = r1, rmax = r2, rs;

		rs = ((rmax-rmin)/rmax)*_rs;

		pf_rx_a.init(rmin, rmax, rs);
		pf_rx_b1.init(rmin, rmax, rs);
		pf_rx_dp.init(rmin, rmax, rs);
		pf_rx_davp.init(rmin, rmax, rs);
		pf_rx_w.init(rmin, rmax, rs);
		pf_rx_wr.init(rmin, rmax, rs);
		pf_rx_q.init(rmin, rmax, rs);
		

		// plot the profiles for the relaxation scenarios
		double psi_0r = 0.0, k_0r = 0.0;

		CRelaxedLoop rx;

		void* param = this;

		do {			
			r = pf_rx_a.get_plot_pos();
			psi_0r = CGSL::integrate(CSmoothLoop::calc_psi, 0.0, r, param);
			k_0r = CGSL::integrate(CSmoothLoop::calc_k, 0.0, r, param);
			
			// Scenario 2.1: same as Scenario 2, except conservation follows expansion radius, see Loop.cpp.
			////////////////////////////////////////////////////////////////////////////////////////////////
			rx = calc_rx(r, psi_0r, r, k_0r/(psi_0r*psi_0r), r);
			rx.profile(_rs);
			//CRelaxedLoop::dump_pfs(r, rx);

			pf_rx_a.plot(rx.get_a1());
			pf_rx_b1.plot(rx.get_b1());
			pf_rx_dp.plot(rx.calc_p(r) - calc_p(r));
			pf_rx_davp.plot(rx.calc_p_av(0.0, r) - rx.calc_p_av(r, r2));
			pf_rx_w.plot(rx.calc_w(r));
			pf_rx_wr.plot(calc_w(r) - rx.calc_w(r));
			pf_rx_q.plot(rx.calc_q(r));
			////////////////////////////////////////////////////////////////////////////////////////////////

		} while (!pf_rx_a.done());

		profiled_rx_scenarios = true;
	} // end of <if (pf_rx_scenarios)> clause
	/////////////////////////////////////////////////////////////////////////////////////////////		
}


/**
 * Return the relaxed state of this loop as a CRelaxedLoop object.
 * Relaxation is performed according to Scenario 2.1.
 *
 * @param _r the relaxation radius; is usually greater than the theshold loop radius
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. invalid relaxation radius.
 */
CRelaxedLoop CSmoothLoop::relax(double _r) {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (_r < r1 || _r > r2) {
		throw CException("Error! Invalid relaxation radius.");
	}

	void* param = this;

	double psi_0r = CGSL::integrate(CSmoothLoop::calc_psi,0,_r,param);
	double k_0r = CGSL::integrate(CSmoothLoop::calc_k,0,_r,param);

	// relax loop according to Scenario 2.1, see Loop.cpp
	return calc_rx(_r, psi_0r, _r, k_0r/(psi_0r*psi_0r), _r);
}

/**
 * Return the relaxed state of this loop as a CRelaxedLoop object.
 * Relaxation is performed according to Scenario 2.1, see Loop.cpp.
 *
 * The radius of the relaxed state is chosen such the loop does
 * not relax beyond q=n. The relaxation radius of the loop is
 * where q reaches its resonant value.
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 */
CRelaxedLoop CSmoothLoop::relax(void) {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	const double R_MIN = r1;
	const double R_MAX = r2;
	const double R_INC = 0.01;

	double r = R_MIN, r_rx = r1;
	CRelaxedLoop rx_loop = relax(r);
	double q = rx_loop.calc_q(r);
	int first_q_int = abs((int) q);

	double prev_r = r;
	double prev_q = q;
	int q_int = 0;

	// find the relaxation radius that is closest to resonance
	r += R_INC;

	while (r <= R_MAX) {
		try {
			rx_loop = relax(r);
			q = rx_loop.calc_q(r);
		}
		catch (CLoopException& x) {
			// can't find relaxed state, so move to next radius
			r += R_INC;
			continue;
		}
		q_int = abs((int) q);

		if (q_int == first_q_int + 1) {
			// relaxation radius has crossed a q=n layer, where n is a whole number
			r_rx = prev_r + (fabs(q_int - prev_q)/fabs(q - prev_q))*R_INC;
			break;
		}

		prev_r = r;
		prev_q = q;

		r += R_INC;
	} // end of <while (r <= R_MAX)> loop

	if (r > R_MAX || r_rx > R_MAX) {
		// unable to find relaxation radius
		// so assume no expansion
		r_rx = r1;
	}

	// relax loop according to Scenario 2.1, see Loop.cpp
	return relax(r_rx);
}


/**
 * Stream the loop's properties to a file
 * (profile date not included).
 *
 * @param _ofs the output file stream.
 */
void CSmoothLoop::out_props(ofstream& _ofs) const {
	void* param = (CSmoothLoop*) this;

	if (dimensions_initialised && field_initialised) {
		_ofs << lam << WS;                                     //  1. lambda
		_ofs << eta << WS;                                     //  2. eta

		_ofs << b1 << WS;                                      //  3. coefficient term for core field
        _ofs << b2 << WS;                                      //  4. coefficient term for envelope field

		_ofs << CGSL::integrate(calc_psi,0,r1,param) << WS;    //  5. axial flux in the loop
		_ofs << CGSL::integrate(calc_psi,0,r2,param) << WS;    //  6. axial flux in the loop and envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS;
	}

	
	if (dimensions_initialised && field_initialised) {
		_ofs << CGSL::integrate(calc_k,0,r1,param) << WS;      //  7. magnetic helicity in loop
		_ofs << CGSL::integrate(calc_k,0,r2,param) << WS;      //  8. magnetic helicity in loop and envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS;
	}

	if (dimensions_initialised && field_initialised) {
		_ofs << calc_w(r1) << WS;                              //  9. magnetic energy in loop
		_ofs << calc_w(r2) << WS;                              // 10. magnetic energy in loop and envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS;
	}

	
    if (profiled) {
		int rsntk_min_mlt = ceil(0.5*l*pf_rsntk.get_min_plot_val(0.0,r2)/PI);
		int rsntk_max_mlt = floor(0.5*l*pf_rsntk.get_max_plot_val(0.0,r2)/PI);
        double rsntk_scr = 0.0;
		
        if (rsntk_max_mlt >= rsntk_min_mlt) {
            unsigned long plot_cnt = pf_rsntk.get_plot_cnt();
			unsigned long i = 0;

            int prv_rsntk = floor(0.5*l*pf_rsntk.get_plot_val((unsigned long) 0)/PI);
            double rsntk = 0.0;

	        for (i = 1; i < plot_cnt; ++i) {
				rsntk = floor(0.5*l*pf_rsntk.get_plot_val(i)/PI);
                if (rsntk != prv_rsntk) {                      
					rsntk_scr += (rsntk > prv_rsntk ? rsntk*pf_rsntk.get_plot_pos(i) : prv_rsntk*pf_rsntk.get_plot_pos(i-1));
                }
                prv_rsntk = rsntk;
            }
	    }
        else {
			rsntk_min_mlt = 0;
 	        rsntk_max_mlt = 0;
        }

		_ofs << rsntk_min_mlt << WS;   // 11. the minimum multiple of 2PI/L within the resonant k profile between 0 and r2.
		_ofs << rsntk_max_mlt << WS;   // 12. the maximum multiple of 2PI/L
		_ofs << rsntk_scr << WS;       // 13. radially weighted resonance score
	}
	else {
		_ofs << 0 << WS << 0 << WS << 0 << WS;
	}


	_ofs << NL;
}

/**
 * Stream the loop's field profile data to a file.
 *
 * @param _ofs the output file stream.
 * @throw CException should any of the following occur:
 *        1. loop not profiled.
 */
void CSmoothLoop::out_field_pfs(ofstream& _ofs) const {
	if (!profiled) {
		throw CException("Error! Loop not profiled.");
	}

	unsigned long plot_cnt = pf_bz.get_plot_cnt();
	unsigned long i = 0;

	for (i = 0; i < plot_cnt; ++i) {
		_ofs << pf_bz.get_plot_pos(i) << WS;    // 1. radial position
		_ofs << pf_bz.get_plot_val(i) << WS;    // 2. axial field
		_ofs << pf_bth.get_plot_val(i) << WS;   // 3. azimuthal field
		_ofs << pf_twist.get_plot_val(i) << WS; // 4. twist
		_ofs << pf_shear.get_plot_val(i) << WS; // 5. shear
		_ofs << pf_alpha.get_plot_val(i) << WS; // 6. alpha (j.B/B^2)
        _ofs << pf_rsntk.get_plot_val(i) << WS; // 7. resonant k

		_ofs << NL;
	}
}

/**
 * Stream the loop's relaxation profile data to a file.
 *
 * @param _ofs the output file stream.
 * @throw CException should any of the following occur:
 *        1. loop not profiled.
 */
void CSmoothLoop::out_rx_pfs(ofstream& _ofs) const {
	if (!profiled) {
		throw CException("Error! Loop not profiled.");
	}

	unsigned long plot_cnt = pf_rx_a.get_plot_cnt();
	unsigned long i = 0;

	for (i = 0; i < plot_cnt; ++i) {
		double rrx = pf_rx_a.get_plot_pos(i);
		_ofs << rrx << WS;                               // 1. relaxation radius

		double arx = pf_rx_a.get_plot_val(i);
		_ofs << arx << WS;                               // 2. relaxation alpha for when k/(psi*psi) is conserved

		if (CRelaxedLoop::KCRIT_BESSEL_ARG/rrx == arx) {
			_ofs << '0' << WS;                           // 3. relaxed state features helical modes
		}
		else {
			_ofs << '1' << WS;			                 // 3. relaxed state is axisymmetric
		}

		_ofs << pf_rx_b1.get_plot_val(i) << WS;          // 4. b1 coefficient
		_ofs << pf_rx_dp.get_plot_val(i) << WS;          // 5. pressure difference
		_ofs << pf_rx_davp.get_plot_val(i) << WS;        // 6. average pressure difference
		_ofs << pf_rx_w.get_plot_val(i) << WS;           // 7. energy of relaxed state
		_ofs << pf_rx_wr.get_plot_val(i) << WS;          // 8. energy release
		_ofs << pf_rx_q.get_plot_val(i) << WS;           // 9. q value at relaxation radius
	
		_ofs << NL;
	}
}

/**
 * Stream the properties of the loop to a file stream.
 *
 * This function is a friend of the CLoop class.
 *
 * @param _ofs the output file stream.
 * @param _lp the loop object.
 * @return reference to output file.
 */
ofstream& operator<<(ofstream& _ofs, const CSmoothLoop& _lp) {
	_lp.out_props(_ofs);
	return _ofs;
}
