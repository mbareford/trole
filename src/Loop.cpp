/**
 * @file Loop.cpp
 *
 * This file defines the methods of the CLoop class.
 */

#include "StdInclude.h"

#include "Exception.h"
#include "GSL.h"
#include "NetBthCalculator.h"
#include "Profile.h"
#include "RelaxedLoop.h"
#include "RelaxationCalculator.h"
#include "Loop.h"


const double PI = CGSL::PI;

/// the number of relaxation regimes
const unsigned short RX_SCENARIO_CNT = 5;


/**
 * All attributes are zeroed.
 */
CLoop::CLoop(void): dimensions_initialised(false), coeffs_initialised(false), coeffs_normalised(false), field_initialised(false),
                    profiled(false), profiled_rx_scenarios(false), profiles_integrated(false), relaxed(false),
                    netcurrent(false), relaxable(false) {
	reset();
}

/**
 * Set attributes to zero.
 */
void CLoop::reset(void) {
	if (!dimensions_initialised) {
		r1 = r2 = r3 = r4 = 0.0;
		l = 0.0;
	}

	if (!field_initialised) {
		a1 = a2 = a3 = 0.0;
		b1 = b2 = b3 = c2 = c3 = 0.0;
	}

	if (!profiled) {
		pf_bz.reset();
		pf_bth.reset();
		pf_twist.reset();
		pf_shear.reset();
        pf_rsntk.reset();

		vpf_rx_a.clear();
		vpf_rx_b1.clear();
		vpf_rx_dp.clear();
		vpf_rx_davp.clear();
		vpf_rx_wr.clear();
		vpf_rx_q.clear();
	}

	if (!profiles_integrated) {
		lbth_0r1 = lbth_r1r2 = lbth_r2r3 = lbth_r3r4 = 0.0;
		rbz_0r1 = rbz_r1r2 = rbz_r2r3 = rbz_r3r4 = 0.0;
		baty_twist_0r1 = baty_twist_r1r2 = baty_twist_r2r3 = baty_twist_r3r4 = 0.0;
		twist_0r1 = twist_r1r2 = twist_r2r3 = twist_r3r4 = 0.0;
		shear_sq_0r1 = shear_sq_r1r2 = shear_sq_r2r3 = shear_sq_r3r4 = 0.0;
		shear_abs_0r1 = shear_abs_r1r2 = shear_abs_r2r3 = shear_abs_r3r4 = 0.0;
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
 * @param _r1 the core radius.
 * @param _r2 the outer layer radius.
 * @param _r3 the neutralisation layer radius.
 * @param _r4 the potential envelope radius.
 * @param _l the length of the loop in units of _r2.
 * @throw CException should any of the following occur:
 *        1. any of the radial boundaries has a zero value;
 *        2. radial boundaries in incorrect numerical order;
 *        3. length is not positive.
 */
void CLoop::init_dimensions(double _r1, double _r2, double _r3, double _r4, double _l) {
	if (0.0 == _r1 || 0.0 == _r2 || 0.0 == _r3 || 0.0 == _r4) {
		throw CException("Error! Zero radial boundary.");
	}

	if (_r1 > _r2 || _r1 > _r3 || _r1 > _r4 || _r2 > _r3 || _r2 > _r4 || _r3 > _r4) {
		throw CException("Error! Radial boundary out of order (must have 0 < r1 <= r2 <= r3 <= r4).");
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

	netcurrent = _r3 == _r2; // is true if there is NO neutralisation layer
	relaxable = false;

	reset();

	r1 = _r1;
	r2 = _r2;
	r3 = _r3;
	r4 = _r4;

	l = _l;

	dimensions_initialised = true;
}

/**
 * This function calculates the azimuthal field at r3 (in units of b1).
 *
 * This function is used in conjunction with the GSL root solver gsl_root_fsolver_brent.
 *
 * @param _a3 the alpha value for the neutralisation layer.
 * @param _param points to a CNetBthCalculator object
 * @return zero if _a3 really does neutralise the azimuthal field, otherwise
 * the azimuthal field at r3  (in units of b1).
 * @throw CException should any of the following occur:
 *        1. null _param;
 *        2. null CNetBthCalculator object.
 */
double CLoop::calc_r3bth(double _a3, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CNetBthCalculator* calc = (CNetBthCalculator*) (_param);
	if (NULL == calc) {
		throw CException("Error! NULL CNetBthCalculator object.");
	}

	return calc->netbth(_a3);
}

/**
 * This method calculates the alpha value (a3) for a loop's current neutralisation layer (r2-r3).
 * The neutralisation layer current must be such that it causes the azimuthal field to vanish at r3.
 * The value of a3 with the smallest absolute value is returned.
 *
 * @return the neutralisation layer alpha.
 * @throw CLoopException should any of the following occur:
 *        1. error finding neutralisation alpha.
 */
double CLoop::calc_nlalpha(void) {
	CNetBthCalculator calc(r1, r2, r3, a1, a2, b2, c2);

	void* param;
	param = &calc;

	double a3 = 0.0;

	try {
		a3 = CGSL::find_root(calc_r3bth, 50.0, 0.01, (a2 <= 0 ? 1 : -1), param);
	}
	catch (CGSLException& x) {
		throw CLoopException(a1, a2, "error finding neutralisation alpha", x);
	}

    return a3;
}

/**
 * This method approximates the alpha value (a3) for a loop's current neutralisation layer (r2-r3).
 * The tan function can be used to calculate a3 assuming that r3-r2 << r2, i.e., the current neutralisation
 * layer is very thin.
 * The neutralisation layer current must be such that it causes the azimuthal field to vanish at r3.
 *
 * This function is private, it is called by init_field_coeffs.
 *
 * @return the alpha value for the neutralisation layer
 * @throw CException should any of the following occur:
 *        1. potential loop;
 *        2. zero r2bth.
 */
double CLoop::approx_nlalpha(void) {
	double r2bz = 0.0, r2bth = 0.0;

	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa2 = (a2 < 0 ? -1.0 : 1.0);
	double aa1 = fabs(a1);
	double aa2 = fabs(a2);

	if (0.0 == a1 && 0.0 == a2) {
		throw CException("Error! Potential loop.");
	}

	if (0.0 == a2) {
		// potential outer layer
		r2bz = b2;
		r2bth = sa1*b1*J1(aa1*r1)*r1/r2;
	}
	else {
		// potential core or non-potential loop
		r2bz = b2*F0(aa2*r2);
		r2bth = sa2*b2*F1(aa2*r2);
	}

	if (0.0 == r2bth) {
		throw CException("Error! Zero r2bth.");
	}

	double a3 = atan2(r2bth, r2bz)/(r3-r2);
	double sa3 = -1.0*(r2bz/r2bth)*tan(a3*(r3-r2));

	if (sa3 < 0.0) {
	   a3 = -1.0*a3;
	}

	return a3;
}

/**
 * Calculates the magnetic field coefficients and the neutralisation alpha, a3
 * (unless netcurrent is true).
 *
 * This function is private, it is called by init_field.
 *
 * @param _b1norm if true b1 is normalised to 1.0, otherwise psi_0r4 is set to 1.0.
 * @param _approx_nlalpha is true if a3 is to be approximated rather than calculated
 */
void CLoop::init_field_coeffs(bool _b1norm, bool _approx_nlalpha) {
	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa2 = (a2 < 0 ? -1.0 : 1.0);
	double sa1a2 = sa1*sa2;

	coeffs_initialised = false;
	coeffs_normalised = false;

	// assume b1 is normalised
	b1 = 1.0;

	// calculate the coefficients for the outer layer
	if (0.0 == a1 && 0.0 == a2) {
		// potential loop
		b2 = 1.0;
		c2 = 0.0;
	}
	else if (0.0 == a1) {
		// potential core
		double delta = 2.0/(PI*aa2*r1);
		b2 = -1.0*Y1(aa2*r1)/delta;
		c2 = J1(aa2*r1)/delta;
	}
	else if (0.0 == a2) {
		// potential outer layer
		b2 = J0(aa1*r1);
		c2 = 0.0;
	}
	else {
		// entirely non-potential loop
		double delta = 2.0/(PI*aa2*r1);
		b2 = (sa1a2*J1(aa1*r1)*Y0(aa2*r1) - J0(aa1*r1)*Y1(aa2*r1))/delta;
		c2 = (J0(aa1*r1)*J1(aa2*r1) - sa1a2*J1(aa1*r1)*J0(aa2*r1))/delta;
	}


	if (netcurrent) {
		// ignore neutralisation; the loop carries net current
		// skip the current neutralisation layer
		a3 = b3 = c3 = 0.0;

		if (0.0 == a1 && 0.0 == a2) {
			// potential loop
			b4 = b1;
			c4 = 0.0;
		}
		else if (0.0 == a1) {
			// potential core
			b4 = b2*F0(aa2*r2);
			c4 = b2*F1(aa2*r2);
		}
		else if (0.0 == a2) {
			// potential outer layer
			b4 = b2;
			c4 = c2;
		}
		else {
			// entirely non-potential loop
			b4 = b2*F0(aa2*r2);
			c4 = b2*F1(aa2*r2);
		}
	}
	else {
		// the loop carries zero net current
		// calculate what a3 needs to be in order for the azimuthal field to vanish at r3
		if (0.0 == a1 && 0.0 == a2) {
			// potential loop
			a3 = 0.0;
		}
		else if (_approx_nlalpha) {
			a3 = approx_nlalpha();
		}
		else {
			a3 = calc_nlalpha();
		}

		double aa3 = fabs(a3);
		double sa3 = (a3 < 0 ? -1.0 : 1.0);
		double sa1a3 = sa1*sa3;
		double sa2a3 = sa2*sa3;

		// calculate the coefficients for the neutralisation layer and potential envelope
		if (0.0 == a1 && 0.0 == a2) {
			// potential loop
			b4 = b3 = b1;
			c4 = c3 = 0.0;
		}
		else if (0.0 == a1) {
			// potential core (a3 must be non-zero)
			double delta = 2.0/(PI*aa3*b2*r2);
			b3 = (sa2a3*F1(aa2*r2)*Y0(aa3*r2) - F0(aa2*r2)*Y1(aa3*r2))/delta;
			c3 = (F0(aa2*r2)*J1(aa3*r2) - sa2a3*F1(aa2*r2)*J0(aa3*r2))/delta;

			b4 = b3*G0(aa3*r3);
			c4 = 0.0;
		}
		else if (0.0 == a2) {
			// potential outer layer (a3 must be non-zero)
			double delta = 2.0/(PI*aa3*r2);
			b3 = (sa1a3*b1*J1(aa1*r1)*(r1/r2)*Y0(aa3*r2) - b2*Y1(aa3*r2))/delta;
			c3 = (b2*J1(aa3*r2) - sa1a3*b1*J1(aa1*r1)*(r1/r2)*J0(aa3*r2))/delta;

			b4 = b3*G0(aa3*r3);
			c4 = 0.0;
		}
		else {
			// core and outer layer are non-potential
			if (0.0 == a3) {
				// potential neutralisation layer (outer layer neutralises core)
				b3 = b2*F0(aa2*r2);
				c3 = 0.0;

				b4 = b3;
				c4 = 0.0;
			}
			else {
				// entirely non-potential loop
				double delta = 2.0/(PI*aa3*b2*r2);
				b3 = (sa2a3*F1(aa2*r2)*Y0(aa3*r2) - F0(aa2*r2)*Y1(aa3*r2))/delta;
				c3 = (F0(aa2*r2)*J1(aa3*r2) - sa2a3*F1(aa2*r2)*J0(aa3*r2))/delta;

				b4 = b3*G0(aa3*r3);
				c4 = 0.0;
			}
		}
	} // end of netcurrent else clause


	// so far all field coefficients are in units of b1
	coeffs_initialised = true;

	// if the total axial flux is normalised rescale all magnetic field coefficients
	if (!_b1norm) {
		b1 = 1.0/calc_psi(r4);
	}

	// set remaining coefficients so that they are no longer in units of b1
	b2 = b2*b1;
	c2 = c2*b1;
	b3 = b3*b1;
	c3 = c3*b1;
	b4 = b4*b1;
	c4 = c4*b1;

	coeffs_normalised = true;
}

/**
 * Initialises the loop object's alpha values.
 *
 * If this function completes successfully the field initialisation flag is
 * set to true.
 *
 * @param _a1 the alpha value within the core.
 * @param _a2 the alpha value within the outer layer.
 * @param _b1norm if true b1 is normalised to 1.0, otherwise psi_0r3 is set to 1.0.
 * @param _approx_nlalpha is true if a3 is to be approximated rather than calculated.
 * @param _relaxable is true if it is appropriate for this loop to be relaxed.
 * @throw CException should any of the following occur:
 *        1. loop dimensions uninitialised.
 */
void CLoop::init_field(double _a1, double _a2, bool _b1norm, bool _approx_nlalpha, bool _relaxable) {
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

	a1 = _a1;
	a2 = _a2;

	init_field_coeffs(_b1norm, _approx_nlalpha);

	field_initialised = true;
}

/**
 * This function permits the field to be reinitialised. The normalisation
 * is kept the same.
 *
 * @param _a1 the alpha value within the core.
 * @param _a2 the alpha value within the outer layer.
 * @param _approx_nlalpha is true if a3 is to be approximated rather than calculated.
 * @param _relaxable is true if it is appropriate for this loop to be relaxed.
 * @throw CException should any of the following occur:
 *        1. loop dimensions uninitialised.
 */
void CLoop::init_field(double _a1, double _a2, bool _approx_nlalpha, bool _relaxable) {
	if (!field_initialised) {
		throw CException("Error! Loop field uninitialised.");
	}

	init_field(_a1, _a2, (1.0 == b1), _approx_nlalpha, _relaxable);
}


/**
 * Calculate axial component of the magnetic field for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @return the axial component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 */
double CLoop::calc_bz(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r4) {
		throw CException("Error! Radial value out of bounds.");
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double aa3 = fabs(a3);
	double bz = 0.0;

	if (0.0 == r || (0.0 == a1 && 0.0 == a2) || (r <= r1 && 0.0 == a1)) {
		// on axis or potential loop or (within core and potential core)
		bz = b1;
	}
	else if (r <= r1) {
		// within core
		bz = b1*J0(aa1*r);
	}
	else if (r <= r2) {
		// within outer layer
		if (0.0 == a2) {
			// potential outer layer
			bz = b2;
		}
		else {
			// potential core or non-potential loop
			bz = b2*J0(aa2*r) + c2*Y0(aa2*r);
		}
	}
	else if (r <= r3) {
		// within neutralisation layer
		if (0.0 == a3) {
			// outer layer neutralises core
			bz = b3;
		}
		else {
			bz = b3*J0(aa3*r) + c3*Y0(aa3*r);
		}
	}
	else {
		// within potential envelope
		bz = b4;
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
double CLoop::calc_bz(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CLoop* loop = (CLoop*) (_param);
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
double CLoop::calc_bth(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r4) {
		throw CException("Error! Radial value out of bounds.");
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double aa3 = fabs(a3);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa2 = (a2 < 0 ? -1.0 : 1.0);
	double sa3 = (a3 < 0 ? -1.0 : 1.0);
	double bth = 0.0;

	if (0.0 == r || (0.0 == a1 && 0.0 == a2) || (r <= r1 && 0.0 == a1)) {
		// on axis or potential loop or (within core and potential core)
		bth = 0.0;
	}
	else if (r <= r1) {
		// within core
		bth = sa1*b1*J1(aa1*r);
	}
	else if (r <= r2) {
		// within outer layer
		if (0.0 == a2) {
			// potential outer layer
			bth = sa1*b1*J1(aa1*r1)*r1/r;
		}
		else {
			// potential core or non-potential loop
			bth = sa2*(b2*J1(aa2*r) + c2*Y1(aa2*r));
		}
	}
	else if (r <= r3) {
		// within neutralisation layer
		if (0.0 == a3) {
			// outer layer neutralises core
			// a1 and a2 are both non-zero
			bth = 0.0;
			// if in fact a3 did not have a neutralising function
			// bth would be sa2*b2*F1(aa2*r2)*r2/r;
		}
		else {
			bth = sa3*(b3*J1(aa3*r) + c3*Y1(aa3*r));
		}
	}
	else {
		// within potential envelope
		if (netcurrent) {
			if (0.0 == a2) {
				// potential outer layer
				bth = sa1*b1*J1(aa1*r)*r1/r;
			}
			else {
				// potential core or entirely non-potential loop
				bth = sa2*b2*F1(aa2*r2)*r2/r;
			}
		}
		else {
			bth = 0.0;
			// if in fact a3 did not have a neutralising function
			// bth would be sa3*b3*G1(aa3*r3)*r3/r;
		}
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
double CLoop::calc_bth(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CLoop* loop = (CLoop*) (_param);
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
double CLoop::calc_p(double _r) const {
	double bz = calc_bz(_r);
	double bth = calc_bth(_r);

	return (bz*bz + bth*bth)/2.0;
}

/**
 * This is the function used by to calculate the twist at a given radius.
 * It is also used by the GSL routine gsl_integration_qags to calculate the twist integrals.
 *
 * @param _r the radial coordinate
 * @return the twist at _r
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 * @throw CSingularityException should any of the following occur:
 *        1. zero bz.
 */
double CLoop::calc_twist(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r4) {
		throw CException("Error! CLoop::calc_twist(), radial value out of bounds.");
	}

	double tw = 0.0;

	if (0.0 == r) {
		double aa1 = fabs(a1);
		double sa1 = (a1 < 0 ? -1.0 : 1.0);

		tw = sa1*l*aa1/2.0;
	}
	else {
		double bz = calc_bz(r);
		double bth = calc_bth(r);

		if (0.0 == bz) {
			throw CSingularityException(r, "Error! CLoop::calc_twist(), zero bz.");
		}
		else {
			tw = (l*bth) / (r*bz);
		}
	}

	return tw;
}

/**
 * This is the function used by to calculate the shear at a given radius.
 * It is also used by the GSL routine gsl_integration_qags to calculate the shear integrals.
 *
 * @param _r the radial coordinate
 * @return the shear at _r
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 * @throw CSingularityException should any of the following occur:
 *        1. zero bz.
 */
double CLoop::calc_shear(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r4) {
		throw CException("Error! CLoop::calc_shear(), radial value out of bounds.");
	}

	double sh = 0.0;
	double bz = calc_bz(r);
	double bth = calc_bth(r);

	if (0.0 == bz) {
		throw CSingularityException(r, "Error! CLoop::calc_shear(), zero bz.");
	}
	else {
		sh = bth/bz;
	}

	return sh;
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
double CLoop::calc_rsntk(double _r) const {
	double r = _r;
	if (r < 0.0 || r > r4) {
	    throw CException("Error! CLoop::calc_rsntk(), radial value out of bounds.");
	}

	double rsntk = 0.0;

    if (0.0 == r) {
		double aa1 = fabs(a1);
		double sa1 = (a1 < 0 ? -1.0 : 1.0);

		rsntk = -sa1*(aa1/2.0);
	}
	else {	  
	    double bz = calc_bz(r);
	    double bth = calc_bth(r);

		if (0.0 == bz) {
	        throw CSingularityException(r, "Error! CLoop::calc_rsntk(), zero bz.");
		}
		else {
	        rsntk = -bth/(r*bz);
		}
	}

	return rsntk;
}

/**
 * This is the function used by the GSL routine gsl_integration_qags to calculate the integral of the
 * the total twist weighted by area.
 *
 * @param _r the radial coordinate
 * @param _param loop object pointer
 * @return the twist at _r
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null loop object.
 */
double CLoop::calc_twist(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CLoop* loop = (CLoop*) (_param);
	if (NULL == loop) {
		throw CException("Error! NULL loop object.");
	}

	return 2*PI*_r*loop->calc_twist(_r);
}

/**
 * This is the function used by the GSL routine gsl_integration_qags to calculate the integral of the
 * the total twist according to Baty 2001.
 *
 * @param _r the radial coordinate
 * @param _param points to a loop object
 * @return the twist at _r
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null loop object.
 */
double CLoop::calc_baty_twist(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CLoop* loop = (CLoop*) (_param);
	if (NULL == loop) {
		throw CException("Error! NULL loop object.");
	}

	return loop->calc_twist(_r);
}

/**
 * This is the function used by the GSL routine gsl_integration_qags to calculate the integral of the
 * the total shear squared weighted by area.
 *
 * @param _r the radial coordinate
 * @param _param points to a loop object
 * @return the shear squared at _r
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null loop object.
 */
double CLoop::calc_shear_sq(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CLoop* loop = (CLoop*) (_param);
	if (NULL == loop) {
		throw CException("Error! NULL loop object.");
	}

	double sh = loop->calc_shear(_r);

	return 2*PI*_r*sh*sh;
}

/**
 * This is the function used by the GSL routine gsl_integration_qags to calculate the integral of the
 * the total absolute shear weighted by area.
 *
 * @param _r the radial coordinate
 * @param _param points to a loop object
 * @return the absolute shear at _r
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null loop object.
 */
double CLoop::calc_shear_abs(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CLoop* loop = (CLoop*) (_param);
	if (NULL == loop) {
		throw CException("Error! NULL loop object.");
	}

	double sh = loop->calc_shear(_r);

	return 2*PI*_r*fabs(sh);
}

/**
 * Return the axial flux between 0 and _r.
 *
 * @param _r the radial coordinate that is the upper bound for the axial flux calculation
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field coefficients uninitialised;
 *        2. b1 is uninitialised;
 *        3. radial value out of bounds.
 */
double CLoop::calc_psi(double _r) const {
	if (!dimensions_initialised || !coeffs_initialised) {
		throw CException("Error! dimensions and/or field coefficients uninitialised.");
	}

	if (0.0 == b1) {
		throw CException("Error! b1 is uninitialised.");
	}

	if (_r < 0.0 || _r > r4) {
		throw CException("Error! Radial value out of bounds.");
	}

	if (0.0 == _r) {
		return 0.0;
	}

	if (0.0 == a1 && 0.0 == a2) {
		return PI*b1*_r*_r;
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double aa3 = fabs(a3);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa2 = (a2 < 0 ? -1.0 : 1.0);
	double sa3 = (a3 < 0 ? -1.0 : 1.0);
	double sa1a2 = sa1*sa2;
	double sa1a3 = sa1*sa3;
	double sa2a3 = sa2*sa3;

	double psi = 0.0;

	do {
		if (0.0 == a1 && 0.0 == a2) {
			// potential loop
			psi = PI*b1*_r*_r;
			break;
		}

		double r = _r < r1 ? _r : r1;

		// core
		if (0.0 == a1) {
			// potential core
			psi = PI*b1*r*r;
		}
		else {
			// non-potential core
			psi = ((2.0*PI*b1*r)/aa1)*J1(aa1*r);
		}

		if (_r <= r1) {
			break;
		}

		r = _r < r2 ? _r : r2;

		// outer layer
		if (0.0 == a2) {
			// potential outer layer
			psi += (PI*b2*(r*r - r1*r1));
		}
		else {
			// non-potential outer layer
			if (0.0 == a1) {
				// potential core
				psi += ((2.0*PI*b2*r)/aa2)*F1(aa2*r);
			}
			else {
				// non-potential core
				psi += ((2.0*PI*b2*r)/aa2)*F1(aa2*r) - 2.0*PI*b1*r1*J1(aa1*r1)*(sa1a2/aa2);
			}
		}

		if (_r <= r2) {
			break;
		}


		if (!netcurrent) {
			r = _r < r3 ? _r : r3;

			// current neutralisation layer
			if (0.0 == a2) {
				// potential outer layer, a3 must be non-zero
				psi += ((2.0*PI*b3*r)/aa3)*G1(aa3*r) - 2.0*PI*b1*J1(aa1*r1)*r1*(sa1a3/aa3);
			}
			else {
				if (0.0 == a1) {
					// potential core, a3 must be non-zero
					psi += ((2.0*PI*b3*r)/aa3)*G1(aa3*r) - 2.0*PI*b2*r2*F1(aa2*r2)*(sa2a3/aa3);
				}
				else {
					// core and outer layer are non-potential
					if (0.0 == a3) {
						// outer layer neutralises core
						psi += (PI*b3*(r*r - r2*r2));
					}
					else {
						// entirely non-potential loop
						psi += (2.0*PI*b3/aa3)*(r*G1(aa3*r) - r2*G1(aa3*r2));
					}
				}
			}

			if (_r <= r3) {
				break;
			}

		} // end of <if (!netcurrent)> clause


		r = _r < r4 ? _r : r4;

		// potential envelope
		psi += netcurrent ? (PI*b4*(r*r - r2*r2)) : (PI*b4*(r*r - r3*r3));

	} while (0);

	return psi;
}

/**
 * Return the helicity between 0 and _r.
 *
 * @param _r the radial coordinate that is the upper bound for the helicity calculation
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. radial value out of bounds.
 */
double CLoop::calc_k(double _r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (_r < 0.0 || _r > r4) {
		throw CException("Error! Radial value out of bounds.");
	}

	if (0.0 == _r || (0.0 == a1 && 0.0 == a2)) {
		return 0.0;
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double aa3 = fabs(a3);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa2 = (a2 < 0 ? -1.0 : 1.0);
	double sa3 = (a3 < 0 ? -1.0 : 1.0);
	double sa1a2 = sa1*sa2;
	double sa2a3 = sa2*sa3;
	double sa1a3 = sa1*sa3;

	double uk = 0.0, k = 0.0;

	do {
		if (0.0 == a1 && 0.0 == a2) {
			// potential loop
			k = 0.0;
			break;
		}

		double r = _r < r1 ? _r : r1;

		// core
		if (0.0 == a1) {
			// potential core
			k = 0.0;
		}
		else {
			// non-potential core
			uk = (2.0*PI*l*b1*b1/aa1)*(r*r*(J0(aa1*r)*J0(aa1*r) + J1(aa1*r)*J1(aa1*r)) - ((2.0*r)/aa1)*J0(aa1*r)*J1(aa1*r));

			k = sa1*uk;
		}

		if (_r <= r1) {
			break;
		}

		r = _r < r2 ? _r : r2;

		// outer layer
		if (0.0 == a2) {
			// potential outer layer
			uk = PI*l*b2*b1*J1(aa1*r1)*r1*(r*r - r1*r1) + 2.0*PI*l*b1*J1(aa1*r1)*r1*r1*(2*(b1*J1(aa1*r1)/aa1) - b2*r1)*log(r/r1);

			k += sa1*uk;
		}
		else {
			// non-potential outer layer
			if (0.0 == a1) {
				// potential core
				uk = (2.0*PI*l*b2*b2/aa2)*(r*r*(F0(aa2*r)*F0(aa2*r) + F1(aa2*r)*F1(aa2*r)) - ((2.0*r)/aa2)*F0(aa2*r)*F1(aa2*r))
					 -
					 ((2.0*PI*l*b1*b2*r1*r1)/aa2)*F0(aa2*r);

				k += sa2*uk;
			}
			else {
				// non-potential core
				uk = (2.0*PI*l*b2*b2/aa2)*(r*r*(F0(aa2*r)*F0(aa2*r) + F1(aa2*r)*F1(aa2*r)) - ((2.0*r)/aa2)*F0(aa2*r)*F1(aa2*r))
					 -
					 (2.0*PI*l*b2*b2/aa2)*(r1*r1*(F0(aa2*r1)*F0(aa2*r1) + F1(aa2*r1)*F1(aa2*r1)) - ((2.0*r1)/aa2)*F0(aa2*r1)*F1(aa2*r1))
					 +
					 (4.0*PI*l*b2/aa2)*(F0(aa2*r1) - F0(aa2*r))*b1*r1*J1(aa1*r1)*((1.0/aa1) - (sa1a2/aa2));

				k += sa2*uk;
			}
		}


		if (_r <= r2) {
			break;
		}

		if (!netcurrent) {
			r = _r < r3 ? _r : r3;

			// current neutralisation layer
			if (0.0 == a2) {
				// potential outer layer, a3 must be non-zero
				uk = (2.0*PI*l*b3*b3/aa3)*(r*r*(G0(aa3*r)*G0(aa3*r) + G1(aa3*r)*G1(aa3*r)) - ((2.0*r)/aa3)*G0(aa3*r)*G1(aa3*r))
					 -
					 (2.0*PI*l/aa3)*(b2*b2*r2*r2 + b1*b1*J1(aa1*r1)*J1(aa1*r1)*r1*r1 - sa1a3*(2.0*b2*b1*J1(aa1*r1)*r1)/aa3)
					 +
					 (2.0*PI*l*b3/aa3)*(b2/b3 - G0(aa3*r))*((2.0*b1*J1(aa1*r1)*r1)/aa1 - (sa1a3*2.0*b1*J1(aa1*r1)*r1)/aa3 + b2*(r2*r2 - r1*r1));
			}
			else {
				// non-potential outer layer
				if (0.0 == a3) {
					// outer layer neutralises core
					// a1 must be non-zero
					uk = 0.0;
				}
				else {
					uk = (2.0*PI*l*b3*b3/aa3)*(r*r*(G0(aa3*r)*G0(aa3*r) + G1(aa3*r)*G1(aa3*r)) - ((2.0*r)/aa3)*G0(aa3*r)*G1(aa3*r))
						 -
						 (2.0*PI*l*b3*b3/aa3)*(r2*r2*(G0(aa3*r2)*G0(aa3*r2) + G1(aa3*r2)*G1(aa3*r2)) - ((2.0*r2)/aa3)*G0(aa3*r2)*G1(aa3*r2));

					if (0.0 == a1) {
						// potential core, a3 must be non-zero
						uk += (2.0*PI*l*b3/aa3)*(G0(aa3*r2) - G0(aa3*r))*(2.0*b2*r2*F1(aa2*r2)*((1.0/aa2) - (sa2a3/aa3)) + b1*r1*r1);
					}
					else {
						// non-potential loop
						uk += (4.0*PI*l*b3/aa3)*(G0(aa3*r2) - G0(aa3*r))*(b2*r2*F1(aa2*r2)*((1.0/aa2) - (sa2a3/aa3))
							    								        + b1*r1*J1(aa1*r1)*((1.0/aa1) - (sa1a2/aa2)));
					}
				}
			}
			k += sa3*uk;

			if (_r <= r3) {
				break;
			}

		} // end of <if (!netcurrent)> clause


		r = _r < r4 ? _r : r4;

		// potential envelope
		if (netcurrent) {
			if (0.0 == a2) {
				// potential outer layer
				uk = (2.0*PI*l*b1*b1)*((2.0*r1*r1*J1(aa1*r1)*J1(aa1*r1)*log(r/r2)/aa1)
									   -
									   (r1*r1*r1*J0(aa1*r1)*J1(aa1*r1)*log(r/r2))
									   +
									   (0.5*r1*J0(aa1*r1)*J1(aa1*r1)*(r*r - r2*r2)));
			}
			else {
				// potential core or entirely non-potential loop
				uk = 2.0*l*c4*r2*((calc_psi(r2) - PI*b4*r2*r2)*log(r/r2) + (PI*b4/2.0)*(r*r - r2*r2));
			}
			k += sa2*uk;
		}
		else {
			uk = 0.0;
			// if a3 did not in fact have a neutralising function
			// uk would be 2.0*l*b3*G1(aa3*r3)*r3*((calc_psi(r3) - PI*b4*r3*r3)*log(r/r3) + (PI*b4/2.0)*(r*r - r3*r3))
			k += sa3*uk;
		}

	} while (0);

	return k;
}

/**
 * Return the energy between 0 and _r.
 *
 * @param _r the radial coordinate that is the upper bound for the energy calculation
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. radial value out of bounds.
 */
double CLoop::calc_w(double _r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (_r < 0.0 || _r > r4) {
		throw CException("Error! Radial value out of bounds.");
	}

	if (0.0 == _r) {
		return 0.0;
	}

	if (0.0 == a1 && 0.0 == a2) {
		return (l*PI*b1*b1*_r*_r)/2.0;
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double aa3 = fabs(a3);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa3 = (a3 < 0 ? -1.0 : 1.0);
	double sa1a3 = sa1*sa3;

	double w = 0.0;

	do {
		if (0.0 == a1 && 0.0 == a2) {
			// potential loop
			w = (b1*b1*_r*_r)/2.0;
			break;
		}

		double r = _r < r1 ? _r : r1;

		// core
		if (0.0 == a1) {
			// potential core
			w = (b1*b1*r*r)/2.0;
		}
		else {
			// non-potential core
			w = b1*b1*(r*r*(J0(aa1*r)*J0(aa1*r) + J1(aa1*r)*J1(aa1*r)) - (r/aa1)*J0(aa1*r)*J1(aa1*r));
		}

		if (_r <= r1) {
			break;
		}

		r = _r < r2 ? _r : r2;

		// outer layer
		if (0.0 == a2) {
			// potential outer layer
			w += (b2*b2/2.0)*(r*r - r1*r1) + (r1*r1*b1*b1*J1(aa1*r1)*J1(aa1*r1)*log(r/r1));
		}
		else {
			// non-potential outer layer
			if (0.0 == a1) {
				// potential core
				w += b2*b2*(r*r*(F0(aa2*r)*F0(aa2*r) + F1(aa2*r)*F1(aa2*r)) - (r/aa2)*F0(aa2*r)*F1(aa2*r)) - b1*b1*r1*r1;
			}
			else {
				// non-potential core
				w += b2*b2*(r*r*(F0(aa2*r)*F0(aa2*r) + F1(aa2*r)*F1(aa2*r)) - (r/aa2)*F0(aa2*r)*F1(aa2*r)
							-
							r1*r1*(F0(aa2*r1)*F0(aa2*r1) + F1(aa2*r1)*F1(aa2*r1)) + (r1/aa2)*F0(aa2*r1)*F1(aa2*r1));
			}
		}

		if (_r <= r2) {
			break;
		}

		if (!netcurrent) {
			r = _r < r3 ? _r : r3;

			// current neutralisation layer
			if (0.0 == a2) {
				// potential outer layer, a3 must be non-zero
				w += b3*b3*(r*r*(G0(aa3*r)*G0(aa3*r) + G1(aa3*r)*G1(aa3*r)) - (r/aa3)*G0(aa3*r)*G1(aa3*r))
					 -
					 (b2*b2*r2*r2 + b1*b1*J1(aa1*r1)*J1(aa1*r1)*r1*r1 - (sa1a3*b2*b1*J1(aa1*r1)*r1)/aa3);
			}
			else {
				// non-potential outer layer
				if (0.0 == a3) {
					// outer layer neutralises core
					w += ((b3*b3/2.0)*(r*r - r2*r2));
				}
				else {
					w += b3*b3*(r*r*(G0(aa3*r)*G0(aa3*r) + G1(aa3*r)*G1(aa3*r)) - (r/aa3)*G0(aa3*r)*G1(aa3*r)
								-
								r2*r2*(G0(aa3*r2)*G0(aa3*r2) + G1(aa3*r2)*G1(aa3*r2)) + (r2/aa3)*G0(aa3*r2)*G1(aa3*r2));
				}
			}

			if (_r <= r3) {
				break;
			}

		} // end of <if (!netcurrent)> clause


		r = _r < r4 ? _r : r4;

		// potential envelope

		if (netcurrent) {
			if (0.0 == a2) {
				// potential outer layer
				w += ((b1*b1*J0(aa1*r1)*J0(aa1*r1/2.0)*(r*r - r2*r2)) + (b1*b1*J1(aa1*r1)*J1(aa1*r1)*r1*r1*log(r/r2)));
			}
			else {
				// potential core or entirely non-potential loop
				w += ((b4*b4/2.0)*(r*r - r2*r2)) + (c4*c4*r2*r2*log(r/r2));
			}
		}
		else {
			w += ((b4*b4/2.0)*(r*r - r3*r3));
			// if in fact a3 did not have a neutralising function
			// w would be ((b4*b4/2.0)*(r*r - r3*r3)) + (b3*b3*G1(aa3*r3)*G1(aa3*r3)*r3*r3*log(r/r3));
		}

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
 * This function is used in conjunction with GSL root solver, gsl_root_fsolver_brent.
 *
 * @param _a1 the prospective relaxed alpha value.
 * @param _param points to a CRelaxationCalculator object.
 * @return zero if _rxa really does represent the relaxed state.
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null CRelaxationCalculator object.
 */
double CLoop::calc_consrv_err(double _a1, void* _param) {
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
 * @param _rnl current neutralisation layer radius; if _rnl == _r the relaxed state has a current neutralisation surface,
 * if _rnl == 0 the relaxed state has net current
 * @param _psi the axial flux that is conserved by the relaxed loop
 * @param _psi_r the radial range of the axial flux that is conserved
 * @param _consrv the value k/(psi*psi) that is conserved by the relaxed loop
 * @param _consrv_r the radial range of the threshold state over which _consrv is calculated
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. invalid relaxation radius;
 *        3. invalid relaxation neutralisation layer radius.
 * @throw CLoopException should any of the following occur:
 *        1. error finding a1 (relaxation).
 */
CRelaxedLoop CLoop::calc_rx(double _r, double _rnl, double _psi, double _psi_r, double _consrv, double _consrv_r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (_r <= 0.0 || _r > r4) {
		throw CException("Error! Invalid relaxation radius.");
	}

	if (0.0 != _rnl && (_rnl < _r || _rnl > r4)) {
		throw CException("Error! Invalid relaxation neutralisation layer radius.");
	}

	CRelaxedLoop rx;

	// initialise the dimensions of the relaxed loop
	if (_rnl > _r) {
		// the relaxed state has a neutralisation layer
		rx.init_dimensions(_r, _rnl, r4, l);
	}
	else {
		rx.init_dimensions(_r, r4, l);
	}

	if (0.0 != _rnl) {
		// the relaxed state has zero net current, which is achieved via
		// neutralisation layer (_rnl > _r) or a neutralisation surface (_rnl == _r)

		// initialise the relaxed state's potential envelope so that it has the
		// same axial field as that of the envelope for the threshold state
		rx.fix_envelope(b4);
	}


	// calculate relaxed state according to conservation parameters
	double rxa = 0.0;

	CRelaxationCalculator calc(_psi, _psi_r, _consrv, _consrv_r, rx);
	void* param = &calc;

	try {
		rxa = CGSL::find_root(CLoop::calc_consrv_err, 50.0, 0.01, 0, param);
	}
	catch (CGSLException& x) {
		throw CLoopException(a1, a2, "error finding a1 (relaxation)", x);
	}

	rx.init_field(_psi, _psi_r, rxa);

	return rx;
}


/**
 * Integrate area-weighted twist numerically.
 *
 * @param _param points to a loop object.
 */
void CLoop::int_twist(void* _param) {
	twist_0r1 = CGSL::integrate(calc_twist, 0.0, r1, _param);
	twist_r1r2 = CGSL::integrate(calc_twist, r1, r2, _param);
	twist_r2r3 = netcurrent ? 0.0 : CGSL::integrate(calc_twist, r2, r3, _param);
	twist_r3r4 = CGSL::integrate(calc_twist, r3, r4, _param);
}

/**
 * Integrate twist numerically according to Baty 2001.
 *
 * @param _param points to a loop object.
 */
void CLoop::int_baty_twist(void* _param) {
	baty_twist_0r1 = CGSL::integrate(calc_baty_twist, 0.0, r1, _param);
	baty_twist_r1r2 = CGSL::integrate(calc_baty_twist, r1, r2, _param);
	baty_twist_r2r3 = netcurrent ? 0.0 : CGSL::integrate(calc_baty_twist, r2, r3, _param);
	baty_twist_r3r4 = CGSL::integrate(calc_baty_twist, r3, r4, _param);
}

/**
 * Integrate twist analytically according to Velli et al. 1990.
 */
void CLoop::int_velli_twist(void) {
	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double aa3 = fabs(a3);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa2 = (a2 < 0 ? -1.0 : 1.0);
	double sa3 = (a3 < 0 ? -1.0 : 1.0);
	double sa1a3 = sa1*sa3;

	if (0.0 == a1 && 0.0 == a2) {
		// potential loop
		lbth_0r1 = 0.0;
		rbz_0r1 = b1;

		lbth_r1r2 = 0.0;
		rbz_r1r2 = b1;

		lbth_r2r3 = 0.0;
		rbz_r2r3 = b1;

		lbth_r3r4 = 0.0;
		rbz_r3r4 = b1;
	}
	else if (0.0 == a1) {
		// potential core
		lbth_0r1 = 0.0;
		rbz_0r1 = (b1*r1*r1/2.0);

		lbth_r1r2 = (sa2*l/aa2)*(b1 - b2*F0(aa2*r2));
		rbz_r1r2 = (b2*r2/aa2)*F1(aa2*r2);

		if (netcurrent) {
			lbth_r2r3 = 0.0;
			rbz_r2r3 = 0.0;

			lbth_r3r4 = sa2*l*b2*F1(aa2*r2)*r2*log(r4/r2);
			rbz_r3r4 = 0.5*b2*F0(aa2*r2)*(r4*r4 - r2*r2);
		}
		else {
			// the loop carries zero net current
			lbth_r2r3 = (l*sa3*b3/aa3)*(G0(aa3*r2) - G0(aa3*r3));
			rbz_r2r3 = (b3/aa3)*(r3*G1(aa3*r3) - r2*G1(aa3*r2));

			lbth_r3r4 = 0.0;
			// if in fact a3 did not have a neutralising function
			// lbth_r3r4 would be sa3*l*b3*G1(aa3*r3)*r3*log(r4/r3);
			rbz_r3r4 = 0.5*b4*(r4*r4 - r3*r3);
		}
	}
	else if (0.0 == a2) {
		// potential outer layer
		lbth_0r1 = (l*sa1*b1/aa1)*(1-J0(aa1*r1));
		rbz_0r1 = (b1*r1/aa1)*(J1(aa1*r1));

		lbth_r1r2 = sa1*l*b1*J1(aa1*r1)*r1*log(r2/r1);
		rbz_r1r2 = 0.5*b2*(r2*r2 - r1*r1);

		if (netcurrent) {
			lbth_r2r3 = 0.0;
			rbz_r2r3 = 0.0;

			lbth_r3r4 = sa1*l*b1*J1(aa1*r1)*r1*log(r4/r2);
			rbz_r3r4 = 0.5*b2*(r4*r4 - r2*r2);
		}
		else {
			// the loop carries zero net current
			lbth_r2r3 = (sa3*l/aa3)*(b2 - b3*G0(aa3*r3));
			rbz_r2r3 = (1.0/aa3)*(b3*r3*G1(aa3*r3) - sa1a3*b1*J1(aa1*r1)*r1);

			lbth_r3r4 = 0.0;
			// if in fact a3 did not have a neutralising function
			// lbth_r3r4 would be sa3*l*b3*G1(aa3*r3)*r3*log(r4/r3);
			rbz_r3r4 = 0.5*b4*(r4*r4 - r3*r3);
		}
	}
	else {
		// non-potential loop
		lbth_0r1 = (l*sa1*b1/aa1)*(1-J0(aa1*r1));
		rbz_0r1 = (b1*r1/aa1)*(J1(aa1*r1));

		lbth_r1r2 = (l*sa2*b2/aa2)*(F0(aa2*r1) - F0(aa2*r2));
		rbz_r1r2 = (b2/aa2)*(r2*F1(aa2*r2) - r1*F1(aa2*r1));

		if (netcurrent) {
			lbth_r2r3 = 0.0;
			rbz_r2r3 = 0.0;

			lbth_r3r4 = sa2*l*r2*(b2*J1(aa2*r2) + c2*Y1(aa2*r2))*log(r4/r2);
			rbz_r3r4 = 0.5*b2*F0(aa2*r2)*(r4*r4 - r2*r2);
		}
		else {
			// the loop carries zero net current
			if (0.0 == a3) {
				// outer layer neutralises core
				lbth_r2r3 = 0.0;
				rbz_r2r3 = 0.5*b3*(r3*r3 - r2*r2);
			}
			else {
				// current neutralising layer with non-zero a3
				lbth_r2r3 = (l*sa3*b3/aa3)*(G0(aa3*r2) - G0(aa3*r3));
				rbz_r2r3 = (b3/aa3)*(r3*G1(aa3*r3) - r2*G1(aa3*r2));
			}

			lbth_r3r4 = 0.0;
			// if in fact a3 did not have a neutralising function
			// lbth_r3r4 would be sa3*l*b3*G1(aa3*r3)*r3*log(r4/r3)
			rbz_r3r4 = 0.5*b4*(r4*r4 - r3*r3);
		}
	} // end of non-potential loop else clause
}

/**
 * Calculates numerically the shear squared integral over the loop's radial length.
 *
 * @param _param points to a loop object.
 */
void CLoop::int_shear_sq(void* _param) {
	shear_sq_0r1 = CGSL::integrate(calc_shear_sq, 0.0, r1, _param);
	shear_sq_r1r2 = CGSL::integrate(calc_shear_sq, r1, r2, _param);
	shear_sq_r2r3 = netcurrent ? 0.0 : CGSL::integrate(calc_shear_sq, r2, r3, _param);
	shear_sq_r3r4 = CGSL::integrate(calc_shear_sq, r3, r4, _param);
}

/**
 * Calculates numerically the absolute shear integral over the loop's radial length.
 *
 * @param _param points to a loop object.
 */
void CLoop::int_shear_abs(void* _param) {
	shear_abs_0r1 = CGSL::integrate(calc_shear_abs, 0.0, r1, _param);
	shear_abs_r1r2 = CGSL::integrate(calc_shear_abs, r1, r2, _param);
	shear_abs_r2r3 = netcurrent ? 0.0 : CGSL::integrate(calc_shear_abs, r2, r3, _param);
	shear_abs_r3r4 = CGSL::integrate(calc_shear_abs, r3, r4, _param);
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
void CLoop::profile(unsigned long _rs, bool pf_rx_scenarios/*=false*/) {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	profiled = false;
	profiled_rx_scenarios = false;
	profiles_integrated = false;

	reset();

	// construct field profiles
	//////////////////////////////////////////////////////////////////////////
	double r, bz, bth, tw, sh, rsntk;

	void* param = this;

	pf_bz.init(0.0, r4, _rs);
	pf_bth.init(0.0, r4, _rs);
	pf_twist.init(0.0, r4, _rs);
	pf_shear.init(0.0, r4, _rs);
	pf_rsntk.init(0.0, r4, _rs);

	do {
		r = pf_bz.get_plot_pos();

		// calculate axial and azimuthal magnetic components
		bz = calc_bz(r);
		bth = calc_bth(r);
		// calculate twist and shear
		tw = calc_twist(r);
		sh = calc_shear(r);
        // calculate the resonant k value for the given radius
        rsntk = calc_rsntk(r);

		pf_bz.plot(bz);
		pf_bth.plot(bth);
		pf_twist.plot(tw);
		pf_shear.plot(sh);
        pf_rsntk.plot(rsntk);

	} while (!pf_bz.done());
	//////////////////////////////////////////////////////////////////////////

	profiled = true;

	/////////////////////////////////////////////////////////////////////////////////////////////
	if (pf_rx_scenarios) {
		// construct arx, dp and wr profiles for the relaxation regimes
		/////////////////////////////////////////////////////////////////////////////////////////////

		// size the profile vectors
		vpf_rx_a.resize(RX_SCENARIO_CNT);
		vpf_rx_b1.resize(RX_SCENARIO_CNT);
		vpf_rx_dp.resize(RX_SCENARIO_CNT);
		vpf_rx_davp.resize(RX_SCENARIO_CNT);
		vpf_rx_wr.resize(RX_SCENARIO_CNT);
		vpf_rx_q.resize(RX_SCENARIO_CNT);

		// initialise the profiles for the relaxation scenarios
		double rmin, rmax, rs;
		for (int i = 0; i < RX_SCENARIO_CNT; i++) {
			rmin = r2;
			rmax = r4;
			rs = ((rmax-rmin)/rmax)*_rs;

			vpf_rx_a[i].init(rmin, rmax, rs);
			vpf_rx_b1[i].init(rmin, rmax, rs);
			vpf_rx_dp[i].init(rmin, rmax, rs);
			vpf_rx_davp[i].init(rmin, rmax, rs);
			vpf_rx_wr[i].init(rmin, rmax, rs);
			vpf_rx_q[i].init(rmin, rmax, rs);
		}

		// plot the profiles for the relaxation scenarios
		double psi_0r4 = calc_psi(r4), psi_0r = 0.0;
		double k_0r4 = calc_k(r4), k_0r = 0.0;

		CRelaxedLoop rx;

		do {
			r = vpf_rx_a[0].get_plot_pos();
			psi_0r = calc_psi(r);
			k_0r = calc_k(r);

			// Scenario 1: loop has net current.
			///////////////////////////////////////////////////////////////
			rx = calc_rx(r, 0.0, psi_0r4, r4, k_0r4/(psi_0r4*psi_0r4), r4);
			rx.profile(_rs);
			CRelaxedLoop::dump_pfs(r, rx);

			vpf_rx_a[0].plot(rx.get_a1());
			vpf_rx_b1[0].plot(rx.get_b1());
			vpf_rx_dp[0].plot(rx.calc_p(r) - calc_p(r));
			vpf_rx_davp[0].plot(rx.calc_p_av(0.0, r) - rx.calc_p_av(r, r4));
			vpf_rx_wr[0].plot(calc_w(r4) - rx.calc_w(r4));
			vpf_rx_q[0].plot(rx.calc_q(r));
			///////////////////////////////////////////////////////////////

			// For following scenarioes, if loop in fact has net current, all that happens is
			// is that the loop relaxes out to a certain radius.

			// Scenario 2: loop has zero net current due to having a neutralisation surface.
			////////////////////////////////////////////////////////////////////////////////////
			rx = calc_rx(r, r, psi_0r4, r4, k_0r4/(psi_0r4*psi_0r4), r4);
			rx.profile(_rs);
			CRelaxedLoop::dump_pfs(r, rx);

			vpf_rx_a[1].plot(rx.get_a1());
			vpf_rx_b1[1].plot(rx.get_b1());
			vpf_rx_dp[1].plot(rx.calc_p(r) - calc_p(r));
			vpf_rx_davp[1].plot(rx.calc_p_av(0.0, r) - rx.calc_p_av(r, r4));
			vpf_rx_wr[1].plot(calc_w(r4) - rx.calc_w(r4));
			vpf_rx_q[1].plot(rx.calc_q(r));
			////////////////////////////////////////////////////////////////////////////////////

			// Scenario 2.1: same as Scenario 2, except conservation follows expansion radius.
			////////////////////////////////////////////////////////////////////////////////////
			rx = calc_rx(r, r, psi_0r, r, k_0r/(psi_0r*psi_0r), r);
			rx.profile(_rs);
			CRelaxedLoop::dump_pfs(r, rx);

			vpf_rx_a[2].plot(rx.get_a1());
			vpf_rx_b1[2].plot(rx.get_b1());
			vpf_rx_dp[2].plot(rx.calc_p(r) - calc_p(r));
			vpf_rx_davp[2].plot(rx.calc_p_av(0.0, r) - rx.calc_p_av(r, r4));
			vpf_rx_wr[2].plot(calc_w(r) - rx.calc_w(r));
			vpf_rx_q[2].plot(rx.calc_q(r));
			////////////////////////////////////////////////////////////////////////////////////

			// Scenario 3: loop has zero net current due to having a neutralisation layer.
			////////////////////////////////////////////////////////////////////////////////////
			rx = calc_rx((r2/r3)*r, r, psi_0r4, r4, k_0r4/(psi_0r4*psi_0r4), r4);
			rx.profile(_rs);
			CRelaxedLoop::dump_pfs(r, rx);

			vpf_rx_a[3].plot(rx.get_a1());
			vpf_rx_b1[3].plot(rx.get_b1());
			vpf_rx_dp[3].plot(rx.calc_p(r) - calc_p(r));
			vpf_rx_davp[3].plot(rx.calc_p_av(0.0, r) - rx.calc_p_av(r, r4));
			vpf_rx_wr[3].plot(calc_w(r4) - rx.calc_w(r4));
			vpf_rx_q[3].plot(rx.calc_q((r2/r3)*r));
			////////////////////////////////////////////////////////////////////////////////////

			// Scenario 3.1: same as Scenario 3, except conservation follows expansion radius.
			////////////////////////////////////////////////////////////////////////////////////
			rx = calc_rx((r2/r3)*r, r, psi_0r, r, k_0r/(psi_0r*psi_0r), r);
			rx.profile(_rs);
			CRelaxedLoop::dump_pfs(r, rx);

			vpf_rx_a[4].plot(rx.get_a1());
			vpf_rx_b1[4].plot(rx.get_b1());
			vpf_rx_dp[4].plot(rx.calc_p(r) - calc_p(r));
			vpf_rx_davp[4].plot(rx.calc_p_av(0.0, r) - rx.calc_p_av(r, r4));
			vpf_rx_wr[4].plot(calc_w(r) - rx.calc_w(r));
			vpf_rx_q[4].plot(rx.calc_q((r2/r3)*r));
			////////////////////////////////////////////////////////////////////////////////////

		} while (!vpf_rx_a[0].done());

		profiled_rx_scenarios = true;
	} // end of <if (pf_rx_scenarios)> clause
	/////////////////////////////////////////////////////////////////////////////////////////////

	// do the twist and shear integrals
	int_twist(param);
	int_baty_twist(param);
	int_velli_twist();
	int_shear_sq(param);
	int_shear_abs(param);

	profiles_integrated = true;
}


/**
 * Return the relaxed state of this loop as a CRelaxedLoop object.
 * Relaxation is performed according to Scenario 1.0 or 2.1 depending on the value of netcurrent.
 *
 * @param _r the relaxation radius; is usually greater than the theshold loop radius
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. invalid relaxation radius.
 */
CRelaxedLoop CLoop::relax(double _r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (_r < r2 || _r > r4) {
		throw CException("Error! Invalid relaxation radius.");
	}

	double psi_0r = calc_psi(_r);
	double k_0r = calc_k(_r);

	if (netcurrent) {
		// relax loop according to Scenario 1.0
        return calc_rx(_r, 0.0, psi_0r, _r, k_0r/(psi_0r*psi_0r), _r);
	}
	else {
		// relax loop according to Scenario 2.1
		return calc_rx(_r, _r, psi_0r, _r, k_0r/(psi_0r*psi_0r), _r);
	}
}

/**
 * Return the relaxed state of this loop as a CRelaxedLoop object.
 *
 * The radius of the relaxed state is chosen such the loop does
 * not relax beyond q=n. The relaxation radius of the loop is
 * where q reaches its resonant value.
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 */
CRelaxedLoop CLoop::relax(void) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	const double R_MIN = r3;
	const double R_MAX = r4;
	const double R_INC = 0.01;

	double r = R_MIN, r_rx = r3;
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
		r_rx = r3;
	}

	// relax loop
	return relax(r_rx);
}


/**
 * Return the twist at radial coordinate _r.
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. loop has not been profiled.
 */
double CLoop::get_twist(double _r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (!profiled) {
		throw CException("Error! Loop has not been profiled.");
	}

	return pf_twist.get_plot_val(_r);
}

/**
 * Return the average magnetic twist weighted by area for a particular concentric section
 * of the loop - only the following sections are allowed: [0,r1], [0,r2], [0,r3],
 * [0,r4], [r1,r2], [r2,r3] and [r3,r4].
 *
 * @param _rmin the innermost radial boundary of the section (e.g., 0, r1 or r2).
 * @param _rmax the outermost radial boundary of the section (e.g., r1, r2 or r3).
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. loop has not been profiled;
 *        3. minimum radius greater than or equal to maximum radius;
 *        4. invalid radial boundaries.
 */
double CLoop::get_twist_av(double _rmin, double _rmax) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (!profiled) {
		throw CException("Error! Loop has not been profiled.");
	}

	if (_rmin > _rmax) {
		throw CException("Error! Minimum radius greater than maximum radius");
	}

	if (_rmin == _rmax) {
		return 0.0;
	}

	double avtw = 0.0;

	if (0.0 == _rmin && r1 == _rmax) {
		// core
		avtw = twist_0r1 / (PI*r1*r1);
	}
	else if (0.0 == _rmin && r2 == _rmax) {
		// core and outer layer
		avtw = (twist_0r1+twist_r1r2) / (PI*r2*r2);
	}
	else if (0.0 == _rmin && r3 == _rmax) {
		// core, outer layer and current neutralisation layer (i.e., loop)
		avtw = (twist_0r1+twist_r1r2+twist_r2r3) / (PI*r3*r3);
	}
	else if (0.0 == _rmin && r4 == _rmax) {
		// loop and envelope
		avtw = (twist_0r1+twist_r1r2+twist_r2r3+twist_r3r4) / (PI*r4*r4);
	}
	else if (r1 == _rmin && r2 == _rmax) {
		// outer layer
		avtw = twist_r1r2 / ((PI*r2*r2)-(PI*r1*r1));
	}
	else if (r2 == _rmin && r3 == _rmax) {
		// current neutralisation layer
		avtw = netcurrent ? 0.0 : twist_r2r3 / ((PI*r3*r3)-(PI*r2*r2));
	}
	else if (r3 == _rmin && r4 == _rmax) {
		// envelope
		avtw = twist_r3r4 / ((PI*r4*r4)-(PI*r3*r3));
	}
	else {
		throw CException("Error! Invalid radial boundaries.");
	}

	return avtw;
}

/**
 * Return the average magnetic twist (according to Baty 2001) for a particular concentric
 * section of the loop - only the following sections are allowed: [0,r1], [0,r2],
 * [0,r3], [0,r4], [r1,r2], [r2,r3] and [r3,r4].
 *
 * @param _rmin the innermost radial boundary of the section (e.g., 0, r1 or r2).
 * @param _rmax the outermost radial boundary of the section (e.g., r1, r2 or r3).
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. loop has not been profiled;
 *        3. minimum radius greater than or equal to maximum radius;
 *        4. invalid radial boundaries.
 */
double CLoop::get_baty_twist_av(double _rmin, double _rmax) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (!profiled) {
		throw CException("Error! Loop has not been profiled.");
	}

	if (_rmin > _rmax) {
		throw CException("Error! Minimum radius greater than maximum radius");
	}

	if (_rmin == _rmax) {
		return 0.0;
	}

	double avtw = 0.0;

	if (0.0 == _rmin && r1 == _rmax) {
		// core
		avtw = baty_twist_0r1 / r1;
	}
	else if (0.0 == _rmin && r2 == _rmax) {
		// core and outer layer
		avtw = (baty_twist_0r1+baty_twist_r1r2) / r2;
	}
	else if (0.0 == _rmin && r3 == _rmax) {
		// core, outer layer and current neutralisation layer (i.e., loop)
		avtw = (baty_twist_0r1+baty_twist_r1r2+baty_twist_r2r3) / r3;
	}
	else if (0.0 == _rmin && r4 == _rmax) {
		// loop and envelope
		avtw = (baty_twist_0r1+baty_twist_r1r2+baty_twist_r2r3+baty_twist_r3r4) / r4;
	}
	else if (r1 == _rmin && r2 == _rmax) {
		// outer layer
		avtw = baty_twist_r1r2 / (r2-r1);
	}
	else if (r2 == _rmin && r3 == _rmax) {
		// current neutralisation layer
		avtw = netcurrent ? 0.0 : baty_twist_r2r3 / (r3-r2);
	}
	else if (r3 == _rmin && r4 == _rmax) {
		// envelope
		avtw = baty_twist_r3r4 / (r4-r3);
	}
	else {
		throw CException("Error! Invalid radial boundaries.");
	}

	return avtw;
}

/**
 * Return the average magnetic twist (according to Velli et al 1990) for a particular
 * concentric section of the loop - only the following sections are allowed: [0,r1],
 * [0,r2], [0,r3], [0,r4], [r1,r2], [r2,r3] and [r3,r4].
 *
 * @param _rmin the innermost radial boundary of the section (e.g., 0, r1 or r2).
 * @param _rmax the outermost radial boundary of the section (e.g., r1, r2 or r3).
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. loop has not been profiled;
 *        3. minimum radius greater than or equal to maximum radius;
 *        4. invalid radial boundaries.
 */
double CLoop::get_velli_twist_av(double _rmin, double _rmax) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (!profiled) {
		throw CException("Error! Loop has not been profiled.");
	}

	if (_rmin > _rmax) {
		throw CException("Error! Minimum radius greater than maximum radius");
	}

	if (_rmin == _rmax) {
		return 0.0;
	}

	double avtw = 0.0;

	if (0.0 == _rmin && r1 == _rmax) {
		// core
		avtw = lbth_0r1/rbz_0r1;
	}
	else if (0.0 == _rmin && r2 == _rmax) {
		// core and outer layer
		avtw = (lbth_0r1 + lbth_r1r2)/(rbz_0r1 + rbz_r1r2);
	}
	else if (0.0 == _rmin && r3 == _rmax) {
		// core, outer layer and current neutralisation layer (i.e., loop)
		avtw = (lbth_0r1 + lbth_r1r2 + lbth_r2r3)/(rbz_0r1 + rbz_r1r2 + rbz_r2r3);
	}
	else if (0.0 == _rmin && r4 == _rmax) {
		// loop and envelope
		avtw = (lbth_0r1 + lbth_r1r2 + lbth_r2r3 + lbth_r3r4)/(rbz_0r1 + rbz_r1r2 + rbz_r2r3 + rbz_r3r4);
	}
	else if (r1 == _rmin && r2 == _rmax) {
		// outer layer
		avtw = lbth_r1r2/rbz_r1r2;
	}
	else if (r2 == _rmin && r3 == _rmax) {
		// current neutralisation layer
		avtw = netcurrent ? 0.0 : lbth_r2r3/rbz_r2r3;
	}
	else if (r3 == _rmin && r4 == _rmax) {
		// envelope
		avtw = lbth_r3r4/rbz_r3r4;
	}
	else {
		throw CException("Error! Invalid radial boundaries.");
	}

	return avtw;
}

/**
 * Return the shear at radial coordinate _r.
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. loop has not been profiled;
 */
double CLoop::get_shear(double _r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (!profiled) {
		throw CException("Error! Loop has not been profiled.");
	}

	return pf_shear.get_plot_val(_r);
}

/**
 * Return the shear rms for a particular concentric section of the loop
 * - only the following sections are allowed: [0,r1], [0,r2], [0,r3], [0,r4],
 * [r1,r2], [r2,r3] and [r2,r3].
 *
 * @param _rmin the innermost radial boundary of the section (e.g., 0, r1 or r2).
 * @param _rmax the outermost radial boundary of the section (e.g., r1, r2 or r3).
 * @param _nmr if true return the average twist determined numerically, otherwise return
 *             a value that has been calculated analytically.
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. loop has not been profiled;
 *        3. minimum radius greater than or equal to maximum radius;
 *        4. invalid radial boundaries.
 */
double CLoop::get_shear_rms(double _rmin, double _rmax) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (!profiled) {
		throw CException("Error! Loop has not been profiled.");
	}

	if (_rmin > _rmax) {
		throw CException("Error! Minimum radius greater than maximum radius");
	}

	if (_rmin == _rmax) {
		return 0.0;
	}

	double shear_rms = 0.0;

	if (0.0 == _rmin && r1 == _rmax) {
		// core
		shear_rms = sqrt(shear_sq_0r1/(PI*r1*r1));
	}
	else if (0.0 == _rmin && r2 == _rmax) {
		// core and outer layer
		shear_rms = sqrt((shear_sq_0r1+shear_sq_r1r2)/(PI*r2*r2));
	}
	else if (0.0 == _rmin && r3 == _rmax) {
		// loop (i.e., core, outer layer and neutralisation layer)
		shear_rms = sqrt((shear_sq_0r1+shear_sq_r1r2+shear_sq_r2r3)/(PI*r3*r3));
	}
	else if (0.0 == _rmin && r4 == _rmax) {
		// loop and envelope
		shear_rms = sqrt((shear_sq_0r1+shear_sq_r1r2+shear_sq_r2r3+shear_sq_r3r4)/(PI*r4*r4));
	}
	else if (r1 == _rmin && r2 == _rmax) {
		// outer layer
		shear_rms = sqrt(shear_sq_r1r2/((PI*r2*r2)-(PI*r1*r1)));
	}
	else if (r2 == _rmin && r3 == _rmax) {
		// neutralisation layer
		shear_rms = netcurrent ? 0.0 : sqrt(shear_sq_r2r3/((PI*r3*r3)-(PI*r2*r2)));
	}
	else if (r3 == _rmin && r4 == _rmax) {
		// envelope
		shear_rms = sqrt(shear_sq_r3r4/((PI*r4*r4)-(PI*r3*r3)));
	}
	else {
		throw CException("Error! Invalid radial boundaries.");
	}

	return shear_rms;
}

/**
 * Return the shear mean absolute for a particular concentric section of the loop
 * - only the following sections are allowed: [0,r1], [0,r2], [0,r3], [0,r4],
 * [r1,r2], [r2,r3] and [r3,r4].
 *
 * @param _rmin the innermost radial boundary of the section (e.g., 0, r1 or r2).
 * @param _rmax the outermost radial boundary of the section (e.g., r1, r2 or r3).
 * @param _nmr if true return the average twist determined numerically, otherwise return
 *             a value that has been calculated analytically.
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. loop has not been profiled;
 *        3. minimum radius greater than or equal to maximum radius;
 *        4. invalid radial boundaries.
 */
double CLoop::get_shear_mabs(double _rmin, double _rmax) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (!profiled) {
		throw CException("Error! Loop has not been profiled.");
	}

	if (_rmin > _rmax) {
		throw CException("Error! Minimum radius greater than maximum radius");
	}

	if (_rmin == _rmax) {
		return 0.0;
	}

	double shear_mabs = 0.0;

	if (0.0 == _rmin && r1 == _rmax) {
		// core
		shear_mabs = shear_abs_0r1/(PI*r1*r1);
	}
	else if (0.0 == _rmin && r2 == _rmax) {
		// core and outer layer
		shear_mabs = (shear_abs_0r1+shear_abs_r1r2)/(PI*r2*r2);
	}
	else if (0.0 == _rmin && r3 == _rmax) {
		// loop (i.e., core, outer layer and neutralisation layer)
		shear_mabs = (shear_abs_0r1+shear_abs_r1r2+shear_abs_r2r3)/(PI*r3*r3);
	}
	else if (0.0 == _rmin && r4 == _rmax) {
		// loop and envelope
		shear_mabs = (shear_abs_0r1+shear_abs_r1r2+shear_abs_r2r3+shear_abs_r3r4)/(PI*r4*r4);
	}
	else if (r1 == _rmin && r2 == _rmax) {
		// outer layer
		shear_mabs = shear_abs_r1r2/((PI*r2*r2)-(PI*r1*r1));
	}
	else if (r2 == _rmin && r3 == _rmax) {
		// neutralisation layer
		shear_mabs = netcurrent ? 0.0 : shear_abs_r2r3/((PI*r3*r3)-(PI*r2*r2));
	}
	else if (r3 == _rmin && r4 == _rmax) {
		// envelope
		shear_mabs = shear_abs_r3r4/((PI*r4*r4)-(PI*r3*r3));
	}
	else {
		throw CException("Error! Invalid radial boundaries.");
	}

	return shear_mabs;
}


/**
 * The J0 Bessel function.
 *
 * @param _param the value passed to J0.
 * @return the result of calculating J0.
 */
double CLoop::J0(double _param) const {
	return CGSL::J0(_param);
}

/**
 * The J1 Bessel function.
 *
 * @param _param the value passed to J1.
 * @return the result of calculating J1.
 */
double CLoop::J1(double _param) const {
	return CGSL::J1(_param);
}

/**
 * The Y0 Bessel function.
 *
 * @param _param the value passed to Y0.
 * @return the result of calculating Y0.
 */
double CLoop::Y0(double _param) const {
	return CGSL::Y0(_param);
}

/**
 * The Y1 Bessel function.
 *
 * @param _param the value passed to Y1.
 * @return the result of calculating Y1.
 */
double CLoop::Y1(double _param) const {
	return CGSL::Y1(_param);
}

/**
 * The F0 Composite Bessel function.
 *
 * @param _param the value passed to the individual bessel functions.
 * @return the result of calculating F0.
 */
double CLoop::F0(double _param) const {
	return CGSL::J0(_param) + (c2/b2)*CGSL::Y0(_param);
}

/**
 * The F1 Composite Bessel function.
 *
 * @param _param the value passed to the individual bessel functions.
 * @return the result of calculating F1.
 */
double CLoop::F1(double _param) const {
	return CGSL::J1(_param) + (c2/b2)*CGSL::Y1(_param);
}

/**
 * The G0 Composite Bessel function.
 *
 * @param _param the value passed to the individual bessel functions.
 * @return the result of calculating G0.
 */
double CLoop::G0(double _param) const {
	return CGSL::J0(_param) + (c3/b3)*CGSL::Y0(_param);
}

/**
 * The G1 Composite Bessel function.
 *
 * @param _param the value passed to the individual bessel functions.
 * @return the result of calculating G1.
 */
double CLoop::G1(double _param) const {
	return CGSL::J1(_param) + (c3/b3)*CGSL::Y1(_param);
}


/**
 * Stream the loop's properties to a file
 * (profile date not included).
 *
 * @param _ofs the output file stream.
 */
void CLoop::out_props(ofstream& _ofs) const {
	if (dimensions_initialised && field_initialised) {
		_ofs << a1 << WS;                                      //  1. core alpha
		_ofs << a2 << WS;                                      //  2. outer layer alpha
		_ofs << a3 << WS;                                      //  3. neutralisation layer alpha
		                                                       //     is zero if netcurrent is true

		_ofs << b1 << WS;                                      //  4. coefficient term for core field
		_ofs << b2 << WS;                                      //  5. coefficient term for outer layer field
		_ofs << b3 << WS;                                      //  6. coefficient term for neutralisation layer field
                                                               //     is zero if netcurrent is true
		_ofs << b4 << WS;                                      //  7. coefficient term for potential envelope field
		_ofs << c2 << WS;                                      //  8. coefficient term for outer layer field
		_ofs << c3 << WS;                                      //  9. coefficient term for neutralisation layer field
		                                                       //     is zero if netcurrent is true
		_ofs << c4 << WS;                                      // 10. coefficient term for potential envelope field

		_ofs << calc_psi(r1) << WS;                            // 11. axial flux in the core
		_ofs << calc_psi(r2) << WS;                            // 12. axial flux in the core and outer layer
		_ofs << calc_psi(r3) << WS;                            // 13. axial flux in the core, outer layer (and neutralisation layer)
		                                                       //     equal to the previous field if netcurrent is true
		_ofs << calc_psi(r4) << WS;                            // 14. axial flux in the core, outer layer (neutralisation layer)
		                                                       //     and potential envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
	}

	if (profiled) {
		_ofs << pf_bz.get_plot_val(r1) << WS;                  // 15. bz at r1
		_ofs << pf_bth.get_plot_val(r1) << WS;                 // 16. btheta at r1
		_ofs << pf_bz.get_plot_val(r2) << WS;                  // 17. bz at r2
		_ofs << pf_bth.get_plot_val(r2) << WS;                 // 18. btheta at r2
		_ofs << pf_bz.get_plot_val(r3) << WS;                  // 19. bz at r3
		_ofs << pf_bth.get_plot_val(r3) << WS;                 // 20. btheta at r3

		_ofs << pf_twist.get_plot_val(0.0) << WS;              // 21. axial twist
		_ofs << pf_twist.get_plot_val(r1) << WS;               // 22. twist at boundary between core and outer layer
		_ofs << pf_twist.get_plot_val(r2) << WS;               // 23. twist at boundary between outer layer and neutralisation layer
		                                                       //     (potential envelope if netcurrent is true)
		_ofs << pf_twist.get_max_plot_val(0.0,r1) << WS;       // 24. the maximum twist within the core
		_ofs << pf_twist.get_max_plot_val(r1,r2) << WS;        // 25. the maximum twist within the outer layer
		_ofs << pf_twist.get_max_plot_val(r2,r3) << WS;        // 26. the maximum twist within the neutralisation layer
                                                               //     (merely part of the outer layer if netcurrent is true)
		_ofs << pf_twist.get_max_plot_val(r3,r4) << WS;        // 27. the maximum twist within the potential envelope

		_ofs << pf_shear.get_plot_val(0.0) << WS;              // 28. axial shear
		_ofs << pf_shear.get_plot_val(r1) << WS;               // 29. shear at boundary between core and outer layer
		_ofs << pf_shear.get_plot_val(r2) << WS;               // 30. shear at boundary between outer layer and neutralisation layer
		                                                       //     (potential envelope if netcurrent is true)
		_ofs << pf_shear.get_max_plot_val(0.0,r1) << WS;       // 31. the maximum shear within the core
		_ofs << pf_shear.get_max_plot_val(r1,r2) << WS;        // 32. the maximum shear within the outer layer
		_ofs << pf_shear.get_max_plot_val(r2,r3) << WS;        // 33. the maximum shear within the neutralisation layer
		                                                       //     (merely part of the outer layer if netcurrent is true)
		_ofs << pf_shear.get_max_plot_val(r3,r4) << WS;        // 34. the maximum shear within the potential envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
	}

	if (profiles_integrated) {
		// weighted area average twists calculated numerically
		_ofs << get_twist_av(0.0, r1) << WS;                   // 35. core average twist
		_ofs << get_twist_av(r1, r2) << WS;                    // 36. outer layer average twist
		_ofs << get_twist_av(r2, r3) << WS;                    // 37. neutralisation layer average twist
		                                                       //     (merely part of the outer layer if netcurrent is true)
		_ofs << get_twist_av(r3, r4) << WS;                    // 38. potential envelope average twist
		_ofs << get_twist_av(0.0, r2) << WS;                   // 39. average twist across core and outer layer
		_ofs << get_twist_av(0.0, r3) << WS;                   // 40. average twist across core, outer layer
		                                                       //     and neutralisation layer (if netcurrent is false)
		_ofs << get_twist_av(0.0, r4) << WS;                   // 41. average twist across core, outer layer
		                                                       //     (neutralisation layer) and potential envelope

		// average twists calculated numerically according to Baty 2001
		_ofs << get_baty_twist_av(0.0, r1) << WS;              // 42. core average twist
		_ofs << get_baty_twist_av(r1, r2) << WS;               // 43. outer layer average twist
		_ofs << get_baty_twist_av(r2, r3) << WS;               // 44. neutralisation layer average twist
		                                                       //     (merely part of the outer layer if netcurrent is true)
		_ofs << get_baty_twist_av(r3, r4) << WS;               // 45. potential envelope average twist
		_ofs << get_baty_twist_av(0.0, r2) << WS;              // 46. average twist across core and outer layer
		_ofs << get_baty_twist_av(0.0, r3) << WS;              // 47. average twist across core, outer layer
								                               //     and neutralisation layer (if netcurrent is false)
		_ofs << get_baty_twist_av(0.0, r4) << WS;              // 48. average twist across core, outer layer
								                               //     (neutralisation layer) and potential envelope

		// average twists calculated analytically according to Velli et al 1990
		_ofs << get_velli_twist_av(0.0, r1) << WS;             // 49. core average twist
		_ofs << get_velli_twist_av(r1, r2) << WS;              // 50. outer layer average twist
		_ofs << get_velli_twist_av(r2, r3) << WS;              // 51. neutralisation layer average twist
		                                                       //     (merely part of the outer layer if netcurrent is true)
		_ofs << get_velli_twist_av(r3, r4) << WS;              // 52. potential envelope average twist
		_ofs << get_velli_twist_av(0.0, r2) << WS;             // 53. average twist across core and outer layer
		_ofs << get_velli_twist_av(0.0, r3) << WS;             // 54. average twist across core, outer layer
								                               //     and neutralisation layer (if netcurrent is false)
		_ofs << get_velli_twist_av(0.0, r4) << WS;             // 55. average twist across core, outer layer
								                               //     (neutralisation layer) and potential envelope

		// shear rms calculated numerically
		_ofs << get_shear_rms(0.0, r1) << WS;                  // 56. core shear rms
		_ofs << get_shear_rms(r1, r2) << WS;                   // 57. outer layer shear rms
		_ofs << get_shear_rms(r2, r3) << WS;                   // 58. neutralisation layer shear rms
		                                                       //     (merely part of the outer layer if netcurrent is true)
		_ofs << get_shear_rms(r3, r4) << WS;                   // 59. potential envelope shear rms
		_ofs << get_shear_rms(0.0, r2) << WS;                  // 60. shear rms across core and outer layer
		_ofs << get_shear_rms(0.0, r3) << WS;                  // 61. shear rms across core, outer layer and
		                                                       //     neutralisation layer (if netcurrent is false)
		_ofs << get_shear_rms(0.0, r4) << WS;                  // 62. shear rms across core, outer layer
		                                                       //     (neutralisation layer) and potential envelope

		// shear mabs calculated numerically
		_ofs << get_shear_mabs(0.0, r1) << WS;                 // 63. core shear mabs
		_ofs << get_shear_mabs(r1, r2) << WS;                  // 64. outer layer shear mabs
		_ofs << get_shear_mabs(r2, r3) << WS;                  // 65. neutralisation layer shear mabs
		                                                       //     (merely part of the outer layer if netcurrent is true)
		_ofs << get_shear_mabs(r3, r4) << WS;                  // 66. potential envelope shear mabs
		_ofs << get_shear_mabs(0.0, r2) << WS;                 // 67. shear mabs across core and outer layer
		_ofs << get_shear_mabs(0.0, r3) << WS;                 // 68. shear mabs across core, outer layer and
								                               //     neutralisation layer (if netcurrent is false)
		_ofs << get_shear_mabs(0.0, r4) << WS;                 // 69. shear mabs across core, outer layer
								                               //     (neutralisation layer) and potential envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
	}

	if (dimensions_initialised && field_initialised) {
		_ofs << calc_k(r1) << WS;                              // 70. magnetic helicity in core
		_ofs << calc_k(r2) << WS;                              // 71. magnetic helicity in core and outer layer
		_ofs << calc_k(r3) << WS;                              // 72. magnetic helicity in core, outer layer
		                                                       //     (and neutralisation layer)
		_ofs << calc_k(r4) << WS;                              // 73. magnetic helicity in loop and envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
	}

	if (dimensions_initialised && field_initialised) {
		_ofs << calc_w(r1) << WS;                              // 74. magnetic energy in core
		_ofs << calc_w(r2) << WS;                              // 75. magnetic energy in core and outer layer
		_ofs << calc_w(r3) << WS;                              // 76. magnetic energy in core, outer layer
		                                                       //     (and neutralisation layer)
		_ofs << calc_w(r4) << WS;                              // 77. magnetic energy in loop and envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
	}

	if (dimensions_initialised && field_initialised && relaxable) {
		// loop is relaxed according to Scenario 2.1 only if it has zero net current, otherwise
		// loop is simply relaxed out to a certain radius

		// relax loop out to r3 (1.0) and according to Scenario 2.1
		try {
		  CRelaxedLoop rx_loop = relax(r3);
		  _ofs << r3 << WS;                                    // 78. minimum relaxation radius
		  _ofs << calc_w(r3) - rx_loop.calc_w(r3) << WS;       // 79. energy release (rx_r = 1.0)
		  _ofs << rx_loop.get_a1() << WS;                      // 80. relaxation alpha (rx_r = 1.0)
		  _ofs << rx_loop.calc_psi(r3) << WS;                  // 81. axial flux (rx_r = 1.0)
		  _ofs << rx_loop.calc_k(r3) << WS;                    // 82. relaxed helicity (rx_r = 1.0)
		}
		catch (CLoopException& x) {
			_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		}

		// relax loop out to 1.5 and according to Scenario 2.1
		try {
		  CRelaxedLoop rx_loop = relax(1.5);
		  _ofs << "1.5" << WS;                                 // 83. maximum relaxation radius
		  _ofs << calc_w(1.5) - rx_loop.calc_w(1.5) << WS;     // 84. energy release (rx_r = 1.5)
		  _ofs << rx_loop.get_a1() << WS;                      // 85. relaxation alpha (rx_r = 1.5)
		  _ofs << rx_loop.calc_psi(1.5) << WS;                 // 86. axial flux (rx_r = 1.5)
		  _ofs << rx_loop.calc_k(1.5) << WS;                   // 87. relaxed helicity (rx_r = 1.5)
		}
		catch (CLoopException& x) {
		  _ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		}

		// relax loop out to 2.0 and according to Scenario 2.1
		try {
		  CRelaxedLoop rx_loop = relax(2.0);
		  _ofs << "2.0" << WS;                                 // 88. maximum relaxation radius
		  _ofs << calc_w(2.0) - rx_loop.calc_w(2.0) << WS;     // 89. energy release (rx_r = 2.0)
		  _ofs << rx_loop.get_a1() << WS;                      // 90. relaxation alpha (rx_r = 2.0)
		  _ofs << rx_loop.calc_psi(2.0) << WS;                 // 91. axial flux (rx_r = 2.0)
		  _ofs << rx_loop.calc_k(2.0) << WS;                   // 92. relaxed helicity (rx_r = 2.0)
		}
		catch (CLoopException& x) {
		  _ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		}

		// relax loop out to 2.5 and according to Scenario 2.1
		try {
		  CRelaxedLoop rx_loop = relax(2.5);
		  _ofs << "2.5" << WS;                                 // 93. maximum relaxation radius
		  _ofs << calc_w(2.5) - rx_loop.calc_w(2.5) << WS;     // 94. energy release (rx_r = 2.5)
		  _ofs << rx_loop.get_a1() << WS;                      // 95. relaxation alpha (rx_r = 2.5)
		  _ofs << rx_loop.calc_psi(2.5) << WS;                 // 96. axial flux (rx_r = 2.5)
		  _ofs << rx_loop.calc_k(2.5) << WS;                   // 97. relaxed helicity (rx_r = 2.5)
		}
		catch (CLoopException& x) {
		  _ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		}

		// relax loop out to r4 (3.0) and according to Scenario 2.1
		try {
		  CRelaxedLoop rx_loop = relax(r4);
		  _ofs << r4 << WS;                                    // 98. maximum relaxation radius
		  _ofs << calc_w(r4) - rx_loop.calc_w(r4) << WS;       // 99. energy release (rx_r = 3.0)
		  _ofs << rx_loop.get_a1() << WS;                      // 100. relaxation alpha (rx_r = 3.0)
		  _ofs << rx_loop.calc_psi(r4) << WS;                  // 101. axial flux (rx_r = 3.0)
		  _ofs << rx_loop.calc_k(r4) << WS;                    // 102. relaxed helicity (rx_r = 3.0)
		}
		catch (CLoopException& x) {
		  _ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		}

		// relax loop out to where q=1 and according to Scenario 2.1
		try {
		  CRelaxedLoop rx_loop = relax();
		  double rx_r = rx_loop.get_r1();
		  _ofs << rx_r << WS;                                  // 103. q=1 relaxation radius
		  _ofs << calc_w(rx_r) - rx_loop.calc_w(rx_r) << WS;   // 104. energy release (q=1 rx_r)
		  _ofs << rx_loop.get_a1() << WS;                      // 105. relaxation alpha (q=1 rx_r)
		  _ofs << rx_loop.calc_psi(rx_r) << WS;                // 106. axial flux (q=1 rx_r)
		  _ofs << rx_loop.calc_k(rx_r) << WS;                  // 107. relaxed helicity (q=1 rx_r)
		}
		catch (CLoopException& x) {
		  _ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		}
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
	}


    if (profiled) {
		int rsntk_min_mlt = ceil(0.5*l*pf_rsntk.get_min_plot_val(0.0,r3)/PI);
		int rsntk_max_mlt = floor(0.5*l*pf_rsntk.get_max_plot_val(0.0,r3)/PI);
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

		_ofs << rsntk_min_mlt << WS;   // 108. the minimum multiple of 2PI/L within the resonant k profile between 0 and r3.
		_ofs << rsntk_max_mlt << WS;   // 109. the maximum multiple of 2PI/L
		_ofs << rsntk_scr << WS;       // 110. radially weighted resonance score
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
void CLoop::out_field_pfs(ofstream& _ofs) const {
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
        _ofs << pf_rsntk.get_plot_val(i) << WS; // 6. resonant k

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
void CLoop::out_rx_pfs(ofstream& _ofs) const {
	if (!profiled) {
		throw CException("Error! Loop not profiled.");
	}

	unsigned long plot_cnt = vpf_rx_a[0].get_plot_cnt();
	unsigned long i = 0;

	for (i = 0; i < plot_cnt; ++i) {
		double rrx = vpf_rx_a[0].get_plot_pos(i);
		_ofs << rrx << WS;                                   // 1. relaxation radius
		
		for (int j = 0; j < RX_SCENARIO_CNT; j++) {
			double arx = vpf_rx_a[j].get_plot_val(i);
			_ofs << arx << WS;                               // 2. relaxation alpha for when k/(psi*psi) is conserved

			if (CRelaxedLoop::KCRIT_BESSEL_ARG/rrx == arx) {
				_ofs << '0' << WS;                           // 3. relaxed state features helical modes
			}
			else {
				_ofs << '1' << WS;			                 // 3. relaxed state is axisymmetric
			}

			_ofs << vpf_rx_b1[j].get_plot_val(i) << WS;      // 4. b1 coefficient
			_ofs << vpf_rx_dp[j].get_plot_val(i) << WS;      // 5. pressure difference
			_ofs << vpf_rx_davp[j].get_plot_val(i) << WS;    // 6. average pressure difference
			_ofs << vpf_rx_wr[j].get_plot_val(i) << WS;      // 7. energy release
			_ofs << vpf_rx_q[j].get_plot_val(i) << WS;       // 8. q value at relaxation radius
		}

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
ofstream& operator<<(ofstream& _ofs, const CLoop& _lp) {
	_lp.out_props(_ofs);
	return _ofs;
}
