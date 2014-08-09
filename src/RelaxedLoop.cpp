/**
 * @file RelaxedLoop.cpp
 *
 * This file defines the methods of the CRelaxedLoop class.
 */

#include "StdInclude.h"

#include "Exception.h"
#include "GSL.h"
#include "NetBthCalculator.h"
#include "Profile.h"
#include "RelaxedLoop.h"


const double PI = CGSL::PI;
const double CRelaxedLoop::KCRIT_BESSEL_ARG = 3.11;
unsigned int CRelaxedLoop::dump_id = 1;


/**
 * All attributes are zeroed.
 */
CRelaxedLoop::CRelaxedLoop(void): dimensions_initialised(false), envelope_fixed(false), envelope_neutralised(false),
		                          coeffs_initialised(false), coeffs_normalised(false),
		                          field_initialised(false), profiled(false) {
	reset();
}

/**
 * Copy Constructor
 */
CRelaxedLoop::CRelaxedLoop(const CRelaxedLoop& _rx): pf_bz(_rx.pf_bz), pf_bth(_rx.pf_bth),
		                                             pf_twist(_rx.pf_twist), pf_shear(_rx.pf_shear) {
	dimensions_initialised = _rx.dimensions_initialised;
	coeffs_initialised = _rx.coeffs_initialised;
	coeffs_normalised = _rx.coeffs_normalised;
	field_initialised = _rx.field_initialised;
	profiled = _rx.profiled;

	has_nl = _rx.has_nl;
	has_pe = _rx.has_pe;
	is_helical = _rx.is_helical;

	r1 = _rx.r1;
	r2 = _rx.r2;
	r3 = _rx.r3;

	l = _rx.l;

	a1 = _rx.a1;
	a2 = _rx.a2;

	b1 = _rx.b1;
	b2 = _rx.b2;
	c2 = _rx.c2;
	b3 = _rx.b3;
	c3 = _rx.c3;

	envelope_fixed = _rx.envelope_fixed;
	envelope_neutralised = _rx.envelope_neutralised;
	bpe = _rx.bpe;
	cpe = _rx.cpe;
}

/**
 * The assignment operator.
 */
CRelaxedLoop& CRelaxedLoop::operator=(const CRelaxedLoop& _rx) {
	if (this != &_rx) {
		// protect against self-assignment
		dimensions_initialised = _rx.dimensions_initialised;
		coeffs_initialised = _rx.coeffs_initialised;
		coeffs_normalised = _rx.coeffs_normalised;
		field_initialised = _rx.field_initialised;
		profiled = _rx.profiled;

		has_nl = _rx.has_nl;
		has_pe = _rx.has_pe;
		is_helical = _rx.is_helical;

		r1 = _rx.r1;
		r2 = _rx.r2;
		r3 = _rx.r3;

		l = _rx.l;

		a1 = _rx.a1;
		a2 = _rx.a2;

		b1 = _rx.b1;
		b2 = _rx.b2;
		c2 = _rx.c2;
		b3 = _rx.b3;
		c3 = _rx.c3;

		envelope_fixed = _rx.envelope_fixed;
		envelope_neutralised = _rx.envelope_neutralised;
		bpe = _rx.bpe;
		cpe = _rx.cpe;

		pf_bz = _rx.pf_bz;
		pf_bth = _rx.pf_bth;
		pf_twist = _rx.pf_twist;
		pf_shear = _rx.pf_shear;
	}

    return *this;
}

/**
 * Set attributes to zero.
 * @param _full if true all attributes are reset
 */
void CRelaxedLoop::reset(bool _full/* = false*/) {
	if (_full) {
		dimensions_initialised = false;
		envelope_fixed = false;
		envelope_neutralised = false;
		coeffs_initialised = false;
		coeffs_normalised = false;
		field_initialised = false;
		profiled = false;
	}

	if (!dimensions_initialised) {
		r1 = r2 = r3 = 0.0;
		l = 0.0;

		has_nl = false;
		has_pe = false;
	}

	if (!field_initialised) {
		a1 = a2 = 0.0;
		b1 = b2 = b3 = c2 = c3 = 0.0;

		is_helical = false;
	}

	if (!envelope_fixed) {
		bpe = cpe = 0.0;
	}

	if (!profiled) {
		pf_bz.reset();
		pf_bth.reset();
		pf_twist.reset();
		pf_shear.reset();
	}
}


/**
 * Initialises the loop's dimensions - the radial boundaries and length.
 * The loop has net current, so this loop will NOT feature a current neutralisation layer.
 *
 * Other attributes are zeroed.
 *
 * If this function completes successfully the dimensions initialisation flag is
 * set to true.
 *
 * @param _r1 the core radius.
 * @param _r2 the potential envelope radius.
 * @param _l the length of the loop in units of 1 radial unit.
 * @throw CException should any of the following occur:
 *        1. any of the radial boundaries has a zero value;
 *        2. radial boundaries in incorrect numerical order;
 *        3. length is not positive.
 */
void CRelaxedLoop::init_dimensions(double _r1, double _r2, double _l) {
	if (0.0 == _r1 || 0.0 == _r2) {
		throw CException("Error! Zero radial boundary.");
	}

	if (_r1 > _r2) {
		throw CException("Error! Radial boundary out of order (must have 0 < r1 < r2).");
	}

	if (_l <= 0.0) {
		throw CException("Error! Length is not positive.");
	}

	dimensions_initialised = false;
	envelope_fixed = false;
	envelope_neutralised = false;
	coeffs_initialised = false;
	coeffs_normalised = false;
	field_initialised = false;
	profiled = false;

	reset();

	r1 = _r1;
	r2 = _r2;

	has_nl = false;
	has_pe = r1 < r2;

	l = _l;

	dimensions_initialised = true;
}

/**
 * Initialises the loop's dimensions - the radial boundaries and length.
 * The loop has zero net current, so this loop will feature a current neutralisation layer.
 *
 * Other attributes are zeroed.
 *
 * If this function completes successfully the dimensions initialisation flag is
 * set to true.
 *
 * @param _r1 the core radius.
 * @param _r2 the neutralisation layer radius.
 * @param _r3 the potential envelope radius.
 * @param _l the length of the loop in units of 1 radial unit.
 * @throw CException should any of the following occur:
 *        1. any of the radial boundaries has a zero value;
 *        2. radial boundaries in incorrect numerical order;
 *        3. length is not positive.
 */
void CRelaxedLoop::init_dimensions(double _r1, double _r2, double _r3, double _l) {
	if (0.0 == _r1 || 0.0 == _r2 || 0.0 == _r3) {
		throw CException("Error! Zero radial boundary.");
	}

	if (_r1 > _r2 || _r1 > _r3 || _r2 > _r3) {
		throw CException("Error! Radial boundary out of order (must have 0 < r1 < r2 < r3).");
	}

	if (_l <= 0.0) {
		throw CException("Error! Length is not positive.");
	}

	dimensions_initialised = false;
	envelope_fixed = false;
	envelope_neutralised = false;
	coeffs_initialised = false;
	coeffs_normalised = false;
	field_initialised = false;
	profiled = false;

	reset();

	r1 = _r1;
	r2 = _r2;
	r3 = _r3;

	has_nl = true;
	has_pe = r2 < r3;

	l = _l;

	dimensions_initialised = true;
}


/**
 * Fixes the magnetic coefficients of the potential envelope, so that they
 * will not be changed by a subsequent call to init_field.
 *
 * @param _b the b coefficient.
 * @param _c the c coefficient, will be zero if pre-relaxed (i.e., threshold) state
 * had neutralisation layer.
 * @throw CException should any of the following occur:
 *        1. loop dimensions uninitialised.
 */
void CRelaxedLoop::fix_envelope(double _b, double _c/* = 0.0*/) {
	if (!dimensions_initialised) {
		throw CException("Error! Loop dimensions uninitialised.");
	}

	envelope_fixed = false;
	envelope_neutralised = false;
	coeffs_initialised = false;
	coeffs_normalised = false;
	field_initialised = false;
	profiled = false;

	reset();

	bpe = _b;
	cpe = _c;

	if (0.0 == cpe) {
		envelope_neutralised = true;
	}

	envelope_fixed = true;
}

/**
 * Zeros the C coefficient of the potential envelope, so that it has no current.
 * This will not be changed by a subsequent call to init_field.
 *
 * @throw CException should any of the following occur:
 *        1. loop dimensions uninitialised.
 */
void CRelaxedLoop::neutralise_envelope(void) {
	if (!dimensions_initialised) {
		throw CException("Error! Loop dimensions uninitialised.");
	}

	envelope_fixed = false;
	envelope_neutralised = false;
	coeffs_initialised = false;
	coeffs_normalised = false;
	field_initialised = false;
	profiled = false;

	reset();

	cpe = 0.0;

	envelope_neutralised = true;
}

/**
 * This function calculates the azimuthal field at r2 (in units of b1).
 *
 * This function is used in conjunction with the GSL root solver gsl_root_fsolver_brent.
 *
 * @param _an the alpha value for the neutralisation layer.
 * @param _param points to a CNetBthCalculator object
 * @return zero if _an really does neutralise the azimuthal field, otherwise
 * the azimuthal field at r2.
 * @throw CException should any of the following occur:
 *        1. null _param parameter;
 *        2. null loop label.
 */
double CRelaxedLoop::calc_r2bth(double _a2, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CNetBthCalculator* calc = (CNetBthCalculator*) (_param);
	if (NULL == calc) {
		throw CException("Error! NULL relaxed CNetBthCalculator object.");
	}

	return calc->netbth(_a2);
}

/**
 * This method calculates the alpha value (a2) for a loop's current neutralisation layer (r1-r2).
 * The neutralisation layer current must be such that it causes the azimuthal field to vanish at r2.
 * The value of a2 with the smallest absolute value is returned.
 *
 * @return the alpha value for the neutralisation layer
 * @throw CLoopException should any of the following occur:
 *        1. cannot find neutralisation alpha.
 */
double CRelaxedLoop::calc_nlalpha(void) {
	CNetBthCalculator calc(r1, r2, a1);

	void* param = &calc;

	double a2 = 0.0;

	try {
		a2 = CGSL::find_root(calc_r2bth, 50.0, 0.01, (a1 <= 0 ? 1 : -1), param);
	}
	catch (CGSLException& x) {
		throw CLoopException(a1, 0.0, "error finding relaxed neutralisation alpha", x);
	}

	return a2;
}

/**
 * Calculates the magnetic field coefficients.
 *
 * This function is private, it is called by init_field.
 *
 * @param _psi the axial flux used to set the b1 coefficient.
 * @param _psi_r the radial extent of axial flux used to set the b1 coefficient.
 */
void CRelaxedLoop::init_field_coeffs(double _psi, double _psi_r) {
	double aa1 = fabs(a1);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);

	coeffs_initialised = false;
	coeffs_normalised = false;

	if (has_nl) {
		// loop has neutralisation layer
		if (0.0 == a1) {
			// potential loop
			// calculate what a2 needs to be in order for the azimuthal field to vanish at r2
			a2 = 0.0;

			b2 = 1.0;
			c2 = 0.0;
		}
		else {
			// non-potential loop
			// calculate what a2 needs to be in order for the azimuthal field to vanish at r2
			// if a1 non-zero then so is a2
			a2 = calc_nlalpha();

			double aa2 = fabs(a2);
			double sa2 = (a2 < 0 ? -1.0 : 1.0);
			double sa1a2 = sa1*sa2;

			double delta = 2.0/(PI*aa2*r1);
			b2 = (sa1a2*J1(aa1*r1)*Y0(aa2*r1) - J0(aa1*r1)*Y1(aa2*r1))/delta;
			c2 = (J0(aa1*r1)*J1(aa2*r1) - sa1a2*J1(aa1*r1)*J0(aa2*r1))/delta;
		}

		if (has_pe) {
			// zero net current loop has potential envelope
			if (envelope_fixed) {
				// fix the envelope so that its field coefficients are the
				// same as they were for the pre-relaxed state
				b3 = bpe;
				c3 = cpe;
			}
			else {
				double aa2 = fabs(a2);
				b3 = 0.0 == aa2 ? 1.0 : b2*F0(aa2*r2);
				c3 = 0.0;
			}
		}
	}
	else if (has_pe) {
		// net current loop has potential envelope
		if (envelope_fixed) {
			// fix the envelope so that its field coefficients are the
			// same as they were for the pre-relaxed state
			b2 = bpe;
			c2 = cpe;
		}
		else {
			b2 = 0.0 == aa1 ? 1.0 : J0(aa1*r1);
			c2 = 0.0;
		}
	}

	coeffs_initialised = true;

	// normalise coefficients such that axial flux is conserved
	b1 = 1.0; // normalise b1 so that calc_psi will return a flux in units of b1

	double rlp = has_nl ? r2 : r1;
	if (envelope_fixed && _psi_r > rlp) {
		double psi_lp = calc_psi(rlp);
		double psi_pe = calc_psi(_psi_r) - psi_lp;
		b1 = (_psi - psi_pe)/psi_lp;
	}
	else {
		b1 = _psi/calc_psi(_psi_r);
	}

	if (has_nl) {
		b2 = b2*b1;
		c2 = c2*b1;

		if (has_pe && !envelope_fixed) {
			b3 = b3*b1;
			c3 = c3*b1;
		}
	}
	else if (has_pe && !envelope_fixed) {
		b2 = b2*b1;
		c2 = c2*b1;
	}

	coeffs_normalised = true;
}

/**
 * Initialises the loop object's alpha values and Bessel function coefficients.
 * The details of any other features of the relaxed loop (e.g., neutralisation
 * layer and/or potential envelope) are also calculated by this function.
 *
 * If this function completes successfully the field initialisation flag is
 * set to true.
 *
 * @param _psi the axial flux used to set the b1 coefficient.
 * @param _psi_r the radial extent of axial flux used to set the b1 coefficient.
 * @param _a1 the alpha value within the core.
 * @param _helical_test if true test if relaxed state features helical modes.
 *
 * @throw CException should any of the following occur:
 *        1. loop dimensions uninitialised;
 *        2. zero _psi_r.
 */
void CRelaxedLoop::init_field(double _psi, double _psi_r, double _a1, bool _helical_test/* = true*/) {
	if (!dimensions_initialised) {
		throw CException("Error! Loop dimensions uninitialised.");
	}

	if (0.0 == _psi_r) {
		throw CException("Error! Zero _psi_r.");
	}

	coeffs_initialised = false;
	coeffs_normalised = false;
	field_initialised = false;
	profiled = false;

	reset();

	a1 = _a1;

	// initialise the field coefficents (e.g., b1, b2, c2)
	init_field_coeffs(_psi, _psi_r);

	if (_helical_test) {
		// determine if relaxed state features helical modes
		CRelaxedLoop helical(*this);
		helical.init_field(_psi, _psi_r, KCRIT_BESSEL_ARG/r1, false);

		is_helical = (calc_k(r1) >= helical.calc_k(r1));
		if (is_helical) {
			a1 = KCRIT_BESSEL_ARG/r1;
			init_field_coeffs(_psi, _psi_r);
		}
	}

	field_initialised = coeffs_initialised && coeffs_normalised;
	if (has_nl) {
		field_initialised = field_initialised && envelope_fixed;
	}
}


/**
 * Calculate axial component of the magnetic field for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @return the axial component of the magnetic field
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds.
 */
double CRelaxedLoop::calc_bz(double _r) const {
	if (!validate_radius(_r)) {
		throw CException("Error! Radial value out of bounds.");
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double bz = 0.0;

	double r = _r;
	if (0.0 == r || 0.0 == aa1) {
		return b1;
	}

	if (r <= r1) {
		// core
		bz = b1*J0(aa1*r);
	}
	else {
		if (has_nl) {
			if (r <= r2) {
				// neutralisation layer
				bz = b2*J0(aa2*r) + c2*Y0(aa2*r);
			}
			else {
				// potential envelope
				bz = b3;
			}
		}
		else {
			// potential envelope
			bz = b2;
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
double CRelaxedLoop::calc_bz(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CRelaxedLoop* loop = (CRelaxedLoop*) (_param);
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
double CRelaxedLoop::calc_bth(double _r) const {
	if (!validate_radius(_r)) {
		throw CException("Error! Radial value out of bounds.");
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa2 = (a2 < 0 ? -1.0 : 1.0);
	double bth = 0.0;

	double r = _r;
	if (0.0 == r || 0.0 == aa1) {
		return 0.0;
	}

	if (r <= r1) {
		// core
		bth = sa1*b1*J1(aa1*r);
	}
	else {
		if (has_nl) {
			if (r <= r2) {
				// neutralisation layer
				bth = sa2*(b2*J1(aa2*r) + c2*Y1(aa2*r));
			}
			else  {
				// potential envelope
				//bth = sa2*b2*F1(aa2*r2)*r2/r;
				bth = 0.0; // loop is neutralised
			}
		}
		else {
			// potential envelope
			if (envelope_neutralised) {
				bth = 0.0;
			}
			else {
				bth = sa1*b1*J1(aa1*r1)*r1/r;
			}
		}
	}

	return bth;
}

/**
 * Calculate q value of the magnetic field for a given radial coordinate.
 *
 * @param _r the radial coordinate
 * @return the q value of the magnetic field
 * @throw CException should any of the following occur:
 *        1. radial value out of bounds;
 *        2. zero bth
 */
double CRelaxedLoop::calc_q(double _r) const {
	if (!validate_radius(_r)) {
		throw CException("Error! Radial value out of bounds.");
	}

	double bz = calc_bz(_r);
	double bth = calc_bth(_r);

	if (0.0 == bth) {
		throw CException("Error! Zero bth.");
	}

	return (2*PI*_r*bz)/(l*bth);
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
double CRelaxedLoop::calc_bth(double _r, void* _param) {
	if (NULL == _param) {
		throw CException("Error! NULL _param parameter.");
	}

	CRelaxedLoop* loop = (CRelaxedLoop*) (_param);
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
double CRelaxedLoop::calc_p(double _r) const {
	double bz = calc_bz(_r);
	double bth = calc_bth(_r);

	return (bz*bz + bth*bth)/2.0;
}

/**
 * Get the average pressure of the magnetic field over the specified radial range.
 *
 * @param _rmin the minimum radial coordinate
 * @param _rmax the maximum radial coordinate
 * @return the average magnetic pressure over _rmin to _rmax
 *
 * @throw CException should any of the following occur:
 *        1. minimum radius greater than or equal to maximum radius.
 */
double CRelaxedLoop::calc_p_av(double _rmin, double _rmax) const {
	if (_rmin > _rmax) {
		throw CException("Error! Minimum radius greater than maximum radius");
	}

	if (_rmin == _rmax) {
		return 0.0;
	}

	return (calc_w(_rmax) - calc_w(_rmin))/(l*PI*(_rmax*_rmax - _rmin*_rmin));
}

/**
 * Returns true if the radial coordinate is within the correct range
 * for this relaxed loop.
 *
 * @param _r the radial coordinate
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field coefficients uninitialised;
 */
bool CRelaxedLoop::validate_radius(double _r) const {
	if (!dimensions_initialised || !coeffs_initialised) {
		throw CException("Error! dimensions and/or field coefficients uninitialised.");
	}

	bool valid = false;

	if (_r >= 0.0) {
		double rmax = r1;

		if (has_nl) {
			rmax = has_pe ? r3 : r2;
		}
		else {
			rmax = has_pe ? r2 : r1;
		}

		valid = (_r <= rmax);
	}

	return valid;
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
double CRelaxedLoop::calc_psi(double _r) const {
	if (!dimensions_initialised || !coeffs_initialised) {
		throw CException("Error! dimensions and/or field coefficients uninitialised.");
	}

	if (!validate_radius(_r)) {
		throw CException("Error! Radial value out of bounds.");
	}

	if (0.0 == b1) {
		throw CException("Error! b1 is uninitialised.");
	}

	if (0.0 == _r) {
		return 0.0;
	}

	if (0.0 == a1) {
		// if potential core then everywhere is potential
		return PI*b1*_r*_r;
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa2 = (a2 < 0 ? -1.0 : 1.0);
	double sa1a2 = sa1*sa2;

	double psi = 0.0;

	do {
		double r = _r < r1 ? _r : r1;

		// core
		psi = 2.0*PI*b1*r*J1(aa1*r)/aa1;

		if (_r <= r1) {
			break;
		}

		r = _r < r2 ? _r : r2;

		if (has_nl) {
			// neutralisation layer
			psi += ((2.0*PI*b2*r)/aa2)*F1(aa2*r) - (2.0*PI*b1*r1*J1(aa1*r1))*(sa1a2/aa2);
		}
		else {
			// potential envelope
			psi += (PI*(r*r - r1*r1))*b2;
		}

		if (_r <= r2) {
			break;
		}

		r = _r < r3 ? _r : r3;

		// has neutralisation layer and potential envelope
		psi += PI*b3*(r*r - r2*r2);

	} while (0);

	return psi;
}

/**
 * Return the helicity between 0 and _r.
 *
 * @param _r the radial coordinate that is the upper bound for the helicity calculation
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field coefficients uninitialised;
 *        2. b1 is uninitialised;
 *        3. radial value out of bounds.
 */
double CRelaxedLoop::calc_k(double _r) const {
	if (!dimensions_initialised || !coeffs_initialised) {
		throw CException("Error! dimensions and/or field coefficients uninitialised.");
	}

	if (0.0 == b1) {
		throw CException("Error! b1 is uninitialised.");
	}

	if (!validate_radius(_r)) {
		throw CException("Error! Radial value out of bounds.");
	}

	if (0.0 == _r || 0.0 == a1) {
		// if potential core everywhere is potential
		return 0.0;
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);
	double sa1 = (a1 < 0 ? -1.0 : 1.0);
	double sa2 = (a2 < 0 ? -1.0 : 1.0);
	double sa1a2 = sa1*sa2;

	double uk = 0.0, k = 0.0;

	do {
		double r = _r < r1 ? _r : r1;

		// core
		uk = (2.0*PI*l*b1*b1/aa1)*(r*r*(J0(aa1*r)*J0(aa1*r) + J1(aa1*r)*J1(aa1*r)) - (2.0*(r/aa1))*J0(aa1*r)*J1(aa1*r));

		k = sa1*uk;

		if (_r <= r1) {
			break;
		}

		r = _r < r2 ? _r : r2;

		if (has_nl) {
			// neutralisation layer
			uk = (2.0*PI*l*b2*b2/aa2)*(r*r*(F0(aa2*r)*F0(aa2*r) + F1(aa2*r)*F1(aa2*r)) - 2.0*(r/aa2)*F0(aa2*r)*F1(aa2*r))
				 -
				 (2.0*PI*l*b2*b2/aa2)*(r1*r1*(F0(aa2*r1)*F0(aa2*r1) + F1(aa2*r1)*F1(aa2*r1)) - 2.0*(r1/aa2)*F0(aa2*r1)*F1(aa2*r1))
				 +
				 (4.0*PI*l*b2/aa2)*(F0(aa2*r1) - F0(aa2*r))*b1*r1*J1(aa1*r1)*((1.0/aa1) - (sa1a2/aa2));

			k += sa2*uk;
		}
		else {
			// potential envelope
			if (envelope_neutralised) {
				uk = 0.0;
			}
			else {
				uk = l*b2*b1*J1(aa1*r1)*r1*PI*(r*r - r1*r1) + 2.0*PI*l*b1*J1(aa1*r1)*r1*r1*(2*((b1*J1(aa1*r1))/aa1) - b2*r1)*log(r/r1);
			}

			k += sa1*uk;
		}

		if (_r <= r2) {
			break;
		}

		r = _r < r3 ? _r : r3;

		// neutralisation layer and potential envelope
		uk = 0.0;
		// if a2 did not in fact have a neutralising function
		// uk would be 2.0*l*b2*F1(aa2*r2)*r2*((calc_psi(r2) - PI*b3*r2*r2)*log(r/r2) + PI*b3*(r*r - r2*r2)/2.0);
		k += sa2*uk;

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
double CRelaxedLoop::calc_w(double _r) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	if (!validate_radius(_r)) {
		throw CException("Error! Radial value out of bounds.");
	}

	if (0.0 == _r) {
		return 0.0;
	}

	if (0.0 == a1) {
		// if potential core then everywhere is potential
		return (l*PI*b1*b1*_r*_r)/2.0;
	}

	double aa1 = fabs(a1);
	double aa2 = fabs(a2);

	double w = 0.0;

	do {
		double r = _r < r1 ? _r : r1;

		// core
		if (is_helical) {
			double psi_0r1 = calc_psi(r1);
			double ak_0r1 = fabs(calc_k(r1));
			w = (aa1/2.0)*(ak_0r1 + (l*J0(aa1*r1)*psi_0r1*psi_0r1)/(2.0*PI*r1*J1(aa1*r1)));
		}
		else {
			w = b1*b1*(r*r*(J0(aa1*r)*J0(aa1*r) + J1(aa1*r)*J1(aa1*r)) - (r/aa1)*J0(aa1*r)*J1(aa1*r));
		}

		if (_r <= r1) {
			break;
		}

		r = _r < r2 ? _r : r2;

		if (has_nl) {
			// neutralisation layer
			w += b2*b2*(r*r*(F0(aa2*r)*F0(aa2*r) + F1(aa2*r)*F1(aa2*r)) - (r/aa2)*F0(aa2*r)*F1(aa2*r)
						-
						r1*r1*(F0(aa2*r1)*F0(aa2*r1) + F1(aa2*r1)*F1(aa2*r1)) + (r1/aa2)*F0(aa2*r1)*F1(aa2*r1));
		}
		else {
			// potential envelope
			if (envelope_neutralised) {
				w += ((b2*b2)/2.0)*(r*r - r1*r1);
			}
			else {
				w += ((b2*b2)/2.0)*(r*r - r1*r1) + (r1*r1*b1*b1*J1(aa1*r1)*J1(aa1*r1)*log(r/r1));
			}
		}

		if (_r <= r2) {
			break;
		}

		r = _r < r3 ? _r : r3;

		// has neutralisation layer and potential envelope
		w += (((b3*b3)/2.0)*(r*r - r2*r2));
		// if in fact a2 did not have a neutralising function
		// w would be (((b3*b3)/2.0)*(r*r - r2*r2)) + (r2*r2*b2*b2*F1(aa2*r2)*F1(aa2*r2)*log(r/r2));

	} while (0);


	return l*PI*w;
}


/**
 * Calculates four radial profiles - magnetic axial field, magnetic azimuthal field,
 * magnetic twist ((l*bth)/(r*bz)) and shear (bth/bz).
 *
 * If this function completes successfully the profiled flag is set to true.
 *
 * @param _rs the profile resolution i.e., the number of points that make up the profile.
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised;
 *        2. divide by zero error.
 */
void CRelaxedLoop::profile(unsigned long _rs) {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	profiled = false;
	reset();

	double r, bz, bth, twist, shear;

	double rmax = r1;
	if (has_pe) {
		rmax = has_nl ? r3 : r2;
	}

	pf_bz.init(0.0, rmax, _rs);
	pf_bth.init(0.0, rmax, _rs);
	pf_twist.init(0.0, rmax, _rs);
	pf_shear.init(0.0, rmax, _rs);

	do {
		r = pf_bz.get_plot_pos();

		// calculate axial and azimuthal magnetic components
		bz = calc_bz(r);
		bth = calc_bth(r);

		// calculate twist and shear
		if (0.0 != r) {
			if (0.0 == bz) {
				throw CException("Error! Divide by zero error.");
			}

			twist = (l*bth) / (r*bz);
			shear = bth / bz;
		}
		else {
			twist = l*fabs(a1)/2.0;
			if (a1 < 0.0) {
				twist = -1.0*twist;
			}
			shear = 0.0;
		}

		pf_bz.plot(bz);
		pf_bth.plot(bth);
		pf_twist.plot(twist);
		pf_shear.plot(shear);

	} while (!pf_bz.done());

	profiled = true;
}


/**
 * Return the relaxation alpha.
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised.
 */
double CRelaxedLoop::get_a1(void) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	return a1;
}

/**
 * Return the relaxation b1.
 *
 * @throw CException should any of the following occur:
 *        1. dimensions and/or field uninitialised.
 */
double CRelaxedLoop::get_b1(void) const {
	if (!dimensions_initialised || !field_initialised) {
		throw CException("Error! Dimensions and/or field uninitialised.");
	}

	return b1;
}

/**
 * Return the relaxation radius.
 *
 * @throw CException should any of the following occur:
 *        1. loop is uninitialised.
 */
double CRelaxedLoop::get_r1(void) const {
	if (!dimensions_initialised) {
		throw CException("Error! Loop uninitialised.");
	}

	return r1;
}

/**
 * Return the neutralisation layer or potential envelope radius.
 *
 * @throw CException should any of the following occur:
 *        1. loop is uninitialised;
 *        2. loop does not have a neutralisation layer or a potential envelope.
 */
double CRelaxedLoop::get_r2(void) const {
	if (!dimensions_initialised) {
		throw CException("Error! Loop uninitialised.");
	}

	if (!has_nl && !has_pe) {
		throw CException("Error! Loop does not have a neutralisation layer or a potential envelope.");
	}

	return r2;
}


/**
 * The J0 Bessel function.
 *
 * @param _param the value passed to J0.
 * @return the result of calculating J0.
 */
double CRelaxedLoop::J0(double _param) const {
	return CGSL::J0(_param);
}

/**
 * The J1 Bessel function.
 *
 * @param _param the value passed to J1.
 * @return the result of calculating J1.
 */
double CRelaxedLoop::J1(double _param) const {
	return CGSL::J1(_param);
}

/**
 * The Y0 Bessel function.
 *
 * @param _param the value passed to Y0.
 * @return the result of calculating Y0.
 */
double CRelaxedLoop::Y0(double _param) const {
	return CGSL::Y0(_param);
}

/**
 * The Y1 Bessel function.
 *
 * @param _param the value passed to Y1.
 * @return the result of calculating Y1.
 */
double CRelaxedLoop::Y1(double _param) const {
	return CGSL::Y1(_param);
}

/**
 * The F0 Composite Bessel function.
 *
 * @param _param the value passed to the individual bessel functions.
 * @return the result of calculating F0.
 */
double CRelaxedLoop::F0(double _param) const {
	return CGSL::J0(_param) + (c2/b2)*CGSL::Y0(_param);
}

/**
 * The F1 Composite Bessel function.
 *
 * @param _param the value passed to the individual bessel functions.
 * @return the result of calculating F1.
 */
double CRelaxedLoop::F1(double _param) const {
	return CGSL::J1(_param) + (c2/b2)*CGSL::Y1(_param);
}


/**
 * Stream the loop's properties to a file
 * (profile date not included).
 *
 * @param _ofs the output file stream.
 */
void CRelaxedLoop::out_props(ofstream& _ofs) const {
	if (dimensions_initialised && field_initialised) {
		_ofs << a1 << WS;                                 //  1. core alpha
		_ofs << a2 << WS;                                 //  2. neutralisation layer alpha

		_ofs << has_nl << WS;                             //  3. net current flag
		_ofs << has_pe << WS;                             //  4. potential envelope flag
		_ofs << is_helical << WS;                         //  5. helical flag

		_ofs << b1 << WS;                                 //  6. coefficient term for core field
		_ofs << b2 << WS;                                 //  7. coefficient term for neutralisation layer field
		_ofs << b3 << WS;                             	  //  8. coefficient term for potential envelope field
		_ofs << c2 << WS;                                 //  9. coefficient term for neutralisation layer field
		_ofs << c3 << WS;                             	  // 10. coefficient term for potential envelope field

		_ofs << calc_psi(r1) << WS;                       // 11. axial flux in the core
		_ofs << calc_psi(r2) << WS;                       // 12. axial flux in the core and neutralisation layer
		_ofs << calc_psi(r3) << WS;                       // 13. axial flux in the core, neutralisation layer and potential envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS << 0.0 << WS;
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS;
	}

	if (dimensions_initialised && field_initialised) {
		_ofs << calc_k(r1) << WS;                          // 14. magnetic helicity in core
		_ofs << calc_k(r2) << WS;                          // 15. magnetic helicity in core and neutralisation layer
		_ofs << calc_k(r3) << WS;                          // 16. magnetic helicity in core, neutralisation layer and potential envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS;
	}

	if (dimensions_initialised && field_initialised) {
		_ofs << calc_w(r1) << WS;                          // 17. magnetic energy in core
		_ofs << calc_w(r2) << WS;                          // 18. magnetic energy in core and neutralisation layer
		_ofs << calc_w(r3) << WS;                          // 19. magnetic energy in core, neutralisation layer and potential envelope
	}
	else {
		_ofs << 0.0 << WS << 0.0 << WS << 0.0 << WS ;
	}

	_ofs << NL;
}

/**
 * Stream the loop's profile data to a file.
 *
 * @param _ofs the output file stream.
 * @throw CException should any of the following occur:
 *        1. loop not profiled.
 */
void CRelaxedLoop::out_pfs(ofstream& _ofs) const {
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
		_ofs << NL;
	}
}


/**
 * Dump the field profiles to a file with a specified name.
 *
 * @throw CException should any of the following occur:
 *        1. cannot open dump file.
 */
void CRelaxedLoop::dump_pfs(double _r, const CRelaxedLoop& _rx) {
	if (!(1.0 == _r || 1.5 == _r || 2.0 == _r || 2.5 == _r || 3.0 == _r)) {
		return;
	}

	stringstream fn(stringstream::in | stringstream::out);

	fn << "dump_pf_" << CRelaxedLoop::dump_id << ".txt";

	ofstream of_dump(fn.str().data());
	if (!of_dump) {
		throw CException("Error! Cannot open dump file.");
	}

	try {
		_rx.out_pfs(of_dump);
	}

	catch (CException& x) {
		of_dump.clear();
		of_dump.close();
		throw x;
	}

	if (of_dump.is_open()) {
		of_dump.clear();
		of_dump.close();
	}

	CRelaxedLoop::dump_id++;
}


/**
 * Stream the properties of the loop to a file stream.
 *
 * This function is a friend of the CRelaxedLoop class.
 *
 * @param _ofs the output file stream.
 * @param _lp the loop object.
 * @return reference to output file.
 */
ofstream& operator<<(ofstream& _ofs, const CRelaxedLoop& _lp) {
	_lp.out_props(_ofs);
	return _ofs;
}
