/**
 * @file NetBthCalculator.cpp
 *
 * This file defines the methods of the CNetBthCalculator class.
 */

#include "StdInclude.h"

#include "Exception.h"
#include "GSL.h"
#include "NetBthCalculator.h"


const double PI = CGSL::PI;


/**
 * Construct a CNetBthCalculator object for a loop that has a single outer layer.
 *
 * @param _r1 the core radius
 * @param _r2 the outer layer radius
 * @param _a1 the alpha value within the core (0-r1)
 */
CNetBthCalculator::CNetBthCalculator(double _r1, double _r2, double _a1) {
	double_outer_layer = false;

	r1 = _r1;
	r2 = _r2;

	a1 = _a1;
}

/**
 * Construct a CNetBthCalculator object for a loop that has a double outer layer.
 *
 * @param _r1 the core radius
 * @param _r2 the first outer layer radius
 * @param _r3 the second outer layer radius
 * @param _a1 the alpha value within the core (0-r1)
 * @param _a2 the alpha value within the first outer layer (r1-r2)
 * @param _b2 the b2 magnetic field coefficient in the first outer layer (r1-r2)
 * @param _c2 the c2 magnetic field coefficient in the first outer layer (r1-r2)
 */
CNetBthCalculator::CNetBthCalculator(double _r1, double _r2, double _r3,
		                             double _a1, double _a2,
		                             double _b2, double _c2) {
	double_outer_layer = true;

	r1 = _r1;
	r2 = _r2;
	r3 = _r3;

	a1 = _a1;
	a2 = _a2;

	// normalise the magnetic field coefficient in the core (0-r1)
	b1 = 1.0;
	b2 = _b2;
	c2 = _c2;
}

/**
 * This function calculates the azimuthal field at r2 if double_outer_layer is false
 * and r3 if double_outer_layer is true.
 *
 * @param _a the alpha value for the neutralisation layer.
 *
 * @return zero if _a really does neutralise the azimuthal field, otherwise
 * the azimuthal field at r3 (double_outer_layer=true) or r2 (double_outer_layer=false).
 *
 * @throw CException should any of the following occur:
 *        1. any of the radial boundaries has a zero value;
 *        2. radial boundaries in incorrect numerical order;
 *        3. b2 is zero.
 */
double CNetBthCalculator::netbth(double _a) const {
	double netbth = 0.0;

	if (double_outer_layer) {

		if (0.0 == r1 || 0.0 == r2 || 0.0 == r3) {
			throw CException("Error! Zero radial boundary.");
		}

		if (r1 > r2 || r1 > r3 || r2 > r3) {
			throw CException("Error! Radial boundary out of order (must have 0 < r1 <= r2 <= r3).");
		}

		if (0.0 == b2) {
			throw CException("Error! Zero b2.");
		}

		double aa1 = fabs(a1);
		double aa2 = fabs(a2);
		double aa3 = fabs(_a);

		double sa1 = (a1 < 0 ? -1.0 : 1.0);
		double sa2 = (a2 < 0 ? -1.0 : 1.0);
		double sa3 = (_a < 0 ? -1.0 : 1.0);
		double sa1a3 = sa1*sa3;
		double sa2a3 = sa2*sa3;

		if (0.0 == a1 && 0.0 == a2) {
			// potential loop
			if (0.0 == _a) {
				netbth = 0.0;
			}
			else {
				// can't use coefficients
				netbth = sa3*0.5*PI*aa3*r2*(J1(aa3*r2)*Y1(aa3*r3) - Y1(aa3*r2)*J1(aa3*r3));
			}
		}
		else if (0.0 == a2) {
			// potential outer layer
			if (0.0 == _a) {
				netbth = sa1*b1*J1(aa1*r1)*r1/r3;
			}
			else {
				netbth = sa3*0.5*PI*aa3*r2*(J1(aa3*r3)*(sa1a3*(r1/r2)*b1*J1(aa1*r1)*Y0(aa3*r2) - b2*Y1(aa3*r2))
										    +
										    Y1(aa3*r3)*(b2*J1(aa3*r2) - sa1a3*(r1/r2)*b1*J1(aa1*r1)*J0(aa3*r2)));
			}
		}
		else {
			// potential core or entirely non-potential loop
			if (0.0 == _a) {
				netbth = sa2*b2*F1(aa2*r2)*r2/r3;
			}
			else {
				netbth = sa3*0.5*PI*aa3*b2*r2*(J1(aa3*r3)*(sa2a3*F1(aa2*r2)*Y0(aa3*r2) - F0(aa2*r2)*Y1(aa3*r2))
											   +
											   Y1(aa3*r3)*(F0(aa2*r2)*J1(aa3*r2) - sa2a3*F1(aa2*r2)*J0(aa3*r2)));
			}
		}
	}
	else {
		// single outer layer
		if (0.0 == r1 || 0.0 == r2) {
			throw CException("Error! Zero radial boundary.");
		}

		if (r1 > r2) {
			throw CException("Error! Radial boundary out of order (must have 0 < r1 <= r2).");
		}

		double aa1 = fabs(a1);
		double aa2 = fabs(_a);
		double sa1 = (a1 < 0 ? -1.0 : 1.0);
		double sa2 = (_a < 0 ? -1.0 : 1.0);
		double sa1a2 = sa1*sa2;

		if (0.0 == a1) {
			// potential loop
			if (0.0 == _a) {
				netbth = 0.0;
			}
			else {
				// can't use coefficients
				netbth = sa2*0.5*PI*aa2*r1*(J1(aa2*r1)*Y1(aa2*r2) - Y1(aa2*r1)*J1(aa2*r2));
			}
		}
		else {
			// non-potential loop
			if (0.0 == _a) {
				netbth = sa1*J1(aa1*r1)*r1/r2;
			}
			else {
				netbth = sa2*0.5*PI*aa2*r1*(J1(aa2*r2)*(sa1a2*Y0(aa2*r1)*J1(aa1*r1) - J0(aa1*r1)*Y1(aa2*r1))
										    +
										    Y1(aa2*r2)*(J0(aa1*r1)*J1(aa2*r1) - sa1a2*J0(aa2*r1)*J1(aa1*r1)));
			}
		}
	} // end of <if (double_outer_layer)> else clause

	return netbth;
}

/**
 * The J0 Bessel function.
 *
 * @param _param the value passed to J0.
 * @return the result of calculating J0.
 */
double CNetBthCalculator::J0(double _param) const {
	return CGSL::J0(_param);
}

/**
 * The J1 Bessel function.
 *
 * @param _param the value passed to J1.
 * @return the result of calculating J1.
 */
double CNetBthCalculator::J1(double _param) const {
	return CGSL::J1(_param);
}

/**
 * The Y0 Bessel function.
 *
 * @param _param the value passed to Y0.
 * @return the result of calculating Y0.
 */
double CNetBthCalculator::Y0(double _param) const {
	return CGSL::Y0(_param);
}

/**
 * The Y1 Bessel function.
 *
 * @param _param the value passed to Y1.
 * @return the result of calculating Y1.
 */
double CNetBthCalculator::Y1(double _param) const {
	return CGSL::Y1(_param);
}

/**
 * The F0 Composite Bessel function.
 *
 * @param _param the value passed to the individual bessel functions.
 * @return the result of calculating F0.
 */
double CNetBthCalculator::F0(double _param) const {
	return CGSL::J0(_param) + (c2/b2)*CGSL::Y0(_param);
}

/**
 * The F1 Composite Bessel function.
 *
 * @param _param the value passed to the individual bessel functions.
 * @return the result of calculating F1.
 */
double CNetBthCalculator::F1(double _param) const {
	return CGSL::J1(_param) + (c2/b2)*CGSL::Y1(_param);
}

