/**
 * @file NetBthCalculator.h
 *
 * This file declares the CNetBthCalculator class.
 */

#ifndef CNETBTHCALCULATOR_H_
#define CNETBTHCALCULATOR_H_

/**
 * The CNetBthCalculator class is used to achieve current neutralisation for
 * loop threshold and loop relaxed states. It is used by GSL root finding routines
 * to find the alpha value that causes bth to drop to zero at a certain radius.
 */
class CNetBthCalculator {

public:
	CNetBthCalculator(double _r1, double _r2, double _a1);
	CNetBthCalculator(double _r1, double _r2, double _r3, double _a1, double _a2, double _b2, double _c2);

	/// empty default destructor
	~CNetBthCalculator(void) {}

	double netbth(double _a) const;

private:
	/**
	 * If true the net bth is calculated at the outer boundary (r3) of a loop volume that features
	 * a core (0-r1) enclosed by two outer layers (r1-r2 and r2-r3). Otherwise the net bth is calculated
	 * at the outer boundary (r2) of a simpler loop volume; one that has 2 concentric regions,
	 * a core (0-r1) and just one outer layer (r1-r2).
	 */
	bool double_outer_layer;

	/// the core radius
	double r1;
	/// the first outer layer radius
	double r2;
	/// the second outer layer radius (only used if double_outer_layer is true)
	double r3;

	/// the alpha value within the core (0-r1)
	double a1;
	/// the alpha value within the first outer layer (r1-r2) (only used if double_outer_layer is true)
	double a2;

	/// the magnetic field coefficients in the core and in the first outer layer (only used if double_outer_layer is true)
	double b1, b2, c2;


	double J0(double _param) const;
	double J1(double _param) const;
	double Y0(double _param) const;
	double Y1(double _param) const;
	double F0(double _param) const;
	double F1(double _param) const;

}; // end of CNetBthCalculator class declaration

#endif /* CNETBTHCALCULATOR_H_ */



