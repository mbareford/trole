/**
 * @file RelaxationCalculator.h
 *
 * This file declares the CRelaxationCalculator class.
 */

#ifndef RELAXATIONCALCULATOR_H_
#define RELAXATIONCALCULATOR_H_

/**
 * The CRelaxationCalculator class is used to determine the relaxation state that conserves k/(psi*psi)
 * In addition, this relaxed state may feature a current neutralisation layer and/or a potential envelope.
 * A further possibility is that the relaxed state may feature helical modes.
 */
class CRelaxationCalculator {

public:
	CRelaxationCalculator(double _psi, double _psi_r, double _consrv, double _consrv_r, CRelaxedLoop& _rx);
	/// empty default destructor
	~CRelaxationCalculator(void) {}

	double consrv_err(double _a);

private:
	/// the axial flux that is conserved
	double psi;
	/// the radial range of the axial flux that is conserved
	double psi_r;
	/// the value k/(psi*psi) that is conserved by the relaxed loop
	double consrv;
	/// the radial range (0-_consrv_r) over which k/(psi*psi) is conserved
	double consrv_r;

	/// the relaxed loop
	CRelaxedLoop rx;

}; // end of CRelaxationCalculator class declaration

#endif /* RELAXATIONCALCULATOR_H_ */



