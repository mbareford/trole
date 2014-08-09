/**
 * @file RelaxationCalculator.cpp
 *
 * This file defines the methods of the CRelaxationCalculator class.
 */

#include "StdInclude.h"

#include "Exception.h"
#include "GSL.h"
#include "Profile.h"
#include "RelaxedLoop.h"
#include "RelaxationCalculator.h"


/**
 * Construct a CRelaxationCalculator object for a loop that has a single outer layer.
 *
 * @param _psi the axial flux that is conserved
 * @param _psi_r the radial range of the axial flux that is conserved
 * @param _consrv the value k/(psi*psi) that is conserved by the relaxed loop
 * @param _consrv_r the radial range (0-_consrv_r) over which k/(psi*psi) is conserved
 * @param _rx the relaxed state - its dimensions must be initialised
 */
CRelaxationCalculator::CRelaxationCalculator(double _psi,
		                                     double _psi_r,
		                                     double _consrv,
		                                     double _consrv_r,
		                                     CRelaxedLoop& _rx)
                                             : psi(_psi),
                                               psi_r(_psi_r),
                                               consrv(_consrv),
                                               consrv_r(_consrv_r),
                                               rx(_rx) {}

/**
 * When a loop undergoes a relaxation the quantity k/(psi^2) is conserved.
 * This function is used to find the value of the relaxation alpha that
 * conserves this quantity. It does this by returning the conservation error:
 * the difference in k/(psi^2) between the threshold and relaxed states.
 * If this quantity is conserved the difference should be zero.
 *
 * Essentially, the rx attribute is tested with a range of a1 values in order to
 * find the a1 value that conserves k/(psi*psi). Whether or not the relaxed
 * state features current neutralisation is already expressed within the rx object.
 *
 * An exception is thrown (from CRelaxedLoop::init_field_coeffs) if psi_r is zero,
 * since this would be also make t_psi zero.
 *
 * @param _a1 the prospective relaxed alpha value.
 * @return zero if _a1 really does represent the relaxed state.
 */
double CRelaxationCalculator::consrv_err(double _a1) {
	rx.init_field(psi, psi_r, _a1);

	double t_psi = rx.calc_psi(consrv_r);

	double t_k = rx.calc_k(consrv_r);

	return consrv - (t_k/(t_psi*t_psi));
}
