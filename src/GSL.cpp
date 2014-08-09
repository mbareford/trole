/**
 * @file GSL.cpp
 *
 * This file defines the methods of the CGSL class.
 */

#include "StdInclude.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <sys/time.h>

#include "Exception.h"
#include "GSL.h"

void* CGSL::m_rng = NULL;
const double CGSL::PI = M_PI;


/**
 * the JO Bessel function
 *
 * @param _param the value passed to JO
 * @return the result of calculating J0
 */
double CGSL::J0(double _param) {
	return gsl_sf_bessel_J0(_param);
}

/**
 * the J1 Bessel function
 *
 * @param _param the value passed to J1
 * @return the result of calculating J1
 */
double CGSL::J1(double _param) {
	return gsl_sf_bessel_J1(_param);
}

/**
 * the YO Bessel function
 *
 * @param _param the value passed to YO
 * @return the result of calculating Y0
 */
double CGSL::Y0(double _param) {
	return gsl_sf_bessel_Y0(_param);
}

/**
 * the Y1 Bessel function
 *
 * @param _param the value passed to Y1
 * @return the result of calculating Y1
 */
double CGSL::Y1(double _param) {
	return gsl_sf_bessel_Y1(_param);
}

/**
 * integrate a function _f from _min to _max
 *
 * @param _f the function to be integrated
 * @param _min the lower bound of the integration
 * @param _max the upper bound of the integration
 * @param _param points to a loop object
 * @return the result of calculating the integral
 * @throw CGSLException should underlying GSL function return an error
 */
double CGSL::integrate(double (*_f)(double, void*), double _min, double _max, void* _param) {
	double epsabs = 0.001;     // absolute accuracy
	double epsrel = 0.01;      // relative accuracy
	int max_num_subint = 100;  // maximum number of intervals in integral
	double abserr = 0.0;       // absolute error on integral result
	double res = 0.0;          // integral result
    int err = 0;               // error code
		
    gsl_function gsl_f;
    gsl_f.function = _f;
    gsl_f.params = _param;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(max_num_subint);
	if (NULL == w) {
		throw CGSLException(GSL_ENOMEM, gsl_strerror(GSL_ENOMEM), _min, _max);
	}
    err = gsl_integration_qags(&gsl_f, _min, _max, epsabs, epsrel, max_num_subint, w, &res, &abserr);
    gsl_integration_workspace_free(w);

    if (GSL_SUCCESS != err) {
		throw CGSLException(err, gsl_strerror(err), _min, _max);
	}

	return res;
}

/**
 * find the root of the specified function between min and max
 * parameter values
 *
 * @param _f the function for which to find a root
 * @param _min the lower bound of the search
 * @param _max the upper bound of the search
 * @param _param points to a loop object
 * @return the parameter value at the root
 * @throw CGSLException should underlying GSL function return an error
 */
double CGSL::root(double (*_f)(double, void*), double _min, double _max, void* _param) {
	const double VAL_TOL = 1e-05;
	
	double prev_val = 0.0;	
	double val = 0.0;
    int err = 0;

    gsl_function gsl_f; 
	gsl_f.function = _f;
  	gsl_f.params = _param;

    gsl_root_fsolver* s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	if (NULL == s) {
		throw CGSLException(GSL_ENOMEM, gsl_strerror(GSL_ENOMEM), _min, _max);
	}
    gsl_root_fsolver_set(s, &gsl_f, _min, _max);

	int iter = 0, max_iter = 100;
	do {
		iter++;

		err = gsl_root_fsolver_iterate(s);
		if (GSL_SUCCESS != err && GSL_CONTINUE != err) {
			gsl_root_fsolver_free(s);
			throw CGSLException(err, gsl_strerror(err), _min, _max);
		}

		val = gsl_root_fsolver_root(s);
                
        if (1 == iter) {
			prev_val = val;
        }
	} while (fabs(val-prev_val) > VAL_TOL && iter < max_iter);

	gsl_root_fsolver_free(s);

	return val;
}


/**
 * Search for a function's lowest root in the range -_arg_lim and +_arg_lim. The lowest root is the one
 * associated with the smallest absolute function argument.
 *
 * @param _f the function for which to find a root
 * @param _ARG_LIM the maximum absolute value for the function argument
 * @param _ARG_STEP_SIZE the amount by which the function argument is incremented between
 *                       0 and _arg_lim during root search
 * @param _ROOT_SIGN_PREF the preferred sign for the root in the advent of roots of both signs
 *                        being found in the same interval (1 means positive, -1 means negative
 *                        and 0 means no preference)
 * @param _param points to a loop object
 * @return the value of the function argument at the root
 * @throw CGSLException should underlying GSL function return an error or if root simply cannot be located
 * or _ARG_STEP_SIZE parameter is zero
 */
double CGSL::find_root(double (*_f)(double, void*),
		               const double _ARG_LIM, const double _ARG_STEP_SIZE, const int _ROOT_SIGN_PREF,
		               void* _param) {
	if (0.0 == _ARG_STEP_SIZE) {
		throw CGSLException(0.0, "Error! Zero _ARG_STEP_SIZE.", 0.0, _ARG_LIM);
	}

	const unsigned long INTERVAL_CNT = (unsigned long) (_ARG_LIM / _ARG_STEP_SIZE);
	unsigned long i = 0;

	double prev_arg = 0.0, arg = 0.0, prev_neg_arg = 0.0, neg_arg = 0.0;
	double prev_f_arg = 0.0, f_arg = 0.0, prev_f_neg_arg = 0.0, f_neg_arg = 0.0;

	bool root_fnd = false;
	double root = 0.0;

	f_arg = _f(arg, _param);
	f_neg_arg = _f(neg_arg, _param);

	do {
		i++;

		// look for positive root
		prev_arg = arg;
		prev_f_arg = f_arg;

		arg = i*_ARG_STEP_SIZE;
		f_arg = _f(arg, _param);

		if (0.0 == prev_f_arg || 0.0 == f_arg
			|| (prev_f_arg < 0.0 && f_arg > 0.0) || (prev_f_arg > 0.0 && f_arg < 0.0)) {
			// positive root found in this interval
			root_fnd = true;

			root = CGSL::root(_f, prev_arg, arg, _param);
		}

		// look for negative root
		prev_neg_arg = neg_arg;
		prev_f_neg_arg = f_neg_arg;

		neg_arg = -1.0*i*_ARG_STEP_SIZE;
		f_neg_arg = _f(neg_arg, _param);

		if (0.0 == prev_f_neg_arg || 0.0 == f_neg_arg
			|| (prev_f_neg_arg < 0.0 && f_neg_arg > 0.0) || (prev_f_neg_arg > 0.0 && f_neg_arg < 0.0)) {
			// negative root found in this interval
			double t = CGSL::root(_f, neg_arg, prev_neg_arg, _param);
			if (root_fnd) {
				// positive root has already been found in this interval
				if ((-1 == _ROOT_SIGN_PREF) || (0 == _ROOT_SIGN_PREF && fabs(t) < fabs(root))) {
					root = t;
				}
			}
			else {
				root_fnd = true;
				root = t;
			}
		}

	} while (!root_fnd && i <= INTERVAL_CNT);

	if (!root_fnd) {
		throw CGSLException(0.0, "Error! Cannot find root.", 0.0, _ARG_LIM);
	}

	return root;
}


/**
 * Round a double to the nearest integer.
 *
 * @param _num the double number
 */
int CGSL::round(double _num) {
	int i_num = (long) _num;
	double rem = _num - ((double) i_num);

	if (rem >= 0.5) {
		i_num++;
	}
	else if (rem <= 0.5) {
		i_num--;
	}

	return i_num;
}


/**
 * Allocate the GSL random number generator.
 *
 * @throw CGSLException random number generator already exists.
 * @throw CGSLException not enough memory to allocate random number generator.
 */
void CGSL::random_alloc(void) {
	if (NULL != CGSL::m_rng) {
		throw CGSLException(GSL_FAILURE, "Random number generator already exists.");
	}

    gsl_rng_env_setup();

	CGSL::m_rng = gsl_rng_alloc(gsl_rng_taus);  
  	if (NULL == CGSL::m_rng) {
		throw CGSLException(GSL_ENOMEM, "Not enough memory to allocate random number generator.");
	}

	time_t tm;
	time(&tm);
	gsl_rng_set((const gsl_rng*) CGSL::m_rng, tm);
}

/**
 * Free the GSL random number generator.
 *
 * @throw CGSLException no random number generator
 */
void CGSL::random_free(void) {
	if (NULL == CGSL::m_rng) {
		throw CGSLException(GSL_ENOMEM, "No random number generator.");
	}
	gsl_rng_free((gsl_rng*) CGSL::m_rng);
}

/**
 * Return a random integer that falls within the specified range.
 *
 * @param _min the lower boundary of the random number range
 * @param _max the upper boundary of the random number range
 *
 * @throw CGSLException no random number generator
 */
double CGSL::random(double _min, double _max) {
	if (NULL == CGSL::m_rng) {
		throw CGSLException(GSL_ENOMEM, "No random number generator.");
	}
	double ran = gsl_rng_uniform((const gsl_rng*) CGSL::m_rng);
    return _min + (_max-_min)*ran;
}

/**
 * Return a random binary bit.
 */
bool CGSL::random_binary(void) {
	return random(0.0, 1.0) >= 0.5;
}

/**
 * Return an number randomly selected according to a normal distribution.
 *
 * @param _mean the mean of the normal distribution.
 * @param _var the variance of the normal distribution.
 *
 * @throw CGSLException no random number generator
 */
double CGSL::random_normal(double _mean, double _var) {
 	if (NULL == CGSL::m_rng) {
		throw CGSLException(GSL_ENOMEM, "No random number generator.");
	}
	
	double ran = gsl_ran_gaussian((const gsl_rng*) CGSL::m_rng, _var);
	return _mean + ran;
}
