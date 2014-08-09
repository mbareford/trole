/**
 * @file GSL.h
 *
 * This file declares the CGSL class.
 */

#ifndef GSL_H_
#define GSL_H_

/**
 * The CGSL class cannot be instantiated - all of its methods are declared static.
 * Access to the GSL C routines is achieved via this class.
 */
class CGSL {

public:
	static double J0(double _param);
	static double J1(double _param);
	static double Y0(double _param);
	static double Y1(double _param);

	static double integrate(double (*_f)(double, void*), double _min, double _max, void* _param);
	static double find_root(double (*_f)(double, void*), double _ARG_LIM, double _ARG_STEP_SIZE, const int _ROOT_SIGN_PREF, void* _param);

	static int round(double _num);

    static void random_alloc(void);
	static void random_free(void);
	static double random(double _min, double _max);
	static bool random_binary(void);
	static double random_normal(double _mean, double _var);

	static const double PI;

private:
	/// pointer to random number generator
	static void* m_rng;

	/// default construction not allowed
	CGSL(void) {};
	/// default destruction not allowed
	~CGSL(void) {};

	static double root(double (*_f)(double, void*), double _min, double _max, void* _param);

}; // end of CGSL class declaration

#endif /* GSL_H_ */
