/**
 * @file Exception.h
 *
 * This file declares the CException class.
 */

#ifndef EXCEPTION_
#define EXCEPTION_

using namespace std;

/// the maximum number of characters that can be used to describe an exception
const unsigned short MAX_MSG_LEN = 1024;

/**
 * The CExcepton class records information (e.g., textual descriptions)
 * of exceptional events.
 */
class CException {

public:
	CException(const char* _msg = NULL);
	CException(const CException& _x);
	~CException(void) {};

	const char* msg(void) const { return m_msg; }

private:
	/// the number of characters in the description of the exception
	unsigned short m_msg_len;

	/// the description of the exception
	char m_msg[MAX_MSG_LEN+1];

	friend ostream& operator<<(ostream& _os, const CException& _x);

}; // end of CException class declaration

/**
 * The CNAGExcepton class records information about errors returned by
 * NAG routines.
 */
class CNAGException : public CException {

public:
	CNAGException(int _code = 0, int _extcode = 0, const char* _msg = NULL, double _min = 0.0, double _max = 0.0);
	CNAGException(const CNAGException& _nagx);
	~CNAGException(void) {};

	int code(void) const { return m_code; }
	int extcode(void) const { return m_extcode; }
	double min(void) const { return m_min; }
	double max(void) const { return m_max; }

private:
	/// NAG error code.
	int m_code;
	/// Additional NAG error code.
	int m_extcode;
	/// lower bound of NAG function if any
	double m_min;
	/// upper bound of NAG function if any
	double m_max;

friend ostream& operator<<(ostream& _os, const CNAGException& _x);

}; // end of CNAGException class declaration


/**
 * The CGSLExcepton class records information about errors returned by
 * GSL routines.
 */
class CGSLException : public CException {

public:
	CGSLException(int _code = 0, const char* _msg = NULL, double _min = 0.0, double _max = 0.0);
	CGSLException(const CGSLException& _gslx);
	~CGSLException(void) {};

	int code(void) const { return m_code; }
	double min(void) const { return m_min; }
	double max(void) const { return m_max; }

private:
	/// GSL error code.
	int m_code;
	/// lower bound of GSL function if any
	double m_min;
	/// upper bound of GSL function if any
	double m_max;

friend ostream& operator<<(ostream& _os, const CGSLException& _x);

}; // end of CGSLException class declaration


/**
 * The CSingularityException class records information whenever
 * an intergrand function attempts a divide by zero.
 */
class CSingularityException : public CException {

public:
	CSingularityException(double _r, const char* _msg = NULL);
	CSingularityException(const CSingularityException& _nagx);
	~CSingularityException(void) {};

	double r(void) const { return m_r; }

private:
	/// radius.
	double m_r;

friend ostream& operator<<(ostream& _os, const CSingularityException& _x);

}; // end of CSingularityException class declaration


/**
 * The CLoopExcepton class records information about errors returned by
 * methods of the CLoop class.
 */
class CLoopException : public CException {

public:
	CLoopException(double _a1, double _a2, const char* _msg);
	CLoopException(double _a1, double _a2, const char* _msg, const CNAGException& _nagx);
	CLoopException(double _a1, double _a2, const char* _msg, const CGSLException& _gslx);
	~CLoopException(void) {};

	double a1(void) const { return m_a1; }
	double a2(void) const { return m_a2; }

private:
	/// inner core alpha value
	double m_a1;
	/// outer core alpha value
	double m_a2;
	/// the NAG error associated with the loop exception
	CNAGException m_nagx;
	/// the GSL error associated with the loop exception
	CGSLException m_gslx;

friend ostream& operator<<(ostream& _os, const CLoopException& _x);

}; // end of CLoopException class declaration

/**
 * The CSmoothLoopExcepton class records information about errors returned by
 * methods of the CSmoothLoop class.
 */
class CSmoothLoopException : public CException {

public:
	CSmoothLoopException(double _lam, double _eta, const char* _msg);
	CSmoothLoopException(double _lam, double _eta, const char* _msg, const CNAGException& _nagx);
	CSmoothLoopException(double _lam, double _eta, const char* _msg, const CGSLException& _gslx);
	~CSmoothLoopException(void) {};

	double lam(void) const { return m_lam; }
	double eta(void) const { return m_eta; }

private:
	/// lambda parameter
	double m_lam;
	/// eta parameter
	double m_eta;
	/// the NAG error associated with the loop exception
	CNAGException m_nagx;
	/// the GSL error associated with the loop exception
	CGSLException m_gslx;

friend ostream& operator<<(ostream& _os, const CSmoothLoopException& _x);

}; // end of CSmoothLoopException class declaration

#endif /* EXCEPTION_ */
