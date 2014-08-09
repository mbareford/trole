/**
 * @file Exception.cpp
 *
 * This file defines the methods of the CException class.
 */

#include "StdInclude.h"

#include "Exception.h"

/**
 * The constructor ensures that the m_msg_len and m_msg attributes
 * are properly initialised.
 *
 * @param _msg the text describing the exception
 */
CException::CException(const char* _msg/*=NULL*/) {
	m_msg_len = 0;
	memset(m_msg, 0, sizeof(m_msg));

	if (_msg) {
		m_msg_len = strlen(_msg);

		if (m_msg_len > sizeof(m_msg)-1) {
			m_msg_len = sizeof(m_msg)-1;
		}

		strncpy(m_msg, _msg, m_msg_len);
	}
}

/**
 * The default copy constructor.
 *
 * @param _x the CException object to be copied
 */
CException::CException(const CException& _x) {
	m_msg_len = 0;
	memset(m_msg, 0, sizeof(m_msg));

	const char* _msg = _x.msg();

	if (_msg) {
		m_msg_len = strlen(_msg);

		if (m_msg_len > sizeof(m_msg)-1) {
			m_msg_len = sizeof(m_msg)-1;
		}

		strncpy(m_msg, _msg, m_msg_len);
	}
}

/**
 * Stream the description of the exception one character at a time
 * to standard output.
 *
 * This function is a friend of the CException class.
 *
 * @param _os the output stream
 * @param _x the exception object
 * @return reference to standard output
 */
ostream& operator<<(ostream& _os, const CException& _x) {
	for (unsigned short i = 0; i < _x.m_msg_len; ++i) {
		_os << _x.m_msg[i];
	}
	_os << NL;

	return _os;
}


/**
 * The constructor ensures that the singularity error code attributes
 * are properly initialised.
 *
 * @param _r the radius
 * @param _msg the text describing the exception
 */
CSingularityException::CSingularityException(double _r, const char* _msg/*=NULL*/) : CException(_msg) {
	m_r = _r;
}

/**
 * The default copy constructor.
 *
 * @param _x the CSingularityException object to be copied
 */
CSingularityException::CSingularityException(const CSingularityException& _sngx) : CException((const CException&) _sngx) {
	m_r = _sngx.r();
}

/**
 * Stream the description of the exception one character at a time
 * to standard output.
 *
 * This function is a friend of the CException class.
 *
 * @param _os the output stream
 * @param _x the exception object
 * @return reference to standard output
 */
ostream& operator<<(ostream& _os, const CSingularityException& _x) {
	_os << "Singularity Error: r=" << _x.m_r << ", Msg=" << (CException) _x;

	return _os;
}


/**
 * The constructor ensures that the NAG error code attributes
 * are properly initialised.
 *
 * @param _code the NAG error code
 * @param _extcode additional NAG error code
 * @param _msg the text describing the exception
 * @param _min the lower bound associated with the NAG function that returned an error
 * @param _max the upper bound associated with the NAG function that returned an error
 */
CNAGException::CNAGException(int _code/*=0*/, int _extcode/*=0*/, const char* _msg/*=NULL*/,
		                     double _min/*=0.0*/, double _max/*=0.0*/) : CException(_msg) {
	m_code = _code;
	m_extcode = _extcode;
	m_min = _min;
	m_max = _max;
}

/**
 * The default copy constructor.
 *
 * @param _x the CNAGException object to be copied
 */
CNAGException::CNAGException(const CNAGException& _nagx) : CException((const CException&) _nagx) {
	m_code = _nagx.code();
	m_extcode = _nagx.extcode();
	m_min = _nagx.min();
	m_max = _nagx.max();
}

/**
 * Stream the description of the exception one character at a time
 * to standard output.
 *
 * This function is a friend of the CException class.
 *
 * @param _os the output stream
 * @param _x the exception object
 * @return reference to standard output
 */
ostream& operator<<(ostream& _os, const CNAGException& _x) {
	if (0.0 == _x.m_min && 0.0 == _x.m_max) {
		_os << "NAG Error: Code=" << _x.m_code << ", ExtCode=" << _x.m_extcode << ", Msg=" << (CException) _x;
	}
	else {
		_os << "NAG Error: Code=" << _x.m_code << ", ExtCode=" << _x.m_extcode << ", Msg=" << (CException) _x
		                          << " (min=" << _x.m_min << ",max=" << _x.m_max << ")." << NL;
	}

	return _os;
}


/**
 * The constructor ensures that the GSL error code attributes
 * are properly initialised.
 *
 * @param _code the GSL error code
 * @param _msg the text describing the exception
 * @param _min the lower bound associated with the GSL function that returned an error
 * @param _max the upper bound associated with the GSL function that returned an error
 */
CGSLException::CGSLException(int _code/*=0*/, const char* _msg/*=NULL*/,
		                     double _min/*=0.0*/, double _max/*=0.0*/) : CException(_msg) {
	m_code = _code;
	m_min = _min;
	m_max = _max;
}

/**
 * The default copy constructor.
 *
 * @param _x the CGSLException object to be copied
 */
CGSLException::CGSLException(const CGSLException& _gslx) : CException((const CException&) _gslx) {
	m_code = _gslx.code();
	m_min = _gslx.min();
	m_max = _gslx.max();
}

/**
 * Stream the description of the exception one character at a time
 * to standard output.
 *
 * This function is a friend of the CException class.
 *
 * @param _os the output stream
 * @param _x the exception object
 * @return reference to standard output
 */
ostream& operator<<(ostream& _os, const CGSLException& _x) {
	if (0.0 == _x.m_min && 0.0 == _x.m_max) {
		_os << "NAG Error: Code=" << _x.m_code << ", Msg=" << (CException) _x;
	}
	else {
		_os << "NAG Error: Code=" << _x.m_code << ", Msg=" << (CException) _x
		                          << " (min=" << _x.m_min << ",max=" << _x.m_max << ")." << NL;
	}

	return _os;
}


/**
 * The constructor ensures that the attributes describing the loop exception
 * are properly initialised.
 *
 * @param _a1 the inner core alpha value
 * @param _a2 the inner core alpha value
 * @param _msg the text describing the exception
 */
CLoopException::CLoopException(double _a1, double _a2, const char* _msg) : CException(_msg) {
	m_a1 = _a1;
	m_a2 = _a2;
}

/**
 * The constructor ensures that the attributes describing the loop exception
 * are properly initialised.
 *
 * @param _a1 the inner core alpha value
 * @param _a2 the inner core alpha value
 * @param _msg the text describing the exception
 * @param _nagx a pointer to an associated NAG exception
 */
CLoopException::CLoopException(double _a1, double _a2, const char* _msg, const CNAGException& _nagx) : CException(_msg), m_nagx(_nagx) {
	m_a1 = _a1;
	m_a2 = _a2;
}

/**
 * The constructor ensures that the attributes describing the loop exception
 * are properly initialised.
 *
 * @param _a1 the inner core alpha value
 * @param _a2 the inner core alpha value
 * @param _msg the text describing the exception
 * @param _gslx a pointer to an associated GSL exception
 */
CLoopException::CLoopException(double _a1, double _a2, const char* _msg, const CGSLException& _gslx) : CException(_msg), m_gslx(_gslx) {
	m_a1 = _a1;
	m_a2 = _a2;
}

/**
 * Stream the description of the exception one character at a time
 * to standard output.
 *
 * This function is a friend of the CException class.
 *
 * @param _os the output stream
 * @param _x the exception object
 * @return reference to standard output
 */
ostream& operator<<(ostream& _os, const CLoopException& _x) {
	if (0 == _x.m_nagx.code() && 0 == _x.m_gslx.code()) {
		_os << "Loop Error: a1=" << _x.m_a1 << ", a2=" << _x.m_a2 << ", Msg=" << (CException) _x;
	}
	else {
		if (0 == _x.m_nagx.code()) {
			_os << "Loop Error: a1=" << _x.m_a1 << ", a2=" << _x.m_a2 << ", Msg=" << (CException) _x << _x.m_gslx;
		}		
		else {
			_os << "Loop Error: a1=" << _x.m_a1 << ", a2=" << _x.m_a2 << ", Msg=" << (CException) _x << _x.m_nagx;
		}
	}

	return _os;
}


/**
 * The constructor ensures that the attributes describing the loop exception
 * are properly initialised.
 *
 * @param _lam the lambda parameter
 * @param _eta the eta parameter
 * @param _msg the text describing the exception
 */
CSmoothLoopException::CSmoothLoopException(double _lam, double _eta, const char* _msg) : CException(_msg) {
	m_lam = _lam;
	m_eta = _eta;
}

/**
 * The constructor ensures that the attributes describing the loop exception
 * are properly initialised.
 *
 * @param _lam the lambda parameter
 * @param _eta the eta parameter
 * @param _msg the text describing the exception
 * @param _nagx a pointer to an associated NAG exception
 */
CSmoothLoopException::CSmoothLoopException(double _lam, double _eta, const char* _msg, const CNAGException& _nagx) : CException(_msg), m_nagx(_nagx) {
	m_lam = _lam;
	m_eta = _eta;
}

/**
 * The constructor ensures that the attributes describing the loop exception
 * are properly initialised.
 *
 * @param _lam the lambda parameter
 * @param _eta the eta parameter
 * @param _msg the text describing the exception
 * @param _gslx a pointer to an associated GSL exception
 */
CSmoothLoopException::CSmoothLoopException(double _lam, double _eta, const char* _msg, const CGSLException& _gslx) : CException(_msg), m_gslx(_gslx) {
	m_lam = _lam;
	m_eta = _eta;
}

/**
 * Stream the description of the exception one character at a time
 * to standard output.
 *
 * This function is a friend of the CException class.
 *
 * @param _os the output stream
 * @param _x the exception object
 * @return reference to standard output
 */
ostream& operator<<(ostream& _os, const CSmoothLoopException& _x) {
	if (0 == _x.m_nagx.code() && 0 == _x.m_gslx.code()) {
		_os << "Loop Error: lam=" << _x.m_lam << ", eta=" << _x.m_eta << ", Msg=" << (CException) _x;
	}
	else {
		if (0 == _x.m_nagx.code()) {
			_os << "Loop Error: lam=" << _x.m_lam << ", eta=" << _x.m_eta << ", Msg=" << (CException) _x << _x.m_gslx;
		}		
		else {
			_os << "Loop Error: lam=" << _x.m_lam << ", eta=" << _x.m_eta << ", Msg=" << (CException) _x << _x.m_nagx;
		}
	}

	return _os;
}
