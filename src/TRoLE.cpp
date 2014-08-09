/**
 * @file TROLE.cpp
 *
 * Taylor Relaxation Of Loop Ensemble
 *
 * The loops are represented as force-free fields with
 * piece-wise constant (or smooth) alpha profiles before onset of
 * linear instability and with linear alpha profiles after
 * helicity-conserving Taylor relaxation.
 *
 * v2.0
 */

#include "StdInclude.h"

#include "Exception.h"
#include "GSL.h"
#include "Profile.h"
#include "Threshold.h"
#include "RelaxedLoop.h"
#include "SmoothLoop.h"
#include "Loop.h"
#include "LoopWalker.h"
#include "ParameterSpace.h"


/// the indexes and codes for the argv arguments
const unsigned short ARGSFN_INDEX         = 1;
const unsigned short MODE_INDEX           = 1;
	const char EXAMINE_MODE         = 'E';
	const char EXAMINE_SMOOTH_MODE  = 'F';
	const char SHEPHERD_MODE        = 'S';

const unsigned short R1_INDEX             = 2;
const unsigned short R2_INDEX             = 3;
const unsigned short R3_INDEX             = 4;
const unsigned short R4_INDEX             = 5;
const unsigned short L_INDEX              = 6;
const unsigned short RELAXABLE_FLAG_INDEX = 7;
const unsigned short APPROX_A3_FLAG_INDEX = 8;

const unsigned short NM_CODE_INDEX        = 9;
	const char NM_CODE_PSI  = 'A';
	const char NM_CODE_B1   = 'B';

const unsigned short PF_RS_INDEX          = 10;

// argument indexes for examine mode
///////////////////////////////////////////////
const unsigned short PF_RXS_FLAG_INDEX    = 11;
const unsigned short AFN_INDEX            = 12;
const unsigned char OFNP_INDEX            = 13;
const unsigned short EFN_INDEX            = 14;
///////////////////////////////////////////////

// argument indexes for shepherd mode
///////////////////////////////////////////////
const unsigned short ITA_FN_INDEX         = 11;
const unsigned short ITP_FN_INDEX         = 12;
const unsigned short MIN_ITWR_INDEX       = 13;
const unsigned short MAX_ITWR_INDEX       = 14;
const unsigned short WR_BIN_CNT_INDEX     = 15;
const unsigned short WRBC_FN_INDEX        = 16;
const unsigned short WRPL_FN_INDEX        = 17;
const unsigned short RX_RAD_INDEX         = 18;
const unsigned short RX_CNT_INDEX         = 19;
const unsigned short RX_FN_INDEX          = 20;
const unsigned short STEP_SIZE_INDEX      = 21;
const unsigned short CORRELATED_INDEX     = 22;
const unsigned short SIGMA_INIT_INDEX     = 23;
const unsigned short SIGMA_INDEX          = 24;
const unsigned short LIFETIME_INDEX       = 25;
const unsigned short OBV_DAT_FN_INDEX     = 26;
///////////////////////////////////////////////

unsigned short ITR_ARG_COUNT = 15;
unsigned short SIM_ARG_COUNT = 27;

/**
 * Construct an output filename based on a prefix and a number.
 *
 * @param fn_prefix the path and name prefix.
 * @param fn_num the number to be appended to prefix.
 * @return constructed filename
 * @throw CException should any of the following occur:
 *        1. NULL file name prefix.
 */
string create_ofn(char* fn_prefix, unsigned int fn_num) {
	stringstream fn(stringstream::in | stringstream::out);

	if (!fn_prefix) {
        throw CException("Error! NULL file name prefix.");
	}

	fn << fn_prefix << fn_num << ".txt";

	return fn.str();
}

/**
 * Construct an output filename based on a prefix, a number and a tag.
 *
 * @param fn_prefix the path and name prefix.
 * @param fn_num the number to be appended to prefix.
 * @param fn_tag an additional tag appended to filename.
 * @return constructed filename
 * @throw CException should any of the following occur:
 *        1. NULL file name prefix;
 *        2. NULL file name tag.
 */
string create_ofn(char* fn_prefix, unsigned int fn_num, char* fn_tag) {
	stringstream fn(stringstream::in | stringstream::out);

	if (!fn_prefix) {
		throw CException("Error! NULL file name prefix.");
	}

	if (!fn_tag) {
		throw CException("Error! NULL file name tag.");
	}

	fn << fn_prefix << fn_num << "_" << fn_tag << ".txt";

	return fn.str();
}

/**
 * Calculate the properties of a series of loops whose alpha values are
 * specified in the indicated file.
 * The loops are identical in length and radial dimensions.
 * Radial profiles of the loops' axial, azimuthal and twist fields are calculated.
 * The axial flux, helicity and magnetic energy are calculated.
 * If a loop is indicated as being unstable the details of the loop's relaxed
 * state (alpha and energy) are also calculated.
 *
 * @throw CException should any of the following occur:
 *        1. invalid number of arguments;
 *        2. invalid normalisation code;
 *        3. cannot open alpha file;
 *        4. cannot open output file;
 *        5. cannot open error file;
 *        6. cannot read alpha file;
 *        7. cannot open output field profile file;
 */
void examine(int argc, char argv[][512]) {
	if (ITR_ARG_COUNT != argc) {
		throw CException("Error! Invalid number of arguments.");
	}

	ifstream af;
	ofstream of, ef;

	try {
		double r1 = atof(argv[R1_INDEX]);
		double r2 = atof(argv[R2_INDEX]);
		double r3 = atof(argv[R3_INDEX]);
		double r4 = atof(argv[R4_INDEX]);

		double l = atof(argv[L_INDEX]);
		bool relaxable = (0 != atoi(argv[RELAXABLE_FLAG_INDEX]));
		bool approx_a3 = (0 != atoi(argv[APPROX_A3_FLAG_INDEX]));

		double a1 = 0.0, a2 = 0.0;

		char nm_code = argv[NM_CODE_INDEX][0];
		if (NM_CODE_PSI != nm_code && NM_CODE_B1 != nm_code) {
			throw CException("Error! Invalid normalisation code.");
		}

		unsigned long rs = atol(argv[PF_RS_INDEX]);
		bool pf_rxs = (0 != atoi(argv[PF_RXS_FLAG_INDEX]));

		CLoop loop;

        af.open(argv[AFN_INDEX]);
		if (!af) {
			throw CException("Error! Cannot open alpha file.");
		}

		of.open(create_ofn(argv[OFNP_INDEX], 0).data());
		if (!of) {
			throw CException("Error! Cannot open output file.");
		}

		ef.open(argv[EFN_INDEX]);
		if (!ef) {
			throw CException("Error! Cannot open error file.");
		}

		// count variable for the number of alpha pairs in the alpha file
		unsigned int count = 0;

		// construct a loop for every alpha pair
		while (!af.eof()) {
			af >> a1 >> a2;
			if (0.0 == a1) {
				int a = 10;
				a++;
			}
			if (af.eof()) {
				break;
			}
			else {
				if (!af.good()) {
					throw CException("Error! Cannot read alpha file.");
				}
			}

			try {
				loop.init_dimensions(r1, r2, r3, r4, l);
				loop.init_field(a1, a2, (NM_CODE_B1 == nm_code), approx_a3, relaxable);
				loop.profile(rs, pf_rxs);
			}
			catch (CSingularityException& x) {
				// a singularity was encountered during one of the twist or shear integral calculations
				// output error, reset loop and carry on
				ef << x << " a1 = " << a1 << " a2 = " << a2 << NL;
				//ef << x.r() << ' ' << a1 << ' ' << a2 << NL;
				loop.reset();
			}
			catch (CGSLException& x) {
				ef << x << " a1 = " << a1 << " a2 = " << a2 << NL;
				//ef << x.code() << ' ' << a1 << ' ' << a2 << NL;
				loop.reset();
			}
			catch (CLoopException& x) {
				ef << x;
				continue;
			}
			catch (CException& x) {
				ef << x;
				continue;
			}

			count++;

			// output loop properties
			of << loop;

			// output loop profile
			///////////////////////////////////////////////////////////////////////
			ofstream of_field_pf(create_ofn(argv[OFNP_INDEX], count).data());
			if (!of_field_pf) {
				throw CException("Error! Cannot open output field profile file.");
			}

			ofstream of_rx_pf;
			if (pf_rxs) {
				of_rx_pf.open(create_ofn(argv[OFNP_INDEX], count, (char*) "rx").data());
				if (!of_rx_pf) {
					throw CException("Error! Cannot open output relaxation profile file.");
				}
			}

			try {
				loop.out_field_pfs(of_field_pf);
				if (pf_rxs) {
					loop.out_rx_pfs(of_rx_pf);
				}
			}
			catch (CException& x) {
				of_field_pf.clear();
				of_field_pf.close();
				if (pf_rxs) {
					of_rx_pf.clear();
					of_rx_pf.close();
				}
				throw x;
			}

			if (of_field_pf.is_open()) {
				of_field_pf.clear();
				of_field_pf.close();
			}
			if (pf_rxs && of_rx_pf.is_open()) {
				of_rx_pf.clear();
				of_rx_pf.close();
			}
			///////////////////////////////////////////////////////////////////////

		} // end of <while (!af.eof())> loop

		ef.clear();
		ef.close();

		of.clear();
		of.close();

		af.clear();
		af.close();

	} // end of try block

	catch (CException& x) {
		if (ef.is_open()) {
			ef.clear();
			ef.close();
		}

		if (of.is_open()) {
			of.clear();
			of.close();
		}

		if (af.is_open()) {
			af.clear();
			af.close();
		}

		throw x;
	}
}

/**
 * Calculate the properties of a series of loops whose alpha values are
 * specified in the indicated file.
 * The loops are identical in length and radial dimensions.
 * Radial profiles of the loops' axial, azimuthal and twist fields are calculated.
 * The axial flux, helicity and magnetic energy are calculated.
 * If a loop is indicated as being unstable the details of the loop's relaxed
 * state (alpha and energy) are also calculated.
 *
 * @throw CException should any of the following occur:
 *        1. invalid number of arguments;
 *        2. invalid normalisation code;
 *        3. cannot open alpha file;
 *        4. cannot open output file;
 *        5. cannot open error file;
 *        6. cannot read alpha file;
 *        7. cannot open output field profile file;
 */
void examine_smooth(int argc, char argv[][512]) {
	if (ITR_ARG_COUNT != argc) {
		throw CException("Error! Invalid number of arguments.");
	}

	ifstream lef;
	ofstream of, ef;

	try {
	    double r1 = atof(argv[R1_INDEX]);
        double r2 = atof(argv[R2_INDEX]);
		double l = atof(argv[L_INDEX]);
	    
		bool relaxable = (0 != atoi(argv[RELAXABLE_FLAG_INDEX]));

		double lam = 0.0, eta = 0.0;

		char nm_code = argv[NM_CODE_INDEX][0];
		if (NM_CODE_PSI != nm_code && NM_CODE_B1 != nm_code) {
			throw CException("Error! Invalid normalisation code.");
		}

		unsigned long rs = atol(argv[PF_RS_INDEX]);
		bool pf_rxs = (0 != atoi(argv[PF_RXS_FLAG_INDEX]));

		CSmoothLoop loop;

		lef.open(argv[AFN_INDEX]);
		if (!lef) {
			throw CException("Error! Cannot open lambda-eta file.");
		}

		of.open(create_ofn(argv[OFNP_INDEX], 0).data());
		if (!of) {
			throw CException("Error! Cannot open output file.");
		}

		ef.open(argv[EFN_INDEX]);
		if (!ef) {
			throw CException("Error! Cannot open error file.");
		}

		// count variable for the number of alpha pairs in the alpha file
		unsigned int count = 0;

		// construct a loop for every alpha pair
		while (!lef.eof()) {
			lef >> lam >> eta;
			
			if (lef.eof()) {
				break;
			}
			else {
				if (!lef.good()) {
					throw CException("Error! Cannot read lambda-eta file.");
				}
			}

			try {
				loop.init_dimensions(r1, r2, l);
				loop.init_field(lam, eta, (NM_CODE_B1 == nm_code), relaxable);
				loop.profile(rs, pf_rxs);
			}
			catch (CSingularityException& x) {
				// a singularity was encountered during one of the twist or shear integral calculations
				// output error, reset loop and carry on
				ef << x << " lam = " << lam << " eta = " << eta << NL;				
				loop.reset();
			}
			catch (CGSLException& x) {
				ef << x << " lam = " << lam << " eta = " << eta << NL;
				loop.reset();
			}
			catch (CSmoothLoopException& x) {
				ef << x;
				continue;
			}
			catch (CException& x) {
				ef << x;
				continue;
			}

			count++;

			// output loop properties
			of << loop;

			// output loop profile
			///////////////////////////////////////////////////////////////////////
			ofstream of_field_pf(create_ofn(argv[OFNP_INDEX], count).data());
			if (!of_field_pf) {
				throw CException("Error! Cannot open output field profile file.");
			}

			ofstream of_rx_pf;
			if (pf_rxs) {
				of_rx_pf.open(create_ofn(argv[OFNP_INDEX], count, (char*) "rx").data());
				if (!of_rx_pf) {
					throw CException("Error! Cannot open output relaxation profile file.");
				}
			}

			try {
				loop.out_field_pfs(of_field_pf);
				if (pf_rxs) {
					loop.out_rx_pfs(of_rx_pf);
				}
			}
			catch (CException& x) {
				of_field_pf.clear();
				of_field_pf.close();
				if (pf_rxs) {
					of_rx_pf.clear();
					of_rx_pf.close();
				}
				throw x;
			}

			if (of_field_pf.is_open()) {
				of_field_pf.clear();
				of_field_pf.close();
			}
			if (pf_rxs && of_rx_pf.is_open()) {
				of_rx_pf.clear();
				of_rx_pf.close();
			}
			///////////////////////////////////////////////////////////////////////

		} // end of <while (!af.eof())> loop

		ef.clear();
		ef.close();

		of.clear();
		of.close();

		lef.clear();
		lef.close();

	} // end of try block

	catch (CException& x) {
		if (ef.is_open()) {
			ef.clear();
			ef.close();
		}

		if (of.is_open()) {
			of.clear();
			of.close();
		}

		if (lef.is_open()) {
			lef.clear();
			lef.close();
		}

		throw x;
	}
}

/**
 * Shepherd a loop as it moves through a stability region demarcated
 * by the indicated threshold.
 *
 * The loop relaxes to the a1=a2 line every time it reaches threshold.
 * Helicity is conserved.
 * The details of the relaxation event are recorded.
 *
 * @throw CException should any of the following occur:
 *        1. invalid number of arguments;
 *        2. invalid normalisation code;
 *        3. invalid relaxation radius;
 *        4. zero relaxation count;
 *        5. negative minimum energy release;
 *        6. negative maximum energy release;
 *        7. minimum energy release >= maximum energy release;
 *        8. zero energy release bin count;
 *        9. cannot open relaxation data output file;
 *        10. step size must be greater than zero;
 *        11. sigma parameters must be greater than 0 if walks are correlated;
 *        12. zero lifetime;
 *        13. cannot open energy release output file;
 *        14. cannot open observational data output file;
 *        15. zero tot_step_cnt.
 */
void shepherd(int argc, char argv[][512]) {
	if (SIM_ARG_COUNT != argc) {
		throw CException("Error! Invalid number of arguments.");
	}
	// output files
	ofstream rxf, wr_bcf, wr_plf, odf;

	// assume that the threshold will be defined in alpha space
	bool tw_threshold = false;

	try {
		double r1 = atof(argv[R1_INDEX]);
		double r2 = atof(argv[R2_INDEX]);
		double r3 = atof(argv[R3_INDEX]);
		double r4 = atof(argv[R4_INDEX]);
		double l = atof(argv[L_INDEX]);
		bool approx_a3 = (0 != atoi(argv[APPROX_A3_FLAG_INDEX]));

		char nm_code = argv[NM_CODE_INDEX][0];
		if (NM_CODE_PSI != nm_code && NM_CODE_B1 != nm_code) {
			throw CException("Error! Invalid normalisation code.");
		}

		unsigned long rs = atol(argv[PF_RS_INDEX]);

		// construct and initialise the parameter space with correct instability threshold
		///////////////////////////////////////////////////////////////////////////////////
		CParameterSpace ps;
		if (0 != strcmp(argv[ITA_FN_INDEX], argv[ITP_FN_INDEX])) {
			// using a threshold defined in twist space
			tw_threshold = true;
			ps.set_threshold(argv[ITA_FN_INDEX], argv[ITP_FN_INDEX]);
		}
		else {
			// using a threshold defined in alpha space
			ps.set_threshold(argv[ITA_FN_INDEX]);
		}
		///////////////////////////////////////////////////////////////////////////////////

		// get the relaxation radius; if zero, relaxation is far as q=1 surface
		double rx_rad = atof(argv[RX_RAD_INDEX]);
		if (0.0 != rx_rad && (rx_rad < r1 || rx_rad > r4)) {
			throw CException("Error! Invalid relaxation radius.");
		}

		// the number of times the loop will undergo relaxation
		unsigned long rx_cnt = atol(argv[RX_CNT_INDEX]);
		if (0 == rx_cnt) {
			throw CException("Error! Zero relaxation count.");
		}

		// declare and initialise the variables used during the simulations of relaxations
		unsigned long step_cnt = 0;
		double tot_step_cnt = 0.0;
		double rx_r = 0.0, rx_a = 0.0;
		pt rx_pt = {0.0, 0.0, 0, 0.0}, crs_pt = {0.0, 0.0, 0, 0.0}, alpha_crs_pt = {0.0, 0.0, 0, 0.0};
		double wr = 0.0, tot_wr = 0.0;
		unsigned long rx_i = 0;
		double half_width = 0.0;

		// receive the minimum and maximum energy releases available from the instability threshold
		double min_wr = atof(argv[MIN_ITWR_INDEX]);
		double max_wr = atof(argv[MAX_ITWR_INDEX]);
		if (min_wr < 0.0) {
			throw CException("Error! Negative minimum energy release.");
		}
		if (max_wr < 0.0) {
			throw CException("Error! Negative maximum energy release.");
		}
		if (min_wr >= max_wr) {
			throw CException("Error! Minimum energy release >= maximum energy release.");
		}

		// receive the number of bins for the energy release profile
		unsigned long wr_bin_cnt = atol(argv[WR_BIN_CNT_INDEX]);
		if (0 == wr_bin_cnt) {
			throw CException("Error! Zero energy release bin count.");
		}

		// construct and initialise the energy release profile
		CProfile pf_wr;
		pf_wr.init(min_wr, max_wr, wr_bin_cnt);

		// open the relaxation data file
		rxf.open(argv[RX_FN_INDEX]);
		if (!rxf) {
			throw CException("Error! Cannot open relaxation data output file.");
		}

		// read in the step size
		double step_size = atof(argv[STEP_SIZE_INDEX]);
		if (step_size <= 0.0) {
			throw CException("Error! Step size must be greater than zero.");
		}

		// is the walk correlated
		bool correlated = (0 != atoi(argv[CORRELATED_INDEX]));
		double sigma_init = atof(argv[SIGMA_INIT_INDEX]);
		double sigma = atof(argv[SIGMA_INDEX]);
		if (correlated && (sigma_init <= 0.0 || sigma <= 0.0)) {
			throw CException("Error! Sigma parameters must be greater than 0 if walks are correlated.");
		}

		// read in the loop lifetime
		long lifetime = atol(argv[LIFETIME_INDEX]), age = 0;
		if (0 == lifetime) {
			throw CException("Error! Zero lifetime.");
		}

		// open the observational data file
		odf.open(argv[OBV_DAT_FN_INDEX]);
		if (!odf) {
			throw CException("Error! Cannot open observational data output file.");
		}
		// construct the loop objects
		CLoop loop;
		loop.init_dimensions(r1, r2, r3, r4, l);
		loop.init_field(0.0, 0.0, (NM_CODE_B1 == nm_code), approx_a3, false);
                CLoopWalker loop_walker;
		loop_walker.init(step_size);
		if (correlated) {
			loop_walker.correlate(sigma_init, sigma);
		}
                CRelaxedLoop rx_loop;
                
		// the loop traverses the stability region and relaxes whenever
		// it crosses the threshold
        //printf("rx_cnt=%ld\n",rx_cnt);
		for (rx_i = 0, age = 0; rx_i < rx_cnt; ++rx_i, ++age) {
			printf("rx_i=%ld\n",rx_i);
            //printf("age=%ld\n",age);

			// initialise variables that hold coordinates of threshold crossing
			memset(&crs_pt, 0, sizeof(crs_pt));
			memset(&alpha_crs_pt, 0, sizeof(alpha_crs_pt));

			if (age == lifetime) {
				// simulation continues with new loop
				age = 0;
				// end the walk to ensure that new starting position is chosen
				// next time ps.traverse_stability_region() is called
				loop_walker.stop();
			}

			// traverse stability region until a valid relaxed state has been achieved
			do {				
				// restart loop if previous relaxation was invalid
				if (loop_walker.started()) {
					loop_walker.restart();
				}

				// traverse the stability region until instability threshold is crossed
				ps.traverse_stability_region(&loop_walker, crs_pt, alpha_crs_pt);

				// initialise the threshold state
				loop.init_field(alpha_crs_pt.x, alpha_crs_pt.y, approx_a3, true);
				
				// get relaxed state
				// relax loop according to Scenario 2.1
				if (rx_rad <= 0.0) {
					// relax loop as far as q=1 surface
					rx_loop = loop.relax();
				}
				else {
					// relax loop as far as rx_rad
					rx_loop = loop.relax(rx_rad);
				}
				// get relaxation radius
				rx_r = rx_loop.get_r1();
				// get relaxed alpha
				rx_a = rx_loop.get_a1();
				// exit from loop if relaxed state is stable
			} while (!ps.within_threshold(rx_a, rx_a, true));

			// add the steps taken during walk to instability
			step_cnt = loop_walker.get_step_cnt();
			tot_step_cnt += step_cnt;

			// output the details of the instability threshold crossing
			rxf << crs_pt.x << WS << crs_pt.y << WS << crs_pt.type << WS << crs_pt.ext1 << WS;
			if (tw_threshold) {
				rxf << alpha_crs_pt.x << WS << alpha_crs_pt.y << WS;
			}
			
			// ensure walk will continue from last relaxation position
			loop_walker.start(rx_a, rx_a);

			wr = loop.calc_w(rx_r) - rx_loop.calc_w(rx_r);
			tot_wr += wr;

			// update the energy release profile
			pf_wr.plot_count(wr);

			// output the relaxation radius
			rxf << rx_r << WS;

			// output the relaxation points in alpha space
			rx_pt.x = rx_a;
			rx_pt.y = rx_a;
			rxf << rx_pt.x << WS << rx_pt.y << WS;

			if (tw_threshold) {
				// ps represents some twist space, so we need to translate
				// alpha relaxation point to the equivalent in twist space
				loop.init_field(rx_pt.x, rx_pt.y, approx_a3, false);
				loop.profile(rs);
				//rx_pt.x = loop.get_twist(r1);
				//rx_pt.y = loop.get_twist(r2);
				rx_pt.x = loop.get_twist_av(0.0, r1);
				rx_pt.y = loop.get_twist_av(r1, r2);

				// output the relaxation points in twist space
				rxf << rx_pt.x << WS << rx_pt.y << WS;
			}

			// output other details of the relaxation event
			rxf << step_cnt << WS;
			rxf << tot_step_cnt << WS;
			rxf << wr << WS;
			rxf << NL;

			// output energy release details in the form of mock
			// observational data
			half_width = (wr-min_wr)/(max_wr-min_wr);
			odf << tot_step_cnt << WS << half_width*(-1.0) << WS << 0.0 << NL;
			odf << tot_step_cnt << WS << 0.0 << WS << wr << NL;
			odf << tot_step_cnt << WS << half_width << WS << 0.0 << NL;

		} // end of <for (rx_i = 0; rx_i < rx_cnt; ++rx_i)>

		// output the last line in the relaxation file
		if (0.0 == tot_step_cnt) {
			throw CException("Error! Zero tot_step_cnt.");
		}

		rxf << rx_cnt << WS;
		rxf << tot_step_cnt << WS << tot_wr << WS;
		rxf << tot_step_cnt/rx_cnt << WS << tot_wr/rx_cnt << WS << tot_wr/tot_step_cnt << NL;
		// close the relaxation data file
		rxf.close();

		// close the observational data file
		odf.close();

		// output the energy release profile in bar chart form
		wr_bcf.open(argv[WRBC_FN_INDEX]);
		if (!wr_bcf) {
			throw CException("Error! Cannot open energy release bar chart output file.");
		}

		double zero_pos = pf_wr.get_plot_pos(0);
		double zero_val = pf_wr.get_plot_val((unsigned long) 0);

		wr_bcf << zero_pos << WS << 0 << NL;
		wr_bcf << zero_pos << WS << zero_val << WS << NL;

		double prev_plot_val = zero_val;
		double plot_pos = 0.0, plot_val = 0.0;
		unsigned long plot_cnt = pf_wr.get_plot_cnt();
		unsigned long i = 0;

		for (i = 1; i < plot_cnt; ++i) {
			plot_pos = pf_wr.get_plot_pos(i);
			plot_val = pf_wr.get_plot_val(i);

			if (prev_plot_val != plot_val) {
				wr_bcf << plot_pos << WS;
				wr_bcf << prev_plot_val << WS;
				wr_bcf << NL;
			}

			wr_bcf << plot_pos << WS;
			wr_bcf << plot_val << WS;
			wr_bcf << NL;

			prev_plot_val = plot_val;
		}

		wr_bcf << pf_wr.get_plot_pos(plot_cnt-1) << WS << 0 << NL;

		wr_bcf.close();

		// output the energy release profile in line graph form
		wr_plf.open(argv[WRPL_FN_INDEX]);
		if (!wr_plf) {
			throw CException("Error! Cannot open energy release plot output file.");
		}

		double x = 0.0, y = 0.0;

		for (i = 0; i < plot_cnt-1; ++i) {
			plot_pos = pf_wr.get_plot_pos(i);

			x = plot_pos + ((pf_wr.get_plot_pos(i+1)-plot_pos)/2.0);
			y = pf_wr.get_plot_val(i);

			wr_plf << x << WS << y << WS << log(x) << WS << log(y) << NL;
		}

		wr_plf.close();
	}

	catch (CException& x) {
		if (rxf.is_open()) {
			rxf.clear();
			rxf.close();
		}

		if (odf.is_open()) {
			odf.clear();
			odf.close();
		}

		if (wr_bcf.is_open()) {
			wr_bcf.clear();
			wr_bcf.close();
		}

		if (wr_plf.is_open()) {
			wr_plf.clear();
			wr_plf.close();
		}

		throw x;
	}
}

/**
 * TROLE is either interrogating the properties of loops (INTERROGATE_MODE) or
 * it is shepherding a loop through a specified stability region (SHEPHERD_MODE).
 *
 * @throw CException should any of the following occur:
 *        1. invalid number of arguments;
 *        2. cannot open arguments file;
 *        3. cannot read arguments file;
 *        4. invalid mode;      
 */
int main(int argc, char** argv) {

	try {
		CGSL::random_alloc();
            
	    if (argc != 2) {
	    	throw CException("Error! Invalid number of arguments.");
		}

		char args[32][512];
		char arg[512];
		unsigned int i = 0;

		memset(args, 0, sizeof(args));
		memset(arg, 0, sizeof(arg));

        FILE *argsf = fopen(argv[ARGSFN_INDEX], "r");
		if (NULL == argsf) {
			throw CException("Error! Cannot open arguments file.");
		}               

		while (!feof(argsf)) {
			memset(arg, 0, sizeof(arg));
			fgets(arg, sizeof(arg)-1, argsf);

			if ('#' == arg[0]) {
				// this is the start of comment block
				do {
					memset(arg, 0, sizeof(arg));
					fgets(arg, sizeof(arg)-1, argsf);
				} while ('#' != arg[0]);
				// have reached end of comment block
				continue;
			}
			else if ('\n' == arg[0] && i > 0) {
				// finished reading the argument set
				char mode = args[1][0];
				if (EXAMINE_MODE == mode) {
					examine(1+i, args);
				}
				else if (EXAMINE_SMOOTH_MODE == mode) {
					examine_smooth(1+i, args);
				}
				else if (SHEPHERD_MODE == mode) {
					shepherd(1+i, args);
				}
				else {
					throw CException("Error! Invalid mode.");
				}
				i = 0;
				memset(args, 0, sizeof(args));
			}
			else {
				if (strlen(arg) > 0 && '\n' != arg[0]) {
					strncpy(&(args[1+i][0]), arg, strlen(arg)-1);
					i++;
				}
			}
		} // end of <while (!feof(argsf))> loop
	
		fclose(argsf);
		CGSL::random_free();
    } // end of try block


	catch (CException& x) {
		cout << x;
	}

	return 0;
}
