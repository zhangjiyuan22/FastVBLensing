// VBBinaryLensing v3.6 (2023)
//
// This code has been developed by Valerio Bozza (University of Salerno) and collaborators.
// Any use of this code for scientific publications should be acknowledged by a citation to:
// V. Bozza, E. Bachelet, F. Bartolic, T.M. Heintz, A.R. Hoag, M. Hundertmark, MNRAS 479 (2018) 5157
// If you use astrometry, user-defined limb darkening or Keplerian orbital motion, please cite
// V. Bozza, E. Khalouei and E. Bachelet (arXiv:2011.04780)
// The original methods present in v1.0 are described in
// V. Bozza, MNRAS 408 (2010) 2188
// Check the repository at http://www.fisica.unisa.it/GravitationAstrophysics/VBBinaryLensing.htm
// for the newest version.
//
// The code relies on the root solving algorithm by Jan Skworon and Andy Gould
// described in Skowron & Gould arXiv:1203.1034.
// Please also cite this paper if specifically relevant in your scientific publication.
// The original Fortran code available on http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/
// has been translated to C++ by Tyler M. Heintz and Ava R. Hoag (2017)
//
// GNU Lesser General Public License applies to all parts of this code.
// Please read the separate LICENSE.txt file for more details.


#ifndef __binlens
#define __binlens
#define __unmanaged

#define _L1 x1-((x1+a/2.0)/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*(x1-a/2.0)/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q) // Used in PlotCrits
#define _L2 x2-(x2/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*x2/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q)
#define _LL (y-z)+coefs[21]/(zc-coefs[20])+coefs[22]/zc //Lens equation test, y/z/zc all take secondary lens as origin, only used in NewImages()
#define _J1c coefs[21]/((zc-coefs[20])*(zc-coefs[20]))+coefs[22]/(zc*zc) //#define _J1 m1/((zc-0.5*a)*(zc-0.5*a))+m2/((zc+0.5*a)*(zc+0.5*a))
#define _J2 -2.0*(coefs[21]/((z-coefs[20])*(z-coefs[20])*(z-coefs[20]))+coefs[22]/(z*z*z))
#define _J3 6.0*(coefs[21]/((z-coefs[20])*(z-coefs[20])*(z-coefs[20])*(z-coefs[20]))+coefs[22]/(z*z*z*z))
#define _skew(p1,p2,q1,q2) p1*q2-p2*q1
#define _NP 200.0
#define __rsize 151
#define __zsize 101

#define _sign(x) ((x>0)? +1 : -1)

#include<stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

class _curve;
class _sols;
class _theta;
class complex;
struct annulus;

#ifndef __unmanaged
namespace VBBinaryLensingLibrary {

	public ref class VBBinaryLensing
#else
	class VBBinaryLensing
#endif
	{
	protected:
		int *ndatasat;
		double **tsat,***possat;
		double Mag0, corrquad, corrquad2, safedist;
		int nim0;
		double e,phi,phip,phi0,Om,inc,t0,d3,v3,GM,flagits;
		double Obj[3],rad[3],tang[3],t0old;
		double Eq2000[3],Quad2000[3],North2000[3];
		double ESPLout[__rsize][__zsize], ESPLin[__rsize][__zsize],ESPLoutastro[__rsize][__zsize], ESPLinastro[__rsize][__zsize];
		double *LDtab,*rCLDtab,*CLDtab;
		double scr2, sscr2;
		int npLD;
		bool ESPLoff, multidark;
		annulus *annlist;

		void ComputeParallax(double, double, double *);
		double LDprofile(double r);
		double rCLDprofile(double tc,annulus *,annulus *);
		double BinaryMagSafe(double s, double q, double y1, double y2, double rho, _sols **images);
		_curve *NewImages(complex,complex  *,_theta *);//, float &);
		void OrderImages(_sols *,_curve *);
		void cmplx_laguerre(complex *, int, complex *, int &, bool &);
		void cmplx_newton_spec(complex *, int, complex *, int &, bool &);
		void cmplx_laguerre2newton(complex *, int, complex *, int &, bool &, int);
		void solve_quadratic_eq(complex &, complex &, complex *);
		void solve_cubic_eq(complex &, complex &, complex &, complex *);

	public: 

		double Tol, RelTol, a1,a2, t0_par; 
		bool astrometry;
		int satellite,parallaxsystem,t0_par_fixed,nsat;
		int minannuli,nannuli,NPS,NPcrit;
		double y_1,y_2,av, therr,astrox1,astrox2;


	// Critical curves and caustic calculation
		_sols *PlotCrit(double a,double q);
		void PrintCau(double a,double q,double y1, double y2, double rho);

	// Initialization for calculations including parallax
		void SetObjectCoordinates(char *Coordinates_file, char *Directory_for_satellite_tables);

	// Magnification calculation functions.

		double BinaryMag0(double s,double q,double y1,double y2, _sols **Images);
		double BinaryMag0(double s, double q, double y1, double y2);
		double BinaryMag(double s,double q,double y1,double y2,double rho,double accuracy, _sols **Images);
		double BinaryMag(double s,double q ,double y1,double y2,double rho,double accuracy);
		double BinaryMag2(double s, double q, double y1, double y2, double rho);
		double BinaryMagDark(double s, double q, double y1, double y2, double rho,double accuracy);
		void BinaryMagMultiDark(double s, double q, double y1, double y2, double rho, double *a1_list, int n_filters, double *mag_list, double accuracy);

	// Limb Darkening control
		enum LDprofiles { LDlinear, LDquadratic, LDsquareroot, LDlog, LDuser};
		void SetLDprofile(double(*UserLDprofile)(double), int tablesampling);
		void SetLDprofile(LDprofiles);

	// ESPL functions
		void LoadESPLTable(char *tablefilename);
		double ESPLMag(double u, double rho);
		double ESPLMag2(double u, double rho);
		double ESPLMagDark(double u, double rho);
		double PSPLMag(double u);


	// New (v2) light curve functions, operating on arrays

		void PSPLLightCurve(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
		void PSPLLightCurveParallax(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
		void ESPLLightCurve(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
		void ESPLLightCurveParallax(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);

		void BinaryLightCurve(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
		void BinaryLightCurveW(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
		void BinaryLightCurveParallax(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
		void BinaryLightCurveOrbital(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, double *sep_array, int np);
		void BinaryLightCurveKepler(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, double *sep_array, int np);

		void BinSourceLightCurve(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
		void BinSourceLightCurveParallax(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
		void BinSourceLightCurveXallarap(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, double *sep_array, int np);
		void BinSourceSingleLensXallarap(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, double *y1_array2, double *y2_array2, int np);
		void BinSourceBinLensXallarap(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
        
        void BinaryMag2_Npoint(double *s, double q,  double rho, \
										double *y1s, double *y2s, \
										int np, \
										double *mags);

	// Old (v1) light curve functions, for a single calculation
		double PSPLLightCurve(double *parameters, double t);
		double PSPLLightCurveParallax(double *parameters, double t);
		double ESPLLightCurve(double *parameters, double t);
		double ESPLLightCurveParallax(double *parameters, double t);

		double BinaryLightCurve(double *parameters,double t);
		double BinaryLightCurveW(double *parameters, double t);
		double BinaryLightCurveParallax(double *parameters, double t);
		double BinaryLightCurveOrbital(double *parameters, double t);
		double BinaryLightCurveKepler(double *parameters, double t);

		double BinSourceLightCurve(double *parameters, double t);
		double BinSourceLightCurveParallax(double *parameters, double t);
		double BinSourceLightCurveXallarap(double *parameters, double t);
		double BinSourceBinLensXallarap(double *parameters, double t);
		double BinSourceSingleLensXallarap(double *parameters, double t);
		double BinSourceBinLensPOX(double* parameters, double t);

	// Skowron & Gould root calculator
		void cmplx_roots_gen(complex *, complex *, int, bool, bool);

	// Constructor and destructor

		VBBinaryLensing();
		~VBBinaryLensing();

		private:
			LDprofiles curLDprofile;
	};

	struct annulus{
		double bin;
		double cum;
		double Mag;
		double err;
		double f;
		int nim;
        double LDastrox1,LDastrox2;
		annulus *prev,*next;
	};


#ifndef __unmanaged
}
#endif


class _theta{
public: 
	double th,maxerr,Mag,errworst,astrox1,astrox2;
	_theta *prev,*next; // prev and next are pointers that point to _theta variable

	_theta(double);     // constructor: assign value to th

};

class _thetas{           /* only used in BinaryMag, not used in NewImages and OrderImages */
public: 
	_theta *first,*last; // pointers that point to _theta variable
	int length;

	_thetas(void);        // constructor: assign 0 to length
	~_thetas(void);       // destructor: delete _theta varialbes from *first to *last
	_theta *insert(double); // method: create an object of class '_theta' on heap, assign input value to attribute 'th', 
	                        //         and insert it to all elements in class '_thetas' according to 'th' value, 
							//         the attributes 'first' and 'last' of class '_thetas' may be updated,
							//         some elements' attributes 'prev' and 'next' may be updated, 
							//         return the newly inserted object's pointer
							//         (O(N) complexity for common case, where N is length)
	void remove(_theta*);   // method: remove the element pointed by the input pointer, 
	                        //         and update that element's prev/next element(if they exist)'s 'next'/'prev'
							//         (loop over the link list, thus O(N) complexity)
							//         (why not just delete 'stheta' as you already have the pointer? why border to use while loop?)


};

/*
class complex{
public:
	double re;
	double im;
	complex(double,double);
	complex(double);
	complex(void);
};
*/
class complex {
public:
	double re;
	double im;
	complex(double a, double b) {re = a; im = b; }
	complex(double a)           {re = a; im = 0; }
	complex(void)               {re = 0; im = 0; }
};

class _point{
public:
	double x1;                                    // this point's coordinate
	double x2;
	double parab,ds,dJ,parabastrox1,parabastrox2; // ds is x'^x" at this point
												  // dJ is Jacobian determinant at this point

	complex d,J2;								  // d is z'(theta) at this point
												  // J2 is Partial^2(zs_c)/Partial(z)^2 at this point

	_theta *theta;								  // pointer
	_point(double ,double,_theta *);              // constructor: assign value to x1, x2, and theta(pointer)
												  //              each _point variable corresponds to a _theta variable
												  //              each _theta variable is pointed by multiple _point variables, 
												  //              because for each _theta, there are multiple images(_point)

	_point *next,*prev;                           // pointers that point to _point variable

	double operator-(_point);                     // method: input a _point class variable, 
												  //         return distance^2 between 'current _point' and 'input _point'
};



class _curve{                                    // a _curve class variable is a linked list of _point variables
public:
	int length;
	_point *first,*last;
	_curve *next,*prev;
	_curve *partneratstart,*partneratend;
	double parabstart,parabastrox1,parabastrox2;

	_curve(_point *);                            // constructor: create a _curve class variable with one member(_point class variable),
												 //              input is pointer(p1) that points to the member, 
												 //              set length=1, first=last=p1, 
												 //              set p1->prev = p1->next = NULL,
												 //              set partner_at_start=partner_at_end=NULL

	_curve(void);                                // constructor: create a _curve class variable without any member(_point class variables), 
											     //              set length=0, first=last=NULL, partner_at_start=partner_at_end=NULL
	
	~_curve(void);								 // destructor:  delete _point varialbes from *first to *last

	_curve *divide(_point *);		// method: used only one time in OrderImages, just next to the hot spot 'while loop';
									//
									//     input a pointer(ref) that points to one member(_point class variable),
									//
									//     divide the original _curve object to two objects:
									//     	   original curve now consists of   first   to ref,                   length = l1; 
									// 		   'nc'(new curve)    consists of ref->next to original curve's last, length = original_length - l1
									//         
									//     return pointer 'nc'
									//     (use a for-loop to determine l1; 
									//      if ref doesn't point to a member, or ref==last, or length==1, or length==0, then error)

	void drop(_point *); 		    // method: input a pointer(ref) to an element of current curve
									//         drop that element out of the current curve, 
									//         but the element isn't deleted, still can be accessed by 'ref'
									//         no return valve
									// (if length==0 or 'ref' doesn't point to an element, then nothing happens)
									// (use for loop to check if 'ref' really points to an element, thus O(N), but loop is from 'last' to 'first')

	void append(double,double);     // method: create a new _point variable on heap, 
									//         and then append it to the end of the linked list (object of class _curve)
									// 
									//     constructor of _point class is called, assign input values to attributes x1, x2, 
									//     and set attribute theta(pointer) = NULL
									// 
									//     has no return value

	void append(_point *);			// method: append an existing _point variable (pointed by input pointer 'pp') 
								    //         to the end of the linked list (object of class _curve)

	void prepend(double,double);    //               method: not used
									//         (very similar to method: append)
									//         create a new _point variable on heap, 
									//         and then append it to the BEGINNING of the linked list (object of class _curve)
//	void prepend(_point *);
	_curve *join(_curve *);			// method: used only one time in OrderImages (actually also used one time in PlotCrit)
									// 		   input is a pointer to 'nc'(new curve) variable
									// 		   join all elements of 'nc' at the end of current curve
									// 		   return the pointer to current curve
									//         (set current curve's attribute 'parner_at_end' to nc's)
									//         (redirect nc->parner_at_end's attribute 'parner_at_end' from 'nc' to current curve)

	_curve *joinbefore(_curve *);   			  // method: not used
	_curve *reverse(void);						  // method: not used
	double closest(_point *,_point **); 		  // method: not used (actually used one time in PlotCrit)
	double closest2(_point *,_point **); 		  // method: not used
	void complement(_point **,int,_point **,int); // method: not used
};

class _sols{                                      //       a  _sols class variable is       a linked list of _curve variables, 
											      // while a _curve class variable is again a linked list of _point variables
public:
	int length;
	_curve *first,*last;

	_sols(void);								  // constructor: create a _sols class variable without any member(_curve class variables), 
											      //              set length=0, first=last=NULL
	
	~_sols(void);								  // destructor:  delete _curve varialbes from *first to *last
												  //              internally will use 'delete scan1 ;' for 'length' times. 
												  //
												  //              'delete scan1 ;' will call destructor of _curve class, 
												  //              which will delete all _point varialbes of the curve pointed by 'scan1'


	void drop(_curve *);						  // method: used only one time in OrderImages (near the while loop and divide)
												  //
												  //         input a pointer(ref) to an element/curve of current sols
												  //         drop that element/curve out of the current sols, 
												  //         but the element/curve isn't deleted, still can be accessed by 'ref'
												  //         no return valve
												  // (if length==0 or 'ref' doesn't point to an element/curve, then nothing happens)
												  // (use for loop to check if 'ref' really points to an element, thus O(N), but loop is from 'last' to 'first')
												  //
												  // totally same to   void _curve::drop(_point *ref)


	void append(_curve *);						  // method: append an existing _curve variable (pointed by input pointer 'cc') 
												  //         to the end of current sols (a linked list, object of class _sols)
												  //       no return value

	void prepend(_curve *);	    			      // method: not used
	void join(_sols *);			    			  // method: not used
};


#define MR 8
#define MT 10
#define MAXIT (MT*MR)
#define MAXM 30

#endif

/*
double abs(complex);
complex conj(complex);
complex sqrt(complex);
double real(complex);
double imag(complex);
complex expcmplx(complex);
complex cbrt(complex);
complex operator+(complex, complex);
complex operator-(complex, complex);
complex operator*(complex, complex);
complex operator/(complex, complex);
complex operator+(complex, double);
complex operator-(complex, double);
complex operator*(complex, double);
complex operator/(complex, double);
complex operator+(double, complex);
complex operator-(double, complex);
complex operator*(double, complex);
complex operator/(double, complex);
complex operator+(int, complex);
complex operator-(int, complex);
complex operator*(int, complex);
complex operator/(int, complex);
complex operator+(complex, int);
complex operator-(complex, int);
complex operator*(complex, int);
complex operator/(complex, int);
complex operator-(complex);
bool operator==(complex, complex);
bool operator!=(complex, complex);
*/

inline double abs(complex z) {
	return sqrt(z.re*z.re + z.im*z.im);
}

inline complex conj(complex z) {
	return complex(z.re, -z.im);
}

inline complex sqrt(complex z) {
	double md = sqrt(z.re*z.re + z.im*z.im);
	return (md>0) ? complex((sqrt((md + z.re) / 2)*((z.im>0) ? 1 : -1)), sqrt((md - z.re) / 2)) : 0.0;
}

inline double real(complex z) {
	return z.re;
}

inline double imag(complex z) {
	return z.im;
}

inline complex operator+(complex p1, complex p2) {
	return complex(p1.re + p2.re, p1.im + p2.im);
}

inline complex operator-(complex p1, complex p2) {
	return complex(p1.re - p2.re, p1.im - p2.im);
}

inline complex operator*(complex p1, complex p2) {
	return complex(p1.re*p2.re - p1.im*p2.im, p1.re*p2.im + p1.im*p2.re);
}

inline complex operator/(complex p1, complex p2) {
	double md = p2.re*p2.re + p2.im*p2.im;
	return complex((p1.re*p2.re + p1.im*p2.im) / md, (p1.im*p2.re - p1.re*p2.im) / md);
}

inline complex operator+(complex z, double a) {
	return complex(z.re + a, z.im);
}

inline complex operator-(complex z, double a) {
	return complex(z.re - a, z.im);
}

inline complex operator*(complex z, double a) {
	return complex(z.re*a, z.im*a);
}

inline complex operator/(complex z, double a) {
	return complex(z.re / a, z.im / a);
}

inline complex operator+(double a, complex z) {
	return complex(z.re + a, z.im);
}

inline complex operator-(double a, complex z) {
	return complex(a - z.re, -z.im);
}

inline complex operator*(double a, complex z) {
	return complex(a*z.re, a*z.im);
}

inline complex operator/(double a, complex z) {
	double md = z.re*z.re + z.im*z.im;
	return complex(a*z.re / md, -a*z.im / md);
}


inline complex operator+(complex z, int a) {
	return complex(z.re + a, z.im);
}

inline complex operator-(complex z, int a) {
	return complex(z.re - a, z.im);
}

inline complex operator*(complex z, int a) {
	return complex(z.re*a, z.im*a);
}

inline complex operator/(complex z, int a) {
	return complex(z.re / a, z.im / a);
}

inline complex operator+(int a, complex z) {
	return complex(z.re + a, z.im);
}

inline complex operator-(int a, complex z) {
	return complex(a - z.re, -z.im);
}

inline complex operator*(int a, complex z) {
	return complex(a*z.re, a*z.im);
}

inline complex operator/(int a, complex z) {
	double md = z.re*z.re + z.im*z.im;
	return complex(a*z.re / md, -a*z.im / md);
}

inline complex operator-(complex z) {
	return complex(-z.re, -z.im);
}

inline bool operator==(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return true;
	return false;
}

inline bool operator!=(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return false;
	return true;
}

inline complex expcmplx(complex p1) {
	double r = exp(p1.re);
	double theta = atan2(p1.im, p1.re);
	return complex(r*cos(theta), r*sin(theta));
}

inline complex cbrt(complex z) {
	complex zout;
	double r, r_cube, theta, theta_cube;
	r = abs(z);
	r_cube = pow(r, 0.333333333333);
	theta = atan2(z.im, z.re);
	theta_cube = theta / 3.;
	return 	complex(r_cube*cos(theta_cube), r_cube*sin(theta_cube));
}