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
//#define _PRINT_ERRORS2
//#define _ERRORS_ANALYTIC

#define _CRT_SECURE_NO_WARNINGS

#ifdef _WIN32
char systemslash = '\\';
#else
char systemslash = '/';
#endif

#include "VBBinaryLensingLibrary_v3p6.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>
//#include <map>

#ifndef __unmanaged
using namespace VBBinaryLensingLibrary;
#endif

/******************************************* changed *******************************************/
class _augmented_priority_queue{
public: 

	struct apq_node{
		double   maxerr ;
		_theta * stheta ;
		int      index  ;
	}; 

	struct sum_tree_node{
		double maxerr ;
		double sumerr ;
	}; 

	std::vector<apq_node> apq_array ; 
	std::vector<sum_tree_node> sum_tree_array ; 

	_augmented_priority_queue(void){						// constructor: avoid frequent reallocations by reserving capacity
		apq_array.reserve(64) ; 
		sum_tree_array.reserve(64) ; 
	}


	void push_augmented_heap(double maxerr_to_push, _theta * stheta_to_push)
	{	
		int hole_index = apq_array.size() ; 				// original_heap_size, a.k.a. new_heap_size - 1
		int sum_tree_current_index = hole_index ; 

		struct apq_node apq_node_to_push = {maxerr_to_push, stheta_to_push, hole_index} ;
		struct sum_tree_node sum_tree_node_to_push = {maxerr_to_push, maxerr_to_push} ;

		apq_array.push_back(apq_node_to_push) ;
		sum_tree_array.push_back(sum_tree_node_to_push) ;  
						
		int parent_index = (hole_index - 1) / 2 ; 
		int sum_tree_parent_index = parent_index ; 

		while (hole_index > 0 && apq_array[parent_index].maxerr < maxerr_to_push)
		{
			apq_array[hole_index] = apq_array[parent_index] ;
			hole_index = parent_index ;
			parent_index = (hole_index - 1) / 2 ;
		}
		
		while (sum_tree_current_index > 0)
		{
			sum_tree_array[sum_tree_parent_index].sumerr += maxerr_to_push ; 
			sum_tree_current_index = sum_tree_parent_index ;
			sum_tree_parent_index = (sum_tree_current_index - 1) / 2 ; 
		}

		apq_array[hole_index] = apq_node_to_push ;
	}


	void pop_then_push_augmented_heap(double maxerr_to_push, _theta * stheta_to_push)
	{
		int sum_tree_replace_index = apq_array[0].index ;

		struct apq_node apq_node_to_push = {maxerr_to_push, stheta_to_push, sum_tree_replace_index} ;
		struct sum_tree_node sum_tree_node_to_push = {maxerr_to_push, maxerr_to_push} ;

		int last_index = apq_array.size() - 1 ; 
		int hole_index = 0 ;

		while(true)
		{
			int max_child_index = 2 * hole_index + 1 ; 		// currently, assume left child is max child
			if (max_child_index > last_index)				// if even left child out of range
				break ;
			int right_child_index = max_child_index + 1 ; 	// right child
			if (right_child_index <= last_index && apq_array[right_child_index].maxerr > apq_array[max_child_index].maxerr)
				max_child_index = right_child_index ; 		// now, max child becomes right child

			if (maxerr_to_push >= apq_array[max_child_index].maxerr)
				break ; 
			apq_array[hole_index] = apq_array[max_child_index] ; 
			hole_index = max_child_index ; 
		}

		apq_array[hole_index] = apq_node_to_push ;


		int sum_tree_left_index   = 2 * sum_tree_replace_index + 1 ; 
		int sum_tree_parent_index = (sum_tree_replace_index - 1) / 2 ; 

		if (sum_tree_left_index <= last_index)
		{
			sum_tree_node_to_push.sumerr += sum_tree_array[sum_tree_left_index].sumerr ;

			if ( (sum_tree_left_index + 1) <= last_index)
			{
				sum_tree_node_to_push.sumerr += sum_tree_array[sum_tree_left_index + 1].sumerr ;
			}
		}

		sum_tree_array[sum_tree_replace_index] = sum_tree_node_to_push ; 

		if ( (sum_tree_replace_index == last_index) && (sum_tree_replace_index % 2 == 1) )
		{
			sum_tree_array[sum_tree_parent_index].sumerr = sum_tree_array[sum_tree_parent_index].maxerr + sum_tree_array[sum_tree_replace_index].sumerr ;
			sum_tree_replace_index = sum_tree_parent_index ;
			sum_tree_parent_index = (sum_tree_replace_index - 1) / 2 ; 
		}

		while(sum_tree_replace_index > 0)
		{
			sum_tree_array[sum_tree_parent_index].sumerr = sum_tree_array[sum_tree_parent_index].maxerr + sum_tree_array[2*sum_tree_parent_index+1].sumerr + sum_tree_array[2*sum_tree_parent_index+2].sumerr ;
			sum_tree_replace_index = sum_tree_parent_index ;
			sum_tree_parent_index = (sum_tree_replace_index - 1) / 2 ; 
		}
	}
};



class _priority_queue{
public: 

	struct pq_node{
		double   maxerr ;
		_theta * stheta ;
	}; 

	std::vector<pq_node> pq_array ; 

	_priority_queue(void){									// constructor: avoid frequent reallocations by reserving capacity
		pq_array.reserve(64) ; 
	}
	

	void push_heap(double maxerr_to_push, _theta * stheta_to_push)
	{	
		struct pq_node node_to_push = {maxerr_to_push, stheta_to_push} ;

		pq_array.push_back(node_to_push) ; 
		
		int hole_index = pq_array.size() - 1 ; 				// new_heap_size - 1, a.k.a. original_heap_size				
		int parent_index = (hole_index - 1) / 2 ; 

		while (hole_index > 0 && pq_array[parent_index].maxerr < maxerr_to_push)
		{
			pq_array[hole_index] = pq_array[parent_index] ;
			hole_index = parent_index ;
			parent_index = (hole_index - 1) / 2 ;
		}
		
		pq_array[hole_index] = node_to_push ;
	}


	void pop_heap()
	{
		int last_index = pq_array.size() - 1 ;
		struct pq_node node_to_adjust = pq_array[last_index] ; 
		double value = node_to_adjust.maxerr ;

		pq_array.pop_back() ;

		last_index-- ;
		int hole_index = 0 ;

		while(true)
		{
			int max_child_index = 2 * hole_index + 1 ; 		// currently, assume left child is max child
			if (max_child_index > last_index)				// if even left child out of range
				break ;
			int right_child_index = max_child_index + 1 ; 	// right child
			if (right_child_index <= last_index && pq_array[right_child_index].maxerr > pq_array[max_child_index].maxerr)
				max_child_index = right_child_index ; 		// now, max child becomes right child

			if (value >= pq_array[max_child_index].maxerr)
				break ; 
			pq_array[hole_index] = pq_array[max_child_index] ; 
			hole_index = max_child_index ; 
		}

		pq_array[hole_index] = node_to_adjust ;
	}


	void pop_then_push_heap(double maxerr_to_push, _theta * stheta_to_push)
	{
		struct pq_node node_to_push = {maxerr_to_push, stheta_to_push} ;

		int last_index = pq_array.size() - 1 ; 
		int hole_index = 0 ;

		while(true)
		{
			int max_child_index = 2 * hole_index + 1 ; 		// currently, assume left child is max child
			if (max_child_index > last_index)				// if even left child out of range
				break ;
			int right_child_index = max_child_index + 1 ; 	// right child
			if (right_child_index <= last_index && pq_array[right_child_index].maxerr > pq_array[max_child_index].maxerr)
				max_child_index = right_child_index ; 		// now, max child becomes right child

			if (maxerr_to_push >= pq_array[max_child_index].maxerr)
				break ; 
			pq_array[hole_index] = pq_array[max_child_index] ; 
			hole_index = max_child_index ; 
		}

		pq_array[hole_index] = node_to_push ;
	}
};
/*******************************************   end   *******************************************/





/******************************************* changed *******************************************/
class _skiplist_curve{							 // a _skiplist_curve class variable is a skip list of _point variables
public:
	_point * first, * last ;

	_point * head ;
	_point * last_array[max_skiplist_level+1] ;
	int Level ;

	int length_notation ;

	_skiplist_curve * next, * prev ;
	_skiplist_curve * partneratstart, * partneratend ;
	double parabstart, parabastrox1, parabastrox2 ;

	
	_skiplist_curve(_point * p1, int new_Level)				
	{											// constructor: used one time in BinaryMag() (and) two times in OrderImages(): if (nprec<npres)
												// 				create a _skiplist_curve class variable with one member(_point class variable),
												//              input is pointer(p1) that points to the member, 
												//              set length_notation=1, first=last=p1, 
												//              set p1->prev = p1->next = NULL,
												//              set partner_at_start=partner_at_end=NULL
		length_notation = 1 ;
		first = last = p1 ;
		p1->prev = p1->next = 0 ;
		partneratstart = partneratend = 0 ;

		head = new _point(0.0, 0.0, 0) ;		// head->next_array[0-max_skiplist_level] set to NULL through constructor
												//   p1->next_array[0-max_skiplist_level] should be set to all NULL previously
		for (int i=0; i<(max_skiplist_level+1); i++)
		{
			last_array[i] = head ;
		}

		// Level = 0 ;
		// while ( Level < max_skiplist_level && (rand() % 4) == 0 ) 
		// {
		// 	Level++ ;
		// } 
		Level = new_Level ;

		for (int i=0; i<(Level+1); i++)
		{
			head->next_array[i] = p1 ;
			last_array[i] = p1 ;
		}
	}


	_skiplist_curve(void)				
	{											// constructor: used one time in find_prev_then_divide()
												//				create a _skiplist_curve class variable without any member(_point class variables), 
											    //              set length_notation=0, first=last=NULL, partner_at_start=partner_at_end=NULL
		length_notation = 0 ;
		first = last = 0 ;
		partneratstart = partneratend = 0 ;

		head = new _point(0.0, 0.0, 0) ;		// head->next_array[] set to NULL through constructor

		for (int i=0; i<(max_skiplist_level+1); i++)
		{
			last_array[i] = head ;
		}

		Level = 0 ;
	}


	~_skiplist_curve(void)
	{											// destructor:  delete _point varialbes from *first to *last
		_point *scan1, *scan2 ;
		scan1 = first ;

		if (length_notation>0)
		{
			while(scan1)
			{
				scan2 = scan1->next ;
				delete scan1 ;
				scan1 = scan2 ;
			}
		}

		delete head ;
	}								 


	_skiplist_curve * join(_skiplist_curve * new_curve) 
	{	     									// method: used one time in OrderImages(): while (npres)
												// input is a pointer (new_curve) to '_skiplist_curve' variable
												// join all elements of 'new_curve' at the end of current curve
												// return the pointer to current curve

		last->next = new_curve->first ;			// as we already know length of current curve and new_curve all >= 1, 
												// we only need to change 'last' and two edge elements' 'next'/'prev'
		new_curve->first->prev = last ;
		last = new_curve->last ;
		
		length_notation = 2 ;					// length_notation=2 means length>=2, as two curves both have length>=1

		partneratend = new_curve->partneratend ;      			// set current curve's attribute 'partner_at_end' to new_curve's
		if (partneratend) partneratend->partneratend = this; 	// redirect new_curve->partner_at_end's attribute 'partner_at_end' from 'new_curve' to current curve
		// 'this' is a pointer pointed to current curve
		//
		// Suppose that you create an object named x of class A, 
		// and class A has a nonstatic member function f(). 
		// If you call the function x.f(), the keyword 'this' in the body of f() stores the address of x.

		for (int i=0; i<(new_curve->Level + 1); i++)
		{
			last_array[i]->next_array[i] = new_curve->head->next_array[i] ; 
			last_array[i]                = new_curve->last_array[i] ;
		}
		if (Level < new_curve->Level)
		{
			Level = new_curve->Level ;
		}

		new_curve->first = 0 ;                         
		new_curve->last  = 0 ;
		new_curve->length_notation = 0 ;
		delete new_curve ;						// Although original new_curve's elements have been linked to current curve, 
												// these elements still can be accessed by 'new_curve'. 
												// So directly delete new_curve will cause the deletion of these elements. 
												// So we need to break the relation between 'new_curve' and these elements before delete new_curve. 
												// Note: destructor is called through delete keyword
		return this ;
	}


	void append(_point * pp, int append_Level) 
	{ 											// method: only used one time in OrderImages()
												// 		   append an existing _point variable (pointed by input pointer 'pp') 
												//         to the end of the skiplist (object of class _skiplist_curve)
												//         (the skiplist must already have at least one element)
		pp->next = last->next ; // set to NULL
		pp->prev = last ;
		last->next = pp ;
		last = pp ;

		length_notation = 2 ;					// length_notation=2 means length>=2, as the skiplist must already have at least one element

		// int append_Level = 0 ;
		// while ( append_Level < max_skiplist_level && (rand() % 4) == 0 ) 
		// {
		// 	append_Level++ ;
		// } 

		for (int i=0; i<(append_Level + 1); i++)// pp->next_array[0-max_skiplist_level] should be set to all NULL previously
		{
			last_array[i]->next_array[i] = pp ; 
			last_array[i]                = pp ;
		}
		if (Level < append_Level)
		{
			Level = append_Level ;
		}
	}


	void append(double x1, double x2, int append_Level) 
	{ 											// method: only used one time in BinaryMag()
												// 		   create a new _point variable on heap, 
												//         and then append it to the end of the skiplist (object of class _skiplist_curve)
												// has no return value
												//         (the skiplist must already have at least one element)
		_point * pp ;
		pp = new _point(x1, x2, 0) ; 			// create an object of class _point on heap, the pointer to that object is assigned to 'pp'
												// constructor is called, assign value to attributes x1, x2, 
												// and set attribute theta(pointer) = NULL
												// (note each _point variable corresponds to a _theta variable, but now has not pointed to its _theta variable)
												//
												// pp->next_array[] set to NULL through constructor
		last->next = pp ;
		pp->prev = last ;
		last = pp ;
		pp->next = 0 ;

		length_notation = 2 ;					// length_notation=2 means length>=2, as the skiplist must already have at least one element

		// int append_Level = 0 ;
		// while ( append_Level < max_skiplist_level && (rand() % 4) == 0 ) 
		// {
		// 	append_Level++ ;
		// } 

		for (int i=0; i<(append_Level + 1); i++)// pp->next_array[0-max_skiplist_level] already set to all NULL through constructor
		{
			last_array[i]->next_array[i] = pp ; 
			last_array[i]                = pp ;
		}
		if (Level < append_Level)
		{
			Level = append_Level ;
		}
	}


	_skiplist_curve * find_prev_then_divide(double th)
	{											// method: only used one time in OrderImages()
												// 		   scurve->first->theta <= 'th' <= scurve->last->theta already statisfied when using this method
												// 		   as 'th' is a unique value, scurve->first->theta < 'th' < scurve->last->theta is thus statisfied
		_point * current = head ;  
		_point * update_array[max_skiplist_level+1] ; 

		_skiplist_curve * new_curve ; 

		for (int i=0; i<(max_skiplist_level+1); i++)
		{
			update_array[i] = head ;
		}        	
		
		for (int i = Level; i >= 0; i--) 		
		{
			while (current->next_array[i] && current->next_array[i]->theta->th < th) 
			{
				current = current->next_array[i] ; 
			}
			update_array[i] = current ;
		}


		new_curve = new _skiplist_curve() ;		// create another object of class _skiplist_curve using 'new' keyword, 
												// memory of the object is allocated on heap and the pointer to that object is assigned to new_curve, 
												// constructor is called (create a _skiplist_curve class variable without any member).
												// (as the object is on heap, we can return it's pointer;
												// if the object is on stack, the returned pointer will be a dangling pointer)

		new_curve->first = current->next ;  	// new_curve consists of current->next to original curve's last
		new_curve->first->prev = 0 ;
		new_curve->last = last ;

		new_curve->length_notation = (new_curve->first == new_curve->last)? 1: 2 ; 	// length_notation=1 means length=1
																					// length_notation=2 means length>=2
		new_curve->partneratend = partneratend ; 
		if (partneratend) partneratend->partneratend = new_curve ; 	// if partner_at_end really points to a _skiplist_curve variable, 
																	// then redirect *partner_at_end's attribute 'partner_at_end' from original curve to 'new_curve'
																	// (seems can direct endlessly, like partner_at_end -> partner_at_end -> ... partner_at_end -> ...)

		last = current ;						// change original curve, original curve now consists of first to current
		current->next = 0 ;		

		length_notation = (first == last)? 1: 2 ; 									// length_notation=1 means length=1
																					// length_notation=2 means length>=2
		partneratend = 0 ;       // set from xxx to NULL


		int i = 0 ; 
		//while (i < () && update_array[i]->next_array[i]) 
		while ( i<(max_skiplist_level+1) && update_array[i]->next_array[i] ) 
		{
			new_curve->head->next_array[i] =  update_array[i]->next_array[i] ;
			new_curve->last_array[i] = last_array[i] ; 
			i++ ;
		}
		new_curve->Level = i-1 ;
		
		for (int k=0; k<(new_curve->Level+1); k++)
		{
			last_array[k] = update_array[k] ;
			update_array[k]->next_array[k] = 0 ;	
		}

		int j = 0 ;
		//while (update_array[j] != head)
		while ( j<(max_skiplist_level+1) && update_array[j] != head)
		{
			j++ ;
		}
		Level = j-1 ;


		return new_curve ;
	}
}; 


class _sols_for_skiplist_curve{                 //       a  _sols_for_skiplist_curve class variable is a linked list of _skiplist_curve variables, 
												// while a _skiplist_curve class variable           is a skip list   of _point variables
public:
	int length ;
	_skiplist_curve * first, * last ;


	_sols_for_skiplist_curve(void)				// constructor: create a _sols class variable without any member(_curve class variables), 
											    //              set length=0, first=last=NULL
	{
		length = 0 ;
		first = last = 0 ;
	}				  
	

	~_sols_for_skiplist_curve(void)				// destructor:  delete _skiplist_curve varialbes from *first to *last
												//              internally will use 'delete scan1 ;' for 'length' times. 
												//
												//              'delete scan1 ;' will call destructor of _skiplist_curve class, 
												//              which will delete all _point varialbes of the curve pointed by 'scan1'
	{
		_skiplist_curve * scan1, * scan2 ;
		scan1 = first;
		while (scan1) 
		{
			scan2 = scan1->next;
			delete scan1;           // 'delete scan1' will call destructor of _skiplist_curve class, 
									// which will delete all _point varialbes of the curve pointed by scan1
			scan1 = scan2;
		}
	}			  


	void drop(_skiplist_curve * ref) 
	{ 											// method: used only one time in OrderImages(): if-if, which near the while loop and divide
												//
												//         input a pointer(ref) to an element/curve of current sols
												//         drop that element/curve out of the current sols, 
												//         but the element/curve isn't deleted, still can be accessed by 'ref'
												//         no return valve
												// (if length==0 or 'ref' doesn't point to an element/curve, then nothing happens)
												// (use for loop to check if 'ref' really points to an element, thus O(N), but loop is from 'last' to 'first')
												//
												// totally same to   void _curve::drop(_point *ref)
		_skiplist_curve * scan  ;
		if (length) {
			for (scan = last; scan && (scan != ref); scan = scan->prev);
			if (scan) {
				if (length == 1) {
					first = last = 0;
				}
				else {
					if (ref->prev) {
						ref->prev->next = ref->next;
						if (ref == last) {
							last = ref->prev;
						}
					}
					if (ref->next) {
						ref->next->prev = ref->prev;
						if (ref == first) {
							first = ref->next;
						}
					}
				}
				length--;
			}
		}
	}						  


	void append(_skiplist_curve * cc) 
	{											// method: used one time  in BinaryMag() (and) 
												//              two times in OrderImages(): if (nprec<npres) (and) 
												//              two times in OrderImages(): if (npres<nfoll)
												//		   append an existing _curve variable (pointed by input pointer 'cc') 
												//         to the end of current sols (a linked list, object of class _sols)
												//       no return value
		if (length == 0) {
			first = cc;
			last = cc;
			cc->prev = 0;
		}
		else {
			last->next = cc;
			cc->prev = last;
			last = cc;
		}
		cc->next = 0;
		length++;
	}					
}; 
/*******************************************   end   *******************************************/




//////////////////////////////
//////////////////////////////
////////Constructor and destructor
//////////////////////////////
//////////////////////////////

VBBinaryLensing::VBBinaryLensing() {
	Obj[0] = -0.0397317;
	Obj[1] = 0.998164;
	Obj[2] = -0.045714;
	// reference is ecliptic with x-axis toward the equinox.
	// axial tilt at J2000 is 23:26:21.406 from JPL fundamental ephemeris
	Eq2000[0] = 1;
	Eq2000[1] = Eq2000[2] = Quad2000[0] = North2000[0] = 0;
	Quad2000[1] = 0.9174820003578725;
	Quad2000[2] = -0.3977772982704228;
	North2000[1] = 0.3977772982704228;
	North2000[2] = 0.9174820003578725;
	t0old = 0.;
	Tol = 1.e-2;
	RelTol = 0;
	tsat = 0;
	possat = 0;
	nsat = 0;
	ndatasat = 0;
	satellite = 0;
	parallaxsystem = 0;
	t0_par_fixed = -1;
	t0_par = 7000;
	minannuli = 1;
	curLDprofile = LDlinear;
	a1 = 0;
	npLD = 0;
	LDtab = rCLDtab = CLDtab=0;
	Mag0 = 0;
	NPcrit = 200;
	ESPLoff = true;
	multidark = false;
    astrometry=false;
}

VBBinaryLensing::~VBBinaryLensing() {
	if (nsat) {
		for (int i = 0; i<nsat; i++) {
			for (int j = 0; j<ndatasat[i]; j++) free(possat[i][j]);
			free(tsat[i]);
			free(possat[i]);
		}
		free(tsat);
		free(possat);
		free(ndatasat);
	}
	if (npLD > 0) {
		free(LDtab);
		free(rCLDtab);
	}
}


//////////////////////////////
//////////////////////////////
////////Critical curves and caustics
//////////////////////////////
//////////////////////////////


_sols *VBBinaryLensing::PlotCrit(double a1, double q1) {
	complex  a, q, ej, zr[4], x1, x2;
	_sols *CriticalCurves;
	_curve *Prov, *Prov2, *isso;
	_point *pisso;
	double SD, MD, CD, centeroffset;

	a = complex(a1, 0.0);
	q = complex(q1, 0.0);
	centeroffset = a1 / 2.0*(1.0 - q1) / (1.0 + q1);

	CriticalCurves = new _sols;
	for (int i = 0; i<4; i++) {
		Prov = new _curve;
		CriticalCurves->append(Prov);
	}

	for (int j = 0; j<NPcrit; j++) {
		ej = complex(cos(2 * j*M_PI / NPcrit), -sin(2 * j*M_PI / NPcrit));
		complex  coefs[5] = { a*a / 16.0*(4.0 - a*a*ej)*(1.0 + q),a*(q - 1.0),(q + 1.0)*(1.0 + a*a*ej / 2.0),0.0,-(1.0 + q)*ej };
		cmplx_roots_gen(zr, coefs, 4, true, true);
		if (j > 0) {
			Prov2 = new _curve();
			for (int i = 0; i < 4; i++) {
				Prov2->append(zr[i].re + centeroffset, zr[i].im);
			}
			for (Prov = CriticalCurves->first; Prov;  Prov = Prov->next) {
				Prov2->closest(Prov->last, &pisso);
				Prov2->drop(pisso);
				Prov->append(pisso);
			}
		}
		else {
			Prov = CriticalCurves->first;
			for (int i = 0; i < 4; i++) {
				Prov->append(zr[i].re + centeroffset, zr[i].im);
				Prov = Prov->next;
			}
		}
	}

	Prov = CriticalCurves->first;
	while (Prov->next) {
		SD = *(Prov->first) - *(Prov->last);
		MD = 1.e100;
		isso = 0;
		for (Prov2 = Prov->next; Prov2; Prov2 = Prov2->next) {
			CD = *(Prov2->first) - *(Prov->last);
			if (CD<MD) {
				MD = CD;
				isso = Prov2;
			}
		}
		if (MD<SD) {
			CriticalCurves->drop(isso);
			Prov->join(isso);
		}
		else {
			Prov = Prov->next;
		}
	}

	// Caustics

	for (Prov = CriticalCurves->last; Prov; Prov = Prov->prev) {
		Prov2 = new _curve;
		for (_point *scanpoint = Prov->first; scanpoint; scanpoint = scanpoint->next) {
			x1 = complex(scanpoint->x1 - centeroffset, 0.0);
			x2 = complex(scanpoint->x2, 0.0);
			Prov2->append(real(_L1) + centeroffset, real(_L2));
		}
		CriticalCurves->append(Prov2);
	}
	return CriticalCurves;
}

void VBBinaryLensing::PrintCau(double a, double q, double y1, double y2, double rho) {
	_sols *CriticalCurves;
	_curve *scancurve;
	_point *scanpoint;
	FILE *f;
	int ncc;

	CriticalCurves = PlotCrit(a, q);
	f = fopen("outcurves.causticdata", "w");
	fprintf(f, "%.16lf %.16lf %.16lf\n", y1,y2,rho);
	ncc = CriticalCurves->length / 2;
	scancurve = CriticalCurves->first;
	for (int i = 0; i<2*ncc; i++) {
	//	scancurve = scancurve->next;
	//}

	//for (int i = 0; i<ncc; i++) {
		fprintf(f, "Curve: %d\n", i + 1);
		for (scanpoint = scancurve->first; scanpoint; scanpoint = scanpoint->next) {
			fprintf(f, "%.16lf %.16lf\n", scanpoint->x1, scanpoint->x2);
		}
		scancurve = scancurve->next;
	}
	fclose(f);
	delete CriticalCurves;
}


//////////////////////////////
//////////////////////////////
////////Parallax settings and computation
//////////////////////////////
//////////////////////////////


void VBBinaryLensing::SetObjectCoordinates(char *modelfile, char *sateltabledir) {
	double RA, Dec, dis, hr, mn, sc, phiprec;
	char filename[256];
	int ic;
	FILE *f;

	if (nsat) {
		for (int i = 0; i<nsat; i++) {
			for (int j = 0; j<ndatasat[i]; j++) free(possat[i][j]);
			free(tsat[i]);
			free(possat[i]);
		}
		free(tsat);
		free(possat);
		free(ndatasat);
	}

	f = fopen(modelfile, "r");
	if (f != 0) {
		fscanf(f, "%lf:%lf:%lf", &hr, &mn, &sc);
		RA = (hr + mn / 60 + sc / 3600)*M_PI / 12,
			fscanf(f, "%lf:%lf:%lf", &hr, &mn, &sc);
		Dec = (fabs(hr) + mn / 60 + sc / 3600)*M_PI / 180;
		if (hr < 0) Dec = -Dec;

		for (int i = 0; i < 3; i++) {
			Obj[i] = (cos(RA)*cos(Dec)*Eq2000[i] + sin(RA)*cos(Dec)*Quad2000[i] + sin(Dec)*North2000[i]);
			rad[i] = Eq2000[i];
			tang[i] = North2000[i];
		}
		fclose(f);


		// Looking for satellite table files in the specified directory
		sprintf(filename, "%s%csatellite*.txt", sateltabledir, systemslash);
		nsat = 0;
		for (unsigned char c = 32; c < 255; c++) {
			filename[strlen(filename) - 5] = c;
			f = fopen(filename, "r");
			if (f != 0) {
				nsat++;
				fclose(f);
			}
		}


		tsat = (double **)malloc(sizeof(double *)*nsat);
		possat = (double ***)malloc(sizeof(double **)*nsat);
		ndatasat = (int *)malloc(sizeof(int)*nsat);

		// Reading satellite table files
		ic = 0;
		for (unsigned char c = 32; c < 255; c++) {
			filename[strlen(filename) - 5] = c;
			f = fopen(filename, "r");
			if (f != 0) {
				int flag2 = 0;
				long startpos = 0;
				char teststring[1000];
				ndatasat[ic] = 1;

				// Finding start of data
				while (!feof(f)) {
					fscanf(f, "%s", teststring);
					if (!feof(f)) {
						fseek(f, 1, SEEK_CUR);
						teststring[5] = 0;
						if (strcmp(teststring, "$$SOE") == 0) {
							flag2 = 1;
							break;
						}
					}
				}
				// Finding end of data
				if (flag2) {
					flag2 = 0;
					startpos = ftell(f);
					while (!feof(f)) {
						fscanf(f, "%[^\n]s", teststring);
						if (!feof(f)) {
							fseek(f, 1, SEEK_CUR);
							teststring[5] = 0;
							if (strcmp(teststring, "$$EOE") == 0) {
								flag2 = 1;
								break;
							}
							else {
								ndatasat[ic]++;
							}
						}
					}
				}

				// Allocating memory according to the length of the table
				tsat[ic] = (double *)malloc(sizeof(double)*ndatasat[ic]);
				possat[ic] = (double **)malloc(sizeof(double *)*ndatasat[ic]);
				for (int j = 0; j < ndatasat[ic]; j++) {
					possat[ic][j] = (double *)malloc(sizeof(double) * 3);
				}
				ndatasat[ic]--;

				// Reading data
				if (f) {
					fseek(f, startpos, SEEK_SET);
					for (int id = 0; id < ndatasat[ic]; id++) {

						if (fscanf(f, "%lf %lf %lf %lf %lf", &(tsat[ic][id]), &RA, &Dec, &dis, &phiprec) == 5) {
							tsat[ic][id] -= 2450000;
							RA *= M_PI / 180;
							Dec *= M_PI / 180;
							for (int i = 0; i < 3; i++) {
								possat[ic][id][i] = dis*(cos(RA)*cos(Dec)*Eq2000[i] + sin(RA)*cos(Dec)*Quad2000[i] + sin(Dec)*North2000[i]);
							}
						}
						else {
							ndatasat[ic] = id;
							break;
						}
					}
					fclose(f);
				}

				ic++;
			}
		}
		if (t0_par_fixed == -1) t0_par_fixed = 0;
	}else {
		printf("\nFile not found!\n");
	}
}

void VBBinaryLensing::ComputeParallax(double t, double t0, double *Et) {
	static double a0 = 1.00000261, adot = 0.00000562; // Ephemeris from JPL website 
	static double e0 = 0.01671123, edot = -0.00004392;
	static double inc0 = -0.00001531, incdot = -0.01294668;
	static double L0 = 100.46457166, Ldot = 35999.37244981;
	static double om0 = 102.93768193, omdot = 0.32327364;
	static double deg = M_PI / 180;
	static double a, e, inc, L, om, M, EE, dE, dM;
	static double x1, y1, vx, vy, Ear[3], vEar[3];
	static double Et0[2], vt0[2], r, sp, ty, Spit;
	int c = 0, ic;

	if (t0_par_fixed == 0) t0_par = t0;
	if (t0_par_fixed == -1) {
		printf("\nUse SetObjectCoordinates to input target coordinates");
	}
	else {

		if (t0_par != t0old) {
			t0old = t0_par;
			ty = (t0_par - 1545) / 36525.0;

			a = a0 + adot*ty;
			e = e0 + edot*ty;
			inc = (inc0 + incdot*ty)*deg;
			L = (L0 + Ldot*ty)*deg;
			om = (om0 + omdot*ty)*deg;

			M = L - om;
			M -= floor((M + M_PI) / (2 * M_PI)) * 2 * M_PI;

			EE = M + e*sin(M);
			dE = 1;
			while (fabs(dE) > 1.e-8) {
				dM = M - (EE - e*sin(EE));
				dE = dM / (1 - e*cos(EE));
				EE += dE;
			}
			x1 = a*(cos(EE) - e);
			y1 = a*sqrt(1 - e*e)*sin(EE);
			//		r=a*(1-e*cos(EE));
			vx = -a / (1 - e*cos(EE))*sin(EE)*Ldot*deg / 36525;
			vy = a / (1 - e*cos(EE))*cos(EE)*sqrt(1 - e*e)*Ldot*deg / 36525;

			Ear[0] = x1*cos(om) - y1*sin(om);
			Ear[1] = x1*sin(om)*cos(inc) + y1*cos(om)*cos(inc);
			Ear[2] = x1*sin(om)*sin(inc) + y1*cos(om)*sin(inc);
			vEar[0] = vx*cos(om) - vy*sin(om);
			vEar[1] = vx*sin(om)*cos(inc) + vy*cos(om)*cos(inc);
			vEar[2] = vx*sin(om)*sin(inc) + vy*cos(om)*sin(inc);

			sp = 0;
			switch (parallaxsystem) {
			case 1:
				for (int i = 0; i < 3; i++) sp += North2000[i] * Obj[i];
				for (int i = 0; i < 3; i++) rad[i] = -North2000[i] + sp*Obj[i];
				break;
			default:
				for (int i = 0; i < 3; i++) sp += Ear[i] * Obj[i];
				for (int i = 0; i < 3; i++) rad[i] = Ear[i] - sp*Obj[i];
				break;
			}

			r = sqrt(rad[0] * rad[0] + rad[1] * rad[1] + rad[2] * rad[2]);
			rad[0] /= r;
			rad[1] /= r;
			rad[2] /= r;
			tang[0] = rad[1] * Obj[2] - rad[2] * Obj[1];
			tang[1] = rad[2] * Obj[0] - rad[0] * Obj[2];
			tang[2] = rad[0] * Obj[1] - rad[1] * Obj[0];

			Et0[0] = Et0[1] = vt0[0] = vt0[1] = 0;
			for (int i = 0; i < 3; i++) {
				Et0[0] += Ear[i] * rad[i];
				Et0[1] += Ear[i] * tang[i];
				vt0[0] += vEar[i] * rad[i];
				vt0[1] += vEar[i] * tang[i];
			}
		}

		ty = (t - 1545) / 36525.0;

		a = a0 + adot*ty;
		e = e0 + edot*ty;
		inc = (inc0 + incdot*ty)*deg;
		L = (L0 + Ldot*ty)*deg;
		om = (om0 + omdot*ty)*deg;

		M = L - om;
		M -= floor((M + M_PI) / (2 * M_PI)) * 2 * M_PI;

		EE = M + e*sin(M);
		dE = 1;
		while (dE > 1.e-8) {
			dM = M - (EE - e*sin(EE));
			dE = dM / (1 - e*cos(EE));
			EE += dE;
		}
		x1 = a*(cos(EE) - e);
		y1 = a*sqrt(1 - e*e)*sin(EE);
		//	r=a*(1-e*cos(EE));

		Ear[0] = x1*cos(om) - y1*sin(om);
		Ear[1] = x1*sin(om)*cos(inc) + y1*cos(om)*cos(inc);
		Ear[2] = x1*sin(om)*sin(inc) + y1*cos(om)*sin(inc);
		Et[0] = Et[1] = 0;
		for (int i = 0; i < 3; i++) {
			Et[0] += Ear[i] * rad[i];
			Et[1] += Ear[i] * tang[i];
		}
		Et[0] += -Et0[0] - vt0[0] * (t - t0_par);
		Et[1] += -Et0[1] - vt0[1] * (t - t0_par);

		if (satellite > 0 && satellite <= nsat) {
			if (ndatasat[satellite - 1] > 2) {
				int left, right;
				if (t < tsat[satellite - 1][0]) {
					ic = 0;
				}
				else {
					if (t > tsat[satellite - 1][ndatasat[satellite - 1] - 1]) {
						ic = ndatasat[satellite - 1] - 2;
					}
					else {
						left = 0;
						right = ndatasat[satellite - 1] - 1;
						while (right - left > 1) {
							ic = (right + left) / 2;
							if (tsat[satellite - 1][ic] > t) {
								right = ic;
							}
							else {
								left = ic;
							}
						}
						ic = left;
					}
				}
				ty = t - tsat[satellite - 1][ic];
				for (int i = 0; i < 3; i++) {
					Spit = possat[satellite - 1][ic][i] * (1 - ty) + possat[satellite - 1][ic + 1][i] * ty;
					Et[0] += Spit*rad[i];
					Et[1] += Spit*tang[i];
				}
			}
		}
	}
}

//////////////////////////////
//////////////////////////////
////////Basic magnification functions
//////////////////////////////
//////////////////////////////


// Static variables in a Function: 
// Even if the function is called multiple times, space for the static variable is allocated only once and live until the program terminates. 
// The value of the variable in the previous call gets carried through the next function call. 
// (Static variables are allocated memory in the data segment, not the stack segment.)
// (Static variables (like global variables) are initialized as 0 if not initialized explicitly. )
//
// Notice static local variable is not global variable! 
// static local variable is only valid in the function where it's declared, although it's value is carried through different function calls. 
// (Even though their lifespan is till the termination of the program, their scope is limited to the block in which they are declared.)
/******************************************* changed *******************************************/
double VBBinaryLensing::BinaryMag0(double a1, double q1, double y1v, double y2v, _sols_for_skiplist_curve **Images) {
//double VBBinaryLensing::BinaryMag0(double a1, double q1, double y1v, double y2v, _sols **Images) {
/*******************************************   end   *******************************************/
	// input 'Images' is the address of a data-segment pointer 'images' which points to _sols class variable
	// thus we can modify the value of pointer 'images' inside the function (i.e. modify the static variable 'images')
	// 
	// or 'images' is copied from actual parameter to formal parameter, and can only modify the formal parameter inside the function
	//
	// notice 'images' is static local variable instead of global variable. 
	// So it is only valid in the function where it's declared, thus we need to explictly pass it/it's address to inner function. 
	// If it's a global variable, we can directly call it in inner function without passing it. 

	// checkpoint3 (jump between BinaryMag0 and BinaryMag by finding this string in VScode)
	static complex a, q, m1, m2, y;
	static double av = -1.0, qv = -1.0;
	static complex coefs[24], d1, d2, dy, dJ, dz;
	static double Mag, Ai;
    
	static _theta *stheta;
	/******************************************* changed *******************************************/
	//static _curve *Prov, *Prov2; // there is also a 'Prov' in NewImages(), but they are different static local variables. 
	static _curve * Prov ;
	static _skiplist_curve * Prov2 ;
	/*******************************************   end   *******************************************/
	static _point *scan1, *scan2;

	//float time ;

	Mag = Ai = -1.0;
	stheta = new _theta(-1.); 	// create a _theta class variable on heap with attribute 'th' = -1.0, 
								// pointer to that object is 'stheta'
	
	if ((a1 != av) || (q1 != qv)) {
		av = a1;
		qv = q1;
		if (q1<1) {             
			a = complex(-a1, 0);
			q = complex(q1, 0);
		}
		else {
			a = complex(a1, 0);
			q = complex(1 / q1, 0);
		}
		m1 = 1.0 / (1.0 + q);
		m2 = q*m1;
		/********************************************************************************/
		// when q1<1, 
		//  m1=1/(1+q1) which is large         m2=q1/(1+q1) which is small
		//    * ---*---------------------------- *
		// (-s,0)  Center of Mass              (0,0)
		//
		// when q1>=1, 
		//  								   m2=1/(1+q1) which is small         m1=q1/(1+q1) which is large
		//    									 * ----------------------------*--- *
		// 									   (0,0)              Center of Mass  (s,0)
		/********************************************************************************/
		// BinaryMag0's input coordinate (y1v,y2v) takes the mass center as origin, 
		// then BinaryMag0 calls function NewImages and still use (y1v,y2v) as input
		//
		// then in function NewImages, y = (y1v,y2v) + coefs[11] becomes the zs used in the coefficient expression,
		//                               = (y1v,y2v) + (-s *  1/(1+q1), 0), when q1<1
		//                               = (y1v,y2v) + ( s * q1/(1+q1), 0), when q1>=1
    	// for example, when q1<1,  the source on second  lens used to be ( s* 1/(1+q1), 0), but now becomes ( 0, 0)
    	//              		    the source on primary lens used to be (-s*q1/(1+q1), 0), but now becomes (-s, 0)
		//              when q1>=1, the source on second  lens used to be ( s* 1/(1+q1), 0), but now becomes ( s, 0)
    	//              		    the source on primary lens used to be (-s*q1/(1+q1), 0), but now becomes ( 0, 0)
		// (here primary and second lens defined as Msecond/Mprimary = q1, no matter q1<1 or q1>=1)
		// (or we can say the left/right lens is primary/second lens)
		// (while m1 is always the massive lens and m2 is the small lens; m2 must be on origin)
		coefs[20] = a;
		coefs[21] = m1;
		coefs[22] = m2;
		coefs[6] = a*a;
		coefs[7] = coefs[6] * a;
		coefs[8] = m2*m2;
		coefs[9] = coefs[6] * coefs[8];
		coefs[10] = a*m2;
		coefs[11] = a*m1;
		coefs[23] = 0;

	}
	y = complex(y1v, y2v);

	// checkpoint1 (jump between BinaryMag0 and NewImages by finding this string in VScode)

	/******************************************* changed *******************************************/
	(*Images) = new _sols_for_skiplist_curve ;
	//(*Images) = new _sols;      // create a _sols class variable without any member (_curve class variables) 
								// and make static pointer 'images'(i.e. *Images) point to that object
	/*******************************************   end   *******************************************/

	corrquad = corrquad2 = 0;   // declaration is in .h (in class definition), 
								// protected thus cannot be accessed from outside the class (i.e. obj.coorquad)

	safedist = 10;				// same as corrquad and corrquad2

	//Prov = NewImages(y, coefs, stheta, time);
	Prov = NewImages(y, coefs, stheta); // pass a static local variable 'stheta' as actual parameter, 
										// then it's copied from actual parameter to formal parameter 'theta' inside NewImages(), 
										// and can only modify the formal parameter inside NewImages()
										// (but obviously we don't need to modify 'stheta' or 'theta', 
										// we only need to check/modify theta->attribute)
										//
										// return a pointer to a _curve variable which is created on heap during NewImages()
										// then assign the returned pointer to static local variable 'Prov'; 
										// this _curve variable is a linked list of good roots (_point variables) at this zs
										// this _curve variable can contain 0, 3, or 5 elements
										// make all elements (_point variables) point to the same _theta variable
										

	if (Prov->length == 0) {			// if this _curve variable contains 0 element, then error
										// only happens when: 
										// (if) one (or more) good roots' dJ is too close to 0 (too close to 0 means -1.0e-7 < dJ < 1.0e-7 )
										// or wrongly polishing a positive parity image to a negative parity image, or vice versa, 
										// causing the summation of good root's parity != -1 ; 
										// and total magnification is not enough close to 1
		delete Prov;
		delete stheta;
		return -1;
	}
	if (q.re < 0.01) {							// planetary test (only needed for q<0.01)
		safedist = y1v + coefs[11].re-1/a.re;	// y1v 				  		   is zs_x centering on mass center, 
												// y1v + coefs[11].re 		   is zs_x centering on the less massive lens, 
												// 1/a.re is (approximately) planet caustic's x coordinate centering on the less massive lens, 
												// ( the accurate expression is 1/a.re/(1+q.re) )										
		safedist *= safedist;
		safedist += y2v*y2v - 36 * q1/(a1*a1);	// safedist is thus distance square between zs and planet caustic - 4 * delta_pc**2, 
												// where delta_pc is defined in Equation(50) of Bozza et al., 2018
												// Question: But why use q1 and a1 as q1 can >1 ?
												//           should be q.re/(a.re*a.re)?
	} 
	Mag = 0.;					
    astrox1=0.;					// declaration is in .h (in class definition), public
    astrox2=0.;
	nim0 = 0;					// declaration is in .h (in class definition), 
								// protected thus cannot be accessed from outside the class (i.e. obj.nim0)
								
	for (scan1 = Prov->first; scan1; scan1 = scan2) {	// loop through all 3 or 5 good roots in _curve variable
		scan2 = scan1->next;		

		/******************************************* changed *******************************************/
		Prov2 = new _skiplist_curve(scan1, 0);						// create an object of class _curve with one member(_point class variable),
												 		// input is pointer(scan1) that points to the member; 
														// pointer to that object is assigned to static local variable 'Prov2'
		/*******************************************   end   *******************************************/

		//Prov2->append(scan1->x1, scan1->x2);
		//Prov2->last->theta = stheta;
		//Prov2->last->d = Prov2->first->d;
		//Prov2->last->dJ = Prov2->first->dJ;
		//Prov2->last->J2 = Prov2->first->J2;
		//Prov2->last->ds = Prov2->first->ds;
		(*Images)->append(Prov2);						// append the newly created _curve variable (pointed by 'Prov2') 
												  		// to the end of sols (a linked list, object of class _sols)
														// which is pointed by *Images (i.e. static local pointer 'images' which is out of the scope of current function)
        
		Ai=fabs(1 / scan1->dJ);							// Ai is magnification of current image, Mag is total magnification
		Mag += Ai;										
		if(astrometry){									// astrox1,astrox2 is magnification weighted SUM of all images' positions
														// centering at the less massive lens
			astrox1 +=scan1->x1*Ai;
			astrox2 +=(scan1->x2)*Ai;
		}
		nim0++;											// number of images

	}													// now we get a _sols variable(pointed by *Images), containing 3 or 5 _curve variables, 
														// and each _curve variable contains only one _point variable
	Prov->length = 0;
	delete Prov;										// delete the _curve variable created during NewImages() after setting length to 0, 
														// delete will first call destructor of _curve class
														// (but as length is set to 0, _point variables won't be deleted), 
														// then free the space of this _curve variable
	delete stheta;
    if(astrometry){										// astrox1,astrox2 now becomes magnification weighted AVERAGE of all images' positions
														// centering at the mass center
		astrox1 /= (Mag);
		astrox1 -=coefs[11].re;
		astrox2 /= (Mag);							
    }
	NPS = 1;											// number of points
	//return time;
	return Mag;				// return the magnification, meanwhile (astrox1,astrox2) can be accessed from outside the class (i.e. obj.astrox1), 
							// corrquad, corrquad2, safedist can be accessed only inside the class, 
							// and we create a _sols class variable and make static pointer 'images'(i.e. *Images) point to that object, 
							// finally this _sols variable(pointed by *Images) contains 3 or 5 _curve variables, 
							// and each _curve variable contains only one _point variable
							// (all _point variables used to point to the same _theta variable, 
							// but the _theta variable has been deleted, makes (_point variable).theta a dangling pointer)
							// (if the _curve variable created during NewImages() contains 0 element, then return -1)
}

double VBBinaryLensing::BinaryMag0(double a1, double q1, double y1v, double y2v) {
	/******************************************* changed *******************************************/
	//static _sols *images;  	// create a pointer 'images' (points to _sols class variable) on data segment
	static _sols_for_skiplist_curve * images ;
	/*******************************************   end   *******************************************/
	static double mag;
	mag = BinaryMag0(a1, q1, y1v, y2v, &images); // pass the address of pointer 'images' as the actual parameter, 
												 // enabling to modify 'images' inside BinaryMag0
												 // 
												 // assign the returned magnification to static local variable 'mag'

	delete images;			// delete the _sols class object pointed by 'images'
							// delete will first call destructor of _sols class, 
							// 				which will delete _curve varialbes from *first to *last
							//              internally will use 'delete scan1 ;' for 'length' times. 
							//              'delete scan1 ;' will call destructor of _curve class, 
							//              which will delete all _point varialbes of the curve pointed by 'scan1'
							// then free the space of this _sols variable
							// 
							// (if we didn't modify 'images' previously, 
							// now we can't delete the _sols class object created inside BinaryMag0, 
							// unless we returned the pointer to that object through BinaryMag0, just like what happens to 'mag')

	return mag;				// return the magnification, meanwhile (astrox1,astrox2) can be accessed from outside the class (i.e. obj.astrox1), 
							// corrquad, corrquad2, safedist can be accessed only inside the class
							// (if the _curve variable created during NewImages() contains 0 element, then return -1)
							// 
							// no dynamically allocated memory left
}

/******************************************* changed *******************************************/
double VBBinaryLensing::BinaryMagSafe(double s, double q, double y1v, double y2v, double RS, _sols_for_skiplist_curve **images) {
//double VBBinaryLensing::BinaryMagSafe(double s, double q, double y1v, double y2v, double RS, _sols **images) {
/*******************************************   end   *******************************************/
	static double Mag, mag1, mag2, RSi, RSo, delta1,delta2;
	static int NPSsafe;
	Mag = BinaryMag(s, q, y1v, y2v, RS,Tol,images);
	RSi = RS;
	RSo = RS;
	NPSsafe = NPS;
	if (Mag < 0) {
		mag1 = -1;
		delta1 = 3.33333333e-8;
		while (mag1 < 0.1 && RSi>=0) {
			delete *images;
			delta1 *= 3.;
			RSi = RS - delta1;
			mag1 = (RSi > 0) ? BinaryMag(s, q, y1v, y2v, RSi, Tol, images) : BinaryMag0(s,q,y1v,y2v,images);
//			printf("\n-safe1 %lf %lf %d", RSi, mag1, NPS);
			NPSsafe += NPS;
		}
		if(mag1<0) mag1=1.0;
		mag2 = -1;
		delta2 = 3.33333333e-8;
		while (mag2 < 0.1) {
			delta2 *= 3.;
			RSo = RS + delta2;
			delete *images;
			mag2 = BinaryMag(s, q, y1v, y2v, RSo, Tol, images);
//			printf("\n-safe2 %lf %lf %d", RSo,mag2,NPS);
			NPSsafe += NPS;
		}
		Mag = (mag1*delta2 + mag2*delta1)/(delta1+delta2);
	}
	NPS = NPSsafe;

	return Mag;
}

/********************************************************/
/* the overloading is the same as BinaryMag0, see above */ 
/********************************************************/
/******************************************* changed *******************************************/
double VBBinaryLensing::BinaryMag(double a1, double q1, double y1v, double y2v, double RSv, double Tol, _sols_for_skiplist_curve **Images) {
//double VBBinaryLensing::BinaryMag(double a1, double q1, double y1v, double y2v, double RSv, double Tol, _sols **Images) {
/*******************************************   end   *******************************************/
	// checkpoint3 (jump between BinaryMag0 and BinaryMag by finding this string in VScode)
	static complex a, q, m1, m2, y0, y, yc, z, zc;
	static double av = -1.0, qv = -1.0;
	static complex coefs[24], d1, d2, dy, dJ, dz;
	static double thoff = 0.01020304,errbuff;
	static double Mag, th;
////////////////////////////  
	static double errimage, maxerr, currerr, Magold;
	static int NPSmax, flag, NPSold,flagbad,flagbadmax=3;
	/******************************************* changed *******************************************/
	//static _curve *Prov, *Prov2;	// there is also a 'Prov' in NewImages(), but they are different static local variables. 
	static _curve * Prov ;
	static _skiplist_curve * Prov2 ;
	/*******************************************   end   *******************************************/
	static _point *scan1, *scan2;
	static _thetas *Thetas;			// 'Thetas' and 'itheta' are unique in BinaryMag, while others also declared in BinaryMag0
	static _theta *stheta, *itheta;
	/******************************************* changed *******************************************/
	static std::minstd_rand engine_start{std::random_device{}()} ;
	/*******************************************   end   *******************************************/
	
	/******************************************* changed *******************************************/
	static _augmented_priority_queue APQ ;

	if (APQ.apq_array.capacity() > 2048)
	{
		APQ.apq_array.resize(2048) ;
		APQ.sum_tree_array.resize(2048) ;

		APQ.apq_array.shrink_to_fit() ;
		APQ.sum_tree_array.shrink_to_fit() ;
	}
	
	APQ.apq_array.clear() ;
	APQ.sum_tree_array.clear() ;
	/*******************************************   end   *******************************************/

	/******************************************* changed *******************************************/
	int new_and_append_Level_start = 0 ;
	//while ( new_and_append_Level_start < max_skiplist_level && (rand() % 4) == 0 ) 
	while ( new_and_append_Level_start < max_skiplist_level && (engine_start() % 4) == 0 ) 
	{
		new_and_append_Level_start++ ;
	} 
	/*******************************************   end   *******************************************/

	//float time;

#ifdef _PRINT_TIMES
	static double tim0, tim1;
#endif

	// Initialization of the equation coefficients
	
	if ((a1 != av) || (q1 != qv)) {	// completely same as BinaryMag0, except for coefs[23] = RSv; instead of coefs[23] = 0;
		av = a1;
		qv = q1;
		if (q1<1) {
			a = complex(-a1, 0);
			q = complex(q1, 0);
		}
		else {
			a = complex(a1, 0);
			q = complex(1 / q1, 0);
		}
		m1 = 1.0 / (1.0 + q);
		m2 = q*m1;

		coefs[20] = a;
		coefs[21] = m1;
		coefs[22] = m2;
		coefs[6] = a*a;
		coefs[7] = coefs[6] * a;
		coefs[8] = m2*m2;
		coefs[9] = coefs[6] * coefs[8];
		coefs[10] = a*m2;
		coefs[11] = a*m1;
		
	}
	coefs[23] = RSv;

	y0 = complex(y1v, y2v);			// coordinate of source center originating at mass center
	NPS = 1;
	if (Tol>1.) {					// if Tol>1, it means max number of points
		errimage = 0.;
		NPSmax = (int)(Tol);
	}
	else {							// if Tol<=1,      Tol is absolute error tolerance in magnification
		errimage = Tol*M_PI*RSv*RSv;// 			  errimage is absolute error tolerance in image area
		NPSmax =10000; // era 32000
	}
	errbuff = 0;

	// Calculation of the images
	/******************************************* changed *******************************************/
	(*Images) = new _sols_for_skiplist_curve ;
	//(*Images) = new _sols;			// create a _sols class variable without any member (_curve class variables) 
									// and make static pointer 'images'(i.e. *Images) point to that object
	/*******************************************   end   *******************************************/

	Thetas = new _thetas;			// create a _thetas class variable without any member (_theta class variables) 
									// and make static pointer 'Thetas' point to that object
	//std::map<double, _theta *> Thetas_BST ; 	
									// BST is not needed, because as we have to find the interval with largest error and 
									// calculate new 'th', we already know between which two _theta variables to insert

	th = thoff;						// initial 'th' = 0.01020304, slightly shifted from 0
	stheta = Thetas->insert(th);	// create a '_theta' variable, assign input value 'th=0.01020304' to attribute 'th', 
	                        		// and insert it to (currently empty) '_thetas' variable pointed by 'Thetas', 
									// assign the returned newly inserted object's pointer to static local pointer 'stheta'
	/******************************************* changed *******************************************/
	//stheta->maxerr = 1.e100;		// assign 1.e100 to attribute 'maxerr'
	stheta->maxerr  = 0. ;
	stheta->Mag     = 0. ;	
	stheta->astrox1 = 0. ;			
    stheta->astrox2 = 0. ;
	/*******************************************   end   *******************************************/			

	y = y0 + complex(RSv*std::cos(thoff), RSv*std::sin(thoff)); // takes the mass center as origin
	
	/******************************************* changed *******************************************/
	itheta = stheta ;				// let itheta point to the first _theta variable. 
									// because stheta will point to the last _theta variable later, 
									// while we need a pointer to the first _theta variable when inserting the third _theta variable
	/*******************************************   end   *******************************************/

#ifdef _PRINT_TIMES
	tim0 = Environment::TickCount;
#endif
	flag = 0;
	flagbad = 0;
	while (flag == 0) {
		//Prov = NewImages(y, coefs, stheta, time);
		Prov = NewImages(y, coefs, stheta);	// pass a static local variable 'stheta' as actual parameter, 
											// then it's copied from actual parameter to formal parameter 'theta' inside NewImages(), 
											// and can only modify the formal parameter inside NewImages()
											// (but obviously we don't need to modify 'stheta' or 'theta', 
											// we only need to check/modify theta->attribute)
											//
											// return a pointer to a _curve variable which is created on heap during NewImages()
											// then assign the returned pointer to static local variable 'Prov'; 
											// this _curve variable is a linked list of good roots (_point variables) at this zs
											// this _curve variable can contain 0, 3, or 5 elements
											// make all elements (_point variables) point to the same _theta variable

		if (Prov->length > 0) {				// try 'th' = 0.01020304, 0.01020304+0.01, 0.01020304+0.02, ... until the _curve variable contains >0 elements
			flag = 1;
		}
		else {								// Prov->length == 0
											// will happens when: 
											// (if) one (or more) good roots' dJ is too close to 0 (too close to 0 means -1.0e-7 < dJ < 1.0e-7 )
											// or wrongly polishing a positive parity image to a negative parity image, or vice versa, 
											// causing the summation of good root's parity != -1 ; 
											// or 
											// good[worst2] is somehow bad but not bad enough

			delete Prov;					// delete the _curve variable itself, 
											// because another _curve variable will be created in next call of NewImages();
											// no _point variable needs to be freed
			stheta->th += 0.01;
			if (stheta->th > 2.0 * M_PI) {	// if the _curve variable always contains 0 element in ~628 trials, then return -1
											// and no dynamically allocated memory left except for the _sols class variable
				delete Thetas;
				return -1;
			}
			y = y0 + complex(RSv*std::cos(stheta->th), RSv*std::sin(stheta->th));
		}
	}										// now the _curve variable pointed by 'Prov' contains 3 or 5 elements 
											// and stheta->th may = 0.01020304, 0.01020304+0.01, 0.01020304+0.02, ...
	
	/******************************************* changed *******************************************/
	APQ.push_augmented_heap(0., itheta) ;	// this element will be first popped out anyway, 
											// so we can just set maxerr_to_push to 0.
	/*******************************************   end   *******************************************/

#ifdef _PRINT_TIMES
	tim1 = Environment::TickCount;
	GM += tim1 - tim0;
#endif
	stheta = Thetas->insert(2.0*M_PI + Thetas->first->th);
	stheta->maxerr = 0.;
	stheta->Mag = 0.;
    stheta->astrox1 = 0.;
    stheta->astrox2 = 0.;
	stheta->errworst = Thetas->first->errworst;
	// the _thetas variable now contains 2 elements, looks like: 
	// 			   (itheta points to below element)						(stheta points to below element)
	// th          0.01020304  										  	0.01020304 + 2*M_PI
	// maxerr      0.	 											  	0.
	// Mag         0.													0.
	// astrox1,2   (0.,0.)				                              	(0.,0.)
	// errworst    distance between two ghost roots (3 good roots)    	distance between two ghost roots (3 good roots)
	//			   -1.e100   						(5 good roots)   	-1.e100   					  	 (5 good roots)

	for (scan1 = Prov->first; scan1; scan1 = scan2) {	// loop through all 3 or 5 good roots in _curve variable
		scan2 = scan1->next;
		/******************************************* changed *******************************************/
		//Prov2 = new _curve(scan1);
		Prov2 = new _skiplist_curve(scan1, new_and_append_Level_start) ;			// create an object of class _curve with one member(_point class variable),
												 		// input is pointer(scan1) that points to the member; 
														// pointer to that object is assigned to static local variable 'Prov2'

		Prov2->append(scan1->x1, scan1->x2, new_and_append_Level_start) ;			// create a new _point variable on heap, 
														// and then append it to the end of the _curve variable pointed by 'Prov2'
														// 
														// constructor of _point class is called, assign input values to attributes x1, x2, 
														// and set attribute theta(pointer) = NULL
		/*******************************************   end   *******************************************/

		Prov2->last->theta = stheta;
		Prov2->last->d = Prov2->first->d;
		Prov2->last->dJ = Prov2->first->dJ;
		Prov2->last->ds = Prov2->first->ds;

		(*Images)->append(Prov2);						// append the newly created _curve variable (pointed by 'Prov2') 
												  		// to the end of sols (a linked list, object of class _sols)
														// which is pointed by *Images (i.e. static local pointer 'images' which is out of the scope of current function)
	}													
	// now we get a _sols variable(pointed by *Images), containing 3 or 5 _curve variables, 
	// and each _curve variable contains 2 _point variables, looks like: 
	// 
	// (x1,x2) 		(x1,x2) at theta=0.01020304 (origin at mass center)	(x1,x2) at theta=0.01020304 (origin at mass center)
	// theta   		points to _theta variable with th=0.01020304  		points to _theta variable with th=0.01020304 + 2*M_PI
	// d			z'(theta) at theta=0.01020304						z'(theta) at theta=0.01020304
	// dJ			Jacobian determinant at theta=0.01020304			Jacobian determinant at theta=0.01020304
	// ds			x'(theta)^x"(theta) at theta=0.01020304				x'(theta)^x"(theta) at theta=0.01020304
	// J2 			Partial^2(zs_c)/Partial(z)^2 at theta=0.01020304	in-determinate state

	Prov->length = 0;
	delete Prov;

	// checkpoint2 (jump between BinaryMag and OrderImages by finding this string in VScode)

	th = M_PI + Thetas->first->th;			// th = 0.01020304 + M_PI
	flag = 0;
	Magold = -1.;
	NPSold = 2;
	

	/******************************************* changed *******************************************/
	Mag = 0. ;								// initialize total Mag to 0. outside the do-loop
											// (notice Mag is a static local variable, so should be initialized every time calls BinaryMag)
											// 
											// before inserting the third element, 
											// principally we should calculate the Mag contributed by the only existing interval as current total Mag, 
											// and the same value will be first reduced in OrderImages()
											// but as the first interval's Mag is set to 0. previously, 
											// we can initialize Mag to 0. then 0. - 0. = 0.
	currerr = 1.e100;
	//currerr = 0. ;
	astrox1 = 0. ;							// declaration is in .h (in class definition), public
	astrox2 = 0. ;
	/*******************************************   end   *******************************************/

	do {
		//stheta = Thetas->insert(th);    	// create a '_theta' variable, assign input value 'th=0.01020304 + M_PI' to attribute 'th', 
											// and insert it to (currently has 2 elements) '_thetas' variable pointed by 'Thetas', 
											// assign the returned newly inserted object's pointer to static local pointer 'stheta'
											/****************************************************************************************/
											// this insert should be optimized. 
											// it costs 3.5 s, while complex_laguerre2newton only costs 4 s
											// 
											// we still need a sorted data structure to store _thetas variable, 
											// because we need to know stheta->prev and next; 
											// 
											// Fortunately, BST is not needed, 
											// because as we have found the interval with largest error and calculated new 'th' before inserting, 
											// we already know between which two _theta variables to insert
											/****************************************************************************************/
		/******************************************* changed *******************************************/
		stheta = Thetas->insert_at_certain_position(itheta, th) ; 
											// this method can only be used when inserting an element in the middle of linked list
											// i.e. *first's 'th' < current 'th' < *last's 'th'
											// but we can safely use this method as we only insert in the middle of linked list inside the do-loop
											// 
											// and the new element is forced to be inserted between itheta and itheta->next, 
											// which means it's the programmer's responsibility to guarantee itheta->th < th < itheta->next->th holds
											//
											// however, we already know new element should be inserted in interval [itheta, itheta->next] and
											// th is calculated to be between itheta->th and itheta->next->th
		/*******************************************   end   *******************************************/

		y = y0 + complex(RSv*std::cos(th), RSv*std::sin(th));
#ifdef _PRINT_TIMES
		tim0 = Environment::TickCount;
#endif
		//if (NPS == 422) {
		//	NPS = NPS;
		//}

		//Prov = NewImages(y, coefs, stheta, time);
		Prov = NewImages(y, coefs, stheta);	// pass a static local variable 'stheta' as actual parameter, 
											// then it's copied from actual parameter to formal parameter 'theta' inside NewImages(), 
											// and can only modify the formal parameter inside NewImages()
											// (but obviously we don't need to modify 'stheta' or 'theta', 
											// we only need to check/modify theta->attribute)
											//
											// return a pointer to a _curve variable which is created on heap during NewImages()
											// then assign the returned pointer to static local variable 'Prov'; 
											// this _curve variable is a linked list of good roots (_point variables) at this zs
											// this _curve variable can contain 0, 3, or 5 elements
											// make all elements (_point variables) point to the same _theta variable

#ifdef _PRINT_TIMES
		tim1 = Environment::TickCount;
		GM += tim1 - tim0;
#endif
		if (Prov->length > 0) {
			flagbad = 0;
			
			/******************************************* changed *******************************************/
			Mag -= stheta->prev->Mag ; 
			//currerr -= theta->prev->maxerr ; 

			if (astrometry){
				astrox1 -= stheta->prev->astrox1 ; 
				astrox2 -= stheta->prev->astrox2 ; 
			}
			/*******************************************   end   *******************************************/

			OrderImages((*Images), Prov); 	// at the beginning, _sols variable pointed by *Images can have 3 _curve variables
											// while _curve variable pointed by 'Prov' have 5 _point variables, or vice versa
											//
											// _curve variable pointed by 'Prov' is a linked list of good roots (_point variables) at this zs
											// can contain 3 or 5 elements
											// all elements (_point variables) point to the same _theta variable (pointed by 'stheta')
											//
											// _curve variable pointed by 'Prov' is deleted inside OrderImages()

			/******************************************* changed *******************************************/
			Mag += stheta->prev->Mag ;
			Mag += stheta->Mag ; 

			//currerr += theta->prev->maxerr ;
			//currerr += theta->maxerr ; 

			if (astrometry){
				astrox1 += stheta->prev->astrox1 ; 
				astrox1 += stheta->astrox1 ; 

				astrox2 += stheta->prev->astrox2 ; 
				astrox2 += stheta->astrox2 ; 
			}
			/*******************************************   end   *******************************************/

			if ((stheta->th - stheta->prev->th)*RSv < 1.e-11/* || stheta->maxerr > jumperrfactor * currerr || stheta->prev->maxerr > jumperrfactor * currerr*/) {
				
				/******************************************* changed *******************************************/
				//currerr -= (stheta->maxerr + stheta->prev->maxerr) ; 
											// reduce the two 'maxerr's from currerr, which are newly added into currerr in OrderImages()
				/*******************************************   end   *******************************************/
				errbuff += stheta->maxerr + stheta->prev->maxerr;
				stheta->maxerr = 0;
				stheta->prev->maxerr = 0;				// stop to insert new theta behind stheta and stheta->prev
			}

			/******************************************* changed *******************************************/
			APQ.pop_then_push_augmented_heap(stheta->prev->maxerr, stheta->prev) ;
			APQ.push_augmented_heap(stheta->maxerr, stheta) ;
			/*******************************************   end   *******************************************/

		} else {							// Prov->length == 0
											// will happens when: 
											// (if) one (or more) good roots' dJ is too close to 0 (too close to 0 means -1.0e-7 < dJ < 1.0e-7 )
											// or wrongly polishing a positive parity image to a negative parity image, or vice versa, 
											// causing the summation of good root's parity != -1 ; 
											// or 
											// good[worst2] is somehow bad but not bad enough

			delete Prov;					// delete the _curve variable itself, 
											// because another _curve variable will be created in next call of NewImages();
											// no _point variable needs to be freed

			flagbad++; 						// flagbad is initialized to 0 around 150 lines ahead

			if (flagbad == flagbadmax) {	// flagbad == 3
				if (NPS < 16) {
					delete Thetas;
					return -1;				// no dynamically allocated memory left except for the _sols class variable
				}
				/******************************************* changed *******************************************/
				//currerr -= stheta->prev->maxerr ; 
											// previously this interval is identified to have the largest 'maxerr', 
											// but three 'th' all lead to Prov->length == 0, 
											// then give up to insert new theta bewteen stheta->prev and stheta->next; 
											// 
											// so, reduce this largest 'maxerr' from currerr
				/*******************************************   end   *******************************************/
				errbuff += stheta->prev->maxerr;
				stheta->prev->maxerr = 0;	// give up to insert new theta behind stheta->prev
				
				/******************************************* changed *******************************************/
				APQ.pop_then_push_augmented_heap(0., stheta->prev) ;
				/*******************************************   end   *******************************************/
				
				NPS--;
				NPSmax--;
			}
			else {							// flagbad == 1 or 2
											// 
											// if current stheta->th is closer to stheta->next->th, then th = (th + prev->th)/2 or th = (th + 2*prev->th)/3
											//							  		  stheta->prev->th       th = (th + next->th)/2 or th = (th + 2*next->th)/3
				th = (th - stheta->prev->th >= stheta->next->th - th) ? (th + flagbad * stheta->prev->th) / (1 + flagbad) : (th + flagbad * stheta->next->th) / (1 + flagbad);
			}
			Thetas->remove(stheta);			// remove and delete the element pointed by the input pointer, 
	                        				// and update that element's prev/next element(if they exist)'s 'next'/'prev'
											// (loop over the linked list, thus O(N) complexity)
		}	
		// initially, th = (stheta->prev->th + stheta->next->th)/2 after insertion, 
		// if this 'th' leads to Prov->length == 0, then flagbad=1, th = (th + prev->th)/2 = (3*prev->th + next->th)/4, 
		// remove and delete the old _theta variable, and insert a _theta variable with new 'th' in the next iteration of do-loop
		// (the stheta are different while the stheta->prev/next are same elements for the old and new _theta variables), 
		// and call NewImages() again;
		//  
		// if the new 'th' again leads to Prov->length == 0, then flagbad=2, th = (th + 2*next->th)/3 = (prev->th + 3*next->th)/4, 
		// and same as above; 
		// 
		// if the new 'th' leads to Prov->length == 0 for the third time, then flagbad=3, 
		// if we only have 1-15 points already inserted except for the current one, then return -1; 
		// or 
		// we give up to insert new theta bewteen stheta->prev and stheta->next, 
		// by setting stheta->prev->maxerr from the largest one to 0 and picking out the original second largest maxerr
		// (but the original stheta->prev->maxerr is accumulated into errbuff)
		// 
		// To summarize: 
		// if the bisection 'th' leads to Prov->length == 0, try the two quadrisection points successively, 
		// if the three 'th' all lead to Prov->length == 0, then give up to insert new theta bewteen stheta->prev and stheta->next
		// (but the original stheta->prev->maxerr is accumulated into errbuff)
		

		if (flagbad == 0 || flagbad == flagbadmax) {
			/******************************************* changed *******************************************/
			 flagbad = 0 ;						// should be a flagbad = 0 ; here? 
												// Otherwise, after flagbad==3, select the original second largest error interval, 
												// and if this interval also leads to Prov->length == 0, then flagbad==4, 
												// thus goes into an unstoppable loop until Prov->length > 0
			//maxerr = 0. ; //currerr = Mag = 0.;
			//astrox1 = astrox2 = 0.;
			
			itheta  = APQ.apq_array[0].stheta ;
			currerr = APQ.sum_tree_array[0].sumerr ;
			/*******************************************   end   *******************************************/

			/******************************************* changed *******************************************/
//			stheta = Thetas->first;
//
//			while (stheta->next) {              /****************************************************************************************/
//												// Don't loop through all elements(except 'last') of Thetas(_thetas variable) when inserting a new _theta variable
//												//
//												// this while loop costs almost all 10 s execution time of BinaryMag in high magnification region, 
//												// where the total execution time is 80 s and complex_laguerre2newton only costs 4 s. 
//												// (in triple code, this problem is alleviated(complex_laguerre2newton = MultiMag = 40 s), 
//												// because complex_laguerre2newton is more time-consuming in triple?)
//												//
//												// currerr and Mag should be calculated through an incremental way
//												// do {
//												//    	currerr = currerr - itheta->maxerr ;
//												//		stheta = Thetas->insert(th);
//												//      NewImages(); OrderImages(); // currently don't understand their behavior, so may be wrong
//												//            						// maybe several _theta variables' maxerr are changed in OrderImages()
//												//      currerr = currerr + stheta->maxerr + stheta->prev->maxerr ; 
//												//      find the new 'itheta' with largest maxerr ;
//												//		th = (itheta->th + itheta->next->th) / 2 ;
//												// }while()
//												// 
//												// we should store the pointer to the _theta variable which has the largest and 2nd-largest maxerr outside the do-loop,
//												// after inserting a new _theta variable in the middle of the largest-maxerr-section, 
//												// only compare the new-got two maxerr and the original 2nd-largest maxerr
//												//
//												// Correction: the above won't work! 
//												// if the new-got two maxerr are both very small, 
//												// the original 2nd-largest maxerr becomes the NEW largest maxerr, and in which you insert a new _theta variable, 
//												// but now you don't know what is the NEW 2nd-largest maxerr, so same procedure cannot continue
//												//
//												// so we need to seek for a new data structure, like Binary Search Tree, 
//												// which can enable finding the largest maxerr and inserting new-got two maxerr
//												// (there will be branch and DRAM access, 
//												// try to use prefetch (prefetch both braches before the if(to mitigate latency)
//												/****************************************************************************************/
//				// currerr += stheta->maxerr;
//				// Mag += stheta->Mag;
//
//				//if (astrometry) {
//				//	astrox1 += stheta->astrox1;
//				//	astrox2 += stheta->astrox2;
//				//}
//
// #ifndef _uniform
// 				if (stheta->maxerr > maxerr) {
// 					maxerr = stheta->maxerr;
// #else
// 				if (stheta->next->th * 0.99999 - stheta->th > maxerr) {
// 					maxerr = stheta->next->th - stheta->th;
// #endif
// 					itheta = stheta;
// 				}
// 				stheta = stheta->next;
// #ifdef _selectimage
// 				if ((NPS == NPSmax - 1) && (fabs(floor(stheta->th / M_PI * _npoints / 2 + 0.5) - stheta->th / M_PI * _npoints / 2) < 1.e-8)) {
// 					printf("%d %.15le\n", (int)floor(stheta->th / M_PI * _npoints / 2 + 0.5), Mag);
// 				}
// #endif
// 			}	// end of while-loop
			/*******************************************   end   *******************************************/

			th = (itheta->th + itheta->next->th) / 2.0 ; 
			NPS++;
//#ifndef _uniform
			if (fabs(Magold - Mag) * 2 < errimage) {
												// errimage is absolute error tolerance in image area
				flag++;							// flag was initialized to 0 outside the do-loop and never used inside the do-loop
			}
			else {
				flag = 0;
				Magold = Mag;
				NPSold = NPS + 8;				// NPSold was initialized to 2 outside the do-loop
			}
// #else
// 			currerr = 2 * errimage;
// 			if (NPS == 2 * NPSold) {
// 				if (fabs(Magold - Mag) * 2 < errimage) {
// 					flag = NPSold;
// 				}
// 				else {
// 					flag = 0;
// 					NPSold = NPS;
// 					Magold = Mag;
// 				}
// 			}
// #endif
#ifdef _PRINT_ERRORS2
			printf("\nNPS= %d Mag = %lf maxerr= %lg currerr =%lg th = %lf", NPS, Mag / (M_PI * RSv * RSv), maxerr / (M_PI * RSv * RSv), currerr / (M_PI * RSv * RSv), th);
#endif
		}
	} while ((currerr > errimage) && (currerr > RelTol * Mag) && (NPS < NPSmax) && ((flag < NPSold)/* || NPS<8 ||(currerr>10*errimage)*/)/*&&(flagits)*/);
    
	if(astrometry){
		astrox1 /= (Mag);
		astrox2 /= (Mag);
    }
	Mag /= (M_PI*RSv*RSv);
	therr = (currerr+errbuff) / (M_PI*RSv*RSv);
 
	delete Thetas;
//	if (NPS == NPSmax) return 1.e100*Tol; // Only for testing

	//return NPS;
	return Mag;
       
}

/********************************************************/
/* the overloading is the same as BinaryMag0, see above */ 
/********************************************************/
double VBBinaryLensing::BinaryMag(double a1, double q1, double y1v, double y2v, double RSv, double Tol) {
	/******************************************* changed *******************************************/
	//static _sols *images;
	static _sols_for_skiplist_curve * images ;
	/*******************************************   end   *******************************************/
	static double mag;
	mag = BinaryMag(a1, q1, y1v, y2v, RSv, Tol, &images);
	delete images;
	return mag;
}


double VBBinaryLensing::BinaryMag2(double s, double q, double y1v, double y2v, double rho) {
	static double Mag, rho2, y2a;//, sms , dy1, dy2;
	//static int c;
	/******************************************* changed *******************************************/
	//static _sols *Images;
	static _sols_for_skiplist_curve *Images ;
	/*******************************************   end   *******************************************/

	//c = 0;


	y2a = fabs(y2v);

	Mag0 = BinaryMag0(s, q, y1v, y2a, &Images);
	delete Images;
	rho2 = rho*rho;
	corrquad *= 6 * (rho2 + 1.e-4*Tol);
	corrquad2 *= (rho+1.e-3);
	if (corrquad<Tol && corrquad2<1 && (/*rho2 * s * s<q || */ safedist>4 * rho2)) {
		Mag = Mag0;
	}
	else {
		Mag = BinaryMagDark(s, q, y1v, y2a, rho, Tol);
	}
	Mag0 = 0;

	if (y2v < 0) {
		y_2 = y2v;
		astrox2 = -astrox2;
	}
	return Mag;
}


double VBBinaryLensing::BinaryMagDark(double a, double q, double y1, double y2, double RSv, double Tolnew) {
	static double Mag, Magold, Tolv;
    static double LDastrox1,LDastrox2;
	static double tc, lc, rc, cb,rb;
	static int c, flag;
	static double currerr, maxerr;
	static annulus *first, *scan, *scan2;
	static int nannold, totNPS;
	/******************************************* changed *******************************************/
	//static _sols *Images;
	static _sols_for_skiplist_curve *Images;
	/*******************************************   end   *******************************************/

	Mag = -1.0;
	Magold = 0.;
	Tolv = Tol;
	LDastrox1 = LDastrox2 = 0.0;
	c = 0;
	totNPS = 1;

	Tol = Tolnew;
	y_1 = y1;
	y_2 = y2;
	while ((Mag<0.9) && (c<3)) {

		first = new annulus;
		first->bin = 0.;
		first->cum = 0.;
		if (Mag0 > 0.5) {
			first->Mag = Mag0;
			first->nim = nim0;
		}
		else {
			first->Mag = BinaryMag0(a, q, y_1, y_2, &Images);
			first->nim = Images->length;
			delete Images;
		}
		if (astrometry) {
			first->LDastrox1 = astrox1 * first->Mag;
			first->LDastrox2 = astrox2 * first->Mag;
		}
		scr2 = sscr2 = 0;
		first->f = LDprofile(0);
		first->err = 0;
		first->prev = 0;


		first->next = new annulus;
		scan = first->next;
		scan->prev = first;
		scan->next = 0;
		scan->bin = 1.;
		scan->cum = 1.;
		scan->Mag = BinaryMagSafe(a, q, y_1, y_2, RSv, &Images);
		if(astrometry){
			scan->LDastrox1 = astrox1*scan->Mag;
			scan->LDastrox2 = astrox2*scan->Mag;
		}
		totNPS += NPS;
		scan->nim = Images->length;
		delete Images;
		scr2 = sscr2 = 1;
		scan->f = LDprofile(0.9999999);
		if (scan->nim == scan->prev->nim) {
			scan->err = fabs((scan->Mag - scan->prev->Mag)*(scan->prev->f - scan->f) / 4);
		}
		else {
			scan->err = fabs((scan->Mag)*(scan->prev->f - scan->f) / 4);
		}

		Magold = Mag = scan->Mag;
		if(astrometry){
			LDastrox1=scan->LDastrox1;
			LDastrox2=scan->LDastrox2;
		}
		//			scan->err+=scan->Mag*Tolv*0.25; //Impose calculation of intermediate annulus at mag>4. Why?
		currerr = scan->err;
		flag = 0;
		nannuli = nannold = 1;
		while (((flag<nannold + 5) && (currerr>Tolv) && (currerr>RelTol*Mag)) || (nannuli<minannuli)) {
			maxerr = 0;
			for (scan2 = first->next; scan2; scan2 = scan2->next) {
#ifdef _PRINT_ERRORS_DARK
				printf("\n%d %lf %le | %lf %le", nannuli, scan2->Mag, scan2->err, Mag, currerr);
#endif
				if (scan2->err>maxerr) {
					maxerr = scan2->err;
					scan = scan2;
				}
			}

			nannuli++;
			Magold = Mag;
			Mag -= (scan->Mag*scan->bin*scan->bin - scan->prev->Mag*scan->prev->bin*scan->prev->bin)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
            if(astrometry){
				LDastrox1 -= (scan->LDastrox1*scan->bin*scan->bin - scan->prev->LDastrox1*scan->prev->bin*scan->prev->bin)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
				LDastrox2 -= (scan->LDastrox2*scan->bin*scan->bin - scan->prev->LDastrox2*scan->prev->bin*scan->prev->bin)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
			}
			currerr -= scan->err;
			lc = scan->prev->cum;
			rc = scan->cum;
			tc = (lc + rc) *0.5;
			cb = rCLDprofile(tc,scan->prev,scan);
			scan->prev->next = new annulus;
			scan->prev->next->prev = scan->prev;
			scan->prev = scan->prev->next;
			scan->prev->next = scan;
			scan->prev->bin = cb;
			scan->prev->cum = tc;
			scan->prev->f = LDprofile(cb);
			scan->prev->Mag = BinaryMagSafe(a, q, y_1, y_2, RSv*cb, &Images);
			if(astrometry){
				scan->prev->LDastrox1=astrox1*scan->prev->Mag;
				scan->prev->LDastrox2=astrox2*scan->prev->Mag;
			}
			totNPS += NPS;
			scan->prev->nim = Images->length;
			if (scan->prev->prev->nim == scan->prev->nim) {
				scan->prev->err = fabs((scan->prev->Mag - scan->prev->prev->Mag)*(scan->prev->prev->f - scan->prev->f)*(scan->prev->bin*scan->prev->bin - scan->prev->prev->bin*scan->prev->prev->bin) / 4);
			}
			else {
				scan->prev->err = fabs((scan->prev->bin*scan->prev->bin*scan->prev->Mag - scan->prev->prev->bin*scan->prev->prev->bin*scan->prev->prev->Mag)*(scan->prev->prev->f - scan->prev->f) / 4);
			}
			if (scan->nim == scan->prev->nim) {
				scan->err = fabs((scan->Mag - scan->prev->Mag)*(scan->prev->f - scan->f)*(scan->bin*scan->bin - scan->prev->bin*scan->prev->bin) / 4);
			}
			else {
				scan->err = fabs((scan->bin*scan->bin*scan->Mag - scan->prev->bin*scan->prev->bin*scan->prev->Mag)*(scan->prev->f - scan->f) / 4);
			}
			rb = (scan->Mag + scan->prev->prev->Mag - 2 * scan->prev->Mag);
			scan->prev->err += fabs(rb*(scan->prev->prev->f - scan->prev->f)*(scan->prev->bin*scan->prev->bin - scan->prev->prev->bin*scan->prev->prev->bin));
			scan->err += fabs(rb*(scan->prev->f - scan->f)*(scan->bin*scan->bin - scan->prev->bin*scan->prev->bin));
#ifdef _PRINT_ERRORS_DARK
			printf("\n%d", Images->length);
#endif
			delete Images;

			Mag += (scan->bin*scan->bin*scan->Mag - cb*cb*scan->prev->Mag)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
			Mag += (cb*cb*scan->prev->Mag - scan->prev->prev->bin*scan->prev->prev->bin*scan->prev->prev->Mag)*(scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin*scan->prev->bin - scan->prev->prev->bin*scan->prev->prev->bin);
			currerr += scan->err + scan->prev->err;
            if(astrometry){
				LDastrox1 += ( scan->bin*scan->bin*scan->LDastrox1 - cb*cb*scan->prev->LDastrox1)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
				LDastrox1 += ( cb*cb*scan->prev->LDastrox1 - scan->prev->prev->bin*scan->prev->prev->bin*scan->prev->prev->LDastrox1)*(scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin*scan->prev->bin - scan->prev->prev->bin*scan->prev->prev->bin);
				LDastrox2 += ( scan->bin*scan->bin*scan->LDastrox2 - cb*cb*scan->prev->LDastrox2)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
				LDastrox2 += ( cb*cb*scan->prev->LDastrox2 - scan->prev->prev->bin*scan->prev->prev->bin*scan->prev->prev->LDastrox2)*(scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin*scan->prev->bin - scan->prev->prev->bin*scan->prev->prev->bin);
			}


			if (fabs(Magold - Mag) * 2<Tolv) {
				flag++;
			}
			else {
				flag = 0;
				nannold = nannuli;
			}

		}

		if (multidark) {
			annlist = first;
		}else{
			while (first) {
				scan = first->next;
				delete first;
				first = scan;
			}
		}

		Tolv /= 10;
		c++;
	}
	NPS = totNPS;
	therr = currerr;
    if(astrometry){
		LDastrox1/=Mag;
		LDastrox2/=Mag;
		astrox1=LDastrox1;
		astrox2=LDastrox2;
    }
	return Mag;
}

void VBBinaryLensing::BinaryMagMultiDark(double a, double q, double y1, double y2, double RSv, double *a1_list, int nfil, double *mag_list, double Tol) {
	annulus *scan;
	int imax = 0;
	double Mag,r2,cr2,scr2,a1;

	multidark = true;

	for (int i = 1; i < nfil; i++) {
		if (a1_list[i] > a1_list[imax]) imax = i;
	}
	a1 = a1_list[imax];
	mag_list[imax] = BinaryMagDark(a, q, y1, y2, RSv, Tol);

	for (int i = 0; i < nfil; i++) {
		if (i != imax) {
			Mag = 0;
			a1 = a1_list[i];
			for (scan = annlist->next; scan; scan = scan->next) {
				r2 = scan->bin*scan->bin;
				cr2 = 1 - r2;
				scr2 = sqrt(cr2);
				scan->cum = (3 * r2*(1 - a1) - 2 * a1*(scr2*cr2 - 1)) / (3 - a1);
				Mag += (scan->bin*scan->bin*scan->Mag - scan->prev->bin*scan->prev->bin*scan->prev->Mag)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
			}
			mag_list[i] = Mag;
		}
	}

	while (annlist) {
		scan = annlist->next;
		delete annlist;
		annlist = scan;
	}

	multidark = false;
}

double VBBinaryLensing::LDprofile(double r) {
	static int ir;
	static double rr,ret;
	switch(curLDprofile){
	case LDuser:
		rr = r * npLD;
		ir = (int)rr;
		rr -= ir;
		ret= LDtab[ir] * (1 - rr) + LDtab[ir + 1] * rr;
		break;
	case LDlinear:
		ret = 3 / (3 - a1)*(1 - a1 * scr2);
		break;
	case LDsquareroot:
		ret = 3 / (3 - a1 - 0.6*a2)*(1 - a1 * scr2 - a2 * sscr2);
        break;
	case LDquadratic:
		ret = 3 / (3 - a1 - 0.5*a2)*(1 - a1 * scr2 - a2 * sscr2);
		break;
	case LDlog:
		ret = 3 / (3 - a1 + 0.666666666666 * a2)*(1 - a1 * scr2 - a2 * sscr2);
		break;
	}
	return ret;
}

double VBBinaryLensing::rCLDprofile(double tc,annulus *left,annulus *right) {
	static int ic;
	static double rc,cb,lc,r2,cr2,cc,lb,rb;

	switch (curLDprofile) {
	case LDuser:
		rc = tc * npLD;
		ic = (int)rc;
		rc -= ic;
		cb = rCLDtab[ic] * (1 - rc) + rCLDtab[ic + 1] * rc;
		break;
	case LDlinear:
		lb = left->bin;
		rb = right->bin;
		lc = left->cum;
		rc = right->cum;
		do {
			cb = rb + (tc - rc)*(rb - lb) / (rc - lc);
			r2 = cb * cb; 
			cr2 = 1 - r2;
			scr2 = 1-sqrt(cr2);
			cc = (3 * r2 - a1 * (r2 - 2 * scr2*cr2)) / (3 - a1);
			if (cc > tc) {
				rb = cb;
				rc = cc;
			}
			else {
				lb = cb;
				lc = cc;
			}
		} while (fabs(cc - tc) > 1.e-5);
		break;
	case LDsquareroot:
		lb = left->bin;
		rb = right->bin;
		lc = left->cum;
		rc = right->cum;
		do {
			cb = rb + (tc - rc)*(rb - lb) / (rc - lc);
			r2 = cb * cb;
			cr2 = 1 - r2;
			scr2 = sqrt(cr2);
			sscr2 = 1 - sqrt(scr2);
			scr2 = 1 - scr2;
			cc = (3 * r2 - a1 * (r2 - 2 * scr2*cr2) - 0.6*a2*(r2 - 4 * sscr2*cr2)) / (3 - a1 - 0.6*a2);
			if (cc > tc) {
				rb = cb;
				rc = cc;
			}
			else {
				lb = cb;
				lc = cc;
			}
		} while (fabs(cc - tc) > 1.e-5);
		break;
	case LDquadratic:
		lb = left->bin;
		rb = right->bin;
		lc = left->cum;
		rc = right->cum;
		do {
			cb = rb + (tc - rc)*(rb - lb) / (rc - lc);
			r2 = cb * cb;
			cr2 = 1 - r2;
			scr2 = 1- sqrt(cr2);
			sscr2 = scr2*scr2;
			cc = (3 * r2 - a1 * (r2 - 2 * scr2*cr2) + a2*(4*scr2-(2+4*scr2)*r2+1.5*r2*r2)) / (3 - a1 - 0.5*a2);
			if (cc > tc) {
				rb = cb;
				rc = cc;
			}
			else {
				lb = cb;
				lc = cc;
			}
		} while (fabs(cc - tc) > 1.e-5);
		break;
	case LDlog:
		lb = left->bin;
		rb = right->bin;
		lc = left->cum;
		rc = right->cum;
		do {
			cb = rb + (tc - rc)*(rb - lb) / (rc - lc);
			r2 = cb * cb;
			cr2 = 1 - r2;
			scr2 = sqrt(cr2);
			sscr2 = scr2*log(scr2);
			scr2 = 1 - scr2;
			cc = (3 * r2 - a1 * (r2 - 2 * scr2*cr2) + 2*a2*(scr2*(1+scr2*(scr2/3-1)) + sscr2*cr2)) / (3 - a1 + 0.6666666666666666*a2);
			if (cc > tc) {
				rb = cb;
				rc = cc;
			}
			else {
				lb = cb;
				lc = cc;
			}
		} while (fabs(cc - tc) > 1.e-5);
		break;
	}

	return cb;
}

void VBBinaryLensing::SetLDprofile(double (*UserLDprofile)(double),int newnpLD) {
	int ic,ir;
	if (npLD > 0) {
		free(LDtab);
		free(rCLDtab);
	}
	if (newnpLD > 0) {
		npLD = newnpLD;
		double npLD2 = npLD * npLD;
		LDtab = (double *)malloc(sizeof(double)*(npLD+1));
		CLDtab = (double *)malloc(sizeof(double)*(npLD + 1));
		rCLDtab = (double *)malloc(sizeof(double)*(npLD+1));

		LDtab[0] = UserLDprofile(0.);
		CLDtab[0] = 0.;
		for (int i = 1; i <= npLD; i++) {
			LDtab[i] = UserLDprofile(((double) i)/ npLD);
			CLDtab[i] = CLDtab[i-1]+(LDtab[i]*i + LDtab[i - 1]*(i-1));
		}
		for (int i = 0; i <= npLD; i++) {
			LDtab[i] *= npLD2/CLDtab[npLD];
			CLDtab[i] /= CLDtab[npLD];
		}
		ic = 1;
		rCLDtab[0] = 0;
		ir = 1;
		while (ic < npLD) {
			while (CLDtab[ir] * npLD < ic && ir<npLD) ir++;
			rCLDtab[ic] = ((CLDtab[ir] - ((double) ic) / npLD)*(ir-1) + (((double)ic) / npLD - CLDtab[ir-1])*ir) / (CLDtab[ir] - CLDtab[ir - 1])/npLD;
			ic++;
		}
		rCLDtab[npLD] = 1;


		//printf("\n\n--------------------");
		//annulus left, right;
		//left.cum = left.bin=0;
		//right.cum = right.bin=1;
		//for (int i = 0; i <= npLD; i++) {
		//	double rl, fl;
		//	rl = ((double)i) / npLD;
		//	scr2 = 1 - sqrt(1 - rl * rl);
		//	fl = LDprofile(rl);
		//	rl = rCLDprofile(((double) i) / npLD,&left,&right);
		//	printf("\n%lf %lf %lf %lf", fl, LDtab[i], rl, rCLDtab[i]);
		//}
		//printf("\n--------------------\n\n");

		free(CLDtab);
		curLDprofile = LDuser;
	}
	else {
		npLD = 0;
		curLDprofile = LDlinear;
	}
}

void VBBinaryLensing::SetLDprofile(LDprofiles LDval) {
	if(npLD > 0) {
		npLD = 0;
		free(LDtab);
		free(rCLDtab);
	}
	curLDprofile = LDval;
}

void VBBinaryLensing::LoadESPLTable(char *filename){
	FILE *f;

	if((f = fopen(filename, "rb"))!=0){
		// ******************************************* changed *******************************************
		fread(ESPLin, sizeof(double), __rsize_ESPL * __zsize_ESPL, f);
		fread(ESPLout, sizeof(double), __rsize_ESPL * __zsize_ESPL, f);
        fread(ESPLinastro, sizeof(double), __rsize_ESPL * __zsize_ESPL, f);
		fread(ESPLoutastro, sizeof(double), __rsize_ESPL * __zsize_ESPL, f);
		// *******************************************   end   *******************************************
		fclose(f);
		ESPLoff=false;
	}else{
		printf("\nESPL table not found !");
	}
}


double VBBinaryLensing::PSPLMag(double u) {
	static double u2,u22;
	u2 = u * u;
	u22 = u2 + 2;
	if (astrometry) {
		astrox1 = u+u/u22;
	}
	return  u22 / sqrt(u2 * (u2 + 4));
}


double VBBinaryLensing::ESPLMag(double u, double RSv) {
	double mag,z,fr,cz,cr,u2;
	int iz, ir;
       
	if (ESPLoff) {
		printf("\nLoad ESPL table first!");
		return 0;
	}
         
	fr = -10.857362047581296* log(0.01* RSv);
	// ******************************************* changed *******************************************
	if (fr > __rsize_ESPL - 1) fr = __rsize_ESPL -1.000001;
	// *******************************************   end   *******************************************
	if (fr < 0) printf("Source too large!");
	ir = (int) floor(fr);
	fr -= ir;
	cr = 1 - fr;

	z = u / RSv;

	if (z < 1) {
		// ******************************************* changed *******************************************
		z *= __zsize_ESPL -1;
		// *******************************************   end   *******************************************
		iz = (int) floor(z);
		z -= iz;
		cz = 1 - z;
		mag = sqrt(1 + 4. / (RSv*RSv));
		mag *= ESPLin[ir][iz] * cr*cz + ESPLin[ir + 1][iz] * fr*cz + ESPLin[ir][iz + 1] * cr*z + ESPLin[ir + 1][iz + 1] * fr*z;
                if (astrometry) {
                	astrox1=(1-1./(4+RSv*RSv))*u;
                	astrox1 *= ESPLinastro[ir][iz] * cr*cz + ESPLinastro[ir + 1][iz] * fr*cz + ESPLinastro[ir][iz + 1] * cr*z + ESPLinastro[ir + 1][iz + 1] * fr*z;
                }
	}
	else {
		z = 0.99999999999999 / z;
		// ******************************************* changed *******************************************
		z *= __zsize_ESPL - 1;
		// *******************************************   end   *******************************************
		iz = (int)floor(z);
		z -= iz;
		cz = 1 - z;

		u2 = u*u;
		mag = (u2 + 2) / sqrt(u2*(u2 + 4));
		mag *= ESPLout[ir][iz] * cr*cz + ESPLout[ir + 1][iz] * fr*cz + ESPLout[ir][iz + 1] * cr*z + ESPLout[ir + 1][iz + 1] * fr*z;
		if (astrometry) {
			astrox1 = u * (u2 + 3) / (u2 + 2);
			astrox1 *= ESPLoutastro[ir][iz] * cr*cz + ESPLoutastro[ir + 1][iz] * fr*cz + ESPLoutastro[ir][iz + 1] * cr*z + ESPLoutastro[ir + 1][iz + 1] * fr*z;
		}
	} 

	return mag;
}

double VBBinaryLensing::ESPLMag2(double u, double rho) {
	double Mag, u2,u6,rho2Tol;
	int c = 0;

	//Tol1_4 = sqrt(2 / Tol);
	//u2 = u*u;
	//u3Tol = u2*u*Tol;

	//if (u2 < Tol1_4) {
	//	rho2 = rho*rho;
	//	if (u3Tol > rho2*(1 + Tol1_4*rho)) {
	
	u2 = u*u;
	rho2Tol = rho*rho/Tol;
	u6 = u2*u2*u2;

	if (u6*(1+0.003*rho2Tol) > 0.027680640625*rho2Tol*rho2Tol) {
		Mag = (u2+2)/(u*sqrt(u2+4));
		if (astrometry) {
			astrox1 = u * (1 + 1 / (u2 + 2));
		}
	}
	else {
		Mag = ESPLMagDark(u, rho);
	}
	Mag0 = 0;
	return Mag;
}

double VBBinaryLensing::ESPLMagDark(double u, double RSv) {
	double Mag = -1.0, Magold = 0., Tolv = Tol;
	double tc, rb, lc, rc, cb,u2;
	int c = 0, flag;
	double currerr, maxerr;
	annulus *first, *scan, *scan2;
	int nannold, totNPS = 1;
        double LDastrox1=0.0;

		while ((Mag < 0.9) && (c < 3)) {

			first = new annulus;
			first->bin = 0.;
			first->cum = 0.;

			u2 = u * u;
			first->Mag = Mag0 = (u2 + 2) / (u*sqrt(u2 + 4));
			first->nim = 2;
			if (astrometry) {
				astrox1 = u * (u2 + 3) / (u2 + 2);
				first->LDastrox1 = astrox1 * first->Mag;
			}

			scr2 = sscr2 = 0;
			first->f = LDprofile(0);
			first->err = 0;
			first->prev = 0;


			first->next = new annulus;
			scan = first->next;
			scan->prev = first;
			scan->next = 0;
			scan->bin = 1.;
			scan->cum = 1.;
			scan->Mag = ESPLMag(u, RSv);//ESPLMag(u, RSv, Tolv, &Images);
			if (astrometry) {
				scan->LDastrox1 = astrox1 * scan->Mag;
			}
			scan->nim = 2;
			scr2 = sscr2 = 1;
			scan->f = LDprofile(0.9999999);
		scan->err = fabs((scan->Mag - scan->prev->Mag)*(scan->prev->f - scan->f) / 4);

		Magold = Mag = scan->Mag;
		if(astrometry){
			LDastrox1=scan->LDastrox1;			 
		}
		//	
		//			scan->err+=scan->Mag*Tolv*0.25; //Impose calculation of intermediate annulus at mag>4. Why?
		currerr = scan->err;
		flag = 0;
		nannuli = nannold = 1;
		while (((flag<nannold + 5) && (currerr>Tolv) && (currerr>RelTol*Mag)) || (nannuli<minannuli)) {
			maxerr = 0;
			for (scan2 = first->next; scan2; scan2 = scan2->next) {
#ifdef _PRINT_ERRORS_DARK
				printf("\n%d %lf %le | %lf %le", nannuli, scan2->Mag, scan2->err, Mag, currerr);
#endif
				if (scan2->err>maxerr) {
					maxerr = scan2->err;
					scan = scan2;
				}
			}

			nannuli++;
			Magold = Mag;
			Mag -= (scan->Mag*scan->bin*scan->bin - scan->prev->Mag*scan->prev->bin*scan->prev->bin)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
	    if(astrometry){
		        LDastrox1 -= (scan->LDastrox1*scan->bin*scan->bin - scan->prev->LDastrox1*scan->prev->bin*scan->prev->bin)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
				 
			}		
                        currerr -= scan->err;
			rc = scan->cum;
			lc = scan->prev->cum;
			tc = (lc + rc) / 2;
			cb = rCLDprofile(tc, scan->prev, scan);

			scan->prev->next = new annulus;
			scan->prev->next->prev = scan->prev;
			scan->prev = scan->prev->next;
			scan->prev->next = scan;
			scan->prev->bin = cb;
			scan->prev->cum = tc;
			scan->prev->f = LDprofile(cb);
			scan->prev->Mag = ESPLMag(u, RSv*cb);
                        if(astrometry){
				scan->prev->LDastrox1=astrox1*scan->prev->Mag;
				 
			}
			scan->prev->nim = 2;
			scan->prev->err = fabs((scan->prev->Mag - scan->prev->prev->Mag)*(scan->prev->prev->f - scan->prev->f)*(scan->prev->bin*scan->prev->bin - scan->prev->prev->bin*scan->prev->prev->bin) / 4);
			scan->err = fabs((scan->Mag - scan->prev->Mag)*(scan->prev->f - scan->f)*(scan->bin*scan->bin - scan->prev->bin*scan->prev->bin) / 4);
			rb = (scan->Mag + scan->prev->prev->Mag - 2 * scan->prev->Mag);
			scan->prev->err += fabs(rb*(scan->prev->prev->f - scan->prev->f)*(scan->prev->bin*scan->prev->bin - scan->prev->prev->bin*scan->prev->prev->bin));
			scan->err += fabs(rb*(scan->prev->f - scan->f)*(scan->bin*scan->bin - scan->prev->bin*scan->prev->bin));
#ifdef _PRINT_ERRORS_DARK
			printf("\n%d", Images->length);
#endif

			Mag += (scan->bin*scan->bin*scan->Mag - cb*cb*scan->prev->Mag)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
			Mag += (cb*cb*scan->prev->Mag - scan->prev->prev->bin*scan->prev->prev->bin*scan->prev->prev->Mag)*(scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin*scan->prev->bin - scan->prev->prev->bin*scan->prev->prev->bin);
                        if(astrometry){
				LDastrox1 += ( scan->bin*scan->bin*scan->LDastrox1 - cb*cb*scan->prev->LDastrox1)*(scan->cum - scan->prev->cum) / (scan->bin*scan->bin - scan->prev->bin*scan->prev->bin);
				LDastrox1 += ( cb*cb*scan->prev->LDastrox1 - scan->prev->prev->bin*scan->prev->prev->bin*scan->prev->prev->LDastrox1)*(scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin*scan->prev->bin - scan->prev->prev->bin*scan->prev->prev->bin);
				 
			}
			currerr += scan->err + scan->prev->err;

			if (fabs(Magold - Mag) * 2<Tolv) {
				flag++;
			}
			else {
				flag = 0;
				nannold = nannuli;
			}

		}

		while (first) {
			scan = first->next;
			delete first;
			first = scan;
		}

		Tolv /= 10;
		c++;
	}
	therr = currerr;
        if(astrometry){
		LDastrox1/=Mag;
		astrox1=LDastrox1;
		 
    }
	return Mag;
}

//////////////////////////////
//////////////////////////////
////////New (v2) light curve functions
//////////////////////////////
//////////////////////////////

void VBBinaryLensing::PSPLLightCurve(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double u0=exp(pr[0]), t0=pr[2],tE_inv=exp(-pr[1]),tn,u;

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - t0) *tE_inv;
		u = tn*tn + u0*u0;

		y1s[i] = -tn;
		y2s[i] = -u0;
		mags[i]= (u + 2) / sqrt(u*(u + 4));

	}
}


void VBBinaryLensing::PSPLLightCurveParallax(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double u0 = pr[0], t0 = pr[2], tE_inv = exp(-pr[1]), tn, u,u1,pai1=pr[3],pai2=pr[4];
	double Et[2];
	t0old = 0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);
		tn = (ts[i] - t0) * tE_inv + pai1*Et[0] + pai2*Et[1];
		u1 = u0 + pai1*Et[1] - pai2*Et[0];
		u = tn*tn + u1*u1;

		y1s[i] = -tn;
		y2s[i] = -u1;
		mags[i] = (u + 2) / sqrt(u*(u + 4));
	}

}


void VBBinaryLensing::ESPLLightCurve(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double u0 = exp(pr[0]), t0 = pr[2], tE_inv = exp(-pr[1]), tn, u,rho=exp(pr[3]);

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - t0) *tE_inv;
		u = sqrt(tn*tn + u0*u0);

		y1s[i] = -tn;
		y2s[i] = -u0;
		mags[i] = ESPLMag2(u,rho);

	}
}

void VBBinaryLensing::ESPLLightCurveParallax(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double u0 = pr[0], t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, u1, rho = exp(pr[3]), pai1 = pr[4], pai2 = pr[5];
	double Et[2];
	t0old = 0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);
		tn = (ts[i] - t0) * tE_inv + pai1*Et[0] + pai2*Et[1];
		u1 = u0 + pai1*Et[1] - pai2*Et[0];
		u = sqrt( tn*tn + u1*u1);

		y1s[i] = -tn;
		y2s[i] = -u1;
		mags[i] = ESPLMag2(u, rho);
	}

}


void VBBinaryLensing::BinaryLightCurve(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]);
	double salpha = sin(pr[3]), calpha = cos(pr[3]);

	//	_sols *Images; double Mag; // For debugging

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - pr[6]) * tE_inv;
		y1s[i] = pr[2] * salpha - tn*calpha;
		y2s[i] = -pr[2] * calpha - tn*salpha;
		mags[i] = BinaryMag2(s, q, y1s[i], y2s[i], rho);

		//Mag=BinaryMag(s, q, y1s[i], y2s[i], rho, Tol,&Images); // For debugging
		//delete Images;
		//mags[i] -= Mag;
		//if (fabs(mags[i]) > Tol) {
		//	printf("\n%lf %lf %lf", y1s[i], y2s[i], mags[i]);
		//}
	}
}


void VBBinaryLensing::BinaryLightCurveW(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]),t0,u0;
	double salpha = sin(pr[3]), calpha = cos(pr[3]),xc;

	xc = (s - 1 / s) / (1 + q);
	if (xc<0) xc = 0.;
	t0 = pr[6] + xc*calpha/tE_inv;
	u0 = pr[2] + xc*salpha;

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - t0) * tE_inv;
		y1s[i] = u0 * salpha - tn*calpha;
		y2s[i] = -u0 * calpha - tn*salpha;
		mags[i] = BinaryMag2(s, q, y1s[i], y2s[i], rho);
	}
}


void VBBinaryLensing::BinaryLightCurveParallax(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tn,u, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];
	t0old = 0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);
		tn = (ts[i] - t0) * tE_inv + pai1*Et[0] + pai2*Et[1];
		u = u0 + pai1*Et[1] - pai2*Et[0];
		y1s[i] = u * salpha - tn*calpha;
		y2s[i] = -u * calpha - tn*salpha;
		mags[i] = BinaryMag2(s, q, y1s[i], y2s[i], rho);
	}
}	


void VBBinaryLensing::BinaryLightCurveOrbital(double *pr, double *ts, double *mags, double *y1s, double *y2s, double *seps, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];
	double w, phi0, inc, phi, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm,s_true;
	double w13, w123, den, den0, u;
	t0old = 0;

	w13 = w1*w1 + w3*w3;
	w123 = sqrt(w13 + w2*w2);
	w13 = sqrt(w13);
	if (w13>1.e-8) {
		w3 = (w3>1.e-8) ? w3 : 1.e-8;
		w = w3*w123 / w13;
		inc = acos(w2*w3 / w13 / w123);
		phi0 = atan2(-w1*w123, w3*w13);
	}
	else {
		w = w2;
		inc = 0.;
		phi0 = 0.;
	}
	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(inc);
	Sinc = sin(inc);
	den0 = sqrt(Cphi0*Cphi0 + Cinc*Cinc*Sphi0*Sphi0);
	s_true = s / den0; // orbital radius
	COm = (Cphi0*calpha + Cinc*salpha*Sphi0) / den0;
	SOm = (Cphi0*salpha - Cinc*calpha*Sphi0) / den0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);

		phi = (ts[i] - t0_par)*w + phi0;
		Cphi = cos(phi);
		Sphi = sin(phi);
		den = sqrt(Cphi*Cphi + Cinc*Cinc*Sphi*Sphi);
		seps[i] = s_true*den; // projected separation at time ts[i]

		u = u0 + pai1*Et[1] - pai2*Et[0];
		tn = (ts[i] - t0) * tE_inv + pai1*Et[0] + pai2*Et[1];
		y1s[i] = (Cphi*(u*SOm - tn*COm) + Cinc*Sphi*(u*COm + tn*SOm)) / den;
		y2s[i] = (-Cphi*(u*COm + tn*SOm) - Cinc*Sphi*(tn*COm - u*SOm)) / den;
		mags[i] = BinaryMag2(seps[i], q, y1s[i], y2s[i], rho);
	}
}


void VBBinaryLensing::BinaryLightCurveKepler(double *pr, double *ts, double *mags, double *y1s, double *y2s, double *seps, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], alpha = pr[3], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11], szs = pr[12], ar = pr[13]+1.e-8;
	double Et[2];
	double u, w22, w11, w33, w12, w23, szs2, ar2, EE, dE;
	double wt2, smix, sqsmix, e, h, snu, co1EE0, co2EE0,cosE,sinE, co1tperi, tperi, EE0, M, a, St, psi, dM, conu, n;
	double arm1, arm2;
	double X[3], Y[3], Z[3],r[2],x[2];
	t0old = 0;

	smix = 1 + szs * szs;
	sqsmix = sqrt(smix);
	w22 = w2 * w2;
	w11 = w1 * w1;
	w33 = w3 * w3;
	w12 = w11 + w22;
	w23 = w22 + w33;
	wt2 = w12 + w33;

	szs2 = szs * szs;
	ar2 = ar * ar;
	arm1 = ar - 1;
	arm2 = 2 * ar - 1;
//	n = sqrt(wt2) / (ar*sqrt(-1 + 2 * ar)*sqrt(smix));
	n = sqrt(wt2 / arm2 / smix) / ar;
	Z[0] = -szs * w2;
	Z[1] = szs * w1 - w3;
	Z[2] = w2;
	h = sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
	for (int i = 0; i < 3; i++) Z[i] /= h;
	X[0] = -ar * w11 + arm1 * w22 - arm2 * szs*w1*w3 + arm1 * w33;
	X[1] = -arm2  *w2*(w1 + szs * w3);
	X[2] = arm1 * szs*w12 - arm2 * w1*w3 - ar * szs*w33;
	e = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	for (int i = 0; i < 3; i++) X[i] /= e;
	e /= ar * sqsmix*wt2;
	Y[0] = Z[1] * X[2] - Z[2] * X[1];
	Y[1] = Z[2] * X[0] - Z[0] * X[2];
	Y[2] = Z[0] * X[1] - Z[1] * X[0];

//	h = sqrt((smix)*w22 + (szs*w1 - w3)*(szs*w1 - w3));
//	co1e = (1 - ar)*(1 - ar) + ar2 * szs2 + (-1 + 2 * ar)*(w11*(1 - szs2) - szs2 * w22 + 2 * szs*w1*w3) / wt2;
//	co1nu = ar2 * szs2 + arm2 * (w11*(1 - szs2) - szs2 * w22 + 2 * szs*w1*w3) / wt2;
//	co1e = arm1*arm1 + co1nu;
//	coe2 = ar2 * smix;
//	e = sqrt(co1e/coe2);
//	co1nu = (-1 + 2 * ar)*sqrt(smix)*(w1 + szs * w3)*(szs2*(w12)-2 * szs*w1*w3 + w23);
//	co2nu = ar * e*h*(smix)*sqrt((smix))*wt2;
	conu = (X[0] + X[2] * szs) / sqsmix;
	co1EE0 = conu + e;
	co2EE0 = 1 + e * conu;
	cosE = co1EE0 / co2EE0;
	EE0 = acos(cosE);
	snu = (Y[0] + Y[2] * szs);
	EE0 *= (snu > 0) ? 1 : -1;
	sinE = sqrt(1 - cosE * cosE)*((snu > 0) ? 1 : -1);
	co1tperi = e * sinE;
	tperi = t0_par - (EE0 - co1tperi) / n;
//	coX = ar * e*sqrt(smix)*wt2;
	//coX1 = -ar * w11 + (-1 + ar)*w22 + (1 - 2 * ar)*szs*w1*w3 + (-1 + ar)*w33;
	//coX2 = (-1 + 2 * ar)*w1*w23 + szs2 * w1*((-1 + ar)*w12 - ar * w33) + szs * w3*((2 - 3 * ar)*w11 + ar * w23);
	//coY1 = -(-1 + 2 * ar)*w2*(w1 + szs * w3);
	//coY2 = w2 * (-szs2 * w12 + 2 * szs*w1*w3 - w23 + ar * (-4 * szs*w1*w3 + szs2 * (w12 - w33) + (-w11 + w23)));
	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);
		M = n * (ts[i] - tperi);
		EE = M + e * sin(M);
		dE = 1;
		while (fabs(dE) > 1.e-8) {
			dM = M - (EE - e * sin(EE));
			dE = dM / (1 - e * cos(EE));
			EE += dE;
		}

		a = ar * s*sqrt(smix);

		r[0] = a * (cos(EE) - e);
		r[1] = a * sqrt(1 - e * e)*sin(EE);
		x[0] = r[0] * X[0] + r[1] * Y[0];  // (coX1*x[1] + coX2 * y[1] / h) / coX;
		x[1] = r[0] * X[1] + r[1] * Y[1];   //(coY1*x[1] + y[1] * coY2 / h) / coX;
		St = sqrt(x[0] * x[0] + x[1] * x[1]);
		psi = atan2(x[1], x[0]);// +((ar > 1) ? 0 : M_PI);
		u = u0 + pai1 * Et[1] - pai2 * Et[0];
		tn = (ts[i] - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
		y1s[i] = -tn * cos(alpha + psi) + u * sin(alpha + psi);
		y2s[i] = -u * cos(alpha + psi) - tn * sin(alpha + psi);
		seps[i] = St;

		mags[i] = BinaryMag2(seps[i], q, y1s[i], y2s[i], rho);

	}
}

void VBBinaryLensing::BinSourceLightCurve(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double u1 = pr[2], u2=pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR=exp(pr[1]), tn, u;

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - t01) * tE_inv;
		u = tn*tn + u1*u1;

		y1s[i] = -tn;
		y2s[i] = -u1;
		mags[i] = (u + 2) / sqrt(u*(u + 4));

		tn = (ts[i] - t02) * tE_inv;
		u = tn*tn + u2*u2;

		mags[i] += FR*(u + 2) / sqrt(u*(u + 4));
		mags[i] /= (1 + FR);

	}

}


void VBBinaryLensing::BinSourceLightCurveParallax(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u, u0, pai1 = pr[6], pai2 = pr[7], w1 = pr[8], w2 = pr[9], w3 = pr[10];
	double Et[2];
	t0old = 0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et); 
		
		tn = (ts[i] - t01) * tE_inv + pai1*Et[0] + pai2*Et[1];
		u0 = u1 + pai1*Et[1] - pai2*Et[0];
		u = tn*tn + u0*u0;

		y1s[i] = -tn;
		y2s[i] = -u0;
		mags[i] = (u + 2) / sqrt(u*(u + 4));

		tn = (ts[i] - t02) * tE_inv + pai1*Et[0] + pai2*Et[1];
		u0 = u2 + pai1*Et[1] - pai2*Et[0];
		u = tn*tn + u0*u0;

		mags[i] += FR*(u + 2) / sqrt(u*(u + 4));
		mags[i] /= (1 + FR);
	}
}


void VBBinaryLensing::BinSourceLightCurveXallarap(double *pr, double *ts, double *mags, double *y1s, double *y2s,double *seps, int np) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u, u0, pai1 = pr[6], pai2 = pr[7], q = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double th,Cth,Sth;
	double Et[2];
	double s,s_true,w, phi0, inc, phi, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm;
	double w13, w123, den, den0,du0,dt0;
	t0old = 0;

	s = sqrt((u1 - u2)*(u1 - u2) + (t01 - t02)*(t01 - t02) * (tE_inv*tE_inv));
	th = atan2((u1 - u2) , (tE_inv*(t01 - t02)));
	Cth = cos(th);
	Sth = sin(th);
	u0 = (u1 + u2*q) / (1 + q);
	t0 = (t01 + t02*q) / (1 + q);

	w13 = w1*w1 + w3*w3;
	w123 = sqrt(w13 + w2*w2);
	w13 = sqrt(w13);
	if (w13>1.e-8) {
		w3 = (w3>1.e-8) ? w3 : 1.e-8;
		w = w3*w123 / w13;
		inc = acos(w2*w3 / w13 / w123);
		phi0 = atan2(-w1*w123, w3*w13);
	}
	else {
		w = w2;
		inc = 0.;
		phi0 = 0.;
	}
	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(inc);
	Sinc = sin(inc);
	den0 = sqrt(Cphi0*Cphi0 + Cinc*Cinc*Sphi0*Sphi0);
	s_true = s / den0;
	COm = (Cphi0*Cth + Cinc*Sth*Sphi0) / den0;
	SOm = (Cphi0*Sth - Cinc*Cth*Sphi0) / den0;


	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);

		phi = (ts[i] - t0_par)*w + phi0;
		Cphi = cos(phi);
		Sphi = sin(phi);
		den = sqrt(Cphi*Cphi + Cinc*Cinc*Sphi*Sphi);
		seps[i] = s_true*den;

		dt0 = s_true*(COm*Cphi-Cinc*SOm*Sphi)/(1+q)*q;  //Position of the primary component with respect to center of mass
		du0 = s_true*(SOm*Cphi + Cinc*COm*Sphi) / (1 + q)*q;

		tn= -((ts[i] - t0_par) * tE_inv + dt0 + pai1*Et[0] + pai2*Et[1]);
		u= -(u0+du0 + pai1*Et[1] - pai2*Et[0]);
		y1s[i] = tn;
		y2s[i] = u;
		u = tn*tn + u*u;

		mags[i] = (u + 2) / sqrt(u*(u + 4));

		tn =-( (ts[i] - t0_par) * tE_inv -dt0/q + pai1*Et[0] + pai2*Et[1]); // Position of the secondary component
		u =-( u0 - du0/q + pai1*Et[1] - pai2*Et[0]);
		u = tn*tn + u*u;

		mags[i] += FR*(u + 2) / sqrt(u*(u + 4));
		mags[i] /= (1 + FR);
	}
}


void VBBinaryLensing::BinSourceBinLensXallarap(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), u0;
	double salpha = sin(pr[3]), calpha = cos(pr[3]), xi1 = pr[7], xi2 = pr[8], omega = pr[9], inc = pr[10], phi = pr[11], qs = exp(pr[12]);

	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u02, rho2, tn2, y1s2, y2s2, qs4;



	if (t0_par_fixed == 0) t0_par = pr[6];


	for (int i = 0; i < np; i++) {

		phit = omega * (ts[i] - t0_par);

		disp[0] = cos(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));

		disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);

		Xal[0] = xi1 * disp[0] + xi2 * disp[1];
		Xal[1] = xi2 * disp[0] - xi1 * disp[1];
		tn = (ts[i] - pr[6]) * tE_inv + Xal[0];
		u0 = pr[2] + Xal[1];
		y1s[i] = u0 * salpha - tn * calpha;
		y2s[i] = -u0 * calpha - tn * salpha;
		Mag = BinaryMag2(s, q, y1s[i], y2s[i], rho);

		disp2[0] = -cos(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));

		disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

		Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
		Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];
		tn2 = (ts[i] - pr[6]) * tE_inv + Xal2[0];
		u02 = pr[2] + Xal2[1];
		y1s2 = u02 * salpha - tn2 * calpha;
		y2s2 = -u02 * calpha - tn2 * salpha;
		rho2 = rho * pow(qs, 0.89);
		Mag2 = BinaryMag2(s, q, y1s2, y2s2, rho2);
		qs4 = pow(qs, 4.0);
		mags[i] = (Mag + qs4 * Mag2) / (1 + qs4);
	}
}

void VBBinaryLensing::BinSourceSingleLensXallarap(double* pr, double* ts, double* mags, double* y1s, double* y2s, double* y1s2, double* y2s2, int np) {
	double  t0 = pr[1], rho = exp(pr[3]), tn, tE_inv = exp(-pr[2]), u0;
	double  xi1 = pr[4], xi2 = pr[5], omega = pr[6], inc = pr[7], phi = pr[8], qs = exp(pr[9]);

	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u02, rho2, tn2, qs4, u, u2;



	if (t0_par_fixed == 0) t0_par = pr[1];


	for (int i = 0; i < np; i++) {

		phit = omega * (ts[i] - t0_par);

		disp[0] = cos(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));

		disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);

		Xal[0] = xi1 * disp[0] + xi2 * disp[1];
		Xal[1] = xi2 * disp[0] - xi1 * disp[1];
		tn = (ts[i] - pr[1]) * tE_inv + Xal[0];
		u0 = pr[0] + Xal[1];
		u = sqrt(tn * tn + u0 * u0);

		y1s[i] = -tn;
		y2s[i] = -u0;
		Mag = ESPLMag2(u, rho);  //If you want only the second source put =0, otherwise replace ESPLMag2(u, rho);


		disp2[0] = -cos(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));

		disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

		Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
		Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];
		tn2 = (ts[i] - pr[1]) * tE_inv + Xal2[0];
		u02 = pr[0] + Xal2[1];
		u2 = sqrt(tn2 * tn2 + u02 * u02);
		y1s2[i] = -tn2;
		y2s2[i] = -u02;
		rho2 = rho * pow(qs, 0.89);
		Mag2 = ESPLMag2(u2, rho2);  //If you want only the second source put =0, otherwise replace ESPLMag2(u2, rho2);
		qs4 = pow(qs, 4.0);
		mags[i] = (Mag + qs4 * Mag2) / (1 + qs4);
	}
}

//////////////////////////////
//////////////////////////////
////////Old (v1) light curve functions
//////////////////////////////
//////////////////////////////


double VBBinaryLensing::PSPLLightCurve(double *pr, double t) {
	double u0 = exp(pr[0]), t0 = pr[2], tE_inv = exp(-pr[1]), tn, u;

	tn = (t - t0) *tE_inv;
	u = tn*tn + u0*u0;

	y_1 = -tn;
	y_2 = -u0;
	return (u + 2) / sqrt(u*(u + 4));

}


double VBBinaryLensing::PSPLLightCurveParallax(double *pr, double t) {
	double u0 = pr[0], t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, u1, pai1 = pr[3], pai2 = pr[4];
	double Et[2];

	ComputeParallax(t, t0, Et);
	tn = (t - t0) * tE_inv + pai1*Et[0] + pai2*Et[1];
	u1 = u0 + pai1*Et[1] - pai2*Et[0];
	u = tn*tn + u1*u1;

	y_1 = -tn;
	y_2 = -u1;
	return (u + 2) / sqrt(u*(u + 4));

}


double VBBinaryLensing::ESPLLightCurve(double *pr, double t) {
	double u0 = exp(pr[0]), t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, rho = exp(pr[3]);

	tn = (t - t0) *tE_inv;
	u = sqrt(tn*tn + u0*u0);

	y_1 = -tn;
	y_2 = -u0;
	return ESPLMag2(u, rho);

}

double VBBinaryLensing::ESPLLightCurveParallax(double *pr, double t) {
	double u0 = pr[0], t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, u1, rho = exp(pr[3]), pai1 = pr[4], pai2 = pr[5];
	double Et[2];

	ComputeParallax(t, t0, Et);
	tn = (t - t0) * tE_inv + pai1*Et[0] + pai2*Et[1];
	u1 = u0 + pai1*Et[1] - pai2*Et[0];
	u = sqrt(tn*tn + u1*u1);

	y_1 = -tn;
	y_2 = -u1;
	return ESPLMag2(u, rho);

}


double VBBinaryLensing::BinaryLightCurve(double *pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]);
	double salpha = sin(pr[3]), calpha = cos(pr[3]);

	tn = (t - pr[6]) * tE_inv;
	y_1 = pr[2] * salpha - tn*calpha;
	y_2 = -pr[2] * calpha - tn*salpha;
	return BinaryMag2(s, q, y_1, y_2, rho);

}


double VBBinaryLensing::BinaryLightCurveW(double *pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0, u0;
	double salpha = sin(pr[3]), calpha = cos(pr[3]), xc;

	xc = (s - 1 / s) / (1 + q);
	if (xc<0) xc = 0.;
	t0 = pr[6] + xc*calpha / tE_inv;
	u0 = pr[2] + xc*salpha;

	tn = (t - t0) * tE_inv;
	y_1 = u0 * salpha - tn*calpha;
	y_2 = -u0 * calpha - tn*salpha;
	return BinaryMag2(s, q, y_1, y_2, rho);

}


double VBBinaryLensing::BinaryLightCurveParallax(double *pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tn, u, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];

	ComputeParallax(t, t0, Et);
	tn = (t - t0) * tE_inv + pai1*Et[0] + pai2*Et[1];
	u = u0 + pai1*Et[1] - pai2*Et[0];
	y_1 = u * salpha - tn*calpha;
	y_2 = -u * calpha - tn*salpha;
	return BinaryMag2(s, q, y_1, y_2, rho);

}




double VBBinaryLensing::BinaryLightCurveOrbital(double *pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];
	double w, phi0, inc, phi, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm, s_true;
	double w13, w123, den, den0, u;

	w13 = w1*w1 + w3*w3;
	w123 = sqrt(w13 + w2*w2);
	w13 = sqrt(w13);
	if (w13>1.e-8) {
		w3 = (w3>1.e-8) ? w3 : 1.e-8;
		w = w3*w123 / w13;
		inc = acos(w2*w3 / w13 / w123);
		phi0 = atan2(-w1*w123, w3*w13);
	}
	else {
		w = w2;
		inc = 0.;
		phi0 = 0.;
	}

	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(inc);
	Sinc = sin(inc);
	den0 = sqrt(Cphi0*Cphi0 + Cinc*Cinc*Sphi0*Sphi0);
	s_true = s / den0;
	COm = (Cphi0*calpha + Cinc*salpha*Sphi0) / den0;
	SOm = (Cphi0*salpha - Cinc*calpha*Sphi0) / den0;

	ComputeParallax(t, t0, Et);

	phi = (t - t0_par)*w + phi0;
	Cphi = cos(phi);
	Sphi = sin(phi);
	den = sqrt(Cphi*Cphi + Cinc*Cinc*Sphi*Sphi);
	av = s_true*den;

	u = u0 + pai1*Et[1] - pai2*Et[0];
	tn = (t - t0) * tE_inv + pai1*Et[0] + pai2*Et[1];
	y_1 = (Cphi*(u*SOm - tn*COm) + Cinc*Sphi*(u*COm + tn*SOm)) / den;
	y_2 = (-Cphi*(u*COm + tn*SOm) - Cinc*Sphi*(tn*COm - u*SOm)) / den;
	return BinaryMag2(av, q, y_1, y_2, rho);
}

double VBBinaryLensing::BinaryLightCurveKepler(double *pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], alpha = pr[3], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11], szs = pr[12], ar = pr[13]+1.e-8;
	double Et[2];
	double u, w22, w11, w33, w12, w23, szs2, ar2, EE, dE;
	double wt2, smix, sqsmix, e, h, snu, co1EE0, co2EE0, cosE, sinE, co1tperi, tperi, EE0, M, a, St, psi, dM, conu, n;
	double arm1, arm2;
	double X[3], Y[3], Z[3], r[2], x[2];
	t0old = 0;

	smix = 1 + szs * szs;
	sqsmix = sqrt(smix);
	w22 = w2 * w2;
	w11 = w1 * w1;
	w33 = w3 * w3;
	w12 = w11 + w22;
	w23 = w22 + w33;
	wt2 = w12 + w33;

	szs2 = szs * szs;
	ar2 = ar * ar;
	arm1 = ar - 1;
	arm2 = 2 * ar - 1;
	n = sqrt(wt2 / arm2 / smix) / ar;
	Z[0] = -szs * w2;
	Z[1] = szs * w1 - w3;
	Z[2] = w2;
	h = sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
	for (int i = 0; i < 3; i++) Z[i] /= h;
	X[0] = -ar * w11 + arm1 * w22 - arm2 * szs*w1*w3 + arm1 * w33;
	X[1] = -arm2 * w2*(w1 + szs * w3);
	X[2] = arm1 * szs*w12 - arm2 * w1*w3 - ar * szs*w33;
	e = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	for (int i = 0; i < 3; i++) X[i] /= e;
	e /= ar * sqsmix*wt2;
	Y[0] = Z[1] * X[2] - Z[2] * X[1];
	Y[1] = Z[2] * X[0] - Z[0] * X[2];
	Y[2] = Z[0] * X[1] - Z[1] * X[0];

	conu = (X[0] + X[2] * szs) / sqsmix;
	co1EE0 = conu + e;
	co2EE0 = 1 + e * conu;
	cosE = co1EE0 / co2EE0;
	EE0 = acos(cosE);
	snu = (Y[0] + Y[2] * szs);
	EE0 *= (snu > 0) ? 1 : -1;
	sinE = sqrt(1 - cosE * cosE)*((snu > 0) ? 1 : -1);
	co1tperi = e * sinE;
	tperi = t0_par - (EE0 - co1tperi) / n;

	ComputeParallax(t, t0, Et);
	M = n * (t - tperi);
	EE = M + e * sin(M);
	dE = 1;
	while (fabs(dE) > 1.e-8) {
		dM = M - (EE - e * sin(EE));
		dE = dM / (1 - e * cos(EE));
		EE += dE;
	}
	
	a = ar * s*sqrt(smix);

	r[0] = a * (cos(EE) - e);
	r[1] = a * sqrt(1 - e * e)*sin(EE);
	x[0] = r[0] * X[0] + r[1] * Y[0];  // (coX1*x[1] + coX2 * y[1] / h) / coX;
	x[1] = r[0] * X[1] + r[1] * Y[1];   //(coY1*x[1] + y[1] * coY2 / h) / coX;
	St = sqrt(x[0] * x[0] + x[1] * x[1]);
	psi = atan2(x[1], x[0]);// +((ar > 1) ? 0 : M_PI);

	u = u0 + pai1 * Et[1] - pai2 * Et[0];
	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	y_1 = -tn * cos(alpha + psi) + u * sin(alpha + psi);
	y_2 = -u * cos(alpha + psi) - tn * sin(alpha + psi);
	
	return BinaryMag2(St, q, y_1, y_2, rho);
	
}

//double VBBinaryLensing::BinaryLightCurveKepler(double *pr, double t) {
//	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], alpha = pr[3], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11], szs = pr[12], ar = pr[13];
//	double Et[2];
//	double u, w22, w11, w33, w12, w23, szs2, ar2, coe2, coX, coX1, coX2, coY1, coY2, EE, dE;
//	double wt2, smix, e, h, co1e, co1nu, co2nu, co1EE0, co2EE0, co1tperi, tperi, EE0, nu, M, a, St, psi, dM, conu, n;
//	double x[3], y[3];
//	t0old = 0;
//
//	wt2 = w1 * w1 + w2 * w2 + w3 * w3;
//	smix = 1 + szs * szs;
//	w22 = w2 * w2;
//	w11 = w1 * w1;
//
//	w33 = w3 * w3;
//	w12 = w11 + w22;
//	w23 = w22 + w33;
//	szs2 = szs * szs;
//	ar2 = ar * ar;
//	n = sqrt(wt2) / (ar*sqrt(-1 + 2 * ar)*sqrt(smix));
//	h = sqrt((smix)*w22 + (szs*w1 - w3)*(szs*w1 - w3));
//	co1e = (1 - ar)*(1 - ar) + ar2 * szs2 + (-1 + 2 * ar)*(w11*(1 - szs2) - szs2 * w22 + 2 * szs*w1*w3) / wt2;
//	coe2 = ar2 * (smix);
//	e = sqrt(co1e) / sqrt(coe2);
//	co1nu = (-1 + 2 * ar)*sqrt(smix)*(w1 + szs * w3)*(szs2*(w12)-2 * szs*w1*w3 + w23);
//	co2nu = ar * e*h*(smix)*sqrt((smix))*wt2; 
//	nu = asin(co1nu / co2nu);
//	conu = cos(nu);
//	co1EE0 = conu + e;
//	co2EE0 = 1 + e * conu;
//	EE0 = acos(co1EE0 / co2EE0);
//	co1tperi = e * sin(EE0);
//	tperi = t0_par - (EE0 - co1tperi) / n;
//	coX = ar * e*sqrt(smix)*wt2;
//	coX1 = -ar * w11 + (-1 + ar)*w22 + (1 - 2 * ar)*szs*w1*w3 + (-1 + ar)*w33;
//	coX2 = (-1 + 2 * ar)*w1*w23 + szs2 * w1*((-1 + ar)*w12 - ar * w33) + szs * w3*((2 - 3 * ar)*w11 + ar * w23);
//	coY1 = -(-1 + 2 * ar)*w2*(w1 + szs * w3);
//	coY2 = w2 * (-szs2 * w12 + 2 * szs*w1*w3 - w23 + ar * (-4 * szs*w1*w3 + szs2 * (w12 - w33) + (-w11 + w23)));
//	
//	ComputeParallax(t, t0, Et);
//	M = n * (t - tperi);
//	EE = M + e * sin(M);
//	dE = 1;
//	while (fabs(dE) > 1.e-8) {
//		dM = M - (EE - e * sin(EE));
//		dE = dM / (1 - e * cos(EE));
//		EE += dE;
//	}
//
//	a = ar * s*sqrt(smix);
//
//	x[1] = a * (cos(EE) - e);
//	y[1] = a * sqrt(1 - e * e)*sin(EE);
//	x[2] = (coX1*x[1] + coX2 * y[1] / h) / coX;
//	y[2] = (coY1*x[1] + y[1] * coY2 / h) / coX;
//	St = sqrt(x[2] * x[2] + y[2] * y[2]);
//	psi = atan2(y[2], x[2])+ ((ar>1)? 0 : M_PI);
//	u = u0 + pai1 * Et[1] - pai2 * Et[0];
//	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
//	y_1 = -tn * cos(alpha + psi) + u * sin(alpha + psi);
//	y_2 = -u * cos(alpha + psi) - tn * sin(alpha + psi);
//
//	return BinaryMag2(St, q, y_1, y_2, rho);
//
//}

double VBBinaryLensing::BinSourceLightCurve(double *pr, double t) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u,mag;

	tn = (t - t01) * tE_inv;
	u = tn*tn + u1*u1;

	y_1 = -tn;
	y_2 = -u1;
	mag = (u + 2) / sqrt(u*(u + 4));

	tn = (t - t02) * tE_inv;
	u = tn*tn + u2*u2;

	mag += FR*(u + 2) / sqrt(u*(u + 4));
	mag /= (1 + FR);

	return mag;

}


double VBBinaryLensing::BinSourceLightCurveParallax(double *pr, double t) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u, u0, pai1 = pr[6], pai2 = pr[7], w1 = pr[8], w2 = pr[9], w3 = pr[10];
	double Et[2],mag;

	ComputeParallax(t, t0, Et);

	tn = (t - t01) * tE_inv + pai1*Et[0] + pai2*Et[1];
	u0 = u1 + pai1*Et[1] - pai2*Et[0];
	u = tn*tn + u0*u0;

	y_1 = -tn;
	y_2 = -u0;
	mag = (u + 2) / sqrt(u*(u + 4));

	tn = (t - t02) * tE_inv + pai1*Et[0] + pai2*Et[1];
	u0 = u2 + pai1*Et[1] - pai2*Et[0];
	u = tn*tn + u0*u0;

	mag += FR*(u + 2) / sqrt(u*(u + 4));
	mag /= (1 + FR);

	return mag;
}


double VBBinaryLensing::BinSourceLightCurveXallarap(double *pr, double t) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u, u0, pai1 = pr[6], pai2 = pr[7], q = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double th, Cth, Sth;
	double Et[2],mag;
	double s, s_true, w, phi0, inc, phi, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm;
	double w13, w123, den, den0, du0, dt0;

	s = sqrt((u1 - u2)*(u1 - u2) + (t01 - t02)*(t01 - t02) * (tE_inv*tE_inv));
	th = atan2((u1 - u2), (tE_inv*(t01 - t02)));
	Cth = cos(th);
	Sth = sin(th);
	u0 = (u1 + u2*q) / (1 + q);
	t0 = (t01 + t02*q) / (1 + q);

	w13 = w1*w1 + w3*w3;
	w123 = sqrt(w13 + w2*w2);
	w13 = sqrt(w13);
	if (w13>1.e-8) {
		w3 = (w3>1.e-8) ? w3 : 1.e-8;
		w = w3*w123 / w13;
		inc = acos(w2*w3 / w13 / w123);
		phi0 = atan2(-w1*w123, w3*w13);
	}
	else {
		w = w2;
		inc = 0.;
		phi0 = 0.;
	}
	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(inc);
	Sinc = sin(inc);
	den0 = sqrt(Cphi0*Cphi0 + Cinc*Cinc*Sphi0*Sphi0);
	s_true = s / den0;
	COm = (Cphi0*Cth + Cinc*Sth*Sphi0) / den0;
	SOm = (Cphi0*Sth - Cinc*Cth*Sphi0) / den0;

	ComputeParallax(t, t0, Et);

	phi = (t - t0_par)*w + phi0;
	Cphi = cos(phi);
	Sphi = sin(phi);
	den = sqrt(Cphi*Cphi + Cinc*Cinc*Sphi*Sphi);
	av = s_true*den;

	dt0 = s_true*(COm*Cphi - Cinc*SOm*Sphi) / (1 + q)*q;  //Position of the primary component with respect to center of mass
	du0 = s_true*(SOm*Cphi + Cinc*COm*Sphi) / (1 + q)*q;

	tn = -((t - t0_par) * tE_inv - dt0 + pai1*Et[0] + pai2*Et[1]);
	u = -(u0 + du0 + pai1*Et[1] - pai2*Et[0]);
	y_1 = tn;
	y_2 = u;
	u = tn*tn + u*u;

	mag = (u + 2) / sqrt(u*(u + 4));

	tn = -((t - t0_par) * tE_inv + dt0 / q + pai1*Et[0] + pai2*Et[1]); // Position of the secondary component
	u = -(u0 - du0 / q + pai1*Et[1] - pai2*Et[0]);
	u = tn*tn + u*u;

	mag += FR*(u + 2) / sqrt(u*(u + 4));
	mag /= (1 + FR);

	return mag;

}


double VBBinaryLensing::BinSourceBinLensXallarap(double* pr, double t) {

	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), u0;
	double salpha = sin(pr[3]), calpha = cos(pr[3]), xi1 = pr[7], xi2 = pr[8], omega = pr[9], inc = pr[10], phi = pr[11], qs = exp(pr[12]);

	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u02, rho2, tn2, y1s2, y2s2, y1s, y2s, mags, qs4;



	if (t0_par_fixed == 0) t0_par = pr[6];




	phit = omega * (t - t0_par);

	disp[0] = sin(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));

	disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);

	Xal[0] = xi1 * disp[0] + xi2 * disp[1];
	Xal[1] = xi2 * disp[0] - xi1 * disp[1];
	tn = (t - pr[6]) * tE_inv + Xal[0];
	u0 = pr[2] + Xal[1];
	y1s = u0 * salpha - tn * calpha;
	y2s = -u0 * calpha - tn * salpha;
	Mag = BinaryMag2(s, q, y1s, y2s, rho);

	disp2[0] = -sin(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));

	disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

	Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
	Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];
	tn2 = (t - pr[6]) * tE_inv + Xal2[0];
	u02 = pr[2] + Xal2[1];
	y1s2 = u02 * salpha - tn2 * calpha;
	y2s2 = -u02 * calpha - tn2 * salpha;
	rho2 = rho * pow(qs, 0.89);
	Mag2 = BinaryMag2(s, q, y1s2, y2s2, rho2);
	qs4 = pow(qs, 4.0);
	mags = (Mag + qs4*Mag2) / (1 + qs4);

	return mags;


}

double VBBinaryLensing::BinSourceSingleLensXallarap(double* pr, double t) {

	double  t0 = pr[1], rho = exp(pr[3]), tn, tE_inv = exp(-pr[2]), u0;
	double  xi1 = pr[4], xi2 = pr[5], omega = pr[6], inc = pr[7], phi = pr[8], qs = exp(pr[9]);

	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u02, rho2, tn2, y1s2, y2s2, qs4, u, y1s, y2s, mags, u2;



	if (t0_par_fixed == 0) t0_par = pr[1];

	phit = omega * (t - t0_par);

	disp[0] = sin(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));

	disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);

	Xal[0] = xi1 * disp[0] + xi2 * disp[1];
	Xal[1] = xi2 * disp[0] - xi1 * disp[1];
	tn = (t - pr[1]) * tE_inv + Xal[0];
	u0 = pr[0] + Xal[1];
	u = sqrt(tn * tn + u0 * u0);

	y1s = -tn;
	y2s = -u0;
	Mag = ESPLMag2(u, rho); //If you want only the second source put =0, otherwise replace ESPLMag2(u, rho);


	disp2[0] = -sin(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));

	disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

	Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
	Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];
	tn2 = (t - pr[1]) * tE_inv + Xal2[0];
	u02 = pr[0] + Xal2[1];
	u2 = sqrt(tn2 * tn2 + u02 * u02);
	y1s2 = -tn2;
	y2s2 = -u02;
	rho2 = rho * pow(qs, 0.89);
	Mag2 = ESPLMag2(u2, rho2); //If you want only the second source put =0, otherwise replace ESPLMag2(u2, rho2);
	qs4 = pow(qs, 4.0);
	mags = (Mag + qs4 * Mag2) / (1 + qs4);
	return mags;
}

double VBBinaryLensing::BinSourceBinLensPOX(double* pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tE_inv = exp(-pr[5]), t0 = pr[6];
	double pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];
	double tn, w, phi0, phil, incl, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm, s_true;
	double w13, w123, den, den0, u;

	double xi1 = pr[12], xi2 = pr[13], omega = pr[14], inc = pr[15], phi = pr[16], qs = exp(pr[17]);
	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u01, u02, rho2, tn1, tn2, mags, qs4;

	w13 = w1 * w1 + w3 * w3;
	w123 = sqrt(w13 + w2 * w2);
	w13 = sqrt(w13);
	if (w13 > 1.e-8) {
		w3 = (w3 > 1.e-8) ? w3 : 1.e-8;
		w = w3 * w123 / w13;
		incl = acos(w2 * w3 / w13 / w123);
		phi0 = atan2(-w1 * w123, w3 * w13);
	}
	else {
		w = w2;
		incl = 0.;
		phi0 = 0.;
	}

	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(incl);
	Sinc = sin(incl);
	den0 = sqrt(Cphi0 * Cphi0 + Cinc * Cinc * Sphi0 * Sphi0);
	s_true = s / den0;
	COm = (Cphi0 * calpha + Cinc * salpha * Sphi0) / den0;
	SOm = (Cphi0 * salpha - Cinc * calpha * Sphi0) / den0;

	ComputeParallax(t, t0, Et);

	phil = (t - t0_par) * w + phi0;
	Cphi = cos(phil);
	Sphi = sin(phil);
	den = sqrt(Cphi * Cphi + Cinc * Cinc * Sphi * Sphi);
	av = s_true * den;
	u = u0 + pai1 * Et[1] - pai2 * Et[0];
	tn = (t-t0)*tE_inv + pai1 * Et[0] + pai2 * Et[1];

	phit = omega * (t - t0_par);

	disp[0] = sin(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));
	disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);
	disp2[0] = -sin(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));
	disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

	Xal[0] = xi1 * disp[0] + xi2 * disp[1];
	Xal[1] = xi2 * disp[0] - xi1 * disp[1];
	Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
	Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];

	tn1 = tn + Xal[0];
	u01 = u + Xal[1];
	tn2 = tn + Xal2[0];
	u02 = u + Xal2[1];
	rho2 = rho * pow(qs, 0.89);
	qs4 = pow(qs, 4.0);


//	y1s = u01 * salpha - tn1 * calpha;
//	y2s = -u01 * calpha - tn1 * salpha;
//	y1s2 = u02 * salpha - tn2 * calpha;
//	y2s2 = -u02 * calpha - tn2 * salpha;

	y_1 = (Cphi * (u02 * SOm - tn2 * COm) + Cinc * Sphi * (u02 * COm + tn2 * SOm)) / den;
	y_2 = (-Cphi * (u02 * COm + tn2 * SOm) - Cinc * Sphi * (tn2 * COm - u02 * SOm)) / den;
	Mag2 = BinaryMag2(av, q, y_1, y_2, rho2);

	y_1 = (Cphi * (u01 * SOm - tn1 * COm) + Cinc * Sphi * (u01 * COm + tn1 * SOm)) / den;
	y_2 = (-Cphi * (u01 * COm + tn1 * SOm) - Cinc * Sphi * (tn1 * COm - u01 * SOm)) / den;
	Mag = BinaryMag2(av, q, y_1, y_2, rho);

	mags = (Mag + qs4 * Mag2) / (1 + qs4);

	return mags;
}



///////////////////////////////////////////////
///////////////////////////////////////////////
////////// Internal private functions
///////////////////////////////////////////////
///////////////////////////////////////////////

// dza=(z-a), where a=complex(-a1, 0)
// J1 = m1/(z-a)**2 + m2/(z-0)**2, which is just the conjugate of m1/(zc-a)**2 + m2/(zc-0)**2
// same in triple, m1/(z-z1)**2 + m2/(z-z2)**2 + m3/(z-z3)**2 is the conjugate of m1/(zc-z1c)**2 + m2/(zc-z2c)**2 + m3/(zc-z3c)**2
// J1 is Partial(zs_c)/Partial(z), J1c is Partial(zs)/Partial(z_c)
// dJ is Jacobian determinant (but dJ is declared as a complex variable)
// J2 is Partial^2(zs_c)/Partial(z)^2, which is the conjugate of Partial^2(zs)/Partial(z_c)^2
#define _Jacobians1 \
	z=zr[i];\
	dza=z-coefs[20];\
	za2 = dza*dza;\
	zb2=z*z;\
	J1= coefs[21]/za2+coefs[22]/zb2;\
	J1c=conj(J1);\
	dJ=1-J1*J1c;\
	J2=-2.*(coefs[21]/(za2*dza)+coefs[22]/(zb2*z));

/* only used when using BinaryMag */
// dy is zs'(theta) = i*rho*exp(i*theta)
// dz is z'(theta)  = { zs'(theta) - Partial(zs)/Partial(z_c) * conj( zs'(theta) ) } / dJ, 
// z'(theta) can be calculated from zs'(theta), z_c(at a certain theta)
// 
// set newly appended _point variable's attribute dJ(double) to dJ(complex).re, which is Jacobian determinant at this point/root
// similarily, attribute d is z'(theta)                    at this root
//						J2 is Partial^2(zs_c)/Partial(z)^2 at this root
// 						ds is x'(theta)^x"(theta)          at this root
// change the origin of z=(x1,x2) from the less massive lens to mass center
#define _Jacobians2\
	dy = complex(-sin(theta->th), cos(theta->th))*coefs[23];\
	dz = (dy - J1c*conj(dy)) / dJ.re;\
	Prov->last->x1 -= coefs[11].re;\
	Prov->last->dJ = dJ.re;\
	Prov->last->d = dz;\
	Prov->last->ds = (imag(dy*dz*dz*J2) + coefs[23].re*coefs[23].re) / dJ.re;

	/******************************************* changed *******************************************/
	//Prov->last->J2 = J2;
	/*******************************************   end   *******************************************/
        

/* Quadrupole test, only used when using BinaryMag0 */
// set newly appended _point variable's attribute dJ(double) to dJ(complex).re, which is Jacobian determinant at this point/root
// J3(original) = 6* [ m1/(z-a)**4 + m2/(z-0)**4 ] = Partial^3(zs_c)/Partial(z)^3, which is just the conjugate of Partial^3(zs)/Partial(z_c)^3
// dJ2(double) is Jacobian determinant**2
// J3(final) = J3(original) * J1c * J1c = Partial^3(zs_c)/Partial(z)^3 * [ Partial(zs)/Partial(z_c) ]**2
// ob2(double) = | Partial^2(zs_c)/Partial(z)^2 |**2 * [6 - 6*dJ(double) + dJ2(double)]
// J2(final) = J2(original)**2 * J1c**3 = [ Partial^2(zs_c)/Partial(z)^2 ]**2 * [ Partial(zs)/Partial(z_c) ]**3
//
// thus, equation(34) of Bozza et al., 2018 reads (correction of finite source from quadrupole term)
//                         mu_QI = (-1) * Re{6*conj(J2) - ob2(double) + 2*dJ(double)*conj(J3)} * rho**2 / dJ(double)**5
//							     = (-1) * {6*J2.re - ob2(double) + 2*dJ(double)*J3.re} * rho**2 / dJ(double)**5
//								 =        {ob2(double) - 6*J2.re - 2*dJ(double)*J3.re} * rho**2 / dJ(double)**5
// notice, each image has one of this quantity
//
// besides, as mu_QI can occasionally vanish (cusp approach by a fixed angle), we can construct a cusp correction term, 
// 						   mu_CI = Im{3*conj(J2)} * rho**2 / dJ(double)**5
// 							     = (-3) * J2.im * rho**2 / dJ(double)**5
// the final criteria is 
// 						   Sigma(I){cQ * (abs_mu_QI + abs_mu_CI)} < delta
//
// thus, we can let cq (correction of quadruple) = 
//						   0.5 * [ fabs{ob2(double) - 6*J2.re - 2*dJ(double)*J3.re} + fabs{(-3) * J2.im} ] / fabs{dJ(double)**5}
// and we have 
//                         abs_mu_QI + abs_mu_CI = 2 * cq * rho**2
// note: fabs() will always return a floating-point number even if the argument is an integer, 
//	     whereas abs() will return a floating-point or an integer depending upon the argument
// note: J2 has been permanently changed once you call _Jacobians3
#define _Jacobians3\
	Prov->last->dJ = dJ.re;\
	J3=6.*(coefs[21]/(za2*za2)+coefs[22]/(zb2*zb2));\
	dJ2=dJ.re*dJ.re;\
	za2=J1c*J1c;\
	J3=J3*za2;\
	ob2=(J2.re*J2.re+J2.im*J2.im)*(6-6*dJ.re+dJ2);\
	J2=J2*J2*za2*J1c;\
	cq= 0.5*(fabs(ob2-6.*J2.re-2.*J3.re*dJ.re)+3*fabs(J2.im))/fabs(dJ.re*dJ2*dJ2); 

/* ghost images test, only used when using BinaryMag0 */
// z_alt_c = yc - f(zG) = yc - (-1)*[ m1/(zG-a) + m2/(zG-0) ], which is just z^ defined in equation (41)
// J_alt_c = m1/(z_alt_c-a)**2 + m2/(z_alt_c-0)**2, which is just f'(z^) shown in equation (43)
// J1c = m1/(zG_c-a)**2 + m2/(zG_c-0)**2, i.e. Partial(zs)/Partial(zG_c) or f'(zG_c)
// JJalt2 = 1 - f'(zG_c) * conj(f'(z^)) = conj(J^), where J^ is defined in equation (44)
// 
// J2 is f"(zG), J1c is f'(zG_c), and JJalt2 is conj(J^), 
// thus J3 is conj[ J^ * f"(zG_c) * f'(zG) ], which is shown in equation (47), we call it conj[ denominator_component ]
// finally, J3(final) = { conj[ J^*f"(zG_c)*f'(zG) ] - J^*f"(zG_c)*f'(zG) * conj(f'(z^)) } / { conj(J^)*conj(J^)*dJ(zG,double) }, 
//          which makes 1/conj[ J3(final) ] just the thing inside abs in equation (47)
// thus, as cq should = abs(J3(final)), we have equation (47) in the form of 0.5 / cq > cG * rho
// 
// both ghost images should pass this test separately
#define _Jacobians4\
	zaltc=yc+coefs[21]/dza+coefs[22]/z;\
	za2=(zaltc-coefs[20]);\
	Jaltc=coefs[21]/(za2*za2)+coefs[22]/(zaltc*zaltc);\
	Jalt = conj(Jaltc);\
	JJalt2=(1-J1c*Jalt);\
	J3=J2*J1c*JJalt2;\
	J3=(J3-conj(J3)*Jalt)/(JJalt2*JJalt2*dJ.re);\
	cq=sqrt(J3.re*J3.re+J3.im*J3.im);
	/******************************************* changed *******************************************/
	// cq=(J3.re*J3.re+J3.im*J3.im);
	/*******************************************   end   *******************************************/
	

_curve* VBBinaryLensing::NewImages(complex yi, complex* coefs, _theta* theta) {//, float & time) {
	static complex  y, yc, z, zc, J1, J1c, dy, dz, dJ, J2, J3, dza, za2, zb2, zaltc, Jalt, Jaltc, JJalt2;
	static complex zr[5] = { 0.,0.,0.,0.,0. };
	static double dlmin = 1.0e-4, dlmax = 1.0e-3, good[5], dJ2, ob2, cq;
	static int worst1, worst2, worst3,  f1, checkJac;//bad,
	//static double av = 0.0, m1v = 0.0, disim, disisso;
	static _curve* Prov;
	static _point* scan, * prin, * fifth, * left, * right, * center; 

#ifdef _PRINT_TIMES
	static double tim0, tim1;
#endif

	// checkpoint1 (jump between BinaryMag0 and NewImages by finding this string in VScode)

	// BinaryMag0's input coordinate (y1v,y2v) takes the mass center as origin, 
    // then BinaryMag0 calls function NewImages and still use (y1v,y2v) as input, 
    // then in function NewImages, y = (y1v,y2v) + (-s*1/(1+q), 0) becomes the zs used in the coefficient expression, 
    // for example, the source on second  lens used to be ( s*1/(1+q), 0), but now becomes ( 0, 0)
    //              the source on primary lens used to be (-s*q/(1+q), 0), but now becomes (-s, 0)
	// (See a complete comment in BinaryMag0, here only takes q1<1 as an example)
	y = yi + coefs[11];
	yc = conj(y);

	// coefs[6]=a*a; coefs[7]=a*a*a; coefs[8]=m2*m2; coefs[9]=a*a*m2*m2; coefs[10]=a*m2; coefs[11]=a*m1; coefs[20]=a; coefs[21]=m1; coefs[22]=m2;
	
	coefs[0] = coefs[9] * y;
	coefs[1] = coefs[10] * (coefs[20] * (coefs[21] + y * (2 * yc - coefs[20])) - 2 * y);
	coefs[2] = y * (1 - coefs[7] * yc) - coefs[20] * (coefs[21] + 2 * y * yc * (1 + coefs[22])) + coefs[6] * (yc * (coefs[21] - coefs[22]) + y * (1 + coefs[22] + yc * yc));
	coefs[3] = 2 * y * yc + coefs[7] * yc + coefs[6] * (yc * (2 * y - yc) - coefs[21]) - coefs[20] * (y + 2 * yc * (yc * y - coefs[22]));
	coefs[4] = yc * (2 * coefs[20] + y);
	coefs[4] = yc * (coefs[4] - 1) - coefs[20] * (coefs[4] - coefs[21]);
	coefs[5] = yc * (coefs[20] - yc);

	//bad = 1;
	//disim = -1.;
	f1 = 0;

#ifdef _PRINT_TIMES
	tim0 = Environment::TickCount;
#endif
	//auto begin = std::chrono::high_resolution_clock::now() ;
										       /***************************************************************************/
	cmplx_roots_gen(zr, coefs, 5, true, true); // Question: polish using original polynomial, so how to guarantee not to polish to same root???
											   // Note: VBBL use previous point's roots as new point's initial guess!!!
											   //       remeber to test the performance with use_roots_as_starting_points=ture/false
											   // Another question is what's the scope of reusing? a source contour or a lightcurve?
											   // Still reuse for multi-call of BinaryMag2 in a C++ test code? what about a Python test code?
											   /***************************************************************************/
	// cmplx_roots_gen(complex *roots, complex *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points)
	// use_roots_as_starting_points=ture will turn off resetting zr[5]={ 0.,0.,0.,0.,0. } at the beginning of cmplx_roots_gen
	// Thus even with "static complex zr[5] = { 0.,0.,0.,0.,0. };", 
	// after the first call of NewImages, the five roots are retained in zr[] and used as initial guess in the second call of NewImages
	// The reason is: zr[] is static!


	//auto end   = std::chrono::high_resolution_clock::now() ;
	//auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin) ;
	//time = (float)(elapsed_time.count() * 1e-9) ; // nanosecond to second
    //printf("cmplx_roots_gen needs %e (second)\n", total_time) ;

#ifdef _PRINT_TIMES
	tim1 = Environment::TickCount;
	inc += tim1 - tim0;
#endif
	// apply lens equation to check if it is really solved
	for (int i = 0; i < 5; i++) {             // only cares about the worst, second-worst, and third-worst roots
											  // ranked by good[i] which is just abs(lens_equation) 
		z = zr[i];
		zc = conj(z);
		good[i] = abs(_LL); // Lens equation check
		// #define _LL (y-z)+coefs[21]/(zc-coefs[20])+coefs[22]/zc 
		// Lens equation test, y/z/zc all take secondary lens as origin, only used in NewImages()
		switch (i) {
		case 0:
			worst1 = i;
			break; 							  // break out of the switch block
		case 1:
			if (good[i] > good[worst1]) {
				worst2 = worst1;
				worst1 = i;
			}
			else worst2 = i;
			break;
		case 2:
			if (good[i] > good[worst1]) {      // good[i] > good[worst1]
				worst3 = worst2;
				worst2 = worst1;
				worst1 = i;
			}
			else if (good[i] > good[worst2]) { //           good[worst1] >= good[i] > good[worst2]
				worst3 = worst2;
				worst2 = i;
			}
			else worst3 = i;                   //                                     good[worst2] >= good[i]
			break;

		default: //case 3, 4
			if (good[i] > good[worst1]) {
				worst3 = worst2;
				worst2 = worst1;
				worst1 = i;
			}
			else if (good[i] > good[worst2]) {
				worst3 = worst2;
				worst2 = i;
			}
			else if (good[i] > good[worst3]) {
				worst3 = i;
			}
		}
	}

	Prov = new _curve; 							// create an object of class _curve on heap
												// pointer to that object is assigned to static local variable 'Prov', 
												// notice: no bracket, so leaves attributes in in-determinate state.
												//         (see https://www.quora.com/What-is-the-difference-between-malloc-new-and-operator-new-in-C-C)
	checkJac = 0;
	//	if (!((good[worst3] < dlmin) && ((good[worst1] < dlmin) || (good[worst2] > dlmax)))) {  // old check for unacceptable roots

	// 3 good roots
	if (good[worst2] * dlmin > (good[worst3]+1.e-12) ) { // good[worst2] > (1.e4*good[worst3] + 1.e-8), so three good roots
														 // Adams 1967 (take a third degree polynomial as an example): 
														 // abs(polynomial) < 1.0e-15 * [ abs(a3x**3+a2x**2+a1x+a0) + abs(a3x**3+a2x**2+a1x) + abs(a3x**3+a2x**2) + abs(a3x**3) ]
		for (int i = 0; i < 5; i++) {
			if ((i != worst1) && (i != worst2)) {        // if current zr[i] is one of the three good roots
				//if((i==worst3)&&(good[i]>dlmax)&&(good[worst2]>1.e2*good[worst3])){
				//	zr[i]=(coefs[21].re<coefs[22].re)? 0.5*coefs[20]+coefs[21]/(0.5*coefs[20]-yc-coefs[22]/coefs[20]) : -0.5*coefs[20]+coefs[22]/(-0.5*coefs[20]-yc+coefs[21]/coefs[20]);
				//}
				Prov->append(zr[i].re, zr[i].im); 		 // create a new _point variable on heap, 
				// just (*Prov).append();				 // and then append it to the END of the linked list 'Prov'(object of class _curve)
														 //
														 // constructor of _point class is called, 
														 // assign input values (the good root zr[i]) to attributes x1, x2, 
														 // and set attribute theta(pointer) = NULL
														 //
														 // has no return value

				_Jacobians1								// J1 is Partial(zs_c)/Partial(z), J1c is Partial(zs)/Partial(z_c)
														// dJ is Jacobian determinant (but dJ is declared as a complex variable)
														// J2 is Partial^2(zs_c)/Partial(z)^2

					if (theta->th >= 0) {				// the formal parameter 'theta' that points to a _theta variable on heap
														//   in BinaryMag0, 
														// input pointer 'theta' points to a _theta variable whose attribute 'th'=-1
														//   in BinaryMag,
														// input pointer 'theta' points to two _theta variables successively, 
														// means starting from two 'th' : 0.01020304 and 0.01020304 + M_PI,
														// while the Thetas(_thetas variable) first consists of three elements: 0.01020304, 0.01020304 + M_PI, 0.01020304 + 2*M_PI,
														//   then new 'th' is calculated from 'th = (itheta->th + itheta->next->th) / 2;' and inserted to Thetas(_thetas variable), 
														// where itheta has the largest 'maxerr' in all elements(except 'last') of Thetas(_thetas variable), 
														// note itheta is acquired through looping all elements(except 'last') 

						_Jacobians2						// BinaryMag
														// change the origin of (attribute) z=(x1,x2) from the less massive lens to mass center
														// assign value to attribute d, which is z'(theta)                    at this root
														//							J2  		 Partial^2(zs_c)/Partial(z)^2 at this root
														// 							ds 			 x'(theta)^x"(theta)          at this root
														//							dJ			 Jacobian determinant         at this root
					}
					else {
						_Jacobians3						// BinaryMag0, thus perform the Quadrupole test
														// calculate cq for this root, where (abs_mu_QI + abs_mu_CI) = 2 * cq * rho**2
														//
							corrquad += cq;				// sum current root's cq to the total cq
														// the passing criteria is cQ * 2 * corrquad * rho**2 < delta
														// corrquad is set to 0 in BinaryMag0, just in front of calling NewImages
					}
				
				checkJac += (fabs(Prov->last->dJ) > 1.e-7) ? _sign(Prov->last->dJ) : 10;
														// if current root's dJ:
														//						     1.0e-7 < dJ, checkJac += 1
														//			  -1.0e-7 < dJ < 1.0e-7,      checkJac += 10
														//       dJ < -1.0e-7, 					  checkJac += -1
														// checkJac is thus the summation of each root's parity 
														// (dJ>0 is positive parity, while dJ<0 is negative parity)
														// 
														// if there are 3 good roots, 
														// there will be one positive parity image and two negative parity images, 
														// so checkJac should be -1
														// 
														// if there are 5 good roots, 
														// the new pair of images will be one positive parity and one negative parity, 
														// so still should have checkJac=-1

				Prov->last->theta = theta;				// set newly appended _point variable's attribute theta(pointer) to formal parameter 'theta'
														// finally make these three _point variables all point to the same _theta variable 

			}
		}

		if (theta->th < 0) {							// BinaryMag0, thus perform the ghost images test
			dz = zr[worst2] - zr[worst1];				// not used

			int i = worst1;
			_Jacobians1									// calculate zr[worst1]'s J1, J1c, dJ(complex), J2
				_Jacobians4								// calculate cq for this ghost image, 
														// where the passing criteria is 0.5 / cq > cG * rho
														// and the two ghost images should pass this criteria separately
				corrquad2 = cq;

			i = worst2;
			_Jacobians1									// calculate zr[worst2]'s J1, J1c, dJ(complex), J2
				_Jacobians4
				if (cq > corrquad2) corrquad2 = cq;		// corrquad2 = max{cq_worst1, cq_worst2}, 
														// thus if 0.5 / corrquad2 > cG * rho, both ghost images will pass the criteria
			//_Jacobians3
			//corrquad2 +=  1/cq;

		}
		else {
			theta->errworst = abs(zr[worst1] - zr[worst2]);	
														// BinaryMag
														// distance between two ghost roots
														// aim is to estimate how far zs is from a caustic
														// see equation (37) of Bozza, 2010
		}

	} 
	else { 
		 if (good[worst2]*dlmax > good[worst3] + 1.e-12 && theta->th>=0) { // Dubious cases. Better exclude them 
														// (1.e3*good[worst3] + 1.e-9) < good[worst2] <= (1.e4*good[worst3] + 1.e-8), 
														// and in BinaryMag (we have lots of thetas so can exclude one), 
														// which means good[worst2] is somehow bad but not bad enough. 
														// thus abort all five real/ghost roots and
														// return the pointer to an empty _curve variable
												
			return Prov;
		}
		else {		// 5 good roots						// good[worst2] <= (1.e3*good[worst3] + 1.e-9) in BinaryMag
														// or 
														// good[worst2] <= (1.e4*good[worst3] + 1.e-8) in BinaryMag0, 
														// treat these two cases both as 5 good roots
			f1 = 0;
			for (int i = 0; i < 5; i++) {
				Prov->append(zr[i].re, zr[i].im);

				_Jacobians1
					if (theta->th >= 0) {
						_Jacobians2						// BinaryMag, same as above (but on 5 good roots)
														// change the origin of (attribute) z=(x1,x2) from the less massive lens to mass center
														// assign value to attribute d, which is z'(theta)                    at this root
														//							J2  		 Partial^2(zs_c)/Partial(z)^2 at this root
														// 							ds 			 x'(theta)^x"(theta)          at this root
														//							dJ			 Jacobian determinant         at this root
					}
					else {
						_Jacobians3						// BinaryMag0, thus perform the Quadrupole test
														// (as 5 good roots, no ghost images test!)
														// calculate cq for this root, where (abs_mu_QI + abs_mu_CI) = 2 * cq * rho**2
														//
							corrquad += cq;				// sum current root's cq to the total cq (all 5 roots)
														// the passing criteria is cQ * 2 * corrquad * rho**2 < delta
														// corrquad is set to 0 in BinaryMag0, just in front of calling NewImages
					}
				checkJac += (fabs(Prov->last->dJ) > 1.e-7) ? _sign(Prov->last->dJ) : 10;
														// same as 3 good roots case, see above

				Prov->last->theta = theta;				// finally make these five _point variables all point to the same _theta variable 

				if (fabs(dJ.re) < 1.e-5) f1 = 1;		// if one (or more) of the five dJ(double) is close to 0, set f1=1
			}

			theta->errworst = -1.e100;					// 5 good roots, both BinaryMag0 and BinaryMag will set theta->errworst to a incredible valve

			// check Jacobians in ambiguous cases
			if (f1) {									// if one (or more) of the five dJ(double) is close to 0 ( -1.0e-5 < dJ < 1.0e-5 )
														// but none of them is too close to 0 ( all have abs(dJ) > 1.0e-7 ), 
														// otherwise checkJac will != -1 and we don't care about anymore
				left = right = center = fifth = 0;
				dJ.re = 0;
				for (scan = Prov->first; scan; scan = scan->next) {	
														// loop through all 5 good roots (scan points to _point variables)
														// 
														// if the root's y coordinate has the same sign as zs's y 
														// (only one root can have this property in binary lensing with both lenses on x-axis), 
														// then use 'prin' to point to this _point variable
														// 
														// find out the root with largest abs(dJ) in the remaining 4 roots (smallest magnification)
														// use 'fifth' to point to this _point variable, and use dJ.re to store this largest abs(dJ)
					if (_sign(scan->x2) == _sign(y.im)) {
						prin = scan;
					}
					else {
						dz.re = fabs(scan->dJ);
						if (dz.re > dJ.re) {
							fifth = scan;
							dJ.re = dz.re;
						}
					}
				}
				for (scan = Prov->first; scan; scan = scan->next) {
					if ((scan != prin) && (scan != fifth)) {
														// loop through the remaining 3 good roots (scan points to _point variables)
														// use 'left'/'center'/'right' to point to these 3 _point variables, 
														// according to the remaining 3 roots' x coordinate
						if (left) {
							if (scan->x1 < left->x1) {
								if (left != right) {
									center = left;
								}
								left = scan;
							}
							else {
								if (scan->x1 > right->x1) {
									if (left != right) {
										center = right;
									}
									right = scan;
								}
								else {
									center = scan;
								}
							}
						}
						else {
							left = right = center = scan;
						}
					}
				}
														// left and right images should have negative parity, 
														// while center image should have positive parity 
														// (see Schneider & Weiss, 1986, Fig. 6)
														// Question: 
														// But why can we manipulate dJ's sign for these 3 images?
														// And checkJac has already been calculated before, so why bother to manipulate dJ's sign now? 
				if (left->dJ > 0) left->dJ = -left->dJ;	
				if (center->dJ < 0) center->dJ = -center->dJ;
				if (right->dJ > 0) right->dJ = -right->dJ;
			}
		}
	}
	if (checkJac != -1) {								// checkJac, which is the summation of each root's parity, != -1
														// might be caused by 
														// one (or more) root's dJ is too close to 0 (too close to 0 means -1.0e-7 < dJ < 1.0e-7 )
														// or can be caused by 
														// wrongly polishing a positive parity image to a negative parity image, or vice versa
//		printf("\ncheckJac!");
		if (theta->th < 0) {							// BinaryMag0
														// dJ.re is total magnification
														// if total magnification is enough close to 1, 
														// keep all elements and return the pointer to the _curve variable
														// (and set corrquad = 0)
														// 
														// or drop and delete all elements of curve 'Prov'
														// thus finally return the pointer to an empty _curve variable
			dJ = 0;
			for (scan = Prov->first; scan; scan = scan->next) {
				dJ = dJ+ 1 / fabs(scan->dJ);
			}
			if (fabs(dJ.re - 1) < Tol) {
				checkJac = -1;
				corrquad = 0;
			}
		}
		if(checkJac!=-1){								// BinaryMag
														// drop and delete all elements of curve 'Prov'
														// thus finally return the pointer to an empty _curve variable
			_point* scan2;
			for (scan = Prov->first; scan; scan = scan2) {
				scan2 = scan->next;
				Prov->drop(scan);
				delete scan;
			}
		}
	}
	return Prov;										// return the pointer to a _curve variable which is created on heap during NewImages()
														// this _curve variable is a linked list of good roots (_point variables) at this zs
														// this _curve variable can contain 0, 3, or 5 elements
														// make all elements (_point variables) point to the same _theta variable
														// 
														//
														//
														// if input actual/formal parameter 'theta' points to a _theta variable with th<0, 
														// means NewImages() is used in BinaryMag0, 
														// 
														// 				   (1.e4*good[worst3] + 1.e-8) < good[worst2], so three good roots
														// 
														// three good roots(_point variables)' attribute Jacobian determinant (at each root) is calculated
														// perform the Quadrupole test on three good roots: 
														// 		calculate cq for each root, where (abs_mu_QI + abs_mu_CI) = 2 * cq * rho**2
														// 		sum three roots' cq to the total cq, corrquad
														// 		the passing criteria is cQ * 2 * corrquad * rho**2 < delta
														// 
														// perform the ghost images test on two ghost roots:
														// 		calculate cq for each ghost image, where the passing criteria is 0.5 / cq > cG * rho
														// 	and the two ghost images should pass this criteria separately;
														// 		then take corrquad2 = max{cq_worst1, cq_worst2}, 
														// 	thus if 0.5 / corrquad2 > cG * rho, both ghost images will pass the criteria
														//
														// (if) one (or more) of the three good roots' dJ is too close to 0 (too close to 0 means -1.0e-7 < dJ < 1.0e-7 )
														// or wrongly polishing a positive parity image to a negative parity image, or vice versa, 
														// will cause the summation of three root's parity != -1, then check the total magnification of three good roots:
														// if total magnification is enough close to 1: 
														// 		keep three good roots and return the pointer to the _curve variable (also set corrquad = 0)
														// or:
														// 		drop and delete all elements of curve 'Prov', thus finally return the pointer to an empty _curve variable
														// 
														//
														//
														// good[worst2] <= (1.e4*good[worst3] + 1.e-8), so five good roots
														// 
														// five good roots(_point variables)' attribute Jacobian determinant (at each root) is calculated
														// perform the Quadrupole test on five good roots, and no ghost images test(corrquad2=0)
														// set theta->errworst to -1.e100
														// 
														// if one (or more) of the five dJ is close to 0 ( -1.0e-5 < dJ < 1.0e-5 ), 
														// manipulate dJ's sign for 3/5 images (but after summing five roots' parity)
														//
														// if the summation of five roots' parity != -1, then check the total magnification of five good roots, same as above
														//
														//
														//
														//
														// otherwise means used in BinaryMag,
														// 												 (1.e4*good[worst3] + 1.e-8) < good[worst2], so three good roots
														// three good roots(_point variables)' attributes 
														// z'(theta), Partial^2(zs_c)/Partial(z)^2, x'(theta)^x"(theta), Jacobian determinant (at each root) is calculated,  
														// origin of (attribute) z=(x1,x2) is changed from the less massive lens to mass center;
														// and distance between two ghost roots is assigned to theta->errworst
														// 
														// (if) one (or more) root's dJ is too close to 0 (too close to 0 means -1.0e-7 < dJ < 1.0e-7 )
														// or wrongly polishing a positive parity image to a negative parity image, or vice versa, 
														// will cause the summation of three root's parity != -1, finally return the pointer to an empty _curve variable
														// 
														// (1.e3*good[worst3] + 1.e-9) < good[worst2] <= (1.e4*good[worst3] + 1.e-8)
														// which means good[worst2] is somehow bad but not bad enough.
														// thus abort all five real/ghost roots and return the pointer to an empty _curve variable
														// 
										// good[worst2] <= (1.e3*good[worst3] + 1.e-9), so five good roots
														// same as three good roots case, but set theta->errworst to -1.e100
														// and if one (or more) of the five dJ is close to 0 ( -1.0e-5 < dJ < 1.0e-5 )
														// but none of them is too close to 0 ( all have abs(dJ) > 1.0e-7 ), manipulate dJ's sign for 3/5 images
														// (otherwise checkJac will != -1 and we don't care about anymore, just return the pointer to an empty _curve variable)
														// 
														// Note: VBBL use previous point's roots as new point's initial guess!!!
														// 		 after the first call of NewImages(), the five roots are retained in zr[] 
														//       and used as initial guess in the second call of NewImages() because zr[] is static
}

/******************************************* changed *******************************************/
//void VBBinaryLensing::OrderImages(_sols *Sols, _curve *Newpts) {
void VBBinaryLensing::OrderImages(_sols_for_skiplist_curve * Sols, _curve * Newpts) {
	static double A[5][5];
	static _skiplist_curve *cprec[5];
	static _skiplist_curve *cpres[5];
	static _skiplist_curve *cfoll[5];
	//static _curve *cprec[5];
	//static _curve *cpres[5];
	//static _curve *cfoll[5];
	static _point *scan, *scan2, *isso[2];//*scan3, 
	static _skiplist_curve *scurve, *scurve2 ;
	//static _curve *scurve, *scurve2;

	static std::minstd_rand engine{std::random_device{}()} ;
/*******************************************   end   *******************************************/

	_theta *theta;
	static double th, mi, cmp, cmp2,cmp_2,dx2,avgx2,avgx1,avg2x1,pref,d2x2,dx1,d2x1,avgwedgex1,avgwedgex2,parab1,parab2;
        
	int nprec = 0, npres, nfoll = 0, issoc[2], ij;
	// int l2 ;

	/******************************************* changed *******************************************/
	int new_and_append_Level = 0 ;
	while ( new_and_append_Level < max_skiplist_level && (engine() % 4) == 0 ) 
	{
		new_and_append_Level++ ;
	} 
	/*******************************************   end   *******************************************/

	// checkpoint2 (jump between BinaryMag and OrderImages by finding this string in VScode)

	theta = Newpts->first->theta; 						// pointer to _theta variable
														// all elements in _curve variable point to the same _theta variable
	th = theta->th;


	theta->Mag = theta->prev->Mag = theta->maxerr = theta->prev->maxerr = 0; 
        theta->astrox1 = theta->prev->astrox1 = theta->astrox2 = theta->prev->astrox2 = 0;
														// two segments need to be updated after inserting one theta
	if (Newpts->length == 3) {
		mi = theta->next->errworst - theta->errworst;  	// distance between two ghost roots was assigned to theta->errworst
														// if theta->next also has only 3 good roots, then mi is g(i+1)-g(i)
														// 												or mi is ~ -1.e100
		if ((mi > theta->errworst) && (theta->prev->errworst > 0.)) {
			theta->prev->maxerr = mi*mi;				
		}												// if 'theta', 'theta->next' both have only 3 good roots, 
														// and g(i+1)-g(i)>g(i),
														// but 'theta->prev' also has only 3 good roots, 
														// which means zs(i-1) outside the caustic
														// 
														// interval [i-1, i] might thus be undersampled, 
														// then [g(i+1)-g(i)]^2 is added to maxerr(i-1)

		mi = theta->prev->errworst - theta->errworst; 	// equation (38) in Bozza, 2010
														// if theta->prev also has only 3 good roots, then mi is g(i-1)-g(i)
														// 												or mi is ~ -1.e100
		if ((mi>theta->errworst) && (theta->next->errworst>0.)) {
			theta->maxerr = mi*mi;
		}												// if 'theta', 'theta->prev' both have only 3 good roots, 
														// and g(i-1)-g(i)>g(i),
														// but 'theta->next' also has only 3 good roots, 
														// which means zs(i+1) outside the caustic
														// 
														// interval [i, i+1] might thus be undersampled, 
														// then [g(i-1)-g(i)]^2 is added to maxerr(i)
	}
	// (i-1)--------------------i---------(i+0.5)-------(i+i)
	//        	[i-1, i]				[i, i+1]
	// 
	// if (i-1) and i both outside the caustic, and g(i-1) - g(i) > g(i), then we suspect (i+1) will be inside the caustic;  
	// while if (i+1) is actually outside the caustic, we then suspect the interval [i, i+1] is undersampled
	// 
	// Symmetrically, 
	// if (i+1) and i both outside the caustic, and g(i+1) - g(i) > g(i), then we suspect (i-1) will be inside the caustic;  
	// while if (i-1) is actually outside the caustic, we then suspect the interval [i-1, i] is undersampled
	//
	// So, whether an interval [i, i+1] is undersampled is determined by (i-1)->i and (i+2)->(i+1)    Question(following): 
	// After inserting (i+0.5), we need to calculate interval [i+0.5, i+1]'s error by     i->(i+0.5), but what about the contribution by (i+2)->(i+1)?
	// 														  [i, i+0.5]		      (i+1)->(i+0.5), but what about the contribution by (i-1)->   i ?
	// 
	//
	//
	// the _thetas variable now contains 3 elements, looks like: 
	// 																	(theta points to the following element)
	// th          0.01020304  										  	0.01020304 + M_PI								0.01020304 + 2*M_PI
	// maxerr      ?													? 											  	0.
	// Mag         0.													0.												0.
	// astrox1,2   (0.,0.)												(0.,0.)			                               	(0.,0.)
	// errworst    distance between two ghost roots (3 good roots)    	distance between two ghost roots (3 good roots)	distance between two ghost roots (3 good roots)
	//			   -1.e100   						(5 good roots)   	-1.e100   					  	 (5 good roots)	-1.e100   					  	 (5 good roots)
	

	// For each image, we find the place to insert the new points
	scurve = Sols->first;
	for (int i = 0; i<Sols->length; i++) 
	{
		if (th < scurve->first->theta->th) 				// can occur when current 'scurve' is the 4th, 5th image contours
		{
			if (th > scurve->first->theta->prev->prev->th) 	
			{
				cfoll[nfoll] = scurve; // image involved at start
				nfoll++;								// the current 'scurve' is only a following curve of new inserted points
														// (not a preceding curve)
				scurve2 = scurve->next;					
				Sols->drop(scurve);						// drop 'scurve' because it's a following curve of new inserted points
				i--;
				scurve = scurve2;
			}											// if        scurve->first->theta->prev->prev < 'th' < scurve->first->theta
														// 
														// if 'th' sits between the last 'th' outside caustic and the first 'th' inside caustic 
														// (last 'th' outside caustic < 'th' < first 'th' inside caustic), 
														// then 'th' has the possibility to either be inside/outside caustic: 
														// 	scurve->first->theta             is the first 'th' inside caustic
														// 	scurve->first->theta->prev       is the newly inserted 'th'
														// 	scurve->first->theta->prev->prev is the original last 'th' outside caustic
														// (as newly inserted 'th' was calculated by two end points, 
														//  scurve->first->theta->prev->prev must exist! )
														// 
														// if 'th' < the last 'th' outside caustic, we will thus have 'th' < the first 'th' inside caustic, 
														// 'th' must be outside caustic if the whole source contour only has one time caustic crossing
														// (but 'th' can be inside caustic if the second caustic crossing occurs) : 
														// but what about scurve->first->theta->prev->prev? 
														// 	scurve->first->theta             is the first 'th' inside caustic
														// 	scurve->first->theta->prev       is the original last 'th' outside caustic
														// 	scurve->first->theta->prev->prev might be the newly inserted 'th' or a 'th' outside caustic
														// so the if (th > scurve->first->theta->prev->prev->th) can't be true
														

			else 
			{
				scurve = scurve->next;
			}											// if 'th' <= scurve->first->theta->prev->prev 
														//
														// just the situation in previous paragraph
		}
		else 
		{
			if (th > scurve->last->theta->th) 			// can occur when current 'scurve' is the 4th, 5th image contours
			{
				if (th < scurve->last->theta->next->next->th) 
				{
					cprec[nprec] = scurve; // image involved at end
					nprec++;							// the current 'scurve' is only a preceding curve of new inserted points
														// (not a following curve)

				}										// if scurve->last->theta < 'th' < scurve->last->theta->next->next
			}


			else 										// can occur when current 'scurve' is 1st to 5th image contours
														// scurve->first->theta <= 'th' <= scurve->last->theta
			{
				// image involved at the centre			// insertion happens between itheta->th and itheta->next->th
														// if we record the 3/5 pointers to corresponding _point variables in _theta variable, 
														// we can try to for(i): scan = itheta->pointer[i] to divide each _curve variable?
														// we cannot divide 5 curves one by one just by knowing 5 pointers from _theta variable, 
														// without knowing which curve each pointer belongs to (see source code of divide)
														//
														// skiplist is the answer, 
														// because 5 _curve variables will be broken and reconnected with each other in an unpredicted form
														// and we can not relabel every _point variables behind 'itheta' after each reconnection
														// skiplist is as easy to divide and join just as normal linked list (I have thinked this carefully), 
														// but dividing and joining a Binary Search Tree is much harder 
														//
														// now the remaining question is:
														// can augmented skiplist know the Nelement behind an element after dividing and joining?
														// or do we really need the length information for the divided curve?
														// 
														// to summarize, the key requirement is: 
														// given an aribitrary theta, for each curve, we need to know at which place to divide the curve
														// if we only know five pointers, we can not know each pointer belongs to which curve
														// (if you insist to do this, 
														// like let each _point know which curve it belongs to through some kind of augmentation, 
														// you will need to update all _point variables after a reconnection, that's O(n))
														// or we just loop through each curve, 
														// you will know for each curve, at which pointer to divide the curve
				/******************************************* changed *******************************************/
				// scan = scurve->last;
				// l2 = 0 ; // l2 is Nmember from scan->next to last (excluding the one pointed by 'scan')
				// while (scan->theta->th > th) {
				// 	scan = scan->prev;
				// 	l2++ ;
				// }
				// cfoll[nfoll] = scurve->divide(scan, l2);
				cfoll[nfoll] = scurve->find_prev_then_divide(th) ;
				/*******************************************   end   *******************************************/
				nfoll++;
				cprec[nprec] = scurve;
				nprec++;								// the current 'scurve' (first to scan)      is a preceding curve of new inserted points
														// the new curve 'nc'   (scan->next to last) is a following curve of new inserted points
			}
			scurve = scurve->next;
		}
	}
	npres = Newpts->length;								// 1. Sols->length=5, but nprec=3, nfoll=3,     npres=3 can occur (if-else or else-if-else)
														// 2.                                       but npres=5 can occur
														// 3. Sols->length=5,     nprec=3, nfoll=5,     npres=5 can occur (if-if)
														// 4.                                           npres=3 can occur 
														// 5. Sols->length=5,     nprec=5, nfoll=3,     npres=5 can occur (else-if-if)
														// 6.                                           npres=3 can occur 
														// 7. Sols->length=5,     nprec=5, nfoll=5,     npres=5 can occur (else-else)
														// 8.                                       but npres=3 can occur
														// 9. Sols->length=3,     nprec=3, nfoll=3,     npres=5 can occur 
														//10.                                           npres=3 can occur

	//if((theta->th>4.7116917419)&&(theta->th<4.711691759)){
	//	theta->th=theta->th;
	//}
	//if((theta->prev->th>4.7116917419)&&(theta->prev->th<4.711691759)){
	//	theta->th=theta->th;
	//}


	// case of creating new images

	if (nprec<npres) {									// previous case 2,3,9
		mi = 1.e100;	
		scan = Newpts->first;
		for (int i = 0; i<Newpts->length - 1; i++) {	// i=0,1,2,3
			scan2 = scan->next;
			for (int j = i + 1; j<Newpts->length; j++) {// i=0, j=1,2,3,4;   i=1, j=2,3,4;   i=2, j=3,4;   i=3, j=4
				cmp = (*scan2) - (*scan);				// operator-: distance^2 between 'current _point' and 'input _point'
				if (cmp<mi) {
					mi = cmp;
					isso[0] = scan;
					isso[1] = scan2;
				}
				scan2 = scan2->next;
			}
			scan = scan->next;
		}												// find the min distance^2 pair in all C_5^2=10 possible pairs
														// 'mi' is the min distance^2, isso[0] and isso[1] point to that pair
														// this pair is regarded as the two new created images
		Newpts->drop(isso[0]);
		Newpts->drop(isso[1]);
		/******************************************* changed *******************************************/
		//scurve = new _curve(isso[0]);					// new _skiplist_curve
		scurve  = new _skiplist_curve(isso[0], new_and_append_Level) ;	
		//isso[0]->prev = isso[0]->next = 0;				// weren't set to NULL when dropping, but already set to NULL when new _curve
		//scurve2 = new _curve(isso[1]);
		scurve2 = new _skiplist_curve(isso[1], new_and_append_Level) ;
		//isso[1]->prev = isso[1]->next = 0;
		/*******************************************   end   *******************************************/
		scurve->partneratstart = scurve2;
		scurve2->partneratstart = scurve;				// make two new curves 'scurve', 'scurve2' become each other's 'partner_at_start'
		Sols->append(scurve);
		Sols->append(scurve2);							// in following case, Sols->length becomes 7: 
														// 	Sols->length=5, but nprec=3, nfoll=3, but npres=5 can occur (if-else or else-if-else)
														// in other two cases, Sols->length still =5
		cpres[3] = scurve;
		cpres[4] = scurve2;
		scan = isso[0];
		scan2 = isso[1];

		cmp2 = fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
		cmp = sqrt(mi / cmp2);   // Delta theta tilde	// equation (29) in Bozza, 2010
 
		cmp_2 = cmp*cmp;
		mi = cmp_2*cmp*0.04166666667;					// 1/24 * (delta theta ~)^3
		parab1 = -(-scan->ds + scan2->ds)*mi; 			// here assume: scan as +, scan2 as -
		parab2=-0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
														// Ac(p2)'s sign not checked
		scurve->parabstart = 0.5 * (parab1 + parab2);	// in if(nprec<npres) case, stored in _curve attribute: scurve->parabstart,             
														//		where scurve=new _curve(isso[0]) (and) isso[0] and isso[1] are just dropped from Newpts
														// in if(nprec>npres) case, stored in _point attribute: scan->parab, 					
														// 		where scan=cprec[issoc[0]]->last (and) issoc[0] and issoc[1] are the closest pair among C_5^2=10 possible cprec[i]->last pairs
														// in nprec==npres case (essentially), stored in _point attribute: scan->parab, 
														// 		where scan=cprec[issoc[0]]->last (and) cprec[issoc[0]]->last and isso[1] is the closest cprec[]->last -- Newpts pair 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: created image:
        if(astrometry){
			avgwedgex1=-(-scan->x1*scan->ds + scan2->x1*scan2->ds)*mi;
			avgwedgex2=-(-scan->x2*scan->ds + scan2->x2*scan2->ds)*mi;
			dx2=-(-scan->d.im + scan2->d.im);
			d2x2=dx2*dx2; 
			dx1=-(-scan->d.re + scan2->d.re);
			d2x1=dx1*dx1; 
			scurve->parabastrox1 =-0.125*d2x1*dx2*mi-avgwedgex1;
			scurve->parabastrox2 =-0.125*d2x2*dx1*mi+avgwedgex2;
		}
#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scurve->parabstart, (scan->ds + scan2->ds)*mi / 2, fabs(scurve->parabstart)*(cmp*cmp) / 10, 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) - 2 * cmp*cmp2)*cmp);
#endif

		mi = fabs((parab1-parab2)*0.5) + fabs(scurve->parabstart)*(cmp_2 *0.1) + 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) - 2 * cmp*cmp2)*cmp;
														// mi = E1new(c) + E3(c) + E2(c)
#ifdef _noparab
		mi = fabs(scurve->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scurve->parabstart = 0.;
#endif
 
#ifdef _selectimage
		if (_selectionimage)
#endif
		pref=(scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5;	// trapezium approximation, have checked sign
        theta->prev->Mag -= ((scan->dJ>0) ? -1 : 1)*(pref + scurve->parabstart);
														// according to the assumption (scan as +, scan2 as -), 
														// 	scan->dJ > 0 thus Mag -= -(A(t)+A(p)), thus correct sign; 
														// if finally scan as -, scan2 as +, 
														// 	scan->dJ < 0 thus Mag -= (-A(t)-A(p)), thus still correct sign

		theta->prev->maxerr += mi;						// add the crossing-critical-curve term to theta->prev

#ifdef _ERRORS_ANALYTIC
		char filnam[32];
		sprintf(filnam, "%02dprev.txt", NPS);
		FILE* f = fopen(filnam, "a+");
		fprintf(f, "%.15le %.15le %.15le\n", -((scan->dJ > 0) ? -1 : 1)* (pref), -((scan->dJ > 0) ? -1 : 1)* (scurve->parabstart), mi);
		fclose(f);
#endif
		scurve2->parabstart = -scurve->parabstart;
               
                if(astrometry){
		        dx2=scan2->x2 - scan->x2;
		        avgx1=scan->x1 + scan2->x1;
		        avg2x1=avgx1*avgx1;
		        avgx2=scan->x2 + scan2->x2;
		        theta->prev->astrox1 += ((scan->dJ>0) ? -1 : 1)*(avg2x1*dx2*0.125+scurve->parabastrox1);
		        theta->prev->astrox2 -= ((scan->dJ>0) ? -1 : 1)*(pref*avgx2*0.25+scurve->parabastrox2);
		        scurve2->parabastrox2=-scurve->parabastrox2;
		        scurve2->parabastrox1=-scurve->parabastrox1;
                }

		 

	}

														// 1. Sols->length=5, but nprec=3, nfoll=3,     npres=3 can occur (if-else or else-if-else)
														// 2.                                       but npres=5 can occur
														// 3. Sols->length=5,     nprec=3, nfoll=5,     npres=5 can occur (if-if)
														// 4.                                           npres=3 can occur 
														// 5. Sols->length=5,     nprec=5, nfoll=3,     npres=5 can occur (else-if-if)
														// 6.                                           npres=3 can occur 
														// 7. Sols->length=5,     nprec=5, nfoll=5,     npres=5 can occur (else-else)
														// 8.                                       but npres=3 can occur
														// 9. Sols->length=3,     nprec=3, nfoll=3,     npres=5 can occur 
														//10.                                           npres=3 can occur
	// case of image destruction 
	if (nprec>npres) {									// previous case 6,8
		mi = 1.e100;
		for (int i = 0; i<nprec - 1; i++) {				// i=0,1,2,3
			for (int j = i + 1; j<nprec; j++) {			// i=0, j=1,2,3,4;   i=1, j=2,3,4;   i=2, j=3,4;   i=3, j=4
				cmp = *(cprec[i]->last) - *(cprec[j]->last);
														// operator-: distance^2 between 'current _point' and 'input _point'
				if (cmp<mi) {
					mi = cmp;
					issoc[0] = i;
					issoc[1] = j;
				}
			}
		}												// find the min distance^2 pair in all C_5^2=10 possible cprec[i]->last pairs
														// 'mi' is the min distance^2, cprec[issoc[0]]->last and cprec[issoc[1]]->last point to that pair
														// this pair is regarded corresponding to the two destructed images
		cprec[issoc[0]]->partneratend = cprec[issoc[1]];
		cprec[issoc[1]]->partneratend = cprec[issoc[0]];// make two curves become each other's 'partner_at_end'

		scan = cprec[issoc[0]]->last;
		scan2 = cprec[issoc[1]]->last;

		cmp2 = fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
		cmp = sqrt(mi / cmp2);							// equation (29) in Bozza, 2010
		cmp_2 = cmp*cmp;
		mi = cmp_2*cmp *0.04166666666667;				// 1/24 * (delta theta ~)^3
		parab1 = -(scan->ds - scan2->ds)*mi;
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);			// in if(nprec<npres) case, stored in _curve attribute: scurve->parabstart
														// in if(nprec>npres) case, stored in _point attribute: scan->parab
														// 
														// in if(nprec<npres) case, stored in _curve attribute: scurve->parabstart,             
														//		where scurve=new _curve(isso[0]) (and) isso[0] and isso[1] are just dropped from Newpts
														// in if(nprec>npres) case, stored in _point attribute: scan->parab, 					
														// 		where scan=cprec[issoc[0]]->last (and) issoc[0] and issoc[1] are the closest pair among C_5^2=10 possible cprec[i]->last pairs
														// in nprec==npres case (essentially), stored in _point attribute: scan->parab, 
														// 		where scan=cprec[issoc[0]]->last (and) cprec[issoc[0]]->last and isso[1] is the closest cprec[]->last -- Newpts pair
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: destructed image:
                if(astrometry){
		        avgwedgex1=-(scan->x1*scan->ds - scan2->x1*scan2->ds)*mi;
		        avgwedgex2=-(scan->x2*scan->ds - scan2->x2*scan2->ds)*mi;
			dx2=-(scan->d.im-scan2->d.im);
		        d2x2=dx2*dx2;
		        dx1=-(scan->d.re-scan2->d.re);
		        d2x1=dx1*dx1; 
		        scan->parabastrox1 =-0.125*d2x1*dx2*mi-avgwedgex1;
		        scan->parabastrox2 =-0.125*d2x2*dx1*mi+avgwedgex2;
		       
               }
#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds + scan2->ds)*mi / 2, fabs(scan->parab)*(cmp*cmp) / 10, 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) + 2 * cmp*cmp2)*cmp);
#endif

		mi = fabs((parab1-parab2)*0.5) + fabs(scan->parab)*(cmp*cmp *0.1) + 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) + 2.0 * cmp*cmp2)*cmp;
														// mi = E1new(c) + E3(c) + E2(c)
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
		pref=(scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5;	// trapezium approximation
		theta->prev->Mag += ((scan->dJ>0) ? -1 : 1)*(pref+ scan->parab);
#ifdef _ERRORS_ANALYTIC
		char filnam[32];
		sprintf(filnam, "%02dprev.txt", NPS);
		FILE* f = fopen(filnam, "a+");
		fprintf(f, "%.15le %.15le %.15le\n", ((scan->dJ > 0) ? -1 : 1)* (pref), ((scan->dJ > 0) ? -1 : 1)* (scan->parab), mi);
		fclose(f);
#endif
		if(astrometry){
			dx2=scan2->x2 - scan->x2;
			avgx1=scan->x1 + scan2->x1;
			avg2x1=avgx1*avgx1;
			avgx2=scan->x2 + scan2->x2;
			theta->prev->astrox1 -= ((scan->dJ>0) ? -1 : 1)*(avg2x1*dx2*0.125+scan->parabastrox1);
			theta->prev->astrox2 += ((scan->dJ>0) ? -1 : 1)*(pref*avgx2*0.25+scan->parabastrox2);
		}
		theta->prev->maxerr += mi;						// add the crossing-critical-curve term to theta->prev
		scan2->parab = -scan->parab;
        if(astrometry){
			scan2->parabastrox2 =-scan->parabastrox2;
			scan2->parabastrox1 =-scan->parabastrox1;
		}

 
		nprec -= 2;										// nprec from 5 to 3
														// originally, cprec[issoc[0]]->last and cprec[issoc[1]]->last correspond to the two destructed images
														// and issoc[0] and issoc[1] can be any index from 0 to 4; 
														// 
														// now, cprec[0], cprec[1], cprec[2] point to the three normal image contours, 
														// while cprec[3] and cprec[4] are meaningless
		ij = 0;
		for (int i = 0; i<nprec; i++) {
			if (i == issoc[0]) ij++;
			if (i == issoc[1] - 1) ij++;
			cprec[i] = cprec[i + ij];
		}
	}


														// 1. Sols->length=5, but nprec=3, nfoll=3,     npres=3 can occur (if-else or else-if-else)
														// 2.                                       but npres=5 can occur
														// 3. Sols->length=5,     nprec=3, nfoll=5,     npres=5 can occur (if-if)
														// 4.                                           npres=3 can occur 
														// 5. Sols->length=5,     nprec=5, nfoll=3,     npres=5 can occur (else-if-if)
														// 6.                                           npres=3 can occur 
														// 7. Sols->length=5,     nprec=5, nfoll=5,     npres=5 can occur (else-else)
														// 8.                                       but npres=3 can occur
														// 9. Sols->length=3,     nprec=3, nfoll=3,     npres=5 can occur 
														//10.                                           npres=3 can occur

	// constructing distance matrix with previous images
														// in if(nprec<npres) case, i.e. case 2,3,9, 
														// 	finally nprec=3, npres(still)=5, Newpts->length from =5 to =3, but cpres[3,4] are already meaningful
														// 
														// in if(nprec>npres) case, i.e. case 6,8, 
														// 	finally nprec from =5 to =3, npres=3, and only cprec[0,1,2] are meaningful, cprec[3,4] are meaningless
														// 
														// in remaining cases, i.e. case 1,4,10
														// nprec=npres=3			
														//                     i.e. case 5,7
														// nprec=npres=5							
	mi = 1.e100;
	for (int i = 0; i<nprec; i++) {	// loop through cprec[i]->last
		cpres[i] = cprec[i];							// notice! very tricky, see comments just in front of if(nfoll)
		scan = Newpts->first;
		for (int j = 0; j<nprec; j++) {	// loop through Newpts
			A[i][j] = (signbit(cprec[i]->last->dJ) == signbit(scan->dJ))? *(cprec[i]->last) - *scan : 100;
			if (A[i][j]<mi) {
				mi = A[i][j];
				issoc[0] = i;
				issoc[1] = j;
				isso[1] = scan;
			}
			scan = scan->next;
		}
	}													// construct 3*3 or 5*5 distance matrix between cprec[]->last (as i) and Newpts (as j)
														// matrix value is distance^2 for same parity, 100 for different parity
														// 
														// 'mi' is the min distance^2, which is also stored in A[issoc[0]][issoc[1]] 
														// cprec[issoc[0]]->last and isso[1] point to that pair
														
	// association with previous images
	while (nprec) {
		scan = cprec[issoc[0]]->last;
		scan2 = isso[1];								// start from the closest pair

		cmp2 = mi / fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
		cmp = (scan->theta->th - scan2->theta->th);		// - delta theta
		cmp_2 = cmp*cmp;
		mi = cmp_2*cmp *0.0416666666666667;  ////// (1/24 cube(delta Teta))
														// - 1/24 * (delta theta)^3
		parab1 = (scan->ds + scan2->ds)*mi; // Old parabolic correction, -A(p1)
		// New parabolic correction
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im - scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re - scan->d.re)) * cmp;
														// - 1/12 * something * (delta theta) = -A(p2), have checked something's sign
		scan->parab = 0.5*(parab1+parab2); 	// -A(p)	// in if(nprec<npres) case, stored in _curve attribute: scurve->parabstart
														// in if(nprec>npres) case, stored in _point attribute: scan->parab
														// in current nprec==npres case (essentially), stored in _point attribute: scan->parab
														// 
														// in if(nprec<npres) case, stored in _curve attribute: scurve->parabstart,             
														//		where scurve=new _curve(isso[0]) (and) isso[0] and isso[1] are just dropped from Newpts
														// in if(nprec>npres) case, stored in _point attribute: scan->parab, 					
														// 		where scan=cprec[issoc[0]]->last (and) issoc[0] and issoc[1] are the closest pair among C_5^2=10 possible cprec[i]->last pairs
														// in nprec==npres case (essentially), stored in _point attribute: scan->parab, 
														// 		where scan=cprec[issoc[0]]->last (and) cprec[issoc[0]]->last and isso[1] is the closest cprec[]->last -- Newpts pair
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: ordinary image:
        if(astrometry){
			avgwedgex1=(scan->x1*scan->ds + scan2->x1*scan2->ds)*mi;
		    avgwedgex2=(scan->x2*scan->ds + scan2->x2*scan2->ds)*mi;
			dx2=scan->d.im+scan2->d.im;
			d2x2=dx2*dx2;
			dx1=scan->d.re+scan2->d.re;
			d2x1=dx1*dx1; 
			scan->parabastrox1 =-0.125*d2x1*dx2*mi-avgwedgex1;
			scan->parabastrox2 =-0.125*d2x2*dx1*mi+avgwedgex2;
		}
#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds - scan2->ds)*mi / 2, fabs(scan->parab)*(cmp2) / 10, fabs(scan->parab)*(1.5*fabs(cmp2 / (cmp*cmp) - 1)));
#endif

		/******************************************* changed *******************************************/
		// Question: typo? should be cmp_2 instead of cmp2
		// original: 
		//mi = fabs((parab1 - parab2) *0.5) + fabs( scan->parab * (cmp2 *0.1 + 1.5*fabs(cmp2 / (cmp_2) - 1)) );
		// should be: 
		mi = fabs((parab1 - parab2) *0.5) + fabs( scan->parab * (cmp_2 *0.1 + 1.5*fabs(cmp2 / (cmp_2) - 1)) );
														// mi = E1new + E3 + E2
		/*******************************************   end   *******************************************/
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
		pref=(scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5;	// trapezium approximation
																	// -A(t)
		theta->prev->Mag += ((scan->dJ>0) ? -1 : 1)*( pref+ scan->parab);	// dJ>0 thus Mag += -(-A(t)-A(p)), thus correct sign
														// add this term to theta->prev
#ifdef _ERRORS_ANALYTIC
				char filnam[32];
				sprintf(filnam, "%02dprev.txt", NPS);
				FILE* f = fopen(filnam, "a+");
				fprintf(f, "%.15le %.15le %.15le\n", ((scan->dJ > 0) ? -1 : 1)* (pref), ((scan->dJ > 0) ? -1 : 1)* (scan->parab),mi);
				fclose(f);
#endif
		if(astrometry){
			dx2=scan2->x2 - scan->x2;
			avgx1=scan->x1 + scan2->x1;
			avg2x1=avgx1*avgx1;
			avgx2=scan->x2 + scan2->x2;
			theta->prev->astrox1 -= ((scan->dJ>0) ? -1 : 1)*(avg2x1*dx2*0.125+scan->parabastrox1);
			theta->prev->astrox2 += ((scan->dJ>0) ? -1 : 1)*(pref*avgx2*0.25+scan->parabastrox2);
        }
		theta->prev->maxerr += mi;						// add this term to theta->prev
                 


		Newpts->drop(isso[1]);
		cprec[issoc[0]]->append(isso[1], new_and_append_Level);
		cprec[issoc[0]]->partneratend = 0;				// append the new point to the end of the curve
														// partner_at_end set to NULL, not correct for case 5, but will be adjusted later

		nprec--;
		for (int i = issoc[0]; i<nprec; i++) {
			cprec[i] = cprec[i + 1];
			for (int j = 0; j<nprec + 1; j++) {
				A[i][j] = A[i + 1][j];
			}
		}												// nprec from 3 to 2 (or 5 to 4)
														// drop cprec[issoc[0]], and now cprec[0], cprec[1] are meaningful while cprec[2] is meaningless
														// 							 (or cprec[0-3])						 (or cprec[4])				
														// drop the corresponding row in A[][], now A[][]'s meaningful range is 2*3 (or 4*5)
																											
		for (int j = issoc[1]; j<nprec; j++) {
			for (int i = 0; i<nprec; i++) {
				A[i][j] = A[i][j + 1];
			}
		}												// drop the corresponding column in A[][], now A[][]'s meaningful range is 2*2 (or 4*4)

		mi = 1.e100;
		for (int i = 0; i<nprec; i++) {
			scan = Newpts->first;
			for (int j = 0; j<nprec; j++) {
				if (A[i][j]<mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
					isso[1] = scan;
				}
				scan = scan->next;
			}
		}												// 'mi' is the 'remaining' min distance^2, which is also stored in A[issoc[0]][issoc[1]] 
														// cprec[issoc[0]]->last and isso[1] point to that pair
	}					
	delete Newpts;										// as 'while(nprec)' terminates, 
														// all 3 (or 5) points in Newpts have been dropped and appended to the end of 3 (or 5) cprec[] curves, 
														// according to: first connect the closest pair in 3*3 (or 5*5) pairs, 
														//               then  connect the closest pair in remaining 2*2 (or 4*4) pairs, ...
														// (which is different from algorithm used in Assignment Problem)
														// and Newpts->length == 0, thus can be deleted directly

#ifdef _PRINT_ERRORS
	printf("\nN");
#endif
														// 1. Sols->length=5, but nprec=3, nfoll=3,     npres=3 can occur (if-else or else-if-else)
														// 2.                                       but npres=5 can occur
														// 3. Sols->length=5,     nprec=3, nfoll=5,     npres=5 can occur (if-if), but the two 'foll' are dropped from Sols
														// 4.                                           npres=3 can occur 		 , but the two 'foll' are dropped from Sols
														// 5. Sols->length=5,     nprec=5, nfoll=3,     npres=5 can occur (else-if-if)
														// 6.                                           npres=3 can occur 
														// 7. Sols->length=5,     nprec=5, nfoll=5,     npres=5 can occur (else-else)
														// 8.                                       but npres=3 can occur			 
														// 9. Sols->length=3,     nprec=3, nfoll=3,     npres=5 can occur 
														//10.                                           npres=3 can occur

														// status of each variable till now: 
														// 
														// in if(nprec<npres) case, i.e. case 2,3,9, 
														// 	now nprec=(from3to)0, npres(still)=5, Newpts->length=(from5to3to)0, 
														// 									   but cpres[3,4] are already meaningful, point to the two new created image contours
														//  cprec[] and A[][] are meaningless, but cpres[0,1,2] store the original cprec[0,1,2] which point to the three normal image contours
														// 
														// in if(nprec>npres) case, i.e. case 6,8, 
														// 	finally nprec=(from5to3to)0, npres=3, 
														//	cprec[0,1,2,3,4] are meaningless, but cpres[0,1,2] store the original cprec[0,1,2] which point to the three normal image contours
														// 
														// in remaining cases, i.e. case 1,4,10
														// 	nprec=(from3to)0, npres=3, 
														// 	cprec[0,1,2] are meaningless, but cpres[0,1,2] store the original cprec[0,1,2]			
														//                     i.e. case 5,7
														// 	nprec=(from5to)0, npres=5, 
														// 	cprec[0,1,2,3,4] are meaningless, but cpres[0,1,2,3,4] store the original cprec[0,1,2,3,4]
	// following images
	if (nfoll) {
		// Case of creating new images
		if (npres<nfoll) {								// previous case 4,8
			mi = 1.e100;
			for (int i = 0; i<nfoll - 1; i++) {			// i=0,1,2,3
				for (int j = i + 1; j<nfoll; j++) {		// i=0, j=1,2,3,4;   i=1, j=2,3,4;   i=2, j=3,4;   i=3, j=4
					cmp = *(cfoll[i]->first) - *(cfoll[j]->first);
					if (cmp<mi) {
						mi = cmp;
						issoc[0] = i;
						issoc[1] = j;
					}
				}
			}											// find the min distance^2 pair in all C_5^2=10 possible cfoll[i]->first pairs
														// 'mi' is the min distance^2, cfoll[issoc[0]]->first and cfoll[issoc[1]]->first point to that pair
														// this pair is regarded corresponding to the two created images
			cfoll[issoc[0]]->partneratstart = cfoll[issoc[1]];
			cfoll[issoc[1]]->partneratstart = cfoll[issoc[0]];	// make two curves become each other's 'partner_at_start'

			scan = cfoll[issoc[0]]->first;
			scan2 = cfoll[issoc[1]]->first;

			cmp2 = fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
			cmp = sqrt(mi / cmp2);						// delta theta ~, equation (29) in Bozza, 2010
			cmp_2 = cmp*cmp;
			mi = cmp_2*cmp *0.04166666666666667;		// 1/24 * (delta theta ~)^3
			parab1 = (scan->ds - scan2->ds)*mi;			// Ac(p1), here assume: scan as +, scan2 as -
			parab2 = -0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
														// - 1/12 * something * (delta theta ~)
														// Ac(p2)'s sign not checked
			cfoll[issoc[0]]->parabstart = 0.5 * (parab1 + parab2);
														// in if(nprec<npres) case, stored in _curve attribute: scurve->parabstart,             
														//		where scurve=new _curve(isso[0]) (and) isso[0] and isso[1] are just dropped from Newpts
														// in if(nprec>npres) case, stored in _point attribute: scan->parab, 					
														// 		where scan=cprec[issoc[0]]->last (and) issoc[0] and issoc[1] are the closest pair among C_5^2=10 possible cprec[i]->last pairs
														// in nprec==npres case (essentially), stored in _point attribute: scan->parab, 
														// 		where scan=cprec[issoc[0]]->last (and) cprec[issoc[0]]->last and isso[1] is the closest cprec[]->last -- Newpts pair
														// 
														// in if(npres<nfoll) case, stored in _curve attribute: cfoll[issoc[0]]->parabstart, 
														//		where cfoll[issoc[0]]->first and cfoll[issoc[1]]->first are the closest pair among C_5^2=10 possible cfoll[i]->first pairs
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: created image:
                if(astrometry){
		        avgwedgex1=(scan->x1*scan->ds - scan2->x1*scan2->ds)*mi;
		        avgwedgex2=(scan->x2*scan->ds - scan2->x2*scan2->ds)*mi;
			dx2=-(-scan->d.im+scan2->d.im);
		        d2x2=dx2*dx2;
		        dx1=-(-scan->d.re+scan2->d.re);
		        d2x1=dx1*dx1;
		        cfoll[issoc[0]]->parabastrox1 =-0.125*d2x1*dx2*mi-avgwedgex1; 
		        cfoll[issoc[0]]->parabastrox2 =-0.125*d2x2*dx1*mi+avgwedgex2;
                }
#ifdef _PRINT_ERRORS
			printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, cfoll[issoc[0]]->parabstart, (scan->ds + scan2->ds)*mi / 2, fabs(cfoll[issoc[0]]->parabstart)*(cmp*cmp) / 10, 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) - 2 * cmp*cmp2)*cmp);
#endif
			mi = fabs((parab1-parab2) *0.5) + fabs(cfoll[issoc[0]]->parabstart)*(cmp*cmp *0.1) + 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) - 2.0 * cmp*cmp2)*cmp;
														// mi = E1new(c) + E3(c) + E2(c), have checked sign
#ifdef _noparab
			mi = fabs(cfoll[issoc[0]]->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
			cfoll[issoc[0]]->parabstart = 0.;
#endif
#ifdef _selectimage
			if (_selectionimage)
#endif
			pref=(scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5;// trapezium approximation, have checked sign
																	// Ac(t)
			theta->Mag -= ((scan->dJ>0) ? -1 : 1)*(pref + cfoll[issoc[0]]->parabstart);
														// according to the assumption (scan as +, scan2 as -), 
														// 	scan->dJ > 0 thus Mag -= -(Ac(t)+Ac(p)), thus correct sign; 
														// if finally scan as -, scan2 as +, 
														// 	scan->dJ < 0 thus Mag -= (-A(t)-A(p)), thus still correct sign
#ifdef _ERRORS_ANALYTIC
			char filnam[32];
			sprintf(filnam, "%02dfoll.txt", NPS);
			FILE* f = fopen(filnam, "a+");
			fprintf(f, "%.15le %.15le %.15le\n", -((scan->dJ > 0) ? -1 : 1)* (pref), -((scan->dJ > 0) ? -1 : 1)* (cfoll[issoc[0]]->parabstart), mi);
			fclose(f);
#endif
            if(astrometry){
				dx2=scan2->x2 - scan->x2;
		      	avgx1=scan->x1 + scan2->x1;
				avg2x1=avgx1*avgx1;
				avgx2=scan->x2 + scan2->x2;
				theta->astrox1 += ((scan->dJ>0) ? -1 : 1)*(avg2x1*dx2*0.125+cfoll[issoc[0]]->parabastrox1);
				theta->astrox2 -= ((scan->dJ>0) ? -1 : 1)*(pref*avgx2*0.25+cfoll[issoc[0]]->parabastrox2);
			}
			theta->maxerr += mi;						// add the crossing-critical-curve term to theta
                        
			cfoll[issoc[1]]->parabstart = -cfoll[issoc[0]]->parabstart; 
			if(astrometry){
				cfoll[issoc[1]]->parabastrox2=-cfoll[issoc[0]]->parabastrox2;
				cfoll[issoc[1]]->parabastrox1=-cfoll[issoc[0]]->parabastrox1;
			}

			Sols->append(cfoll[issoc[0]]);				
			Sols->append(cfoll[issoc[1]]);				// in case 4, the two 'foll' curves are previously dropped from Sols
														// now append them back to Sols with no other thing happens
														// (see comments just in front of if(nfoll))
														// in case 8, two new 'foll' curves are appended to Sols, causing Sols->length=(from5to)7
			nfoll -= 2;
			ij = 0;
			for (int i = 0; i<nfoll; i++) {
				if (i == issoc[0]) ij++;
				if (i == issoc[1] - 1) ij++;
				cfoll[i] = cfoll[i + ij];
			}											// nfoll from 5 to 3
														// originally, cfoll[issoc[0]]->first and cfoll[issoc[1]]->first correspond to the two created images
														// and issoc[0] and issoc[1] can be any index from 0 to 4; 
														// 
														// now, cfoll[0], cfoll[1], cfoll[2] point to the three normal image contours, 
														// while cfoll[3] and cfoll[4] are meaningless
		}



														// 1. Sols->length=5, but nprec=3, nfoll=3,     npres=3 can occur (if-else or else-if-else)
														// 2.                                       but npres=5 can occur
														// 3. Sols->length=5,     nprec=3, nfoll=5,     npres=5 can occur (if-if), but the two 'foll' are dropped from Sols
														// 4.                                           npres=3 can occur 		 , but the two 'foll' are dropped from Sols
														// 5. Sols->length=5,     nprec=5, nfoll=3,     npres=5 can occur (else-if-if)
														// 6.                                           npres=3 can occur 
														// 7. Sols->length=5,     nprec=5, nfoll=5,     npres=5 can occur (else-else)
														// 8.                                       but npres=3 can occur			 
														// 9. Sols->length=3,     nprec=3, nfoll=3,     npres=5 can occur 
														//10.                                           npres=3 can occur

														// status of each variable till now: 
														// 
														// in if(nprec<npres) case, i.e. case 2,3,9, 
														// 	now nprec=(from3to)0, npres(still)=5, Newpts->length=(from5to3to)0, 
														// 									   but cpres[3,4] are already meaningful, point to the two new created image contours
														//  cprec[] and A[][] are meaningless, but cpres[0,1,2] store the original cprec[0,1,2] which point to the three normal image contours
														// 
														// in if(nprec>npres) case, i.e. case 6,8, 
														// 	finally nprec=(from5to3to)0, npres=3, 
														//	cprec[0,1,2,3,4] are meaningless, but cpres[0,1,2] store the original cprec[0,1,2] which point to the three normal image contours
														// 
														// in remaining cases, i.e. case 1,4,10
														// 	nprec=(from3to)0, npres=3, 
														// 	cprec[0,1,2] are meaningless, but cpres[0,1,2] store the original cprec[0,1,2]			
														//                     i.e. case 5,7
														// 	nprec=(from5to)0, npres=5, 
														// 	cprec[0,1,2,3,4] are meaningless, but cpres[0,1,2,3,4] store the original cprec[0,1,2,3,4]
														//
														// then
														//
														// in if(npres<nfoll) case, i.e. case 4,8, 
														// 	finally nfoll=(from5to)3, npres=3, 
														//	cfoll[3,4] are meaningless, but cfoll[0,1,2] point to the three normal image contours, 
														// 	in case 4, two 'foll' curves (previously dropped from Sols) are now appended back to Sols with no other thing happens
		// case of image destruction 
		if (npres>nfoll) {								// previous case 2,5,9
			mi = 1.e100;
			for (int i = 0; i<npres - 1; i++) {			// i=0,1,2,3
				for (int j = i + 1; j<npres; j++) {		// i=0, j=1,2,3,4;   i=1, j=2,3,4;   i=2, j=3,4;   i=3, j=4
					cmp = *(cpres[i]->last) - *(cpres[j]->last);
					if (cmp<mi) {
						mi = cmp;
						issoc[0] = i;
						issoc[1] = j;
					}
				}
			}											// find the min distance^2 pair in all C_5^2=10 possible cpres[i]->last pairs
														// 'mi' is the min distance^2, cpres[issoc[0]]->last and cpres[issoc[1]]->last point to that pair
														// this pair is regarded corresponding to the two destructed images
														// (for case 2,9, no need to do this, because already know cpres[3,4] correspond to the two destructed images;
														// but for case 5, this procedure is necessary)
			cpres[issoc[0]]->partneratend = cpres[issoc[1]];
			cpres[issoc[1]]->partneratend = cpres[issoc[0]];// make two curves become each other's 'partner_at_end'

			scan = cpres[issoc[0]]->last;
			scan2 = cpres[issoc[1]]->last;

			cmp2 = fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
			cmp = sqrt(mi / cmp2);						// delta theta ~, equation (29) in Bozza, 2010
			cmp_2 = cmp*cmp;
			mi = cmp_2*cmp *0.0416666666667;			// 1/24 * (delta theta ~)^3
			parab1 = -(scan->ds - scan2->ds)*mi;
			parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
			scan->parab = 0.5 * (parab1 + parab2);		// in if(nprec<npres) case, stored in _curve attribute: scurve->parabstart,             
														//		where scurve=new _curve(isso[0]) (and) isso[0] and isso[1] are just dropped from Newpts
														// in if(nprec>npres) case, stored in _point attribute: scan->parab, 					
														// 		where scan=cprec[issoc[0]]->last (and) issoc[0] and issoc[1] are the closest pair among C_5^2=10 possible cprec[i]->last pairs
														// in nprec==npres case (essentially), stored in _point attribute: scan->parab, 
														// 		where scan=cprec[issoc[0]]->last (and) cprec[issoc[0]]->last and isso[1] is the closest cprec[]->last -- Newpts pair
														// 
														// in if(npres<nfoll) case, stored in _curve attribute: cfoll[issoc[0]]->parabstart, 
														//		where cfoll[issoc[0]]->first and cfoll[issoc[1]]->first are the closest pair among C_5^2=10 possible cfoll[i]->first pairs
														// in if(npres>nfoll) case, stored in _point attribute: scan->parab, 					
														// 		where scan=cpres[issoc[0]]->last (and) issoc[0] and issoc[1] are the closest pair among C_5^2=10 possible cpres[i]->last pairs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: destructed image:
            if(astrometry){
		        avgwedgex1=-(scan->x1*scan->ds - scan2->x1*scan2->ds)*mi;
		        avgwedgex2=-(scan->x2*scan->ds - scan2->x2*scan2->ds)*mi;
				dx2=-(scan->d.im-scan2->d.im);
		        d2x2=dx2*dx2;
		        dx1=-(scan->d.re-scan2->d.re);
		        d2x1=dx1*dx1; 
		        scan->parabastrox1 =-0.125*d2x1*dx2*mi-avgwedgex1;
		        scan->parabastrox2 =-0.125*d2x2*dx1*mi+avgwedgex2;
            }

#ifdef _PRINT_ERRORS
			printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds + scan2->ds)*mi / 2, fabs(scan->parab)*(cmp*cmp) / 10, 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) + 2 * cmp*cmp2)*cmp);
#endif

			mi = fabs((parab1-parab2) *0.5) + fabs(scan->parab)*(cmp*cmp *0.1) + 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) + 2.0 * cmp*cmp2)*cmp;
#ifdef _noparab
			mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
			scan->parab = 0.;
#endif
#ifdef _selectimage
			if (_selectionimage)
#endif
			pref=(scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5;
			theta->Mag += ((scan->dJ>0) ? -1 : 1)*( pref + scan->parab);
#ifdef _ERRORS_ANALYTIC
			char filnam[32];
			sprintf(filnam, "%02dfoll.txt", NPS);
			FILE* f = fopen(filnam, "a+");
			fprintf(f, "%.15le %.15le %.15le\n", ((scan->dJ > 0) ? -1 : 1)* (pref), ((scan->dJ > 0) ? -1 : 1)* (scan->parab), mi);
			fclose(f);
#endif
			if(astrometry){
				dx2=scan2->x2 - scan->x2;
				avgx1=scan->x1 + scan2->x1;
				avg2x1=avgx1*avgx1;
				avgx2=scan->x2 + scan2->x2;
				theta->astrox1 -= ((scan->dJ>0) ? -1 : 1)*(avg2x1*dx2*0.125+scan->parabastrox1);
				theta->astrox2 += ((scan->dJ>0) ? -1 : 1)*(pref*avgx2*0.25+scan->parabastrox2);
			}
			theta->maxerr += mi;						// add the crossing-critical-curve term to theta
			scan2->parab = -scan->parab;
			if(astrometry){
				scan2->parabastrox2=-scan->parabastrox2;
				scan2->parabastrox1=-scan->parabastrox1;                        
			}

			npres -= 2;
			ij = 0;
			for (int i = 0; i<npres; i++) {
				if (i == issoc[0]) ij++;
				if (i == issoc[1] - 1) ij++;
				cpres[i] = cpres[i + ij];
			}											// npres from 5 to 3
														// originally, cpres[issoc[0]]->last and cpres[issoc[1]]->last correspond to the two distructed images
														// and issoc[0] and issoc[1] can be any index from 0 to 4 (in case 5); 
														// 
														// now, cpres[0], cpres[1], cpres[2] point to the three normal image contours, 
														// while cpres[3] and cpres[4] are meaningless
		}





														// 1. Sols->length=5, but nprec=3, nfoll=3,     npres=3 can occur (if-else or else-if-else)
														// 2.                                       but npres=5 can occur
														// 3. Sols->length=5,     nprec=3, nfoll=5,     npres=5 can occur (if-if), but the two 'foll' are dropped from Sols
														// 4.                                           npres=3 can occur 		 , but the two 'foll' are dropped from Sols
														// 5. Sols->length=5,     nprec=5, nfoll=3,     npres=5 can occur (else-if-if)
														// 6.                                           npres=3 can occur 
														// 7. Sols->length=5,     nprec=5, nfoll=5,     npres=5 can occur (else-else)
														// 8.                                       but npres=3 can occur			 
														// 9. Sols->length=3,     nprec=3, nfoll=3,     npres=5 can occur 
														//10.                                           npres=3 can occur

														// status of each variable till now: 
														// 
														// in if(nprec<npres) case, i.e. case 2,3,9, 
														// 	now nprec=(from3to)0, npres(still)=5, Newpts->length=(from5to3to)0, 
														// 									   but cpres[3,4] are already meaningful, point to the two new created image contours
														//  cprec[] and A[][] are meaningless, but cpres[0,1,2] store the original cprec[0,1,2] which point to the three normal image contours
														// 
														// in if(nprec>npres) case, i.e. case 6,8, 
														// 	finally nprec=(from5to3to)0, npres=3, 
														//	cprec[0,1,2,3,4] are meaningless, but cpres[0,1,2] store the original cprec[0,1,2] which point to the three normal image contours
														// 
														// in remaining cases, i.e. case 1,4,10
														// 	nprec=(from3to)0, npres=3, 
														// 	cprec[0,1,2] are meaningless, but cpres[0,1,2] store the original cprec[0,1,2]			
														//                     i.e. case 5,7
														// 	nprec=(from5to)0, npres=5, 
														// 	cprec[0,1,2,3,4] are meaningless, but cpres[0,1,2,3,4] store the original cprec[0,1,2,3,4]
														//
														// then
														//
														// in if(npres<nfoll) case, i.e. case 4,8, 
														// 	finally nfoll=(from5to)3, npres=3, 
														//	cfoll[3,4] are meaningless, but cfoll[0,1,2] point to the three normal image contours, 
														// 	in case 4, two 'foll' curves (previously dropped from Sols) are now appended back to Sols with no other thing happens
														// 	in case 8, two new 'foll' curves are appended to Sols, causing Sols->length=(from5to)7
														// 
														// in if(npres>nfoll) case, i.e. case 2,5,9, 
														// 	finally npres=(from5to)3, nfoll=3,
														// cpres[3,4] are meaningless, but cpres[0,1,2] point to the three normal image contours	
														//
														// in remaining cases, i.e. case 1,6,10
														// 	npres=3, nfoll=3, 
														// 	cpres[0,1,2] store the original cprec[0,1,2] which point to the three normal image contours			
														//                     i.e. case 3,7
														// 	npres=5, nfoll=5, 
														//	cpres[0,1,2,3,4] store the original cprec[0,1,2,3,4] in case 7
														//  	or  
														//  cpres[0,1,2] point to the three normal image contours + cpres[3,4] point to the two new created image contours 									
		// constructing distance matrix with subsequent images
 
		mi = 1.e100;
		for (int i = 0; i<npres; i++) {	// loop through cpres[i]->last
			for (int j = 0; j<npres; j++) {	// loop through cfoll[j]->first
				A[i][j] = signbit(cpres[i]->last->dJ) == signbit(cfoll[j]->first->dJ)? *(cpres[i]->last) - *(cfoll[j]->first) : 100;
				if (A[i][j]<mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}	
		}												// construct 3*3: case 2,9,6,8,1,4,10,5 
														//        or 5*5: case 3,7
														// distance matrix between cpres[]->last (as i) and cfoll[]->first (as j)
														// matrix value is distance^2 for same parity, 100 for different parity
														// 
														// 'mi' is the min distance^2, which is also stored in A[issoc[0]][issoc[1]] 
														// cpres[issoc[0]]->last and cfoll[issoc[1]]->first point to that pair

		// association with the images below
		while (npres) {
			scan = cpres[issoc[0]]->last;
			scan2 = cfoll[issoc[1]]->first;				// start from the closest pair

			cmp2 = mi / fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
			cmp = (scan->theta->th - scan2->theta->th);	// - delta theta
			cmp_2 = cmp*cmp;
			mi = cmp_2*cmp *0.041666666667;				// - 1/24 * (delta theta)^3
			parab1 = (scan->ds + scan2->ds)*mi; // Old parabolic correction, -A(p1)
			// New parabolic correction
			parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im - scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re - scan->d.re)) * cmp;
														// - 1/12 * something * (delta theta) = -A(p2), have checked something's sign

			scan->parab = 0.5*(parab1+parab2);	// -A(p)// in if(nprec<npres) case, stored in _curve attribute: scurve->parabstart,             
														//		where scurve=new _curve(isso[0]) (and) isso[0] and isso[1] are just dropped from Newpts
														// in if(nprec>npres) case, stored in _point attribute: scan->parab, 					
														// 		where scan=cprec[issoc[0]]->last (and) issoc[0] and issoc[1] are the closest pair among C_5^2=10 possible cprec[i]->last pairs
														// in nprec==npres case (essentially), stored in _point attribute: scan->parab, 
														// 		where scan=cprec[issoc[0]]->last (and) cprec[issoc[0]]->last and isso[1] is the closest cprec[]->last -- Newpts pair
														// 
														// in if(npres<nfoll) case, stored in _curve attribute: cfoll[issoc[0]]->parabstart, 
														//		where cfoll[issoc[0]]->first and cfoll[issoc[1]]->first are the closest pair among C_5^2=10 possible cfoll[i]->first pairs
														// in if(npres>nfoll) case, stored in _point attribute: scan->parab, 					
														// 		where scan=cpres[issoc[0]]->last (and) issoc[0] and issoc[1] are the closest pair among C_5^2=10 possible cpres[i]->last pairs
														// in npres==nfoll case, stored in _point attribute: scan->parab, 
														// 		where scan=cpres[issoc[0]]->last (and) cpres[issoc[0]]->last and cfoll[issoc[1]]->first is the closest cpres[]->last -- cfoll[]->first pair
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: ordinary image:
            if(astrometry){
		        avgwedgex1=(scan->x1*scan->ds + scan2->x1*scan2->ds)*mi;
		        avgwedgex2=(scan->x2*scan->ds + scan2->x2*scan2->ds)*mi;
				dx2=scan->d.im+scan2->d.im;
		        d2x2=dx2*dx2;
		        dx1=scan->d.re+scan2->d.re;
		        d2x1=dx1*dx1; 
		        scan->parabastrox1 =-0.125*d2x1*dx2*mi-avgwedgex1;
		        scan->parabastrox2 =-0.125*d2x2*dx1*mi+avgwedgex2;
            }    

#ifdef _PRINT_ERRORS
			printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds - scan2->ds)*mi / 2, fabs(scan->parab)*(cmp2) / 10, fabs(scan->parab)*(1.5*fabs(cmp2 / (cmp*cmp) - 1)));
#endif
			/******************************************* changed *******************************************/
			// Question: typo? should be cmp_2 instead of cmp2
			// original: 
			//mi = fabs((parab1-parab2) *0.5) + fabs( scan->parab * (cmp2 *0.1 + 1.5*fabs(cmp2 / (cmp_2) - 1)) );
			// should be: 
			mi = fabs((parab1-parab2) *0.5) + fabs( scan->parab * (cmp_2*0.1 + 1.5*fabs(cmp2 / (cmp_2) - 1)) );
															// mi = E1new + E3 + E2
			/*******************************************   end   *******************************************/
#ifdef _noparab
			mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
			scan->parab = 0.;
#endif
#ifdef _selectimage
			if (_selectionimage)
#endif						
            pref=(scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5;	// trapezium approximation, -A(t)
			theta->Mag += ((scan->dJ>0) ? -1 : 1)*(pref + scan->parab);	// dJ>0 thus Mag += -(-A(t)-A(p)), thus correct sign
																		// add this term to theta
			if(astrometry){
				dx2=scan2->x2 - scan->x2;
				avgx1=scan->x1 + scan2->x1;
				avg2x1 = avgx1 * avgx1;
				avgx2=scan->x2 + scan2->x2;
				theta->astrox1 -= ((scan->dJ>0) ? -1 : 1)*(avg2x1*dx2*0.125+scan->parabastrox1);
				theta->astrox2 += ((scan->dJ>0) ? -1 : 1)*(pref*avgx2*0.25+scan->parabastrox2);
			}
			theta->maxerr += mi;										// add this term to theta
#ifdef _ERRORS_ANALYTIC
				char filnam[32];
				sprintf(filnam, "%02dfoll.txt", NPS);
				FILE* f = fopen(filnam, "a+");
				fprintf(f, "%.15le %.15le %.15le\n", ((scan->dJ > 0) ? -1 : 1)* (pref), ((scan->dJ > 0) ? -1 : 1)* (scan->parab), mi);
				fclose(f);
#endif
			cpres[issoc[0]]->join(cfoll[issoc[1]]);		// join 'cfoll' curve to 'cpres' curve

			npres--;
			for (int i = issoc[0]; i<npres; i++) {
				cpres[i] = cpres[i + 1];
				for (int j = 0; j<npres + 1; j++) {
					A[i][j] = A[i + 1][j];
				}
			}											// npres from 3 to 2 (or 5 to 4)
														// drop cpres[issoc[0]], and now cpres[0], cpres[1] are meaningful while cpres[2] is meaningless
														// 							 (or cpres[0-3])						 (or cpres[4])				
														// drop the corresponding row in A[][], now A[][]'s meaningful range is 2*3 (or 4*5)

			for (int j = issoc[1]; j<npres; j++) {
				cfoll[j] = cfoll[j + 1];
				for (int i = 0; i<npres; i++) {
					A[i][j] = A[i][j + 1];
				}
			}											// drop cfoll[issoc[1]], and now cfoll[0], cfoll[1] are meaningful while cfoll[2] is meaningless
														// 							 (or cfoll[0-3])						 (or cfoll[4])			
														// drop the corresponding column in A[][], now A[][]'s meaningful range is 2*2 (or 4*4)

			mi = 1.e100;
			for (int i = 0; i<npres; i++) {
				for (int j = 0; j<npres; j++) {
					if (A[i][j]<mi) {
						mi = A[i][j];
						issoc[0] = i;
						issoc[1] = j;
					}
				}
			}
		}
	}

}

//////////////////////////////
//////////////////////////////
////////_point methods
//////////////////////////////
//////////////////////////////

/******************************************* changed *******************************************/
// inline added
inline _point::_point(double x, double y, _theta *theta1) {
	x1 = x;
	x2 = y;
	theta = theta1;
	next = 0 ;
	prev = 0 ;
	for (int i=0; i<(max_skiplist_level+1); i++)
	{
		next_array[i] = 0 ;
	}
}

inline double _point::operator-(_point p2) {
	return (x1 - p2.x1)*(x1 - p2.x1) + (x2 - p2.x2)*(x2 - p2.x2);
}
/*******************************************   end   *******************************************/


//////////////////////////////
//////////////////////////////
////////_curve methods
//////////////////////////////
//////////////////////////////

inline _curve::_curve(void) {
	length = 0;
	first = last = 0;
	partneratstart = partneratend = 0;
}

inline _curve::_curve(_point *p1) {
	length = 1;
	first = last = p1;
	p1->prev = p1->next = 0;
	partneratstart = partneratend = 0;
}

inline _curve::~_curve(void) {
	_point *scan1, *scan2;
	scan1 = first;
	for (int i = 0; i<length; i++) {
		scan2 = scan1->next;
		delete scan1;
		scan1 = scan2;
	}
}

/******************************************* changed *******************************************/
_curve *_curve::divide(_point *ref, int length2) {	// length2 is Nmember from ref->next to last (excluding the one pointed by 'ref')
	_point *scan;
	_curve *nc;
	//int l1;

	//l1 = 1;
	//for (scan = first; scan != ref; scan = scan->next) l1++; // loop from first, if ref==first, then l1=1;
															 // l1 is Nmember from first to ref(including the one pointed by 'ref')

	nc = new _curve();		// create another object of class _curve using 'new' keyword, 
	                        // memory of the object is allocated on heap and the pointer to that object is assigned to 'nc'(new curve), 
							// constructor is called (create a _curve class variable without any member).
							// (as the object is on heap, we can return it's pointer;
							// if the object is on stack, the returned pointer will be a dangling pointer) 

	nc->first = ref->next;  // 'nc'(new curve) consists of ref->next to original curve's last, length = original_length - l1
	nc->first->prev = 0;
	nc->last = last;
	//nc->length = length - l1;
	nc->length = length2 ;
	nc->partneratend = partneratend; 
	if (partneratend) partneratend->partneratend = nc; // if partner_at_end really points to a _curve variable, 
													   // then redirect *partner_at_end's attribute 'partner_at_end' from original curve to 'nc'
													   // (seems can direct endlessly, like partner_at_end -> partner_at_end -> ... partner_at_end -> ...)

							// change original curve 
	//length = l1;			// original curve now consists of first to ref, length = l1
	length = length - length2 ;
	last = ref;
	ref->next = 0;
	partneratend = 0;       // set from xxx to NULL
	return nc;
}
/*******************************************   end   *******************************************/

void _curve::append(double x1, double x2) { // method: create a new _point variable on heap, 
											//         and then append it to the end of the linked list (object of class _curve)
											// has no return value
	_point *pp;
	pp = new _point(x1, x2, 0); // create an object of class _point on heap, the pointer to that object is assigned to 'pp'
								// constructor is called, assign value to attributes x1, x2, 
								// and set attribute theta(pointer) = NULL
								// (note each _point variable corresponds to a _theta variable, but now has not pointed to its _theta variable)
	if (length == 0) {
		first = pp;
		last = pp;
		pp->prev = 0;
	}
	else {
		last->next = pp;
		pp->prev = last;
		last = pp;
	}
	pp->next = 0;
	length++;
}

void _curve::append(_point *pp) { // method: append an existing _point variable (pointed by input pointer 'pp') 
								  //         to the end of the linked list (object of class _curve)
								  //         (the linked list must already have at least one element)

	pp->next = last->next; // set to NULL
	pp->prev = last;
	last->next = pp;
	last = pp;
	length++;
}

void _curve::prepend(double x1, double x2) { // same as _curve::append, just append it to the beginning of the linked list (object of class _curve)
											 // not used
	_point *pp;
	pp = new _point(x1, x2, 0);
	if (length == 0) {
		first = pp;
		last = pp;
		pp->next = 0;
	}
	else {
		first->prev = pp;
		pp->next = first;
		first = pp;
	}
	pp->prev = 0;
	length++;
}

_curve *_curve::join(_curve *nc) {	     // method: used once in OrderImages (actually also used one time in PlotCrit
										 // input is a pointer to 'nc'(new curve) variable
										 // join all elements of 'nc' at the end of current curve
										 // return the pointer to current curve

	if (length>0) {						 // deal with four conditions: length > 0 && nc->length > 0, needs to change 'last' and two edge elements' 'next'/'prev'
										 //                            length ==0 && nc->length > 0, needs to change 'first'/'last'
										 //                            length > 0 && nc->length ==0, needs to change nothing
										 //                            length ==0 && nc->length ==0, needs to change nothing
		last->next = nc->first;
	}
	else {
		first = nc->first;
	};
	if (nc->length>0) {
		nc->first->prev = last;
		last = nc->last;
	}
	length += nc->length;

	partneratend = nc->partneratend;      				 // set current curve's attribute 'partner_at_end' to nc's
	if (partneratend) partneratend->partneratend = this; // redirect nc->partner_at_end's attribute 'partner_at_end' from 'nc' to current curve
	// 'this' is a pointer pointed to current curve
	//
	// Suppose that you create an object named x of class A, 
	// and class A has a nonstatic member function f(). 
	// If you call the function x.f(), the keyword 'this' in the body of f() stores the address of x.

	nc->first = 0;                         
	nc->last = 0;
	nc->length = 0;
	delete nc;							  // Although original 'nc''s elements have been linked to current curve, 
										  // these elements still can be accessed by 'nc'. 
									      // So directly delete nc will cause the deletion of these elements. 
										  // So we need to break the relation between 'nc' and these elements before delete nc. 
										  // Note: destructor is called through delete keyword
	return this;
}

_curve *_curve::joinbefore(_curve *nc) { // method: not used
	if (length>0) {
		first->prev = nc->last;
	}
	else {
		last = nc->last;
	};
	if (nc->length>0) {
		nc->last->next = first;
		first = nc->first;
	}
	length += nc->length;
	nc->first = 0;
	nc->last = 0;
	nc->length = 0;
	delete nc;
	return this;
}

_curve *_curve::reverse(void) {		     // method: not used
	_point *scan1, *scan2, *scambio;
	if (length>1) {
		scan1 = first;
		while (scan1) {
			scan2 = scan1->next;
			scambio = scan1->next;
			scan1->next = scan1->prev;
			scan1->prev = scambio;
			scan1 = scan2;
		}
		scambio = first;
		first = last;
		last = scambio;
	}
	return this;
}

void _curve::drop(_point *ref) {         // method: input a pointer(ref) to an element of current curve
									     //         drop that element out of the current curve, 
										 //         but the element isn't deleted, still can be accessed by 'ref'
										 //         no return valve
										 // (if length==0 or 'ref' doesn't point to an element, then nothing happens)
										 // (use for loop to check if 'ref' really points to an element, thus O(N), but loop is from 'last' to 'first')
	_point *scan;
	if (length) {
		for (scan = last; scan && (scan != ref); scan = scan->prev); // if 'ref' really points to an element, then scan == ref
															         // 									  else scan == first->prev == NULL

															// still four conditions: ref == first and ref == last, just length==1 case
															//						  ref == first and ref != last, reset 'first' and ref->next's 'prev'
															//						  ref != first and ref == last
															//						  ref != first and ref != last, reset ref->prev's 'next' and ref->next's 'prev'
		if (scan) {
			if (length == 1) {
				first = last = 0;
			}
			else {
				if (ref->prev) {         // if ref != first
					ref->prev->next = ref->next;
					if (ref == last) {
						last = ref->prev;
					}
				}
				if (ref->next) {         // if ref != last
					ref->next->prev = ref->prev;
					if (ref == first) {
						first = ref->next;
					}
				}
			}
			length--;
		}
	}
}

double _curve::closest2(_point *ref, _point **clos2) { // method: not used
	double mi = 1.e100, mi2 = 1.e100, FP;
	_point *scan, *clos;
	if (length>1) {
		clos = *clos2 = first;
		for (scan = first; scan != 0; scan = scan->next) {
			FP = *scan - *ref;
			if (FP<mi) {
				mi2 = mi;
				mi = FP;
				*clos2 = clos;
				clos = scan;
			}
			else if (FP<mi2) {
				mi2 = FP;
				*clos2 = scan;
			}
		}
	}
	else {
		*clos2 = 0;
	}
	return (**clos2 - *ref);
}

double _curve::closest(_point *ref, _point **clos) { // method: not used (actually used one time in PlotCrit)
	double mi = 1.e100, FP;
	_point *scan;
	for (scan = first; scan != 0; scan = scan->next) {
		FP = *scan - *ref;
		if (FP<mi) {
			mi = FP;
			*clos = scan;
		}
	}
	return mi;
}

void _curve::complement(_point **sott, int lensott, _point **res, int lenres) { // method: not used
	int flag, i;
	_point *scan;
	i = 0;
	for (scan = first; scan != 0; scan = scan->next) {
		flag = 0;
		for (int j = 0; (j<lensott) && (!flag); j++) {
			if (scan == sott[j]) {
				flag = 1;
			}
		}
		if ((!flag) && (i<lenres)) {
			res[i] = scan;
			i++;
		}
	}
}

//////////////////////////////
//////////////////////////////
////////_sols methods
//////////////////////////////
//////////////////////////////


_sols::_sols(void) {
	length = 0;
	first = last = 0;
}

_sols::~_sols(void) {
	_curve *scan1, *scan2;
	scan1 = first;
	while (scan1) {
		scan2 = scan1->next;
		delete scan1;           // 'delete scan1' will call destructor of _curve class, 
								// which will delete all _point varialbes of the curve pointed by scan1
		scan1 = scan2;
	}
}

void _sols::append(_curve *cc) {// method: append an existing _curve variable (pointed by input pointer 'cc') 
								//         to the end of the linked list (object of class _sols)
	if (length == 0) {
		first = cc;
		last = cc;
		cc->prev = 0;
	}
	else {
		last->next = cc;
		cc->prev = last;
		last = cc;
	}
	cc->next = 0;
	length++;
}

void _sols::prepend(_curve *cc) { // not used
	if (length == 0) {
		first = cc;
		last = cc;
		cc->next = 0;
	}
	else {
		first->prev = cc;
		cc->next = first;
		first = cc;
	}
	cc->prev = 0;
	length++;
}

void _sols::drop(_curve *ref) { // method: used only one time in OrderImages (near the while loop and divide)
								//
								//         input a pointer(ref) to an element/curve of current sols
								//         drop that element/curve out of the current sols, 
								//         but the element/curve isn't deleted, still can be accessed by 'ref'
								//         no return valve
								// (if length==0 or 'ref' doesn't point to an element/curve, then nothing happens)
								// (use for loop to check if 'ref' really points to an element, thus O(N), but loop is from 'last' to 'first')
								//
								// totally same to   void _curve::drop(_point *ref)
	_curve *scan;
	if (length) {
		for (scan = last; scan && (scan != ref); scan = scan->prev);
		if (scan) {
			if (length == 1) {
				first = last = 0;
			}
			else {
				if (ref->prev) {
					ref->prev->next = ref->next;
					if (ref == last) {
						last = ref->prev;
					}
				}
				if (ref->next) {
					ref->next->prev = ref->prev;
					if (ref == first) {
						first = ref->next;
					}
				}
			}
			length--;
		}
	}
}

void _sols::join(_sols *nc) {                   // method: not used
	if (length>0) {
		last->next = nc->first;
	}
	else {
		first = nc->first;
	};
	if (nc->length>0) {
		nc->first->prev = last;
		last = nc->last;
	}
	length += nc->length;
	nc->first = 0;
	nc->last = 0;
	nc->length = 0;
	delete nc;
}

//////////////////////////////
//////////////////////////////
////////_theta methods
//////////////////////////////
//////////////////////////////

_theta::_theta(double th1) {
	th = th1;
}
_thetas::_thetas(void) {
	length = 0;
}

_thetas::~_thetas(void) {
	_theta *scan, *scan2;    // pointers that point to _theta variable
	scan = first;            // start from attribute 'first', which is a pointer that point to _theta variable
	while (scan) {             // if scan != NULL
		scan2 = scan->next;  // visit the object (class _theta)'s attribute 'next' pointed by scan
		delete scan;         // delete the object (class _theta) pointed by scan (pointed by first at beginning), 
		                     // causing the delete of it's prev and next attributes
		scan = scan2;       
		                     // my guess: if scan == last, then scan2 = scan->next = NULL, then loop terminates
							 // Question: when will ~_thetas() be called? I don't see explicit call. 
	}
}

_theta *_thetas::insert(double th) { // it's return value is a pointer that points to _theta variable
	_theta *scan, *scan2;

	scan2 = new _theta(th); // create an object of class _theta using 'new' keyword, 
	                        // memory of the object is allocated on heap and the pointer to that object is assigned to scan2, 
							// constructor is called.
							// (as the object is on heap, we can return it's pointer;
							// if the object is on stack, the returned pointer will be a dangling pointer) 
	if (length) {   // if length != 0
		if (th<first->th) {           // if current 'th'  <  the object(pointed by first)'s attribute 'th'
		                              // means the newly inserted element will become the new first

			first->prev = scan2;          // change the object(pointed by first)'s attribute 'prev' from NULL to scan2
			scan2->next = first;      	  // set the newly inserted element's attribute 'next' to first(the original first)
			scan2->prev = 0;          	  // set                                        'prev' to NULL
			first = scan2;                // change the attribute 'first' to newly inserted one
		}
		else {
			if (th>last->th) {        // the newly inserted element will become the new last
				last->next = scan2;       // same as above, from NULL to scan2
				scan2->prev = last;       //                set to last(the original last)
				scan2->next = 0;          //                set to NULL
				last = scan2;             // change the attribute 'last' to newly inserted one
			}
			else {                    // *first's 'th' < current 'th' < *last's 'th', 
									  // thus the newly inserted element must be in the middle of link list
				scan = first;
				while (scan->th<th) scan = scan->next; // loop from the element pointed by 'first', 
													   // if it's 'th' < current th, then jump to the next element, 
													   // till some element's 'th' > current th (now scan points to that element)
				
				scan2->next = scan;                    // now newly inserted element should be placed just in front of scan
				scan2->prev = scan->prev;              // thus scan and scan->prev will be newly inserted element's 'next' and 'prev'

				scan->prev->next = scan2;              // scan2 will be scan->prev's 'next' and scan's 'prev'
				scan->prev = scan2;
			}
		}
	}
	else {         // if length== 0
		first = scan2;    // becomes the only element
		last = scan2;
		scan2->next = 0;  // set the newly inserted object's attributes 'next' and 'prev' to NULL
		scan2->prev = 0;
	}
	length++;             // class '_thetas's attribute 'length' +=1
	//	scan2->maxerr=0.;
	return scan2;         // return the newly inserted object's pointer (object on heap)
}

/******************************************* changed *******************************************/
_theta * _thetas::insert_at_certain_position(_theta * itheta, double th)
{
	_theta *scan2;

	scan2 = new _theta(th); // create an object of class _theta using 'new' keyword, 
	                        // memory of the object is allocated on heap and the pointer to that object is assigned to scan2, 
							// constructor is called.
							// (as the object is on heap, we can return it's pointer;
							// if the object is on stack, the returned pointer will be a dangling pointer) 
	scan2->prev = itheta ;
	scan2->next = itheta->next ;

	itheta->next->prev = scan2 ;
	itheta->next = scan2 ;

	length++ ; 

	return scan2 ;			// this method can only be used when inserting an element in the middle of linked list
							// i.e. *first's 'th' < current 'th' < *last's 'th'
							// and the new element is forced to be inserted between itheta and itheta->next, 
							// which means it's the programmer's responsibility to guarantee itheta->th < th < itheta->next->th holds
}
/*******************************************   end   *******************************************/

void _thetas::remove(_theta* stheta) { // has no return value
	_theta *scan;
	scan = first;
	while (scan!=0) {              // only last's 'next'==0
		if (scan == stheta) {          // scan == input pointer
			if (scan != first) scan->prev->next = stheta->next; // update scan's prev/next elements(if they exist)' 'next'/'prev'
			if (scan != last)  scan->next->prev = stheta->prev; // 
			delete stheta;                                      // delete the element pointed by scan(or stheta), 
		                                                        // causing the delete of it's prev and next attributes and so on
			length--;							// (if the original length==1, means both first and last, just delete the element)
			break;      										// terminates the while loop
		}
		scan = scan->next;         // when scan==last, then scan = last->next = 0, then while loop terminates
	}
								   // why not just delete stheta? why border to use while loop?
								   // Question: why not update 'first'/'last' when stheta == 'first'/'last'? (although we know we won't delete 0 and 2*PI)
}

//////////////////////////////
//////////////////////////////
////////complex methods and operators
//////////////////////////////
//////////////////////////////

/*
complex::complex(double a, double b) {
	re = a;
	im = b;
}

complex::complex(double a) {
	re = a;
	im = 0;
}

complex::complex(void) {
	re = 0;
	im = 0;
}

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
*/

//////////////////////////////
//////////////////////////////
////////_Skowron & Gould functions, translated by Tyler M. Heintz and Ava R. Hoag
//////////////////////////////
//////////////////////////////
// See copyright notice for these functions


void VBBinaryLensing::cmplx_roots_gen(complex *roots, complex *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points) {
	//roots - array which will hold all roots that had been found.
	//If the flag 'use_roots_as_starting_points' is set to
	//.true., then instead of point(0, 0) we use value from
	//this array as starting point for cmplx_laguerre

	//poly - is an array of polynomial cooefs, length = degree + 1,
	//poly[0] x ^ 0 + poly[1] x ^ 1 + poly[2] x ^ 2 + ...

	//degree - degree of the polynomial and size of 'roots' array

	//polish_roots_after - after all roots have been found by dividing
	//original polynomial by each root found,
	//you can opt in to polish all roots using full
	//polynomial

	//use_roots_as_starting_points - usually we start Laguerre's 
	//method from point(0, 0), but you can decide to use the
	//values of 'roots' array as starting point for each new
	//root that is searched for.This is useful if you have
	//very rough idea where some of the roots can be.
	//

	complex poly2[MAXM];
	static int i, j, n, iter;
	bool success;
	complex coef, prev;

	if (!use_roots_as_starting_points) {
		for (int jj = 0; jj < degree; jj++) {
			roots[jj] = complex(0, 0);
		}
	}

	for (j = 0; j <= degree; j++) poly2[j] = poly[j];

	// Don't do Laguerre's for small degree polynomials
	if (degree <= 1) {
		if (degree == 1) roots[0] = -poly[0] / poly[1];
		return;
	}

	// cmplx_laguerre2newton
	// when Laguerre, calculate criterion for each iteration, may go to SG or Newton;
	// (and the implemented Laguerre is 10% faster than zroots due to no module is calculated)
	// when SG,       calculate criterion once per 10 iterations, may go to Newton; 
	// when Newton,   calculate criterion only one time at the beginning, and do at most 10 iterations, 
	//				  if still not returned after 10 iterations, then go to Laguerre
	// claims to be 1.3 times faster than full Laguerre method

	// using 3 times root solving + 3 divisions + solve_quadratic_eq leads to 10% speedup in solving 5th degree polynomial

	// polishing usually takes only one step, thus using Laguerre or Newton does not matter the overall speed

	// cmplx_root_gen is thus 1.1*1.3*1.1 = 1.6-1.7 times faster than zroots for 5th degree polynomial

	for (n = degree; n >= 3; n--) {
		cmplx_laguerre2newton(poly2, n, &roots[n - 1], iter, success, 2);
		if (!success) {
			roots[n - 1] = complex(0, 0);
			cmplx_laguerre(poly2, n, &roots[n - 1], iter, success);
		}

		// Divide by root
		coef = poly2[n];
		for (i = n - 1; i >= 0; i--) {
			prev = poly2[i];
			poly2[i] = coef;
			coef = prev + roots[n - 1] * coef;
		}
	}


	//Find the to last 2 roots
	solve_quadratic_eq(roots[1], roots[0], poly2);
	//cmplx_laguerre2newton(poly2, 2, &roots[1], iter, success, 2);
	//if (!success) {
	//	solve_quadratic_eq(roots[1], roots[0], poly2);
	//}
	//else {
	//	roots[0] = -(roots[1] + poly2[1] / poly2[2]); // Viete's Formula for the last root
	//}



	if (polish_roots_after) {
		/******************************************* changed *******************************************/
		// no need to polish roots[4]
		//for (n = 0; n < degree; n++) {
		for (n = 0; n < degree-1; n++) {
		/*******************************************   end   *******************************************/
			cmplx_newton_spec(poly, degree, &roots[n], iter, success); // Polish roots with full polynomial
		}
	}

	return;
}

void VBBinaryLensing::solve_quadratic_eq(complex &x0, complex &x1, complex *poly) {
	complex a, b, c, b2, delta;
	a = poly[2];
	b = poly[1];
	c = poly[0];
	b2 = b*b;
	delta = sqrt(b2 - 4 * a*c);
	if (real(conj(b)*delta) >= 0) {
		x0 = -0.5*(b + delta);
	}
	else {
		x0 = -0.5*(b - delta);
	}
	if (x0 == complex(0., 0.)) {
		x1 = complex(0., 0.);
	}
	else { //Viete's formula
		x1 = c / x0;
		x0 = x0 / a;
	}
	return;

}

/*
void VBBinaryLensing::solve_cubic_eq(complex &x0, complex &x1, complex &x2, complex *poly) {
	//Cubic equation solver for comples polynomial (degree=3)
	//http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
	// poly is an array of polynomial cooefs, length = degree+1, poly[0] is constant
	//	0				1				2			3
	//poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + poly[3] x^3
	complex zeta = complex(-0.5, 0.8660254037844386);
	complex zeta2 = complex(-0.5, -0.8660254037844386);
	double third = 0.3333333333333333;
	complex s0, s1, s2;
	complex E1; //x0+x1+x2
	complex E2; //x0*x1+x1*x2+x2*x0
	complex E3; //x0*x1*x2
	complex A, B, a_1, E12, delta, A2;

	complex val, x;
	a_1 = 1 / poly[3];
	E1 = -poly[2] * a_1;
	E2 = poly[1] * a_1;
	E3 = -poly[0] * a_1;

	s0 = E1;
	E12 = E1*E1;
	A = 2.0 * E1 * E12 - 9.0 * E1 * E2 + 27.0 * E3;
	B = E12 - 3.0 * E2;
	//quadratic equation z^2 - A * z + B^3 where roots are equal to s1^3 and s2^3
	A2 = A * A;
	delta = sqrt(A2 - 4.0 * (B * B * B));
	if (real(conj(A) * delta) >= 0.0) { // scalar product to decide the sign yielding bigger magnitude
		s1 = cbrt(0.5 * (A + delta));
	}
	else
	{
		s1 = cbrt(0.5 * (A - delta));
	}
	if (s1.re == 0.0 && s1.im == 0.0) {
		s2 = complex(0, 0);
	}
	else {
		s2 = B / s1;
	}

	x0 = third * (s0 + s1 + s2);
	x1 = third * (s0 + s1 * zeta2 + s2 * zeta);
	x2 = third * (s0 + s1 * zeta + s2 * zeta2);

	return;

}
*/

void VBBinaryLensing::cmplx_laguerre(complex *poly, int degree, complex *root, int &iter, bool &success) {
	//Subroutine finds one root of a complex polynomial using
	//Laguerre's method. In every loop it calculates simplified 
	//Adams' stopping criterion for the value of the polynomial.
	//
	//Uses 'root' value as a starting point(!!!!!)
	//Remember to initialize 'root' to some initial guess or to
	//point(0, 0) if you have no prior knowledge.
	//
	//poly - is an array of polynomial cooefs
	//
	//length = degree + 1, poly(1) is constant
	//	1              2				3
	//poly(1) x ^ 0 + poly(2) x ^ 1 + poly(3) x ^ 2 + ...
	//
	//degree - a degree of the polynomial
	//
	//root - input: guess for the value of a root
	//output : a root of the polynomial
	//iter - number of iterations performed(the number of polynomial
	//evaluations and stopping criterion evaluation)
	//
	//success - is false if routine reaches maximum number of iterations
	//
	//For a summary of the method go to :
	//http://en.wikipedia.org/wiki/Laguerre's_method
	//
	static int FRAC_JUMP_EVERY = 10;
	const int FRAC_JUMP_LEN = 10;
	double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297,
		0.91577881, 0.25921289, 0.50487203,
		0.08177045, 0.13653241, 0.306162,
		0.37794326, 0.04618805, 0.75132137 }; // some random numbers

	double faq; //jump length
	double FRAC_ERR = 2.0e-15; //Fractional Error for double precision
								// 2*pi + omega + sigma = 1.9984e-15
	complex p, dp, d2p_half; //value of polynomial, 1st derivative, and 2nd derivative
	static int i, k;//j, 
	bool good_to_go;
	complex denom, denom_sqrt, dx, newroot;
	double ek, absroot, abs2p;
	complex fac_newton, fac_extra, F_half, c_one_nth;
	double one_nth, n_1_nth, two_n_div_n_1;
	complex c_one = complex(1, 0);
	complex zero = complex(0, 0);
	double stopping_crit2;

	//--------------------------------------------------------------------------------------------

	//EXTREME FAILSAFE! not usually needed but kept here just to be on the safe side. Takes care of first coefficient being 0
	if (false) {
		if (degree < 0) {
			printf("Error: cmplx_laguerre: degree<0");
			return;
		}
		if (poly[degree] == complex(0, 0)) {
			if (degree == 0) return;
			cmplx_laguerre(poly, degree - 1, root, iter, success);
		}
		if (degree <= 1) {
			if (degree == 0) {
				success = false; // we just checked if poly[0] is zero and it isnt
				printf("Warning: cmplx_laguerre: degree = 0 and poly[0] does not equal zero, no roots");
				return;
			}
			else {
				*root = -poly[0] / poly[1];
				return;
			}
		}
	} // End of EXTREME failsafe

	good_to_go = false;
	one_nth = 1.0 / degree;
	n_1_nth = (degree - 1.0)*one_nth;
	two_n_div_n_1 = 2.0 / n_1_nth;
	c_one_nth = complex(one_nth, 0.0);
	for (i = 1; i <= MAXIT; i++) {
		ek = abs(poly[degree]); 	// Preparing stopping criterion
									// |a5|
									// Question: why the (2*pi+omega)/(2*pi+omega+sigma)=3.5/4.5 is not multiplied?
									// 			 ek will thus be overestimated
		absroot = abs(*root);		// |z|
		// Calculate the values of polynomial and its first and second derivatives
		p = poly[degree];			// a5 
		dp = zero;
		d2p_half = zero;
		for (k = degree - 1; k >= 0; k--) {
			d2p_half = dp + d2p_half*(*root);
			dp = p + dp * *root;
			p = poly[k] + p*(*root); // b_k
									// a4 + a5 * z
									// a3 + z * (a4 + a5 * z)
									// ...
									// a0 + z * ()

									 //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
									 //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
									 //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
									 //Eq 8.
			ek = absroot*ek + abs(p);	// k=4, |z| * |a5| + |a4 + a5 * z|
										// k=3, |z| * (|z| * |a5| + |a4 + a5 * z|)  +  |a3 + z * (a4 + a5 * z)|
										// ...
										// k=0, |z| * ()  +  |a0 + z * ()|
		}
		iter += 1;

		abs2p = real(conj(p)*p);
		if (abs2p == 0) return;
		stopping_crit2 = pow(FRAC_ERR*ek, 2.0);
		if (abs2p < stopping_crit2) {
			//(simplified a little Eq. 10 of Adams 1967)
			//do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01*stopping_crit2) {
				return; // we are at a good place!
			}
			else {
				good_to_go = true;
			}
		}
		else {
			good_to_go = false;
		}

		faq = 1.0;
		denom = zero;
		if (dp != zero) {
			fac_newton = p / dp;							// -delta_Newton
			fac_extra = d2p_half / dp;
			F_half = fac_newton*fac_extra;					// 0.5*F(z)
			denom_sqrt = sqrt(c_one - two_n_div_n_1*F_half);// sqrt( 1 - 2*n/(n-1) * 0.5*F(z) ) = sqrt( 1 - n/(n-1) * F(z) )
															// the old sqrt implementation always has complex angle (0, pi)
															// but the Skowron&Gould requires (-0.5*pi, 0.5*pi), which means positive real part
															// 
															// thus change the sqrt implementation in .h and comment the following check
			/******************************************* changed *******************************************/
			//NEXT LINE PROBABLY CAN BE COMMENTED OUT. Check if compiler outputs positive real
			//if (real(denom_sqrt) >= 0.0) {
			//	denom = c_one_nth + n_1_nth*denom_sqrt;		// denom = 1/n + (n-1)/n * sqrt( 1 - n/(n-1) * F(z) ), 
															// if sqrt returns the result with positive real part
			//}
			//else {
			//	denom = c_one_nth - n_1_nth*denom_sqrt;
			//}

			denom = c_one_nth + n_1_nth*denom_sqrt ;		// denom = 1/n + (n-1)/n * sqrt( 1 - n/(n-1) * F(z) ), 
															// where NEW sqrt returns the result with positive real part
			/*******************************************   end   *******************************************/
			
		}

		if (denom == 0) {
			dx = (absroot + 1.0)*expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI));
		}
		else {
			dx = fac_newton / denom;						// dx = -delta_Newton / denom = -delta_Laguerre, 
															// where delta_Newton and delta_Laguerre are defined in paper
		}

		newroot = *root - dx;								// because dx = -delta_Laguerre, and x_next = x_current + delta_Laguerre

		if (newroot == *root) return; //nothing changes so return
		if (good_to_go) {
			*root = newroot;
			return;
		}
		if (i % FRAC_JUMP_EVERY == 0) { //decide whether to do a jump of modified length (to break cycles)
			faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
			newroot = *root - faq*dx; // do jump of semi-random length
		}
		*root = newroot;
	}
	success = false; // too many iterations here
	return;
}

void VBBinaryLensing::cmplx_newton_spec(complex *poly, int degree, complex *root, int &iter, bool &success) {
	//Subroutine finds one root of a complex polynomial
	//Newton's method. It calculates simplified Adams' stopping 
	//criterion for the value of the polynomial once per 10 iterations (!),
	//after initial iteration. This is done to speed up calculations
	//when polishing roots that are known preety well, and stopping
	// criterion does significantly change in their neighborhood.

	//Uses 'root' value as a starting point (!!!!!)
	//Remember to initialize 'root' to some initial guess.
	//Do not initilize 'root' to point (0,0) if the polynomial 
	//coefficients are strictly real, because it will make going 
	//to imaginary roots impossible.

	// poly - is an array of polynomial cooefs
	//	length = degree+1, poly(1) is constant 
	//0					1				2
	//poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + ...
	//degree - a degree of the polynomial
	// root - input: guess for the value of a root
	//		  output: a root of the polynomial
	//iter - number of iterations performed (the number of polynomial evaluations)
	//success - is false if routine reaches maximum number of iterations

	//For a summary of the method go to: 
	//http://en.wikipedia.org/wiki/Newton's_method

	int FRAC_JUMP_EVERY = 10;
	const int FRAC_JUMP_LEN = 10;
	double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers
	double faq; //jump length
	double FRAC_ERR = 2e-15;
	complex p; //value of polynomial
	complex dp; //value of 1st derivative
	int i, k;
	bool good_to_go;
	complex dx, newroot;
	double ek, absroot, abs2p;
	complex zero = complex(0, 0);
	double stopping_crit2;

	iter = 0;
	success = true;

	//the next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version
	if (false) { //change false to true if you would like to use caustion about haveing first coefficient == 0
		if (degree < 0) {
			printf("Error: cmplx_newton_spec: degree<0");
			return;
		}
		if (poly[degree] == zero) {
			if (degree == 0) return;
			cmplx_newton_spec(poly, degree, root, iter, success);
			return;
		}
		if (degree <= 1) {
			if (degree == 0) {
				success = false;
				printf("Warning: cmplx_newton_spec: degree=0 and poly[0]!=0, no roots");
				return;
			}
			else {
				*root = -poly[0] / poly[1];
				return;
			}
		}
	}
	//end EXTREME Failsafe
	good_to_go = false;

	stopping_crit2 = 0.0; //value not important, will be initialized anyway on the first loop
	for (i = 1; i <= MAXIT; i++) {
		faq = 1.0;
		//prepare stoping criterion
		//calculate value of polynomial and its first two derivatives
		p = poly[degree];
		dp = zero;
		if (i % 10 == 1) { //calculate stopping criterion every tenth iteration
			ek = abs(poly[degree]);
			absroot = abs(*root);
			for (k = degree - 1; k >= 0; k--) {
				dp = p + dp * (*root);
				p = poly[k] + p * (*root); //b_k
										   //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
										   //Communications of ACM, Volume 10 Issue 10, Oct. 1967, p. 655
										   //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
										   //Eq. 8
				ek = absroot * ek + abs(p);
			}
			stopping_crit2 = pow(FRAC_ERR * ek, 2);
		}
		else { // calculate just the value and derivative
			for (k = degree - 1; k >= 0; k--) { //Horner Scheme, see for eg. Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
				dp = p + dp * (*root);
				p = poly[k] + p * (*root);
			}
		}

		iter = iter + 1;

		abs2p = real(conj(p) * p);
		if (abs2p == 0.0) return;
		if (abs2p < stopping_crit2) { //simplified a little Eq. 10 of Adams 1967
			if (dp == zero) return; //if we have problem with zero, but we are close to the root, just accept
									//do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01 * stopping_crit2) return; //return immediatley because we are at very good place
			else {
				good_to_go = true; //do one iteration more
			}
		}

		else {
			good_to_go = false; //reset if we are outside the zone of the root
		}
		if (dp == zero) {
			//problem with zero
			dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i% FRAC_JUMP_LEN] * 2 * M_PI));
		}
		else {
			dx = p / dp; // Newton method, see http://en.wikipedia.org/wiki/Newton's_method
		}
		newroot = *root - dx;
		if (newroot == *root) return; //nothing changes -> return
		if (good_to_go) {//this was jump already after stopping criterion was met
			*root = newroot;
			return;
		}
		if (i % FRAC_JUMP_EVERY == 0) { // decide whether to do a jump of modified length (to break cycles)
			faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
			newroot = *root - faq * dx;
		}
		*root = newroot;
	}
	success = false;
	return;
	//too many iterations here
}

void VBBinaryLensing::cmplx_laguerre2newton(complex *poly, int degree, complex *root, int &iter, bool &success, int starting_mode) {
	//Subroutine finds one root of a complex polynomial using
	//Laguerre's method, Second-order General method and Newton's
	//method - depending on the value of function F, which is a 
	//combination of second derivative, first derivative and
	//value of polynomial [F=-(p"*p)/(p'p')].

	//Subroutine has 3 modes of operation. It starts with mode=2
	//which is the Laguerre's method, and continues until F
	//becames F<0.50, at which point, it switches to mode=1,
	//i.e., SG method (see paper). While in the first two
	//modes, routine calculates stopping criterion once per every
	//iteration. Switch to the last mode, Newton's method, (mode=0)
	//happens when becomes F<0.05. In this mode, routine calculates
	//stopping criterion only once, at the beginning, under an
	//assumption that we are already very close to the root.
	//If there are more than 10 iterations in Newton's mode,
	//it means that in fact we were far from the root, and
	//routine goes back to Laguerre's method (mode=2).

	//Uses 'root' value as a starting point (!!!!!)
	//Remember to initialize 'root' to some initial guess or to 
	//point (0,0) if you have no prior knowledge.

	//poly - is an array of polynomial cooefs
	//	0					1				2
	//	poly[0] x^0 + poly[1] x^1 + poly[2] x^2
	//degree - a degree of the polynomial
	//root - input: guess for the value of a root
	//		output: a root of the polynomial
	//iter - number of iterations performed (the number of polynomial
	//		 evaluations and stopping criterion evaluation)
	//success - is false if routine reaches maximum number of iterations
	//starting_mode - this should be by default = 2. However if you  
	//				  choose to start with SG method put 1 instead.
	//				  Zero will cause the routine to
	//				  start with Newton for first 10 iterations, and
	//				  then go back to mode 2.

	//For a summary of the method see the paper: Skowron & Gould (2012)

	int FRAC_JUMP_EVERY = 10;
	const int FRAC_JUMP_LEN = 10;
	double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers

	double faq; //jump length
	double FRAC_ERR = 2.0e-15;

	complex p; //value of polynomial
	complex dp; //value of 1st derivative
	complex d2p_half; //value of 2nd derivative
	int i, j, k;
	bool good_to_go;
	//complex G, H, G2;
	complex denom, denom_sqrt, dx, newroot;
	double ek, absroot, abs2p, abs2_F_half;
	complex fac_netwon, fac_extra, F_half, c_one_nth;
	double one_nth, n_1_nth, two_n_div_n_1;
	int mode;
	complex c_one = complex(1, 0);
	complex zero = complex(0, 0);
	double stopping_crit2;

	iter = 0;
	success = true;
	stopping_crit2 = 0; //value not important, will be initialized anyway on the first loop

						//next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
	if (false) {//change false to true if you would like to use caution about having first coefficent == 0
		if (degree < 0) {
			printf("Error: cmplx_laguerre2newton: degree < 0");
			return;
		}
		if (poly[degree] == zero) {
			if (degree == 0) return;
			cmplx_laguerre2newton(poly, degree, root, iter, success, starting_mode);
			return;
		}
		if (degree <= 1) {
			if (degree == 0) {//// we know from previous check that poly[0] not equal zero
				success = false;
				printf("Warning: cmplx_laguerre2newton: degree = 0 and poly[0] = 0, no roots");
				return;
			}
			else {
				*root = -poly[0] / poly[1];
				return;
			}
		}
	}
	//end EXTREME failsafe

	j = 1;
	good_to_go = false;

	mode = starting_mode; // mode = 2 full laguerre, mode = 1 SG, mode = 0 newton

	for (;;) { 	//infinite loop, just to be able to come back from newton, if more than 10 iteration there

				// when Laguerre, calculate criterion for each iteration, may go to SG or Newton;
				// (and the implemented Laguerre is 10% faster than zroots due to no module is calculated)
				// when SG,       calculate criterion once per 10 iterations, may go to Newton; 
				// when Newton,   calculate criterion only one time at the beginning, and do at most 10 iterations, 
				//				  if still not returned after 10 iterations, then go to Laguerre
				// claims to be 1.3 times faster than full Laguerre method

			   ////////////
			   ///mode 2///
			   ////////////

		if (mode >= 2) {//Laguerre's method
			one_nth = 1.0 / (degree); ///
			n_1_nth = (degree - 1) * one_nth; ////
			two_n_div_n_1 = 2.0 / n_1_nth;
			c_one_nth = complex(one_nth, 0.0);

			for (i = 1; i <= MAXIT; i++) {
				faq = 1.0;

				//prepare stoping criterion
				ek = abs(poly[degree]);
				absroot = abs(*root);
				//calculate value of polynomial and its first two derivative
				p = poly[degree];
				dp = zero;
				d2p_half = zero;
				for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
					d2p_half = dp + d2p_half * (*root);
					dp = p + dp * (*root);
					p = poly[k - 1] + p * (*root); // b_k
												   //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
												   //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
												   //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
												   //Eq 8.
					ek = absroot * ek + abs(p);
				}
				abs2p = real(conj(p) * p); // abs(p)
				iter = iter + 1;
				if (abs2p == 0) return;

				stopping_crit2 = pow(FRAC_ERR * ek, 2);
				if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
											 //do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01*stopping_crit2) return; // ten times better than stopping criterion
															 //return immediately, because we are at very good place
					else {
						good_to_go = true; //do one iteration more
					}
				}
				else {
					good_to_go = false; //reset if we are outside the zone of the root
				}

				denom = zero;
				if (dp != zero) {
					fac_netwon = p / dp;
					fac_extra = d2p_half / dp;
					F_half = fac_netwon * fac_extra;

					abs2_F_half = real(conj(F_half) * F_half);
					if (abs2_F_half <= 0.0625) {//F<0.50, F/2<0.25
												//go to SG method
						if (abs2_F_half <= 0.000625) {//F<0.05, F/2<0.02
							mode = 0; //go to Newton's
						}
						else {
							mode = 1; //go to SG
						}
					}

					denom_sqrt = sqrt(c_one - two_n_div_n_1*F_half);// sqrt( 1 - 2*n/(n-1) * 0.5*F(z) ) = sqrt( 1 - n/(n-1) * F(z) )
																	// the old sqrt implementation always has complex angle (0, pi)
																	// but the Skowron&Gould requires (-0.5*pi, 0.5*pi), which means positive real part
																	// 
																	// thus change the sqrt implementation in .h and comment the following check

					/******************************************* changed *******************************************/
					//NEXT LINE PROBABLY CAN BE COMMENTED OUT 
					//if (real(denom_sqrt) > 0.0) {
					//	//real part of a square root is positive for probably all compilers. You can \F9
					//	//test this on your compiler and if so, you can omit this check
					//	denom = c_one_nth + n_1_nth * denom_sqrt;	// denom = 1/n + (n-1)/n * sqrt( 1 - n/(n-1) * F(z) ), 
					//												// if sqrt returns the result with positive real part
					//}
					//else {
					//	denom = c_one_nth - n_1_nth * denom_sqrt;
					//}

					denom = c_one_nth + n_1_nth * denom_sqrt ;		// denom = 1/n + (n-1)/n * sqrt( 1 - n/(n-1) * F(z) ), 
																	// where NEW sqrt returns the result with positive real part
					/*******************************************   end   *******************************************/
				}
				if (denom == zero) {//test if demoninators are > 0.0 not to divide by zero
					dx = (abs(*root) + 1.0) + expcmplx(complex(0.0, FRAC_JUMPS[i% FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
				}
				else {
					dx = fac_netwon / denom;
				}
				newroot = *root - dx;
				if (newroot == *root) return; // nothing changes -> return
				if (good_to_go) {//this was jump already after stopping criterion was met
					*root = newroot;
					return;
				}
				if (mode != 2) {
					*root = newroot;
					j = i + 1; //remember iteration index
					break; //go to Newton's or SG
				}
				if ((i% FRAC_JUMP_EVERY) == 0) {//decide whether to do a jump of modified length (to break cycles)
					faq = FRAC_JUMPS[((i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN)];
					newroot = *root - faq * dx; // do jump of some semi-random length (0 < faq < 1)
				}
				*root = newroot;
			} //do mode 2

			if (i >= MAXIT) {
				success = false;
				return;
			}
		}

		////////////
		///mode 1///
		////////////

		if (mode == 1) {//SECOND-ORDER GENERAL METHOD (SG)

			for (i = j; i <= MAXIT; i++) {
				faq = 1.0;
				//calculate value of polynomial and its first two derivatives
				p = poly[degree];
				dp = zero;
				d2p_half = zero;
				if ((i - j) % 10 == 0) {
					//prepare stopping criterion
					ek = abs(poly[degree]);
					absroot = abs(*root);
					for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
						d2p_half = dp + d2p_half * (*root);
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root); //b_k
													   //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
													   //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
													   //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
													   //Eq 8.
						ek = absroot * ek + abs(p);
					}
					stopping_crit2 = pow(FRAC_ERR*ek, 2);
				}
				else {
					for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
						d2p_half = dp + d2p_half * (*root);
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root); //b_k
					}
				}
				abs2p = real(conj(p) * p); //abs(p)**2
				iter = iter + 1;
				if (abs2p == 0.0) return;

				if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
					if (dp == zero) return;
					//do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01*stopping_crit2) return; //ten times better than stopping criterion
															 //ten times better than stopping criterion
					else {
						good_to_go = true; //do one iteration more
					}
				}
				else {
					good_to_go = false; //reset if we are outside the zone of the root
				}
				if (dp == zero) {//test if denominators are > 0.0 not to divide by zero
					dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i% FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
				}
				else {
					fac_netwon = p / dp;
					fac_extra = d2p_half / dp;
					F_half = fac_netwon * fac_extra;

					abs2_F_half = real(conj(F_half) * F_half);
					if (abs2_F_half <= 0.000625) {//F<0.05, F/2<0.025
						mode = 0; //set Newton's, go there after jump
					}
					dx = fac_netwon * (c_one + F_half); //SG
				}
				newroot = *root - dx;
				if (newroot == *root) return; //nothing changes -> return
				if (good_to_go) {
					*root = newroot; //this was jump already after stopping criterion was met
					return;
				}
				if (mode != 1) {
					*root = newroot;
					j = i + 1; //remember iteration number
					break; //go to Newton's
				}
				if ((i% FRAC_JUMP_EVERY) == 0) {// decide whether to do a jump of modified length (to break cycles)
					faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
					newroot = *root - faq * dx; //do jump of some semi random lenth (0 < faq < 1)		
				}
				*root = newroot;
			}
			if (i >= MAXIT) {
				success = false;
				return;
			}

		}
		//------------------------------------------------------------------------------- mode 0
		if (mode == 0) { // Newton's Method

			for (i = j; i <= j + 10; i++) { // Do only 10 iterations the most then go back to Laguerre
				faq = 1.0;

				//calc polynomial and first two derivatives
				p = poly[degree];
				dp = zero;
				if (i == j) { // Calculating stopping criterion only at the beginning
					ek = abs(poly[degree]);
					absroot = abs(*root);
					for (k = degree; k >= 1; k--) {
						dp = p + dp*(*root);
						p = poly[k - 1] + p*(*root);
						ek = absroot*ek + abs(p);
					}
					stopping_crit2 = pow(FRAC_ERR*ek, 2.0);
				}
				else {
					for (k = degree; k >= 1; k--) {
						dp = p + dp*(*root);
						p = poly[k - 1] + p*(*root);
					}
				}
				abs2p = real(conj(p)*p);
				iter = iter + 1;
				if (abs2p == 0.0) return;

				if (abs2p < stopping_crit2) {
					if (dp == zero) return;
					// do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01*stopping_crit2) {
						return; // return immediately since we are at a good place
					}
					else {
						good_to_go = true; // do one more iteration
					}
				}
				else {
					good_to_go = false;
				}

				if (dp == zero) {
					dx = (abs(*root) + 1.0)*expcmplx(complex(0.0, 2 * M_PI*FRAC_JUMPS[i % FRAC_JUMP_LEN])); // make a random jump
				}
				else {
					dx = p / dp;
				}

				newroot = *root - dx;
				if (newroot == *root) return;
				if (good_to_go) {
					*root = newroot;
					return;
				}
				*root = newroot;
			}
			if (iter >= MAXIT) {
				//too many iterations
				success = false;
				return;
			}
			mode = 2; //go back to Laguerre's. Happens when could not converge with 10 steps of Newton
		}

	}/// end of infinite loop
}






void VBBinaryLensing::BinaryMag2_Npoint(double *s, double q,  double rho, \
										double *y1s, double *y2s, \
										int np, \
										double *mags) 
{ 
	//double mags[np] ;

	for (int i = 0; i < np; i++) 
	{
		mags[i] = BinaryMag2(s[i], q, y1s[i], y2s[i], rho) ;
	}

	//return mags ;
}





extern "C"
{

void * wrapBinaryMag2(double s, double q, double x, double y, double rho, double Gamma, double EPSILON, double *Mag)
{
		VBBinaryLensing VBBL;
		VBBL.a1 = Gamma;
		VBBL.Tol = EPSILON;

		//std::ofstream outFile2("test_VBBL.dat");
		//outFile2<<Gamma<<" ";
		//outFile2<<EPSILON<<" ";
		//outFile2.close();

        *Mag = VBBL.BinaryMag2(s, q, x, y, rho);

		return 0;
}

void * wrapBinaryMag2_Npoint(double *s, double q, double rho, \
							 double *x, double *y, \
							 double Gamma, double EPSILON, \
							 int np, \
							 double *mags)
{
		VBBinaryLensing VBBL;
		VBBL.a1 = Gamma;
		VBBL.Tol = EPSILON;

        VBBL.BinaryMag2_Npoint(s, q, rho, x, y, np, mags); 

		return 0;
}

}
