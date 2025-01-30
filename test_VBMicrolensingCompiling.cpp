/**************************************************************************************/
// this code is calling VBMicrolensing to calculate a series of points' magnification, 
// to test the speed of VBMicrolensing code at different magnification
/**************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <chrono>

#include"VBMicrolensing_lib_compiling_optimization/VBMicrolensingLibrary.h"

int main(int argc, char *argv[]) 
{
    //printf("Number of command-line arguments: %d\n", argc);
    if (argc != 4)
    {
        printf("wrong number of argument!") ;
        exit(0) ;
    }

    //printf("s   = %s\n", argv[1]);
    //printf("q   = %s\n", argv[2]);
    //printf("rho = %s\n", argv[3]);

    printf("y_min = %s\n", argv[1]);
    printf("x_min = %s\n", argv[2]);

    printf("file number = %s\n", argv[3]);
    

    //double s   = atof( argv[1] ) ;
    //double q   = atof( argv[2] ) ;
    //double rho = atof( argv[3] ) ;

    //int Npoint_y = atoi( argv[4] ) ;
    //int Npoint_x = atoi( argv[5] ) ;
    double y_min = atof( argv[1] ) ;
    double x_min = atof( argv[2] ) ;

    char file_number = argv[3][0] ;


    // set by yourself. mass center is the origin here. 
    //double y_min = -0.1 ;
    //double x_min = -0.1 ;
    int Npoint_y = 251 ;
    int Npoint_x = 251 ;

    double y_max = -y_min ;
    double x_max = -x_min ;

    
    // shift origin from mass center to magnification center
    /*
    double shift_y = 0. ;
    double shift_x ;

    if(s <= 1)
    {
        shift_x = 0. ;
    }
    else // s > 1
    {
        shift_x = (-1.)*(s-1./s)*q/(1.+q) ;
    }

    y_min += shift_y ;
    y_max += shift_y ;

    x_min += shift_x ;
    x_max += shift_x ;
    */
    
    // calculate dy (and dx) according to y(x)_min, y(x)_max, Npoint_y(x)
    double dy = (y_max - y_min) / (double)Npoint_y;
    double dx = (x_max - x_min) / (double)Npoint_x;


    double  q_array[3] = { 1., 0.001, 0.0001} ;
    complex s_array[3] = { complex(0.,0.), complex(1.,0.), complex(0.,0.9)} ;
    double  rho = 0.001 ;


    // first declare an instance to the VBMicrolensing class
	VBMicrolensing VBML;
    // Choose the algorithm you want to use (Nopoly, Multipoly, Singlepoly).
	VBML.SetMethod(VBMicrolensing::Method::Multipoly);

    // Call the function setLensGeometry, 
    // which initializes the arrays describing the geometric configuration of the lens system 
    // based on the chosen algorithm and configuration.
	VBML.SetLensGeometry(3, q_array, s_array);

    VBML.a1  = 0.;
    VBML.Tol = 0.001;
    VBML.RelTol = 0.0001;

    
    // use a variable to store a point's magnification
    double mag ;

    // declare arrays to store all points' magnification and computation time
    float * magnification_array ;
    magnification_array    = (float *)malloc( Npoint_y * Npoint_x * 4 ) ;
    float * computation_time_array ;
    computation_time_array = (float *)malloc( Npoint_y * Npoint_x * 4 ) ;

    // measure the total time
    auto begin_total = std::chrono::high_resolution_clock::now() ;

    // double loop to calculate ( Npoint_y * Npoint_x ) points' magnification using BinaryMag2
    for(int y_arg=0; y_arg < Npoint_y; y_arg++)
    {   
        double y_current = y_min + y_arg * dy ;

        for(int x_arg=0; x_arg < Npoint_x; x_arg++)
        {
            double x_current = x_min + x_arg * dx ; 

            int total_arg = y_arg * Npoint_x + x_arg ;

            auto begin = std::chrono::high_resolution_clock::now();
            
            mag = VBML.MultiMag( complex(x_current,y_current), rho, 0.001) ;
            
            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

            magnification_array[total_arg]    = (float)mag ;
            computation_time_array[total_arg] = (float)(elapsed_time.count() * 1e-9) ; // nanosecond to second
            //printf("CPU code                   needs %e (second)\n", cpu_code_time) ; 
        }
    }
    
    auto end_total   = std::chrono::high_resolution_clock::now() ;
    auto elapsed_time_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end_total - begin_total) ;
    float total_time = (float)(elapsed_time_total.count() * 1e-9) ; // nanosecond to second
    printf("total needs %e (second)\n", total_time) ;

    //printf("%f\n", magnification_array[0]) ;
    //printf("%e\n", computation_time_array[0]) ;
    
    // write arrays into a file
    FILE *file_pointer ;

    
    // optimized version
    char file_name[] = "../result/test_VBMicrolensing_result_compiling_optimization_1_with_RelTol_1eminus4.txt" ;
    file_name[60] = file_number ;
    

    file_pointer = fopen(file_name,"w") ; 
    if(file_pointer == NULL) 
    {
        printf("The file is not opened. The program will exit now") ;
        exit(0) ;
    }
    else
    {
        printf("file:%s is now opened.\n", file_name) ;

        for(int arg=0; arg < Npoint_y * Npoint_x; arg++)
        {
            //fprintf(file_pointer, "%f\n", magnification_array[arg]) ; 
            fprintf(file_pointer, "%f %e\n", magnification_array[arg], computation_time_array[arg]) ; 
        }
        
        fclose(file_pointer) ; 

        printf("Data successfully written in file.\n") ;
        printf("The file is now closed.\n") ;
    }
    

    // free the allocated memory space
    free(magnification_array) ;
    free(computation_time_array) ;

    return 0;
}
