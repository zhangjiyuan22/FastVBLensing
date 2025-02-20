# FastVBLensing
Speed Up VBBinaryLensing and VBMicroLensing Through Compiling and Algorithmic Optimization

Search for * changed *  (which is star space changed space star) to find all changes

Rich comments to VBBL are added in VBBL_lib_algorithmic_compiling_optimization/ 

## How to run?
1. git clone https://github.com/zhangjiyuan22/FastVBLensing.git
2. cd FastVBLensing
3. chmod a+x compile_all_VBBL.sh
4. chmod a+x compile_all_VBMicrolensing.sh
5. ./compile_all_VBBL.sh
6. ./compile_all_VBMicrolensing.sh
7. export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:bin/
8. cd bin
### VBBL
10. ./test_VBBLAlgorithmicCompilingOptimization.out 1.0 0.001 0.001 -1.0 -1.0 1
    <br>(which means s=1.0, q=0.001, rho=0.001, x_range/y_range=np.linspace(-1.0, 1.0, 251), resulting file named as '1', <br>default Tol=1e-3, RelTol=1e-4, No Limb-Darkening, using BinaryMag2)<br>
    (x_range/y_range centers on Primary Lens)<br>
    (actually, it's no strictly np.linspace(-1.0, 1.0, 251), the ending point is not included, but just use this notation)
11. ./test_VBBLAlgorithmicCompilingOptimization.out 1.0 0.001 0.001 -0.1 -0.1 2
    <br>(which means s=1.0, q=0.001, rho=0.001, x_range/y_range=np.linspace(-0.1, 0.1, 251), resulting file named as '2', <br>default Tol=1e-3, RelTol=1e-4, No Limb-Darkening, using BinaryMag2)
12. ./test_VBBLAlgorithmicCompilingOptimization.out 1.0 0.001 0.001 -0.01 -0.01 3
    <br>(which means s=1.0, q=0.001, rho=0.001, x_range/y_range=np.linspace(-0.01, 0.01, 251), resulting file named as '3', <br>default Tol=1e-3, RelTol=1e-4, No Limb-Darkening, using BinaryMag2)
    <br>(Purpose of running three times with different ranges: get enough points in low/medium/high magnification region, respectively)
13. cd ../result
14. ls <br>(now there are three files generated. in each file, first column is magnification, second column is computation time in second.)
15. repeat step 9-11 with other two versions, then you can compare the performance between No Optimization, Compiling Optimization, and Algorithmic Compiling Optimization.
16. you can compare the correctness of results between different versions using e.g. diff <(awk '{print $1}' test_VBBL_result_algorithmic_compiling_optimization_3_with_RelTol_1eminus4.txt) <(awk '{print $1}' test_VBBL_result_compiling_optimization_3_with_RelTol_1eminus4.txt)
#### typical time used for VBBL
|                                    | x_range/y_range= np.linspace(-1.0,1.0,251) | x_range/y_range= np.linspace(-0.1,0.1,251) | x_range/y_range= np.linspace(-0.01,0.01,251) |
|------------------------------------|--------------------------------------------|--------------------------------------------|----------------------------------------------|
| Algorithmic Compiling Optimization | 0.199 s                                    | 4.837 s                                    | 12.654 s                                     |
| Compiling Optimization             | 0.202 s                                    | 5.171 s                                    | 18.400 s                                     |
| No Optimization                    | 0.801 s                                    | 12.154 s                                   | 34.111 s                                     |
### VBMicrolensing
17. ./test_VBMicrolensingAlgorithmicCompilingOptimization.out -1.0 -1.0 1
    <br>(which means x_range/y_range=np.linspace(-1.0, 1.0, 251), resulting file named as '1', <br>default s2=1.0, q2=0.001, s3=0.9, q3=0.0001, psi=90 degree, rho=0.001, Tol=1e-3, RelTol=1e-4, No Limb-Darkening, using MultiMag with Multipoly method)
    <br>(x_range/y_range centers on Primary Lens)
19. run remaining two ranges and again other two versions like in VBBL, and compare results   
#### typical time used for VBMicrolensing
|                                    | x_range/y_range= np.linspace(-1.0,1.0,251) | x_range/y_range= np.linspace(-0.1,0.1,251) | x_range/y_range= np.linspace(-0.01,0.01,251) |
|------------------------------------|--------------------------------------------|--------------------------------------------|----------------------------------------------|
| Algorithmic Compiling Optimization | 13.591 s                                   | 22.230 s                                   | 94.480 s                                     |
| Compiling Optimization             | 13.492 s                                   | 23.129 s                                   | 137.472 s                                    |
| No Optimization                    | 43.300 s                                   | 69.079 s                                   | 309.948 s                                    |
