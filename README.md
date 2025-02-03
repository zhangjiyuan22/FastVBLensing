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
9. ./test_VBBLAlgorithmicCompilingOptimization.out 1.0 0.001 0.001 -1.0 -1.0 1
    (which means s=1.0, q=0.001, rho=0.001, x_range/y_range=np.linspace(-1.0, 1.0, 251), resulting file named as '1', default Tol=1e-3, RelTol=1e-4, No Limb-Darkening)
10. ./test_VBBLAlgorithmicCompilingOptimization.out 1.0 0.001 0.001 -0.1 -0.1 2
    (which means s=1.0, q=0.001, rho=0.001, x_range/y_range=np.linspace(-0.1, 0.1, 251), resulting file named as '2', default Tol=1e-3, RelTol=1e-4, No Limb-Darkening)
11. ./test_VBBLAlgorithmicCompilingOptimization.out 1.0 0.001 0.001 -0.01 -0.01 3
    (which means s=1.0, q=0.001, rho=0.001, x_range/y_range=np.linspace(-0.01, 0.01, 251), resulting file named as '3', default Tol=1e-3, RelTol=1e-4, No Limb-Darkening)
    (Purpose of running three times with different ranges: get enough points in low/medium/high magnification region, respectively)
   
