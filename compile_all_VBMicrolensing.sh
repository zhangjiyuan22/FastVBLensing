#Check for the existence of the dynamic library path environment variable(LD_LIBRARY_PATH)
echo $LD_LIBRARY_PATH
#If there is nothing to be displayed, add a default path value (or not if you wish to)
LD_LIBRARY_PATH=/usr/local/lib
#add the desired path (need to change the following path according to the location of this folder in your system)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/zhangjiyuan/VBMicrolensing_and_VBBL/bin/






##############################################
### No Optimization version VBMicrolensing ###
##############################################

### build a dynamic library(.a is static link, .so is dynamic/runtime link)
rm -rf bin/lib_VBMicrolensingLibraryNoOptimization.so
g++ -fPIC -O3 -g -flto -Wall -Wextra -shared -march=native -o bin/lib_VBMicrolensingLibraryNoOptimization.so VBMicrolensing_lib_no_optimization/VBBinaryMag.cpp VBMicrolensing_lib_no_optimization/VBclasses.cpp VBMicrolensing_lib_no_optimization/VBESPL.cpp VBMicrolensing_lib_no_optimization/VBLightCurves.cpp VBMicrolensing_lib_no_optimization/VBLimbDarkening.cpp VBMicrolensing_lib_no_optimization/VBMicrolensingLibrary.cpp VBMicrolensing_lib_no_optimization/VBMultiMag.cpp VBMicrolensing_lib_no_optimization/VBparallaxcaustics.cpp VBMicrolensing_lib_no_optimization/VBpolynomials.cpp VBMicrolensing_lib_no_optimization/VBSkowronGould.cpp
chmod -x bin/lib_VBMicrolensingLibraryNoOptimization.so


### build the test code which calls the dynamic library
rm -rf bin/test_VBMicrolensingNoOptimization.out
#the source file should be in front of the dynamic library
g++ -O3 -g -Wall -Wextra -march=native test_VBMicrolensingNoOptimization.cpp -Lbin -l_VBMicrolensingLibraryNoOptimization -o bin/test_VBMicrolensingNoOptimization.out






#####################################################
### Compiling Optimization version VBMicrolensing ###
#####################################################

### build a dynamic library(.a is static link, .so is dynamic/runtime link)
rm -rf bin/lib_VBMicrolensingLibraryCompilingOptimization.so
g++ -fPIC -O3 -g -flto -Wall -Wextra -shared -march=native -o bin/lib_VBMicrolensingLibraryCompilingOptimization.so VBMicrolensing_lib_compiling_optimization/VBMicrolensingLibrary.cpp
chmod -x bin/lib_VBMicrolensingLibraryCompilingOptimization.so


### build the test code which calls the dynamic library
rm -rf bin/test_VBMicrolensingCompilingOptimization.out
#the source file should be in front of the dynamic library
g++ -O3 -g -Wall -Wextra -march=native test_VBMicrolensingCompilingOptimization.cpp -Lbin -l_VBMicrolensingLibraryCompilingOptimization -o bin/test_VBMicrolensingCompilingOptimization.out







#################################################################
### Algorithmic Compiling Optimization version VBMicrolensing ###
#################################################################

### build a dynamic library(.a is static link, .so is dynamic/runtime link)
rm -rf bin/lib_VBMicrolensingLibraryAlgorithmicCompilingOptimization.so
g++ -fPIC -O3 -g -flto -Wall -Wextra -shared -march=native -o bin/lib_VBMicrolensingLibraryAlgorithmicCompilingOptimization.so VBMicrolensing_lib_algorithmic_compiling_optimization/VBMicrolensingLibrary.cpp
chmod -x bin/lib_VBMicrolensingLibraryAlgorithmicCompilingOptimization.so


### build the test code which calls the dynamic library
rm -rf bin/test_VBMicrolensingAlgorithmicCompilingOptimization.out
#the source file should be in front of the dynamic library
g++ -O3 -g -Wall -Wextra -march=native test_VBMicrolensingAlgorithmicCompilingOptimization.cpp -Lbin -l_VBMicrolensingLibraryAlgorithmicCompilingOptimization -o bin/test_VBMicrolensingAlgorithmicCompilingOptimization.out




