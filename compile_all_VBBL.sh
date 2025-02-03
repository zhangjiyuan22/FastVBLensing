###Check for the existence of the dynamic library path environment variable(LD_LIBRARY_PATH)
#echo $LD_LIBRARY_PATH
###If there is nothing to be displayed, add a default path value (or not if you wish to)
#LD_LIBRARY_PATH=/usr/local/lib
###add the desired path (need to change the following path according to the location of this folder in your system)
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/zhangjiyuan/FastVBLensing/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:bin/






####################################
### No Optimization version VBBL ###
####################################

### build a dynamic library(.a is static link, .so is dynamic/runtime link)
rm -rf bin/lib_VBBinaryLensingLibraryNoOptimization.so
g++ -fPIC -O3 -g -flto -Wall -Wextra -shared -march=native -o bin/lib_VBBinaryLensingLibraryNoOptimization.so VBBL_lib_no_optimization/VBBinaryLensingLibrary_v3p6.cpp
#g++ -pthread -Wno-unused-result -Wsign-compare -DDYNAMIC_ANNOTATIONS_ENABLED=1 -DNDEBUG -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -fexceptions -fstack-protector-strong -grecord-gcc-switches -m64 -mtune=generic -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection -D_GNU_SOURCE -fPIC -fwrapv -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -fexceptions -fstack-protector-strong -grecord-gcc-switches -m64 -mtune=generic -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection -D_GNU_SOURCE -fPIC -fwrapv -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -fexceptions -fstack-protector-strong -grecord-gcc-switches -m64 -mtune=generic -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection -D_GNU_SOURCE -fPIC -fwrapv -fPIC -I/usr/include/python3.9 -c VBBL_lib_original/VBBinaryLensingLibrary_v3p6.cpp -o bin/lib_VBBinaryLensingLibraryOriginal.so -std=c++17
chmod -x bin/lib_VBBinaryLensingLibraryNoOptimization.so


### build the test code which calls the dynamic library
rm -rf bin/test_VBBLNoOptimization.out
#the source file should be in front of the dynamic library
g++ -O3 -g -Wall -Wextra -march=native test_VBBLNoOptimization.cpp -Lbin -l_VBBinaryLensingLibraryNoOptimization -o bin/test_VBBLNoOptimization.out






###########################################
### Compiling Optimization version VBBL ###
###########################################

### build a dynamic library(.a is static link, .so is dynamic/runtime link)
rm -rf bin/lib_VBBinaryLensingLibraryCompilingOptimization.so
g++ -fPIC -O3 -g -flto -Wall -Wextra -shared -march=native -o bin/lib_VBBinaryLensingLibraryCompilingOptimization.so VBBL_lib_compiling_optimization/VBBinaryLensingLibrary_v3p6.cpp
chmod -x bin/lib_VBBinaryLensingLibraryCompilingOptimization.so


### build the test code which calls the dynamic library
rm -rf bin/test_VBBLCompilingOptimization.out
#the source file should be in front of the dynamic library
g++ -O3 -g -Wall -Wextra -march=native test_VBBLCompilingOptimization.cpp -Lbin -l_VBBinaryLensingLibraryCompilingOptimization -o bin/test_VBBLCompilingOptimization.out






#######################################################
### Algorithmic Compiling Optimization version VBBL ###
#######################################################

### build a dynamic library(.a is static link, .so is dynamic/runtime link)
rm -rf bin/lib_VBBinaryLensingLibraryAlgorithmicCompilingOptimization.so
g++ -fPIC -O3 -g -flto -Wall -Wextra -shared -march=native -o bin/lib_VBBinaryLensingLibraryAlgorithmicCompilingOptimization.so VBBL_lib_algorithmic_compiling_optimization/VBBinaryLensingLibrary_v3p6.cpp
chmod -x bin/lib_VBBinaryLensingLibraryAlgorithmicCompilingOptimization.so


### build the test code which calls the dynamic library
rm -rf bin/test_VBBLAlgorithmicCompilingOptimization.out
#the source file should be in front of the dynamic library
g++ -O3 -g -Wall -Wextra -march=native test_VBBLAlgorithmicCompilingOptimization.cpp -Lbin -l_VBBinaryLensingLibraryAlgorithmicCompilingOptimization -o bin/test_VBBLAlgorithmicCompilingOptimization.out

