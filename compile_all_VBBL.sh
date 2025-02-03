#Check for the existence of the dynamic library path environment variable(LD_LIBRARY_PATH)
echo $LD_LIBRARY_PATH
#If there is nothing to be displayed, add a default path value (or not if you wish to)
LD_LIBRARY_PATH=/usr/local/lib
#add the desired path (need to change the following path according to the location of this folder in your system)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/zhangjiyuan/VBMicrolensing_and_VBBL/bin/






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




##############################
### optimized version VBBL ###
##############################

### build a dynamic library(.a is static link, .so is dynamic/runtime link)
rm -rf bin/lib_VBBinaryLensingLibraryOptimized.so
g++ -fPIC -O3 -g -flto -Wall -Wextra -shared -march=skylake-avx512 -o bin/lib_VBBinaryLensingLibraryOptimized.so VBBL_lib_optimized/VBBinaryLensingLibrary_v3p6.cpp
chmod -x bin/lib_VBBinaryLensingLibraryOptimized.so


### build the test code which calls the dynamic library
rm -rf bin/test_VBBLOptimized.out
#the source file should be in front of the dynamic library
g++ -O3 -g -Wall -Wextra -march=skylake-avx512 test_VBBLOptimized.cpp -Lbin -l_VBBinaryLensingLibraryOptimized -o bin/test_VBBLOptimized.out










##############################
### 4.4 version VBBL ###
##############################

### build a dynamic library(.a is static link, .so is dynamic/runtime link)
rm -rf bin/lib_VBBinaryLensingLibraryV4p4.so
g++ -fPIC -O3 -g -flto -Wall -Wextra -shared -march=skylake-avx512 -o bin/lib_VBBinaryLensingLibraryV4p4.so VBBL_lib_v4p4/VBBinaryLensingLibrary_v3p6.cpp
#g++ -fPIC -O1 -g -flto -Wall -Wextra -shared -o bin/lib_VBBinaryLensingLibraryV4p3.so VBBL_lib_v4p3/VBBinaryLensingLibrary_v3p6.cpp
chmod -x bin/lib_VBBinaryLensingLibraryV4p4.so


### build the test code which calls the dynamic library
rm -rf bin/test_VBBLV4p4.out
#the source file should be in front of the dynamic library
g++ -O3 -g -Wall -Wextra -march=skylake-avx512 test_VBBLV4p4.cpp -Lbin -l_VBBinaryLensingLibraryV4p4 -o bin/test_VBBLV4p4.out
#g++ -O1 -g -Wall -Wextra test_VBBLV4p3.cpp -Lbin -l_VBBinaryLensingLibraryV4p3 -o bin/test_VBBLV4p3.out

