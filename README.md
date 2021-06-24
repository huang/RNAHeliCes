# RNAHeliCes and HiKinetics
Abstract folding space analysis based on helices and RNA-kinetics based on folding space abstraction

Copyright (C) 2011-20 Jiabin Huang, Björn Voss.
Send comments/bug reports to: J. Huang <j.huang@uke.de>.

## Dependencies
RNAliHiKinetics needs RNAliHelices and RNAliHipath whose installation is described below. Additional dependencies are:

compile time:
* C++ compiler (for example GCC g++)
* C compiler (for example GCC)
* GNU make >= 3.81

compile time and runtime:
* Boost Libraries (>1.58): program_options, date_time
* LAPACK 

## Installation
```sh
git clone https://github.com/huang/RNAHeliCes
./configure CFLAGS="-fno-stack-protector" CPPFLAGS="-std=c++98" CXXFLAGS="-std=c++98 -fno-stack-protector"
make
sudo make install
```
Notes:
  - If there are any linking problems after installing, please add '/usr/local/lib' to LD_LIBRARY_PATH.
```sh
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export LD_LIBRARY_PATH
```
  - ./src/libs/libRNA.a contains all Vienna package routines, this may be system dependent. If so, please reload this from Vienna RNA Package.
  - If without SUPER privilege, use 
  ./configure --with-boost-include-path="/your_path/boost_1_58_installed/include" --with-boost-lib-path="/your_path/boost_1_58_installed/lib" CFLAGS='-g -O2 -fno-stack-protector' CPPFLAGS='-std=c++98 -I/your_path/boost_1_58_installed/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2' CXXFLAGS='-std=c++98 -g -O2 -w -fno-stack-protector'
     
## Test run
```sh
RNAHeliCes examples/collosoma_slrna.seq
HiPath -f examples/switches_4.faa -k 40 -P ./librna/vienna/rna_turner1999.par
HiTed examples/riboswitches.fas -t 1 -r 1
```

## HiKinetics
See ./HiKinetics/README for HiKinetics.
    
## Citations
  [1] Huang, J., & Voß, B. (2011, September). RNAHeliCes—Folding Space Analysis Based on Position Aware Structure Abstraction. In German Conference on Bioinformatics, Weihenstephan, Germany (Vol. 79).
  
  [2] Huang, J., Backofen, R., & Voß, B. (2012). Abstract folding space analysis based on helices. RNA, 18(12), 2135-2147.
  
  [3] Huang, J., & Voß, B. (2012, September). Reducing the search space in RNA helix based folding. European Conference on Computational Biology 2012. European Conference on Computational Biology (ECCB) 2012.
  
  [4] Huang, J., & Voß, B. (2014). Analysing RNA-kinetics based on folding space abstraction. BMC bioinformatics, 15(1), 60.
