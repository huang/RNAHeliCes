# RNAHeliCes and HiKinetics
Abstract folding space analysis based on helices and RNA-kinetics based on folding space abstraction

Copyright (C) 2011-20 Jiabin Huang, Björn Voss.
Send comments/bug reports to: J. Huang <j.huang@uke.de>.


## Dependencies
1. (Optional) if without SUPER privilege, install libboost locally
#install boost without root access
#https://programmer.ink/think/how-to-compile-and-install-boost-libraries-on-linux.html
#under conda-env cc_env  # gcc5.4
#NEED python2.7: downgrade the python from 3.8 --> 2.7
conda install python=2.7
./bootstrap.sh --with-libraries=all /home/jhuang/anaconda3/envs/cc_env/bin/x86_64-conda_cos6-linux-gnu-gcc
#./bootstrap.sh --with-libraries=date_time,program_options --with-toolset=gcc
./b2 install --prefix=/home/jhuang/boost_1_58_installed    #--with-toolset=/home/jhuang/anaconda3/envs/cc_env/bin/gcc
#ldconfig  # requires the root permission!

## Installation
```sh
git clone https://github.com/huang/RNAHeliCes
./configure CFLAGS="-fno-stack-protector" CPPFLAGS="-std=c++98" CXXFLAGS="-std=c++98 -fno-stack-protector"
make
sudo make install
```
   Notes:
     - If there are any linking problems after installing, please check if /usr/local/lib is contained in the environmental variable LD_LIBRARY_PATH. 
     - ./src/libs/libRNA.a contains all Vienna package routines, this may be system dependent. If so, please reload this from Vienna RNA Package.
     - If without SUPER privilege, use 
     ./configure --with-boost-include-path="/home/jhuang/boost_1_58_installed/include" --with-boost-lib-path="/home/jhuang/boost_1_58_installed/lib" CFLAGS='-g -O2 -fno-stack-protector' CPPFLAGS='-std=c++98 -I/home/jhuang/boost_1_58_installed/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2' CXXFLAGS='-std=c++98 -g -O2 -w -fno-stack-protector'
     
## Test run
```sh
RNAHeliCes examples/collosoma_slrna.seq'  #run an example
HiPath -h                                 #get help for calculating energy barriers of an energy landscape
HiPath -f examples/switches_4.faa -k 40 -P ./librna/vienna/rna_turner1999.par  #run an example
HiTed -h                                  #get help for calculating minimum Hishape based Tree edit distance
HiTed examples/riboswitches.fas -t 1 -r 1 #run an example.
```

## HiKinetics
See ./HiKinetics/README for HiKinetics.

## INSTALL and COPYING
See "INSTALL"        for detailed installation instructions, and
    "COPYING"        for disclaimer and copyright.
    
## Citations
[1] Huang, J., & Voß, B. (2011, September). RNAHeliCes—Folding Space Analysis Based on Position Aware Structure Abstraction. In German Conference on Bioinformatics, Weihenstephan, Germany (Vol. 79).
[2] Huang, J., Backofen, R., & Voß, B. (2012). Abstract folding space analysis based on helices. RNA, 18(12), 2135-2147.
[3] Huang, J., & Voß, B. (2014). Analysing RNA-kinetics based on folding space abstraction. BMC bioinformatics, 15(1), 60.
