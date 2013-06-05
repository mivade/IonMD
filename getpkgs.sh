#!/bin/bash
#
# This is just a simple script to grab all the dependencies you need
# on Debian wheezy. As with **NOTE THAT THIS USES SUDO, SO READ THIS
# AND USE AT YOUR OWN RISK**
#
# I wanted to make this general, but if you know how to apt-pin,
# jessie has packages for nlopt

sudo aptitude install build-essential libgsl0-dbg libgsl0-dev git python-numpy python-scipy mayavi2 python-imaging python-matplotlib python-ctypeslib
wget http://ab-initio.mit.edu/nlopt/nlopt-2.3.tar.gz
tar -xvzf nlopt-2.3.tar.gz
cd nlopt-2.3
./configure --enable-debug --with-cxx
make
sudo make install

