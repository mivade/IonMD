#!/bin/bash
#
# This is just a simple script to grab all the dependencies you need
# on Debian wheezy. **NOTE THAT THIS USES SUDO, SO READ THIS AND USE
# AT YOUR OWN RISK**
#
# I wanted to make this general, but if you know how to apt-pin,
# jessie has packages for nlopt
#
# This is a work in progress... I'm adding to this as I find new
# problems. I'll get around to fully testing it at some point in the
# near future.

sudo aptitude install build-essential libgsl0-dbg libgsl0-dev git python-numpy python-scipy mayavi2 python-imaging python-ctypeslib python-pip
sudo pip install matplotlib # wheezy's version is out of date
wget http://ab-initio.mit.edu/nlopt/nlopt-2.3.tar.gz
tar -xvzf nlopt-2.3.tar.gz
cd nlopt-2.3
./configure --enable-debug --with-cxx
make
sudo make install
mkdir data
mkdir images
