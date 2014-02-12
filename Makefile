# C++11 is required for some convenience issues with Armadillo, namely
# using initializers.
CXX = g++
CXXFLAGS = -Wall -O3 -fPIC -flto -fopenmp -march=native -g -std=c++11
LIBS = -lm -lgsl -lgslcblas
SRCS = *.cpp
HDRS = *.hpp
OBJS = minimize.o
RESTFLAGS = --math-output="MathJax" --cloak-email-addresses

ionmd.so: $(SRCS) $(HDRS) params.py README.html minimize.o
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(LIBS) ionmd.cpp -oionmd.so

minimize.o: minimize.cpp minimize.hpp
	$(CXX) $(CXXFLAGS) -c minimize.cpp -ominimize.o

params.py: params.hpp
	h2xml $(PWD)/params.hpp -o params.xml
	xml2py params.xml -s Params -o params.py

README.html: README.rst
	rst2html $(RESTFLAGS) README.rst README.html

clean:
	rm -f *.o
	rm *.so
	rm params.py
	rm params.pyc
