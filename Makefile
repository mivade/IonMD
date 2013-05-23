CXX = g++
CXXFLAGS = -Wall -O3 -fPIC -flto -fopenmp -march=native -g
LIBS = -lm -lgsl -lgslcblas -lnlopt
SRCS = *.cpp
HDRS = *.hpp
RESTFLAGS = --math-output="MathJax" --cloak-email-addresses

ionmd.so: $(SRCS) $(HDRS) params.py README.html
	$(CXX) $(CXXFLAGS) -shared $(LIBS) $(OBJS) ionmd.cpp -oionmd.so

# apparently this doesn't work right on 64 bit machines
# set Params._pack_ = 4 to Params._pack_ = 8 in params.py
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