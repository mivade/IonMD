CXX = g++
CXXFLAGS = -Wall -O3 -fPIC -flto -fopenmp -march=native -g
LIBS = -lm -lgsl -lgslcblas
SRCS = *.cpp
HDRS = *.hpp

ionmd.so: $(SRCS) $(HDRS) params.py
	$(CXX) $(CXXFLAGS) -shared $(LIBS) $(OBJS) ionmd.cpp -oionmd.so

params.py: params.hpp
	h2xml $(PWD)/params.hpp -o params.xml
	xml2py params.xml -s Params -o params.py

clean:
	rm -f *.o
	rm *.so
