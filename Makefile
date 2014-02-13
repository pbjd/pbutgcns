include pbi.mk

CXXFLAGS := -std=c++0x -Wall -Wuninitialized -Wno-div-by-zero -fpermissive \
			-pedantic -c -fmessage-length=0 -MMD -MP

SRCS := $(wildcard src/*.cpp)
OBJS := $(SRCS:.cpp=.o)
LIBS := -lpbdata -lblasr -llog4cpp -lboost_thread -lboost_program_options -lpthread 

all: OPTIMIZE = -O3

all: pbutgcns

debug: OPTIMIZE = -g -fno-inline

debug: pbutgcns

pbutgcns: $(OBJS)
	$(CXX) $(LIBDIRS) -static -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OPTIMIZE) $(INCDIRS) -Isrc -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"

clean:
	rm -f src/*.d
	rm -f src/*.o
	rm -f pbutgcns

.PHONY: all clean
