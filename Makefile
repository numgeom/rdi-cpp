# Commands
CXX = c++

USE_OPENMP = 1
USE_METIS = 1
WLSLIB_HOME = ../../wlslib
MOMP2CPP_HOME = ../../momp2cpp
CPPFLAGS = -Isrc -I${WLSLIB_HOME}/cpp/src -I${WLSLIB_HOME}/cpp/src/extern

ifeq (${DEBUG},)
OPTFLAGS = -O3 -DNDEBUG
else ifeq (${DEBUG},1)
OPTFLAGS = -O2
else
OPTFLAGS = -O0 -fno-omit-frame-pointer -fsanitize=address
endif

LIBRDI_CXXFLAGS = -fPIC -g -Wall -Wno-comment ${OPTFLAGS} ${CXXFLAGS}

ifeq (${LAPACK_LIBS},)
LAPACK_LIBS = -llapack -lblas
endif

ifneq (${LAPACK_LIB_PATH},)
LIBRDI_LIBS = -L${LAPACK_LIB_PATH} ${LAPACK_LIBS}
else
LIBRDI_LIBS = ${LAPACK_LIBS}
endif

ifeq (${USE_OPENMP},1)
LIBRDI_CXXFLAGS += -fopenmp
LIBRDI_LIBS += -fopenmp
endif

ifeq (${USE_METIS},1)
CPPFLAGS += -I${MOMP2CPP_HOME}/include -DRDI_USE_METIS
LIBRDI_LIBS += -L${MOMP2CPP_HOME}/lib/glnxa64 -lmetis
endif

LIBRDI_LIBS += ${LIBS}

VPATH = src src/codegen_src

librdi.so: rdi.o
	${CXX} -shared -o $@ $< ${LIBRDI_LIBS} ${LDFLAGS}

rdi.o: rdi.cpp rdi.hpp $(shell ls src/codegen_src/)
	@if [ -z "${WLSLIB_HOME}" ]; then \
		echo "\nERROR: WLSLIB_HOME=/path/to/wlslib must be passed in\n" 1>&2; \
		exit 1; \
	fi
	${CXX} -c ${CPPFLAGS} ${LIBRDI_CXXFLAGS} $<

.PHONY: clean

clean:
	rm -f *.so *.o
