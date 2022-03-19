# Commands
CXX = c++

USE_OPENMP = 1
USE_METIS = 1
WLSLIB_HOME =
CPPFLAGS = -I${WLSLIB_HOME}/cpp/extern -I${WLSLIB_HOME}/cpp/src -I.

ifeq (${DEBUG},)
OPTFLAGS = -O3 -DNDEBUG
else ifeq (${DEBUG},1)
OPTFLAGS = -O2
else
OPTFLAGS = -O3 -DNDEBUG -fno-omit-frame-pointer -fsanitize=address
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
CPPFLAGS += -DRDI_USE_METIS
LIBRDI_LIBS += -lmetis
endif

LIBRDI_LIBS += ${LIBS}

VPATH = codegen_src/

librdi.so: rdi.o
	${CXX} -shared -o $@ $< ${LIBRDI_LIBS} ${LDFLAGS}

rdi.o: rdi.cpp rdi.hpp $(shell ls codegen_src/)
	@if [ -z "${WLSLIB_HOME}" ]; then \
		echo "\nERROR: WLSLIB_HOME=/path/to/wlslib must be passed in\n" 1>&2; \
		exit 1; \
	fi
	${CXX} -c ${CPPFLAGS} ${LIBRDI_CXXFLAGS} $<

.PHONY: clean

clean:
	rm -f *.so *.o
