
WIN = 0

ifeq ($(WIN),1)
	CUDA_ROOT = $(shell echo $$CUDA_PATH)
	NVCC = "$(CUDA_ROOT)/bin/nvcc"
	CUDA_LIBS = "$(CUDA_ROOT)/lib/win32/cuda.lib"
	EXT = dll
	NVCCFLAGS = -Xptxas -v -Xcudafe -\# -shared -Xptxas -abi=no
else
	NVCC = "$(shell which nvcc)"
	CUDA_ROOT = $(shell dirname $(NVCC))/../
	CUDA_LIBS = -lcuda
	EXT = so
	NVCCFLAGS = -Xptxas -v -Xcudafe -\# -shared -Xptxas -abi=no \
			-Xcompiler -fPIC -Xcompiler -fvisibility=hidden
endif

GEN_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\" 
GEN_SM13 = -gencode=arch=compute_13,code=\"sm_13,compute_13\" 
GEN_SM10 = -gencode=arch=compute_10,code=\"sm_10,compute_10\" 

INC = -I"$(CUDA_ROOT)/include" -I. 

DEPS = ./Makefile \
	sort_engine.cu \
	sort_engine.h \
	$(wildcard b40c/util/*.cuh) \
	$(wildcard b40c/util/**/*.cuh) \
	$(wildcard b40c/radix_sort/*.cuh) \
	$(wildcard b40c/radix_sort/**/*.cuh) 

LIBNAME = sort_engine

all: $(LIBNAME)_sm10.$(EXT) $(LIBNAME)_sm13.$(EXT) $(LIBNAME)_sm20.$(EXT)
	touch built

clean :
	rm -f  *.$(EXT) *.lib *.exp *.dll built

$(LIBNAME)_sm10.$(EXT) : $(DEPS)
	$(NVCC) $(GEN_SM10) -o $@ sort_engine.cu $(NVCCFLAGS) $(INC) -O3  

$(LIBNAME)_sm13.$(EXT) : $(DEPS)
	$(NVCC) $(GEN_SM13) -o $@ sort_engine.cu $(NVCCFLAGS) $(INC) -O3  

$(LIBNAME)_sm20.$(EXT) : $(DEPS)
	$(NVCC) $(GEN_SM20) -o $@ sort_engine.cu $(NVCCFLAGS) $(INC) -O3  

