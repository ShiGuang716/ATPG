CODE = -DBYTERLE


#compilers
#PCC = g++
PCC = mpicxx
PCC_MPI = mpicxx
# PCFLAGS = -fcilkplus -lcilkrts -std=c++14 -DCILK -g -O3  $(INTT) $(INTE) $(CODE) $(PD) $(MEM) -openmp -DMPI_DEBUG
#PCFLAGS = -fcilkplus -lcilkrts -std=c++14 -DCILK -g -O3 $(INTT) $(INTE) $(CODE) $(PD) $(MEM) -openmp -DSORT
PCFLAGS = -fcilkplus -lcilkrts -std=c++14 -DCILK -g -O3 $(INTT) $(INTE) $(CODE) $(PD) $(MEM) -openmp

COMMON= ligra.h graph.h compressedVertex.h vertex.h utils.h IO.h parallel.h gettime.h index_map.h maybe.h sequence.h edgeMap_utils.h binary_search.h quickSort.h blockRadixSort.h transpose.h parseCommandLine.h byte.h byteRLE.h nibble.h byte-pd.h byteRLE-pd.h nibble-pd.h vertexSubset.h  encoder.C decoder.C 


ALL =  ATPG_FAN_MPI8
# ALL=  ATPG_FAN GOODSIM ATPG_FAN_MPI5 ATPG_FAN_MPI6 
# ALL= ATPG_FAN_MPI ATPG_FAN_MPI4 GOODSIM ATPG_FAN ATPG_FAN2


all: $(ALL)

# FAN: FAN.C $(COMMON) preprocess.h define.h FAN1.h implication.h
# 	$(PCC) $(PCFLAGS) -o $@ $<
# FAN: FAN.C $(COMMON) preprocess.h define.h FAN.h implication.h
# 	$(PCC) $(PCFLAGS) -o $@ $<

% : %.C $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

$(COMMON):
	ln -s ../ligra/$@ .

.PHONY : clean

clean :
	rm -f *.o $(ALL)

cleansrc :
	rm -f *.o $(ALL)
	rm $(COMMON)
