# make trae FC=f77 por defecto; se reemplaza por gfortran salvo que el usuario defina FC
ifeq ($(origin FC),default)
FC = gfortran
endif
FFLAGS ?= -O2

SRC := $(wildcard src/*.f90)
BIN := $(patsubst src/%.f90,bin/%,$(SRC))

.PHONY: all clean

all: $(BIN)

bin/%: src/%.f90 | bin
	$(FC) $(FFLAGS) -o $@ $<

bin:
	mkdir -p bin

clean:
	rm -rf bin
