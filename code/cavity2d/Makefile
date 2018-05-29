BIN = CAVITY_2D_STAGGERED

.PHONY: default
default: CAVITY_2D_STAGGERED

CAVITY_2D_STAGGERED: CAVITY_2D_STAGGERED.F90 Makefile
	gfortran -O3 $< -o $@

.PHONY: clean
clean:
	rm -rf $(BIN) *.mod

.PHONY: clean-run
clean-run:
	rm -rf $(BIN) *.DAT *.log
