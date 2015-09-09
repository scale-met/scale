#This is created by "tools/depend.pl" at Wed Jan  8 09:02:23 2014

test_gadg.o:	test_gadg.F90 gadg_algorithm.o gadg_base.o gadg_macros.inc gadg_reduction.o

gadg_algorithm.o:	gadg_algorithm.F90 gadg_algorithm_accum.inc gadg_base.o gadg_config.inc gadg_config_platform.inc gadg_macros.inc

gadg_base.o:	gadg_base.F90 gadg_base.inc gadg_config.inc gadg_config_platform.inc

gadg_reduction.o:	gadg_reduction.F90 gadg_base.o gadg_reduction_detail.inc

gadg_constants.F90:	gadg_constants.F90_fpx 

gadg_constants.o:	gadg_constants.F90 

