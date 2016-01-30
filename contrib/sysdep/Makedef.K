#
# ------ FOR K computer -----
#

##### additional compile option for all SDM modules

ADDITIONAL_FFLAGS_SDM =

##### additional compile option for specific SDM modules

# sdm_condensation_water.f90
#FFLAGS_SDM1 = -Kprefetch_indirect,simd=2
FFLAGS_SDM1 = -Kprefetch_indirect

# sdm_sd2fluid.f90
#FFLAGS_SDM2 = -Ksimd=2
FFLAGS_SDM2 = 

# sdm_coalescence.f90
#FFLAGS_SDM3 = -Ksimd=2
FFLAGS_SDM3 = 

# sdm_motion.f90
#FFLAGS_SDM4 = -Ksimd=2
FFLAGS_SDM4 = 

##### additional compile option for the MT (Mersenne Twister) module
FFLAGS_MT1 = 

