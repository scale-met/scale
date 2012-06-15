mkpara.f90 is program to make "micpara.dat", which is needed to run the Spectral Bin Microphysical Model(SBM) implemented in SCALE3 model.
In the "micpara.dat" file, information of cloud microphysical parameters (e.g. radius of each cloud bin, capacitance, collision kernel or so) are described, and the SBM read the file before main calculation of SCALE3.

How to use?

1, Modify the Makefile
 Please modify the FC1, FC2, and FFLAGS in Makefile for your distribution

2, Setup the Number of hydrometeor's bin, and collection kernel
 Please open the mkpara.f90 by editor (vi or emacs), and modify the "nbin", "ncld", "ndrz", and "kphase". The meaning of these parameters area shonw below

 nbin : total number of hydrometeors bin (33 in default)
 ncld : number of cloud bin (it should be set to nbin/3) (11 in default)
 ndrz : number of drizzle bin (it should be set to 2*nbin/3) (22 in default)
 kphase : flag to select the collision kernel
          1 -> Hydro-dynamic kernel (shown in Iguchi et al. (2008) JGR) (default)
          2 -> Long kernel (Long 1973, JAS)
          3 -> Golovin Kernel (Golovin 1963)


3, Input "make" command

4, Run the binary named "mkpara".
 If the ascii file named micpara.dat is created, you succeed to create "micpara.dat"

5, The "micpara.dat"
