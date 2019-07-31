FC = mpiifort
QV_OPT =  -c -cpp -O2 -mcmodel=medium -D NCAR -D SINGLE -D TOROIDAL -D JUPITER -D PULSE -D RUN_FLUID # -D B1 -D PULSE # -D GYRO  -fPIC -g
QV_LOPT = 
FFLAGS = $(QV_OPT)


OBJECT = indeces.o param.o mesh.o fields.o parparam.o df.o distgen.o e1terms.o main.o \
	 grid_dipolar.o f1.o f2.o bc1.o bc2.o pred.o corr.o \
         Va.o Va_2.o density1.o density2.o initial_pulse.o initial_epsilon5.o \
	 lemks.o current.o initial_flr.o field_arrays.o initial_arrays.o \
         field_update.o particle_update.o rho_i.o e1_poloidal.o \
         epar1.o epar2.o e1.o eperp_gyro.o eperp.o ec_calc.o \
	 ec_calc_poloidal.o e1_term_write.o \
         full_write.o full_read.o input_read.o input_out.o \
         pcurrent0.o pcurrent1.o moment_bc.o \
         fdistrb.o fdistr2D.o df_write.o padistr2.o mudistr.o pusher.o \
         fluid_write.o particle_write.o \
         divb.o divj.o \
         fl_normalize.o p_normalize.o scaling.o \
         p_initialize.o sp_init.o v_initialize3.o pa_init.o \
         simp.o tridag.o filterd.o getname2.o getname_mpi.o \
         distgen_vx.o distgen_vy.o  distgen_vz.o distgen_vz2.o \
         j1delta_det.o  distgen_extra.o

OBJECT2 = indeces.o input_read.o df.o df_extract.o getname2.o param.o \
          parparam.o

%.o: %.f90 
	$(FC) $(FFLAGS) $<

mhd     : $(OBJECT)
	$(FC)  $(OBJECT) $(QV_LOPT) -o gke_f90.exe

df      : $(OBJECT2)
	$(FC)  $(OBJECT2) $(QV_LOPT) -o a.out

datareallyclean :
	rm d_* *.dat s_* t_*


dataclean	:
	rm d_u*  d_b* d_m* d_e* d_j* d_s* d_dn* d_rh* d_de* d_ge* d_dj* d_db* s_*

clean	:
	rm *.o *.mod input.out gke_f90.exe

