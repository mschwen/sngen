#########################################################################
#  Makefile for SNO Supernova Generator: 'sngen'
#########################################################################

####################### Compilers #######################

# linux/Mac OS X computer
 FC =            /usr/bin/g77
# SUN Solaris machines on site ('surf', 'solstice')
# FC =            /opt/SUNWspro/bin/f77

################### Library/Code Locations #############

SNO_CODE =      /home/fermi/my_snoman/snoman/10_0000/code

##################### Compiler Flags ####################

#
# Mac OS X (chose either production or debugging)
#
# production flags
#FFLAGS =        -Wall -Os -fno-second-underscore -fno-f2c -I$(SNO_CODE)
# debugging flags
#FFLAGS =        -g -Wall -fno-second-underscore -fno-f2c -I$(SNO_CODE)

#
# Linux (chose either production or debugging; chose X86 flags in either case)
#
# production flags
FFLAGS =        -Wall -Os -mieee-fp -ffloat-store -fno-second-underscore -fno-f2c -I$(SNO_CODE)
# degugging flags
#FFLAGS =        -g -Wall -mieee-fp -ffloat-store -fno-second-underscore -fno-f2c -I$(SNO_CODE)
FFLAGS_X86 =    -malign-double -march=i686 

#
# SUN Solaris
#
# production flags
#FFLAGS =        -xO4 -dalign -native -xlibmil -I$(SNO_CODE)

################### Libraries ###################

# SNOMAN libraries
SNOLIB  =       $(SNO_CODE)/libsnolib.a

# CERNLIB libraries
#CERNLIB =       /usr/local/lib/libpacklib.a \
                /usr/local/lib/libmathlib.a

# Mac OS X linking errors (saveFP, rest FP, prebinding)
#MACOSX =        -noprebind /usr/lib/libgcc.a

#LIBS =          $(SNOLIB) $(CERNLIB) $(MACOSX)
# Linux CERNLIB libraries
CERNLIB =       $(CERNDIR)/libpacklib.a \
                $(CERNDIR)/libmathlib.a
LIBS =          $(SNOLIB) $(CERNLIB) -lnsl -L/usr/lib/gcc/i386-redhat-linux/4.1.1/ -lgfortran
# SUN Solaris CERNLIB libraries
#CERNLIB =       /usr/local/cern/pro/lib/libpacklib.a \
#                /usr/local/cern/pro/lib/libmathlib.a
#LIBS =          $(SNOLIB) $(CERNLIB) -lsocket -lnsl

############################################################################

SRCS =  sngen.f  \
           read_user.f  \
           sn_interactions.f  \
           flux_integral.f  \
               sn_flux.f  \
                 get_nu_type.f  \
                 flavour_swap.f  \
                   E_swap.f  \
                 get_flux_index.f  \
                 sn_num.f  \
                     sn_lum.f  \
                     sn_Eavg.f  \
                 fd_spectrum.f  \
                    fd_norm.f  \
           count_file_open.f  \
           count_sum.f  \
              read_file_flux.f \
              det_flux.f  \
                 get_int_type.f  \
                 get_xsect.f  \
                    file_open.f  \
                    get_target.f  \
                    xsect_nu_e_es.f  \
                    xsect_nu_p_es.f  \
                    get_threshold.f  \
                    xsect_nu_p_cc.f  \
                    xsect_nu_d_cc.f  \
                       aspline.f  \
                       aspline_eval.f  \
                    xsect_nu_12C_cc.f  \
                    xsect_nu_16O_cc.f  \
                    xsect_nu_208Pb_cc_1n.f  \
                    xsect_nu_208Pb_cc_2n.f  \
                    xsect_nu_d_nc.f  \
                    xsect_nu_16O_nc.f  \
                    xsect_nu_12C_nc.f  \
                    xsect_nu_208Pb_nc_1n.f  \
                    xsect_nu_208Pb_nc_2n.f  \
              central.f  \
              time_delay.f \
           poisson.f  \
              poidev.f  \
                 gammln.f  \
           sample.f  \
              file_check.f  \
              file_increment.f  \
              sample_nu.f  \
                 ran_val.f  \
                    ran3.f  \
              sample_part.f  \
                 sn_direction.f  \
                    ephem.f  \
                    celest2equator.f  \
                    equator2horizon.f  \
                    horizon2det.f  \
                 get_num_part.f  \
                 en_nu_p_es.f  \
                 dir_nu_p_es.f  \
                 en_nu_e_es.f  \
                    get_esdif_rad.f  \
                 dir_nu_e_es.f  \
                    combine_angles.f  \
                 part_cc.f  \
                    read_titles.f  \
                       sdb_open_all.f  \
                       sdb_read.f  \
                       brent2.f  \
                    en_nu_d_cc.f  \
                       xsect_nu_d_diff.f  \
                          f_xs_diff.f  \
                          f1dim.f  \
                    dir_nu_d_cc.f  \
                    en_nu_d_n_cc.f  \
                    dir_nu_d_n_cc.f  \
                    en_nu_p_cc.f  \
                       xsect_nu_p_diff.f  \
                    dir_nu_p_cc.f  \
                    en_nu_Pb_cc_1n.f  \
                    en_nu_Pb_cc_2n.f  \
                    en_nu_12C_cc.f  \
                    en_nu_16O_cc.f  \
                 part_nc.f  \
                    en_nu_d_nc.f  \
                    en_nu_12C_nc.f  \
                    en_nu_Pb_nc_1n.f  \
                    en_nu_Pb_nc_2n.f  \
                    en_nu_16O_nc_15O.f  \
                    en_nu_16O_nc_15N.f  \
                    dir_nu_d_nc.f  \
                    dir_ran.f \
           count_file_close.f  \
           summary.f  \
              group_sum1d.f  \
              group_isum1d.f  \
           stats.f  \
              stats_src.f  \
              group_sum2d.f \
              group_sum2d_file.f \
           sn_daq.f  \
           sn_output.f  \
              create_mcpl.f  \
              create_cmd.f  \
                get_medium.f  \
                file_cmd_type.f

OBJS = $(SRCS:%.f=%.o)

all : sngen

clean :
	rm -f $(OBJS)

nuclean :
	rm -f $ sam/nu*.out

partclean :
	rm -f $ sam/part*.out

samclean :
	rm -f $ sam/nu*.out
	rm -f $ sam/part*.out

snoclean :
	rm -f $ snoman/new/mcpl*.out
	rm -f $ snoman/new/sn*.cmd

outclean :
	rm -f $ lumen_nu_e.out
	rm -f $ lumen_nu_eb.out
	rm -f $ lumen_nu_mu.out
	rm -f $ fd_nu_e.out
	rm -f $ fd_nu_eb.out
	rm -f $ fd_nu_mu.out
	rm -f $ model_stats.out
	rm -f $ xsect/*
	rm -f $ sum_total.out
	rm -f $ countrate_all.out
	rm -f $ countrate_es.out
	rm -f $ countrate_cc.out
	rm -f $ countrate_nc.out
	rm -f $ countrate_nu_e.out
	rm -f $ countrate_nu_eb.out
	rm -f $ countrate_nu_mu.out
	rm -f $ countrate_d2o.out
	rm -f $ countrate_h2o.out
	rm -f $ countrate_av.out
	rm -f $ sam/*.out
	rm -f $ time_overlap.out
	rm -f $ snoman/new/mcpl*.out
	rm -f $ snoman/new/sn*.cmd

allclean :
	rm -f $ lumen_fit.out
	rm -f $ fd_shape_burst.out
	rm -f $ fd_shape_cooling.out
	rm -f $ model_stats.out
	rm -f $ xsect/*
	rm -f $ sum_total.out
	rm -f $ countrate_nu*.out
	rm -f $ sam/*
	rm -f $ time_overlap.out
	rm -f $ snoman/new/mcpl*.out
	rm -f $ snoman/new/sn*.cmd
	rm -f $ snoman/new/sn*.zdab

sngen : $(OBJS)
	$(FC) -o $@ $(OBJS) $(FFLAGS) $(FFLAGS_X86) $(LIBS)

$(OBJS) : %.o : %.f
	$(FC) $(FFLAGS) $(FFLAGS_X86) -c $*.f

# Include file (*.inc) dependencies

sngen.o : interaction_total.inc version.inc
read_user.o : file_flux_size.inc \
              physics_constants.inc trig_constants.inc convert_units.inc
sn_interactions.o : interaction_total.inc interaction_names.inc sn_classes.inc
count_sum.o : interaction_total.inc file_flux_size.inc convert_units.inc \
              time_size.inc
read_file_flux.o : file_flux_size.inc convert_units.inc
sn_flux.o : flux_size.inc file_flux_size.inc id_nu.inc convert_units.inc
get_nu_type.o : nu_class_size.inc id_nu.inc
get_flux_index.o : file_flux_size.inc
sn_num.o : id_nu.inc convert_units.inc
sn_lum.o : id_nu.inc convert_units.inc lumen_size.inc aspline_size.inc
sn_Eavg.o : id_nu.inc lumen_size.inc aspline_size.inc
fd_spectrum.o : id_nu.inc convert_units.inc
det_flux.o : trig_constants.inc convert_units.inc id_int.inc
get_int_type.o : sn_classes.inc id_int.inc
get_target.o : targets.inc
xsect_nu_e_es.o : energy_size.inc id_nu.inc \
                   physics_constants.inc trig_constants.inc convert_units.inc
det_flux_cc.o : xsect_size.inc
get_threshold.o : threshold.inc
xsect_nu_p_cc.o : energy_size.inc convert_units.inc
xsect_nu_d_cc.o : energy_size.inc xsect_size.inc aspline_size.inc \
                   convert_units.inc
aspline.o : aspline_size.inc
xsect_nu_12C_cc.o : energy_size.inc xsect_size.inc aspline_size.inc \
                     convert_units.inc
xsect_nu_O_cc.o : energy_size.inc convert_units.inc
det_flux_nc.o : xsect_size.inc
xsect_nu_d_nc.o : energy_size.inc xsect_size.inc aspline_size.inc \
                   convert_units.inc
xsect_nu_16O_nc.o : energy_size.inc convert_units.inc
xsect_nu_12C_nc.o : energy_size.inc xsect_size.inc aspline_size.inc \
                     convert_units.inc
time_delay.o : physics_constants.inc convert_units.inc id_nu.inc
poisson.o : random_seed.inc
sample_nu.o : interaction_total.inc sample_count_size.inc file_flux_size.inc
ran_val.o : random_seed.inc
sample_part.o : sample_count_size.inc id_int.inc
sn_direction.o : sno_coordinates.inc trig_constants.inc
celest2equator.o : trig_constants.inc
equator2horizon.o : trig_constants.inc
horizon2det.o : trig_constants.inc
get_num_part.o : id_int.inc id_part.inc
en_nu_e_es.o : id_nu.inc physics_constants.inc
get_esdif_rad.o : physics_constants.inc trig_constants.inc id_nu.inc
dir_nu_e_es.o : physics_constants.inc trig_constants.inc
part_cc.o : id_part.inc
#read_titles.o : z.inc mt.inc memory_size.inc sdb.inc sdb_ipath.inc \
#                sdb_options.inc bank_mctk.inc id_particles.inc \
#                id_en_types.inc titles_ncxs.inc titles_ncxb.inc \
#                titles_ccxs.inc titles_ccxb.inc titles_ccxt.inc \
#                titles_ncxt.inc
#sdb_open_all.o : db.inc db_io.inc db_misc.inc mzca.inc
#sdb_read.o : memory_size.inc zunit.inc mt.inc sdb.inc sdb_options.inc \
#             sdb_ipath.inc id_errors.inc sdb_contents_type.inc \
#             titles_dbhdr.inc
en_nu_d_cc.o : physics_constants.inc
xsect_nu_d_diff.o : physics_constants.inc
dir_nu_d_cc.o : trig_constants.inc
en_nu_d_n_cc.o : physics_constants.inc
dir_nu_d_n_cc.o : trig_constants.inc physics_constants.inc
en_nu_p_cc.o : id_part.inc physics_constants.inc
xsect_nu_p_diff.o : physics_constants.inc trig_constants.inc
dir_nu_p_cc.o : id_part.inc physics_constants.inc trig_constants.inc
part_nc.o : id_part.inc
en_nu_d_nc.o : physics_constants.inc
dir_nu_d_nc.o : trig_constants.inc
summary.o : interaction_total.inc \
            sn_classes.inc
stats.o : interaction_total.inc id_nu.inc time_size.inc
stats_src.o : id_nu.inc convert_units.inc
group_sum1d.o : interaction_total.inc
group_isum1d.o : interaction_total.inc
group_sum2d.o : interaction_total.inc sample_count_size.inc
group_sum2d_file.o : interaction_total.inc time_size.inc
sn_daq.o : interaction_total.inc sn_count_size.inc
sn_output.o : interaction_total.inc
create_mcpl.o : interaction_total.inc sample_count_size.inc convert_units.inc
create_cmd.o : interaction_total.inc id_part.inc id_cmd_type.inc
file_cmd_type.o : id_cmd_type.inc

# DO NOT DELETE
