#C data file for simple example
#C file created using an r4ss function
#C file write time: 2026-04-28  10:37:28
#
1 #_styr
30 #_endyr
1 #_nseas
12 #_months_per_seas
2 #_Nsubseasons
1 #_spawn_month
2 #_Nsexes
32 #_Nages
1 #_N_areas
1 #_Nfleets
#_fleetinfo
#_type	surveytiming	area	units	need_catch_mult	fleetname
1	-1	1	1	0	Fishery	#_1
#_Catch data
#_year	season	fleet	catch	catch_se
 -999	1	1	1e-20	0.05	#_1         
    1	1	1	  100	0.05	#_2         
    2	1	1	  100	0.05	#_3         
    3	1	1	  100	0.05	#_4         
    4	1	1	  100	0.05	#_5         
    5	1	1	  100	0.05	#_6         
    6	1	1	  100	0.05	#_7         
    7	1	1	  100	0.05	#_8         
    8	1	1	  100	0.05	#_9         
    9	1	1	  100	0.05	#_10        
   10	1	1	  100	0.05	#_11        
   11	1	1	  100	0.05	#_12        
   12	1	1	  100	0.05	#_13        
   13	1	1	  100	0.05	#_14        
   14	1	1	  100	0.05	#_15        
   15	1	1	  100	0.05	#_16        
   16	1	1	  100	0.05	#_17        
   17	1	1	  100	0.05	#_18        
   18	1	1	  100	0.05	#_19        
   19	1	1	  100	0.05	#_20        
   20	1	1	  100	0.05	#_21        
   21	1	1	  100	0.05	#_22        
   22	1	1	  100	0.05	#_23        
   23	1	1	  100	0.05	#_24        
   24	1	1	  100	0.05	#_25        
   25	1	1	  100	0.05	#_26        
   26	1	1	  100	0.05	#_27        
   27	1	1	  100	0.05	#_28        
   28	1	1	  100	0.05	#_29        
   29	1	1	  100	0.05	#_30        
   30	1	1	  100	0.05	#_31        
-9999	0	0	    0	   0	#_terminator
#_CPUE_and_surveyabundance_observations
#_Units:  0=numbers; 1=biomass; 2=F; >=30 for special types
#_Errtype:  -1=normal; 0=lognormal; >0=T
#_SD_Report: 0=no sdreport; 1=enable sdreport
#_fleet	units	errtype	SD_report
1	1	0	0	#_Fishery
#
#_CPUE_data
#_year	month	index	obs	se_log
    1	6	1	50	0.3	#_1         
    2	6	1	50	0.3	#_2         
    3	6	1	50	0.3	#_3         
    4	6	1	50	0.3	#_4         
    5	6	1	50	0.3	#_5         
    6	6	1	50	0.3	#_6         
    7	6	1	50	0.3	#_7         
    8	6	1	50	0.3	#_8         
    9	6	1	50	0.3	#_9         
   10	6	1	50	0.3	#_10        
   11	6	1	50	0.3	#_11        
   12	6	1	50	0.3	#_12        
   13	6	1	50	0.3	#_13        
   14	6	1	50	0.3	#_14        
   15	6	1	50	0.3	#_15        
   16	6	1	50	0.3	#_16        
   17	6	1	50	0.3	#_17        
   18	6	1	50	0.3	#_18        
   19	6	1	50	0.3	#_19        
   20	6	1	50	0.3	#_20        
   21	6	1	50	0.3	#_21        
   22	6	1	50	0.3	#_22        
   23	6	1	50	0.3	#_23        
   24	6	1	50	0.3	#_24        
   25	6	1	50	0.3	#_25        
   26	6	1	50	0.3	#_26        
   27	6	1	50	0.3	#_27        
   28	6	1	50	0.3	#_28        
   29	6	1	50	0.3	#_29        
   30	6	1	50	0.3	#_30        
-9999	0	0	 0	  0	#_terminator
0 #_N_discard_fleets
#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)
#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal
#
#_discard_fleet_info
#
#_discard_data
#
#_meanbodywt
0 #_use_meanbodywt
 #_DF_for_meanbodywt_T-distribution_like
#
#_population_length_bins
3 # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector
43 #_N_lbinspop
#_lbin_vector_pop
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 #_lbin_vector_pop
1 #_use_lencomp
#
#_len_info
#_mintailcomp	addtocomp	combine_M_F	CompressBins	CompError	ParmSelect	minsamplesize
-1	0.001	0	0	0	0	0.001	#_Fishery
21 #_N_lbins
#_lbin_vector
12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 #_lbin_vector
#
#_lencomp
#_year	month	fleet	sex	part	Nsamp	f12	f13	f14	f15	f16	f17	f18	f19	f20	f21	f22	f23	f24	f25	f26	f27	f28	f29	f30	f31	f32	m12	m13	m14	m15	m16	m17	m18	m19	m20	m21	m22	m23	m24	m25	m26	m27	m28	m29	m30	m31	m32
    1	1	1	1	0	10	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	#_1         
-9999	0	0	0	0	 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	#_terminator
32 #_N_agebins
#
#_agebin_vector
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 #_agebin_vector
#
#_ageing_error
1 #_N_ageerror_definitions
#_	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1	   -1
0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001
#
#_age_info
#_mintailcomp	addtocomp	combine_M_F	CompressBins	CompError	ParmSelect	minsamplesize
-1	0.001	0	0	0	0	0.001	#_Fishery
1 #_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths
 #_combine males into females at or below this bin number
#_X.9999	X0	X0.1	X0.2	X0.3	X0.4	X0.5	X0.6	X0.7	X0.8	X0.9	X0.10	X0.11	X0.12	X0.13	X0.14	X0.15	X0.16	X0.17	X0.18	X0.19	X0.20	X0.21	X0.22	X0.23	X0.24	X0.25	X0.26	X0.27	X0.28	X0.29	X0.30	X0.31	X0.32	X0.33	X0.34	X0.35	X0.36	X0.37	X0.38	X0.39	X0.40	X0.41	X0.42	X0.43	X0.44	X0.45	X0.46	X0.47	X0.48	X0.49	X0.50	X0.51	X0.52	X0.53	X0.54	X0.55	X0.56	X0.57	X0.58	X0.59	X0.60	X0.61	X0.62	X0.63	X0.64	X0.65	X0.66	X0.67	X0.68	X0.69	X0.70	X0.71
-9999	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	#_terminator
#
#_MeanSize_at_Age_obs
0 #_use_MeanSize_at_Age_obs
0 #_N_environ_variables
0 #_N_sizefreq_methods
0 #_do_tags
0 #_morphcomp_data
0 #_use_selectivity_priors
#
999
