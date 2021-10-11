! KED Fitting program July 17, 2021
! 	
!
!	Modified to read attenuation parameters from file "meked_attn_coeffs.csv"
!	Attenuation coefficient bias correction parametes now read from attn parameters file
!	Adjusted U and Pu attn coef to account
!	Added subroutine CALL AMU_VALUES() and modified function attn_fact_7_co to use the new data source
!	KED_FIT no longer specific to U, Pu, Np, Am and Cm - now accommodates 5 arbitrary elements defined
!	        from MEKED.exe  
!
!	Added interference correction for Pu and U X-rays
!	Added interference corrections for Bi and Pb X-rays
! 
!	modified to accommodate 4k spectra   							4/7/2021
!	revised HPGE detection efficiency    							4/12/2021
!	revised SS attenuation parameters    							4/12/2021
!	read misc attenuation parameters from file						4/13/2021
!	expanded container, solution and misc attenuation parameters from 8 to 10 coeffients	4/15/2021
!	Added trap for NaN in function FCHISQR1							4/21/2021
!	Fixed date truncation for *.;spe input files						4/29/2021
!	Corrected calculation for weig(i) in curfit						5/20/2021
!	Limit reference peak step tail, adjust W X-ray data					5/29/2021
!	Corrected reference peak escape width							6/01/2021
!
!	This software is intended as a research and development tool only.
!	
!	commandline call:   ked_fit_4k_v1p1 [nchans] [filename]
!	
!	nchans >2048 assumes a 4k spectrum, filename is the name of the analysis information file
!	the program will also require that the following files are present in the same directory
!		kfit_parms
!		KED_FIT_CONSTANTS  
!		MEKED_attn_coeffs
!
!	Questions: contact R.D.McElroy, Oak Ridge National Laboratory, Oak Ridge Tennessee
!
!
program KEDFIT_1
    integer re_i
    DOUBLE PRECISION temp_array(101), tails(4), e_escape(2), i_escape(2)
    DOUBLE PRECISION Y
!
    DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
    DOUBLE PRECISION kfit_errors(36), zdat(4100)
    DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
    DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100), BKG(4100)
    DOUBLE PRECISION YFIT(4100), YDAT(4100), ERR2(4100)
    DOUBLE PRECISION TAIL_AREA, TAIL_DECAY, FLAMDA, GATE
    DOUBLE PRECISION array(36, 36), BDUM(36), WEIGHT(4100)
    DOUBLE PRECISION deltaa(36), a1(36), AA(36), SIGMAA(36)
    DOUBLE PRECISION RAND_SPEC(4100), STEP_SIZE, ang_min
    DOUBLE PRECISION ENER_SCAT, ENER_ESC, ENERGY, b_scat, chi1, chi2
    DOUBLE PRECISION chisq_init, chisqr, D_E, deltae, dltchi, e0, e1
    DOUBLE PRECISION e_0, e_l1, e_l2, e_zero, i_first
    DOUBLE PRECISION t_real, rand_count, rand_scale, raw_count
    DOUBLE PRECISION sum1, sum_width, spec_fold2, x_chan, z
    DOUBLE PRECISION COVAR_ARRAY(36,36)
    DOUBLE PRECISION COVAR_TEMP1(36),covar_out1(36,36)
    DOUBLE PRECISION COVAR_TEMP2(36),covar_out2(36,36)
    DOUBLE PRECISION back_scat, functn2, FCHISQ1
    Double PRECISION cd_init, raw_start, cd_bkg
    Double PRECISION gaussian, convergence
    DOUBLE PRECISION amu_lib(10,8), cmu_lib(2,10), sol_lib(2,10), edges(10), misc_lib(6,10), misc_edges(6)
    DOUBLE PRECISION sol_edges(2), cont_edges(2), lor_widths(10), act_norm_factor(10)
    INTEGER i_kfit(36), idex_a(36), free_params(36), NTERMS
    INTEGER ichn, i_scat_hi, i_scat_lo, ichn_start, ichn_stop
    INTEGER idone, ichannels, ij, ispit, itmax, it_flag, j, k, l1, ijj
    INTEGER mtrms, mxp,i_cd_min, i_cd_max, n_cd_chan
    INTEGER ifil, lcmd, namflag
    INTEGER nchannels, nchan, nchn, niter, nmax, nmax2, nmin2, nt_max
    REAL E_LO1,E_LO2,E_HI1,E_HI2
!
    character string1*8
    character FILNAM1*60,SPECTEMP*60, newparname*60, Sample_ID*60, FILNAM2*60
    character inputfilename*60,  outfilname*100, specout2*100, results_rpt*100
    CHARACTER FILNAME*60, SPECOUT*60, PARMSOUT*60, named_spec_file*100
    CHARACTER arg*60 ,  seq_char*60
    CHARACTER(10) acq_date_mm, acq_date_dd, acq_date_yy
    CHARACTER(10) acq_time

    INTEGER i_conc(7), i_attn(7), chan_num, chan_num_in, seq_num
    INTEGER NPTS, num_NFREE, Nmin, MODE, nfree
    INTEGER NLOW, NHIGH, I, ILIVE, j_dex, istt,istp
    COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
    COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
    COMMON/KEDDAT3/YDAT, WEIGHT
    COMMON/attn_parms/amu_lib, cmu_lib, sol_lib, edges, sol_edges, cont_edges, lor_widths, act_norm_factor, misc_lib, misc_edges
	COMMON/filesize/mxp       
!
!     get max number of channels from the command line
!
	mxp=2048

  	CALL get_command_argument(1, arg)
!    	print *,' arg= ', arg
	READ(arg,*)chan_num_in               		  ! converts from character to a number
	if(chan_num_in.EQ.4096) mxp=4096

	CALL get_command_argument(2, arg)		  ! file name for input spectrum
	filnam2=trim(arg)

!	print *, "filname2 = " , filnam2, "=channels =  ", mxp
!
!
!	Fitting function free parameter list
!	paramters to be fit are designated by the free_parameters array 
!		free_parameter(n)=1 means the nth parameter is free
!		free_parameter(n)=0 means the nth parameter is fixed
!
!        kfit_parms(1) =Element 1  Concentration  (g/mL)
!        kfit_parms(2) =Element 2  Concentration  (g/mL)
!        kfit_parms(3) =Element 3  Concentration  (g/mL)
!        kfit_parms(4) =Element 4  Concentration  (g/mL)
!        kfit_parms(5) =Element 5  Concentration  (g/mL)
!        kfit_parms(6) =Matrix concentration  (g/mL)
!        kfit_parms(7) =Escape peak scaling factor
!        kfit_parms(8) =e_0 - X-ray generator end pt energy (keV)
!        kfit_parms(9) =x-ray generator intensity
!        kfit_parms(10) =Kramer distribution shape parameter (exponent)
!        kfit_parms(11) =detector back scatter intensity
!        kfit_parms(12) =Step Bkg Scale Factor
!        kfit_parms(13) = spectrum energy offset term (keV)
!        kfit_parms(14) =delta_e, spectrum energy cal (keV/channel)
!        kfit_parms(15) =detector resolution (keV) in terms of gaussian width (not FWHM) 
!        kfit_parms(16) =lorenztian width of Element 1 K-edge transition (keV)
!        kfit_parms(17) =lorenztian width of Element 2 K-edge transition (keV)
!        kfit_parms(18) =thickness of stainless steel beam filter (cm)
!        kfit_parms(19) =reference peak (Cd109) intensity
!        kfit_parms(20) =reference peak (Cd109) energy (keV)
!        kfit_parms(21) =Bi X-ray intensity
!        kfit_parms(22) =1st peak tail intenisty 
!        kfit_parms(23) =1st peak tail decay constant
!        kfit_parms(24) =2nd peak tail intenisty
!        kfit_parms(25) =2nd peak tail decay constant
!        kfit_parms(26) =Tungsten X-ray peak intensity
!        kfit_parms(27) =small angle scatter intensity 
!        kfit_parms(28) =small angle scatter bounding angle
!        kfit_parms(29) =Cd peak step background increment
!        kfit_parms(30) =U K-X-ray peak intensity
!        kfit_parms(31) =Pu  K-X-ray peak intensity
!        kfit_parms(32) =X-ray generator peak intensity
!        kfit_parms(33) =X-ray generator gaussian width (keV)
!        kfit_parms(34) =X-ray generator peak tail intensity
!        kfit_parms(35) =X-ray generator peak tail decay const
!        kfit_parms(36) =Pb X-ray Peak Intensity
!
!	KFIT contants array definitions
!        kfit_consts(1)	= cell length or vial diameter (cm)
!        kfit_consts(2)	= sample vial single wall thickness (cm) 
!        kfit_consts(3)	= X-ray generator Be window thickness (cm)
!        kfit_consts(4)	= cd filter thickness (cm)
!        kfit_consts(5)	= HPGe efficiency correction reference energy (keV)
!        kfit_consts(6)	= Kramer distribution shape parameter(quadratic exponent)
!        kfit_consts(7)	= X-ray Generator Incidence Angle (degrees)
!        kfit_consts(8)	= X-ray Generator Takeoff Angle (degrees)
!        kfit_consts(9) = Tungsten mass attenuation coefficenct modifier 1st term
!        kfit_consts(10) = TUngsten mass attenuation coefficenct modifier 2nd term
!        kfit_consts(11) = detector backscatter mininum scatter angle (degrees)
!        kfit_consts(12) = use fitted background correction (1=yes, 0 =no)
!        kfit_consts(13) = Include X-ray generator enpoint peak
!        kfit_consts(14) = generator off-set (keV)
!        kfit_consts(15) = not used
!        kfit_consts(16) = Vial type (glass = 1, poly = anything else)	
!        kfit_consts(17) = not used	
!        kfit_consts(18) = not used	
!        kfit_consts(19) = not used	
!        kfit_consts(20) = not used	
!        kfit_consts(21) = Fitting Region lower energy limit (keV)
!        kfit_consts(22) = Fitting Region upper energy limit (keV)
!        kfit_consts(23) = Fitting Regio lower energy limit (keV)	
!        kfit_consts(24) = Maximum number of iterations for fit
!        kfit_consts(25) = allow fit to continue afer max iterations	
!        kfit_consts(26) = override kfit_constants path length with value from data file	
!        kfit_consts(27) = not used	
!        kfit_consts(28) = not used	
!        kfit_consts(29) = not used	
!        kfit_consts(30) = Background Low Energy ROI lower limit (keV)
!        kfit_consts(31) = Background Low Energy ROI upper limit (keV)
!        kfit_consts(32) = Background High Energy ROI lower limit (keV)
!        kfit_consts(33) = Background High Energy ROI upper limit (keV)
!        kfit_consts(34) = Background Low Energy ROI Cd correction factor
!        kfit_consts(35) = Background Low Energy ROI Cd Scaling factor 
!        kfit_consts(36) = fractional change in fitting parameters for fderiv()
!
!        
!    	read in attenuation coefficients from file "MEKED_attn_coeffs.txt"

	CALL read_amu_values()
!
!
!	--------------------------------------------------------------------------------------------
!
!	Gettrial KED parameters from FIT_PARMS.TXT
!
!	--------------------------------------------------------------------------------------------
!
!
!    FILNAME='kfit_parms_##.txt'								! kfit_parms.txt contains initial seed values and use flag for each
!
    FILNAME="KFIT_PARMS_" // FILNAM2(10:len(TRIM(FILNAM2))) // ".txt"				! grab index number from spectrum file to get correct seed parameters
!
!	print *,  'test ', "KFIT_PARMS_" // FILNAM2(10:len(TRIM(FILNAM2))) // ".txt", filname
!
    print *, "KED Fit 4k Ver. 1.1 (Build July 17, 2021)"
!
    open (unit=10,file=FILNAME)
!
    read (10,*)FILNAM1										! filnam1 = root name for the output files read from 1st line of kfit_parms.txt
    do 1 l1 = 1, 36
    read (10,*)i,kfit_parms(l1), free_params(l1)
1   continue
    print *,' end read '
    close (unit=10)

    FILNAM1=FILNAM2         									!  spectrum file name from the command line
    print *, ' file name = ', FILNAM1
!
!	Create file names for use during fit
!	All files except for the results file (*.res) are required to be within the same directory as the executable
!	so no prefixes are used
!
    SPECOUT=TRIM(TRIM(FILNAM1)//".SPC")			!  final output text file containing raw spec and fitted spectrum
    SPECTEMP=TRIM(TRIM(FILNAM1)//".TMP")		!  temp output text file containing raw spec and fitted spectrum
    PARMSOUT=TRIM(TRIM(FILNAM1)//".PAR")		!  output file including fit parameters uncertainties and covariances after iteratation
    FILNAME=TRIM(TRIM(FILNAM1)// ".txt")		!  input spectrum text file
    newparname=TRIM(TRIM(FILNAM1)// ".rpt")		!  output text file containing final fit parameters uncertainties and covariance matrix
    PRINT*, FILNAME, SPECTEMP
!
!	--------------------------------------------------------------------------------------------
!
!		READ IN FITTING CONSTANTS FROM KED_FIT_CONSTANTS.DAT
!
!	--------------------------------------------------------------------------------------------
!
    open (unit=10,file="KED_FIT_CONSTANTS.txt")
!
    do 999 l1 = 1, 36
    read (10,*)i,kfit_consts(l1)
999 continue
    print *,' end read constants'
    close (unit=10)
!	
	convergence= kfit_consts(23)
!
!	DEFINE FIT REGION
!
    E_L1=kfit_consts(21)
    E_L2=kfit_consts(22)
    ichn_start=int((e_l1-kfit_parms(13))/kfit_parms(14))
    ichn_stop=int((e_l2-kfit_parms(13))/kfit_parms(14))
    NMIN=ichn_start
    NMAX=ICHN_STOP
!
!    mxp=4096
    MTRMS=36
    nterms = 0
!	nimn=ichn_start
    nmax=ichn_stop
    nchannels = ichn_stop - ichn_start + 1
    nt_max = 36
    DO 2 i = 1 , 36
    i_kfit(i) = 0
    idex_a(i) = 0
    a1(i) = 0
    kfit_errors(i) = 0
2   CONTINUE
!
    j = 0
    DO 3 i = 1 , nt_max
    If (free_params(i).EQ.1) j = j + 1
    If (free_params(i).EQ.1) nterms = nterms + 1
    If (free_params(i).EQ.1) i_kfit(i) = j
!	 position of kfit in a e.g. U only  nitric conc kfit_param (6) -> a1(2) so i_kfit(6)=2
    If (free_params(i).EQ.1) idex_a(j) = i
!	 position of a in kfit e.g. U only  nitric conc a1(2)-> kfit_param(6) so idex_a(2)=6
3   CONTINUE
    DO 4 I = 1 , NTERMS
    PRINT *, i, IDEX_A(I)
4   CONTINUE
!
!	SETUP DATA ARRAYS
!
    DO 5 I = 1,7
    I_CONC(I)=1
    I_ATTN(I)=1
5   CONTINUE
!
    DO 11 I = 1 , MTRMS
    A1(I)=0.
    AA(I)=0.
    DELTAA(I)=0.
    SIGMAA(I)=0.
    BDUM(I)=0.
    DO 10 J= 1, MTRMS
    ARRAY(I,J)=0.0
10  CONTINUE
11  CONTINUE
!
!	********************************************
!
    DO 15 I=1,NT_MAX
    PRINT *,I, KFIT_PARMS(I)
15  CONTINUE
!	************ Initialize arrays  ************
    DO 20 I = 1,4096
    SCAT_DAT(I)=0.
    RAW_DATA(I)=0.
    RAND_TEMP=0.
    XDAT(I)=0.
    weight(i)=0.
    ERROR1(I)=1.
20  CONTINUE
!
    E_0=kfit_parms(13)
    D_E=kfit_parms(14)

!
!	--------------------------------------------------------------------------------------------
!
!		READ IN SPECTRUM FROM DATA FILE:  SPEC_FIT_##.TXT
!
!	--------------------------------------------------------------------------------------------
!
!  	inputfilename = source of input data - for bookkeeping purposes
!	outfilname = full name for results file  e.g. c:\ked_temp\results\meked_file1.res
!
    open (unit=10,file=FILNAME)	
    read (10,*)inputfilename, outfilname, Sample_ID, acq_date_mm, acq_time	
    read (10,*)ichn, ilive, t_real, vial_dia, sample_temp, U_enrich, pu_wgt
    print *, ichn, ilive,t_real,  vial_dia, sample_temp, U_enrich, pu_wgt
    print *, ' data file name = ', inputfilename
    print *, ' Sample Id = ', Sample_ID
    print *, ' acq date = ' , acq_date_mm 
    print *, ' acq time = ' , acq_time
    print *, ' ***************************************** '


!       may want to use e0 and e1 as starting values for e_cal
!
    read (10,*)i,x_chan
    do 25 l1 = 1, mxp-1
    read (10,*)i,x_chan
!	IF (L1.GT.ICHN_STOP) GOTO 25
!	IF (L1.LT.ICHN_START) GOTO 25
    raw_data(l1)=x_chan
    IF (RAW_DATA(l1).GT.0) ERROR1(l1)=SQRT(RAW_DATA(l1))
    If (error1(l1).NE.0) weight(i) = 1.0 /error1(l1)**2
25  continue
    print *,' end read data file'
    close (unit=10)
    NCHN=ICHN_STOP-ICHN_START+1
    NMAX=NCHN

    if((vial_dia.GE.0).and.(kfit_consts(26).EQ.1)) kfit_consts(1)=vial_dia
!
!       ********** end data file read **********
!
!	*******   CORRECT FOR RANDOM SUMMING   *******
!
!	THESE PARAMETERS TO GO INTO A SETUP FILE
    NLOW=50
    NHIGH=2000
    GATE= 0.00000078813904

	istt=int((kfit_parms(8)+ 4-kfit_parms(13))/kfit_parms(14))     ! start rand correction background 4 keV above X-ray HV setting

!	istt=1689
	istp=2000
	if(mxp.GT. 2050) istt=int(1689*mxp/2048)
	if(mxp.GT. 2050) istp=int(2000*mxp/2048)
	if(mxp.GT. 2050) NHIGH=4000

!
    DO 50 I = 2, istp
    RAND_TEMP(I)=spec_fold2(NLOW, NHIGH, I, GATE, ILIVE, Raw_data)
50  CONTINUE
!
    tail_area=0.0
    tail_decay=1.0
    sum_width=0.5
    DO 60 I = 2 , istp
    ENERGY=E0+E1*FLOAT(I)
!   XDAT(I)=spec_tail(ENERGY,tail_area,tail_decay,E1)+RAND_TEMP(I)
    XDAT(I)=RAND_TEMP(I)
60  CONTINUE
!
    RAW_COUNT=0.0
    RAND_COUNT=0.0

    DO 65 I = istt , istp
    RAW_COUNT=RAW_COUNT+RAW_DATA(I)	!  sum of counts in rand corr ROI 
    RAND_COUNT=RAND_COUNT+XDAT(I)	! sum of estimated counts in rand corr ROI
65  CONTINUE
    RAND_SCALE=(RAW_COUNT-kfit_consts(29)*(istp-istt+1))/RAND_COUNT			! scale random bkg over the ROI including a fixed background component
    DO 70 I = 1, istp
    XDAT(I)=RAW_DATA(I)-XDAT(I)*RAND_SCALE-kfit_consts(29)				! subtract  random rate and a fixed background (from constants file)
    ERROR1(I)=(ERROR1(I)**2+(RAND_TEMP(I)*RAND_SCALE)**2/ABS(RAND_COUNT))**0.5
    If (error1(i).NE.0) weight(i) = 1.0 /error1(i)**2
70  CONTINUE
!
!   XDAT is now the random summing corrected spectrum
!
!	********** STRIP OUT STEP BACKGROUND   ******
!   estimate Cd bkg correction
!   assume Cd109 peak is near channel 983

	i_cd_min=975
	if(mxp.gt.2050) i_cd_min=int(975*(mxp/2048))
	i_cd_max=i_cd_min
	raw_start=0

    DO 72 i=i_cd_min, i_cd_min + int(20*(mxp/2048))
        if (XDAT(i).gt. raw_start) i_cd_max=i
        if (XDAT(i).gt. raw_start) raw_start=XDAT(i)
72  continue
    n_cd_chan=int(.5+2.5*2.35*kfit_parms(15)/kfit_parms(14))    ! +/- 2.5 FHWM from max channel
    cd_init=0

    do 73 i= i_cd_max - n_cd_chan, i_cd_max + n_cd_chan
    cd_init=cd_init+xdat(i)
73  continue
    cd_bkg=xdat(i_cd_max - n_cd_chan+1)+xdat(i_cd_max + n_cd_chan+1)
    cd_init=cd_init-cd_bkg*n_cd_chan

    Print *, 'Cd initial chan = ',i_cd_max, ' count =  ', XDAT(i_cd_max), ' mxp =  ',mxp
    Print *, 'Cd initial area = ',cd_init
!
    E_LO1=KFIT_CONSTS(30)
    E_LO2=KFIT_CONSTS(31)
    E_HI1=KFIT_CONSTS(32)
    E_HI2=KFIT_CONSTS(33)
!
	do 770 I = 1, mxp
	BKG(I)=0.0
770     continue
        If (kfit_consts(12).EQ.1) goto 774
    CALL BKG_STRIP(KFIT_PARMS, BKG, cd_init)
!	CALL BKG_STRIP2(step_size, BKG)
    DO 74 I=1 , mxp-1
    XDAT(I)=XDAT(I)-BKG(I)
74  CONTINUE
774       continue
!
!	XDAT IS NOW THE RANDOM AND BKG CORRECTED SPECTRUM
!
!	******** create back scatter array  ************
!
    ang_min = kfit_consts(11)
    i_scat_lo = Int((60 - kfit_parms(13))/kfit_parms(14))
    i_scat_hi = Int((130 - kfit_parms(13))/kfit_parms(14))
!
    do 75 i = i_scat_lo , i_scat_hi
    ENER_SCAT = kfit_parms(13) + kfit_parms(14) * i
    scat_dat(i)=back_scat(ener_SCAT,i,kfit_parms)
!   IF(ENER_00.LT.100) PRINT*, I, ENER_SCAT, SCAT_DAT(I)
75  continue
!
!    ******************************************
!
    PRINT *, 'REACHED CONDENSE'
!	******* CONDENSE MATRIX  **********
    DO 100 i = 1 , nterms
    a1(i) = kfit_parms(idex_a(i))
100 CONTINUE
!
    DO 110 i = 1 , 4
    tails(i) = abs(kfit_parms(21 + i))
110 CONTINUE
    e_zero = kfit_parms(13)
    deltae = kfit_parms(14)
!
    DO 120 i = 1 , nterms
    deltaa(i) = a1(i) * kfit_consts(36)
    If (deltaa(i).EQ.0) deltaa(i) = 0.0001
120 CONTINUE
!
    nmin2 = ICHN_START
    nmax2 = ICHN_STOP
!
    npts = 1 + nmax2 - nmin2
    Sum1 = 0.
    DO 130 i = nmin2 , nmax2
    Sum1 = Sum1 + xdat(i)
130 CONTINUE
    DO 140 i = 1 , mxp
    ERR2(i) = error1(i)
    ydat(i) = xdat(i)
140 CONTINUE
!
    itmax = INT(KFIT_CONSTS(24))
    it_flag=KFIT_CONSTS(25)
    i_first=0
190 continue
    chi1 = 0.
    niter = 0
    MODE = 1
!
!	START FITTING  ****************
!
    FLAMDA = 0.001
    PRINT *, 'REACHED FITTING'
!
!
    do 200 i = 1 , mxp
    ZDAT(i) = e_zero + float(i) * deltae
200 continue
!
    do 201 i=ichn_start, ichn_stop
    yfit(i)=functn2(i,kfit_parms)
!   PRINT *, I, YFIT(I)
201 continue
    num_free = npts - nterms
!  ydat is the stripped data array, yfit is the calculated response
    chisq_init = FCHISQ1(ydat, weight, npts, num_free, NMIN, yfit)
    print *, 'chisq_init = ', chisq_init
!
    open (unit=20,file=SPECTEMP,action="write",status="replace")
    write (20,*) "file:", FILNAME
    write (20,*) "CHAN  RAW_DATA    RAW_ERROR   Corrected   fit result"
    do 202 I = 1, mxp
    ENERGY=E_0+D_E*I
    z=functn2(i,kfit_parms)
!       TEMP=functn_kfit(energy, KFIT_PARMS)
    write (20,*) i,RAW_datA(i),ERROR1(I),XDAT(I), Z
202 continue
    close (unit=20)
!
210 continue
    do 5555 i = 1 , itmax
    do 220 ij = 1 , nterms
    deltaa(i) = a1(i) * kfit_consts(36)
220 continue
    niter = 1 + niter
!
    print *, FILNAME, SPECTEMP
    print *, 'iteration = ', niter
    print *, 'nterms = ', nterms
!
!
    FLAMDA = 0.001
    MODE=1
    PRINT *,'REACHED CURFIT'
!
    Call curfit(ydat, ERR2, npts, nterms, MODE, a1, kfit_parms,deltaa, SIGMAA, &
    FLAMDA, yfit, CHISQR, nmin2, free_params,Idex_a, COVAR_ARRAY)
!
    Call a_to_kfit(SIGMAA, kfit_errors, idex_a, nterms)
    PRINT *, 'Nterms after= ', nterms
	do 245 ijj=15, 17
	kfit_parms(ijj)=abs(kfit_parms(ijj))
245 continue
	do 246 ijj= 22, 25
	kfit_parms(ijj)=abs(kfit_parms(ijj))
246 continue
    CHI2 = CHISQR
    PRINT *, 'chi2 = ', CHI2, NMIN2
    DO 250 ispit = 1 , nt_max
    PRINT *, 'a1(',ispit ,')=',kfit_parms(ispit),kfit_errors(ispit)
250 CONTINUE
    DLTCHI = Abs(CHI2 - chi1) / CHI2
    chi1 = CHI2
    PRINT *, 'chi_square = ', CHI1

	if(convergence.LE.0) convergence=0.0001

    If (DLTCHI.LT.convergence) GoTo 5600
5555    continue
    STRING1="N"
    print*,' CHI-SQUARE DID NOT CONVERGE AFTER ',ITMAX,' ITERATIONS'
        if(it_flag.NE.0) print *,' CONTINUE (Y,N)?'
        if(it_flag.NE.0) read *, STRING1
	if(it_flag.EQ.0) STRING1="N"
    IF(STRING1.NE. 'Y')GOTO 5600
    IF(I_FIRST.NE.0) GOTO 5600
    print *,' MAX NUMBER OF ITERATIONS?'
    read *,ITMAX
    GoTo 210
!
5600    continue
    DO 5650 i = 1 , nterms
    AA(i) = a1(i)
5650    CONTINUE
    DO 5700 i = 1 , nterms
    SIGMAA(i + 1) = SIGMAA(i + 1)
5700    CONTINUE
!
    PRINT *, 'chi_square = ', CHI2
!
1001    continue
!
!		****************
!
!	EXPAND COVARIANCE ARRAY
!
    DO 6010 I = 1 , 36
    DO 6000 J = 1 , 36
    COVAR_TEMP1(J)=COVAR_ARRAY(I,J)
6000    CONTINUE
!
    Call a_to_kfit(COVAR_TEMP1, COVAR_TEMP2, idex_a, nterms)
    DO 6010 K=1, 36
    COVAR_OUT1(I,K)=COVAR_TEMP2(K)
6010    CONTINUE
!
    DO 6030 I = 1 , 36
    DO 6020 J = 1 , 36
    COVAR_TEMP1(J)=COVAR_OUT1(J,I)
6020    CONTINUE
    Call a_to_kfit(COVAR_TEMP1, COVAR_TEMP2, idex_a, nterms)
!
    DO 6030 K=1, 36
    COVAR_OUT1(K,I)=COVAR_TEMP2(K)
6030    CONTINUE
!
!      write parameters to par file
!
    open (unit=20,file=PARMSOUT,action="write",status="replace")
    write (20,*) filname
    write (20,*) CHI2
    do 6100 I = 1, 36
    write (20,*) i,KFIT_PARMS(I),kfit_errors(I)
6100    continue
    DO 6110 I = 1, 36
    DO 6110 J = 1, 36
    write (20,*) i,J, COVAR_OUT1(I,J)
6110    continue
    close (unit=20)
!
!	***********************************
!
    open (unit=20,file=SPECOUT,action="write",status="replace")
!
!    write (20,*) "file:", FILNAME, FILNAME, inputfilename


    write (20,*) "file: ", TRIM(FILNAME), " ", TRIM(inputfilename) , " ", TRIM(Sample_ID)
    write (20,*) "CHAN RAW_DATA RAW_ERROR Corrected fit"
    do 325 I = 1, mxp
!
    ENERGY=E_0+D_E*I
    B_SCAT=SCAT_DAT(I)
!	z=functn2(i,kfit_parms)+B_SCAT
    z=functn2(i,kfit_parms)
    write (20,*) i,RAW_DATA(I),ERROR1(I), XDAT(I), Z
325 continue
    close (unit=20)
    i_first = i_first+1
    print *,'looping back'
    itmax=1
!	if (I_first .eq. 1) goto 190
!
!      calc and write new error parameters to results file
!
        call Est_error(kfit_parms, kfit_errors, j_dex)
!
!	***********************************
!	***********************************
!
	SPECOUT2=TRIM(TRIM(outfilname)//".res")

    open (unit=20,file=specout2,action="write",status="replace")
!

    write (20,*) TRIM(FILNAME), " ", TRIM(inputfilename) , " ", TRIM(Sample_ID), " ", TRIM(acq_date_mm) , " ", acq_time
    write (20,*) CHI2, vial_dia, sample_temp, U_enrich, pu_wgt, ilive, t_real
    do 455 i = 1, 36
    write (20,*) i,kfit_parms(i), kfit_errors(i), free_params(i)
455 continue
    DO 460 I = 1, 36
    DO 460 J = 1, 36
    COVAR_OUT1(I,J)= COVAR_OUT1(I,J)*(float(J_dex))**2
    write (20,*) i,J, COVAR_OUT1(I,J)
460    continue
   write (20,*) "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
!
    write (20,*) "file:", FILNAME, FILNAME, inputfilename
    write (20,*) "CHAN RAW_DATA RAW_ERROR Corrected fit"
    do 326 I = 1, mxp
!
    ENERGY=E_0+D_E*I
    B_SCAT=SCAT_DAT(I)
!	z=functn2(i,kfit_parms)+B_SCAT
    z=functn2(i,kfit_parms)
    write (20,*) i,RAW_DATA(I),ERROR1(I), XDAT(I), Z
326 continue
    close (unit=20)

!	
!
!	***********************************
!	Write final parameters to file
!
!
    open (unit=20,file=newparNAME,action="write",status="replace")
!
    write (20,*) FILNAME, inputfilename
    write (20,*) CHI2, vial_dia, sample_temp, U_enrich, pu_wgt, j_dex
    do 425 i = 1, 36
    write (20,*) i,kfit_parms(i), kfit_errors(i)
425 continue
    DO 430 I = 1, 36
    DO 430 J = 1, 36
    COVAR_OUT1(I,J)= COVAR_OUT1(I,J)*(float(J_dex))**2
    write (20,*) i,J, COVAR_OUT1(I,J)
430    continue
    close (unit=20)
!
!	***********************************
!
        goto 441
    SPECOUT2=TRIM(TRIM(outfilname)//".rpt")
!
	    open (unit=20,file=specout2,action="write",status="replace")
!
	    write (20,*) FILNAME, inputfilename
	    write (20,*) CHI2, vial_dia, sample_temp, U_enrich, pu_wgt, j_dex
	    do 435 i = 1, 36
	    write (20,*) i,kfit_parms(i), kfit_errors(i)
435 continue
	    DO 440 I = 1, 36
	    DO 440 J = 1, 36
!	    COVAR_OUT1(I,J)= COVAR_OUT1(I,J)*(float(J_dex))**2
    write (20,*) i,J, COVAR_OUT1(I,J)
440    continue
    close (unit=20)
441      continue

    END
!
!	*********************************
!	*********************************
!	*********************************
!
    SUBROUTINE a_to_kfit(a1, kfit_parms, idex_a, nterms)
    DOUBLE PRECISION kfit_parms(36), A1(36)
    INTEGER idex_a(36),  NTERMS
!	Stuffs variable parameters from a1(i) back into kfit_parms(i)
    DO 10 i = 1 , nterms
    kfit_parms(idex_a(i)) = a1(i)
10  CONTINUE
    RETURN
    END
!
!	*********************************
!
    Function FCHISQ1(y, weight, npts, nfree, NMN1, yfit)
    DOUBLE PRECISION YFIT(4100), Y(4100), WEIGHT(4100), fchisq1
    INTEGER NPTS, NFREE, NMN1
!	COMMON/CALCDAT/YFIT
    AFREE = FLOAT(nfree)
    NMN2 = NMN1 + npts - 1
    CHISQR = 0.
!    If (nfree.LE.0) Return
    DO 30 i = NMN1 , NMN2
!	PRINT *, I, Y(I), YFIT(I), WEIGHT(I)
    CHISQR = CHISQR + weight(i) * (y(i) - yfit(i))**2
30  CONTINUE
    FCHISQ1 = CHISQR / AFREE
    print *, 'fchis', npts, nfree, nmn1
!	PRINT *, FREE, WEIGHT(1200), FCHISQ1
    RETURN
    END
!
!	**********************************************************
!	FUNCTIONS THAT DECRIBE THE KED RESPONSE
!	**********************************************************
!
    Function functn2(chan_num, kfit_parms)
!       * * * KEDFUNC * * *
!
!	KED Response function
!
!	Ported from KED VBA routine
!	April 8, 2016
!
    DOUBLE PRECISION functn2, UPUXrays, Pb_Xrays, Bi_xrays
    DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
    DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
    DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100)
    INTEGER i_conc(7), i_attn(7), chan_num, CHAN_1
    COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
    COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
!
    DOUBLE PRECISION e_start, deltae, gauss_width
    DOUBLE PRECISION p_i, p_e, ref_e, escape_a, p_w
    DOUBLE PRECISION temp_array(201),e_escape(2)
    DOUBLE PRECISION ener_00, y, ener_0,i_escape(2)
    DOUBLE PRECISION gaussian,functn_kfit,E_T, GAUS, FUNC_K, ESC_1, ESC_2
    DOUBLE PRECISION bent_kramers4, Y_J, T_B, ENER_11, TEMP, TEMP_VAL
    DOUBLE PRECISION W_1,W_2, SCAT_ANGLE, E_MIN, CD_PEAK
    DOUBLE PRECISION scat_p1, scat_p2, scat_i, ANGLES
    DOUBLE PRECISION SCAT_E, SS_THICK, SCATTER, ANG_MIN
    DOUBLE PRECISION ENER10, ENER11, TAIL1, TAIL2
    DOUBLE PRECISION lo_tail_area,lo_tail_decay
    DOUBLE PRECISION fast_tail_area,fast_tail_decay
    DOUBLE PRECISION STEP_SIZE, BKG(4100), ENERGY
    DOUBLE PRECISION UPU_K_alpha, Pb_K_alpha, Bi_K_alpha
    DOUBLE PRECISION x_AREA, E_0, TAILCD, CD_TAIL_I, CD_TAIL_T
    DOUBLE PRECISION gen_width, gen_e, gen_dum, gen_int, gen_area
    DOUBLE PRECISION ener1, HPGE_eff	
!
!	calls   UPUXrays, Pb_K_alpha, Bi_K_alpha, gaussian
!		FUNCTN_KFIT, bent_kramers4, g_smooth, BKG_STRIP2, W_Voigt
!

    chan_1 = chan_num
!    E_0=kfit_parms(8)
    n=51
    e_start = kfit_parms(13)
    deltae = kfit_parms(14)
    Gauss_width = kfit_parms(15)
    p_i = kfit_parms(19)
    p_e = kfit_parms(20)
!	tail_long = kfit_parms(21)
    ref_e = kfit_consts(5)
    p_w = Gauss_width * (p_e / ref_e)**0.5
    escape_a = kfit_parms(7)
    e_escape(1) = 9.88
    e_escape(2) = 10.99
    i_escape(1) = 1.0
    i_escape(2) = 1.0
!
!	Calculate counts from X-ray peaks from U, Pu, Pb and Bi peaks 
!	includes K-alpha and K-beta lines
!	intended as an inteference correction i.e. non quantitative
!
        UPU_K_alpha =UPUXrays(chan_num, kfit_parms)    !  calc characteristic X-ray peak area
        Pb_K_alpha = Pb_Xrays(chan_num, kfit_parms)    !  calc characteristic X-ray peak area
        Bi_K_alpha = Bi_Xrays(chan_num, kfit_parms)    !  calc characteristic X-ray peak area
!
    slo_tail_area = abs(kfit_parms(22))
    slo_tail_decay = abs(kfit_parms(23))
    fast_tail_area = abs(kfit_parms(24))
    fast_tail_decay = abs(kfit_parms(25))
!
    ener_00 = e_start + FLOAT(chan_1) * deltae
    y = 0
!
!	create bremsstrahlung array, TEMP_ARRAY, for smoothing
!
    do 10  j = 1 , n
    ener_0 = ener_00 + float(j - (int(n/2)+1)) * deltae
!	escape peaks from continuum and 109Cd peak
!
    E_T = ener_0 + e_escape(1) + 0.10
    GAUS=gaussian(E_T, p_e, p_i, p_w*0.1, deltae)
    FUNC_K=FUNCTN_KFIT(E_T, KFIT_PARMS)
    esc_1=escape_area(E_T, 1)*escape_a*(FUNC_K+3.0*GAUS)
!
    E_T = ener_0 + e_escape(2)
    GAUS=gaussian(E_T, p_e, p_i, p_w*0.1, deltae)
    FUNC_K=FUNCTN_KFIT(E_T, KFIT_PARMS)
    esc_2=escape_area(E_T, 2)*escape_a*(FUNC_K+3.0*GAUS)
!
    T_KFIT=FUNCTN_KFIT(ener_0, KFIT_PARMS)
    T_B= bent_kramers4(1, ener_0,1.0, KFIT_PARMS)
    temp_array(j) = T_KFIT + (esc_1+esc_2)+T_B
!
10  CONTINUE
!
!	add tailing
!
    DO 20 j = 1 , 100
    ener10 = deltae * float(j) + ener_00
    ener11 = deltae * float(j)
    TEMP = FUNCTN_KFIT(ener10, KFIT_PARMS)
    TEMP_1=gaussian(ener10, p_e, p_i, p_w, deltae)
!
    TAIL1=slo_tail_area*EXP(-slo_tail_decay*ener11)
    TAIL2=fast_tail_area*EXP(-fast_tail_decay*ener11)
!    TAILCD=CD_tail_i*EXP(-CD_tail_T*ener11)        
    TAILCD=0.
!
    y = y + temp * (TAIL1 + TAIL2)+TEMP_1*TAILCD
!
20  CONTINUE
!
    temp_val=g_smooth(ener_00,deltae,Gauss_width,temp_array,n)
!
    ref_e = kfit_consts(5)
    w_1 = kfit_consts(9)
    w_2 = kfit_consts(10)
!    scat_angle = kfit_consts(11)
!    e_min = kfit_consts(12)
    scat_p1 = kfit_consts(13)
    scat_p2 = kfit_consts(14)
    scat_i = abs(kfit_parms(11))
!    angles = kfit_parms(12)
	scat_e = kfit_parms(8)
    ss_thick = kfit_parms(18)
!
    scatter = 0.0
    ang_min = kfit_consts(11)
!
!	scatter = scat_i * back_scat(ener00, chan_num, kfit_parms)
    scatter = scat_i * scat_dat(chan_num)
!	cdpeak = 0
    cdpeak = gaussian(ener_00, p_e, p_i, p_w, deltae)
	If ((ener_00.LT.p_e).AND.(ener_00.GT.(P_e-22.))) cdpeak = cdpeak + kfit_parms(29) * p_i
	If ((ener_00.LT.p_e).AND.(ener_00.LE.(P_e-22.))) cdpeak = cdpeak + kfit_parms(29) * p_i*0.55
    functn2 = temp_val + scatter + y + cdpeak + UPU_K_alpha + Pb_K_alpha+ Bi_K_alpha
!
!       Use fitted background or adjusted Ottmar background correction
        If (kfit_consts(12).NE.1) goto 73
          ENERGY=e_start+deltae*FLOAT(chan_num)
          step_size=KFIT_PARMS(12)

 	   FUNCTN2=FUNCTN2+BKG_STRIP2(CHAN_NUM, step_size)

 73	continue
!
!     generator endpoint peak
!
	if(kfit_consts(13).NE.1) GOTO 90
        energy=e_start +float(chan_num)*deltae
	gen_width=kfit_parms(33)
	gen_e=kfit_parms(8)+kfit_consts(14)				! endpoint energy + offset from constants
	gen_int=kfit_parms(32)*kfit_parms(9)				! peak fraction * X-ray intensity
	gen_area=gen_int*gaussian(energy, gen_e, 1.0D0, gen_width, deltae)
	FUNCTN2=FUNCTN2+gen_area
!	Calculate geneator peak tail area
	gen_width=kfit_parms(35)
	gen_int=kfit_parms(34)*kfit_parms(9)				! tail fraction * X-ray intensity
	gen_area=gen_int*DEXP(-DABS((energy-gen_e)/gen_width))		! double sided exponential
	FUNCTN2=FUNCTN2+gen_area
90        continue
!
!	add tungsten X-ray peaks 
!
!	to provide future capability for lower Z elements the tungsten lines will need to be included
!
	ENERGY=e_start+deltae*FLOAT(chan_num)
	ENER1=e_start+deltae*FLOAT(chan_num)
	X_AREA=KFIT_parms(26)*HPGE_eff(ENER1)/HPGE_eff(59.2D0)*deltae
	X_AREA=X_AREA*W_VOIGT(energy,Gauss_width,lo_tail_area,lo_tail_decay,fast_tail_area,fast_tail_decay,deltae)
	FUNCTN2=FUNCTN2+X_AREA
!
	if(isnan(functn2)) test1=1
!	if(isnan(functn2)) print *, chan_num, temp_val, scatter, y, cd_peak
!	print *, temp_val, scatter, y, cd_peak
	RETURN
	End
!
!	*********************************
!
	Function g_smooth(energy, deltae, Gauss_width, x,n)
!	Dim e0(2), i0(2), lor_w(2),
	DOUBLE PRECISION x(101), tails(4), ENERGY, E0, XDUM
	DOUBLE PRECISION gaussian, Gauss_width, Gauss_width2, DELTAE
!
	npks = 1
	XDUM=1.0
!
!	Gauss_width = 0.219
!	step_bkg_inc = 0.00043
!
	Gauss_width2 = Gauss_width*(energy/115.6)**0.5
	y = 0
!
	DO 10  i = 1 , n
	e0 = energy + float(int(n/2) - i) * deltae
	y = y + x(i) * gaussian(energy, e0, XDUM, Gauss_width2, deltae)
10	CONTINUE
!
	g_smooth = y
!
	RETURN
	End
!
!	*********************************
!
	Function bent_kramers4(num, energy, delta_angle, KFIT_PARMS)
	DOUBLE PRECISION  i_0
!
	DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
	DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
	DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100)
	INTEGER i_conc(7), i_attn(7), chan_num
        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
        COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
	DOUBLE PRECISION FUNCTN_KFIT, ENERGY_UNSCAT, ENERGY, y, bent_kramers4
!
	Pi = 3.14159265358979
	m_e = 510.99906
	bent_kramers4 = 0.
	num_div = 1
	d_a = delta_angle / num_div
!
	d_e = 0.05
	e_0 = kfit_parms(8)
	i_0 = kfit_parms(27)
	angle = kfit_parms(28)

!
	ALPHA = energy / m_e
	ALPHA2 = (energy + 0.05) / m_e
	rad_th = angle * Pi / 180
	cos_th = Cos(angle * Pi / 180)
!
	fF_1 = 1 / (1 - ALPHA * (1 - cos_th))
	fF_2 = 1 / (1 - ALPHA2 * (1 - cos_th))
!	d_scat_de = (fF_2 - fF_1) / 0.05
	d_scat_de = 1
	energy_unscat = fF_1 * energy
!
	y = FUNCTN_KFIT(energy_unscat, KFIT_PARMS)
!
	if (ENERGY_UNSCAT.LT.E_0) bent_kramers4 = i_0*y/d_scat_de
	RETURN
	End
!
!	*********************************
!
	Function functn_kfit(ener, KFIT_PARMS)
!
	DOUBLE PRECISION  temp_array(51), tails(4), e_escape(2)
	DOUBLE PRECISION ENER, ebel, ebel1, functn_kfit
	INTEGER i_cd,iconc(6)
!
	DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
	DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
	DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100), xray_ref_energy, energy
	DOUBLE PRECISION x, x_intensity, attn, hpge_eff
	INTEGER i_conc(7), i_attn(7), chan_num, vial_type
        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
        COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
!
	DIMENSION e_array(51), concent(7)
!
	DO 10 i = 1 , 7
	concent(i) = kfit_parms(i)
10	CONTINUE
	vial_type=kfit_consts(16)
!
	e_0 = kfit_parms(8)
	x_intensity = kfit_parms(9)
	Par_x = kfit_parms(10)
	par_y = kfit_consts(6)
	ang_1 = kfit_consts(7)
	ang_2 = kfit_consts(8)
	w_1 = kfit_consts(9)
	w_2 = kfit_consts(10)
!
	deltae = kfit_parms(14)
	gaus_width = kfit_parms(15)
!
	widthU = kfit_parms(16)
	widthpu = kfit_parms(17)
	ss_thick = kfit_parms(18)
	vial_d = kfit_consts(1)
       	vial_t = kfit_consts(2)
	be_thick = kfit_consts(3)
	cd_thick = kfit_consts(4)
	xray_ref_energy = kfit_consts(5)
!
	i_ss = i_attn(1)
	I_w = i_attn(2)
	icd = i_attn(3)
	i_poly = i_attn(4)
	i_be = i_attn(5)
!
	cd_rho = 8.65
	be_rho = 1.85
	rho_vial = 0.95
	if(vial_type.EQ.1)rho_vial=2.5
	ss_rho = 7.9
!
	energy = ener
	ebel1=ebel(energy,e_0,Par_x,par_y,ang_1,ang_2,w_1,w_2)
	x=x_intensity*ebel1
!
	attn=Exp(-ss_thick*ss_rho*attn_fact_ss_7(energy,i_ss))
	ATTN_TEMP=attn_fact_7_co(energy,deltae,KFIT_PARMS)
	attn=attn*Exp(-vial_d*ATTN_TEMP)
	attn=attn*Exp(-2.0*vial_t*rho_vial*attn_fact_vial(energy, vial_type))
	attn=attn*Exp(-be_thick*be_rho*attn_fact_be_co(energy))
	attn=attn*Exp(-cd_thick*cd_rho*attn_fact_cd_7(energy,icd))
!
	rel_eff=HPGE_eff(energy)/HPGE_eff(xray_ref_energy)
	functn_kfit=x*attn*rel_eff
!	print *, x_intensity, x, attn, ebel1
!	print *, energy, xray_ref_energy, functn_kfit
!	print *, '************************'
!
	RETURN
	END
!
!	*********************************
!
Subroutine read_amu_values()
!
!	Reads in attenuation coefficient parameters from file MEKED_attn_coeffs.txt
!
	DOUBLE PRECISION amu_lib(10,8), cmu_lib(2,10), sol_lib(2,10), edges(10), misc_lib(6,10), misc_edges(6)
	DOUBLE PRECISION sol_edges(2), cont_edges(2), lor_widths(10), act_norm_factor(10), energy
	CHARACTER FILNAME*60, FILNAM1*60 
	COMMON/attn_parms/amu_lib, cmu_lib, sol_lib, edges, sol_edges, cont_edges, lor_widths, act_norm_factor, misc_lib, misc_edges
	integer i, j, l1,l2
	REAL a1, b1, c1

    FILNAME='MEKED_attn_coeffs.txt'
!
    open (unit=10,file=FILNAME)
!
    read (10,*)FILNAM1

!	Get actinide attenuation coefficient parameters, amu_lib
	do 1111 l1 = 1, 10
	do 1111 l2 = 1, 8
	read (10,*) a1, b1, c1
	amu_lib(l1, l2)=c1
1111	continue

!	Get actinide K-Edge energies, edges  
	do 1112 l1 = 1, 10
	read (10,*) a1, edges(l1)
1112	continue

!	Get lorentzian widths to broaden the K-edge transistion, lor_widths
	do 1113 l1 = 1, 10
	read (10,*) a1, lor_widths(l1)
1113	continue

!	Get container attenuation coefficient parameters, cmu_lib
	do 1114 l1 = 1, 2
	do 1114 l2 = 1, 11
	read (10,*) a1, b1, c1
	if(l2.EQ.11) cont_edges(l1)=c1
	if(l2.LT.11) cmu_lib(l1,l2)=c1
1114	continue

!	Get solution attenuation coefficient parameters, sol_lib
	do 1115 l1 = 1, 2
	do 1115 l2 = 1, 11
	read (10,*) a1, b1, c1
	if(l2.EQ.11) sol_edges(l1)=c1
	if(l2.LT.11) sol_lib(l1,l2)=c1

1115	continue

!	Get attenuation factor correction factors, act_norm_factor
	do 1116 l1=1, 10
	read (10,*) a1, c1
	act_norm_factor(l1)=c1
1116	continue

!	Read attenuation factor correction factors errors, act_norm_factor
	do 1117 l1=1, 10
	read (10,*) a1, c1
1117	continue
!
!	Get misc attenuation coefficient parameters, misc_lib
	do 1118 l1 = 1, 3
	do 1118 l2 = 1, 11
	read (10,*) a1, b1, c1
	if(l2.EQ.11) misc_edges(l1)=c1
	if(l2.LT.11) misc_lib(l1,l2)=c1
	
1118	continue

    close (unit=10)

	RETURN
	END


!
!	*********************************
!
!
!	*********************************
!
	Function attn_fact_7_co(energy, deltae,KFIT_PARMS)
!
!	calculates 
!
	DOUBLE PRECISION  concent(7), LOG_E, energy
	DOUBLE PRECISION  MU_BELOW, MU_ABOVE,HNO3_attn_fact
!
	DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
	DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
	DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100)
	INTEGER i_conc(7), i_attn(7), chan_num

	DOUBLE PRECISION LOR_U, LOR_PU, LOR_AM, LOR_NP, LOR_CM
	DOUBLE PRECISION th_edge, U_EDGE, PU_EDGE, AM_EDGE, NP_EDGE
	DOUBLE PRECISION CM_EDGE, HNO3_edge, PI, ALG_E,Y
	DOUBLE PRECISION amu_lib(10,8), cmu_lib(2,10), sol_lib(2,10), edges(10), misc_lib(6,10), misc_edges(6)
	DOUBLE PRECISION sol_edges(2), cont_edges(2), lor_widths(10), act_norm_factor(10)
	DOUBLE PRECISION lorentz_w(5), act_attn_fact(6)

        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
        COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
	COMMON/attn_parms/amu_lib, cmu_lib, sol_lib, edges, sol_edges, cont_edges, lor_widths, act_norm_factor, misc_lib, misc_edges
!
	Pi = 3.1415926535898
!
	do 1 i = 1 , 7
	concent(i)=kfit_parms(i)
1	continue
	lorentz_w(1) = abs(KFIT_PARMS(16))			
	lorentz_w(2) = abs(KFIT_PARMS(17))
	lorentz_w(3)= lor_widths(3)
	lorentz_w(4) = lor_widths(4)
	lorentz_w(5) = lor_widths(5)

	HNO3_edge=sol_edges(1)
!
	attn_fact_7_co = 0.0
!
	attn1 = 0
!
!
	alg_e = Log(energy)
!
	i_dex=1
	amu_below = exp(amu_lib(i_dex,4)*alg_e**3+amu_lib(i_dex,3)*alg_e**2+amu_lib(i_dex,2)*alg_e + amu_lib(i_dex,1))
	amu_above = exp(amu_lib(i_dex,8)*alg_e**3+amu_lib(i_dex,7)*alg_e**2+amu_lib(i_dex,6)*alg_e + amu_lib(i_dex,5))
	y = 2.0*(energy-edges(i_dex))/lorentz_w(i_dex)
	atany = atan(y)
	act_attn_fact(i_dex)=(amu_below+(amu_above-amu_below)*(0.5+atany/Pi))
	act_attn_fact(i_dex)=act_attn_fact(i_dex)*concent(i_dex)
!	comment out the next line to broaden above the K-edge
	If (energy.GT.edges(i_dex)) act_attn_fact(i_dex) = amu_above*concent(i_dex)
!
!	Element 2
	i_dex=2
	amu_below = exp(amu_lib(i_dex,4)*alg_e**3+amu_lib(i_dex,3)*alg_e**2+amu_lib(i_dex,2)*alg_e + amu_lib(i_dex,1))
	amu_above = exp(amu_lib(i_dex,8)*alg_e**3+amu_lib(i_dex,7)*alg_e**2+amu_lib(i_dex,6)*alg_e + amu_lib(i_dex,5))
	y = 2.0*(energy-edges(i_dex))/lorentz_w(i_dex)
	atany = atan(y)
	act_attn_fact(i_dex)=(amu_below+(amu_above-amu_below)*(0.5+atany/Pi))
	act_attn_fact(i_dex)=act_attn_fact(i_dex)*concent(i_dex)
!	comment out the next line to broaden above the K-edge
	If (energy.GT.edges(i_dex)) act_attn_fact(i_dex) = amu_above*concent(i_dex)
!
!	Elements 3 to 5 are broadened above the K-edge (yes I could have used a loop)
!
!	Element 3
	i_dex=3
	amu_below = exp(amu_lib(i_dex,4)*alg_e**3+amu_lib(i_dex,3)*alg_e**2+amu_lib(i_dex,2)*alg_e + amu_lib(i_dex,1))
	amu_above = exp(amu_lib(i_dex,8)*alg_e**3+amu_lib(i_dex,7)*alg_e**2+amu_lib(i_dex,6)*alg_e + amu_lib(i_dex,5))
	y = 2.0*(energy-edges(i_dex))/lorentz_w(i_dex)
	atany = atan(y)
	act_attn_fact(i_dex)=(amu_below+(amu_above-amu_below)*(0.5+atany/Pi))
	act_attn_fact(i_dex)=act_attn_fact(i_dex)*concent(i_dex)
!
!	Element  4
	i_dex=4
	amu_below = exp(amu_lib(i_dex,4)*alg_e**3+amu_lib(i_dex,3)*alg_e**2+amu_lib(i_dex,2)*alg_e + amu_lib(i_dex,1))
	amu_above = exp(amu_lib(i_dex,8)*alg_e**3+amu_lib(i_dex,7)*alg_e**2+amu_lib(i_dex,6)*alg_e + amu_lib(i_dex,5))
	y = 2.0*(energy-edges(i_dex))/lorentz_w(i_dex)
	atany = atan(y)
	act_attn_fact(i_dex)=(amu_below+(amu_above-amu_below)*(0.5+atany/Pi))
	act_attn_fact(i_dex)=act_attn_fact(i_dex)*concent(i_dex)
!
!	Element  5
	i_dex=5
	amu_below = exp(amu_lib(i_dex,4)*alg_e**3+amu_lib(i_dex,3)*alg_e**2+amu_lib(i_dex,2)*alg_e + amu_lib(i_dex,1))
	amu_above = exp(amu_lib(i_dex,8)*alg_e**3+amu_lib(i_dex,7)*alg_e**2+amu_lib(i_dex,6)*alg_e + amu_lib(i_dex,5))
	y = 2.0*(energy-edges(i_dex))/lorentz_w(i_dex)
	atany = atan(y)
	act_attn_fact(i_dex)=(amu_below+(amu_above-amu_below)*(0.5+atany/Pi))
	act_attn_fact(i_dex)=act_attn_fact(i_dex)*concent(i_dex)
!
!	Nitric Acid
	i_dex=1
	amu_below = exp(sol_lib(i_dex,4)*alg_e**3+sol_lib(i_dex,3)*alg_e**2+sol_lib(i_dex,2)*alg_e + sol_lib(i_dex,1))
	amu_above = exp(sol_lib(i_dex,9)*alg_e**3+sol_lib(i_dex,8)*alg_e**2+sol_lib(i_dex,7)*alg_e + sol_lib(i_dex,6))
	HNO3_attn_fact=amu_below*concent(6)
	If (energy.GT.HNO3_edge) HNO3_attn_fact=amu_above*concent(6)
!
	attn_fact_7_co=0
	do 20 i_dex = 1, 5
	attn_fact_7_co =attn_fact_7_co + act_attn_fact(i_dex)*act_norm_factor(i_dex) 
20	continue
	attn_fact_7_co=attn_fact_7_co+HNO3_attn_fact
	RETURN
	END
!
!	*********************************
!
	Function ebel(energy, e_0, Par_x, par_y, ang_1, ang_2, w_1, w_2)
!
	DOUBLE PRECISION ebel, MU_Z, F_ABS, I_A, PAR_NU, DN_DE, energy, alnz, par_z
	DOUBLE PRECISION U_0
	ebel = 0.0
	If (energy.GT.e_0) GoTo 100
	d_e = 0.09
	cal_c = 1.0
	tim = 1.0
	I_A = 15.0
	omega = (0.4 / 10)**2

!	tungsten target
	par_z = 74
	par_a = 183.84
!
	Pi = 3.14159265358979
	phi = (90.0 - ang_1) / 180.0 * Pi
	eps = ang_2 / 180.0 * Pi
!
	par_m = 0.1382 - 0.9211 / par_z**0.5
	par_j = 0.0135 * par_z
	U_0 = e_0 / energy
!
	alnz=dlog(par_z)
	par_nu=(e_0**par_m)*(0.1904-0.2236*alnz+0.1292*alnz**2-0.0149*alnz**3)
!
	rho_zm = par_a / par_z*(0.00000787 * par_j**0.5*e_0**1.5 + 0.000000735 * e_0**2)
!
	rho_zbar=rho_zm*(0.49269-1.0987*par_nu+0.78557*par_nu**2)*dlog(U_0)/(0.70256-1.09865*par_nu+1.0046*par_nu**2+dlog(U_0))
!
	mu_z = attn_fact_W(energy) * (w_1 + w_2 * energy)
	f_abs=(1-Exp(-mu_z*2*rho_zbar*Sin(phi)/Sin(eps)))/(mu_z*2*rho_zbar*Sin(phi)/Sin(eps))
!
	dn_de=cal_c*omega*I_A*tim*par_z*(e_0/energy-1)**Par_x*f_abs*d_e
!
	ebel = dn_de * (e_0 / energy)**par_y
100	continue
!
	RETURN
	End
!
!	*********************************
!
	Function attn_fact_W_PE2(energy, par_1, par_2)
!	attn coeff Tungsten Photoelectric
!	no longer used in KED_FIT
!
	DOUBLE PRECISION energy
!
!	change over point
	edge = 69.53
!
	aL1 = 232603
	aL2 = -2.734
	aH1 = 938475
	aH2 = -2.679
!
	attn_fact_W_PE2 = 1
	attn_fact_W_PE2 = aL1*energy**aL2
	If (energy.GT.edge) attn_fact_W_PE2 = aH1*energy**aH2
!
	attn_fact_W_PE2 = (par_1 + par_2 * energy) * attn_fact_W_PE2
!
	RETURN
	END
!
!	*********************************
!
	Function attn_fact_W(energy)
!	attn coeff Tungsten 
!	used in function ebel 
!
	DOUBLE PRECISION energy
	edge = 69.53
!
	aL1 = 160974
	aL2 = -2.608
	aH1 = 493198
	aH2 = -2.523
!
	attn_fact_W = 1
	attn_fact_W = aL1 * energy**aL2
	If (energy.GT.edge) attn_fact_W = aH1 * energy**aH2
!
	RETURN
	End
!
!	*********************************
!
	Function attn_fact_ss_7(energy, i_co)
	DOUBLE PRECISION energy
	DOUBLE PRECISION amu_lib(10,8), cmu_lib(2,10), sol_lib(2,10), edges(10), misc_lib(6,10), misc_edges(6)
	DOUBLE PRECISION sol_edges(2), cont_edges(2), lor_widths(10), act_norm_factor(10)
	COMMON/attn_parms/amu_lib, cmu_lib, sol_lib, edges, sol_edges, cont_edges, lor_widths, act_norm_factor, misc_lib, misc_edges
!
!	attn coeff for SS
!	parameters read in from attenuation factors data file
!
	edge = misc_edges(1)
!	attn parameters with coherent scattering

	Amu_ss_L1 = misc_lib(1,1)
	Amu_ss_L2 = misc_lib(1,2)
	Amu_ss_L3 = misc_lib(1,3)
	Amu_ss_L4 = misc_lib(1,4)
	Amu_ss_L5 = misc_lib(1,5)

	Amu_ss_H1 = misc_lib(1,6)
	Amu_ss_H2 = misc_lib(1,7)
	Amu_ss_H3 = misc_lib(1,8)
	Amu_ss_H4 = misc_lib(1,9)
	Amu_ss_H5 = misc_lib(1,10)
!
100	continue

	al_e= log(energy)

	Amu_below = exp(Amu_ss_L5*al_e**4+Amu_ss_L4*al_e**3+Amu_ss_L3*al_e**2+Amu_ss_L2*al_e + Amu_ss_L1)
	Amu_above = exp(Amu_ss_H5*al_e**4+Amu_ss_H4*al_e**3+Amu_ss_H3*al_e**2+Amu_ss_H2*al_e + Amu_ss_H1)
!
	attn_fact_ss_7 = Amu_below
	If (energy.GT.edge) attn_fact_ss_7 = Amu_above
!
	RETURN
	End
!
!	*********************************
!
	Function attn_fact_ss_7_old(energy, i_co)
	DOUBLE PRECISION energy
!
!	attn coeff for SS
!	use of the more accurate reprentation in attn_fact_ss_7 seems to have adversely impacted the matrix density determination
!	older function has been left in case we want to explore the impact on performannce
!
	edge = 60
!	attn parameters with coherent scattering
	Amu_ss_L1 = -0.1964
	Amu_ss_L2= -2.2074
	Amu_ss_L3 = 4.6049
	Amu_ss_H1 = 1.8729
	Amu_ss_H2 = -9.3876
	Amu_ss_H3 = 10.8492
!
	If (i_co.EQ.1) GoTo 100
!
!	attn parameters without coherent scattering
	Amu_ss_L1 = -0.2392
	Amu_ss_L2 = -2.1287
	Amu_ss_L3 = 4.5662
	Amu_ss_H1 = 2.0803
	Amu_ss_H2 = -10.227
	Amu_ss_H3 = 11.6519
!
100	continue
	al_e = Log10(energy)
	Amu_below = 10**(Amu_ss_L1*al_e**2+Amu_ss_L2*al_e+Amu_ss_L3)
	Amu_above = 10**(Amu_ss_H1*al_e**2+Amu_ss_H2*al_e+Amu_ss_H3)
!
	attn_fact_ss_7_old = Amu_below
	If (energy.GT.edge) attn_fact_ss_7_old = Amu_above
!
	RETURN
	End
!
!	*********************************
!
	Function attn_fact_Ge_7(energy, i_co)
	DOUBLE PRECISION energy
!
!	attn coeff for Ge
! 	used to calculate detector backscatter in function fack_scat
!
	edge = 11.1
!	attn parameters with coherent scattering
	amu_ge_L1 = -0.3273
	amu_ge_L2 = -2.2106
	amu_ge_L3 = 4.1075
	amu_ge_H1 = 0.4545
	amu_ge_H2 = -4.105
	amu_ge_H3 = 6.1549
!
	If (i_co.EQ.1) GoTo 100
!
!	attn parameters without coherent scattering
	amu_ge_L1 = -0.3674
	amu_ge_L2 = -2.1865
	amu_ge_L3 = 4.1042
	amu_ge_H1 = 0.4528
	amu_ge_H2 = -4.1419
	amu_ge_H3 = 6.1966
!
100	continue
!
	alog_e = Log10(energy)
	amu_below = 10**(amu_ge_L1*alog_e**2+amu_ge_L2*alog_e+amu_ge_L3)
	amu_above = 10**(amu_ge_H1*alog_e**2+amu_ge_H2*alog_e+amu_ge_H3)
!
	attn_fact_Ge_7 = amu_below
!
	If (energy.GT.edge) attn_fact_Ge_7 = amu_above
	RETURN
	END
!
!	*********************************
!
	Function attn_fact_cd_7(energy, i_co)
	DOUBLE PRECISION energy
	DOUBLE PRECISION amu_lib(10,8), cmu_lib(2,10), sol_lib(2,10), edges(10), misc_lib(6,10), misc_edges(6)
	DOUBLE PRECISION sol_edges(2), cont_edges(2), lor_widths(10), act_norm_factor(10)
	COMMON/attn_parms/amu_lib, cmu_lib, sol_lib, edges, sol_edges, cont_edges, lor_widths, act_norm_factor, misc_lib, misc_edges
!
!	attn coeff for cadmium
!	
!	attn parameters with coherent scattering
!
	al_e2=log(energy)
	edge=misc_edges(3)

	amu_cd_L1 = misc_lib(3,1)
	amu_cd_L2 = misc_lib(3,2)
	amu_cd_L3 = misc_lib(3,3)
	amu_cd_L4 = misc_lib(3,4)
	amu_cd_L5 = misc_lib(3,5)

	amu_cd_H1 = misc_lib(3,6)
	amu_cd_H2 = misc_lib(3,7)
	amu_cd_H3 = misc_lib(3,8)
	amu_cd_H4 = misc_lib(3,9)
	amu_cd_H5 = misc_lib(3,10)

	Amu_below = exp(amu_cd_L5*al_e2**4+amu_cd_L4*al_e2**3+amu_cd_L3*al_e2**2+amu_cd_L2*al_e2 + amu_cd_L1)
	Amu_above = exp(amu_cd_H5*al_e2**4+amu_cd_H4*al_e2**3+amu_cd_H3*al_e2**2+amu_cd_H2*al_e2 + amu_cd_H1)
!
	attn_fact_cd_7 = amu_below
	If (energy.GT.edge) attn_fact_cd_7 = amu_above

	If (i_co.EQ.1) GoTo 100
!
!	attn parameters without coherent scattering
!
	edge = 102
	amu_cd_L1 = 0.03
	amu_cd_L2 = -2.8373
	amu_cd_L3 = 5.7089
	amu_cd_H1 = 1.2242
	amu_cd_H2 = -7.6469
	amu_cd_H3 = 10.555
!
	alog_e = Log10(energy)
	amu_below = 10**(amu_cd_L1*alog_e**2+amu_cd_L2*alog_e+amu_cd_L3)
	amu_above = 10**(amu_cd_H1*alog_e**2+amu_cd_H2*alog_e+amu_cd_H3)
!
	attn_fact_cd_7 = amu_below
	If (energy.GT.edge) attn_fact_cd_7 = amu_above

100	continue

	RETURN
	END
!
!	*********************************
!
	Function attn_fact_vial(energy, vial_type)
	DOUBLE PRECISION energy, alg_e
!	DOUBLE PRECISION amu_POLY_L0, amu_POLY_L1, amu_POLY_L2, amu_POLY_L3
!	DOUBLE PRECISION amu_SiO2_L0, amu_SiO2_L1, amu_SiO2_L2, amu_SiO2_L3
	DOUBLE PRECISION amu_L0, amu_L1, amu_L2, amu_L3
	INTEGER vial_type
!
	DOUBLE PRECISION amu_lib(10,8), cmu_lib(2,10), sol_lib(2,10), edges(10), misc_lib(6,10), misc_edges(6)
	DOUBLE PRECISION sol_edges(2), cont_edges(2), lor_widths(10), act_norm_factor(10)
	COMMON/attn_parms/amu_lib, cmu_lib, sol_lib, edges, sol_edges, cont_edges, lor_widths, act_norm_factor, misc_lib, misc_edges
!
!	attn coeff for polyethylene
!
	amu_L0=cmu_lib(1,1)
	amu_L1=cmu_lib(1,2)
	amu_L2=cmu_lib(1,3)
	amu_L3=cmu_lib(1,4)

	if(vial_type.NE.1) GOTO 10

	amu_L0=cmu_lib(2,1)
	amu_L1=cmu_lib(2,2)
	amu_L2=cmu_lib(2,3)
	amu_L3=cmu_lib(2,4)
!
10	continue
	alg_e=LOG(energy)
	attn_fact_vial = exp(amu_L3*alg_e**3+amu_L2*alg_e**2+amu_L1*alg_e+amu_L0)

	RETURN
	END
!
!	*********************************
!
	Function attn_fact_be_co(energy)
	DOUBLE PRECISION energy
	DOUBLE PRECISION amu_lib(10,8), cmu_lib(2,10), sol_lib(2,10), edges(10), misc_lib(6,10), misc_edges(6)
	DOUBLE PRECISION sol_edges(2), cont_edges(2), lor_widths(10), act_norm_factor(10)
	COMMON/attn_parms/amu_lib, cmu_lib, sol_lib, edges, sol_edges, cont_edges, lor_widths, act_norm_factor, misc_lib, misc_edges
!
!	attn coeff for Beryllium
!
	aH1 = misc_lib(2, 1)
	aH2 = misc_lib(2, 2)
!
	attn_fact_be_co = aH1 * energy**aH2
!
	RETURN
	END
!
!	*********************************
!
	Function HPGE_eff(energy)
!
!	ISOCS EFFICIENT CALCULATION FOR GL0210 HPGE DETECTOR
!	Detector is located in a tungsten colimator
!	aperature is 0.8 mm ID
!	material is uniformly distributed in a low density matrix (i.e. air)
!	No attenuators are used,the vial material is set to air.
!	vial ID = 1.413 cm
!
!	the intent is use the results to determine relative - not absolute - efficiency with energy
!
	DOUBLE PRECISION energy, HPGE_eff, eff
	REAL eff_p(7)

	eff_p(1) = -15.66
	eff_p(2) = -0.2712
	eff_p(3) = 1.082
	eff_p(4) = -0.09303
	eff_p(5) = -0.1436
	eff_p(6) = 0.02888
	eff_p(7) = 1020.0

!
	x = Log(eff_p(7) / energy)
	eff = 0.0
!
	DO 10 i = 1 , 6
	eff = eff + eff_p(i) * x**(i - 1)
10	CONTINUE
!
	eff = DExp(eff)
	HPGE_eff = eff
	RETURN
	END
!
!	*********************************
!
	Function gaussian(energy, e0, P0, SIGMA, deltae)
!
	DOUBLE PRECISION ENERGY, E0, P0, SIGMA, DELTAE
	DOUBLE PRECISION Pi, y1, y, gaussian
	Pi = 3.14159265358979
	y1 = ((energy - e0)**2)/2.0/ SIGMA**2
	y =deltae*P0*Exp(-y1)/(2.0 * Pi)**0.5/SIGMA
	gaussian = y
	RETURN
	END
!
!	*********************************
!
	Function escape_area(energy, np)
!
!	np is escape peak 1 or 2 (9.88 or 10.99 keV)
	DOUBLE PRECISION ENERGY
	DIMENSION esc_e(2), esc_mu(2)
	esc_e(1) = 9.88
	esc_e(2) = 10.99
	esc_mu(1) = 37.42
	esc_mu(2) = 37.42
!
	edge = 11.1
!
	frac = 0.467
	If (np.EQ.2) frac = 0.0603
!
	A1 = 77700
	A2 = -2.523
!
	attn_fact_ge = 2.81
!
	If (energy.GT.edge)attn_fact_ge = A1*energy**A2
!
	escape_area = frac*attn_fact_ge/esc_mu(np)/2.5
	RETURN
	END
!
!	*********************************
!
	Function spec_fold2(nlow, nhigh, npt, gate, live_time, XRAW)
!
!	computes the random coincidence summing function
!
!	RAW_DATA contains the raw spectrum in counts
!	nlow, nhigh are the start and stop channels over which the convolution is performed
!	npt is the channel number of the resultant spectrum
!	gate is the coincidence gate time in seconds (~0.5 E-6)
!	Note: This is a first order correction. However, for summed energies of less than
!	150 kev the second order terms should be negligible.
!
!
	DOUBLE PRECISION XRAW(4100), GATE, spec_fold2
	INTEGER nlow, nhigh, npt, live_time
!
	alive_time=float(live_time)
	temp = 0.0
	do 10, i = nlow , nhigh
	If (i.GT.npt) GoTo 20
	If ((npt - i).LT. 1) GoTo 10
	temp=temp+XRAW(i)*XRAW(npt-i)*gate/alive_time**2
10 	continue
20 	continue
!
	spec_fold2 = temp
	RETURN
	END
!
!	*********************************
!
    Function spec_tail(energy,lo_tail_area,lo_tail_decay,delta_e)
	DOUBLE PRECISION LO_TAIL_AREA, lo_tail_decay
	DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
	DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
	DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100)
	INTEGER i_conc(7), i_attn(7), chan_num
        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
        COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
	DOUBLE PRECISION ENERGY
!
!	TAILING IS A REACH AHEAD APPROACH, NREACH CHANNELS
!
	fast_tail_area = 0.0
	fast_tail_decay = 1.0
	NREACH=50
	temp = 0.0
	do 10 i = 1 , NREACH
	e_i = i * delta_e
	temp=temp+RAND_TEMP(i)*lo_tail_area*Exp(-e_i*lo_tail_decay)
!
10	continue
	spec_tail = temp
	RETURN
	END
!
!	*********************************
!
    Function back_scat(energy_scat1,i_10, KFIT_PARMS)
!
    DOUBLE PRECISION ENERGY_SCAT1, back_scat, e_1
    INTEGER I, i_10, i_1, I_START
    DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
    DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
    DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100)
    INTEGER i_conc(7), i_attn(7), chan_num
    COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
    COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
!
    DOUBLE PRECISION F_S, M_E, R_C, E_0
    DOUBLE PRECISION X_INTENSITY,KRAM_POWER,SCAT_I,SCAT_E,SPEC_0
    DOUBLE PRECISION DELTA_E, TEMP_SUM,angle_min, E1
    DOUBLE PRECISION COS_SCAT_DAT, ang_scat, P_E
    DOUBLE PRECISION ATTN_TEMP, REL_COR, TEMP1
!
    energy_scat = energy_scat1 + 1.8
    i_0 = i_10 + 20
    f_s = 1.0 / 137.04
    M_E = 511.0
    r_c = 0.38616
    e_0 = kfit_parms(8)
    x_intensity = kfit_parms(9)
    Kram_power = kfit_parms(10)
    scat_i = Abs(kfit_parms(11))
    spec_0 = kfit_parms(13)
    delta_e = kfit_parms(14)
    temp_sum = 0.0
    angle_min = (KFIT_CONSTS(11)*3.1415926/180)
!
!	i_1 = int((energy - spec_0) / delta_e)
!	I_start = Max(i_0, i_1)
!
    DO 50 I = I_0 , I_0+1000
    e_1 = spec_0 + delta_e * (i)
!
    If (e_1.GT.e_0) GoTo 100
    cos_ang_scat = 0.0
!
    TEMP2=Abs(1.0-(e_1/energy_scat-1.0)*m_e/e_1)
    If (TEMP2.GT.1) GoTo 50
    cos_ang_scat=(1.0-(e_1/energy_scat-1.0)*m_e/e_1)
    ang_scat = Acos(cos_ang_scat)
!
    If (ang_scat.LT.angle_min) GoTo 50
    P_E=1./(1+e_1/m_e*(1.-cos_ang_scat))
    ATTN_TEMP=Exp(-attn_fact_Ge_7(e_1, 1)*5.323* 1.0)
    rel_cor=ATTN_TEMP/(1-ATTN_TEMP)
    TEMP1=xdat(i)*(f_s*r_c*P_E)**2 * (P_E+1/P_E-1+cos_ang_scat**2)/2*rel_cor
!
    If (cos_ang_scat.NE.0) temp_sum = temp_sum + TEMP1
50  continue
100 continue
!
    back_scat = temp_sum
    RETURN
    END
!
!	*********************************
!
    SUBROUTINE BKG_STRIP(KFIT_PARMS, BKG, cd_init)
!
    INTEGER LO_CHAN1, LO_CHAN2,HI_CHAN1, HI_CHAN2
    INTEGER LO_AVG_CHAN, HI_AVG_CHAN
    DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
    DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
    DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100), BKG(4100)
    DOUBLE PRECISION BKG_ADJUST, cd_init
    INTEGER i_conc(7), i_attn(7), chan_num, mxp
    COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
    COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
	COMMON/filesize/mxp
    REAL E_LO1,E_LO2,E_HI1,EHI2
!
    PRINT *, '-------BKG SUB----------------'
	E_LO1=kfit_consts(30)
	E_LO2=kfit_consts(31)
	E_HI1=kfit_consts(32)
	E_HI2=kfit_consts(33)
	BKG_ADJUST=kfit_consts(34)*cd_init
    E0=KFIT_PARMS(13)
    E1=KFIT_PARMS(14)
    LO_CHAN1=INT(E_LO1-E0)/E1
    LO_CHAN2=INT(E_LO2-E0)/E1
    HI_CHAN1=INT(E_HI1-E0)/E1
    HI_CHAN2=INT(E_HI2-E0)/E1
    DIFF1=LO_CHAN2-LO_CHAN1+1
    DIFF2=HI_CHAN2-HI_CHAN1+1
    LO_AVG_CHAN=INT((LO_CHAN1+LO_CHAN2)/2+.5)
    HI_AVG_CHAN=INT((HI_CHAN1+HI_CHAN2)/2+.5)
!
    SUM_LO=0.
    SUM_HI=0.
    SUM_CENTER=0.
!
    DO 10 I= LO_CHAN1, LO_CHAN2
    SUM_LO=SUM_LO+XDAT(I)
10  CONTINUE
!
    DO 20 I= HI_CHAN1, HI_CHAN2
    SUM_HI=SUM_HI+XDAT(I)
20  CONTINUE
!
    DO 30 I= LO_CHAN2+1, HI_CHAN1-1
    SUM_CENTER=SUM_CENTER+XDAT(I)
30  CONTINUE
!
    AVG_LO=(SUM_LO/DIFF1-BKG_ADJUST)*kfit_consts(35)
    AVG_HI=SUM_HI/DIFF2
!
    DO 40 I=1, LO_AVG_CHAN
    XDAT(I)=XDAT(I)-AVG_LO
40  CONTINUE
    DO 50 I=HI_AVG_CHAN, mxp-1
    XDAT(I)=XDAT(I)-AVG_HI
50  CONTINUE
!
    RUN_SUM=0.0
    PRINT *, cd_init, BKG_ADJUST, SUM_LO
    PRINT *, LO_CHAN1, LO_CHAN2,HI_CHAN1, HI_CHAN2
    PRINT *, LO_AVG_CHAN, HI_AVG_CHAN
    PRINT *, DIFF1, DIFF2
    PRINT *, AVG_LO, AVG_HI, SUM_CENTER
    DO 60 I = LO_AVG_CHAN+1, HI_AVG_CHAN-1
    RUN_SUM = RUN_SUM + XDAT(I)
    TEMPBKG=AVG_HI+(1.0-RUN_SUM/SUM_center)*(AVG_LO-AVG_HI)
    BKG(I)=TEMPBKG
60  CONTINUE
    RETURN
    END
!
!	*********************************
!
	FUNCTION BKG_STRIP2(I_chn, step_size)
!
	DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
	DOUBLE PRECISION STEP_SIZE,SUM0
	DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
	DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100)
	INTEGER i_conc(7), i_attn(7), chan_num, I_CHN, mxp
        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
        COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
	COMMON/filesize/mxp

!
!	PRINT *, '-------BKG SUB2----------------'
	SUM0=0.0
!
	SUM0=0.0
	DO 20 J= I_CHN, mxp-1
	SUM0=SUM0+STEP_SIZE*XDAT(J)
20	CONTINUE
	BKG_STRIP2=SUM0
!	print *, I_chn, step_size, sum0
	RETURN
	END
!
!	*********************************
!
	SUBROUTINE BKG_STRIP3(step_size, BKG)
!
	DOUBLE PRECISION kfit_parms(36), kfit_consts(36)
	DOUBLE PRECISION STEP_SIZE,SUM0
	DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
	DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100), BKG(4100)
	INTEGER i_conc(7), i_attn(7), chan_num, I_CHN, mxp
        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
        COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
	COMMON/filesize/mxp

!
!	PRINT *, '-------BKG SUB2----------------'
	SUM0=0.0
!
	DO 10 I= 1, mxp
	bkg(i)=0.0
10	CONTINUE
!
	do 25 I=1, mxp-1
	SUM0=0.0
	DO 20 J= I, mxp
	SUM0=SUM0+STEP_SIZE*XDAT(J)
20	CONTINUE
	BKG(I)=SUM0
!	PRINT *, I, STEP_SIZE, SUM0
25	CONTINUE
	RETURN
	END
!
!	*********************************
!
  SUBROUTINE curfit(y,ERROR_m,npts, nterms, MODE,A1, kfit_parms, &
  deltaa, SIGMAA, FLAMDA, yfit, CHISQR, NMIN, &
  free_params, idex_a, COVAR_ARRAY)
!
	DOUBLE PRECISION kfit_parms(36), kfit_consts(36), functn2
	DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
	DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100)
	DOUBLE PRECISION YFIT(4100), Y(4100), ERROR_M(4100)
	DOUBLE PRECISION CHISQ1, FCHISQ1
	INTEGER i_conc(7), i_attn(7), chan_num, IDEX_A(36)
!
	INTEGER NPTS, NTERMS, MODE, NMIN, FREE_PARAMS(36)
	dOUBLE PRECISION SIGMAA(36), FLAMDA,DELTAA(36)
	DOUBLE PRECISION CHISQR
        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
        COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
!	COMMON/CALCDAT/YFIT
!
	DOUBLE PRECISION array_m(36, 36), BDUM(36), kfit_temp(36)
	DOUBLE PRECISION WEIG(4100), ALPHA(36, 36), BETA(36), deriv(36)
	DOUBLE PRECISION A1(36), b(36)
	INTEGER JX, IX
	DOUBLE PRECISION ALPHA_TEMP, TEMPOUT, COVAR_ARRAY(36,36)
!
!	*********************************************************************************
!	*    Based in large part on Bevington's CURFIT routine in                       *
!	*   "Data Reduction and Error Analysis fro the Physical Sciences" (1st Edition) *
!	*********************************************************************************
!	
	nfree = npts - nterms
	NMAX = NMIN + npts - 1
	matloops = 0
11	nfree = npts - nterms
!	----------------------  EVALUATE WEIGHTS ----------------
20	continue
        DO 29  i = NMIN , NMAX
	WEIG(i)=0.
        If (ERROR_m(i).EQ.0.) ERROR_m(i) = Abs(y(i)) ** 0.5
!    	If (ERROR_m(i).NE.0.) ERROR_m(i) = (Abs(y(i)) ** 0.5)
29      WEIG(i) = 1 / ERROR_m(i) ** 2
!
!	--------------- EVALUATE ALPHA AND BETA MATRICES -----------
!
31	DO 34 J = 1, nterms
        BETA(j) = 0
        DO 34 k = 1 , j
34      ALPHA(j, k) = 0
!
	Call a_to_kfit(a1, kfit_parms, idex_a, nterms)
41	DO 50 i = NMIN , NMAX
	Call fderiv2(i,a1,nterms,deltaa,deriv,free_params,idex_a,y,kfit_parms)
!
	DO 46 j = 1 , nterms
	BETA(j)=BETA(j)+WEIG(i)*(y(i)-functn2(i,kfit_parms))*deriv(j)
	DO 46  k = 1 , j
46	ALPHA(j,k)=ALPHA(j,k)+WEIG(i)*deriv(j)*deriv(k)
!
50	CONTINUE
51      DO 53 j = 1 , nterms
        DO 53 k = 1 , j
53      ALPHA(k, j) = ALPHA(j, k)
!
!	----------   EVALUATE CHI SQUARE AT STARTING POINT    ---------
!
61	DO 62 i = NMIN , NMAX
62      yfit(i) = functn2(i,kfit_parms)
63	CHISQ1 = FCHISQ1(y, WEIG, npts, nfree, NMIN, yfit)
!
!	 --------  INVERT MODIFIED CURVATURE MATRIX TO FIND NEW PARAMETERS
!
71	DO 74 j = 1 , nterms
	DO 73 k = 1 , nterms
        F1 = (ALPHA(j, j) * ALPHA(k, k)) ** 0.5
73      If (F1 .NE. 0) array_m(j, k) = (ALPHA(j, k) / F1)
74	Array_m(j, j) = 1.0 + (FLAMDA)
!
	matloops = matloops + 1
80	Call MATInv(nterms, det, array_m)
81	DO 84 j = 1 , nterms
	b(j) = a1(j)
	DO 84 k = 1 , nterms
    	F1 = (ALPHA(j, j) * ALPHA(k, k)) ** 0.5
84  	If (F1.NE.0.0) b(j)=b(j)+BETA(k)*array_m(j,k)/F1
!
!	 ------------- IF CHI SQUARE INCREASED, INCREASE FLAMDA AND TRY AGAIN
!
        DO 87 i_101 = 1 , 36
        kfit_temp(i_101) = kfit_parms(i_101)
87	CONTINUE
!
        Call a_to_kfit(b, kfit_temp, idex_a, nterms)
!
91	DO 92 i = NMIN , NMAX
92	yfit(i) = functn2(i,kfit_temp)
93	CHISQR = FCHISQ1(y, WEIG, npts, nfree, NMIN, yfit)
!
	if (isnan(CHISQR)) GOTO 101		
!		by chance a trial parameter has caused the valuation of functn2 to result in NaN (not a number)
!		in this case the loop is bypassed and the next iteration is begun
!
	If (CHISQ1.GE.CHISQR) GOTO 101
	if (matloops.GE.50) GOTO 101
	PRINT *, 'loops ', matloops , CHISQR
95	FLAMDA = 10.0 * FLAMDA
	GoTo 71
!
!	-----------  EVALUATE NEW PARAMETERS
!
101	CONTINUE
	if (isnan(CHISQR)) PRINT *, " **** Error in minimization - loop bypassed **** "
	DO 103 j = 1 , nterms
	if (isnan(CHISQR)) goto 102
	a1(J) = b(J)
102	alpha_temp = ALPHA(j, j)
	If (alpha_temp .EQ. 0.0) SIGMAA(j) = 0.0
    	If (alpha_temp .NE.0.0)SIGMAA(j)=((1.+FLAMDA)*array_m(j,j)/alpha_temp)**0.5
103	CONTINUE
        FLAMDA = FLAMDA / 10.0
	Call a_to_kfit(a1, kfit_parms, idex_a, nterms)
110	continue
!
	PRINT *, I, KFIT_PARMS(I), SIGMAA(I)
!
!	CALCULATE COVARIANCE MATRIX
!
	DO 150 JX = 1 , NTERMS
	DO 150 IX = 1 , NTERMS
	alpha_temp = ALPHA(IX, JX)
	If (alpha_temp .EQ. 0) tempout = 0.0
	If (alpha_temp.NE.0)tempout=((1.0+FLAMDA)*array_m(IX, JX)/alpha_temp)
	COVAR_ARRAY(IX,JX) = tempout
150	CONTINUE
  	PRINT *, 'IN CURFIT CHISQR =', CHISQR, CHISQ1
	RETURN
	END
!
!	*********************************
!
	SUBROUTINE fderiv2(chan_num, a1, nterms, deltaa, deriv,free_params, idex_a, ydat, kfit_parms)
!
!
	DOUBLE PRECISION A1(36), b(36), kfit_temp(36)
!
	DOUBLE PRECISION kfit_parms(36), kfit_consts(36), functn2
	DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
	DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100), YDAT(4100)
	DOUBLE PRECISION YFIT(4100), DERIV(36), DELTAA(36)
	INTEGER i_conc(7), i_attn(7), chan_num, IDEX_A(36)
	INTEGER NPTS, num_NFREE, Nmin, MODE, NTERMS
	INTEGER FREE_PARAMS(36)
        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
        COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
!
	orig_val = functn2(chan_num,kfit_parms)
	DO 10 i = 1 , 36
	kfit_temp(i) = kfit_parms(i)
10	CONTINUE
!
	DO 20 i = 1 , nterms
	b(i) = a1(i)
20	CONTINUE
!
	DO 30 i = 1 , nterms
	b(i) = a1(i) - deltaa(i)
	Call a_to_kfit(b, kfit_temp, idex_a, nterms)
!
	low_val = functn2(chan_num, kfit_temp)
!
	b(i) = a1(i) + deltaa(i)
	Call a_to_kfit(b, kfit_temp, idex_a, nterms)
!
	high_val = functn2(chan_num, kfit_temp)
!
	DERIV(I)=0.
	If (deltaa(i).NE.0)deriv(i) = (high_val - low_val) / 2 / deltaa(i)
	b(i) = a1(i)
30	CONTINUE
!
	RETURN
	END
!
!	*********************************
!
!       M A T I N V
!
!       PURPOSE
!     INVERT A SYMETRIC MATRIX AND CALCULATE ITS DETERMINANT
!
!       USAGE
!     CALL MATINV (NORDER,DET)
!
!       DESCRIPTION OF PARAMETERS
!     ARRAY  -  INPUT MATRIX WHICH IS REPLACED BY ITS INVERSE
!     NORDER -  DEGREE OF MATRIX (ORDER OF DETERMINANT)
!     DET    -  DETERMINANT OF INPUT MATRIX
!
!       SUBROUTINES USED
!     NONE
!
!       COMMENTS
!     DIMENSION STATEMENTS VALID FOR NORDER UP TO 10
!     THIS ROUTINE CONVERTS THE SINGLE PRECISION ARRAY TO A
!     DOUBLE PRECISION ARRAY FOR THE COMPUTATIONS. THE ARRAY RETURNED
!     IS SINGLE PRECISION
!
!       THIS IS A MODIFICATION OF THE ROUTINE IN BEVINGTON PP.302-303.
!           FEB.  1978,   J.R.BOYCE
!                   MODIFIED FOR OVERLAY PROGS JAN. 1984   R.D.M.
!
        SUBROUTINE MATINV(NORDR,DET, ARRAY)
        DOUBLE PRECISION ARRAY(36,36),AMAX,SAVE,A,B,BDUM(36)
        DIMENSION IK(36),JK(36)
!    COMMON/POLQ/BDUM,ARRAY
!       TYPE *,' $'
        DET=1.
        DO 100 K=1,NORDR
!-------FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX
        AMAX=0.D0
21      DO 30 I=K,NORDR
        DO 30 J=K,NORDR
        A=DABS(AMAX)
        B=DABS(ARRAY(I,J))
        IF(A-B) 24,24,30
24      AMAX=ARRAY(I,J)
        IK(K)=I
        JK(K)=J
30      CONTINUE
!-------INTERCHANGE ROWS AND COLLUMNS TO PUT AMAX IN ARRAY(K,K)
        IF(AMAX) 41,32,41
32      DET=0.
        RETURN
41      I=IK(K)
        IF (I-K) 21,51,43
43      DO 50 J=1,NORDR
        SAVE=ARRAY(K,J)
        ARRAY(K,J)=ARRAY(I,J)
50      ARRAY(I,J)=-SAVE
51      J=JK(K)
        IF(J-K) 21,61,53
53      DO 60 I=1,NORDR
        SAVE = ARRAY(I,K)
        ARRAY(I,K)=ARRAY(I,J)
60      ARRAY(I,J)=-SAVE
!-------ACCUMULATE ELEMENTS OF INVERSE MATRIX
61      DO 70 I=1,NORDR
        IF(I.EQ.K) GO TO 70
        ARRAY(I,K)=-ARRAY(I,K)/AMAX
70      CONTINUE
        DO 80 I=1,NORDR
        DO 80 J=1,NORDR
        IF(I.EQ.K) GO TO 80
        IF(J.EQ.K) GO TO 80
        ARRAY(I,J)=ARRAY(I,J) + ARRAY(I,K)*ARRAY(K,J)
80      CONTINUE
        DO 90 J=1,NORDR
        IF(J.EQ.K) GO TO 90
        ARRAY(K,J)= ARRAY(K,J)/AMAX
90      CONTINUE
        ARRAY(K,K)=1.D0/AMAX
100     DET=DET*SNGL(AMAX)
!     TYPE *,' $$'
!-------RESTORE ORDERING OF MATRIX
        DO 130 L=1,NORDR
        K=NORDR-L+1
        J=IK(K)
        IF(J.LE.K) GO TO 111
        DO 110 I=1,NORDR
        SAVE=ARRAY(I,K)
        ARRAY(I,K)=-ARRAY(I,J)
110     ARRAY(I,J)=SAVE
111     I=JK(K)
        IF(I.LE.K) GO TO 130
        DO 120 J=1,NORDR
        SAVE=ARRAY(K,J)
        ARRAY(K,J)=-ARRAY(I,J)
120     ARRAY(I,J)=SAVE
130     CONTINUE
        RETURN
        END
!
!	*********************************
!
	Function W_VOIGT(energy,Gauss_width,lo_tail_a,lo_tail_decay,fast_tail_a,fast_tail_decay,deltae)
!
!   Creates tungsten X-ray peaks
!	The relative peak intensities will not be very accurate because the source distribution/attenuation in not defined.
!	Estimation of the true coincidence summing peaks will be further complicated by the source distribution
!
	DOUBLE PRECISION ENERGY, GAUSS_WIDTH, E0(12), I0(12),lor_w(12)
	DOUBLE PRECISION lo_tail_a,lo_tail_decay,fast_tail_a, czosnyka
	DOUBLE PRECISION deltae,gauss_width_0, Y_TEMP, LOR_WIDTH
	DOUBLE PRECISION TAIL, E_TEMP,fast_tail_decay, E_1, char_x1, tail_x
!
 	DATA E0/59.318, 57.981, 57.425, 67.244, 69.067, 66.950, 69.273, 67.685, 69.484, 68.82, 80.0, 80.1/
! 	DATA E0/59.318, 57.981, 57.425, 67.235, 69.088, 75.543, 77.042, 78.846,  0, 0, 0, 0/
	DATA lor_w/0.0432,0.0374,0.043,0.0486,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05/
!   
!
! 	 DATA I0/0.3769, 0.2197, 0.0001, 0.0826, 0.0287, 0.0429, 0.0004, 0.0019, 0.004, 0.0265, 0.0634, 0.0305/
!	 DATA I0/0.47,0.274,0.000206,0.103,0.0358,0.0535,0.0006,0.00241,0.0051,0.0, 0.0, 0.0/
! 	 DATA I0/0.37699, 0.21978, 0.00017, 0.23239, 0.07365, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0/
	 DATA I0/0.47, 0.25818, 0.00018, 0.1785625, 0.0665375, 0.0916875, 0.0011125, 0.0042375, 0.009625, 0.0773 , 0.0, 0.0/

	W_VOIGT = 0.0
	y_temp = 0.0
	tail = 0.0
	npks=10
!
	do 10 j = 1 , npks
	gauss_width_0 = (e0(j) / 88.0)**0.5 * Gauss_width
	e_temp = e0(j)
!	
	lor_width = lor_w(j)
	y_temp = i0(j) * czosnyka(energy,e_temp,lor_width,gauss_width_0)
!
	tail = 0
	DO 5 i = 16 , 16
	e_1 = e0(j) + FLOAT(16 - i) * deltae
	tail=tail+tail_x(e_1,energy,lo_tail_a,lo_tail_decay,fast_tail_a,fast_tail_decay,gauss_width_0)
5	CONTINUE
	tail = i0(j) * tail
!
	W_VOIGT = W_VOIGT + y_temp + tail
!
10	continue
!
	return
	end
!
!	*********************************
!
	Function czosnyka(energy, e0, lor_width, gauss_sig)
!
!	X-ray peak shape based on Czonyka et.al.
!
	DOUBLE PRECISION ENERGY, E0, LOR_WIDTH, GAUSS_SIG, czosnyka
	DOUBLE PRECISION ka(3), la(3), ma(3), na(3)
	DOUBLE PRECISION A, V, V_SUM, X, Y, TOT
!
	ka(1) = 0.0
	la(1) = 0.0
	ma(1) = 1.32272
	na(1) = 0.081905
!
	ka(2) = 1.09148
	la(2) = 0.090227
	ma(2) = 1.29081
	na(2) = 0.0093116
!
	ka(3) = 2.30556
	la(3) = 0.0035776
	ma(3) = 1.17417
	na(3) = -0.0116099
!
	a = lor_width / 2 / gauss_sig
	v = (energy - e0) / gauss_sig
	v_sum = 0.0
!
	DO 10 i = 1 ,3
	x = (ka(i) * la(i) + na(i) * (a + ma(i)))
	y = (ka(i)**2 + (a + ma(i))**2 + v**2)
	Z = (ka(i)**2+(a+ma(i))**2 + v**2)**2-4*ka(i)**2*v**2
	tot = (x * y - 2 * ka(i) * la(i) * v**2) / Z
	v_sum = v_sum + tot
10	CONTINUE
!
	czosnyka = v_sum / a / gauss_sig**2
	return
	end
!
!	*********************************
!
	Function tail_x(energy,e0,lo_tail_a,lo_tail_decay,fast_tail_a,fast_tail_decay,Gauss_width)
!
!	calculates peak tail shape based on Gunnink
!
	DOUBLE PRECISION lo_tail_a,lo_tail_decay,fast_tail_a,fast_tail_decay
	DOUBLE PRECISION TAIL_X, ENERGY, E0, T_X, LO_LAMBA, high_lamba
	DOUBLE PRECISION Gauss_width
!
	tail_x = 0.0D0
	If (energy .GT. e0) t_x = 0
	If (energy .GT. e0) GoTo 10
	lo_lambda = Abs(lo_tail_decay)
	high_lamba = Abs(fast_tail_decay)
	ALPHA = -1 / 2 * Gauss_width**2
	t_x = (lo_tail_a*Exp(-lo_lambda*(e0 - energy))+fast_tail_a*Exp(-high_lamba*(e0-energy)))
	t_x = t_x * (1 - Exp(0.4 * ALPHA * (e0 - energy)))
10	continue
	tail_x = t_x
101	continue
	return
	end
!
!	*********************************
!
    function char_x1(E_0, k)
!	Generatres characteristic W X-rays
    double precision R, S_inv, f_abs, rho_z_bar, rho_zm, nu, U0, JAE, m_z
    double precision par_z, E_0, EX0(9), I0(9), par_A, z_k, b_k, NUM
    double precision omega, I_A, t, omega_k, p_k, const1, char_x1
    double precision kfit_consts(36)
    DOUBLE PRECISION xdat(4100), scat_dat(4100), exsp
    INTEGER i_conc(7), i_attn(7)
    INTEGER k
    COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
    DATA EX0/59.318,57.981,57.425,67.244,69.067,66.950,62.273,67.685,69.484/
	DATA I0/0.47,0.274,0.000206,0.103,0.0358,0.0535,0.0006,0.00241,0.0051/
    par_z=74.0
    par_A=183.0
    omega_k=1.0
    omega=1.0
    const1=1.0
    I_A=15.0
    t=1.0
!
    z_k=2.0
    b_k=0.35
    JAE=0.135*par_z
    m_z=0.132-0.9211/par_z**0.5
    U0=E_0/EX0(k)
    nu=E_0**m_z*(0.1904-0.2236*log(par_z)+0.1292*(log(par_z))**2-0.0149*(log(par_z))**3)
    rho_zm=par_a/par_z*(0.787E-5*JAE**0.5*E_0**1.5 + 0.735E-6*E_0**2)
    rho_z_bar=rho_zm*log(U0)*(0.49269-0.10987*nu+0.78857*nu**2)/(0.70256-1.09865*nu+1.0046*nu**2+log(U0))
    S_inv=z_k*b_k/par_z*(U0*log(U0)+1-U0)*(1+16.05*sqrt(JAE/EX0(k))*(sqrt(U0)*log(U0)+2.0*(1-sqrt(U0)))/(U0*log(U0)+1.0-U0))
    R=1-0.0081517*par_z+3.613E-5*par_z**2+0.009583*par_z*exp(-U0)+0.001141*E_0
!
	Pi = 3.14159265358979
	phi = (90.0 - KFIT_CONSTS(7)) / 180.0 * Pi
	eps = KFIT_CONSTS(8) / 180.0 * Pi
	exsp=ex0(k)
    xtemp=1.0-EXP(-attn_fact_W(exsp)*2*RHO_Z_BAR*SIN(PHI)/SIN(EPS))
    f_abs=xtemp/(attn_fact_W(exsp)*2.0*rho_z_bar*SIN(PHI)/SIN(EPS))
!
    Num=const1*omega*I_A*t*R*S_inv*omega_k*I0(k)*f_abs
    char_x1=num
    return

    end
!
!	*********************************
!
        SUBROUTINE Est_error(kfit_parms, kfit_errors, m_dex)
!
!	Adjust uncertainties to compensate for collapse in curfit
!
    DOUBLE PRECISION functn2, FCHISQ1, chi_init
    DOUBLE PRECISION conc_err_array(6), kfit_errors(36), yfit1(4100), weight(4100)
    DOUBLE PRECISION kfit_parms(36), kfit_consts(36), yfit(4100)
    DOUBLE PRECISION xdat(4100), ERROR1(4100), scat_dat(4100)
    DOUBLE PRECISION raw_data(4100),RAND_TEMP(4100), ydat(4100)
    INTEGER i_conc(7), i_attn(7), chan_num, CHAN_1
    COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat
    COMMON/KEDDAT2/RAW_DATA, RAND_TEMP, ERROR1
    COMMON/KEDDAT3/YDAT, weight
!
!
    DOUBLE PRECISION e_start, deltae, gauss_width
    DOUBLE PRECISION p_i, p_e, ref_e, escape_a, p_w
    DOUBLE PRECISION temp_array(201),e_escape(2)
    DOUBLE PRECISION ener_00, y, ener_0,i_escape(2)
    DOUBLE PRECISION E_T, GAUS, FUNC_K, ESC_1, ESC_2
    DOUBLE PRECISION Y_J, T_B, ENER_11, TEMP, TEMP_VAL
    DOUBLE PRECISION W_1,W_2, SCAT_ANGLE, E_MIN, CD_PEAK
    DOUBLE PRECISION scat_p1, scat_p2, scat_i, ANGLES
    DOUBLE PRECISION SCAT_E, SS_THICK, SCATTER, ANG_MIN
    DOUBLE PRECISION ENER10, ENER11, TAIL1, TAIL2
    DOUBLE PRECISION lo_tail_area,lo_tail_decay
    DOUBLE PRECISION fast_tail_area,fast_tail_decay
    DOUBLE PRECISION STEP_SIZE, BKG(4100), ENERGY
    DOUBLE PRECISION x_AREA, E_0, TAILCD, CD_TAIL_I, CD_TAIL_T
    DOUBLE PRECISION chi_new, chi_target
    INTEGER NMIN, NMAX, i_dex, j, num_free, j_dex, iiii, m_dex
!
	NMIN=INT((kfit_consts(21)-kfit_parms(13))/kfit_parms(14))
	NMAX=INT((kfit_consts(22)-kfit_parms(13))/kfit_parms(14))
	npts=1+NMAX-NMIN
!	num_free=npts-20

!	NMIN=972
!	NMAX=1642
!	npts=NMAX-NMIN+1
	num_free=npts-20

        DO 5 i=1,36
	TEMP_ARRAY(i)=KFIT_PARMS(i)
5        continue
	
  	DO 10 i = NMIN , NMAX
        yfit(i) = functn2(i, kfit_parms)
	jjjj=1
10	continue

	chi_init = FCHISQ1(ydat, weight, npts, num_free, NMIN, yfit)
	chi_target=(1.0+chi_init**2)**0.5

	j=1
	j_dex=1
	DO 100 j=1,30

	do 30 i_dex=1, 36
	kfit_parms(i_dex)=TEMP_ARRAY(i_dex)+kfit_errors(i_dex)*float(j)
30	continue

	j_dex=j
  	DO 50 k = NMIN , NMAX
50      yfit1(k) = functn2(k,kfit_parms)

        chi_new=FCHISQ1(ydat, weight, npts, num_free, NMIN, yfit1)

	if(chi_new .GE. chi_target) GOTO 110
100	continue

110	continue


!        DO 120 iiii=1,36
!	if(j_dex.GT.1) kfit_errors(i_dex)=kfit_errors(i_dex)*float(j_dex-1)
!     120	KFIT_PARMS(iiii)=TEMP_ARRAY(iiii)
!     150	continue


        DO 1000 i=1,36
	m_dex=j_dex-1
	if (j_dex .LE. 1) m_dex=1
	kfit_errors(i)=kfit_errors(i)*float(m_dex)
1000	KFIT_PARMS(i)=TEMP_ARRAY(i)


	return
	end
!
!
!        -------------------------------------------------------------
!
        FUNCTION UPUXrays(chan_num, kfit_parms)
	DOUBLE PRECISION UPUXrays, VOIGT
        DOUBLE PRECISION kfit_parms(36), kfit_consts(36), xdat(4100), scat_dat(4100)
	DOUBLE PRECISION xrf_e_U(10), xrf_e_Pu(10), xrf_absU(10), xrf_absPu(10)
        DOUBLE PRECISION lor_w_U(10), lor_w_Pu(10)
	DOUBLE PRECISION ener_00, e_start, deltae, Gauss_width, tails(4)
	DOUBLE PRECISION HPGE_eff, u_pk, pb_pk
        INTEGER i_conc(7), i_attn(7), chan_num, CHAN_1

        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat

! 	Generates U and Pu X-ray peaks 

!        U xrf energies
        xrf_e_U(1) = 98.43158    !  Ka1
        xrf_e_U(2) = 94.65084    !  Ka2
        xrf_e_U(3) = 93.83986    !  Ka3
        xrf_e_U(4) = 111.29508   !  kb1
        xrf_e_U(5) = 114.507     !  kb2
        xrf_e_U(6) = 110.4167    !  Kb3
        xrf_e_U(7) = 115.011     !  kb4
        xrf_e_U(8) = 112.009     !  kb5
        xrf_e_U(9) = 115.377     !  KO23
        xrf_e_U(10) = 115.58     !  KP23

!        Pu xrf energies 
        xrf_e_Pu(1) = 103.734    !  Ka1
        xrf_e_Pu(2) = 99.5232    !  Ka2
        xrf_e_Pu(3) = 98.6815    !  Ka3
        xrf_e_Pu(4) = 117.2437   !  Kb1
        xrf_e_Pu(5) = 120.4376   !  kb2
        xrf_e_Pu(6) = 116.2518   !  kb3
        xrf_e_Pu(7) = 120.9437   !  kb4
        xrf_e_Pu(8) = 117.91665  !  kb5
        xrf_e_Pu(9) = 121.543    !  KO23
        xrf_e_Pu(10) = 121.7668  !  KP23

!         U xrf yields
        xrf_absU(1) = 45.1
        xrf_absU(2) = 28.2
        xrf_absU(3) = 0.099
        xrf_absU(4) = 10.7
        xrf_absU(5) = 4.15
        xrf_absU(6) = 5.65
        xrf_absU(7) = 0.12
        xrf_absU(8) = 0.397
        xrf_absU(9) = 0.95
        xrf_absU(10) = 0.159

!        pu xrf yields
        xrf_absPu(1) = 45.1
        xrf_absPu(2) = 28.4
        xrf_absPu(3) = 0.114
        xrf_absPu(4) = 10.7
        xrf_absPu(5) = 4.18
        xrf_absPu(6) = 5.44
        xrf_absPu(7) = 0.13
        xrf_absPu(8) = 0.413
        xrf_absPu(9) = 0.99
        xrf_absPu(10) = 0.157

!        U lorenztian widths
        lor_w_U(1) = 104.5 / 1000
        lor_w_U(2) = 106.3 / 1000
        lor_w_U(3) = 112.3 / 1000
        lor_w_U(4) = 104.2 / 1000
        lor_w_U(5) = 105.1 / 1000
        lor_w_U(6) = 110.4 / 1000
        lor_w_U(7) = 100.7 / 1000
        lor_w_U(8) = 99.8 / 1000
        lor_w_U(9) = 106.3 / 1000
        lor_w_U(10) = 106.3 / 1000

!        Pu lorentz widths
        lor_w_Pu(1) = 113.06 / 1000
        lor_w_Pu(2) = 115.9 / 1000
        lor_w_Pu(3) = 122.1 / 1000
        lor_w_Pu(4) = 112.2 / 1000
        lor_w_Pu(5) = 113.9 / 1000
        lor_w_Pu(6) = 119.4 / 1000
        lor_w_Pu(7) = 108.95 / 1000
        lor_w_Pu(8) = 108.12 / 1000
        lor_w_Pu(9) = 114.4 / 1000
        lor_w_Pu(10) = 114.4 / 1000
!
	e_start = kfit_parms(13)
	deltae = kfit_parms(14)
	Gauss_width= kfit_parms(15)
	tails(1) = Abs(kfit_parms(22))
	tails(2) = Abs(kfit_parms(23))
	tails(3) = Abs(kfit_parms(24))
	tails(4) = Abs(kfit_parms(25))

	ener_00 = e_start + chan_num * deltae


        u_pk = kfit_parms(30) * VOIGT(ener_00, 10, Gauss_width, xrf_e_U, xrf_absU, tails, deltae, lor_w_U)
        u_pk = u_pk * HPGE_eff(ener_00)/HPGE_eff(xrf_e_U(1))
        pu_pk = kfit_parms(31) * VOIGT(ener_00, 10, Gauss_width, xrf_e_Pu, xrf_absPu, tails, deltae, lor_w_Pu)
        pu_pk = pu_pk * HPGE_eff(ener_00)/HPGE_eff(xrf_e_Pu(1))

	UPUXrays=u_pk + pu_pk
	return
	end
!
!        -------------------------------------------------------------
!
!
!
!        -------------------------------------------------------------
!
        FUNCTION Pb_Xrays(chan_num, kfit_parms)
	DOUBLE PRECISION Pb_Xrays, VOIGT
        DOUBLE PRECISION kfit_parms(36), kfit_consts(36), xdat(4100), scat_dat(4100)
	DOUBLE PRECISION xrf_e_pb(10), xrf_abs_pb(10)
        DOUBLE PRECISION lor_w_pb(10), x_num
	DOUBLE PRECISION ener_00, e_start, deltae, Gauss_width, tails(4), HPGE_eff
        INTEGER i_conc(7), i_attn(7), chan_num, CHAN_1

        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat

! 	Generates Pb X-ray peaks 

!        Pb xrf energies
        xrf_e_pb(1) = 74.969    !  Ka1
        xrf_e_pb(2) = 72.805    !  Ka2
        xrf_e_pb(3) = 72.144    !  Ka3
        xrf_e_pb(4) = 84.938    !  kb1
        xrf_e_pb(5) = 87.300    !  kb2
        xrf_e_pb(6) = 84.450    !  Kb3
        xrf_e_pb(7) = 87.580    !  kb4
        xrf_e_pb(8) = 85.470    !  kb5
        xrf_e_pb(9) = 87.911    !  KO23
        xrf_e_pb(10) = 88.003   !  KP23

!         Pb xrf yields
        xrf_abs_pb(1) = 46.2
        xrf_abs_pb(2) = 27.7
        xrf_abs_pb(3) = 0.0428
        xrf_abs_pb(4) = 10.7
        xrf_abs_pb(5) = 3.91
        xrf_abs_pb(6) = 5.58
        xrf_abs_pb(7) = 0.09
        xrf_abs_pb(8) = 0.312
        xrf_abs_pb(9) = 0.70
        xrf_abs_pb(10) = 0.0165

!        Pb lorenztian widths
        lor_w_pb(1) = 30.0 / 1000
        lor_w_pb(2) = 30.0 / 1000
        lor_w_pb(3) = 30.0 / 1000
        lor_w_pb(4) = 30.0 / 1000
        lor_w_pb(5) = 30.0 / 1000
        lor_w_pb(6) = 30.0 / 1000
        lor_w_pb(7) = 30.0 / 1000
        lor_w_pb(8) = 30.0 / 1000
        lor_w_pb(9) = 30.0 / 1000
        lor_w_pb(10) = 30.0 / 1000
!
	e_start = kfit_parms(13)
	deltae = kfit_parms(14)
	Gauss_width= kfit_parms(15)
	tails(1) = Abs(kfit_parms(22))
	tails(2) = Abs(kfit_parms(23))
	tails(3) = Abs(kfit_parms(24))
	tails(4) = Abs(kfit_parms(25))

	x_num=chan_num
	ener_00 = e_start + x_num * deltae

        pb_pk = kfit_parms(36) * VOIGT(ener_00, 10, Gauss_width, xrf_e_pb, xrf_abs_pb, tails, deltae, lor_w_pb)
        pb_pk = pb_pk * (HPGE_eff(ener_00)/HPGE_eff(xrf_e_pb(1)))*deltae

	Pb_Xrays=pb_pk
	return
	end
!
!
!        -------------------------------------------------------------
!
        FUNCTION Bi_Xrays(chan_num, kfit_parms)
	DOUBLE PRECISION Bi_Xrays, VOIGT
        DOUBLE PRECISION kfit_parms(36), kfit_consts(36), xdat(4100), scat_dat(4100)
	DOUBLE PRECISION xrf_e_bi(10), xrf_abs_bi(10)
        DOUBLE PRECISION lor_w_bi(10), HPGE_eff
	DOUBLE PRECISION ener_00, e_start, deltae, Gauss_width, tails(4)
        INTEGER i_conc(7), i_attn(7), chan_num, CHAN_1

        COMMON/KEDDAT/kfit_consts,i_conc,i_attn,xdat,scat_dat

! 	Generates Bi X-ray peaks

!        Bi xrf energies
        xrf_e_bi(1) = 77.107    !  Ka1
        xrf_e_bi(2) = 74.815    !  Ka2
        xrf_e_bi(3) = 74.138    !  Ka3
        xrf_e_bi(4) = 87.349    !  kb1
        xrf_e_bi(5) = 89.784    !  kb2
        xrf_e_bi(6) = 86.830    !  Kb3
        xrf_e_bi(7) = 90.074    !  kb4
        xrf_e_bi(8) = 87.892    !  kb5
        xrf_e_bi(9) = 90.421    !  KO23
        xrf_e_bi(10) = 90.522   !  KP23

!         Bi xrf yields
        xrf_abs_bi(1) = 46.2
        xrf_abs_bi(2) = 27.7
        xrf_abs_bi(3) = 0.0474
        xrf_abs_bi(4) = 10.7
        xrf_abs_bi(5) = 3.93
        xrf_abs_bi(6) = 5.59
        xrf_abs_bi(7) = 0.09
        xrf_abs_bi(8) = 0.321
        xrf_abs_bi(9) = 0.73
        xrf_abs_bi(10) = 0.031

!        Bi lorenztian widths
        lor_w_bi(1) = 30.0 / 1000
        lor_w_bi(2) = 30.0 / 1000
        lor_w_bi(3) = 30.0 / 1000
        lor_w_bi(4) = 30.0 / 1000
        lor_w_bi(5) = 30.0 / 1000
        lor_w_bi(6) = 30.0 / 1000
        lor_w_bi(7) = 30.0 / 1000
        lor_w_bi(8) = 30.0 / 1000
        lor_w_bi(9) = 30.0 / 1000
        lor_w_bi(10) = 30.0 / 1000
!
	e_start = kfit_parms(13)
	deltae = kfit_parms(14)
	Gauss_width= kfit_parms(15)
	tails(1) = Abs(kfit_parms(22))
	tails(2) = Abs(kfit_parms(23))
	tails(3) = Abs(kfit_parms(24))
	tails(4) = Abs(kfit_parms(25))

	x_num=chan_num
	ener_00 = e_start + x_num * deltae

        Bi_pk = kfit_parms(21) * VOIGT(ener_00, 10, Gauss_width, xrf_e_bi, xrf_abs_bi, tails, deltae, lor_w_bi)
        Bi_pk = Bi_pk * HPGE_eff(ener_00)/HPGE_eff(xrf_e_bi(1))*deltae

	Bi_Xrays=Bi_pk
	return
	end
!
!        -------------------------------------------------------------
!
	Function VOIGT(energy, npks, Gauss_width, e0, i0, tails, deltae, lor_w)
	DOUBLE PRECISION VOIGT, czosnyka, lor_width, gauss_width_0, e_temp, tail_x
	DOUBLE PRECISION energy, Gauss_width, e0(10), i0(10), tails(4), deltae, lor_w(10)
	DOUBLE PRECISION y_temp, tail, e_1, lo_tail_area, lo_tail_decay, fast_tail_area, fast_tail_decay
	INTEGER npks, i, j

!	Caulculates Voigt Peak shape using method devloped by Czosnyka with Gunnink's tailing function

	VOIGT = 0.0
	y_temp = 0.0
	tail = 0.0

	lo_tail_area = tails(1)
	lo_tail_decay = tails(2)
	fast_tail_area = tails(3)
	fast_tail_decay = tails(4)

	DO 20 j = 1,  npks
	  gauss_width_0 = Gauss_width * (e0(j) / 88.0) ** 0.5 
	  e_temp = e0(j)
	  lor_width = lor_w(j)
	  y_temp = i0(j) * czosnyka(energy, e_temp, lor_width, gauss_width_0)

	tail = 0
	DO 10  i = 16, 16				! This loop allows finer detail for tailing - to use open the loop to 1, 30 
	e_1 = e0(j) + FLOAT(16 - i) * deltae
	tail = tail + tail_x(e_1, energy, lo_tail_area, lo_tail_decay, fast_tail_area, fast_tail_decay, gauss_width_0)
10	continue	
	tail = i0(j) * tail

	VOIGT = VOIGT + (y_temp) + tail

20	continue

	return
	end
!
!        -------------------------------------------------------------
!
	FUNCTION read_int(str)
	  CHARACTER(*), INTENT(in) :: str
	  CHARACTER(:), ALLOCATABLE :: instr
 

	  instr = adjustl(str)
	  instr = instr(1:VERIFY(instr,'0123456789', .TRUE.)-1)
	  ! if the string doesn't have a leading digit instr will be empty, return a guard value
 	 IF(instr=='') instr = '-999'
 	 READ(instr,*) read_int
	END FUNCTION read_int
