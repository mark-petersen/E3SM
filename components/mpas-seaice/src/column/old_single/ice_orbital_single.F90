!  SVN:$Id: ice_orbital.F90 1175 2017-03-02 19:53:26Z akt $
!=======================================================================

! Orbital parameters computed from date
! author:  Bruce P. Briegleb, NCAR 
!
! 2006 ECH: Converted to free source form (F90)
! 2014 ECH: Moved routines from csm_share/shr_orb_mod.F90

      module ice_orbital

      use ice_kinds_mod
      use ice_constants_colpkg, only: c2, p5, pi, secday
      use ice_warnings_single, only: add_warning

      implicit none
      private
#ifdef CCSMCOUPLED
      public :: compute_coszen
#else
      public :: shr_orb_params, compute_coszen
#endif

!=======================================================================
 
      contains
 
!=======================================================================

! Uses orbital and lat/lon info to compute cosine solar zenith angle
! for the specified date.
!
! author:  Bruce P. Briegleb, NCAR 

      subroutine compute_coszen (tlat,          tlon,     &
                                 calendar_type, days_per_year, &
                                 nextsw_cday,   yday,  sec, &
                                 coszen,        dt)

      use ice_constants_colpkg, only: eccen, mvelpp, lambm0, obliqr, decln, eccf
#ifdef CCSMCOUPLED
      use shr_orb_mod, only: shr_orb_decl
#endif
 
      real (kind=real_kind), intent(in) :: &
         tlat, tlon          ! latitude and longitude (radians)

      character (len=char_len), intent(in) :: &
         calendar_type       ! differentiates Gregorian from other calendars

      integer (kind=int_kind), intent(in) :: &
         days_per_year, &    ! number of days in one year
         sec                 ! elapsed seconds into date

      real (kind=real_kind), intent(in) :: &
         nextsw_cday     , & ! julian day of next shortwave calculation
         yday                ! day of the year

      real (kind=real_kind), intent(inout) :: &
         coszen              ! cosine solar zenith angle 
                             ! negative for sun below horizon
 
      real (kind=real_kind), intent(in) :: &
         dt                  ! thermodynamic time step
      real (kind=8)    :: ydayp1_R8
      ! local variables
      real (kind=real_kind) :: ydayp1 ! day of year plus one time step
      real (kind=8)    :: eccen_R8
      real (kind=8)    :: mvelpp_R8
      real (kind=8)    :: lambm0_R8
      real (kind=8)    :: obliqr_R8
      real (kind=8)    :: decln_R8
      real (kind=8)    :: eccf_R8
! Solar declination for next time step
 
#ifdef CCSMCOUPLED
      if (calendar_type == "GREGORIAN") then
         ydayp1 = min(nextsw_cday, real(days_per_year,kind=real_kind))
      else
         ydayp1 = nextsw_cday
      endif

      !--- update coszen when nextsw_cday valid
      if (ydayp1 > -0.5_real_kind) then
#else
      ydayp1 = yday + sec/secday
#endif
      ydayp1_R8 = real(ydayp1, kind = 8)
      eccen_R8 = real(eccen, kind = 8)
      mvelpp_R8 = real(mvelpp, kind = 8)
      lambm0_R8  = real(lambm0, kind =8)
      obliqr_R8  = real(obliqr, kind = 8)
      decln_R8  = real(decln, kind = 8)
      eccf_R8  = real(eccf, kind = 8)
      call shr_orb_decl(ydayp1_R8, eccen_R8, mvelpp_R8, & 
                        lambm0_R8, obliqr_R8, decln_R8, & 
                        eccf_R8) 
      decln = real(decln_R8, kind=real_kind)
      eccf = real(eccf_R8, kind=real_kind)
      coszen = sin(tlat)*sin(decln) &
             + cos(tlat)*cos(decln) &
             *cos((sec/secday-p5)*c2*pi + tlon) !cos(hour angle)
 
#ifdef CCSMCOUPLED
      endif
#endif

      end subroutine compute_coszen
 
!===============================================================================

#ifndef CCSMCOUPLED
SUBROUTINE shr_orb_params( iyear_AD , eccen , obliq , mvelp    , &
           &               obliqr   , lambm0, mvelpp, log_print, &
                           l_stop, stop_label)

!-------------------------------------------------------------------------------
!
! Calculate earths orbital parameters using Dave Threshers formula which 
! came from Berger, Andre.  1978  "A Simple Algorithm to Compute Long-Term 
! Variations of Daily Insolation".  Contribution 18, Institute of Astronomy 
! and Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
!
!------------------------------Code history-------------------------------------
!
! Original Author: Erik Kluzek
! Date:            Oct/97
!
!-------------------------------------------------------------------------------

   !----------------------------- Arguments ------------------------------------
   integer(int_kind),intent(in)    :: iyear_AD  ! Year to calculate orbit for
   real   (real_kind),intent(inout) :: eccen     ! orbital eccentricity
   real   (real_kind),intent(inout) :: obliq     ! obliquity in degrees
   real   (real_kind),intent(inout) :: mvelp     ! moving vernal equinox long
   real   (real_kind),intent(out)   :: obliqr    ! Earths obliquity in rad
   real   (real_kind),intent(out)   :: lambm0    ! Mean long of perihelion at
                                                   ! vernal equinox (radians)
   real   (real_kind),intent(out)   :: mvelpp    ! moving vernal equinox long
                                                   ! of perihelion plus pi (rad)
   logical(log_kind),intent(in)    :: log_print ! Flags print of status/error

   logical(log_kind),intent(out)   :: l_stop    ! if true, abort model
   character (len=char_len), intent(out) :: stop_label

   !------------------------------ Parameters ----------------------------------
   real   (real_kind),parameter :: SHR_ORB_UNDEF_REAL = 1.e36_real_kind ! undefined real 
   integer(int_kind),parameter :: SHR_ORB_UNDEF_INT  = 2000000000        ! undefined int
   integer(int_kind),parameter :: poblen =47 ! # of elements in series wrt obliquity
   integer(int_kind),parameter :: pecclen=19 ! # of elements in series wrt eccentricity
   integer(int_kind),parameter :: pmvelen=78 ! # of elements in series wrt vernal equinox
   real   (real_kind),parameter :: psecdeg = 1.0_real_kind/3600.0_real_kind ! arc sec to deg conversion

   real   (real_kind) :: degrad = pi/180._real_kind   ! degree to radian conversion factor
   real   (real_kind) :: yb4_1950AD         ! number of years before 1950 AD

   real   (real_kind),parameter :: SHR_ORB_ECCEN_MIN  =   0.0_real_kind ! min value for eccen
   real   (real_kind),parameter :: SHR_ORB_ECCEN_MAX  =   0.1_real_kind ! max value for eccen
   real   (real_kind),parameter :: SHR_ORB_OBLIQ_MIN  = -90.0_real_kind ! min value for obliq
   real   (real_kind),parameter :: SHR_ORB_OBLIQ_MAX  = +90.0_real_kind ! max value for obliq
   real   (real_kind),parameter :: SHR_ORB_MVELP_MIN  =   0.0_real_kind ! min value for mvelp
   real   (real_kind),parameter :: SHR_ORB_MVELP_MAX  = 360.0_real_kind ! max value for mvelp

   character(len=*),parameter :: subname = '(shr_orb_params)'
 
   ! Cosine series data for computation of obliquity: amplitude (arc seconds),
   ! rate (arc seconds/year), phase (degrees).
 
   real   (real_kind), parameter :: obamp(poblen) =  & ! amplitudes for obliquity cos series
   &      (/   -2462.2214466_real_kind, -857.3232075_real_kind, -629.3231835_real_kind,   &
   &            -414.2804924_real_kind, -311.7632587_real_kind,  308.9408604_real_kind,   &
   &            -162.5533601_real_kind, -116.1077911_real_kind,  101.1189923_real_kind,   &
   &             -67.6856209_real_kind,   24.9079067_real_kind,   22.5811241_real_kind,   &
   &             -21.1648355_real_kind,  -15.6549876_real_kind,   15.3936813_real_kind,   &
   &              14.6660938_real_kind,  -11.7273029_real_kind,   10.2742696_real_kind,   &
   &               6.4914588_real_kind,    5.8539148_real_kind,   -5.4872205_real_kind,   &
   &              -5.4290191_real_kind,    5.1609570_real_kind,    5.0786314_real_kind,   &
   &              -4.0735782_real_kind,    3.7227167_real_kind,    3.3971932_real_kind,   &
   &              -2.8347004_real_kind,   -2.6550721_real_kind,   -2.5717867_real_kind,   &
   &              -2.4712188_real_kind,    2.4625410_real_kind,    2.2464112_real_kind,   &
   &              -2.0755511_real_kind,   -1.9713669_real_kind,   -1.8813061_real_kind,   &
   &              -1.8468785_real_kind,    1.8186742_real_kind,    1.7601888_real_kind,   &
   &              -1.5428851_real_kind,    1.4738838_real_kind,   -1.4593669_real_kind,   &
   &               1.4192259_real_kind,   -1.1818980_real_kind,    1.1756474_real_kind,   &
   &              -1.1316126_real_kind,    1.0896928_real_kind/)
 
   real   (real_kind), parameter :: obrate(poblen) = & ! rates for obliquity cosine series
   &        (/  31.609974_real_kind, 32.620504_real_kind, 24.172203_real_kind,   &
   &            31.983787_real_kind, 44.828336_real_kind, 30.973257_real_kind,   &
   &            43.668246_real_kind, 32.246691_real_kind, 30.599444_real_kind,   &
   &            42.681324_real_kind, 43.836462_real_kind, 47.439436_real_kind,   &
   &            63.219948_real_kind, 64.230478_real_kind,  1.010530_real_kind,   &
   &             7.437771_real_kind, 55.782177_real_kind,  0.373813_real_kind,   &
   &            13.218362_real_kind, 62.583231_real_kind, 63.593761_real_kind,   &
   &            76.438310_real_kind, 45.815258_real_kind,  8.448301_real_kind,   &
   &            56.792707_real_kind, 49.747842_real_kind, 12.058272_real_kind,   &
   &            75.278220_real_kind, 65.241008_real_kind, 64.604291_real_kind,   &
   &             1.647247_real_kind,  7.811584_real_kind, 12.207832_real_kind,   &
   &            63.856665_real_kind, 56.155990_real_kind, 77.448840_real_kind,   &
   &             6.801054_real_kind, 62.209418_real_kind, 20.656133_real_kind,   &
   &            48.344406_real_kind, 55.145460_real_kind, 69.000539_real_kind,   &
   &            11.071350_real_kind, 74.291298_real_kind, 11.047742_real_kind,   &
   &             0.636717_real_kind, 12.844549_real_kind/)
 
   real   (real_kind), parameter :: obphas(poblen) = & ! phases for obliquity cosine series
   &      (/    251.9025_real_kind, 280.8325_real_kind, 128.3057_real_kind,   &
   &            292.7252_real_kind,  15.3747_real_kind, 263.7951_real_kind,   &
   &            308.4258_real_kind, 240.0099_real_kind, 222.9725_real_kind,   &
   &            268.7809_real_kind, 316.7998_real_kind, 319.6024_real_kind,   &
   &            143.8050_real_kind, 172.7351_real_kind,  28.9300_real_kind,   &
   &            123.5968_real_kind,  20.2082_real_kind,  40.8226_real_kind,   &
   &            123.4722_real_kind, 155.6977_real_kind, 184.6277_real_kind,   &
   &            267.2772_real_kind,  55.0196_real_kind, 152.5268_real_kind,   &
   &             49.1382_real_kind, 204.6609_real_kind,  56.5233_real_kind,   &
   &            200.3284_real_kind, 201.6651_real_kind, 213.5577_real_kind,   &
   &             17.0374_real_kind, 164.4194_real_kind,  94.5422_real_kind,   &
   &            131.9124_real_kind,  61.0309_real_kind, 296.2073_real_kind,   &
   &            135.4894_real_kind, 114.8750_real_kind, 247.0691_real_kind,   &
   &            256.6114_real_kind,  32.1008_real_kind, 143.6804_real_kind,   &
   &             16.8784_real_kind, 160.6835_real_kind,  27.5932_real_kind,   &
   &            348.1074_real_kind,  82.6496_real_kind/)
 
   ! Cosine/sine series data for computation of eccentricity and fixed vernal 
   ! equinox longitude of perihelion (fvelp): amplitude, 
   ! rate (arc seconds/year), phase (degrees).
 
   real   (real_kind), parameter :: ecamp (pecclen) = & ! ampl for eccen/fvelp cos/sin series
   &      (/   0.01860798_real_kind,  0.01627522_real_kind, -0.01300660_real_kind,   &
   &           0.00988829_real_kind, -0.00336700_real_kind,  0.00333077_real_kind,   &
   &          -0.00235400_real_kind,  0.00140015_real_kind,  0.00100700_real_kind,   &
   &           0.00085700_real_kind,  0.00064990_real_kind,  0.00059900_real_kind,   &
   &           0.00037800_real_kind, -0.00033700_real_kind,  0.00027600_real_kind,   &
   &           0.00018200_real_kind, -0.00017400_real_kind, -0.00012400_real_kind,   &
   &           0.00001250_real_kind/)
 
   real   (real_kind), parameter :: ecrate(pecclen) = & ! rates for eccen/fvelp cos/sin series
   &      (/    4.2072050_real_kind,  7.3460910_real_kind, 17.8572630_real_kind,  &
   &           17.2205460_real_kind, 16.8467330_real_kind,  5.1990790_real_kind,  &
   &           18.2310760_real_kind, 26.2167580_real_kind,  6.3591690_real_kind,  &
   &           16.2100160_real_kind,  3.0651810_real_kind, 16.5838290_real_kind,  &
   &           18.4939800_real_kind,  6.1909530_real_kind, 18.8677930_real_kind,  &
   &           17.4255670_real_kind,  6.1860010_real_kind, 18.4174410_real_kind,  &
   &            0.6678630_real_kind/)
 
   real   (real_kind), parameter :: ecphas(pecclen) = & ! phases for eccen/fvelp cos/sin series
   &      (/    28.620089_real_kind, 193.788772_real_kind, 308.307024_real_kind,  &
   &           320.199637_real_kind, 279.376984_real_kind,  87.195000_real_kind,  &
   &           349.129677_real_kind, 128.443387_real_kind, 154.143880_real_kind,  &
   &           291.269597_real_kind, 114.860583_real_kind, 332.092251_real_kind,  &
   &           296.414411_real_kind, 145.769910_real_kind, 337.237063_real_kind,  &
   &           152.092288_real_kind, 126.839891_real_kind, 210.667199_real_kind,  &
   &            72.108838_real_kind/)
 
   ! Sine series data for computation of moving vernal equinox longitude of 
   ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).      
 
   real   (real_kind), parameter :: mvamp (pmvelen) = & ! amplitudes for mvelp sine series 
   &      (/   7391.0225890_real_kind, 2555.1526947_real_kind, 2022.7629188_real_kind,  &
   &          -1973.6517951_real_kind, 1240.2321818_real_kind,  953.8679112_real_kind,  &
   &           -931.7537108_real_kind,  872.3795383_real_kind,  606.3544732_real_kind,  &
   &           -496.0274038_real_kind,  456.9608039_real_kind,  346.9462320_real_kind,  &
   &           -305.8412902_real_kind,  249.6173246_real_kind, -199.1027200_real_kind,  &
   &            191.0560889_real_kind, -175.2936572_real_kind,  165.9068833_real_kind,  &
   &            161.1285917_real_kind,  139.7878093_real_kind, -133.5228399_real_kind,  &
   &            117.0673811_real_kind,  104.6907281_real_kind,   95.3227476_real_kind,  &
   &             86.7824524_real_kind,   86.0857729_real_kind,   70.5893698_real_kind,  &
   &            -69.9719343_real_kind,  -62.5817473_real_kind,   61.5450059_real_kind,  &
   &            -57.9364011_real_kind,   57.1899832_real_kind,  -57.0236109_real_kind,  &
   &            -54.2119253_real_kind,   53.2834147_real_kind,   52.1223575_real_kind,  &
   &            -49.0059908_real_kind,  -48.3118757_real_kind,  -45.4191685_real_kind,  &
   &            -42.2357920_real_kind,  -34.7971099_real_kind,   34.4623613_real_kind,  &
   &            -33.8356643_real_kind,   33.6689362_real_kind,  -31.2521586_real_kind,  &
   &            -30.8798701_real_kind,   28.4640769_real_kind,  -27.1960802_real_kind,  &
   &             27.0860736_real_kind,  -26.3437456_real_kind,   24.7253740_real_kind,  &
   &             24.6732126_real_kind,   24.4272733_real_kind,   24.0127327_real_kind,  &
   &             21.7150294_real_kind,  -21.5375347_real_kind,   18.1148363_real_kind,  &
   &            -16.9603104_real_kind,  -16.1765215_real_kind,   15.5567653_real_kind,  &
   &             15.4846529_real_kind,   15.2150632_real_kind,   14.5047426_real_kind,  &
   &            -14.3873316_real_kind,   13.1351419_real_kind,   12.8776311_real_kind,  &
   &             11.9867234_real_kind,   11.9385578_real_kind,   11.7030822_real_kind,  &
   &             11.6018181_real_kind,  -11.2617293_real_kind,  -10.4664199_real_kind,  &
   &             10.4333970_real_kind,  -10.2377466_real_kind,   10.1934446_real_kind,  &
   &            -10.1280191_real_kind,   10.0289441_real_kind,  -10.0034259_real_kind/)
 
   real   (real_kind), parameter :: mvrate(pmvelen) = & ! rates for mvelp sine series 
   &      (/    31.609974_real_kind, 32.620504_real_kind, 24.172203_real_kind,   &
   &             0.636717_real_kind, 31.983787_real_kind,  3.138886_real_kind,   &
   &            30.973257_real_kind, 44.828336_real_kind,  0.991874_real_kind,   &
   &             0.373813_real_kind, 43.668246_real_kind, 32.246691_real_kind,   &
   &            30.599444_real_kind,  2.147012_real_kind, 10.511172_real_kind,   &
   &            42.681324_real_kind, 13.650058_real_kind,  0.986922_real_kind,   &
   &             9.874455_real_kind, 13.013341_real_kind,  0.262904_real_kind,   &
   &             0.004952_real_kind,  1.142024_real_kind, 63.219948_real_kind,   &
   &             0.205021_real_kind,  2.151964_real_kind, 64.230478_real_kind,   &
   &            43.836462_real_kind, 47.439436_real_kind,  1.384343_real_kind,   &
   &             7.437771_real_kind, 18.829299_real_kind,  9.500642_real_kind,   &
   &             0.431696_real_kind,  1.160090_real_kind, 55.782177_real_kind,   &
   &            12.639528_real_kind,  1.155138_real_kind,  0.168216_real_kind,   &
   &             1.647247_real_kind, 10.884985_real_kind,  5.610937_real_kind,   &
   &            12.658184_real_kind,  1.010530_real_kind,  1.983748_real_kind,   &
   &            14.023871_real_kind,  0.560178_real_kind,  1.273434_real_kind,   &
   &            12.021467_real_kind, 62.583231_real_kind, 63.593761_real_kind,   &
   &            76.438310_real_kind,  4.280910_real_kind, 13.218362_real_kind,   &
   &            17.818769_real_kind,  8.359495_real_kind, 56.792707_real_kind,   &
   &            8.448301_real_kind,  1.978796_real_kind,  8.863925_real_kind,   &
   &             0.186365_real_kind,  8.996212_real_kind,  6.771027_real_kind,   &
   &            45.815258_real_kind, 12.002811_real_kind, 75.278220_real_kind,   &
   &            65.241008_real_kind, 18.870667_real_kind, 22.009553_real_kind,   &
   &            64.604291_real_kind, 11.498094_real_kind,  0.578834_real_kind,   &
   &             9.237738_real_kind, 49.747842_real_kind,  2.147012_real_kind,   &
   &             1.196895_real_kind,  2.133898_real_kind,  0.173168_real_kind/)

   real   (real_kind), parameter :: mvphas(pmvelen) = & ! phases for mvelp sine series
   &      (/    251.9025_real_kind, 280.8325_real_kind, 128.3057_real_kind,   &
   &            348.1074_real_kind, 292.7252_real_kind, 165.1686_real_kind,   &
   &            263.7951_real_kind,  15.3747_real_kind,  58.5749_real_kind,   &
   &             40.8226_real_kind, 308.4258_real_kind, 240.0099_real_kind,   &
   &            222.9725_real_kind, 106.5937_real_kind, 114.5182_real_kind,   &
   &            268.7809_real_kind, 279.6869_real_kind,  39.6448_real_kind,   &
   &            126.4108_real_kind, 291.5795_real_kind, 307.2848_real_kind,   &
   &             18.9300_real_kind, 273.7596_real_kind, 143.8050_real_kind,   &
   &            191.8927_real_kind, 125.5237_real_kind, 172.7351_real_kind,   &
   &            316.7998_real_kind, 319.6024_real_kind,  69.7526_real_kind,   &
   &            123.5968_real_kind, 217.6432_real_kind,  85.5882_real_kind,   &
   &            156.2147_real_kind,  66.9489_real_kind,  20.2082_real_kind,   &
   &            250.7568_real_kind,  48.0188_real_kind,   8.3739_real_kind,   &
   &             17.0374_real_kind, 155.3409_real_kind,  94.1709_real_kind,   &
   &            221.1120_real_kind,  28.9300_real_kind, 117.1498_real_kind,   &
   &            320.5095_real_kind, 262.3602_real_kind, 336.2148_real_kind,   &
   &            233.0046_real_kind, 155.6977_real_kind, 184.6277_real_kind,   &
   &            267.2772_real_kind,  78.9281_real_kind, 123.4722_real_kind,   &
   &            188.7132_real_kind, 180.1364_real_kind,  49.1382_real_kind,   &
   &            152.5268_real_kind,  98.2198_real_kind,  97.4808_real_kind,   &
   &            221.5376_real_kind, 168.2438_real_kind, 161.1199_real_kind,   &
   &             55.0196_real_kind, 262.6495_real_kind, 200.3284_real_kind,   &
   &            201.6651_real_kind, 294.6547_real_kind,  99.8233_real_kind,   &
   &            213.5577_real_kind, 154.1631_real_kind, 232.7153_real_kind,   &
   &            138.3034_real_kind, 204.6609_real_kind, 106.5938_real_kind,   &
   &            250.4676_real_kind, 332.3345_real_kind,  27.3039_real_kind/)
 
   !---------------------------Local variables----------------------------------
   integer(int_kind) :: i       ! Index for series summations
   real   (real_kind) :: obsum   ! Obliquity series summation
   real   (real_kind) :: cossum  ! Cos series summation for eccentricity/fvelp
   real   (real_kind) :: sinsum  ! Sin series summation for eccentricity/fvelp
   real   (real_kind) :: fvelp   ! Fixed vernal equinox long of perihelion
   real   (real_kind) :: mvsum   ! mvelp series summation
   real   (real_kind) :: beta    ! Intermediate argument for lambm0
   real   (real_kind) :: years   ! Years to time of interest ( pos <=> future)
   real   (real_kind) :: eccen2  ! eccentricity squared
   real   (real_kind) :: eccen3  ! eccentricity cubed
   integer (int_kind), parameter :: s_loglev    = 0         
   character(len=char_len_long) :: warning ! warning message

   !-------------------------- Formats -----------------------------------------
   character(*),parameter :: svnID  = "SVN " // &
   "$Id: ice_orbital.F90 1175 2017-03-02 19:53:26Z akt $"
   character(*),parameter :: svnURL = "SVN <unknown URL>" 
!  character(*),parameter :: svnURL = "SVN " // &
!  "$URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_121022/shr/shr_orb_mod.F90 $"
   character(len=*),parameter :: F00 = "('(shr_orb_params) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_orb_params) ',a,i9)"
   character(len=*),parameter :: F02 = "('(shr_orb_params) ',a,f6.3)"
   character(len=*),parameter :: F03 = "('(shr_orb_params) ',a,es14.6)"

   !----------------------------------------------------------------------------
   ! radinp and algorithms below will need a degree to radian conversion factor

   l_stop = .false.
   stop_label = ' '
 
   if ( log_print .and. s_loglev > 0 ) then
     write(warning,F00) 'Calculate characteristics of the orbit:'
     call add_warning(warning)
     write(warning,F00) svnID
     call add_warning(warning)
!    write(warning,F00) svnURL
!    call add_warning(warning)
   end if
 
   ! Check for flag to use input orbit parameters
 
   IF ( iyear_AD == SHR_ORB_UNDEF_INT ) THEN

      ! Check input obliq, eccen, and mvelp to ensure reasonable
 
      if( obliq == SHR_ORB_UNDEF_REAL )then
         write(warning,F00) trim(subname)//' Have to specify orbital parameters:'
         call add_warning(warning)
         write(warning,F00) 'Either set: iyear_AD, OR [obliq, eccen, and mvelp]:'
         call add_warning(warning)
         write(warning,F00) 'iyear_AD is the year to simulate orbit for (ie. 1950): '
         call add_warning(warning)
         write(warning,F00) 'obliq, eccen, mvelp specify the orbit directly:'
         call add_warning(warning)
         write(warning,F00) 'The AMIP II settings (for a 1995 orbit) are: '
         call add_warning(warning)
         write(warning,F00) ' obliq =  23.4441'
         call add_warning(warning)
         write(warning,F00) ' eccen =   0.016715'
         call add_warning(warning)
         write(warning,F00) ' mvelp = 102.7'
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable oblip'
      else if ( log_print ) then
         write(warning,F00) 'Use input orbital parameters: '
         call add_warning(warning)
      end if
      if( (obliq < SHR_ORB_OBLIQ_MIN).or.(obliq > SHR_ORB_OBLIQ_MAX) ) then
         write(warning,F03) 'Input obliquity unreasonable: ', obliq
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable obliq'
      end if
      if( (eccen < SHR_ORB_ECCEN_MIN).or.(eccen > SHR_ORB_ECCEN_MAX) ) then
         write(warning,F03) 'Input eccentricity unreasonable: ', eccen
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable eccen'
      end if
      if( (mvelp < SHR_ORB_MVELP_MIN).or.(mvelp > SHR_ORB_MVELP_MAX) ) then
         write(warning,F03) 'Input mvelp unreasonable: ' , mvelp
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable mvelp'
      end if
      eccen2 = eccen*eccen
      eccen3 = eccen2*eccen

   ELSE  ! Otherwise calculate based on years before present
 
      if ( log_print .and. s_loglev > 0) then
         write(warning,F01) 'Calculate orbit for year: ' , iyear_AD
         call add_warning(warning)
      end if
      yb4_1950AD = 1950.0_real_kind - real(iyear_AD,real_kind)
      if ( abs(yb4_1950AD) .gt. 1000000.0_real_kind )then
         write(warning,F00) 'orbit only valid for years+-1000000'
         call add_warning(warning)
         write(warning,F00) 'Relative to 1950 AD'
         call add_warning(warning)
         write(warning,F03) '# of years before 1950: ',yb4_1950AD
         call add_warning(warning)
         write(warning,F01) 'Year to simulate was  : ',iyear_AD
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable year'
      end if
 
      ! The following calculates the earths obliquity, orbital eccentricity
      ! (and various powers of it) and vernal equinox mean longitude of
      ! perihelion for years in the past (future = negative of years past),
      ! using constants (see parameter section) given in the program of:
      !
      ! Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term Variations
      ! of Daily Insolation.  Contribution 18, Institute of Astronomy and
      ! Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium.
      !
      ! and formulas given in the paper (where less precise constants are also
      ! given):
      !
      ! Berger, Andre.  1978.  Long-Term Variations of Daily Insolation and
      ! Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367
      !
      ! The algorithm is valid only to 1,000,000 years past or hence.
      ! For a solution valid to 5-10 million years past see the above author.
      ! Algorithm below is better for years closer to present than is the
      ! 5-10 million year solution.
      !
      ! Years to time of interest must be negative of years before present
      ! (1950) in formulas that follow. 
 
      years = - yb4_1950AD
 
      ! In the summations below, cosine or sine arguments, which end up in
      ! degrees, must be converted to radians via multiplication by degrad.
      !
      ! Summation of cosine series for obliquity (epsilon in Berger 1978) in
      ! degrees. Convert the amplitudes and rates, which are in arc secs, into
      ! degrees via multiplication by psecdeg (arc seconds to degrees conversion
      ! factor).  For obliq, first term is Berger 1978 epsilon star; second
      ! term is series summation in degrees.
  
      obsum = 0.0_real_kind
      do i = 1, poblen
         obsum = obsum + obamp(i)*psecdeg*cos((obrate(i)*psecdeg*years + &
         &       obphas(i))*degrad)
      end do
      obliq = 23.320556_real_kind + obsum
 
      ! Summation of cosine and sine series for computation of eccentricity 
      ! (eccen; e in Berger 1978) and fixed vernal equinox longitude of 
      ! perihelion (fvelp; pi in Berger 1978), which is used for computation 
      ! of moving vernal equinox longitude of perihelion.  Convert the rates, 
      ! which are in arc seconds, into degrees via multiplication by psecdeg.
 
      cossum = 0.0_real_kind
      do i = 1, pecclen
        cossum = cossum+ecamp(i)*cos((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do
 
      sinsum = 0.0_real_kind
      do i = 1, pecclen
        sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do
 
      ! Use summations to calculate eccentricity
 
      eccen2 = cossum*cossum + sinsum*sinsum
      eccen  = sqrt(eccen2)
      eccen3 = eccen2*eccen
 
      ! A series of cases for fvelp, which is in radians.
         
      if (abs(cossum) .le. 1.0E-8_real_kind) then
        if (sinsum .eq. 0.0_real_kind) then
          fvelp = 0.0_real_kind
        else if (sinsum .lt. 0.0_real_kind) then
          fvelp = 1.5_real_kind*pi
        else if (sinsum .gt. 0.0_real_kind) then
          fvelp = .5_real_kind*pi
        endif
      else if (cossum .lt. 0.0_real_kind) then
        fvelp = atan(sinsum/cossum) + pi
      else if (cossum .gt. 0.0_real_kind) then
        if (sinsum .lt. 0.0_real_kind) then
          fvelp = atan(sinsum/cossum) + 2.0_real_kind*pi
        else
          fvelp = atan(sinsum/cossum)
        endif
      endif
 
      ! Summation of sin series for computation of moving vernal equinox long
      ! of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For mvelp,
      ! first term is fvelp in degrees; second term is Berger 1978 psi bar 
      ! times years and in degrees; third term is Berger 1978 zeta; fourth 
      ! term is series summation in degrees.  Convert the amplitudes and rates,
      ! which are in arc seconds, into degrees via multiplication by psecdeg.  
      ! Series summation plus second and third terms constitute Berger 1978
      ! psi, which is the general precession.
 
      mvsum = 0.0_real_kind
      do i = 1, pmvelen
        mvsum = mvsum + mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*years + &
        &       mvphas(i))*degrad)
      end do
      mvelp = fvelp/degrad + 50.439273_real_kind*psecdeg*years + 3.392506_real_kind + mvsum
 
      ! Cases to make sure mvelp is between 0 and 360.
 
      do while (mvelp .lt. 0.0_real_kind)
        mvelp = mvelp + 360.0_real_kind
      end do
      do while (mvelp .ge. 360.0_real_kind)
        mvelp = mvelp - 360.0_real_kind
      end do

   END IF  ! end of test on whether to calculate or use input orbital params
 
   ! Orbit needs the obliquity in radians
 
   obliqr = obliq*degrad
 
   ! 180 degrees must be added to mvelp since observations are made from the
   ! earth and the sun is considered (wrongly for the algorithm) to go around
   ! the earth. For a more graphic explanation see Appendix B in:
   !
   ! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth Orbital
   ! Periods.  J. of Geophysical Research 98:10,341-10,362.
   !
   ! Additionally, orbit will need this value in radians. So mvelp becomes
   ! mvelpp (mvelp plus pi)
 
   mvelpp = (mvelp + 180._real_kind)*degrad
 
   ! Set up an argument used several times in lambm0 calculation ahead.
 
   beta = sqrt(1._real_kind - eccen2)
 
   ! The mean longitude at the vernal equinox (lambda m nought in Berger
   ! 1978; in radians) is calculated from the following formula given in 
   ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
   ! 1978) is 0.

   lambm0 = 2._real_kind*((.5_real_kind*eccen + .125_real_kind*eccen3)*(1._real_kind + beta)*sin(mvelpp)  &
   &      - .250_real_kind*eccen2*(.5_real_kind    + beta)*sin(2._real_kind*mvelpp)            &
   &      + .125_real_kind*eccen3*(1._real_kind/3._real_kind + beta)*sin(3._real_kind*mvelpp))
 
   if ( log_print ) then
     write(warning,F03) '------ Computed Orbital Parameters ------'
     call add_warning(warning)
     write(warning,F03) 'Eccentricity      = ',eccen
     call add_warning(warning)
     write(warning,F03) 'Obliquity (deg)   = ',obliq
     call add_warning(warning)
     write(warning,F03) 'Obliquity (rad)   = ',obliqr
     call add_warning(warning)
     write(warning,F03) 'Long of perh(deg) = ',mvelp
     call add_warning(warning)
     write(warning,F03) 'Long of perh(rad) = ',mvelpp
     call add_warning(warning)
     write(warning,F03) 'Long at v.e.(rad) = ',lambm0
     call add_warning(warning)
     write(warning,F03) '-----------------------------------------'
     call add_warning(warning)
   end if
 
END SUBROUTINE shr_orb_params

!===============================================================================

SUBROUTINE shr_orb_decl(calday ,eccen ,mvelpp ,lambm0 ,obliqr ,delta ,eccf)

!-------------------------------------------------------------------------------
!
! Compute earth/orbit parameters using formula suggested by
! Duane Thresher.
!
!---------------------------Code history----------------------------------------
!
! Original version:  Erik Kluzek
! Date:              Oct/1997
!
!-------------------------------------------------------------------------------

   !------------------------------Arguments--------------------------------
   real   (real_kind),intent(in)  :: calday ! Calendar day, including fraction
   real   (real_kind),intent(in)  :: eccen  ! Eccentricity
   real   (real_kind),intent(in)  :: obliqr ! Earths obliquity in radians
   real   (real_kind),intent(in)  :: lambm0 ! Mean long of perihelion at the 
                                              ! vernal equinox (radians)
   real   (real_kind),intent(in)  :: mvelpp ! moving vernal equinox longitude
                                              ! of perihelion plus pi (radians)
   real   (real_kind),intent(out) :: delta  ! Solar declination angle in rad
   real   (real_kind),intent(out) :: eccf   ! Earth-sun distance factor (ie. (1/r)**2)
 
   !---------------------------Local variables-----------------------------
   real   (real_kind),parameter :: dayspy = 365.0_real_kind  ! days per year
   real   (real_kind),parameter :: ve     = 80.5_real_kind   ! Calday of vernal equinox
                                                     ! assumes Jan 1 = calday 1
 
   real   (real_kind) ::   lambm  ! Lambda m, mean long of perihelion (rad)
   real   (real_kind) ::   lmm    ! Intermediate argument involving lambm
   real   (real_kind) ::   lamb   ! Lambda, the earths long of perihelion
   real   (real_kind) ::   invrho ! Inverse normalized sun/earth distance
   real   (real_kind) ::   sinl   ! Sine of lmm
 
   ! Compute eccentricity factor and solar declination using
   ! day value where a round day (such as 213.0) refers to 0z at
   ! Greenwich longitude.
   !
   ! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
   ! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
   ! 35:2362-2367.
   !
   ! To get the earths true longitude (position in orbit; lambda in Berger 
   ! 1978) which is necessary to find the eccentricity factor and declination,
   ! must first calculate the mean longitude (lambda m in Berger 1978) at
   ! the present day.  This is done by adding to lambm0 (the mean longitude
   ! at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
   ! an increment (delta lambda m in Berger 1978) that is the number of
   ! days past or before (a negative increment) the vernal equinox divided by
   ! the days in a model year times the 2*pi radians in a complete orbit.
 
   lambm = lambm0 + (calday - ve)*2._real_kind*pi/dayspy
   lmm   = lambm  - mvelpp
 
   ! The earths true longitude, in radians, is then found from
   ! the formula in Berger 1978:
 
   sinl  = sin(lmm)
   lamb  = lambm  + eccen*(2._real_kind*sinl + eccen*(1.25_real_kind*sin(2._real_kind*lmm)  &
   &     + eccen*((13.0_real_kind/12.0_real_kind)*sin(3._real_kind*lmm) - 0.25_real_kind*sinl)))
 
   ! Using the obliquity, eccentricity, moving vernal equinox longitude of
   ! perihelion (plus), and earths true longitude, the declination (delta)
   ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
   ! rho will be used), and thus the eccentricity factor (eccf), can be 
   ! calculated from formulas given in Berger 1978.
 
   invrho = (1._real_kind + eccen*cos(lamb - mvelpp)) / (1._real_kind - eccen*eccen)
 
   ! Set solar declination and eccentricity factor
 
   delta  = asin(sin(obliqr)*sin(lamb))
   eccf   = invrho*invrho
 
   return
 
END SUBROUTINE shr_orb_decl
#endif

!=======================================================================
 
      end module ice_orbital
 
!=======================================================================
