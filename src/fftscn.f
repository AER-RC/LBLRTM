C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      Subroutine FFTSCN (IFILE,JFILE)
C**********************************************************************
C     FFTSCN and its associated subroutines convolve an LBLRTM spectrum
C     with an instrument response function (scanning function) using
C     Fourier Transforms.  The program transforms the spectrum to the
C     Fourier domain, multiplies it by the Fourier transform of the
C     scanning function, and transforms it back to the frequency domain
C     Compared to convolving in the spectral domain, this technique
C     has the advantages that it is more accurate and more efficient 
C     for large spectral regions.
C
C     For further details, see:
C         FFTSCN: A Program for Spectral Scanning Using Fourier
C                  Transforms
C         By: William Gallery
C             S. A. Clough
C             AER, Inc., Cambridge, MA
C         Technical Report, Phillips Lab, Geophysics Directorate,
C            PL-TR-92-2131
C
C         (Phillips Lab was formerly the Air Force Geophysics Lab [AFGL])
C
C     Created December, 1989 by:
C                               William Gallery
C                               OptiMetrics, Inc.
C     Current address:          AER, Inc.
C                               Cambridge, MA
C
C     Version 1.0: Dec, 1992
C
C     Version 1.1: June, 1993
C         Changes in Version 1.1:
C         1. No longer is the minimum size of an FFT equal to LPTSMX, 
C            rather it is the smallest power of 2 greater than the 
C            number of data points.  LPTSMX (set in file fftparm.inc)
C            can now be set to a large value consistent with real (not
C            virtual) memory.
C         2. Fixed a bug in HARM1D.  The variables NT, NTV2, and MT were
C            set in one call to HARM1D but were not preserved for subsequent
C            calls.  The problem only occured on some platforms, notably
C            Silicon Graphics workstations.  The problem was solved with a 
C            SAVE statement.
C         3. Added Norton-Beer, Brault, and Kaiser-Bessel Scanning Functions
C
C     Version 1.2: November, 1993
C         Changes in Version 1.2
C         1. Added option to interpolate the final result onto a user-
C            specified grid. 
C         2. If the DV of the input spectrum is greater than HWHM/6, then
C            the input spectrum is first interpolated onto a frequency grid
C            where DV <= HWHM/6.
C         3. Fixed bug in Norton-Beer functions: it always returned the weak
C            function.
C         4. Fixed bug: IFILE is now closed after each scanning function 
C            request.  On some machines (Convex), not doing this and requesting
C            another operation on IFILE generated an error.
C         5. Fixed bug: the output spectral grid was sometimes off by
C            one dv.
C
C     Program instructions:
C     The program commands are contained on a single input record plus
C     up to 3 additional records, depending upon the case.  Multiple 
C     commands may be contained on successive records.  A zero or 
C     negative number in the first field terminates the sequence.
C
C     The format of a command is as follows:
C
C     Field:     HWHM      V1      V2   JEMIT   JFNIN  MRATIN   DVOUT 
C     Column:	 1-10	11-20   21-30   31-35	36-40   41-45   46-55 
C     Field:	IUNIT    IFIL    NFIL   JUNIT     IVX   NOFIX
C     Column:   56-60   61-65   66-70   71-75   76-78   79-80
C
C     Format(3F10.3,3I5,F10.3,4I5,I3,I2)
C
C     HWHM    Half Width at Half Maximum of the scanning function, or 
C             if JFN < 0, then maximum optical path difference of an 
C             equivalent interferometer. If HWHW <= 0, then exit FFTSCN
C     V1      Initial wavenumer for the scanned result
C     V2      Final wavenumber for the scanned result
C     JEMIT   = 0:  convolve with transmittance
C             = 1:  convolve with radiance
C     JFNIN   Selects the Scanning Function 
C             = 0: Boxcar.  Halfwidth is truncated to M du/2, 
C                  where M is an integer
C             = 1: Triangle
C             = 2: Gauss
C             = 3: Sinc2
C             = 4: Sinc
C             = 5: Beer
C             = 6: Hamming
C             = 7: Hanning
C             = 8: Norton-Beer, Weak
C             = 9: Norton-Beer, Moderate
C             =10: norton-Beer, Strong
C             =11: Brault (needs input parameter PARM)
C             =12: Kaiser-Bessel (needs input parameter PARM)
C             =13: Kiruna (asymetric, c1*sinc(u)+c2*sinc(u-u1)
C             If JFN < 0, then HWHM is the maximum optical path 
C             difference of an equivalent interfereometer, apodized to
C             give the scanning function given by |JFN|. 
C     MRAT    For prescanning with a boxcar, the ratio of HWHM of the
C             scanning function to the halfwidth of the boxcar,
C             default = MRATDF (=12). If MRAT < 0, no boxcaring is
C             performed.
C     DVOUT   If DVOUT > 0, then the scanned spectral file is interpolated 
C             onto the grid defined by V1, V2, and DV
C     IUNIT   Unit number of the file containing the spectrum to be 
C             scanned. 
C             = 0: use file on UNIT = IFILE (from LBLRTM, default=12)
C             > 0: use file on UNIT = IUNIT
C             < 0: read a filename from the next record, 20 characters max,
C                  and open this file on UNIT = -IUNIT
C     IFILST  Sequential number of the first LBLRTM file on IUNIT to
C             be scanned
C     NIFILS  Number of LBLRTM files on IUNIT to be scanned, 
C             beginning with IFILST
C     JUNIT   Unit number of the file containing the output spectrum.  
C             = 0: UNIT = IFILE (from LBLRTM default = 11),
C                   filename=TAPExx, where xx = IFILE
C             > 0: UNIT = IUNIT, filename = TAPExx, where xx = IUNIT
C             < 0: read in a filename from the next record, 60 characters max,
C                  and open this file on UNIT = -JUNIT
C     IVX     = -1:  Scanning function is calculated as the FFT of the
C                    Apodization function
C             =  0:  Program decides how to calculate the scanning 
C                    function using the value of CR(JFN).
C             =  1:  Scanning function is calculated analytically
C     
C     NOFIX   For prescanning with a boxcar: if non-zero, then do not
C             deconvolve the scanned spectrum with the boxcar
C
C     Record 2 (For JFN > 10, )
C         PARM (F10.3) parameter defining the scanning function
C
C     Record 3 (for IUNIT < 0)
C         INFILE (A60) Name (including path) of the input spectral file
C
C     Record 4 (for JUNIT < 0)
C         OUTFILE (A60) Name (including path) of the spectral output file
C**********************************************************************

C***********************************************************************
C     LPTSMX is the size of a data block for the FFT.  If the number of
C     data points is LPTSMX or less, the FFT is done in memory.  If it 
C     is larger, then a disk based FFT is performed, and LPTSMX is the 
C     size in words of each record of the file. The in-memory routine is 
C     about 20 percent more efficient than the disk-based routine, for 
C     same number of points.  For efficiency, LPTSMX should be made
C     as large a possible so that the computation can be done in memory.
C     However, on virtual memory machines (e.g. VAX), if LPTSMX is too
C     large, then the computation will be done in VIRTUAL memory.  Since
C     the in-memory routine uses widely scattered points, operating 
C     system will spend a great deal of time swapping data in and out 
C     (thrashing) and the efficiency will be very poor.
C
C     The following problem has been corrected (May, 1993).  The minimum
C     size of an FFT is now the smallest power of 2 greater than the 
C     of data points.  LPTSMX should be set to the larges value consistant
C     with real memory.
C
C*    However, LPTSMX is also the minimum size of an FFT (this is a 
C*    design flaw of this program).  If LPTSMX is much larger than the 
C*    number of points in the spectrum, computational time is wasted.
C*
C*    Therefore, LPTSMX should be set to a value somewhat larger than
C*    the size of the smallest typical spectrum, but no larger than the 
C*    largest value possible without thrashing.
C
C     IBLKSZ is the record length in the OPEN statement, and may be 
C     in bytes or words,  depending on the operating system.  For VAX
C     and CDC, it is in words.  For MicroSoft FORTRAN, it is in bytes.
C***********************************************************************
      PARAMETER (LSIZE = 65536)
C*****Following line for computers where the blocksize is measured
C*****in words, e.g. VAX, CDC/NOS/VE
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX,LPTSM8=LPTSMX/8)

C*****Following line for computers where the blocksize is measured 
C***** in bytes, e.g. MS FORTRAN, SUN, Alliant
      PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*4,LPTSM8=LPTSMX/8)

C*****Following line for computers with 64 bit words where the blocksize is 
C*****measured in bytes, e.g. CRAY
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*8,LPTSM8=LPTSMX/8)

      Parameter (JFNMAX = 13)
C*****Computers with 32 bit words need the Double Precision Statements
C*****Computers with 64 bit words (e.g. Cyber) do not.
C*****Frequency variables start with V
      Implicit Real*8           (V)
      Character*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ

C*****Blank Common carries the spectral data
      COMMON S(2450),R1(2650),XF(251)

C*****HVRFFT carries the module CVSversion numbers
      COMMON /CVRFFT/ HVRFFT

C*****SCNHRD carries the header information for the scanned file
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF 
      DIMENSION FILHDR(2),IFSDID(17)
      EQUIVALENCE (FILHDR(1),XID(1))
      EQUIVALENCE (FSCDID(1),IFSDID(1),IHIRAC),(FSCDID(2),ILBLF4), 
     C (FSCDID(3),IXSCNT),(FSCDID(4 ),IAERSL ),(FSCDID(5),IEMIT), 
     C (FSCDID(6),ISCHDR ),(FSCDID(7 ),IPLOT  ),(FSCDID(8),IPATHL),   
     C (FSCDID(9),JRAD  ),(FSCDID(10),ITEST  ),(FSCDID(11),IMRG),  
     C (FSCDID(12),XSCID),(FSCDID(13),XHWHM  ),(FSCDID(14),IDABS),   
     C (FSCDID(15),IATM ),(FSCDID(16),LAYR1  ),(FSCDID(17),NLAYFS),  
     C (YID(1)    ,HDATE),(YID(2),      HTIME),(YI1,IMULT)        

C*****PANL carries the information from the panel header
      COMMON /PANL/ V1P,V2P,DVP,NP
      DIMENSION PNLHDR(1)
      EQUIVALENCE (PNLHDR(1),V1P)


C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 


      Character*15 HVRFFT
      Character*16 SFNAME,ANAMES(0:JFNMAX)
      Dimension C(0:JFNMAX),CRATIO(0:JFNMAX),CLIMIT(0:JFNMAX)

C*****FUNCT1 and FUNCT2 are used to store the spectrum, the Fourier
C*****transforms, the scanning function and the results
      Dimension FUNCT1(LPTSMX),FUNCT2(LPTSMX)

C*****JFNMAX is number of scanning functions currently defined
C*****ANAMES are their names
C*****SFNAME is the name of the scanning function
      Data ANAMES/'BOXCAR','TRIANGLE','GAUSSIAN','SINC**2','SINC',
     1           'BEER','HAMMING','HANNING','NB - WEAK','NB - MEDIUM',
     2           'NB - STRONG','BRAULT','KAISER-BESSEL','Kurina'/

C*****Ci's are constants = A/ HWHM, where A is the "natural" width
C*****parameter for each function, A = 1/L, where L is the maximum
C*****optical path difference for an equivalent interferometer.
C*****If C < 0, Abs(C) is approximate.
C*****Function -1 is special.  It is a very broad sinc which
C*****is the apodization function corresponding a narrow rectangular
C*****scanning function used to pre-scan the spectrum.
      Data C/     0.0,     2.0,0.849322,2.257609,3.314800,2.100669,
     1       2.195676,     2.0,2.570274,2.367714,2.071759, -3.3, -2.5,
     2       3.314800 /

C*****CRATIO is the critical value of the ratio of the frequency range
C*****to the half width at half maximum of the scanning function.
C*****If this ratio is greater than CRATIO, then the apodization 
C*****function is calculated analytically, otherwise it is calculated 
C*****as the FFT of the scanning function.  The values given here are
C*****educated guesses and may need to be revised. If CRATIO < 0, 
c*****apodization is always calculated as the FFT of the scanning 
C*****function
      Data CRATIO/0., 40., 10., 40., 160., 20., 20. ,20.,
     1   40., 40., 20., 100., 10., -1./

C*****CLIMIT: the limits of the scanned spectrum are expanded by
C*****HWHM*CLIMIT(JFN) to allow for the wrap around effect at V1 and
C*****V2 (effectively, V2 wraps around to V1).  Like CRATIO, these
C*****values are educated guesses (except for the triangle, where the
C*****bound of the scanning function is exactly 2*HWHM.)
      Data CLIMIT/0., 2., 3., 40., 160., 20., 20., 20.,
     1   40., 40., 20., 100., 10., 160/

C*****Note: the values of C, CRATIO, and CLIMIT for JFN=11 (Brault) 
C*****corespond to a value of PARM of about .9, in which case the 
C*****scanning function is very near a sinc.

C*****MRATIO is the minimum ratio of the HWHM of the scanning function
C*****to the width of the boxcar used to prescan the spectrum.  MRATDF
C*****is the default value of this parameter.
      Data MRATDF/12/

C*****MDVHW  is the minimum allowed value of the input dv to the hwhm
C*****If hwhm/dv is less than MDVHW/2, then the spectrum is first 
C*****interpolated onto a grid with dv' = 2*hwhm/MDVHW.  Note: MDVHW 
C*****and MRATIO are coupled: if MDVHW > MRATIO, then the program
C*****will interpolate onto a fine grid but then boxcar back onto a
C*****coarser grid
      Data MDVHW /12/

C*****Assign CVS version number to module fftscn.f

      HVRFFT = '$Revision$' 

      Write(IPR,'(''1FFTSCN: SPECTRAL SMOOTHING IN THE '',
     1            ''FOURIER DOMAIN*****'')')

C*****Read in scan commands until HWHM .le. 0.
  100 Continue
      PARM1 = 0.0
      PARM2 = 0.0
      PARM3 = 0.0

      Read(IRD,10,End=120) HWHM, V1, V2, JEMIT, JFNIN, MRATIN,
     1     DVOUT,IUNIT,IFILST,NIFILS, JUNIT, IVX, NOFIX
   10 Format(3F10.3,3I5,F10.5,4I5,I3,I2)

      Write(IPR,15)  HWHM, V1, V2, JEMIT, JFNIN, MRATIN, DVOUT, IUNIT,
     1     IFILST, NIFILS, JUNIT, IVX, NOFIX
   15 Format(/,'       HWHM        V1        V2   JEMIT   JFNIN',
     1   '  MRATIN', /,1X,3F10.4,3I8,12X,//,
     2   '      DVOUT   IUNIT  IFILST  NIFILS   JUNIT     IVX   NOFIX',
     3    /,(1X,F10.4,6I8))

      If(HWHM .LE. 0) Goto 115

      If(ABS(JFNIN) .ge. 11) Then
C*****    Functions 11 and above require further parameters
          Read(IRD,'(3F10.4)',End=120,Err=122) PARM1,PARM2,PARM3
          WRITE(IPR,16) PARM1,PARM2,PARM3
  16      Format(/,'     PARM1     PARM2     PARM3',/,' ',3F10.4)
      Endif

C*****Check whether JEMIT, JFNIN are within bounds
      If(JEMIT .LT. 0 .OR. JEMIT .GT. 1) THEN
          WRITE(IPR,'('' JEMIT is out of range: STOP'')')
          Goto 125
      Endif

      JFN = ABS(JFNIN)
      If(JFN .GT. JFNMAX) THEN
         Write(IPR,'('' JFN is out of range, max is '',I5,
     1               '':STOP'')') JFNMAX
         Goto 125
      Endif

C*****Calculate the parameter A = 1/L, where L is the maximum optical 
C*****length of an equivalent interferometer (JFNIN > 0) or HWHM = 
C*****half width at half height of the scanning function (JFNIIN < 0)
      SFNAME = ANAMES(JFN)
      If (JFNIN .GT. 0) Then
          If (C(JFN) .GT. 0) Then
              A = C(JFN)*HWHM
              PATHL = 1.0/A
          Else
              Write(IPR,*) ' Error: for this scanning function, you ',
     1            'must specify the interferometer optical ',
     2            'path difference'
              Goto 125
          Endif
      Elseif (JFNIN .LT. 0) Then
          If (C(JFN) .LT. 0) Then
              Write(IPR,*) ' Note: for this scanning function, the ',
     1            'calculated halfwidth is only approximate'
          Endif
          PATHL = HWHM
          A = 1./PATHL
          HWHM = A/ABS(C(JFN))
      Elseif (JFNIN .EQ. 0) Then
          A = 0.
          PATHL = 0.
      Endif

      Write(IPR,17) SFNAME,HWHM,PATHL,A
  17  FORMAT(/,' Scanning function selected is:           ',A,/,
     1         ' Halfwidth at Half Height is:           ',G14.6,/,
     2         ' Length of equivalent interferometer is:',G14.6,/,
     3         ' Parameter A = 1/Length is:             ',G14.6)

C*****Decode boxcaring parameter MRATIN
      If(MRATIN .LT. 0) Then
         NOBOX = 1
      Elseif( MRATIN .EQ. 0) Then
         NOBOX = 0
         MRATIO = MRATDF
      Else
         NOBOX = 0
         MRATIO = MRATIN
      Endif

C*****Check the status of the input and output spectra files, and read
C*****in the filenames if necessary
      Call Ckfile(IFILE,IUNIT,0,IERR)
      If(IERR .NE. 0) Goto 125
      
      Call Ckfile(JFILE,JUNIT,1,IERR)
      If(IERR .NE. 0) Goto 125

      If(IFILST .LE. 0) IFILST = 1
      If(NIFILS .EQ. 0) NIFILS = 1
      If(NIFILS .LT. 0) NIFILS = 99

C*****Loop over NIFILS on IFILE
      Rewind IFILE
      If(IFILST .GT. 1) Then
          Call Skipfl(IFILST-1, IFILE, IEOF)
          If(IEOF .EQ. 0) Then
              WRITE(IPR,'('' EOF encountered skipping files: STOP'')')
              Goto 125
          Endif
      Endif

      Do 105 I=1,NIFILS

C*****Read in file header
      Call Gethdr(IFILE,1,JDATA,IEOFSC)
      If(IEOFSC .EQ. 0) Then
          Write(IPR,'('' EOF encountered reading file: STOP'')')
          Goto 125
      Endif

C*****Check the data in the header against the scan request and
C*****determine the data conversion.
C*****JCONVT   action
C*****     0   single panel, no conversion (scanned trans or rad)
C*****     1   single panel, optical depth to transmittance (mono.)
C*****     2   two panels, get first (mono. radiance)
C*****     3   two panels, set second (mono. transmittance)

      JCONVT = -1
      If(JEMIT .EQ. 0) Then
          If(JDATA .EQ. 0) JCONVT = 1
          If(JDATA .EQ. 1) JCONVT = 0
          If(JDATA .EQ. 2) JCONVT = 3
      Elseif(JEMIT .EQ. 1) Then
          If(JDATA .EQ. 2) JCONVT = 2
          If(JDATA .EQ. 3) JCONVT = 0
      Endif

      If(JCONVT .EQ. -1) Then
          Write(IPR,*) ' Data on file not compatible with scan: STOP'
          Goto 125
      Endif


C*****Check frequency limits and adjust, if necessary.  V1 and V2 are
C*****the requested limits, V1C and V2C are the limits of the input  
C*****data, V1Z and V2Z are the actual limits of the of the scanned 
C*****data and V1S and V2S are the limits of the data as output,  
C*****adjusted from V1 and V2 for several effects:
C*****1.  gridding: V1S and V2S must be on the frequency grid of the 
C*****    output spectrum 
C*****2.  edge effects: V1Z and V2Z are expanded from V1 and V2 
C*****    by HWHM*CLIMIT(JFN) so that edge effects do not contaminate 
C*****    the scanned spectrum
C*****3.  regridding: prescanning with the boxcar resamples onto a new
C*****    frequency grid
C*****4.  interpolation: V2S is adjusted to fall on the interpolated
C*****    grid 

      If(V1 .ge. V2) Then
          Write(IPR,*) ' FFTSCN - input error: Initial V >= final V:',
     1      'STOP'
          Goto 125
      Elseif(V1 .ge. V2C  .or.  V2 .le. V1C) Then
          Write(IPR,*) ' FFTSCN -input error: V1 to V2 is outside the',
     1      ' data range: STOP'
          Goto 125
      Endif

      V1Z = V1-HWHM*CLIMIT(JFN)
      V2Z = V2+HWHM*CLIMIT(JFN)
C*****Adjust V1Z, V2Z to fall on current frequency grid
      V1Z = V1C+DBLE(INT((V1Z-V1C)/DV+1.01))*DV
      V2Z = V1C+DBLE(INT((V2Z-V1C)/DV+1.01))*DV

      If (V1Z .LT. V1C) Then
          V1Z = V1C
          Write(IPR,*) ' First data point used: beginning of spectrum',
     1         ' may suffer from wraparound effects'
      Endif
      If (V2Z .GT. V2C) Then
          V2Z = V2C
          Write(IPR,*) ' Last data point used: end of spectrum may', 
     1         ' suffer from wraparound effects'
      Endif

      Write(IPR,'(/,A,F12.5,A,F12.5)')
     1   ' Frequency limits of internal scanning calculation:',
     2   V1Z,' to',V2Z

      If ((V2Z-V1Z)/HWHM .lt. CLIMIT(JFN)) Then
          Write(IPR,*) ' FFTSCN - error: Insufficient Data for an',
     1         'accurate calculation: STOP'
          Goto 125
      Endif

C*****If the frequency grid is too coarse compared to the HWHM,
C*****then interpolate the spectrum onto a new, finer frequency grid.
C*****The new DV must be <= 2*HWHM/MDVHW, where the parameter MDVHW
C*****is nominally 12 (educated guess.)  
      DVN = 2.*HWHM/MDVHW
      If (DV .gt. DVN) Then
          Write(IPR,*) ' FFTSCN: input grid too coarse. ',
     1        'Interpolating onto a grid with dv = ',DVN

C*****    The interpolated file will be on UNIT = KFILE
          Call Getunt(KFILE)
          Open(KFILE,Status='Scratch',Form='Unformatted')

C*****    Interpolate onto new grid.
          Call Intpdr(IFILE,KFILE,V1Z,V2Z,DVN,IFILST,JEMIT,IERR)
          If (IERR .ne. 0) Then
              Write(IPR,*) 'FFTSCN - error: Error in interpolation:',
     1             ' STOP'
              Goto 125
          Endif

C*****    Need to reposition KFILE to just after the file header
          Rewind KFILE

          Call Gethdr(KFILE,0,JDATA,IEOFSC)
C*****    Following statement correct ??? (have interpolated to the 
C*****    desired parameter??)
          JCONVT = 0
      Else
          KFILE = IFILE
      Endif
 
C*****Decide whether to calculate the apodization function 
C*****analytically (IVX = 1) or as the fft of the scanning
C*****function (IVX = -1). The input value of IVX (if not 0)
C*****overides.
      If( IVX .LT. -1  .OR.  IVX .GT. 1) Then
          Write(IPR,*) ' FFTSCN - error: IVX is out of range, = ',
     1       IVX,': STOP'
          Goto 125
      Endif
      If(IVX .EQ. 0) Then
          If( (V2Z-V1Z)/HWHM .GT. CRATIO(JFN)) Then
              IVX = 1
          Else
              IVX = -1
          Endif
      Endif
      If(JFN .EQ. 0) IVX = 0

      If(IVX .EQ. 1) Then
          Write(IPR,'(/,A,A)') ' Apodization function will be ',
     1         'calculated analytically'
      Elseif (IVX .EQ. -1) Then
          Write(IPR,'(/,A,A)') ' Apodization function will be ',
     1         'calculated as the fft of the scanning function'
      Endif

C*****Parmameters checked, are OK.

C*****Load spectrum from V1Z to V2Z into FUNCT1 (or file on LFILE1
C*****if there are more than LPTSMX points). V1Z and V2Z will be
C*****adjusted to fall on the input frequency grid. Total points
C*****is LTOTAL, the total records (including zeroed records) is LREC
C*****Both LTOTAL and LREC must be a power of 2 (LREC can be 1).
C*****If LREC = 1, then a memory based FFT is used.  If LREC>1, then
C*****a disk based FFT is used.

C*****Get free file unit number for LFILE1 
      Call Getunt(LFILE1)

      Call Loadsp(KFILE,LFILE1,JCONVT,JEMIT,V1Z,V2Z,LREC,LTOTAL,
     1            FUNCT1,IERROR)
      If(IERROR .NE. 0) Then
           Write(IPR,'(2A)') ' LOADSP detected an error: STOP'
           Goto 125
      Endif

      M = 1
      If(JFN .EQ. 0) Then
C*****     Rectangular scanning function (boxcar)
           M = INT(2.0*HWHM/DV)
           If(M .LT. 2) Then
               Write(IPR,*) ' Error - width of rectangular scanning ',
     1               'function is less than 2*DV: STOP'
               Goto 125
           Endif
      Elseif (NOBOX .EQ. 0) Then
C*****     If the scanning function is wide enough compared to the DV,
C*****     it is possible to pre-scan the data with a rectanglar 
C*****     function (boxcar) before the convolution.  This procedure
C*****     reduces the number of points in the spectrum, saving
C*****     computer time and reducing storage space.  The effect of 
C*****     the rectangular smoothing can be eliminated by deconvolving
C*****     after smoothing with the regular scanning function.  This 
C*****     procedure is disabled if MRATIN < 0.

C*****     The value of M, the number of points to average, is based 
C*****     on the criterion that the full width of the boxcar (M*DV) be
C*****     no more than 2*HWHM/MRATIO.
           M = 2.0*HWHM/(MRATIO*DV)
           If(M .lt. 2) Write(IPR,*) ' Boxcaring not possible--M lt 2'
      Endif
      If((M .GE. 2) .AND. (NOBOX .EQ. 0)) Then
C*****    Adjust V1, V2, and DV
          V1Z = V1Z+DV*(M-1)/2.0
          KTOTAL = LTOTAL/M
          DV  = M*DV
          V2Z = V1Z+DV*(KTOTAL-1)
          Write(IPR,25) ' Prescanning the spectrum with the BOXCAR ',
     1       'instrument function',
     2       ' Ratio of boxcar to HWHM (MRATIO) = ',MRATIO, 
     3       ' Shrink ratio (M) = ',M,
     4       ' New DV = ',DV
   25     Format(/,A,A,/,A,I4,/,A,I4,/,A,F12.5)

C*****    Perform the averaging
          Call Boxcar(M,LFILE1,FUNCT1,LREC,LTOTAL,JEMIT)
      Endif

C*****If JFN = 0, then the smoothing (boxcar) is complete.
      If(JFN .EQ. 0) Goto 95

C*****If the FFT can be performed in memory (LREC = 1), then find the
C*****smallest sized FFT possible = LPTFFT = smallest power of 2 > LTOTAL 
C*****and < LPTSMX
      If(LREC .eq. 1) Then
          POWER =  LOG( REAL(LTOTAL))/ LOG(2.0)
          If(POWER-INT(POWER) .EQ. 0) Then
              LPTFFT = 2**INT(POWER)
          Else
              LPTFFT = 2**(INT(POWER)+1)
          Endif
         
          Write(IPR,*) ' In-memory FFT: points, FFT points = ',
     1       LTOTAL,', ',LPTFFT
      Else
          LPTFFT = -1
          Write(IPR,*) ' Disk-based FFT: points,  FFT records = ',
     1       LTOTAL,', ',LREC
      Endif

C*****Take FFT of spectrum, store in FUNCT1 if LTOTAL > LPTSMX, else
C*****put on file on unit LFILE1
      Call Fourtr(LREC,LPTFFT,LFILE1,FUNCT1,1)

C*****Calculate apodization and store in FUNCT2 or on LFILE2
      Call Getunt(LFILE2)
      Call Scntrn(JFN,A,IVX,DV,LREC,LPTFFT,LFILE2,FUNCT2,
     1            PARM1,PARM2,PARM3)

C*****Multiply FUNCT1 and FUNCT2 and store the result in FUNCT1 or
C*****on LFILE1
      Call Multrn(LREC,LFILE2,LFILE1,FUNCT2,FUNCT1)

C*****If the spectrum has been pre-scanned with a rectangle, remove the
C*****effects of the rectangle by deconvolving, i.e. divide the 
C*****transform of the scanned spectrum by the transform of the 
C*****rectangle.  This later function is a very broad sinc.
      If(M .GE. 2 .AND. NOBOX.EQ.0 .AND. NOFIX .EQ. 0) Then

C*****    Calculate transform of the boxcar which is a rectangle of
C*****    HWHM = DV/2, where DV is the new DV. The parameter A =< HWHM
          Write(IPR,*) ' Deconvolving'
          Call Scntrn(-1,DV/2.0,1,DV,LREC,LPTFFT,LFILE2,FUNCT2,
     *         0.0,0.0,0.0)
          Call Multrn(LREC,LFILE2,LFILE1,FUNCT2,FUNCT1)
      Endif

C*****Transform FUNCT1 back to spectral domain, store in FUNCT1 or on 
C*****LFILE1
      Call Fourtr(LREC,LPTFFT,LFILE1,FUNCT1,-1)

   95 Continue

C*****Adjust V1S and V2S to fit the current spectral grid 
C*****Actually, let V1S be one DV less than the largest VZ <= V1
C*****and let V2S be one DV larger than the smallest VZ >= V2
C*****This procedure gives two points beyond V1 and V2, as required for
C*****4 point interpolation 
      N1 = INT((V1-V1Z)/DV)
      V1S = V1Z+(N1-1)*DV
      If (V1S .LT. V1Z) Then
          N1 = 1
          V1S = V1Z
      Endif

      N2 = INT((V2-V1Z)/DV+3.01)
      V2S = V1Z+(N2-1)*DV
      If (V2S .GT. V2Z) Then 
          N2 = (V2Z-V1Z)/DV+1.01
          V2S = V1Z+(N2-1)*DV
      Endif

      NTOTAL = N2-N1+1

      Write(IPR,'(/,A,/A,F12.5,/,A,F12.5,/,A,F12.5,/,A,I7)')
     1    ' Adjusted Limits of Scanned Spectrum:',
     2    ' V1 = ',V1S,' V2 = ',V2S,' DV = ',DV,' N  =',NTOTAL

C*****Write FUNCT1 = scanned spectrum out to a file on JUNIT in LBLRTM
C*****format.  First, modify the header.
      V1C = V1S
      V2C = V2S
      JVAR = 0
      XSCID = JVAR + 10*(JFN + 10*(JEMIT))
      XHWHM = HWHM

      ISCHDR = ISCHDR+1

C*****If DVOUT > 0, then write the spectral file to a temporary file, 
C*****and interpolate it to JFILE
      If (DVOUT .gt. 0) Then
          Call Getunt(NFILE)
          Open(NFILE,Status='Scratch',Form='Unformatted')
          Rewind IFILE
      Else
          NFILE = JFILE
      Endif

      Call BUFOUT(NFILE,FILHDR(1),NFHDRF)
      Call Wrtspc(N1,NTOTAL,LREC,LFILE1,FUNCT1,NFILE)

C*****Interpolate the spectrum onto the desired grid
      If (DVOUT .gt. 0) Then

          Write(IPR,'(/,A,2F12.4,F12.6)') ' Interpolating the '//
     1        'spectrum to the grid definded by:', V1,V2,DVOUT

          Call Intpdr(NFILE,JFILE,V1,V2,DVOUT,1,JEMIT,IERR)
          If (IERR .ne. 0) Then
              Write(IPR,*) 'FFTSCN - error: Error in interpolation:',
     1             ' STOP'
              Goto 125
          Endif
          Close(NFILE)
      Endif

C*****End of convolution for this scan request

      Close(LFILE1)
      If (LFILE2.NE.0) Close(LFILE2)

C*****End of file loop for this scan request
  105 Continue

C*****Read next scan request
      Goto 100

  115 Continue
      Write(IPR,*) ' End of scan function requests'
      Return

  120 Continue
      Write(IPR,*) ' FFTSCN: End of file detected on unit ',IRD
      Stop 'Stopped in FFTSCN'

  122 Continue
   
      Write(IPR,*) ' FFTSCN: Error in reading scan fn parameters'
      Stop 'Stopped in FFTSCN'

  125 Continue
      Write(IPR,*) ' Stopped in FFTSCN due to error'
      Stop ' Stopped in FFTSCN: See TAPE6'
      End
      Subroutine Boxcar(M,LFILE,FUNCT,LREC,LTOTAL,JEMIT)
C************************************************************************
C     This subroutine performs the "BOXCAR" smoothing of the spectrum.
C     This smoothing simply averages M adjacent points together, reducing
C     the number of output points by a factor of M. (It is not a running
C     average where the number of input and output points are the same.)
C     The spectrum is stored either on the file LFILE, if the number of
C     records LREC > 1, or in the array FUNCT if LREC = 1.  LTOTAL is
C     the total number of output points after smoothing.
C     The smoothed spectrum is stored in the array FUNCT if LTOTAL <
C     LPTSMX, or is written back to the file on unit LFILE, in which
C     LREC becomes the new number of records. (The old records for L >
C     LREC will still exist but will be ignored.)
C     LREC must still be a power of 2, and zeroed records will be added
C     as necessary. JEMIT flags transmittance (0) or radiance (1)
C************************************************************************

C***********************************************************************
C     LPTSMX is the size of a data block for the FFT.  If the number of
C     data points is LPTSMX or less, the FFT is done in memory.  If it 
C     is larger, then a disk based FFT is performed, and LPTSMX is the 
C     size in words of each record of the file. The in-memory routine is 
C     about 20 percent more efficient than the disk-based routine, for 
C     same number of points.  For efficiency, LPTSMX should be made
C     as large a possible so that the computation can be done in memory.
C     However, on virtual memory machines (e.g. VAX), if LPTSMX is too
C     large, then the computation will be done in VIRTUAL memory.  Since
C     the in-memory routine uses widely scattered points, operating 
C     system will spend a great deal of time swapping data in and out 
C     (thrashing) and the efficiency will be very poor.
C
C     The following problem has been corrected (May, 1993).  The minimum
C     size of an FFT is now the smallest power of 2 greater than the 
C     of data points.  LPTSMX should be set to the larges value consistant
C     with real memory.
C
C*    However, LPTSMX is also the minimum size of an FFT (this is a 
C*    design flaw of this program).  If LPTSMX is much larger than the 
C*    number of points in the spectrum, computational time is wasted.
C*
C*    Therefore, LPTSMX should be set to a value somewhat larger than
C*    the size of the smallest typical spectrum, but no larger than the 
C*    largest value possible without thrashing.
C
C     IBLKSZ is the record length in the OPEN statement, and may be 
C     in bytes or words,  depending on the operating system.  For VAX
C     and CDC, it is in words.  For MicroSoft FORTRAN, it is in bytes.
C***********************************************************************
      PARAMETER (LSIZE = 65536)
C*****Following line for computers where the blocksize is measured
C*****in words, e.g. VAX, CDC/NOS/VE
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX,LPTSM8=LPTSMX/8)

C*****Following line for computers where the blocksize is measured 
C***** in bytes, e.g. MS FORTRAN, SUN, Alliant
      PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*4,LPTSM8=LPTSMX/8)

C*****Following line for computers with 64 bit words where the blocksize is 
C*****measured in bytes, e.g. CRAY
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*8,LPTSM8=LPTSMX/8)
C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 
      Dimension FUNCT(LPTSMX),BOX(LPTSMX)

C*****Calculate the number of points and records after smoothing
C*****KRDATA is the number of output points (including partial 
C*****records) which include real data (not just FILL).
      KTOTAL = LTOTAL/M
      ARDATA =  REAL(KTOTAL)/ REAL(LPTSMX)
      If((ARDATA-INT(ARDATA)) .EQ. 0.0) Then
          KRDATA = INT(ARDATA)
      Else
          KRDATA = INT(ARDATA)+1
      Endif
C*****Find KREC = smallest power of 2 .GE. KRDATA
      POWER =  LOG( REAL(KRDATA))/ LOG(2.0)
      If((POWER-INT(POWER)) .EQ. 0) Then
          KREC = 2**INT(POWER)
      Else
          KREC = 2**INT(POWER+1)
      Endif

C*****Get the first record from LFILE, if necessary
      If(LREC .GT. 1) Then
          IREC = 1
          Read(LFILE,rec=IREC,err=900) (FUNCT(I),I=1,LPTSMX)
      Endif

C*****I is index for input points, J is index for output points within
C*****a record. IREC is index for input records, JREC is for
C*****output records.
      I = 0
      J = 0
      JREC = 0

C*****Loop over output points
      Do 210 K=1,KTOTAL

C*****    Loop over the M elements for each output point
          SUM = 0.0
          DO 200 MM=1,M
              If(I .EQ. LPTSMX) Then
C*****            Get another record from LFILE
                  IREC = IREC+1
                  Read(LFILE,rec=IREC,err=900) (FUNCT(II),II=1,LPTSMX)
                  I = 0
              Endif
              I = I+1
              SUM = SUM + FUNCT(I)
  200     Continue

          J = J+1
          BOX(J) = SUM/ REAL(M)

          If(J .EQ. LPTSMX) Then
C*****        Write out a record
              JREC = JREC+1
              Write(LFILE,rec=JREC,err=910) (BOX(JJ),JJ=1,LPTSMX)
              J = 0
          Endif

  210 Continue

C*****If necessary, fill out the end of the last record with
C*****0's (radiance) or 1's (transmittance) 
      If(J .GT. 0) Then
          If(JEMIT .EQ. 0) THEN
              FILL = 1.0
          Else
              FILL = 0.0
          Endif
          Do 220 JJ=J+1,LPTSMX
              BOX(JJ) = FILL
  220     Continue
      Endif

C*****If only one record, copy BOX to FUNCT
      If (KREC .EQ. 1) Then
          Do 300 JJ=1,LPTSMX
              FUNCT(JJ) = BOX(JJ)
  300     Continue
C*****Else if there is an unwritten partial record, write it
      Elseif (J .GT. 0) Then
          JREC = JREC+1
          Write(LFILE,rec=JREC,err=910) (BOX(J),J=1,LPTSMX)
      Else
      Endif

C*****Write empty records if necessary.
      If(KREC .GT. KRDATA) Then
          Do 310 J=1,LPTSMX
              BOX(J) = FILL
  310     Continue
          Do 320 K=KRDATA+1,KREC
              JREC = JREC+1
              Write(LFILE,rec=JREC,err=910) (BOX(J),J=1,LPTSMX)
  320     Continue
      Endif

C*****At this point JREC should equal KREC
      If(KREC .GT. 1  .AND.  JREC .NE. KREC) Then
          Write(IPR,'(/,2A,I4,A,I4)') ' BOXCAR-ERROR: JREC .NE.',
     1       ' KREC,  JREC = ',JREC,' AND KREC = ',KREC
          Stop 'Stopped in BOXCAR'
      Endif

C*****Reset the total number of records
      LREC = KREC
      LTOTAL = KTOTAL

      Return

  900 Continue
      Write(IPR,*) ' BOXCAR: error reading data from unit ',LFILE,
     1    '   record = ',IREC
      Stop 'Stopped in BOXCAR'

  910 Continue
      Write(IPR,*) ' BOXCAR: error writing data to unit ',LFILE,
     1    '   record = ',KREC
      Stop 'Stopped in BOXCAR'

      End
      Subroutine Ckfile(IFILE,IUNIT,IDIR,IERR)
C***********************************************************************
C     This subroutine checks the status of a file and opens one if 
C     if necessary.  
C     IFILE is the default unit number. 
C     The user specifies a file on IUNIT.  If IUNIT = 0, then the file
C     on IFILE is used.  If IUNIT < 0, then a file name (20 characters
C     max) is read in, and IFILE = -IUNIT
C     IDIR specifies the direction: 0 = Input, 1 = Output.
C     IERR is an error flag: 0 = no error, 1 = error
C     
C     Input:
C     If IUNIT GE 0, Then
C        If IUNIT = 0 Then IFILE = IFILE Else IFILE = IUNIT
C        If a file is not open on IFILE, Then
C            Look for a file of the name TAPExx, where xx = IFILE
C            If TAPExx does not exist, then ERR = 1, Return
C            Open TAPExx on IFILE
C        Return
C     If IUNIT LT 0, Then
C        IFILE = -IUNIT
C        Read in FILENAME
C        If FILENAME does not exist, Then IERR = 1, Return
C        Open FILENAME on IFILE
C        Return
C
C     Output:
C     If IUNIT GE 0, Then
C         If IUNIT = 0 Then IFILE = IFILE Else IFILE = IUNIT
C         IF a file is not open on IFILE, Then
C             Set FILENAME to 'TAPExx' where xx is IFILE
C             If FILENAME exists Then
C                overwrite FILENAME
C             Else
C                Open FILENAME
C         Return
C
C     If IUNIT < 0, Then
C         IFILE = -IUNIT
C         Read in FILENAME
C         IF FILENAME exists, Then IERR = 1, Return
C         Open FILENAME on IFILE
C     Return
C***********************************************************************

C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 
      
      LOGICAL OP,EX
      Character FILNAM*60,CTAPE*4

C*****CTAPEdefines the default prefix for LBLRTM FILENAMEs, e.g. TAPE12
      Data FILNAM/' '/,CTAPE/'TAPE'/
      
      IERR = 0
      If(IDIR .EQ. 0) Then
C*****    File for input
          If(IUNIT .GE. 0) Then
              If(IUNIT .GT. 0) IFILE = IUNIT
              Inquire(UNIT=IFILE,OPENED=OP)
              If(.NOT.OP) Then
                  Write(FILNAM,'(A,I2.2)') CTAPE,IFILE
                  Write(IPR,'(/,A,A)') ' Input file name is: ',FILNAM
                  Inquire(FILE=FILNAM,EXIST=EX)
                  If (EX) Then
                      Open(IFILE,FILE=FILNAM,STATUS='OLD',
     1                    FORM='UNFORMATTED')
                  Else
                      Write(IPR,'(/,A,A)') ' Error: Input file does ',
     1                    'not exist'
                      IERR=1
                      Return
                  Endif
              Endif
          Else
              Read(IRD,'(A)') FILNAM
              Write(IPR,'(A,A)') ' Input  file name is: ',FILNAM
              Inquire(FILE=FILNAM,EXIST=EX)
              If(.NOT. EX) Then
                  Write(IPR,'(3A)') ' Error: input file does not ',
     1                'exist'
                  IERR = 1
                  Return
              Endif
C*****        Get a free file unit number
              Call Getunt(IFILE)
              Open(UNIT=IFILE,FILE=FILNAM,STATUS='OLD',
     1            FORM='UNFORMATTED')
          Endif

      Else
C*****    File for output
          If(IUNIT .GE. 0) Then
              If(IUNIT .GT. 0) IFILE = IUNIT
C*****        Use this file even if it already exists
              Inquire(UNIT=IFILE,OPENED=OP)
              If(OP) Then
                  Rewind IFILE
              Else
                  Write(FILNAM,'(A,I2.2)') CTAPE,IFILE
                  Write(IPR,'(/,A,A)') ' Output file name is: ',FILNAM
                  Open(IFILE,FILE=FILNAM,STATUS='UNKNOWN',
     1                FORM='UNFORMATTED')
              Endif
          Else
              Read(IRD,'(A)') FILNAM
              Write(IPR,'(2A)') ' Output file name is: ',FILNAM
C*****        Get a free file unit number
              Call Getunt(IFILE)
              OPEN(UNIT=IFILE,FILE=FILNAM,STATUS='UNKNOWN',
     1            FORM='UNFORMATTED')
          Endif
      Endif

      Return
      End
      Subroutine Fourtr(LREC,LPTFFT,LFILE,SPECT,IDIR)
C***********************************************************************
C     Fourtr computes the Fourier transform of a data set.
C     If LREC = 1,then the input data is in the array SPECT and the
C     result is returned in SPECT. LPTFFT is the size of the FFT, and
C     must be less or equal to than LPTSMX.
C     If LREC > 1, then the input data is on file LFILE in blocks of 
C     LPTSMX and the result is written to LFILE.  IDIR is 0 for going 
C     from the time domain to the spectral domain, 1 for the other 
C     direction.
C***********************************************************************

C***********************************************************************
C     LPTSMX is the size of a data block for the FFT.  If the number of
C     data points is LPTSMX or less, the FFT is done in memory.  If it 
C     is larger, then a disk based FFT is performed, and LPTSMX is the 
C     size in words of each record of the file. The in-memory routine is 
C     about 20 percent more efficient than the disk-based routine, for 
C     same number of points.  For efficiency, LPTSMX should be made
C     as large a possible so that the computation can be done in memory.
C     However, on virtual memory machines (e.g. VAX), if LPTSMX is too
C     large, then the computation will be done in VIRTUAL memory.  Since
C     the in-memory routine uses widely scattered points, operating 
C     system will spend a great deal of time swapping data in and out 
C     (thrashing) and the efficiency will be very poor.
C
C     The following problem has been corrected (May, 1993).  The minimum
C     size of an FFT is now the smallest power of 2 greater than the 
C     of data points.  LPTSMX should be set to the larges value consistant
C     with real memory.
C
C*    However, LPTSMX is also the minimum size of an FFT (this is a 
C*    design flaw of this program).  If LPTSMX is much larger than the 
C*    number of points in the spectrum, computational time is wasted.
C*
C*    Therefore, LPTSMX should be set to a value somewhat larger than
C*    the size of the smallest typical spectrum, but no larger than the 
C*    largest value possible without thrashing.
C
C     IBLKSZ is the record length in the OPEN statement, and may be 
C     in bytes or words,  depending on the operating system.  For VAX
C     and CDC, it is in words.  For MicroSoft FORTRAN, it is in bytes.
C***********************************************************************
      PARAMETER (LSIZE = 65536)
C*****Following line for computers where the blocksize is measured
C*****in words, e.g. VAX, CDC/NOS/VE
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX,LPTSM8=LPTSMX/8)

C*****Following line for computers where the blocksize is measured 
C***** in bytes, e.g. MS FORTRAN, SUN, Alliant
      PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*4,LPTSM8=LPTSMX/8)

C*****Following line for computers with 64 bit words where the blocksize is 
C*****measured in bytes, e.g. CRAY
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*8,LPTSM8=LPTSMX/8)
C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 

C*****For LREC = 1, SPECT contains the function to be transformed.
      DIMENSION SPECT(*)

C*****A, B, and IWORK are work spaces used for REALBG. If REALBG
C*****is used, then SPECT is also used as work space.
      DIMENSION A(LPTSMX),B(LPTSMX),IWORK(LPTSM8)

C*****SPECT and IWORK are never needed at the same time
C     EQUIVALENCE (SPECT,IWORK)

      If(LREC .EQ. 1) Then
C*****    Use FFT from Numerical Recipies.  REALFT(DATA,N,IDIR)  expects
C*****    DATA to be dimensioned 2*N.
          If(LPTFFT .gt. LPTSMX) Goto 900              
          Call REALFT(SPECT,LPTFFT/2,IDIR)
      Else
C*****    Use file based FFT from M. Esplin

          Call REALBG(LFILE,IDIR,LREC,LPTSMX,A,B,SPECT,IWORK)

      Endif

      Return

  900 Continue
      Write(IPR,*) ' Fourtr: Error-LPRFFT > LPRSMX'
      Write(IPR,*) ' LPTFFT = ',LPTFFT, ' LPTSMX = ',LPTSMX
      Stop

      End
      Subroutine Gethdr (IFILE,IPRNT,JDATA,IEOFSC)
C************************************************************************
C     This subroutine reads the first record an LBLRTM file (the file
C     header) from unit IFILE.  It determines what data is on the file
C     according to the following scheme:
C     Determine what quantity is on the file (i.e., Optical depth,
C     Transmittance, Transmittance and Radiance, Radiance)
C     JDATA    IFILE contains
C         0    optical depth  (monochromatic)
C         1    transmittance  (previously scanned)
C         2    radiance and transmittance  (monochromatic)
C         3    radiance  (previously scanned)
C
C     If IPRNT = 1, then the contents of the file header are printed.
C************************************************************************

C*****Computers with 32 bit words need the Double Precision Statements
C*****Computers with 64 bit words (e.g. Cyber) do not.
C*****Frequency variables start with V
      Implicit Real*8           (V)
      Character*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ

C*****Blank Common carries the spectral data
      COMMON S(2450),R1(2650),XF(251)

C*****SCNHRD carries the header information for the scanned file
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF 
      DIMENSION FILHDR(2),IFSDID(17)
      EQUIVALENCE (FILHDR(1),XID(1))
      EQUIVALENCE (FSCDID(1),IFSDID(1),IHIRAC),(FSCDID(2),ILBLF4), 
     C (FSCDID(3),IXSCNT),(FSCDID(4 ),IAERSL ),(FSCDID(5),IEMIT), 
     C (FSCDID(6),ISCHDR ),(FSCDID(7 ),IPLOT  ),(FSCDID(8),IPATHL),   
     C (FSCDID(9),JRAD  ),(FSCDID(10),ITEST  ),(FSCDID(11),IMRG),  
     C (FSCDID(12),XSCID),(FSCDID(13),XHWHM  ),(FSCDID(14),IDABS),   
     C (FSCDID(15),IATM ),(FSCDID(16),LAYR1  ),(FSCDID(17),NLAYFS),  
     C (YID(1)    ,HDATE),(YID(2),      HTIME),(YI1,IMULT)        

C*****PANL carries the information from the panel header
      COMMON /PANL/ V1P,V2P,DVP,NP
      DIMENSION PNLHDR(4)
      EQUIVALENCE (PNLHDR(1),V1P)


C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 

C*****Read in the file header
      Call BUFIN(IFILE,IEOFSC,FILHDR(1),NFHDRF)

      If(IEOFSC .le. 0) Then
          Write(IPR,*) ' GETHDR: EOF on File ',IFILE,
     1       '  IEOFSC = ',IEOFSC 
          Return
      Endif

C*****Determine what data is on the file
      JDATA = -1
      If(XSCID .le. 0) Then
          ISCAN = 0
      Else
          ISCAN = ISCHDR
      Endif
      If(IEMIT .EQ. 0) Then
          If(ISCAN .LE. 0) Then
              JDATA = 0
          Else
              JDATA = 1
          Endif
      Else
          If(ISCAN .LE. 0) Then
              JDATA = 2
          Else
              JDATA = 3
          Endif
      Endif

      If(JDATA .EQ. -1) Then
          Write(IPR,10) ' GETHDR: error, JDATA = -1, type of data '//
     c       'cannot be determined'
      Endif

      If(IPRNT .EQ. 1 .OR. JDATA .EQ. -1) Then

      Write(IPR,10) IFILE
   10 FORMAT(/,' File Header for file on unit',I3)

      WRITE(IPR,20) XID,(YID(M),M=1,2)
   20 FORMAT(' ',10A8,2X,2(1X,A8,1X))

      WRITE(IPR,30) LAYR1,NLAYFS
   30 FORMAT(/,' Initial Layer = ',I4,',  Final Layer =',I4)

      WRITE(IPR,40) SECANT,PAVE,TAVE,DV,V1C,V2C,JDATA
   40 FORMAT(/,' Secant     =',F12.5,/ ' Press (mb) =',F12.5 /,
     1         ' Temp (K)   =',F9.2 ,/,' DV (cm-1)  =',F15.8,/,
     2         ' V1 (cm-1)  =',F13.6,/ ' V2 (cm-1)  =',F13.6,//,
     3         ' JDATA      =',I6)

      WRITE(IPR,50) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)
   50 FORMAT(/,' Column Density (molecules/cm**2)',
     1   /, '      Wbroad = ',1PE10.3,
     2   /,(A12,       ' = ',1PE10.3))

      Endif

      Return
      End
      Subroutine Getpnl(IFILE,JCONVT,IEOFSC)
C***********************************************************************
C     This subroutine gets one panel of data from the file on IFILE in
C     LBLRTM format. A panel consists of a header plus one or two data
C     records.  The header contains the  the initial and final
C     wavenumbers, wavenumber increment, and number of points in the
C     record. The data contained in the record(s) depends upon the 
C     initial calculation and whether the file has been scanned
C     previously.
C         Initial Calculation   record
C            transmittance      one record, optical depth
C            radiance           two records, radiance then transmittance
C     If the file has been scanned previously, then there is only one
C     record, either transmittance or radiance.
C
C     JCONVT determins which data is to be extracted and whether it is
C     to be converted.
C
C     JCONVT   action
C          0   single record, no conversion (scanned trans or rad)
C          1   single record  optical depth to transmittance (mono.)
C          2   two records, get first (mono. radiance)
C          3   two records, get second (mono. transmittance)
C
C     IEOFSC equals 0 if a hardware end of file is encountered, equals
C      -99 for the end of an LBLRTM file, otherwise it is 1.
C
C     V1P, V2P, DVP, NP are the initial and final wavenumbers, the
C     wavenumber increment and the number of points in the record(s).
C     They are carried in COMMON /PANL/
C
C     The spectral data is returned in S carried in blank common.
C
C***********************************************************************

C*****Computers with 32 bit words need the Double Precision Statements
C*****Computers with 64 bit words (e.g. Cyber) do not.
C*****Frequency variables start with V
      Implicit Real*8           (V)
      Character*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ

C*****Blank Common carries the spectral data
      COMMON S(2450),R1(2650),XF(251)

C*****SCNHRD carries the header information for the scanned file
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF 
      DIMENSION FILHDR(2),IFSDID(17)
      EQUIVALENCE (FILHDR(1),XID(1))
      EQUIVALENCE (FSCDID(1),IFSDID(1),IHIRAC),(FSCDID(2),ILBLF4), 
     C (FSCDID(3),IXSCNT),(FSCDID(4 ),IAERSL ),(FSCDID(5),IEMIT), 
     C (FSCDID(6),ISCHDR ),(FSCDID(7 ),IPLOT  ),(FSCDID(8),IPATHL),   
     C (FSCDID(9),JRAD  ),(FSCDID(10),ITEST  ),(FSCDID(11),IMRG),  
     C (FSCDID(12),XSCID),(FSCDID(13),XHWHM  ),(FSCDID(14),IDABS),   
     C (FSCDID(15),IATM ),(FSCDID(16),LAYR1  ),(FSCDID(17),NLAYFS),  
     C (YID(1)    ,HDATE),(YID(2),      HTIME),(YI1,IMULT)        

C*****PANL carries the information from the panel header
      COMMON /PANL/ V1P,V2P,DVP,NP
      DIMENSION PNLHDR(4)
      EQUIVALENCE (PNLHDR(1),V1P)


C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 

      DIMENSION DUMMY(1)

C*****Read the panel header
      Call BUFIN(IFILE,IEOFSC,PNLHDR(1),NPHDRF)
      If(IEOFSC .LE. 0) Go to 200

C*****Read in the spectral data and convert if necessary

      Call BUFIN(IFILE,IEOFSC,S(1),NP)
      If(IEOFSC .LE. 0) Go to 210
      If(JCONVT .EQ. 3) Call BUFIN(IFILE,IEOFSC,S(1),NP)
      If(IEOFSC .LE. 0) Go to 210
      If(JCONVT .EQ. 1) Then
          Do 100 N=1,NP
              If(S(N) .LT. 0) Then
                 S(N) = 1.0
              Elseif(S(N) .GE. ARGMIN) Then
                 S(N) = EXPMIN
              Else
                 S(N) = EXP(-S(N))
              Endif
  100     Continue
      Endif
      If(JCONVT .EQ. 2) Then
          Call BUFIN(IFILE,IEOFSC,DUMMY(1),1)
          If(IEOFSC .LE. 0) Go to 210

      Endif

  200 Continue

      Return

  210 Continue
      Write(IPR,10) 'GETPNL: error, end of file encountered after'//
     c   ' panel header'
   10 Format(1x,A)
      Stop 'Stopped in GETPNL'
      End
      Subroutine Getunt(IFILE)
C**********************************************************************
C     This subroutine gets the first free file unit number > 60
C     If there are none, stop with error message
C**********************************************************************

C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 

      LOGICAL OP

      Do 100 I=61,99
          Inquire(UNIT=i,OPENED=OP)
          If(.NOT. OP) Then
              IFILE = I
              Return
          Endif
  100 Continue
      
      Write(IPR,'(2A)') ' Error: no free file unit numbers 61 to 99, ',
     1   'Stop'
      Stop
      End

              
      Subroutine Loadsp(IFILE,LFILE,JCONVT,JEMIT,V1,V2,LREC,LTOTAL,
     1                  SPECT,IERROR)
C***********************************************************************
C     This subroutine loads the spectral data between V1 and V2 from
C     the file on IFILE into the array SPECT, if the total number of
C     points LTOTAL < LPTSMX, (ie 1 record) or else writes it to the
C     direct access file on unit LFILE in records of LPTSMX.  LREC is 
C     the number of records.  Both LPTSMX and LREC must be powers of 2
C     (LREC may be 1).  Partial records will be completed with zero's
C     (radiance) or 1's (transmittance) and records of all zero's or
C     1's  will be added as needed.  JCONVT indicates how the data 
C     is to be extracted from the LBLRTM file. JEMIT flags transmittance
C     (0) or radiance(1).
C***********************************************************************

C*****Computers with 32 bit words need the Double Precision Statements
C*****Computers with 64 bit words (e.g. Cyber) do not.
C*****Frequency variables start with V
      Implicit Real*8           (V)
      Character*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ

C*****Blank Common carries the spectral data
      COMMON S(2450),R1(2650),XF(251)

C*****SCNHRD carries the header information for the scanned file
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF 
      DIMENSION FILHDR(2),IFSDID(17)
      EQUIVALENCE (FILHDR(1),XID(1))
      EQUIVALENCE (FSCDID(1),IFSDID(1),IHIRAC),(FSCDID(2),ILBLF4), 
     C (FSCDID(3),IXSCNT),(FSCDID(4 ),IAERSL ),(FSCDID(5),IEMIT), 
     C (FSCDID(6),ISCHDR ),(FSCDID(7 ),IPLOT  ),(FSCDID(8),IPATHL),   
     C (FSCDID(9),JRAD  ),(FSCDID(10),ITEST  ),(FSCDID(11),IMRG),  
     C (FSCDID(12),XSCID),(FSCDID(13),XHWHM  ),(FSCDID(14),IDABS),   
     C (FSCDID(15),IATM ),(FSCDID(16),LAYR1  ),(FSCDID(17),NLAYFS),  
     C (YID(1)    ,HDATE),(YID(2),      HTIME),(YI1,IMULT)        

C*****PANL carries the information from the panel header
      COMMON /PANL/ V1P,V2P,DVP,NP
      DIMENSION PNLHDR(4)
      EQUIVALENCE (PNLHDR(1),V1P)


C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 
C***********************************************************************
C     LPTSMX is the size of a data block for the FFT.  If the number of
C     data points is LPTSMX or less, the FFT is done in memory.  If it 
C     is larger, then a disk based FFT is performed, and LPTSMX is the 
C     size in words of each record of the file. The in-memory routine is 
C     about 20 percent more efficient than the disk-based routine, for 
C     same number of points.  For efficiency, LPTSMX should be made
C     as large a possible so that the computation can be done in memory.
C     However, on virtual memory machines (e.g. VAX), if LPTSMX is too
C     large, then the computation will be done in VIRTUAL memory.  Since
C     the in-memory routine uses widely scattered points, operating 
C     system will spend a great deal of time swapping data in and out 
C     (thrashing) and the efficiency will be very poor.
C
C     The following problem has been corrected (May, 1993).  The minimum
C     size of an FFT is now the smallest power of 2 greater than the 
C     of data points.  LPTSMX should be set to the larges value consistant
C     with real memory.
C
C*    However, LPTSMX is also the minimum size of an FFT (this is a 
C*    design flaw of this program).  If LPTSMX is much larger than the 
C*    number of points in the spectrum, computational time is wasted.
C*
C*    Therefore, LPTSMX should be set to a value somewhat larger than
C*    the size of the smallest typical spectrum, but no larger than the 
C*    largest value possible without thrashing.
C
C     IBLKSZ is the record length in the OPEN statement, and may be 
C     in bytes or words,  depending on the operating system.  For VAX
C     and CDC, it is in words.  For MicroSoft FORTRAN, it is in bytes.
C***********************************************************************
      PARAMETER (LSIZE = 65536)
C*****Following line for computers where the blocksize is measured
C*****in words, e.g. VAX, CDC/NOS/VE
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX,LPTSM8=LPTSMX/8)

C*****Following line for computers where the blocksize is measured 
C***** in bytes, e.g. MS FORTRAN, SUN, Alliant
      PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*4,LPTSM8=LPTSMX/8)

C*****Following line for computers with 64 bit words where the blocksize is 
C*****measured in bytes, e.g. CRAY
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*8,LPTSM8=LPTSMX/8)

      DIMENSION SPECT(LPTSMX)
      IERROR = 0

C*****Skip over panels until the V1 is reached
  100 Continue
      Call GETPNL(IFILE,JCONVT,IEOFSC)
C*****Check for EOF
      If(IEOFSC .lt. 0) Then
          Write(IPR,10) '  LOADSP: end of file encountered on file'//
     c       'LFILE = ',LFILE,', skipping to next scan request'
   10     Format(A,I4,A)
          IERROR = 1
          Return
      Endif

C*****If V1P .GT. V1, then error
      If(V1P .GT. V1) Then
          Write(IPR,20) '  LOADSP: error, first data freq = ',V1P,
     c       ' is greater than requested V1 = ',V1
   20     Format(A,F12.4,A,F12.4)
          IERROR = 1
          Return
      Endif

C*****If V2P .LT. V1 then V1 has not been reached yet.  Get another panel.
      If(V2P .LT. V1) Go to 100

C*****Begin to accumulate the spectrum in SPECT
C*****Variables beginning with L refer to SPECT, with N refer to
C*****the input data in the array S
C*****LTOTAL = total number of points needed from input spectrum
C*****LRDATA = number of records containing real spectral data,
C*****         including partial records (which will be 0 or 1 filled.)
C*****LREC   = smallest power of 2 .GE. LRDATA = total number of
C*****         records on LFILE. Records in excess of LRDATA are all
C*****         0's or 1's
C*****LSLAST  = number of valid points in record LRDATA.
C*****Note: there is a potential roundoff problem using n = (v2-v1)/dv
C*****when n reaches 7 digits, the limit of single precision.  
      LTOTAL = (V2-V1)/DV+1.1
C*****Recompute V2, in case of round off errors
      V2 = V1+DV*(LTOTAL-1)

      ARDATA =  REAL(LTOTAL)/ REAL(LPTSMX)
      If((ARDATA-INT(ARDATA)) .EQ. 0.0) Then
          LRDATA = INT(ARDATA)
      Else
          LRDATA = INT(ARDATA+1)
      Endif

      POWER =  LOG( REAL(LRDATA))/ LOG(2.0)
      If(POWER-INT(POWER) .EQ. 0) Then
          LREC = 2**INT(POWER)
      Else
          LREC = 2**(INT(POWER)+1)
      Endif

      LSLAST = MOD(LTOTAL,LPTSMX)
      If(LSLAST .EQ. 0) Then
          LSLAST = LPTSMX
      Endif

C*****LTOTAL = cululative total points
C*****LS = cumulative total points currently in SPECT
      LTOTAL = 0
      LS = 0

C*****N1 = index of next point in S to go into SPECT
C*****From first panel, N1 will in general not be 1
C*****NN = number of unused points available from S = NP-N1+1
C*****This is a correction from the original fftscn for a
C*****potential roundoff problem that occurs when the total
C*****number of points needed reaches 7 digits, the limit
C*****of single precision
      N1 = (V1 - V1P)/DV+1.00001
      NN = NP-N1+1

C*****Load S into SPECT for a total of LRDATA blocks
      DO 200 L = 1,LRDATA

C*****    LRECMX = total points needed for this record
          If(L .EQ. LRDATA) Then
              LRECMX = LSLAST
          Else
              LRECMX = LPTSMX
          Endif

C*****    While SPECT not full, put S into SPECT
  110     Continue
          If((LRECMX-LS) .EQ. 0) GO TO 180

C*****    If more points are needed from S, get them
          If(N1 .GT. NP) Then
              Call GETPNL(IFILE,JCONVT,IEOFSC)
              If(IEOFSC .LE. 0) Then
                  Write(IPR,*) '  LOADSP: error, ran out of data',
     1          ' before V1, skipping to next scan request'
              IERROR = 1
              Return
              Endif
              N1 = 1
              NN = NP
          Endif

C*****    SPECT has LS points already.  LMAX more points are
C*****    needed.
          LMAX = LRECMX - LS

C*****    NMAX is lesser of LMAX = space available in SPECT
C*****    and NN = points available from S
          NMAX = MIN(LMAX,NN)

C*****    Put points NMAX points from S, starting at N1,
C*****    into SPECT, starting at LS+1
          Do 120 N = N1,N1+NMAX-1
              LS = LS+1
              SPECT(LS) = S(N)
  120     Continue    
          LTOTAL = LTOTAL+NMAX
          N1 = N1+NMAX
          NN = NN-NMAX

          Go To 110

C*****    SPECT is full
  180     Continue

          If(LS .GT. LPTSMX) Then
              Write(IPR,*) 'LOADSP: logic error, overfilled SPECT; ',
     1        ' LS, LPTSMX = ', LS,LPTSMX
              Stop 'Stopped in LOADSP'
          Endif

C*****    If this is the last record and it is not full, then fill
C*****        with 0's if radiance, or 1's if transmittance
          If(JEMIT .EQ. 0) Then
             FILL = 1.0
          Else
             FILL = 0.0
          Endif
          If((L .EQ. LRDATA) .AND. (LS .LT. LPTSMX)) Then
              Do 190 I=LS+1,LPTSMX
                 SPECT(I) = FILL
  190         Continue
          Endif

C*****    If LREC = 1, break out
          If(LREC .EQ. 1) Go To 300

C*****    If LREC > 1 and this is the first record
          If(L .EQ. 1) Then
              Open(LFILE,access='DIRECT',err=420,recl=IBLKSZ,
     1            status='SCRATCH')
          Endif
          WRITE(LFILE,rec=L) SPECT

C*****    End of write block.  Start new SPECT.
          LS = 0

C*****End of loop over data records
  200 Continue

C*****Add dummy records as needed
      If(LREC .GT. LRDATA) Then
          DO 210 I=1,LPTSMX
              SPECT(I) = FILL
  210     Continue

          Do 220 L=LRDATA+1,LREC
              WRITE(LFILE,rec=L,err=430) SPECT
  220     Continue
      Endif

C*****All done.
  300 Continue
      Return

  420 Continue
      Write(IPR,40) 'LOADSP: error in opening direct access file LFILE',
     1   LFILE
   40 Format(1X,A,I4)
      Stop 'Stopped in LOADSP'

  430 Write(IPR,50) ' LOADSP: error in writing to LFILE, record = ',
     1   LREC
   50 Format(A,I5)
      Stop 'Stopped in LOADSP'

      End
      Subroutine Multrn(LREC,LFILE1,LFILE2,FNCT1,FNCT2)
C***********************************************************************
C     This subroutine multiplies two functions FNCT1 and FNCT2 point
C     by point.  If LREC = 1 then the functions are entirely
C     contained in the arrays FNCT1 and FNCT2.  If LREC > 1
C     then the two functions stored on the files on units LFILE1 and
C     LFILE2. The result is stored either in FNCT2 or on file LFILE2.
C     It is assumed that the data is in full blocks of LPTSMX.
C***********************************************************************
C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 
C***********************************************************************
C     LPTSMX is the size of a data block for the FFT.  If the number of
C     data points is LPTSMX or less, the FFT is done in memory.  If it 
C     is larger, then a disk based FFT is performed, and LPTSMX is the 
C     size in words of each record of the file. The in-memory routine is 
C     about 20 percent more efficient than the disk-based routine, for 
C     same number of points.  For efficiency, LPTSMX should be made
C     as large a possible so that the computation can be done in memory.
C     However, on virtual memory machines (e.g. VAX), if LPTSMX is too
C     large, then the computation will be done in VIRTUAL memory.  Since
C     the in-memory routine uses widely scattered points, operating 
C     system will spend a great deal of time swapping data in and out 
C     (thrashing) and the efficiency will be very poor.
C
C     The following problem has been corrected (May, 1993).  The minimum
C     size of an FFT is now the smallest power of 2 greater than the 
C     of data points.  LPTSMX should be set to the larges value consistant
C     with real memory.
C
C*    However, LPTSMX is also the minimum size of an FFT (this is a 
C*    design flaw of this program).  If LPTSMX is much larger than the 
C*    number of points in the spectrum, computational time is wasted.
C*
C*    Therefore, LPTSMX should be set to a value somewhat larger than
C*    the size of the smallest typical spectrum, but no larger than the 
C*    largest value possible without thrashing.
C
C     IBLKSZ is the record length in the OPEN statement, and may be 
C     in bytes or words,  depending on the operating system.  For VAX
C     and CDC, it is in words.  For MicroSoft FORTRAN, it is in bytes.
C***********************************************************************
      PARAMETER (LSIZE = 65536)
C*****Following line for computers where the blocksize is measured
C*****in words, e.g. VAX, CDC/NOS/VE
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX,LPTSM8=LPTSMX/8)

C*****Following line for computers where the blocksize is measured 
C***** in bytes, e.g. MS FORTRAN, SUN, Alliant
      PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*4,LPTSM8=LPTSMX/8)

C*****Following line for computers with 64 bit words where the blocksize is 
C*****measured in bytes, e.g. CRAY
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*8,LPTSM8=LPTSMX/8)

      COMPLEX FNCT1(LPTSMX/2),FNCT2(LPTSMX/2)

C*****FNCT1 and FNCT2  are actually arrays of complex numbers of length
C*****LPTSMX/2, stored as arrays of real numbers of length LPTSMX.
C*****Since FNCTn is the FFT of a real array, only the positive
C*****frequencies are stored.  The imaginary part of FNCTn(1) is zero.
C*****For the in-memory routines, the real part of FNCTn(LPTSMX/2+1)
C*****(=V_MAX) is stored there.  Since the first element of FNCTn
C*****actually contains two REAL numbers, complex multiplication
C*****cannot be used.


      If(LREC .EQ. 1) Then
          FNCT2(1) = REAL(FNCT1(1))*REAL(FNCT2(1))+
     1               SQRT(CMPLX(-1.0))*AIMAG(FNCT1(1))*AIMAG(FNCT2(1))
          Do 100 I = 2,LPTSMX/2
              FNCT2(I) = FNCT1(I)*FNCT2(I)
  100     Continue

      Else
          I1 = 2
          Do 300 J = 1,LREC
               Read(LFILE1,rec=J,err=900) FNCT1
               Read(LFILE2,rec=J,err=910) FNCT2

               If(J .EQ. 1) THEN
                   FNCT2(1) = REAL(FNCT1(1))*REAL(FNCT2(1))+
     1                 SQRT(CMPLX(-1.0))*AIMAG(FNCT1(1))*AIMAG(FNCT2(1))
               Endif

               Do 200 I = I1,LPTSMX/2
                  FNCT2(I) = FNCT1(I)*FNCT2(I)
  200          Continue

               Write(LFILE2,rec=J,err=920) FNCT2
               I1 = 1
  300     Continue

       Endif

       Return

  900 Write(IPR,10) 'MULTRN: error in reading file LFILE1 = ',LFILE1
   10 Format(1X,A,I3)
      Stop 'Stopped in MULTRN'

  910 Write(IPR,10) 'MULTRN: error in reading file LFILE2 = ',LFILE2
      Stop 'Stopped in MULTRN'

  920 Write(IPR,20) 'MULTRAN: error in writing file LFILE2 = ',LFILE2
   20 Format(1X,'A',I3)
      Stop 'Stopped in MULTRN'

      End
      Subroutine Scnfnt(JFN,A,IAPSC,X,DX,LPTS,FUNCT,PARM1,PARM2,PARM3)
C************************************************************************
C     This function calculates the spectral scanning function (IAPSC=1) or
C     the equivalent apodization function (IAPSC=-1) corresponding to JFN,
C     characterized by the parameter A = 1/K, where K is the length of an 
C     equivalent interferometer. PARM1, PARM2, and PARM3 are parameters 
C     required to define some of the scanning functions.
C
C     The scanning functions and their corresponding apodization
C     functions are listed here. The constants Ci convert the half
C     width at half height to the "natural" measure of the width = a.
C     L = 1/a is typically the extent of the apodization function
C     and corresponds to the maximum path length difference of an
C     interferometer.
C
C     JFN  Scanning function
C      1   triangle= 1 - v/a, v < a,  a = C1*hwhm
C      2   Gaussian = exp(-.5*(v/a)**2), a = C2*hwhm
C      3   Sinc**2 = (sin(u)/u)**2, u = Pi*v/a, a = C3*hwhm
C      4   Sinc = sin(u)/u, u = 2*Pi*v/a, a = C4*hwhm
C      5   Beer = J(5/2,u)/u**(5/2) = ((3-u**2)sin(u)-3*u*cos(u))*u**-5
C          J(n,u) is Bessel fn of order n,augument u = 2*Pi*v/a, a=C5*hwhm
C      6   Hamming =sinc(u)+.428752(sinc(u+Pi)+sinc(u-Pi)),
C                     u = 2*Pi*v/a, a = C6*hwhm
C      7   Hanning = sinc(u)+.5*(sinc(u+Pi)+sinc(u-Pi)),
C                     u = 2*Pi*v/a, a = C7*hwhm
C      8   Norton-Beer: weak
C      9   Norton-Beer: moderate
c     10   Norton-Beer: strong
C     11   Brault
C     12   Kaiser-Bessel
C     13   Kiruna: c1*sinc(u)+c2*sinc(u-2*Pi*v_offset/a), u=2*Pi*v/a, a=1/L
C              This is an asymetric scanning function, the corresponding 
C              apodization function is complex.
C
C     JFN  Apodization Function
C      1   Sinc**2 = (sin(z)/z)**2, z = Pi*x*a (a = C1*hwhm)
C      2   Gaussian = exp(-2*Pi*(a*x)**2)
C      3   Triangle, 1 @x=0, 0 @x= 1/a
C      4   Rectangle = 1, x:0 to 1/a; 0, x > 1/a
C      5   Beer = (1-(x*a)**2)**2
C      6   Hamming=(1+.857504*cos(x*a))/1.857504
C      7   Hanning=(1+cos(x*a))/2
C      8   Norton-Beer: weak
C      9   Norton-Beer: moderate
c     10   Norton-Beer: strong
C     11   Brault = 1, 0 < x < p1/a
C                 = (cos(0.5*PI*(x-p1/a)(1/a-p1/a)))**2, p1 <= x < 1/a
C     12   Kaiser-Bessel
C                 = I0(PI*p1*sqrt((1.-(x*a)**2))/I0(Pi*p1), where IO is 
C                   zero-order modified Bessel function of the first kind,
C                   and p1 is a parameter (PARM1) from about 2 to 4
c     13   Kiruna: not implemented
C************************************************************************

C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 

      Dimension FUNCT(LPTS)
C*****The functions are: triangle, gauss, sinc**2, sinc, Beer, Hamming,
C*****and Hanning.  
C*****Function -1 is special.  It is a very broad sinc (only out to about
C*****to the first zero) which is the apodization function corresponding
C*****to a narrow rectangular scanning function used to pre-scan the spectrum.

      PI =2.0*ASIN(1.0)

      X0 = X

      If (JFN .EQ. -1) Go to 200
      Go to (10,20,30,40,50,60,70,80,80,80,90,100,110) JFN

      Write(IPR,*) 'Scnfnt - error: JFN out of range, = ', JFN
      Stop 'Stopped in Scnfnt'

   10 Continue
C*****Triangle
      If(IAPSC .EQ. -1) Then
C*****    Scanning function
          Do 15  L=1,LPTS
              If(X .GT. A) Then
                  FUNCT(L) = 0.0
              Else
                  FUNCT(L) = (1.0-X/A)
              Endif
              X = X0+L*DX
   15     Continue
      Else
C*****    Apodization function
          Do 16 L=1,LPTS-1,2
              If(X .EQ. 0) Then
                  FUNCT(L) = 1.0
              Else
                  FUNCT(L) = (SIN(Pi*X*A)/(Pi*X*A))**2
              Endif
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2.
   16     Continue
      Endif
      Return

   20 Continue
C*****Gaussian,
C*****ARGMIN is largest argument x to exp(-x)
      If(IAPSC .EQ. -1) Then
C*****    Scanning function
          Do 25 L=1,LPTS
              FUNCT(L) = exp( -MIN( 0.5*(X/A)**2, ARGMIN))
              X = X0+L*DX
   25     Continue
      Else
C*****    Apodization function
          Do 26 L=1,LPTS-1,2
              FUNCT(L) = exp( -MIN( 2.0*(Pi*X*A)**2, ARGMIN))
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2
   26     Continue
      Endif
      Return

   30 Continue
C*****Sinc**2
      If(IAPSC .EQ. -1) Then
C*****    Scanning function
          Do 35 L=1,LPTS
              If(X .NE. 0) Then
                  FUNCT(L) = (SIN(Pi*X/A)/(Pi*X/A))**2
              Else
                  FUNCT(L) = 1.0
              Endif
              X = X0+L*DX
   35 Continue
      Else
C*****    Apodization function
          Do 36 L=1,LPTS-1,2
              If(X .LT. 1/A) Then
                  FUNCT(L) = 1.0-X*A
              Else
                  FUNCT(L) = 0.0
              Endif
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2
   36     Continue
      Endif
      Return

   40 Continue
C*****Sinc
      If(IAPSC .EQ. -1) Then
C*****    Scanning function
          Do 45 L=1,LPTS
              If(X .NE. 0) Then
                  FUNCT(L) = SIN(2.0*Pi*X/A)/(2.0*Pi*X/A)
              Else
                  FUNCT(l) = 1.0
              Endif
              X = X0+L*DX
   45     Continue
      Else
C*****    Apodization function
C*****    This function is a rectangle of length L = 1/A.  Because the 
C*****    function is defined on a discreet grid, the length cannot be
C*****    exactly L.  To truncate the function at N = FIX(1/A)/DX would
C*****    lead to a discreet set of rectangles.  To avoid this problem, 
C*****    the last non-zero point is interpolated between 1 and 0 so that
C*****    the area under the rectangle is preserved. This procedure is at 
C*****    best an approximation but at least it provides a continious set
C*****    of functions.  
C*****    (By the way, the actual apodization function is the convolution of
C*****    a rectangle of length 1/A with a sinc(2*Pi*x*vmax), where vmax = 
C*****    is (V2 - V1), or the frequency range of the scanning calculation.  
C*****    The effect of the convolution with the sinc is seen only at the 
C*****    step at x = L as ringing.)
          XSTART = X
          ALEN = 1.0/A
          Do 46 L=1,LPTS-1,2
              If(X .LT. ALEN) Then
                  FUNCT(L) = 1.0
              Else
                  FUNCT(L) = 0.0
              Endif
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2.
   46     Continue
C*****    If the endpoint is in this block, correct it
C*****    ALEN = DX*(N+R), where 0.5 .LE. R .LT 1.5 and N is an integer
C*****    FUNCT(N+1) = NP-0.5
          N = Int((ALEN-XSTART)/DX)+1
          R = (ALEN-XSTART)/DX+1-N
          If (R .LE. 0.5) Then
              R = R+1.0   
              N = N-1
          Endif
          If(N .GE. 1  .AND.  N .LE. LPTS/2-1) Then
              FUNCT(2*N+1) = R-0.5
          Endif
      Endif
      Return

   50 Continue
C*****Beer
      If(IAPSC .EQ. -1) Then
C*****    Scanning function
          Do 55 L=1,LPTS
              U = 2.0*Pi*X/A
              If(U .GE. 0.5) Then
                  FUNCT(L) = 15.*((3.0-U**2)*SIN(U) - 3.0*U*COS(U))
     1                       /(U**5)
              Else
C*****        small angle approximation
                  FUNCT(L) = (1.0-15.0/210*U**2)
              Endif
              X = X0+L*DX
   55     Continue
      Else
C*****    Apodization function
          Do 56 L=1,LPTS-1,2
              If(X .LT. 1/A) Then
                  FUNCT(L) = (1.0 - (X*A)**2)**2
              Else
                  FUNCT(L) = 0.0
              Endif
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2.
   56     Continue
      Endif
      Return

   60 Continue
C*****Hamming
      If(IAPSC .EQ. -1) Then
C*****    Scanning function
          Do 65 L=1,LPTS
              U = 2.0*Pi*X/A
              If(U .NE. 0 .AND. ABS(U) .NE. Pi) Then
                  FUNCT(L) = (SIN(U)/U + 0.428752*
     1                       (SIN(U+Pi)/(U+Pi)+SIN(U-Pi)/(U-Pi)))
              Else If (U .EQ. 0) Then
                  FUNCT(L) = 1.0
              Else If (Abs(U) .EQ. Pi) Then
                  FUNCT(L) = 0.428752
              Endif
              X = X0+L*DX
   65     Continue
      Else
C*****    Apodization function
          Do 66 L=1,LPTS-1,2
              If(X .LT. 1/A) Then
                  FUNCT(L) = 0.53835685*(1.0+0.857504*COS(Pi*X*A))
              Else
                  FUNCT(L) = 0.0
              Endif
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2.
   66     Continue
      Endif

      Return

   70 Continue
C*****Hanning
      If(IAPSC .EQ. -1) Then
C*****    Scanning function
          Do 75 L=1,LPTS
              U = 2.0*Pi*X/A
              If(U .NE. 0 .AND. ABS(U) .NE. Pi) Then
                  FUNCT(L) = (SIN(U)/U + 0.5*
     1                       (SIN(U+Pi)/(U+Pi)+SIN(U-Pi)/(U-Pi)))
              Else If (U .EQ. 0) Then
                  FUNCT(L) = 1.0
              Else If (ABS(U) .EQ. Pi) Then
                  FUNCT(L) = 0.5
               Endif
              X = X0+L*DX
   75     Continue
      Else
C*****    Apodization function
          Do 76 L=1,LPTS-1,2
              If(X .LT. 1./A) Then
                  FUNCT(L) = 0.5*(1.0+COS(Pi*X*A))
              Else
                  FUNCT(L) = 0.0
              Endif
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2.
   76     Continue
      Endif
      Return

  80  Continue
C*****Norton-Beer
C     This is the generalized Norton-Beer function as described in:
C     R. H. Norton and R. Beer "New Apodizing Funcitons for Fourier
C     Spectroscopy", J. Opt. Soc. Am., #66, p259-264 (1976) (Corrected)

      If(IAPSC .EQ. -1) Then
C****     Scanning function
          Write(IPR,*) ' Scnfnt - error: Norton-Beer apodization', 
     c        '  not yet implemented in spectral domain'
          Stop ' Stopped in Scnfnt'
          Do 85 l=1,LPTS
              
              X = X0+L*DX
  85      Continue
      Else
C****     Apodization Function
          Do 89 L=1,LPTS-1,2
              U = 1.0-(X*A)**2
              If(X .LT. 1./A) Then
                  Goto (86,87,88) JFN-7
  86              Continue
C                 Weak Apodization
                  FUNCT(L) = 0.384093-0.087577*U+0.703484*U**2
                  Goto 189

  87              Continue
C                 Medium Apodization
                  FUNCT(L) = 0.152442-0.136176*U+0.983734*U**2
                  Goto 189

  88              Continue
C                 Strong Apodization
                  FUNCT(L) = 0.045355+0.554883*U**2+0.399782*U**4
                  Goto 189
              Else
                  FUNCT(L) = 0.0
              Endif
  189         Continue
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2.
  89      Continue
      Endif
      Return

  90  Continue
C*****Brault
C*****This function requires the parameter P and is defined by:
C*****    X = 0 to P*L, F = 1
C*****    X = P*L to L, F = ((cos(2*PI*(X-P*L)/(L-P*L))+1)/2)**2
C*****    X > L,        F = 0
C*****The valid range of P is [0,1].  A typical value is 0.9.
C*****Note: a value of P = 0 gives cos**2 apodization.

      If(PARM1 .LT. 0.0 .OR. PARM1 .GE. 1.0) Then
          Write(IPR,*) ' SCNFNT - Error: Brault Apodization, ',
     1        'P = ',PARM1,'  Valid range of P is [0,1)'
          Stop 'Stopped in SCNFNT'
      ENDIF

      If(IAPSC .EQ. -1) Then
C****     Scanning function
          Write(IPR,*) ' SCNFNT - Error: Brault Apodization', 
     c        '  not yet implemented in spectral domain'
          Stop ' Stopped in SCNFNT'

          Do 95 l=1,LPTS
C****         Insert code for scanning function here
              X = X0+L*DX
  95      Continue
      Else
C****     Apodization Function
          X1 = PARM1/A
          X2 = 1./A
          Do 99 L=1,LPTS-1,2
              If(X .LT. X1) Then
                  FUNCT(L) = 1.0
              Elseif (X .LT. X2) Then
                  FUNCT(L) = (COS(0.5*PI*(X-X1)/(X2-X1)))**2
              Else
                  FUNCT(L) = 0.0
              Endif
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2.
  99      Continue
      Endif
      Return

 100  Continue
C*****Kaiser-Bessel
C*****This function requires the parameter P and is defined by:
C*****    Y = PI*P**Sqrt(1-(x*A)**2)
C*****    F = Bessel(0,y), where Bessel(0,y) is the zero-order 
C*****        modified Bessel function of the first kind
C*****The valid range of P is [2,4].  

      If(PARM1 .LT. 2.0 .OR. PARM1 .GT. 4.0) Then
          Write(IPR,*) ' SCNFNT - Error: Kaiser-Bessel Apodization, ',
     1        'P = ',PARM1,'  Valid range of P is [2,4]'
          Stop 'Stopped in SCNFNT'
      ENDIF

      If(IAPSC .EQ. -1) Then
C****     Scanning function
          Write(IPR,*) ' Scnfnt - error: Kaiser-Bessel Apodization', 
     c        '  not yet implemented in spectral domain'
          Stop ' Stopped in Scnfnt'

          Do 105 l=1,LPTS
              
              X = X0+L*DX
 105      Continue
      Else
C****     Apodization Function
C****     Expansion for the Bessel function is from 'Numerical Recipies'
C****     Normalization factor: IO(X=0)
          FACTOR = BESSI0(PI*PARM1)
          Do 109 L=1,LPTS-1,2
              IF( X .LE. 1./A) Then
                  Y = PI*PARM1*SQRT(1.0-(X*A)**2)
                  FUNCT(L) = BESSI0(Y)/FACTOR
              Else
                  FUNCT(L) = 0.0
              Endif
              FUNCT(L+1) = 0.0
              X = X0+(L+1)*DX/2.
 109      Continue
      Endif
      Return

  110 Continue
C*****Kiruna (?)
C*****This is an asymetric scanning function consisting of a sinc at 0 
C*****and another sinc at v_offset= PARM1. The magnitudes of the two sinc's
C*****are c1 = PARM2 and c2 = PARM3.

      IF(IAPSC .EQ. -1) Then
C*****    Scanning Function
          Do 115 L=1,LPTS
              U1 = 2.0*Pi*X/A
              If (U1 .EQ. 0.0) THEN 
                  F1=1.0
              ELSE 
                  F1 = SIN(U1)/U1
              ENDIF 
                  
              U2 = 2.0*PI*(X-PARM1)/A
              IF (U2 .EQ. 0.0) THEN 
                  F2 = 1.0
              ELSE 
                  F2 = SIN(U2)/U2
              ENDIF

              FUNCT(L) = PARM2*F1+PARM3*F2
              X = X0+L*DX
  115     Continue
          
      Else
C*****    Apodization Function
C*****    The apodization function for this scanning function is complex
C*****    and is not implemented here.  
          Write(IPR,*) ' SCNFNT - ERROR: Kiruna Function: ',
     1        ' Apodization Function not defined'
              Stop 'Stopped in SCNFNT'

       Endif
       Return

  200 Continue
C*****Rectangle (Boxcar) and Sinc.
C*****This function corresponds to a broad sinc used to deconvolve 
C*****the spectrum from the effect of pre-scanning with a boxcar.  
C*****The transform of the smoothed, pre-scanned spectrum is divided 
C*****by this function which is the transform of the narrow rectangle.
C*****The width of this rectangle is constrained so that at the
C*****value of X at which the apodization function reaches zero 
C*****(i.e. X = length of an equivalent interferometer), the
C*****value of the sinc is 0.9.  The inverse of the sinc is returned 
C*****since this function will divide the transform of the spectrum.

      Do 206 L=1,LPTS-1,2
          Z = 2.0*Pi*X*A
          If(Z .NE. 0) Then
              FUNCT(L) = Z/SIN(Z)
          Else
              FUNCT(L) = 1.0
          Endif
          FUNCT(L+1) = 0.0
          X = X0+(L+1)*DX/2.
  206 Continue
      If(Z .GT. Pi) Then
          Write(IPR,*) ' Scnfnt - error: Function -1, Z is too large',
     1    ' for this approximation.  Prescanning rectangle is too ',
     2    ' wide?'
          Write(IPR,*) ' Z = ',Z
c          Stop 'Stopped in Scnfnt'
      Endif
      Return

      End
 
      FUNCTION BESSI0(X)
      REAL*8           P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     *    1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
      Subroutine Scntrn(JFN,A,IVX,DV,LREC,LPTFFT,LFILE,FUNCT,
     1    PARM1,PARM2,PARM3)
C***********************************************************************
C     This subroutine generates the Fourier transform of the spectral
C     scanning function (a.k.a. instrument response function,
C     apparatus function).  By analogy with Fourier transform spectro-
C     spopy, this function will be refered to as the apodization
C     function.  The apodization function can be generated in two ways:
C         a.  by first generating the spectral scanning function, and
C             taking the Fourier transform, or
C         b.  generating the apodization directly.
C     Approach b. is more efficient but there may be cases where a. is
C     preferable.  This subroutine is written to provide for both
C     approaches; the approach used is determined by IVX.
C
C     The parameters are as follows:
C     JFN specifies the shape of the spectral scanning function while
C     the parameter A = 1/K, where K is the length of an equivalent
C     interferometer, and determines the resolution.  A is related to
C     the halfwidth at half maximum by A = C(JFN)*HWHM.  If IVX = 1, 
C     then the apodization function is calculated directly, if -1, the 
C     scanning function is calculated then transformed. DV is the wavenumber
C     increment of the calculated spectrum, while LREC is the number of
C     records of length LPTSMX in the calculated spectrum. If LREC = 1,
C     then the apodization function is returned in the array spect,
C     otherwise it is written to the direct access file on LFILE in
C     blocks of LPSTMX.  PARM1, PARM2, and PARM3 are parameters used by
C     some of the scanning functions.
C
C     The organization of the arrays needs some explanation.  The scanning
C     function is a real function in the frequency domain, and is symmetric
C     around 0. Let it extend from V = -VMAX+DV to +VMAX in intervals
C     of DV, for a total of N = LREC*LPTSMX points, where VMAX = DV*N/2.
C     Since LPSTMX is a power of 2, N is even.  The positive and negative
C     frequencies of the scanning function are stored in the array FUNCT(I)
C     in the following manner.
C
C      I      V 
C
C      1      0  
C      2      DV
C      3      2*DV
C      .       .
C      N/2    (N/2-1)*DV
C      N/2+1  (N/2)*DV
C      N/2+2  -(N/2-1)*DV
C      N/2+3  -(N/2-3)*DV
C      .       .
C      N-1    -2*DV
C      N      -DV
C
C     Note that the scanning function is symmetric in V around 0.
C     FUNCT is symmetric around N/2+1 so that FUNCT(N+2-I) = FUNCT(I),
C     I = 2,N/2.  I = 1 corresponds to zero frequency, I=N/2+1 to VMAX.
C
C     The apodization function is the Fourier transform of the scanning
C     function so the it is a complex but Hermitan function in the space
C     domain.  Assume the space domain x extends from -Xmax to +Xmax.
C     The apodization function for negative x is the complex
C     conjugate of the apodization function for positive x.
C     Therefore, only the function for positive x need be stored, and
C     it is stored in a real array APOD(I) with the real parts in the
C     odd I's and the imaginary parts in the even.  APOD(0) is the real
C     value for x =0 and APOD(1) is the REAL value for Xmax (there is
C     no imaginary part for x = 0 or Xmax.)  If there are LPTSMX
C     real points in scanning function, then LPTSMX/2 complex points
C     are stored, so the total storage is the same.  The step size in
C     the space domain dX is 1/Vmax = 1/(LREC*LPTSMX*DV),
C     Xmax = 1/(2*DV).
C
C     The subroutine SCNFNT does the actual calculation of either the
C     scanning function or the apodization function, depending on the
C     value of IVX.
C
C     Note: 5/20/96
C     The capability of handling an asymetric scanning function has been
C     added, either by specifying scanning function 13 or by reading in
C     a scanning function.  The corresponding apodization function is 
C     complex but Hermitan (for a symetric scanning function, it is 
C     real.) The flag IASYM designates an asymetric scanning function.
C************************************************************************

C***********************************************************************
C     LPTSMX is the size of a data block for the FFT.  If the number of
C     data points is LPTSMX or less, the FFT is done in memory.  If it 
C     is larger, then a disk based FFT is performed, and LPTSMX is the 
C     size in words of each record of the file. The in-memory routine is 
C     about 20 percent more efficient than the disk-based routine, for 
C     same number of points.  For efficiency, LPTSMX should be made
C     as large a possible so that the computation can be done in memory.
C     However, on virtual memory machines (e.g. VAX), if LPTSMX is too
C     large, then the computation will be done in VIRTUAL memory.  Since
C     the in-memory routine uses widely scattered points, operating 
C     system will spend a great deal of time swapping data in and out 
C     (thrashing) and the efficiency will be very poor.
C
C     The following problem has been corrected (May, 1993).  The minimum
C     size of an FFT is now the smallest power of 2 greater than the 
C     of data points.  LPTSMX should be set to the larges value consistant
C     with real memory.
C
C*    However, LPTSMX is also the minimum size of an FFT (this is a 
C*    design flaw of this program).  If LPTSMX is much larger than the 
C*    number of points in the spectrum, computational time is wasted.
C*
C*    Therefore, LPTSMX should be set to a value somewhat larger than
C*    the size of the smallest typical spectrum, but no larger than the 
C*    largest value possible without thrashing.
C
C     IBLKSZ is the record length in the OPEN statement, and may be 
C     in bytes or words,  depending on the operating system.  For VAX
C     and CDC, it is in words.  For MicroSoft FORTRAN, it is in bytes.
C***********************************************************************
      PARAMETER (LSIZE = 65536)
C*****Following line for computers where the blocksize is measured
C*****in words, e.g. VAX, CDC/NOS/VE
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX,LPTSM8=LPTSMX/8)

C*****Following line for computers where the blocksize is measured 
C***** in bytes, e.g. MS FORTRAN, SUN, Alliant
      PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*4,LPTSM8=LPTSMX/8)

C*****Following line for computers with 64 bit words where the blocksize is 
C*****measured in bytes, e.g. CRAY
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*8,LPTSM8=LPTSMX/8)

C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 

      DIMENSION FUNCT(LPTSMX),FNEXT(1)

C*****IASYM flags asymetric scanning functions(=1)
      DIMENSION IASYM(13)
      DATA IASYM/12*0,1/

C*****Set LPTS equal to the  number of records per block.  For
C*****LREC = 1, this is LPTFFT which may be less than LPTSMX
      If(LREC .EQ. 1) Then
          LPTS = LPTFFT
      Else
          LPTS = LPTSMX
      Endif

      If(LREC .GT.1) Open(LFILE,access='DIRECT',err=900,recl=IBLKSZ,
     1    status='SCRATCH')

      If(IVX .EQ. -1) Then
C*****    Calculate the scanning function, then transform into the
C*****    apodization function.
C*****
C*****    Fill in FUNCT with the scanning function.
C*****    Loop over blocks.  FUNCT is symetric as described above.
C*****    If LREC is odd, then the point of symmetry is at the middle
C*****    of the middle block. If it is even, then the center of symmetry
C*****    is at the end of block LREC/2. If there is a middle block, it
C*****    must be treated separately.

C*****    First loop over full (not middle) blocks, JMAX of them. If JMAX
C*****    is 0 (ie, LREC = 1), skip this section.
          JMAX = LREC/2
          V = 0.0
          Do 120 J = 1,JMAX
              
              VSAVE = V
              Call Scnfnt(JFN,A,-1,V,DV,LPTS,FUNCT,PARM1,PARM2,PARM3)

C*****        Write positive frequencies.
              IREC = J
              Write(LFILE,rec=IREC,err=910) (FUNCT(I),I=1,LPTS)

C*****        Write negative frequencies.
              IF (IASYM(JFN) .EQ. 0) Then 
C*****            FNEXT is next value of FUNCT for positive frequencies
C*****            and is needed for the block with the negative components
                  VV = V
                  Call Scnfnt(JFN,A,-1,VV,DV,1,FNEXT,PARM1,PARM2,PARM3)
                  IREC = LREC+1-J
                  Write(LFILE,rec=IREC,err=910)
     1                 FNEXT(1),(FUNCT(I),I=LPTS,2,-1)
              Else
C*****            Asymetric Scanning function, negative frequncies
C*****            Start at -(VSAVE+DV) 
                  VNEG = -(VSAVE+DV)
                  Call Scnfnt(JFN,A,-1,-VNEG,-DV,LPTS,FUNCT,PARM1,
     1                PARM2,PARM3)
                  IREC = LREC+1-J
                  Write(LFILE,rec=IREC,err=910)
     1                (FUNCT(I),I=LPTS,1,-1)
 
                  V = VV
              Endif
  120     Continue

C*****    Calculate FUNCT for the middle block, if there is one (LREC
C*****    is odd.) Includes LREC = 1.          
          If(MOD(LREC,2) .EQ. 1) Then
              VSAVE = V
C*****        Calculate the positive frequencies.
              Call Scnfnt(JFN,A,-1,V,DV,LPTS/2+1,FUNCT,
     1                    PARM1,PARM2,PARM3)
              
C*****        Fill in the negative frequencies.
              If (IASYM(JFN) .EQ. 0) Then
                  Do 140 I=2,LPTS/2
                      FUNCT(LPTS+2-I) = FUNCT(I)
  140             Continue
                  IREC = JMAX+1
                  If(LREC .GT. 1) Write(LFILE,rec=IREC,err=910)
     1                (FUNCT(I),I=1,LPTS)
              Else
C*****        Need some fancy bookkeeping here. Fill in FUNCT from 
C*****        LPTS/2+2 to LPTS, starting at V = -(LPTS/2-1)*DV
                  V1 = -VSAVE-DV*(LPTS/2-1)
                  Call Scnfnt(JFN,A,-1,V1,DV,LPTS/2-1,FUNCT(LPTS/2+2),
     1                        PARM1,PARM2,PARM3)
                  
              Endif
          Endif

C*****    Take the Fourier transform of the scanning function to
C*****    get the apodization function.

          Call Fourtr(LREC,LPTFFT,LFILE,FUNCT,1)

C*****    Must normalize the scanning function to an area of 1.0.
C*****    The first coefficient of the apodization function is the
C*****    area of the scanning function.
          If(LREC .EQ. 1) Then
              FACTOR = FUNCT(1)
              Do 150 L=1,LPTS
                  FUNCT(L) = FUNCT(L)/FACTOR
  150         Continue
          Else
              Do 170 J=1,LREC
                  Read(LFILE,rec=J,err=920) (FUNCT(L),L=1,LPTS)
                  If(J .eq. 1) FACTOR = FUNCT(1)
                  Do 160 L=1,LPTS
                      FUNCT(L) = FUNCT(L)/FACTOR
  160             Continue
                  Write(LFILE,rec=J,err=910) (FUNCT(L),L=1,LPTS)
  170         Continue
          Endif
          
      Else If (IVX .EQ. 1) Then

C*****    Calculate the apodization function directly
C*****    X is in the space domain, Xmax = 1/dV, dX = 1/Vmax
          X= 0.0

          DX = 1.0/(LPTS*LREC*DV)
C*****    I think this is correct, but maybe it is the following?
C*****    DX = 1.0/((LPTS*LREC-1)*DV)

          Do 200 J=1,LREC
              Call Scnfnt(JFN,A,1,X,DX,LPTS,FUNCT,PARM1,PARM2,PARM3)
              IREC = J
              If(LREC .GT.1) Write(LFILE,rec=IREC,err=910)
     1                            (FUNCT(I),I=1,LPTS)
  200     Continue

      Else
          Write(IPR,*) 'Scntrn - error: IVX .NE. 1 or -1, = ',
     1                  IVX
          Stop 'Stopped in Scntrn'
      Endif
      Return

  900 Continue
      Write(IPR,*) 'Scntrn - error: error opening file: ',LFILE
      Stop 'Stopped in Scntrn'

  910 Continue
      Write(IPR,*) 'Scntrn - error: error writing to file ',LFILE,
     1             ', record number = ',IREC
      Stop 'Stopped in Scntrn'

  920 Continue
      Write(IPR,*) 'Scntrn - error: error reading file ',LFILE,
     1             ', record number = ',IREC
      Stop 'Stopped in Scntrn'


      End
      Subroutine Wrtspc(NV1,NTOTAL,LREC,LFILE,SPECT,JFILE)
C************************************************************************
C     This subroutine rewrites a spectral file from the direct access
C     format used by fftscn to the sequential LBLRTM format. 
C     The input data is on file LFILE, which contains LREC
C     records in blocks of LPTSMX points, or in the array SPECT, if
C     LREC = 1. The output file starts at point NV1 of the input 
c     file (corresponding to V1S) and contains a total of NTOTAL
C     points. NTOTAL may be less than LREC*LPTSMX because the 
C     limits of the scanned spectrum have been expanded to allow 
C     for edge effects and zero fill. If LREC = 1,  then the spectrum
C     is contained entirely in SPECT.  JFILE is the file number to 
C     be written to.
C************************************************************************

C*****Computers with 32 bit words need the Double Precision Statements
C*****Computers with 64 bit words (e.g. Cyber) do not.
C*****Frequency variables start with V
      Implicit Real*8           (V)
      Character*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ

C*****Blank Common carries the spectral data
      COMMON S(2450),R1(2650),XF(251)

C*****SCNHRD carries the header information for the scanned file
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF 
      DIMENSION FILHDR(2),IFSDID(17)
      EQUIVALENCE (FILHDR(1),XID(1))
      EQUIVALENCE (FSCDID(1),IFSDID(1),IHIRAC),(FSCDID(2),ILBLF4), 
     C (FSCDID(3),IXSCNT),(FSCDID(4 ),IAERSL ),(FSCDID(5),IEMIT), 
     C (FSCDID(6),ISCHDR ),(FSCDID(7 ),IPLOT  ),(FSCDID(8),IPATHL),   
     C (FSCDID(9),JRAD  ),(FSCDID(10),ITEST  ),(FSCDID(11),IMRG),  
     C (FSCDID(12),XSCID),(FSCDID(13),XHWHM  ),(FSCDID(14),IDABS),   
     C (FSCDID(15),IATM ),(FSCDID(16),LAYR1  ),(FSCDID(17),NLAYFS),  
     C (YID(1)    ,HDATE),(YID(2),      HTIME),(YI1,IMULT)        

C*****PANL carries the information from the panel header
      COMMON /PANL/ V1P,V2P,DVP,NP
      DIMENSION PNLHDR(4)
      EQUIVALENCE (PNLHDR(1),V1P)


C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 
C***********************************************************************
C     LPTSMX is the size of a data block for the FFT.  If the number of
C     data points is LPTSMX or less, the FFT is done in memory.  If it 
C     is larger, then a disk based FFT is performed, and LPTSMX is the 
C     size in words of each record of the file. The in-memory routine is 
C     about 20 percent more efficient than the disk-based routine, for 
C     same number of points.  For efficiency, LPTSMX should be made
C     as large a possible so that the computation can be done in memory.
C     However, on virtual memory machines (e.g. VAX), if LPTSMX is too
C     large, then the computation will be done in VIRTUAL memory.  Since
C     the in-memory routine uses widely scattered points, operating 
C     system will spend a great deal of time swapping data in and out 
C     (thrashing) and the efficiency will be very poor.
C
C     The following problem has been corrected (May, 1993).  The minimum
C     size of an FFT is now the smallest power of 2 greater than the 
C     of data points.  LPTSMX should be set to the larges value consistant
C     with real memory.
C
C*    However, LPTSMX is also the minimum size of an FFT (this is a 
C*    design flaw of this program).  If LPTSMX is much larger than the 
C*    number of points in the spectrum, computational time is wasted.
C*
C*    Therefore, LPTSMX should be set to a value somewhat larger than
C*    the size of the smallest typical spectrum, but no larger than the 
C*    largest value possible without thrashing.
C
C     IBLKSZ is the record length in the OPEN statement, and may be 
C     in bytes or words,  depending on the operating system.  For VAX
C     and CDC, it is in words.  For MicroSoft FORTRAN, it is in bytes.
C***********************************************************************
      PARAMETER (LSIZE = 65536)
C*****Following line for computers where the blocksize is measured
C*****in words, e.g. VAX, CDC/NOS/VE
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX,LPTSM8=LPTSMX/8)

C*****Following line for computers where the blocksize is measured 
C***** in bytes, e.g. MS FORTRAN, SUN, Alliant
      PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*4,LPTSM8=LPTSMX/8)

C*****Following line for computers with 64 bit words where the blocksize is 
C*****measured in bytes, e.g. CRAY
C     PARAMETER (LPTSMX=LSIZE,IBLKSZ=LPTSMX*8,LPTSM8=LPTSMX/8)

      Dimension SPECT(LPTSMX)
C*****NMAX is the maximum size of an LBLRTM panel
      Data NMAX/2400/

C*****Find the input record containing V1      
C*****If LREC > 1, then read to record containing V1
      NREC = NV1/LPTSMX+1
      If(LREC .GT. 1) Then
           Read(LFILE,rec=NREC,err=900) (SPECT(L),L=1,LPTSMX)
      Endif

C*****LLAST was the last point taken from SPECT
C*****LMORE points are available from SPECT
      LLAST = NV1-(NREC-1)*LPTSMX -1
      LMORE = LPTSMX-LLAST
      NREC = NREC+1

C*****NP is the number of points to be written to this output panel S.
C*****NLAST is the last point currently in S
C*****NDONE is the total number of points written already
      NP = MIN(NMAX,NTOTAL)
      NLAST = 0
      NDONE = 0

C*****While more points remain to be written, loop over output panels.
  120 Continue
      If(NDONE .LT. NTOTAL) Then

C*****    NP is the points needed for the next panel.
          NP = MIN(NMAX,NTOTAL-NDONE)

C*****    While S not full, fill it up
  150     Continue
          If(NLAST .LT. NP) Then

C*****        If more points are needed from SPECT, get them.
              If(LMORE .EQ. 0) Then
                  Read(LFILE,rec=NREC,err=900) (SPECT(L),L=1,LPTSMX)
                  LMORE = LPTSMX
                  LLAST = 0
                  NREC = NREC+1
              Endif

C*****        Add available points from SPECT to S until S is full
              NN = MIN(LMORE,NMAX-NLAST)
              N1 = NLAST + 1
              N2 = N1 + NN - 1
              L = LLAST
              DO 200 N=N1,N2
                  L = L+1
                  S(N) = SPECT(L)
  200         Continue
              NLAST = N2
              LLAST = L
              LMORE = LMORE-NN
              Go To 150
          Endif

C*****    Write out a panel
          If(NDONE .EQ. 0) Then
              V1P = V1C
              DVP = DV
          Else
              V1P = V2P+DVP
          Endif

          V2P = V1P+DVP*(NP-1)

          Call Bufout(JFILE,V1P,NPHDRF)
          Call Bufout(JFILE,S(1),NP)

          NDONE = NDONE+NP
          NLAST = 0
          Go to 120
      Endif

      CALL ENDFIL(JFILE)

      Return

  900 Continue
      Write(IPR,*) 'Wrtspc - err: error reading LFILE = ',LFILE
      Stop 'Stopped in Wrtspc'

      End
      Subroutine INTPDR(IFILE,JFILE,VV1,VV2,DVV,IFILST,JEMIT,IERR)
C**********************************************************************
C     This subroutine is a driver for the subroutine INTERP, which 
C     interpolates the spectral data on IFILE onto the wavenumber grid
C     VV1, VV2, and DVV and writes the result to JFILE.
C     JEMIT selects absorption (-1), transmittance (0), or radiance (1)
C     IERR is non-zero if an error occurs
C**********************************************************************

C*****Computers with 32 bit words need the Double Precision Statements
C*****Computers with 64 bit words (e.g. Cyber) do not.
C*****Frequency variables start with V
      Implicit Real*8           (V)
      Character*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ

C*****Blank Common carries the spectral data
      COMMON S(2450),R1(2650),XF(251)

C*****SCNHRD carries the header information for the scanned file
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF 
      DIMENSION FILHDR(2),IFSDID(17)
      EQUIVALENCE (FILHDR(1),XID(1))
      EQUIVALENCE (FSCDID(1),IFSDID(1),IHIRAC),(FSCDID(2),ILBLF4), 
     C (FSCDID(3),IXSCNT),(FSCDID(4 ),IAERSL ),(FSCDID(5),IEMIT), 
     C (FSCDID(6),ISCHDR ),(FSCDID(7 ),IPLOT  ),(FSCDID(8),IPATHL),   
     C (FSCDID(9),JRAD  ),(FSCDID(10),ITEST  ),(FSCDID(11),IMRG),  
     C (FSCDID(12),XSCID),(FSCDID(13),XHWHM  ),(FSCDID(14),IDABS),   
     C (FSCDID(15),IATM ),(FSCDID(16),LAYR1  ),(FSCDID(17),NLAYFS),  
     C (YID(1)    ,HDATE),(YID(2),      HTIME),(YI1,IMULT)        

C*****PANL carries the information from the panel header
      COMMON /PANL/ V1P,V2P,DVP,NP
      DIMENSION PNLHDR(4)
      EQUIVALENCE (PNLHDR(1),V1P)


C*****IFIL carries file information
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, 
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4

C*****LAMCHN carries hardware specific parameters
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 

C**********************************************************************
C*****THE INPUT DATA WILL BE PUT INTO S(5) = SS(1) WITH THE LAST
C*****4 POINTS OF THE PREVIOUS PANEL PUT INTO S(1 TO 4)
C*****THIS SCHEME PERMITS 6 POINT INTERPOLATION

C*****SS IS NOMINALLY 2401 POINTS BUT MAY NEED TO BE EXTENDED BY
C*****2 POINTS TO PERMIT 4 POINT INTERPOLATION UP TO THE LAST
C*****DATA POINT.

      DIMENSION SS(2406)
      EQUIVALENCE (S(5),SS(1))

C
C     Dimension RSTAT for use in INTERP
C
      DIMENSION RSTAT(3)


      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   J04110
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       J04120
     *               DVSC,XDUM,V1SHFT                                     J04130
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 J04170
      COMMON /STIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            J04180
      COMMON /INPNL/ V1I,V2I,DVI,NNI                                      J04190
      COMMON /OUTPNL/ V1J,V2J,DVJ,NNJ                                     J04200

      LOGICAL OP


      IERR = 0
      IBOUND = 4
      NPTS = 0

      Inquire(IFILE,OPENED=OP)
            If (.NOT. OP) Then
          Write(IPR,*) 'INTPDR: error: IFILE not open. IFILE = ',IFILE
          IERR = 1
          Return
      Endif

      Rewind IFILE
      If (IFILST .gt. 1) Then
          Call SKIPFL(IFILST-1,IFILE,IEOF)

          If (IEOF .eq. 0) Then
              Write(IPR,*) ' INTPDR: error: EOF skipping files'
              IERR = 1
          Endif
      Endif

      CALL GETHDR(IFILE,0,JDATA,IEOFSC)

C*****Ensure that the requested interpolation limits VV1 and VV2 are within
C*****the limits of the spectral data.
      If (VV1 .lt. V1C) Then
          AN = (V1C-VV1)/DVV
          If (AN .eq. INT(AN)) Then
              VV1 = V1+DVV*INT(AN)
          Else
              VV1 = V1+DVV*(INT(AN)+1)
          Endif
      Endif

      If (VV2 .gt. V2C) Then
          AN = (VV2-V2C)/DVV
          If (AN .eq. INT(AN))Then
              VV2 = VV2-DVV*INT(AN)
          Else
              VV2 = VV2-DVV*(INT(AN)+1)
          Endif
      Endif

C*****Reset frequency parameters in the file header
      V1C = VV1
      V2C = VV2
      DV = DVV
C*****Set V1 and V2, although I am not sure of all the implications
      V1 = VV1
      V2 = VV2
      DVO = DVV

C*****The following code was lifted from scnint.f.
C*****It sets various parameters in the file header and initializes
C*****the spectral buffer S. Note: the variables have been renamed from
C*****scnint: T -> S, and S -> SS
      ISCAN = ISCHDR                                                      I07370
      IF (ISCAN.LE.0.OR.XSCID.EQ.-99.) ISCAN = 0                          I07380
      IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR                   I07390
      ISCHDR = ISCAN+10                                                   I07400
C     V1C = V1                                                            I07410
C     V2C = V2                                                            I07420
C     DV = DVO                                                            I07430
C                                                                         I07440
      SCNID = 100*JEMIT                                                   I07450
      XSCID = SCNID+0.01                                                  I07460
C                                                                         I07470
      Call BUFOUT (JFILE,FILHDR(1),NFHDRF)                                I07480
C                                                                         I07490
      JTREM = -1                                                          I07500
      IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0                        I07510
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2                        I07520
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.2)) JTREM = 2                        I07530
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1                        I07540
      ISCANT = MOD(ISCAN,1000)                                            I07550
      IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2                       I07560
      IF (JTREM.LT.0) THEN                                                I07570
         WRITE(IPR,*) ' JTREM.LT.0 AT I07570'                             I07572
         STOP         ' JTREM.LT.0 AT I07570'
      ENDIF                                                                     
C     WRITE (IPR,910) IFILE,IEMIT,JEMIT,JTREM,JABS                        I07580
C                                                                         I07590
      IDATA = -1                                                          I07600
C                                                                         I07610
C     NEED TO SAVE LAST IBOUND POINTS OF EACH PANEL TO ATTACH TO NEXT     I07620
C                                                                         I07630
      IBOUND = 4                                                          I07640
C                                                                         I07650
C     VBOT IS LOWEST NEEDED WAVENUMBER, VTOP IS HIGHEST                   I07660
C                                                                         I07670
      BOUND =  REAL(IBOUND)*DV                                            I07680

C*****What if VBOT and VTOP are outside the limits of the data???

      VBOT = V1-BOUND                                                     I07690
      VTOP = V2+BOUND                                                     I07700

C                                                                         I07710
C     IF (JEMIT.EQ.0.AND.IDABS.EQ.0) BCD = HTRANS                         I07720
C     IF (JEMIT.EQ.0.AND.IDABS.EQ.-1) BCD = HABSRB                        I07730
C     IF (JEMIT.EQ.1) BCD = HRADIA                                        I07740
C     IF (NPTS.GT.0) WRITE (IPR,915) BCD                                  I07750
C                                                                         I07760
C     ZERO OUT T(1 TO IBOUND)                                             I07770
C                                                                         I07780
      DO 10 II = 1, IBOUND                                                I07790
         S(II) = 0.0                                                      I07800
   10 CONTINUE                                                            I07810
C                                                                         I07820
C     READ FROM IFILE UNTIL THE FIRST REQUIRED POINT IS REACHED           I07830
C     AND LOAD DATA INTO SS                                               I07840
C                                                                         I07850
      CALL RDPANL (SS,JTREM,IFILE,ISCAN,JEMIT,ICNVRT)                     I07860
      IF (IEOFSC.LE.0) GO TO 20                                           I07870
C                                                                         I07880
C     DO INTERPOLATION                                                    I07890
C     SET 4 POINT INTERPOLATION FLAG TO 1
C                                                                         I07900
      I4PT = 1
      CALL INTERP (IFILE,JFILE,I4PT,IBOUND,NPTS,JTREM,ISCAN,JEMIT,        I07910
     *             RSTAT,ICNVRT)                                          I07920
C                                                                         I07930
      Call ENDFIL(JFILE)

      Return

 20   Continue

      Write(IPR,*) 'INTPDR: error - EOF on input unit ',IFILE,
     1    ' before V1 was reached'
      IERR = 1

      Return

      End
      SUBROUTINE REALFT(DATA,N,ISIGN)
      REAL*8           WR,WI,WPR,WPI,WTEMP,THETA
      PARAMETER (PI2=6.28318530717959D0)
C*****PARAMETER (PI2=6.283185)
      DIMENSION DATA(*)
      
      THETA=PI2/(2.0*N)
      C1=0.5
      IF (ISIGN.EQ.1) THEN
          C2=-0.5
          CALL FOUR1(DATA,N,+1)
      ELSE
          C2=0.5
          THETA=-THETA
      ENDIF
      WPR=-2.0*SIN(0.5*THETA)**2
      WPI=SIN(THETA)
      WR=1.0+WPR
      WI=WPI
      N2P3=2*N+3
      DO 11 I=2,N/2+1
          I1=2*I-1
          I2=I1+1
          I3=N2P3-I2
          I4=I3+1
          WRS=WR
          WIS=WI
          H1R=C1*(DATA(I1)+DATA(I3))
          H1I=C1*(DATA(I2)-DATA(I4))
          H2R=-C2*(DATA(I2)+DATA(I4))
          H2I=C2*(DATA(I1)-DATA(I3))
          DATA(I1)=H1R+WRS*H2R-WIS*H2I
          DATA(I2)=H1I+WRS*H2I+WIS*H2R
          DATA(I3)=H1R-WRS*H2R+WIS*H2I
          DATA(I4)=-H1I+WRS*H2I+WIS*H2R
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
  11  CONTINUE
      IF (ISIGN.EQ.1) THEN
          H1R=DATA(1)
          DATA(1)=H1R+DATA(2)
          DATA(2)=H1R-DATA(2)
      ELSE
          H1R=DATA(1)
          DATA(1)=C1*(H1R+DATA(2))
          DATA(2)=C1*(H1R-DATA(2))
          CALL FOUR1(DATA,N,-1)
C*****Normalize
          DO 20 I=1,2*N
              DATA(I) = DATA(I)/ REAL(N)
  20      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      REAL*8           WR,WI,WPR,WPI,WTEMP,THETA
      PARAMETER (PI2=6.28318530717959D0)
C*****PARAMETER (PI2=6.2831853)
      DIMENSION DATA(*)
      
      N=2*NN
      J=1
      DO 11 I=1,N,2
          IF(J.GT.I)THEN
              TEMPR=DATA(J)
              TEMPI=DATA(J+1)
              DATA(J)=DATA(I)
              DATA(J+1)=DATA(I+1)
              DATA(I)=TEMPR
              DATA(I+1)=TEMPI
          ENDIF
          M=N/2
  1       IF ((M.GE.2).AND.(J.GT.M)) THEN
              J=J-M
              M=M/2
              GO TO 1
          ENDIF
          J=J+M
  11  CONTINUE
      MMAX=2
  2   IF (N.GT.MMAX) THEN
          ISTEP=2*MMAX
          THETA=PI2/(ISIGN*MMAX)
          WPR=-2.*SIN(0.5*THETA)**2
          WPI=SIN(THETA)
          WR=1.
          WI=0.
          DO 13 M=1,MMAX,2
              DO 12 I=M,N,ISTEP
                  J=I+MMAX
                  TEMPR=(WR)*DATA(J)-(WI)*DATA(J+1)
                  TEMPI=(WR)*DATA(J+1)+(WI)*DATA(J)
                  DATA(J)=DATA(I)-TEMPR
                  DATA(J+1)=DATA(I+1)-TEMPI
                  DATA(I)=DATA(I)+TEMPR
                  DATA(I+1)=DATA(I+1)+TEMPI
  12          CONTINUE
              WTEMP=WR
              WR=WR*WPR-WI*WPI+WR
              WI=WI*WPR+WTEMP*WPI+WI
  13      CONTINUE
          MMAX=ISTEP
          GO TO 2
      ENDIF
      RETURN
      END
      SUBROUTINE REALBG(LFILE,ISIGN,NBLK,IBLKSZ,A,B,S,IWORK)
C***********************************************************************
C     Calculates the FFT of real data.  Also does the inverse transform.
C     The algorithm this routine uses is the same as program REALFT in
C     the book Numerical Recipes The Art of Scientific Computing except
C     that it uses disk swapping.  Some of the actual code from REALFT
C     is used.  Like the Numerical Recipes routine
C     the last point is in data(2).  Unlike REALFT the inverse routine
C     takes care of normalizing by dividing by 1/N.
C     S and IWORK are not used in this routine they are just passed
C     through to BIGFFT.
C     Note of interest: Subroutine REALFT in my Numerical Recipes book i
C     not quite the same as it is on my Numerical Recipes disk.
C     
C     Mark P. Esplin
C     Stewart Radiance Lab.
C     139 The Great Rd.
C     Bedford Ma, 01730
C***********************************************************************
      Real*8           WR,WI,WPR,WPI,THETA
      DIMENSION A(*),B(*),S(*),IWORK(*)
      THETA=6.28318530717959D0/(NBLK*IBLKSZ)
      C1=0.5
      IF (ISIGN.EQ.1) THEN
          C2=-0.5
          CALL BIGFFT(LFILE,1,NBLK,IBLKSZ,A,B,S,IWORK)
      ELSE
          C2=0.5
          THETA=-THETA
      ENDIF
      IALK=1
      READ(LFILE,REC=IALK) (A(IOUT),IOUT=1,IBLKSZ)
C     First point is a special case.
      IF (ISIGN.EQ.1) THEN
          H1R=A(1)
          A(1)=H1R+A(2)
          A(2)=H1R-A(2)
      ELSE
          H1R=A(1)
          A(1)=C1*(H1R+A(2))
          A(2)=C1*(H1R-A(2))
      ENDIF
C     Set complex coefficients
      WPR=-2.0D0*SIN(0.5D0*THETA)**2
      WPI=SIN(THETA)
      WR=1.0D0+WPR
      WI=WPI
C     Do the actual work now.
      LSTBLK=NBLK/2+1
      DO 20 IBLK=NBLK,LSTBLK,-1
          READ(LFILE,REC=IBLK) (B(IOUT),IOUT=1,IBLKSZ)
          DO 10 IA=3,IBLKSZ,2
              IB=IBLKSZ-IA+2
              CALL PROSPN(A(IA),B(IB),WR,WI,WPR,WPI,C1,C2)
  10      CONTINUE
          WRITE(LFILE,REC=IALK) (A(IOUT),IOUT=1,IBLKSZ)
C     The last(first) point in each block is different
          IF(IBLK.EQ.LSTBLK) THEN
              CALL PROSPN(B(1),B(1),WR,WI,WPR,WPI,C1,C2)
          ELSE
C     Finish block B with the first data from next block A.
              IALK=IALK+1
              READ(LFILE,REC=IALK) (A(IOUT),IOUT=1,IBLKSZ)
              CALL PROSPN(A(1),B(1),WR,WI,WPR,WPI,C1,C2)
          ENDIF
          WRITE(LFILE,REC=IBLK) (B(IOUT),IOUT=1,IBLKSZ)
  20  CONTINUE
      IF(ISIGN.EQ.-1) THEN
          CALL BIGFFT(LFILE,-1,NBLK,IBLKSZ,A,B,S,IWORK)
      ENDIF
      RETURN
      END
C     
C     
      SUBROUTINE BIGFFT(LFILE,ISIGN,NBLK,IBLKSZ,A,B,S,IWORK)
C***********************************************************************
C     This routine preforms an in-place disk swapping complex FFT
C     on the data stored in a random access file on unit 4.
C     ISIGN = 1 for forward transform and ISIGN = -1 for inverse.
C     The data is stored in NBLK blocks of size IBLKSZ.
C     The data is stored such that the imaginary part of each complex
C     data point is followed by the real part.
C     The number of blocks of storage used on unit 4 = NBLK + NBLK/8.
C     The extra NBLK/8 blocks is used to store cosine sine tables.
C     A, B, S, and IWORK are work areas of size IBLKSZ.  S and IWORK
C     are not used at the same time, so they can share the same memory s
C     with an equivalence statement.
C     
C     Mark P. Esplin
C     Stewart Radiance Lab.
C     139 The Great Rd.
C     Bedford Ma, 01730
C***********************************************************************
C     
      DIMENSION A(*),B(*),S(*),IWORK(*)
      CALL SORTFF(LFILE,A,B,NBLK,IBLKSZ)
      CALL FFTBLK(LFILE,A,B,IWORK,ISIGN,NBLK,IBLKSZ)
C     Set sine-cosine table.
      CALL SCTABL(S,IBLKSZ)
      WRITE(LFILE,REC=NBLK+1) (S(I),I=1,IBLKSZ)
      CALL RECONS(LFILE,A,B,S,ISIGN,NBLK,IBLKSZ)
C     Renormalize inverse transform
      IF(ISIGN.EQ.-1) THEN
          DO 20 IBLK=1,NBLK
              READ(LFILE,REC=IBLK) (A(IOUT),IOUT=1,IBLKSZ)
              DO 10 I=1,IBLKSZ
                  A(I)=A(I)/NBLK
  10          CONTINUE
              WRITE(LFILE,REC=IBLK) (A(IOUT),IOUT=1,IBLKSZ)
  20      CONTINUE
      ENDIF
      RETURN
      END
C     
C     
C     
C     
C     
      
      SUBROUTINE PROSPN(A,B,WR,WI,WPR,WPI,C1,C2)
C     This routine processes a pair of complex points
C     and calculates the next complex coefficient.
      DIMENSION A(2),B(2)
      Real*8           WR,WI,WPR,WPI,WTEMP
      H1R=C1*(A(1)+B(1))
      H1I=C1*(A(2)-B(2))
      H2R=-C2*(A(2)+B(2))
      H2I=C2*(A(1)-B(1))
      WRS=WR
      WIS=WI
      A(1)=H1R+WRS*H2R-WIS*H2I
      A(2)=H1I+WRS*H2I+WIS*H2R
      B(1)=H1R-WRS*H2R+WIS*H2I
      B(2)=-H1I+WRS*H2I+WIS*H2R
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
      RETURN
      END
C     
C     
C     
C     
      SUBROUTINE SORTFF(LFILE,A,B,NBLK,IBLKSZ)
      COMPLEX A(1),B(1),H
C     THIS SUBROUTINE SORTS DATA TO GET READY FOR AN FFT
C     EVERY NBLK'TH POINT IS SORTED OUT INTO A SEPARATE FILE THEN
C     THE FILES ARE PUT IN BIT INVERTED ORDER
C     NOTE MEMORY REQUIRED FOR A COMPLEX ARRAY IS TWICE THAT OF A REAL
C     ARRAY.
C     INPUT IS IN THE FIRST NBLK BLOCKS OF TAPE4 A RANDOM ACCESS FILE
C     WITH A TOTAL NUMBER OF BLOCKS = NBLK+NBLK/8+1
C     THE DATA IS ALWAYS REWRITTEN IN PLACE
C     IBLKSZ = NUMBER OF STORAGE LOCATIONS IN A BLOCK
C     ICOMSZ = IBLKSZ/2 = NUMBER OF COMPLEX POINTS
C     NBLK = NUMBER OF BLOCKS
C     ICSUBZ = SIZE OF SUB-BLOCKS (COMPLEX) WHICH IS ICOMSZ/NBLK AT STAR
C     NPASS = NUMBER OF PASSES OVER THE DATA
C     THIS SECTION OF CODE SORTS OUT EVER NBLK'TH POINT IN A BLOCK AND
C     AND PUTS THE SUB-BLOCKS IN BIT INVERTED ORDER
      ICOMSZ=IBLKSZ/2
      ICSUBZ=ICOMSZ/NBLK
      NPASS=LLOG2(NBLK)
      DO 20 IBLK=1,NBLK
          READ(LFILE,REC=IBLK) (A(IOUT),IOUT=1,ICOMSZ)
          DO 10 ISUB=1,NBLK
              IP=ISUB-NBLK
              IB=INVER((ISUB-1),NPASS)*ICSUBZ
              DO 10 I=1,ICSUBZ
                  IP=NBLK+IP
  10          B(IB+I)=A(IP)
  20      WRITE(LFILE,REC=IBLK) (B(IOUT),IOUT=1,ICOMSZ)
C     THIS SECTION OF CODE SORTS THE BLOCKS
C     NSUB = NUMBER OF SUB-BLOCKS
C     NBB = NUMBER OF BLOCKS IN EACH BUTTERFLY
C     NBUTTE = NUMBER OF BUTTERFLIES
C     IBU = WHICH BLOCK IN BUTTERFLY
C     IBUTT = WHICH BUTTERFLY
C     IBA = BLOCK A
C     IBB = BLOCK B
      NSUB=NBLK
      NBB=1
      DO 50 IPAS=1,NPASS
          NBUTTE=NSUB/2
          NBB2=NBB*2
          ICSUB2=ICSUBZ*2
          DO 40 IBU=1,NBB
              IBA=IBU-NBB2
              DO 40 IBUTT=1,NBUTTE
                  IBA=IBA+NBB2
                  IBB=IBA+NBB
C     OPERATE ON BLOCK A AND BLOCK B
                  READ(LFILE,REC=IBA) (A(IOUT),IOUT=1,ICOMSZ)
                  READ(LFILE,REC=IBB) (B(IOUT),IOUT=1,ICOMSZ)
                  J=-ICSUB2
                  DO 30 ISUB=1,NSUB,2
                      J=J+ICSUB2
                      JI=J
                      DO 30 I=1,ICSUBZ
                          JI=JI+1
                          JI2=JI+ICSUBZ
                          H=A(JI2)
                          A(JI2)=B(JI)
  30              B(JI)=H
                  WRITE(LFILE,REC=IBA) (A(IOUT),IOUT=1,ICOMSZ)
  40      WRITE(LFILE,REC=IBB) (B(IOUT),IOUT=1,ICOMSZ)
          NBB=NBB2
          ICSUBZ=ICSUB2
  50  NSUB=NSUB/2
      RETURN
      END
C     
C     
C     
                  
      FUNCTION INVER(IVER,N)
C     THIS FUNCTION TAKES "IVER" AND INVERTS THE BITS
C     N IS THE NUMBER OF BITS IN THE WORD
C     THE VALUES OF IVER AND N ARE UNCHANGED
C     THIS FUNCTION WILL NOT WORK FOR NEGATIVE NUMBERS
C     This version uses standard FORTRAN 77.  It would be easier to
C     use a SHIFT and OR function, but they are not standard FORTRAN 77.
      IVE=IVER
      INV=0
      DO 10 I=1,N
          INV=2*INV+MOD(IVE,2)
          IVE=IVE/2
  10  CONTINUE
      INVER=INV
      RETURN
      END
C
C
      FUNCTION LLOG2(N)
C     FINDS INTEGER LOG TO BASE TWO
      LLOG2= LOG( REAL(N))/.69314718+.5
      RETURN
      END
C     
C     
C     
      SUBROUTINE FFTBLK(LFILE,A,S,INV,ISIGN,NBLK,IBLKSZ)
C     THIS ROUTINE DOES AN FFT OF EACH BLOCK WITH THE USE OF HARM1D
      DIMENSION A(*),S(*),INV(*)
      M=LLOG2(IBLKSZ)-1
C     SET UP SIN AND INV TABLES
      CALL HARM1D(A,M,INV,S,0,IFERR)
      IS=2*ISIGN
      DO 10 I=1,NBLK
          READ(LFILE,REC=I) (A(IOUT),IOUT=1,IBLKSZ)
          CALL HARM1D(A,M,INV,S,IS,IERR)
          WRITE(LFILE,REC=I) (A(IOUT),IOUT=1,IBLKSZ)
  10  CONTINUE
      RETURN
      END
C     
C     
C     
C     
C     
      SUBROUTINE SCTABL(S,IBLKSZ)
      DIMENSION S(*)
C     THIS ROUTINE GENERATES A SIN-COS TABLE
C     FROM THETA EQUAL ZERO TO 180 DEGREES
      IBLKS2=IBLKSZ/2
      IBLKS4=IBLKSZ/4
      DELT=4*ATAN(1.)/IBLKS2
      ROOT2=1./SQRT(2.)
      IAD2=IBLKS2+3
      IAD3=IBLKSZ+3
      S(1)=1.
      S(2)=0.
      S(IBLKS2+1)=0.
      S(IBLKS2+2)=1.
      DO 10 I=3,IBLKS4,2
          THA=DELT*(I/2)
          I2=I+1
          COST=COS(THA)
          SINT=SIN(THA)
          S(I)=COST
          S(I2)=SINT
          S(IAD2-I2)=SINT
          S(IAD2-I)=COST
          S(I+IBLKS2)=-SINT
          S(I2+IBLKS2)=COST
          S(IAD3-I2)=-COST
          S(IAD3-I)=SINT
  10  CONTINUE
      S(IBLKS4+1)=ROOT2
      S(IBLKS4+2)=ROOT2
      IBLKS34=IBLKS2+IBLKS4
      S(IBLKS34+1)=-ROOT2
      S(IBLKS34+2)=ROOT2
      RETURN
      END
C     
C     
C     
C     
C     
      SUBROUTINE RECONS(LFILE,A,B,S,ISIGN,NBLK,IBLKSZ)
      COMPLEX A(*),B(*),S(*)
C     THIS ROUTINE RECONSTRUCTS THE BLOCKS THAT HAVE ALREADY BEEN TRANSF
C     SEE SUBROUTINE SORTFF FOR DEFINITIONS OF TERMS
C     INPUT IS IN THE FIRST "NBLK" BLOCKS IN TAPE4
C     BLOCK NBLK+1 MUST CONTAIN A COS-SIN TABLE
C     NSIN= NUMBER OF TIMES THE SIGN TABLE NEEDS TO BE READ
C     NBLSIN= NUMBER OF BLOCKS IN THE COS-SIN TABLE WHEN IT GETS FILLED
C     ISIN= WHICH SIN BLOCK TO READ
C     ISRE= ISIN-NBLK
      ICOMSZ=IBLKSZ/2
      NPASS=LLOG2(NBLK)
      DELTHA=4.*ATAN(1.)/ICOMSZ
      DELT=DELTHA
      NBB=1
      NBUTTE=NBLK/2
      NBLSIN=NBLK/8
      IF(NBLSIN.LE.0) NBLSIN=1
      DO 50 IPASS=1,NPASS
C     TAKE CARE OF COS-SIN TABLE
          NSIN=NBB/4
          IF(NSIN.LE.0) NSIN=1
          NOFSIN=NBLSIN/NSIN
          NOFSI2=NOFSIN/2
          ISIN=NBLK+1-NOFSIN
          DO 40 ISRE=1,NSIN
              ISIN=ISIN+NOFSIN
              READ(LFILE,REC=ISIN) (S(IOUT),IOUT=1,ICOMSZ)
              IF(IPASS.EQ.NPASS) GO TO 10
C     GET TABLE FOR NEXT PASS
              DELTHA=DELTHA/2
              CALL EXPAN(A,B,S,DELTHA,IBLKSZ)
              WRITE(LFILE,REC=ISIN) (A(IOUT),IOUT=1,ICOMSZ)
C     FOR THE FIRST TWO PASSES YOU CAN PUT ALL THE INFORMATION IN THE
C     COS-SIN TABLE IN ONE BLOCK
              IF(IPASS.LE.2) GO TO 10
              ISIN2=ISIN+NOFSI2
              WRITE(LFILE,REC=ISIN2) (B(IOUT),IOUT=1,ICOMSZ)
C     DETERMINE WHICH BLOCKS TO WORK ON
C     THIS SHUFFLE IS TO MAKE IT SO ONLY ONE COS-SIN TABLE BLOCK IS
C     NEEDED TO WORK ON 8 DATA BLOCKS
  10         NBB2=NBB/2
             IBU=ISRE
             CALL HERFFT(LFILE,A,B,S,ISIGN,IBLKSZ,ICOMSZ,IBU,NBB,NBUTTE)
             IF(IPASS.LE.1) GO TO 40
             CALL ADDNIN(S,IBLKSZ)
             IBU=IBU+NBB2
             CALL HERFFT(LFILE,A,B,S,ISIGN,IBLKSZ,ICOMSZ,IBU,NBB,NBUTTE)
             IF(IPASS.LE.2) GO TO 40
             CALL COMPNIN(S,IBLKSZ,DELT)
             IBU=NBB2-ISRE+1
             CALL HERFFT(LFILE,A,B,S,ISIGN,IBLKSZ,ICOMSZ,IBU,NBB,NBUTTE)
             CALL ADDNIN(S,IBLKSZ)
             IBU=NBB-ISRE+1
             CALL HERFFT(LFILE,A,B,S,ISIGN,IBLKSZ,ICOMSZ,IBU,NBB,NBUTTE)
  40      CONTINUE
          NBB=NBB*2
          DELT=DELT/2
  50  NBUTTE=NBUTTE/2
      RETURN
      END
C     
      SUBROUTINE HERFFT(LFILE,A,B,S,ISIGN,IBLKSZ,ICOMSZ,IBU,NBB,
     1    NBUTTE)
C     THIS IS THE HEART OF THE EXTENDED FFT
C     SEE SORTFF FOR DEFINITION OF TERMS
      COMPLEX A(*),B(*),S(*),H,T
      NBB2=NBB*2
      IBA=IBU-NBB2
      DO 40 IBUTT=1,NBUTTE
          IBA=IBA+NBB2
          IBB=IBA+NBB
C     THIS SECTION OF CODE OPERATES ON TWO BLOCKS
          READ(LFILE,REC=IBA) (A(IOUT),IOUT=1,ICOMSZ)
          READ(LFILE,REC=IBB) (B(IOUT),IOUT=1,ICOMSZ)
          IF(ISIGN.EQ.1) THEN
C     The inverse was put on after the fact.
              DO 30 I=1,ICOMSZ
                  H=A(I)
                  T=B(I)*S(I)
                  A(I)=H+T
  30          B(I)=H-T
          ELSE
              DO 35 I=1,ICOMSZ
                  H=A(I)
                  T=B(I)*CONJG(S(I))
                  A(I)=H+T
  35          B(I)=H-T
          ENDIF
          WRITE(LFILE,REC=IBA) (A(IOUT),IOUT=1,ICOMSZ)
          WRITE(LFILE,REC=IBB) (B(IOUT),IOUT=1,ICOMSZ)
  40  CONTINUE
      RETURN
      END
C     
      SUBROUTINE EXPAN(A,B,S,DELTHA,IBLKSZ)
      DIMENSION A(*),B(*),S(*)
C     THIS ROUTINE EXPANDS THE SIN-COS TABLE BY TWO.
C     THE EXPANDED TABLE GOES INTO A AND B FROM THE ORIGINAL IN S
      COSD=COS(DELTHA)
      SIND=SIN(DELTHA)
C     DO FIRST HALF
      J=-1
      DO 10 I=1,IBLKSZ,4
          J=J+2
          J2=J+1
          COSJ=S(J)
          SINJ=S(J2)
          A(I)=COSJ
          A(I+1)=SINJ
          A(I+2)=COSJ*COSD-SINJ*SIND
          A(I+3)=COSJ*SIND+SINJ*COSD
  10  CONTINUE
C     DO SECOND PART
      DO 20 I=1,IBLKSZ,4
          J=J+2
          J2=J+1
          COSJ=S(J)
          SINJ=S(J2)
          B(I)=COSJ
          B(I+1)=SINJ
          B(I+2)=COSJ*COSD-SINJ*SIND
  20  B(I+3)=COSJ*SIND+SINJ*COSD
      RETURN
      END
C     
      SUBROUTINE ADDNIN(S,IBLKSZ)
      DIMENSION S(*)
C     THIS ROUTINE ADDS 90 TO THE ANGLE IN THE COS-SIN BLOCK
      DO 10 I=1,IBLKSZ,2
          I2=I+1
          T=S(I)
          S(I)=-S(I2)
  10  S(I2)=T
      RETURN
      END
C     
      SUBROUTINE COMPNIN(S,IBLKSZ,DELTHA)
      DIMENSION S(*)
C     THIS ROUTINE SUBTRACTS 90 FROM THE ANGLE THEN TAKES THE COS-SIN TA
C     FROM THA TO 90 MINUS THA.
      IBSZ2=IBLKSZ/2
      IEND=IBSZ2+1
      J2=IBLKSZ
      DO 10 I=3,IBSZ2,2
          I2=I+1
          J=J2-1
          T=S(I)
          S(I)=-S(J)
          S(J)=-T
          T=S(I2)
          S(I2)=S(J2)
          S(J2)=T
  10  J2=J2-2
      S(IEND)=-S(IEND)
C     FIND COS-SIN FOR THE FIRST POINT (THE ONE NOT IN OLD BLOCK)
      COSD=COS(DELTHA)
      SIND=SIN(DELTHA)
      S(1)=S(3)*COSD+S(4)*SIND
      S(2)=S(4)*COSD-S(3)*SIND
      RETURN
      END
      
      SUBROUTINE HARM1D(A,M,INV,S,IFSET, IFERR)
C     
C     ..................................................................
C     
C     SUBROUTINE HARM1D
C     
C     PURPOSE
C     PERFORMS DISCRETE COMPLEX FOURIER TRANSFORMS ON A COMPLEX
C     THREE DIMENSIONAL ARRAY
C     
C     USAGE
C     CALL HARM1D(A,M,INV,S,IFSET,IFERR)
C     
C     DESCRIPTION OF PARAMETERS
C     A     - AS INPUT, A CONTAINS THE COMPLEX, 1-DIMENSIONAL
C     ARRAY TO BE TRANSFORMED.  THE REAL PART OF
C     A(I1) IS STORED IN VECTOR FASHION IN A CELL
C     WITH INDEX 2*(I1) + 1 WHERE
C     NI = 2**MI AND I1 = 0,1,...,N1-1 ETC.
C     THE IMAGINARY PART IS IN THE CELL IMMEDIATELY
C     FOLLOWING.
C     THE NUMBER OF CORE LOCATIONS OF
C     ARRAY A IS 2*(N1)
C     M     - A ONE CELL VECTOR WHICH DETERMINES THE SIZES     
C     OF THE DIMENSIONS OF THE ARRAY A.   THE SIZE,     
C     NI, OF THE DIMENSION OF A IS 2**MI
C     INV   - A VECTOR WORK AREA FOR BIT AND INDEX MANIPULATION
C     OF DIMENSION ONE EIGHTH THE NUMBER OF CORE
C     LOCATIONS OF A, VIZ., (1/8)*2*N1
C     S     - A VECTOR WORK AREA FOR SINE TABLES WITH DIMENSION
C     THE SAME AS INV
C     IFSET - AN OPTION PARAMETER WITH THE FOLLOWING SETTINGS
C     0    SET UP SINE AND INV TABLES ONLY
C     1    SET UP SINE AND INV TABLES ONLY AND
C     CALCULATE FOURIER TRANSFORM
C     -1    SET UP SINE AND INV TABLES ONLY AND
C     CALCULATE INVERSE FOURIER TRANSFORM (FOR
C     THE MEANING OF INVERSE SEE THE EQUATIONS
C     UNDER METHOD BELOW)
C     2    CALCULATE FOURIER TRANSFORM ONLY (ASSUME
C     SINE AND INV TABLES EXIST)
C     -2    CALCULATE INVERSE FOURIER TRANSFORM ONLY
C     (ASSUME SINE AND INV TABLES EXIST)
C     IFERR - ERROR INDICATOR.   WHEN IFSET IS 0,+1,-1,
C     IFERR = 1 MEANS THE MAXIMUM M(I) IS LESS THAN 3
C     OR GREATER THAN 20, I=1,2,3  WHEN IFSET IS
C     +2,-2, IFERR = 1 MEANS THAT THE SINE AND INV
C     TABLES ARE NOT LARGE ENOUGH OR HAVE NOT BEEN
C     COMPUTED.  IF ON RETURN IFERR =0 THEN NONE OF
C     THE ABOVE CONDITIONS ARE PRESENT
C     
C     REMARKS
C     THIS SUBROUTINE IS TO BE USED FOR COMPLEX, 1-DIMENSIONAL
C     ARRAYS IN WHICH EACH DIMENSION IS A POWER OF 2.  THE
C     MAXIMUM MI MUST NOT BE LESS THAN 3 OR GREATER THAN 20.    
C     
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C     NONE
C     
C     THIS IS A 1-DIMENSIONAL MODIFICATION OF THE ORIGINAL 
C     3-DIMENSIONAL PROGRAM MODIFIED BY MARK P. ESPLIN.
C     
C     
C     SEE J.W. COOLEY AND J.W. TUKEY, "AN ALGORITHM FOR THE
C     MACHINE CALCULATION OF COMPLEX FOURIER SERIES",
C     MATHEMATICS OF COMPUTATIONS, VOL. 19 (APR. 1965), P. 297.
C     
C     ..................................................................
C     
C     The following SAVE statement needed on some machines,
C     specificlly Silicon Graphics. 
      SAVE NT,NTV2, MT
      DIMENSION A(1),INV(1),S(1),W(2),W2(2),W3(2)
  10  IF( IABS(IFSET) - 1) 900,900,12
  12  MTT=M-2
      ROOT2 = SQRT(2.)
      IF (MTT-MT ) 14,14,13
  13  IFERR=1
      RETURN
  14  IFERR=0
      M1=M
      N1=2**M1
  16  IF(IFSET) 18,18,20
  18  NX= N1
      FN = NX
      DO 19 I = 1,NX
          A(2*I-1) = A(2*I-1)/FN
  19  A(2*I) = -A(2*I)/FN
  20  NP=N1*2
      IL=0
      IL1=1
      MI=M1
  30  IDIF=NP
      KBIT=NP
      MEV = 2*(MI/2)
      IF (MI - MEV )60,60,40
C     
C     M IS ODD. DO L=1 CASE
  40  KBIT=KBIT/2
      KLAST=KBIT-1
      DO 50 K=1,KLAST,2
          KD=K+KBIT
C     
C     DO ONE STEP WITH L=1,J=0
C     A(K)=A(K)+A(KD)
C     A(KD)=A(K)-A(KD)
C     
          T=A(KD)
          A(KD)=A(K)-T
          A(K)=A(K)+T
          T=A(KD+1)
          A(KD+1)=A(K+1)-T
  50  A(K+1)=A(K+1)+T
  52  LFIRST =3
C     
C     DEF - JLAST = 2**(L-2) -1
      JLAST=1
      GO TO 70
C     
C     M IS EVEN
  60  LFIRST = 2
      JLAST=0
  70  DO 240 L=LFIRST,MI,2
          JJDIF=KBIT
          KBIT=KBIT/4
          KL=KBIT-2
C     
C     DO FOR J=0
          DO 80 I=1,IL1,IDIF
              KLAST=I+KL
              DO 80 K=I,KLAST,2
                  K1=K+KBIT
                  K2=K1+KBIT
                  K3=K2+KBIT
C     
                  T=A(K2)
                  A(K2)=A(K)-T
                  A(K)=A(K)+T
                  T=A(K2+1)
                  A(K2+1)=A(K+1)-T
                  A(K+1)=A(K+1)+T
C     
                  T=A(K3)
                  A(K3)=A(K1)-T
                  A(K1)=A(K1)+T
                  T=A(K3+1)
                  A(K3+1)=A(K1+1)-T
                  A(K1+1)=A(K1+1)+T
C     
                  T=A(K1)
                  A(K1)=A(K)-T
                  A(K)=A(K)+T
                  T=A(K1+1)
                  A(K1+1)=A(K+1)-T
                  A(K+1)=A(K+1)+T
C     
                  R=-A(K3+1)
                  T = A(K3)
                  A(K3)=A(K2)-R
                  A(K2)=A(K2)+R
                  A(K3+1)=A(K2+1)-T
  80      A(K2+1)=A(K2+1)+T
          IF (JLAST) 235,235,82
  82      JJ=JJDIF   +1
C     
          ILAST= IL +JJ
          DO 85 I = JJ,ILAST,IDIF
              KLAST = KL+I
              DO 85 K=I,KLAST,2
                  K1 = K+KBIT
                  K2 = K1+KBIT
                  K3 = K2+KBIT
C     
                  R =-A(K2+1)
                  T = A(K2)
                  A(K2) = A(K)-R
                  A(K) = A(K)+R
                  A(K2+1)=A(K+1)-T
                  A(K+1)=A(K+1)+T
C     
                  AWR=A(K1)-A(K1+1)
                  AWI = A(K1+1)+A(K1)
                  R=-A(K3)-A(K3+1)
                  T=A(K3)-A(K3+1)
                  A(K3)=(AWR-R)/ROOT2
                  A(K3+1)=(AWI-T)/ROOT2
                  A(K1)=(AWR+R)/ROOT2
                  A(K1+1)=(AWI+T)/ROOT2
                  T= A(K1)
                  A(K1)=A(K)-T
                  A(K)=A(K)+T
                  T=A(K1+1)
                  A(K1+1)=A(K+1)-T
                  A(K+1)=A(K+1)+T
                  R=-A(K3+1)
                  T=A(K3)
                  A(K3)=A(K2)-R
                  A(K2)=A(K2)+R
                  A(K3+1)=A(K2+1)-T
  85          A(K2+1)=A(K2+1)+T
              IF(JLAST-1) 235,235,90
  90          JJ= JJ + JJDIF
C     
C     NOW DO THE REMAINING J"S
              DO 230 J=2,JLAST
C     
C     FETCH W"S
C     DEF- W=W**INV(J), W2=W**2, W3=W**3
  96              I=INV(J+1)
  98              IC=NT-I
                  W(1)=S(IC)
                  W(2)=S(I)
                  I2=2*I
                  I2C=NT-I2
                  IF(I2C)120,110,100
C     
C     2*I IS IN FIRST QUADRANT
  100             W2(1)=S(I2C)
                  W2(2)=S(I2)
                  GO TO 130
  110             W2(1)=0.
                  W2(2)=1.
                  GO TO 130
C     
C     2*I IS IN SECOND QUADRANT
  120             I2CC = I2C+NT
                  I2C=-I2C
                  W2(1)=-S(I2C)
                  W2(2)=S(I2CC)
  130             I3=I+I2
                  I3C=NT-I3
                  IF(I3C)160,150,140
C     
C     I3 IN FIRST QUADRANT
  140             W3(1)=S(I3C)
                  W3(2)=S(I3)
                  GO TO 200
  150             W3(1)=0.
                  W3(2)=1.
                  GO TO 200
C     
  160             I3CC=I3C+NT
                  IF(I3CC)190,180,170
C     
C     I3 IN SECOND QUADRANT
  170             I3C=-I3C
                  W3(1)=-S(I3C)
                  W3(2)=S(I3CC)
                  GO TO 200
  180             W3(1)=-1.
                  W3(2)=0.
                  GO TO 200
C     
C     3*I IN THIRD QUADRANT
  190             I3CCC=NT+I3CC
                  I3CC = -I3CC
                  W3(1)=-S(I3CCC)
                  W3(2)=-S(I3CC)
  200             ILAST=IL+JJ
                  DO 220 I=JJ,ILAST,IDIF
                      KLAST=KL+I
                      DO 220 K=I,KLAST,2
                          K1=K+KBIT
                          K2=K1+KBIT
                          K3=K2+KBIT
C     
                          R=A(K2)*W2(1)-A(K2+1)*W2(2)
                          T=A(K2)*W2(2)+A(K2+1)*W2(1)
                          A(K2)=A(K)-R
                          A(K)=A(K)+R
                          A(K2+1)=A(K+1)-T
                          A(K+1)=A(K+1)+T
C     
                          R=A(K3)*W3(1)-A(K3+1)*W3(2)
                          T=A(K3)*W3(2)+A(K3+1)*W3(1)
                          AWR=A(K1)*W(1)-A(K1+1)*W(2)
                          AWI=A(K1)*W(2)+A(K1+1)*W(1)
                          A(K3)=AWR-R
                          A(K3+1)=AWI-T
                          A(K1)=AWR+R
                          A(K1+1)=AWI+T
                          T=A(K1)
                          A(K1)=A(K)-T
                          A(K)=A(K)+T
                          T=A(K1+1)
                          A(K1+1)=A(K+1)-T
                          A(K+1)=A(K+1)+T
                          R=-A(K3+1)
                          T=A(K3)
                          A(K3)=A(K2)-R
                          A(K2)=A(K2)+R
                          A(K3+1)=A(K2+1)-T
  220                 A(K2+1)=A(K2+1)+T
C     END OF I AND K LOOPS
C     
  230             JJ=JJDIF+JJ
C     END OF J-LOOP
C     
  235             JLAST=4*JLAST+3
  240         CONTINUE
C     END OF  L  LOOP
C     
C     WE NOW HAVE THE COMPLEX FOURIER SUMS BUT THEIR ADDRESSES ARE
C     BIT-REVERSED.  THE FOLLOWING ROUTINE PUTS THEM IN ORDER
              N1VNT=N1/NT
              JJD1=N1/(N1VNT*N1VNT)
              J=1
  800         JJ1=1
              DO 860 JPP1=1,N1VNT
                  IPP1=INV(JJ1)
                  DO 850 JP1=1,NT
  810                 IP1=INV(JP1)*N1VNT
  830                 I=2*(IPP1+IP1)+1
                      IF (J-I) 840,845,845
  840                 T=A(I)
                      A(I)=A(J)
                      A(J)=T
                      T=A(I+1)
                      A(I+1)=A(J+1)
                      A(J+1)=T
  845                 CONTINUE
  850             J=J+2
  860         JJ1=JJ1+JJD1
C     END OF JPP1 AND JP2
C     
C     
  890         IF(IFSET)891,895,895
  891         DO 892 I = 1,NX
  892         A(2*I) = -A(2*I)
  895         RETURN
C     
C     THE FOLLOWING PROGRAM COMPUTES THE SIN AND INV TABLES.
C     
  900         MT=M-2
              MT = MAX0(2,MT)
  904         IF (MT-20)906,906,905
  905         IFERR = 1
              GO TO 895
  906         IFERR=0
              NT=2**MT
              NTV2=NT/2
C     
C     SET UP SIN TABLE
C     THETA=PIE/2**(L+1) FOR L=1
  910         THETA=.7853981634
C     
C     JSTEP=2**(MT-L+1) FOR L=1
              JSTEP=NT
C     
C     JDIF=2**(MT-L) FOR L=1
C     JDIF=2**(MT-L) FOR L=1
              JDIF=NTV2
              S(JDIF)=SIN(THETA)
              DO 950 L=2,MT
                  THETA=THETA/2.
                  JSTEP2=JSTEP
                  JSTEP=JDIF
                  JDIF=JSTEP/2
                  S(JDIF)=SIN(THETA)
                  JC1=NT-JDIF
                  S(JC1)=COS(THETA)
                  JLAST=NT-JSTEP2
                  IF(JLAST - JSTEP) 950,920,920
  920             DO 940 J=JSTEP,JLAST,JSTEP
                      JC=NT-J
                      JD=J+JDIF
  940             S(JD)=S(J)*S(JC1)+S(JDIF)*S(JC)
  950         CONTINUE
C     
C     SET UP INV(J) TABLE
C     
  960         MTLEXP=NTV2
C     
C     MTLEXP=2**(MT-L). FOR L=1
              LM1EXP=1
C     
C     LM1EXP=2**(L-1). FOR L=1
              INV(1)=0
                  DO 980 L=1,MT
                      INV(LM1EXP+1) = MTLEXP
                      DO 970 J=2,LM1EXP
                          JJ=J+LM1EXP
  970                 INV(JJ)=INV(J)+MTLEXP
                      MTLEXP=MTLEXP/2
  980             LM1EXP=LM1EXP*2
  982             IF(IFSET)12,895,12
                  END
