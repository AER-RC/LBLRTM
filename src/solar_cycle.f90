!  --------------------------------------------------------------------------
! |  Copyright ©, Atmospheric and Environmental Research, Inc., 2015         |
! |                                                                          |
! |  All rights reserved. This source code is part of the LBLRTM software    |
! |  and is designed for scientific and research purposes. Atmospheric and   |
! |  Environmental Research, Inc. (AER) grants USER the right to download,   |
! |  install, use and copy this software for scientific and research         |
! |  purposes only. This software may be redistributed as long as this       |
! |  copyright notice is reproduced on any copy made and appropriate         |
! |  acknowledgment is given to AER. This software or any modified version   |
! |  of this software may not be incorporated into proprietary software or   |
! |  commercial software offered for sale without the express written        |
! |  consent of AER.                                                         |
! |                                                                          |
! |  This software is provided as is without any express or implied          |
! |  warranties.                                                             |
! |                       (http://www.rtweb.aer.com/)                        |
!  --------------------------------------------------------------------------
!
MODULE solar_cycle   ! Constants needed to to scale Kurucz and NRLSSI2 extraterrestrial solar irradiance

!
!
   real, parameter  :: scon_kurucz = 1368.22     ! W/m2
   !should be 1360.95
   real, parameter  :: scon_nrlssi2 = 1360.85     ! W/m2
   real, parameter  :: solcyc_min = 0.01890       !Solar cycle fraction at solar minimum
   real, parameter  :: solcyc_max = 0.3750        !Solar cycle fraction at solar maximum
!
! Mean quiet sun, facular brightening, and sunspot dimming coefficient terms (NRLSSI2, 100-50000 cm-1),
! spectrally integrated from hi-res values
   real, parameter :: Iint = 1360.375     ! Solar quiet term integrated
   real, parameter :: Fint = 0.996047     ! Solar facular brightening term (index-offset), integrated
   real, parameter :: Sint = -0.511590    ! Solar sunspot dimming term (index-offset), integrated
   real, parameter :: Foffset = 0.14959542    ! Solar facular offset
   real, parameter :: Soffset = 0.00066696    ! Solar sunspot offset

! Mg and SB indices for average solar cycle integrated over solar cycle
   real, parameter :: svar_f_avg = 0.1567652     ! Solar variability NRLSSI2 Mg "Bremen" index
   !  time-averaged over Solar Cycles 13-24
   !  and averaged over solar cycle
   real, parameter :: svar_s_avg = 902.71260     ! Solar variability NRLSSI2 SB "SPOT67" index
   !  time-averaged over Solar Cycles 13-24
   !  and averaged over solar cycle


! Mg and SB index look-up tables for average solar cycle as a function of solar cycle

   integer, parameter :: nsolfrac=134        ! Number of elements in solar arrays (12 months
   !  per year over 11-year solar cycle)
   real :: mgavgcyc(nsolfrac)               ! Facular index from NRLSSI2 Mg "Bremen" index
   !  time-averaged over Solar Cycles 13-24
   real :: sbavgcyc(nsolfrac)               ! Sunspot index from NRLSSI2 SB "SPOT67" index
   !  time-averaged over Solar Cycles 13-24
   data mgavgcyc /    0.150737,                                                         &
   &   0.150737,  0.150733,  0.150718,  0.150725,  0.150762,  0.150828, &
   &   0.150918,  0.151017,  0.151113,  0.151201,  0.151292,  0.151403, &
   &   0.151557,  0.151766,  0.152023,  0.152322,  0.152646,  0.152969, &
   &   0.153277,  0.153579,  0.153899,  0.154252,  0.154651,  0.155104, &
   &   0.155608,  0.156144,  0.156681,  0.157178,  0.157605,  0.157971, &
   &   0.158320,  0.158702,  0.159133,  0.159583,  0.160018,  0.160408, &
   &   0.160725,  0.160960,  0.161131,  0.161280,  0.161454,  0.161701, &
   &   0.162034,  0.162411,  0.162801,  0.163186,  0.163545,  0.163844, &
   &   0.164029,  0.164054,  0.163910,  0.163621,  0.163239,  0.162842, &
   &   0.162525,  0.162344,  0.162275,  0.162288,  0.162369,  0.162500, &
   &   0.162671,  0.162878,  0.163091,  0.163251,  0.163320,  0.163287, &
   &   0.163153,  0.162927,  0.162630,  0.162328,  0.162083,  0.161906, &
   &   0.161766,  0.161622,  0.161458,  0.161266,  0.161014,  0.160666, &
   &   0.160213,  0.159690,  0.159190,  0.158831,  0.158664,  0.158634, &
   &   0.158605,  0.158460,  0.158152,  0.157691,  0.157152,  0.156631, &
   &   0.156180,  0.155827,  0.155575,  0.155406,  0.155280,  0.155145, &
   &   0.154972,  0.154762,  0.154554,  0.154388,  0.154267,  0.154152, &
   &   0.154002,  0.153800,  0.153567,  0.153348,  0.153175,  0.153044, &
   &   0.152923,  0.152793,  0.152652,  0.152510,  0.152384,  0.152282, &
   &   0.152194,  0.152099,  0.151980,  0.151844,  0.151706,  0.151585, &
   &   0.151496,  0.151437,  0.151390,  0.151347,  0.151295,  0.151220, &
   &   0.151115,  0.150993,  0.150883,  0.150802,  0.150752,  0.150737, &
   &   0.150737 /
   data sbavgcyc  /     50.3550,                                                         &
   &    50.3550,   52.0179,   59.2231,   66.3702,   71.7545,   76.8671, &
   &    83.4723,   91.1574,   98.4915,  105.3173,  115.1791,  130.9432, &
   &   155.0483,  186.5379,  221.5456,  256.9212,  291.5276,  325.2953, &
   &   356.4789,  387.2470,  422.8557,  466.1698,  521.5139,  593.2833, &
   &   676.6234,  763.6930,  849.1200,  928.4259,  994.9705, 1044.2605, &
   &  1087.5703, 1145.0623, 1224.3491, 1320.6497, 1413.0979, 1472.1591, &
   &  1485.7531, 1464.1610, 1439.1617, 1446.2449, 1496.4323, 1577.8394, &
   &  1669.5933, 1753.0408, 1821.9296, 1873.2789, 1906.5240, 1920.4482, &
   &  1904.6881, 1861.8397, 1802.7661, 1734.0215, 1665.0562, 1608.8999, &
   &  1584.8208, 1594.0162, 1616.1486, 1646.6031, 1687.1962, 1736.4778, &
   &  1787.2419, 1824.9084, 1835.5236, 1810.2161, 1768.6124, 1745.1085, &
   &  1748.7762, 1756.1239, 1738.9929, 1700.0656, 1658.2209, 1629.2925, &
   &  1620.9709, 1622.5157, 1623.4703, 1612.3083, 1577.3031, 1516.7953, &
   &  1430.0403, 1331.5112, 1255.5171, 1226.7653, 1241.4419, 1264.6549, &
   &  1255.5559, 1203.0286, 1120.2747, 1025.5101,  935.4602,  855.0434, &
   &   781.0189,  718.0328,  678.5850,  670.4219,  684.1906,  697.0376, &
   &   694.8083,  674.1456,  638.8199,  602.3454,  577.6292,  565.6213, &
   &   553.7846,  531.7452,  503.9732,  476.9708,  452.4296,  426.2826, &
   &   394.6636,  360.1086,  324.9731,  297.2957,  286.1536,  287.4195, &
   &   288.9029,  282.7594,  267.7211,  246.6594,  224.7318,  209.2318, &
   &   204.5217,  204.1653,  200.0440,  191.0689,  175.7699,  153.9869, &
   &   128.4389,  103.8445,   85.6083,   73.6264,   64.4393,   50.3550, &
   &    50.3550/



end module solar_cycle
