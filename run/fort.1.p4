XCtype 'GGAPBE96'                 / XCtype
NBAND 576                         / Nband
NGRID 90 90 90                    / Ngrid(1:3)
PP 2 'Si-PBE-TM.psv'              / file_ps(1)
PP 2 'O_PBE_RC1.5.psv'            / file_ps(2)
scfconv 1.d-15                    / scf_conv
PROCS  2 1 1 2 1 1                / np_2d(1:6)

#IDIAG "eigen_s"                  /
IMATE   1                         / [ 0 - 3 ]
IROTV   2                         / [ 0 - 2 ]
#
#IOCTRL 3                         / parallel IO for WF
IC      0                         / IC
OC      0                         / OC
SWSCF   1                         / iswitch_scf
DITER   50                        / Diter
NSWEEP  0                         / Nsweep
SWOPT   0                         / cpmd
MBD     8                         /
ETLIMIT 3600.0                    /

GS      11
NBLK    144
SCL     2 2 144 144 2 2

