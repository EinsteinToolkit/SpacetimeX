# gnuplot command file to produce postscript file for 3patch figure

set term postscript eps color solid
set size 1.5, 2.25
set output 'AEIDevelopment_AHFinderDirect_3patch.eps'

unset key
set view 60, 111
set border 0

unset xtics
unset ytics
unset ztics

splot '-' w lines -1,							\
      '-' w lines -1,							\
      '-' w lines -1,							\
      '<../src/misc/select.patch +z <3patch.h.t0.ah1.gp'		\
						using 4:5:6 w lines 1,	\
      '<../src/misc/select.patch +x <3patch.h.t0.ah1.gp'		\
						using 4:5:6 w lines 2,	\
      '<../src/misc/select.patch +y <3patch.h.t0.ah1.gp'		\
						using 4:5:6 w lines 3,	\
      "<../src/misc/select.patch +z <3patch.h.t0.ah1.gp			\
        | gawk '($4 == 0) && ($5 >= 0.0)' -"				\
						using 4:5:6 w lines -1,	\
      "<../src/misc/select.patch +z <3patch.h.t0.ah1.gp			\
        | gawk '($4 >= 0) && ($5 == 0.0)' -"				\
						using 4:5:6 w lines -1,	\
      "<../src/misc/select.patch +x <3patch.h.t0.ah1.gp			\
        | gawk '($5 == 0) && ($6 >= 0.0)' -"				\
						using 4:5:6 w lines -1,	\
      "<../src/misc/select.patch +x <3patch.h.t0.ah1.gp			\
        | gawk '($5 >= 0) && ($6 == 0.0)' -"				\
						using 4:5:6 w lines -1,	\
      "<../src/misc/select.patch +y <3patch.h.t0.ah1.gp			\
        | gawk '($4 == 0) && ($6 >= 0.0)' -"				\
						using 4:5:6 w lines -1,	\
      "<../src/misc/select.patch +y <3patch.h.t0.ah1.gp			\
        | gawk '($4 >= 0) && ($6 == 0.0)' -"				\
						using 4:5:6 w lines -1
0.0	0.0	0.0
1.3	0.0	0.0
eof
0.0	0.0	0.0
0.0	1.3	0.0
eof
0.0	0.0	0.0
0.0	0.0	1.3
eof

set output
