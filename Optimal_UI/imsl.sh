#!/bin/bash
source /opt/imsl//fnl700/macab115e64/bin/fnlsetup.sh
export IMSL_F90FLAGS=$F90FLAGS
/Applications/Absoft16.0//bin/amake2  -k -f makefile.amake
