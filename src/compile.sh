#!/bin/bash
#
REFD=$(pwd)
REFD="$REFD"/../bin/

#
fpc -Mobjfpc BSPabstand.p # bin
mv BSPabstand $REFD/BSPabstand
#
fpc -Mobjfpc BSScherung.p # bin
mv BSScherung $REFD/BSScherung
#
fpc -Mobjfpc PPScherung.p  #bin
mv PPScherung $REFD/PPScherung
#
fpc -Mobjfpc best_plane.p  #bin 
mv best_plane $REFD/best_plane
#
fpc -Mobjfpc dehedral.p  #bin
mv dehedral $REFD/dehedral
#
fpc -Mobjfpc evaluate.p #bin
mv evaluate $REFD/evaluate

# removes object files after compiling with the `free pascal compiler`
rm *.o
