#
#    C1'---C1'
#!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#!!!!!!                DO NOT TOUCH BELOW                   !!!!!!
#!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
set NBP=`wc -l bp | awk '{print $1 }' ` # Anzahl Basenpaare
#
#!  PDB File Vorbereiten
#@@@@@@@@@@@@@@@@@@@@@@

echo $PDB > pdb
echo " "  >> pdb

nawk -v CMDFILE="pdb" -f "$REFD"prepPDB.awk > inn.pdb
rm pdb
#
sed s/\'/\*/g inn.pdb > in.pdb
rm inn.pdb
#
#!  BK Atomen Extrahiern
#@@@@@@@@@@@@@@@@@@@@@@
echo "in.pdb" > pdb
echo "$BK"    >> pdb
echo " "      >> pdb

nawk -v CMDFILE="pdb" -f "$REFD"getC.awk > P.pdb
rm pdb

sed -n  1,"$P1"p  P.pdb > P1.pdb  # P Atome vom Strang 1
sed     1,"$P1"d  P.pdb > P2.pdb  # P Atome vom Strang 2

#!  Berechne Abstand P-Atome BSpline
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo " Abstand Backbone-Atome BSpline Kurve " > ppp
echo " " >> ppp
#____ 1 Strang _____
cp  P1.pdb input1
"$REFD"BSPabstand << eop >> ppp
$KK
eop
rm input1
echo " " >> ppp
#____ 2 Strang _____
cp  P2.pdb input1
"$REFD"BSPabstand  << eop >> ppp
$KK
eop
rm input1
#
#!  Basenpaare Atome extrahieren
#@@@@@@@@@@@@@@@@@@@@@@
#
set count = 0
set loop  = $NBP
#
#
Begin:
#
@ count = $count + 1
#
set A=`sed -n "$count"p bp | awk '{print $1 }' `
set B=`sed -n "$count"p bp | awk '{print $2 }' `
set C=`sed -n "$count"p bp | awk '{print $3 }' `
set D=`sed -n "$count"p bp | awk '{print $4 }' `
#
echo "in.pdb" >  jjj
echo $A       >> jjj
echo $B       >> jjj
echo $C       >> jjj
echo $D       >> jjj

nawk -v CMDFILE="jjj" -f "$REFD"getBP.awk > BP_$A$B-$C$D.pdb

nawk -v CMDFILE="jjj" -f "$REFD"getB.awk > BP_$A$B.pdb

sed 2,3d jjj > jj

nawk -v CMDFILE="jj" -f "$REFD"getB.awk > BP_$C$D.pdb

#
echo "$A$B      " 
echo "$C$D      " 
#
#!  Berechne Beste Ebene
#@@@@@@@@@@@@@@@@@@@@@@@
cp  BP_$A$B-$C$D.pdb input1
"$REFD"best_plane
rm input1
cp  BP_$A$B.pdb input1
"$REFD"best_plane
rm input1
cp  BP_$C$D.pdb input1
"$REFD"best_plane
rm input1

cp  BP_$A$B.pdb input1
cp  BP_$C$D.pdb input2
"$REFD"dehedral 
rm input1
rm input2
rm jjj
rm jj

cp  P1.pdb           input1
cp  BP_$A$B-$C$D.pdb input2
"$REFD"BSScherung << eop
$KK
eop
rm input1
rm input2

cp  P1.pdb      input1
cp  BP_$A$B.pdb input2
"$REFD"BSScherung << eop
$KK
eop
rm input1
rm input2

cp  P1.pdb      input1
cp  BP_$C$D.pdb input2
"$REFD"BSScherung<< eop
$KK
eop
rm input1
rm input2

cp  P2.pdb           input1
cp  BP_$A$B-$C$D.pdb input2
"$REFD"BSScherung << eop
$KK
eop
rm input1
rm input2

cp  P2.pdb      input1
cp  BP_$A$B.pdb input2
"$REFD"BSScherung << eop
$KK
eop
rm input1
rm input2

cp  P2.pdb      input1
cp  BP_$C$D.pdb input2
"$REFD"BSScherung << eop
$KK
eop
rm input1
rm input2

#
# Strand 1
#
if (($count > 0) && ($count <= $loop)) then
  cat BP_$A$B.pdb         > input1
  echo P1.pdb  >  pdb
  echo "$B"   >>  pdb
  nawk -v CMDFILE="pdb" -f "$REFD"getPP.awk >> input1
  rm pdb
  set ref = $B
@ ref = $ref + 1
  echo P1.pdb   >  pdb
  echo "$ref"   >>  pdb
  nawk -v CMDFILE="pdb" -f "$REFD"getPP.awk >> input1
  rm pdb
#
  "$REFD"PPScherung 
  rm input1
 else
  echo " # # "
  echo " # # "
endif

if (($count > 0) && ($count <= $loop)) then
  cat BP_$A$B-$C$D.pdb    > input1
  echo P1.pdb  >  pdb
  echo "$B"   >>  pdb
  nawk -v CMDFILE="pdb" -f "$REFD"getPP.awk >> input1
  rm pdb
  set ref = $B
@ ref = $ref + 1
  echo P1.pdb   >  pdb
  echo "$ref"   >>  pdb
  nawk -v CMDFILE="pdb" -f "$REFD"getPP.awk >> input1
  rm pdb
#
  "$REFD"PPScherung 
  rm input1
 else
  echo " # # "
  echo " # # "
endif
#
# Strand 2
#
if (($count > 0) && ($count <= $loop)) then
  cat BP_$C$D.pdb         > input1
  echo P2.pdb  >  pdb
  echo "$D"   >>  pdb
  nawk -v CMDFILE="pdb" -f "$REFD"getPP.awk >> input1
  rm pdb
  set ref = $D
@ ref = $ref + 1
  echo P2.pdb   >  pdb
  echo "$ref"   >>  pdb
  nawk -v CMDFILE="pdb" -f "$REFD"getPP.awk >> input1
  rm pdb
#
  "$REFD"PPScherung 
  rm input1
 else
  echo " # # "
  echo " # # "
endif
#
if (($count > 0) && ($count <= $loop)) then
  cat BP_$A$B-$C$D.pdb    > input1
  echo P2.pdb  >  pdb
  echo "$D"   >>  pdb
  nawk -v CMDFILE="pdb" -f "$REFD"getPP.awk >> input1
  rm pdb
  set ref = $D
@ ref = $ref + 1
  echo P2.pdb   >  pdb
  echo "$ref"   >>  pdb
  nawk -v CMDFILE="pdb" -f "$REFD"getPP.awk >> input1
  rm pdb
#
  "$REFD"PPScherung 
  rm input1
 else
  echo " # # "
  echo " # # "
endif
#
if ($count < $loop ) goto Begin
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
rm in.pdb
#
rm bp
exit
#

