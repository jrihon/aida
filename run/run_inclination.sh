#!/bin/csh -f
#
#
########################
# INPUT 
#
set PDB="$1"          # PDB File Name
set P1="$2"           # Anzahl P Atomen im Strang 1
set BK="$3"           # Backbone atoms 
#
#
#! create a file called bp with base pair information
#! one on each line like this:  A 1  T 16 ...
#
cp "$1".bp bp
#
set REFD=`pwd`/../bin #set bin/ directory
set KK="3"            # BSline Konstante
#

#set REFD="bin"
#
#!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#!!!!!!                DO NOT TOUCH BELOW                   !!!!!!
#!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
cat << eof > jobscherung
#!/bin/csh -f
set PDB=$PDB
set P1=$P1
set BK="$BK"
set KK=$KK
set REFD=$REFD/
eof
#
cat $REFD/scherung.sh >> jobscherung
#
chmod +x jobscherung
#
./jobscherung > q
#
#
cat << eof > jobanalyse
#!/bin/csh -f
set PDB=$PDB
set BK="$BK"
set KK=$KK
set REFD=$REFD/
eof
#
cat $REFD/analyse.sh >> jobanalyse
#
chmod +x jobanalyse
#
./jobanalyse 
#
rm job*
rm q
rm ppp
rm *.pdb
#
#
mv results_$PDB ../output/
echo "created results_$PDB in output/"
#cat output/results_$PDB
#
#rm $1
#rm $1.bp

