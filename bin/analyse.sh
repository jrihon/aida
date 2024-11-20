#
#
########################
# setting up the output ...
#
echo "*** S C H E R U N G ***" > results_$PDB
echo " " >> results_$PDB
echo " " >> results_$PDB
echo " " >> results_$PDB
echo `date` >> results_$PDB
echo " " >> results_$PDB
echo "INPUT File        :  $PDB " >> results_$PDB
echo " " >> results_$PDB
echo " " >> results_$PDB
echo "Backbone Atome    :  $BK " >> results_$PDB
echo " " >> results_$PDB
echo " " >> results_$PDB
echo "BSpline Konstante :  $KK " >> results_$PDB
echo " " >> results_$PDB
echo " " >> results_$PDB
#
cp  q input1
#
#echo "$REFD evaluate >> results_$PDB"
"$REFD"evaluate >> results_$PDB
rm input1
#
echo " " >> results_$PDB
echo " " >> results_$PDB
#
cat ppp  >> results_$PDB
#

