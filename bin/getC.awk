#
# getC.awk Sortiere alle BK Atome aus einem FILE mit PDB Format
#
BEGIN { 
    getline < CMDFILE
    file1=$0
    getline < CMDFILE
    BK=$0
	readmol(file1,BK)			# read file 1
}
function putatom(a){
			  ATOM    = "ATOM"
			  Nr      = substr(a,7,5)
			  Name    = substr(a,13,4)
			  ResName = substr(a,18,3)
			  ResNr   = substr(a,23,4)
			  X       = substr(a,31,8)
			  Y       = substr(a,39,8)
			  Z       = substr(a,47,8)
			  O       = substr(a,55,6)
			  T       = substr(a,61,6)
		   	     printf "%-4s%7d %-4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	                  	ATOM,Nr,Name,ResName,ResNr,X,Y,Z,O,T
}
function readmol(m,c){
	while(getline < m > 0 ) {
		if ($1 == "ATOM"){
			if ($3 == c) {
               putatom($0)
    			}
		}
    }
}

