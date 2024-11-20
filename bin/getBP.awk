#
# getBP.awk Get Bassenpaar Atome
#
BEGIN { 
    getline < CMDFILE
    file1=$0 
    getline < CMDFILE
    A=$0
    getline < CMDFILE
    B=$0
    getline < CMDFILE
    C=$0
    getline < CMDFILE
    D=$0
	readmol(file1,A,B,C,D)			# read file 1
}
#
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
#
function readmol(m,A,B,C,D){
	while(getline < m > 0 ) {
		if ($1 == "ATOM"){
			if ((substr($4,1,1) == A) && ($5 == B)) { 
			     if ( A ~ "G" ) { guanin($0)  } 
			     if ( A ~ "C" ) { cytosin($0) } 
			     if ( A ~ "A" ) { adenin($0)  } 
			     if ( A ~ "T" ) { thymin($0)  } 
			     if ( A ~ "U" ) { uracil($0)  }
			     if ( A ~ "I" ) { inosin($0)  }
					}
			if ((substr($4,1,1) == C) && ($5 == D)) { 
			     if ( C ~ "G" ) { guanin($0)  } 
			     if ( C ~ "C" ) { cytosin($0) } 
			     if ( C ~ "A" ) { adenin($0)  } 
			     if ( C ~ "T" ) { thymin($0)  } 
			     if ( C ~ "U" ) { uracil($0)  }
			     if ( C ~ "I" ) { inosin($0)  }
                    }
    			}
		}
}
function guanin (a){
	    if (($3 == "C1'") || ($3 == "C1*"))
	      { putatom(a) }
	    if ($3 == "N9") 	
	      { putatom(a) }    
	    if ($3 == "C4") 	
	      { putatom(a) }    
	    if ($3 == "N3") 	
	      { putatom(a) }    
	    if ($3 == "C2") 	
	      { putatom(a) }    
 	    if ($3 == "N2") 	
	      { putatom(a) }    
	    if ($3 == "N1") 	
	      { putatom(a) }    
	    if ($3 == "C6") 	
	      { putatom(a) }    
 	    if ($3 == "O6") 	
	      { putatom(a) }    
	    if ($3 == "C5") 	
	      { putatom(a) }    
 	    if ($3 == "N7") 	
	      { putatom(a) }    
 	    if ($3 == "C8") 	
	      { putatom(a) }    
}
#
function adenin (a){
	    if (($3 == "C1'") || ($3 == "C1*"))
	      { putatom(a) }
	    if ($3 == "N9") 	
	      { putatom(a) }    
	    if ($3 == "C4") 	
	      { putatom(a) }    
	    if ($3 == "N3") 	
	      { putatom(a) }    
	    if ($3 == "C2") 	
	      { putatom(a) }    
	    if ($3 == "N1") 	
	      { putatom(a) }    
	    if ($3 == "C6") 	
	      { putatom(a) }    
	    if ($3 == "N6") 	
	      { putatom(a) }    
	    if ($3 == "C5") 	
	      { putatom(a) }    
	    if ($3 == "N7") 	
	      { putatom(a) }    
	    if ($3 == "C8") 	
	      { putatom(a) }    
}
#
function cytosin (a){
	    if (($3 == "C1'") || ($3 == "C1*"))
	      { putatom(a) }
	    if ($3 == "N1") 	
	      { putatom(a) }    
	    if ($3 == "C2") 	
	      { putatom(a) }    
	    if ($3 == "O2") 	
	      { putatom(a) }    
	    if ($3 == "N3") 	
	      { putatom(a) }    
	    if ($3 == "C4") 	
	      { putatom(a) }    
	    if ($3 == "N4") 	
	      { putatom(a) }    
	    if ($3 == "C5") 	
	      { putatom(a) }    
	    if ($3 == "C6") 	
	      { putatom(a) }    
}
#
function thymin (a){
	    if (($3 == "C1'") || ($3 == "C1*"))
	      { putatom(a) }
	    if ($3 == "N1") 	
	      { putatom(a) }    
	    if ($3 == "C2") 	
	      { putatom(a) }    
	    if ($3 == "O2") 	
	      { putatom(a) }    
	    if ($3 == "N3") 	
	      { putatom(a) }    
	    if ($3 == "C4") 	
	      { putatom(a) }    
	    if ($3 == "O4") 	
	      { putatom(a) }    
	    if ($3 == "C5") 	
	      { putatom(a) }
	    if ($3 == "C5M") 	
	      { putatom(a) }
	    if ($3 == "C5A") 	
	      { putatom(a) }
	    if ($3 == "C6") 	
	      { putatom(a) }    
}
#
function uracil (a){
	    if (($3 == "C1'") || ($3 == "C1*"))
	      { putatom(a) }
	    if ($3 == "N1") 	
	      { putatom(a) }    
	    if ($3 == "C2") 	
	      { putatom(a) }    
	    if ($3 == "O2") 	
	      { putatom(a) }    
	    if ($3 == "N3") 	
	      { putatom(a) }    
	    if ($3 == "C4") 	
	      { putatom(a) }    
	    if ($3 == "O4") 	
	      { putatom(a) }    
	    if ($3 == "C5") 	
	    if ($3 == "C6") 	
	      { putatom(a) }    
}
function inosin (a){
	    if (($3 == "C1'") || ($3 == "C1*"))
	      { putatom(a) }
	    if ($3 == "N9") 	
	      { putatom(a) }    
	    if ($3 == "C4") 	
	      { putatom(a) }    
	    if ($3 == "N3") 	
	      { putatom(a) }    
	    if ($3 == "C2") 	
	      { putatom(a) }    
	    if ($3 == "N1") 	
	      { putatom(a) }    
	    if ($3 == "C6") 	
	      { putatom(a) }    
	    if ($3 == "O6") 	
	      { putatom(a) }    
	    if ($3 == "C5") 	
	      { putatom(a) }    
	    if ($3 == "N7") 	
	      { putatom(a) }    
	    if ($3 == "C8") 	
	      { putatom(a) }    
}

