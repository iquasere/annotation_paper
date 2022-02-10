#! /bin/bash
IDS=$1
UNIPROT=$2
OUTFILE=$3

awk 'BEGIN {
  i=0;
  while (getline < "'"$IDS"'"){
    i+=1;
    ids[$0]=i;
    }
    close("'"$IDS"'");
  }
  {
    if ($0 in ids){
	    print $0" in IDs!"
      getline;
      while ($0 !~ /^>/){
        getline;
      }
	    if (!($0 in ids)){
		    print $0 > "'"$OUTFILE"'";
	    }
    } else {
      print $0 > "'"$OUTFILE"'";
    }
  }' $UNIPROT
