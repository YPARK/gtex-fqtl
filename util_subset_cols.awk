#!/bin/awk -f

BEGIN{
    split(COLS, cols, ",");
    for(j=1; j<=length(cols); ++j) {
	col2pos[j] = cols[j];
	pos2col[cols[j]] = j;
    }
    IGNORECASE = 1;
    if(length(lb) == 0) lb = -10000;
    if(length(ub) == 0) ub = 10000;
}
{
    c = col2pos[1];
    v = $c;
    if($c ~ /na/ || length($c) == 0){
	v = "NA"
    } else {
	if($c <= lb) v = lb
	if($c >= ub) v = ub
    }
    printf v

    for(j=2; j<=length(cols); ++j){
	c = col2pos[j];
	v = $c;
	if($c ~ /na/){
	    v = "NA"
	} else {
	    if($c <= lb) v = lb
	    if($c >= ub) v = ub
	}
	printf FS v;
    }
    printf "\n";
}

