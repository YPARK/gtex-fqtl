#!/usr/bin/awk -f

BEGIN{
  split(ROWS, rows, ",");
  mtot = length(rows);
  for(j=1; j<=mtot; ++j) {
    row2pos[j] = rows[j];
    pos2row[rows[j]] = j;
  }
  IGNORECASE = 1;
}
(NR in pos2row){
  row = pos2row[NR]
  row2data[row] = $0
}
END{
  for(j = 1; j<= mtot; ++j)
    print row2data[j]
}
