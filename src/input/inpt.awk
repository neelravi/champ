#!/usr/bin/gawk -f 
# look for command spec. 
BEGIN{
 if(p2file=="") p2file="auto.p2";
 flg=0;
 csub="";
 ff="";
 printf "# generated by inpt.awk -- do not edit\n" > p2file;
}
FILENAME != ff {
  ff=FILENAME; echof=1;
  if(logf) printf "file %s \n",ff >> logf;
}
/^      [ ]*subroutine/ || /^      [ ]*SUBROUTINE/{
         if((ix=index($2,"("))) csub=substr($2,1,ix-1); else csub=$2;
#	 if(logf) printf "subroutine %s\n" , csub >> logf;
         next;
}
/^C\$INPUT/ { 
  if(echof){
    printf "# %s\n",FILENAME >> p2file; 
    echof=0; }
  if(NF>1){ 
    printf "%s ",$2 >> p2file;
    printf "%s ",csub >> p2file;
    if(NF>2)
      for(i=3;i<=NF;i++)
	printf "%s ",$i >> p2file; 
    printf "\n" >> p2file;
  } else {
    printf "inpt.awk: %s : %d : keyword missing \n",FILENAME,FNR > "/dev/stderr";
  }
  if(logf) printf "key= %s  routine=%s\n",$2,csub >> logf;
}

 