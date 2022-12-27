#!/usr/bin/bash
for file in `ls /home/anubratadas/Documents/bioinformatics_masters_uob/Data_analytics/project_proposal/source`;
do
  if [ -f $file ] && [[ $file == *.txt ]]
  then
 	  echo $file
	  cp $file /home/anubratadas/Documents/bioinformatics_masters_uob/Data_analytics/project_proposal/dump

  fi
done

echo `ls /home/anubratadas/Documents/bioinformatics_masters_uob/Data_analytics/project_proposal/dump| wc -l` files in folder
