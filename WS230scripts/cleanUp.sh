#! /bin/bash
#last modified 20151202

myWD=`pwd ../`
parentDir="$(dirname "$myWD")"
echo 'Are you sure you want to clean up' $parentDir 'directory (y/n)? ' 
read answer
case ${answer:0:1} in
    y|Y )
    rm -r ../samFiles/
	rm -r ../tempData/
	rm `ls -l ../bamFiles/* | grep -v .sorted.bam | awk  '{print $9}'`
	echo samFiles tempData directories and unsorted bam files deleted
    ;;
    * )
        echo Nothing deleted
    ;;
esac


