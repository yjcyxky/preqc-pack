#!/bin/bash

printf 'filesize\tgz_filesize\tnrows' > estimate.tsv

function stat() {
    filesize=`du -s $1 | cut -f 1`
    wcount=`zcat $1 | wc -cl`
    fq_filesize=`echo $wcount | cut -d ' ' -f 2`
    nrow=`echo $wcount | cut -d ' ' -f 1`

    printf "$filesize\t$fq_filesize\t$nrow" >> estimate.tsv
}

ls *.fq.gz > files.txt

while read line
do
    stat $line
done < files.txt