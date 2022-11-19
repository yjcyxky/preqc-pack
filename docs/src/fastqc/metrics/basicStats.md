# [Basic Statistics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/1%20Basic%20Statistics.html)

## Summary

The Basic Statistics module generates some simple composition statistics for the file analysed.

- **File Name**：The original filename of the file which was analysed
- **File type**：Says whether the file appeared to contain actual base calls or colorspace data which had to be converted to base calls
- **Phred Encoding**：Says which ASCII encoding of quality values was found in this file.
- **Total Reads**：A count of the total number of sequences processed. 
- **Total Bases**：A count of the total number of bases in all sequences processed. 
- **Total T Bases**：A count of the total number of base 'T' in all sequences processed. 
- **Total C Bases**：A count of the total number of base 'C' in all sequences processed. 
- **Total G Bases**：A count of the total number of base 'G' in all sequences processed. 
- **Total A Bases**：A count of the total number of base 'A' in all sequences processed. 
- **Total N Bases**： A count of the total number of base 'N' in all sequences processed , where 'N' is the base that cannot be recognized by sequencing. 
- **%GC**：The overall %GC of all bases in all sequences
- **Min Length**：The length of the shortest sequence in the set.
- **Max Length**：The length of the highest sequence in the set.
- **Lowest Char:** The lowest quality char in all sequences processed. 
- **Highest Char:** The highest quality char in all sequences processed. 

## Example

```
"basic_stats"：{
            "file_name"："test.fastq.gz",
            "file_type"："",
            "phred"：{
                "name"："Sanger / Illumina 1.9",
                "offset"：33
            },
            "total_reads"：250000,
            "total_bases"：37500000,
            "t_count"：8473918,
            "c_count"：10050719,
            "g_count"：9888266,
            "a_count"：9086549,
            "n_count"：548,
            "gc_percentage"：0.5317062666666666,
            "min_length"：150,
            "max_length"：150,
            "lowest_char"：35,
            "highest_char"：70
        },
```

