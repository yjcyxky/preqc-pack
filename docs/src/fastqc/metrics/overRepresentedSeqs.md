#  [ Overrepresented Sequences](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html)

## Summary

A normal high-throughput library will contain a diverse set of sequences, with no individual sequence making up a tiny fraction of the whole. Finding that a single sequence is very overrepresented in the set either means that it is highly biologically significant, or indicates that the library is contaminated, or not as diverse as you expected.

This module lists all of the sequence which make up more than 0.1% of the total. To conserve memory only sequences which appear in the first 100,000 sequences are tracked to the end of the file. It is therefore possible that a sequence which is overrepresented but doesn't appear at the start of the file for some reason could be missed by this module.

For each overrepresented sequence the program will look for matches in a database of common contaminants and will report the best hit it finds. Hits must be at least 20bp in length and have no more than 1 mismatch. Finding a hit doesn't necessarily mean that this is the source of the contamination, but may point you in the right direction. It's also worth pointing out that many adapter sequences are very similar to each other so you may get a hit reported which isn't technically correct, but which has very similar sequence to the actual match.

Because the duplication detection requires an exact sequence match over the whole length of the sequence any reads over 75bp in length are truncated to 50bp for the purposes of this analysis. Even so, longer reads are more likely to contain sequencing errors which will artificially increase the observed diversity and will tend to underrepresent highly duplicated sequences.

+ **count**：Total number of  sequences in the file.
+ **Obervation cut off**：The number of unique sequences we want to track in overrepresented module.
+ **Unique sequence count**：The number of unique sequences at present.
+ **Count at unique limit**： The corresponding `count` that has been processed when `unique sequence count`  reaches the `obervation cut off`.
+ **Overrepresented sequences** : It's a set of overrepresented  sequence with attributes "seq"，“count”，"percentage" and contaminant hit", where "seq" is stand for a sequence of bases， “count” is stand for  the number of sentences containing such base sequence, "percentage" is stand for the proportion of "count" in total count, and contaminant hit is stand for the hit results in file of comtaminant .

## Example

```
"overrepresented_seqs": {
            "count": 250000,
            "observation_cut_off": 100000,
            "unique_seq_count": 100000,
            "count_at_unique_limit": 140053,
            "overrepresented_seqs": [
                {
                    "seq": "GTGGCTATTCACAGGCGCGATCCCACTACTGATCAGCACGGGAGTTTTGA",
                    "count": 1941,
                    "percentage": 0.7764,
                    "contaminant_hit": null
                },
                {
                    "seq": "GCAGTGGCTATTCACAGGCGCGATCCCACTACTGATCAGCACGGGAGTTT",
                    "count": 1347,
                    "percentage": 0.5388000000000001,
                    "contaminant_hit": null
                },
                {
                    "seq": "AGTGGCTATTCACAGGCGCGATCCCACTACTGATCAGCACGGGAGTTTTG",
                    "count": 1291,
                    "percentage": 0.5164,
                    "contaminant_hit": null
                },
                {
                    "seq": "GTGCAGTGGCTATTCACAGGCGCGATCCCACTACTGATCAGCACGGGAGT",
                    "count": 916,
                    "percentage": 0.3664,
                    "contaminant_hit": null
                },
                ...
            ]
        },
```

