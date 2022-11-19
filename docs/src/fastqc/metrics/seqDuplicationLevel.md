#  [ Duplicate Sequence](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)

## Summary

In a diverse library most sequences will occur only once in the final set. A low level of duplication may indicate a very high level of coverage of the target sequence, but a high level of duplication is more likely to indicate some kind of enrichment bias (eg PCR over amplification).

This module counts the degree of duplication for every sequence in a library and shows the relative number of sequences with different degrees of duplication.

To cut down on the memory requirements for this module only sequences which first appear in the first 100,000 sequences in each file are analysed, but this should be enough to get a good impression for the duplication levels in the whole file. Each sequence is tracked to the end of the file to give a representative count of the overall duplication level. To cut down on the amount of information any sequences with more than 10 duplicates are placed into grouped bins to give a clear impression of the overall duplication level without having to show each individual duplication value.

Because the duplication detection requires an exact sequence match over the whole length of the sequence, any reads over 75bp in length are truncated to 50bp for the purposes of this analysis. Even so, longer reads are more likely to contain sequencing errors which will artificially increase the observed diversity and will tend to underrepresent highly duplicated sequences.

The module also calculates an expected overall loss of sequence were the library to be deduplicated.

+ **Total percentages**：It takes the full sequence set and shows how its duplication levels are distributed. 
+ **Deduplicated percentages**：The sequences are de-duplicated and the proportions  are the proportions of the deduplicated set which come from different duplication levels in the original data.
+ **Percentage diff**： Calculates an expected overall loss of sequence were the library to be deduplicated.
+ **labels** : It's a set which shows different degrees of duplication.

## Example

```
"dedup_percentages": [
                88.7270625897991,
                8.234690198626954,
                1.3275601587718149,
                0.4567248391318406,
                0.23253294096467203,
                0.16711084299112353,
                0.10936917951231077,
                0.08366047968337169,
                0.060472876913210516,
                0.4965347445786828,
                0.05984315938476192,
                0.041475456999339945,
                0.0011850130571239984,
                0.0017775195856859977,
                0,
                0
            ],
            "total_percentages": [
                62.5830541835812,
                11.616569913245959,
                2.809157665879907,
                1.2885915309298914,
                0.8200779570004981,
                0.7072232509214026,
                0.5399998558899416,
                0.47207498413461896,
                0.383887339423591,
                7.222186751326263,
                2.92962348084703,
                6.041878268559988,
                0.6720163419689051,
                1.9136584762908062,
                0,
                0
            ],
            "percent_diff_seq": 70.53434697022905,
            "labels": [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                ">10",
                ">50",
                ">100",
                ">500",
                ">1k",
                ">5k",
                ">10k"
            ]
```

