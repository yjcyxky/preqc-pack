#  [Kmer Content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html)

## Summary

The analysis of overrepresented sequences will spot an increase in any exactly duplicated sequences, but there are a different subset of problems where it will not work.

- If you have very long sequences with poor sequence quality then random sequencing errors will dramatically reduce the counts for exactly duplicated sequences.
- If you have a partial sequence which is appearing at a variety of places within your sequence then this won't be seen either by the per base content module or the duplicate sequence analysis.

The Kmer module starts from the assumption that any small fragment of sequence should not have a positional bias in its apearance within a diverse library. There may be biological reasons why certain Kmers are enriched or depleted overall, but these biases should affect all positions within a sequence equally. This module therefore measures the number of each 7-mer at each position in your library and then uses a binomial test to look for significant deviations from an even coverage at all positions. Any Kmers with positionally biased enrichment are reported. The top 6 most biased Kmer are additionally reported to show their distribution.

To allow this module to run in a reasonable time only 2% of the whole library is analysed and the results are extrapolated to the rest of the library. Sequences longer than 500bp are truncated to 500bp for this analysis.

+ **Enriched Kmers**： This is the full set of Kmers to be reported.
+ **Enrichments**：Significant deviations from an even coverage at all positions for each kmer in x-labels.
+ **X-categories**：It divides all positions into several groups  which contain a certain digit of position or the range of positions. 
+ **X-labels**：A set of the top 6 most biased Kmer are additionally plotted to show their distribution.

## Example

```
        "kmer_content": {
            "enriched_kmers": [
                {
                    "sequence": "CGCATTT",
                    "count": 7,
                    "lowest_pvalue": 0.0036627573317673523
                },
                {
                    "sequence": "CTCGCTA",
                    "count": 56,
                    "lowest_pvalue": 0
                },
                {
                    "sequence": "CCCCTAT",
                    "count": 14,
                    "lowest_pvalue": 0.0010821496325661428
                },
            	...
            ],
            "enrichments": [
                [
                    0,
                    0,
                    61.712828571428574,
                    0,
                    0,
                    0,
                    0,
                    ...
                ],
                [
                    43.792079314194126,
                    23.179397750686814,
                    10.28547142857143,
                    7.714103571428572,
                    7.714103571428572,
                    7.714103571428572,
                    ...
                ],
                [
                    41.2160746486533,
                    10.301954555860807,
                    0,
                    0,
                    0,
                    0,
                    ...
                ],
                [
                    0,
                    0,
                    20.57094285714286,
                    0,
                    30.856414285714287,
                    41.14188571428572,
                    ...
                ],
                [
                    18.032032658785813,
                    40.56394606370192,
                    20.249521875,
                    8.9997875,
                    6.749840625,
                    6.749840625,
                    ...
                ],
                [
                    6.272011359577675,
                    16.72201319212189,
                    39.65123768115942,
                    18.782165217391306,
                    8.347628985507248,
                    6.260721739130435,
                    ...
                ]
            ],
            "x_categories": [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10-14",
                "15-19",
                "20-24",
                ...
            ],
            "x_labels": [
                "CGCATTT",
                "CTCGCTA",
                "CCCCTAT",
                "CGATTTG",
                "TCGCTAT",
                "CGCTATG"
            ]
        },
```

