# [ Per Base Sequence Content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html)

## Summary

Per Base Sequence Content module shows the proportion of each base position in a file for which each of the four normal DNA bases has been called.

In a random library you would expect that there would be little to no difference between the different bases of a sequence run, so the lines in this plot should run parallel with each other. The relative amount of each base should reflect the overall amount of these bases in your genome, but in any case they should not be hugely imbalanced from each other.

It's worth noting that some types of library will always produce biased sequence composition, normally at the start of the read. Libraries produced by priming using random hexamers (including nearly all RNA-Seq libraries) and those which were fragmented using transposases inherit an intrinsic bias in the positions at which reads start. This bias does not concern an absolute sequence, but instead provides enrichement of a number of different K-mers at the 5' end of the reads. Whilst this is a true technical bias, it isn't something which can be corrected by trimming and in most cases doesn't seem to adversely affect the downstream analysis. It will however produce a warning or error in this module.

+ **X-category**：It divides all base positions into several groups  which contain a certain digit of position or the range of positions. 
+ **G counts**：The counts of each base position in a file for which ’G‘ base has been called.
+ **C counts**：The counts of each base position in a file for which ’C‘ base has been called.
+ **A counts**：The counts of each base position in a file for which ’A‘ base has been called.
+ **T counts**：The counts of each base position in a file for which ’T‘ base has been called.
+ **Percentages**：The proportion of each base position in a file for which each of the four normal DNA bases has been called.Its order in  storage is `[t_percent_array, c_percent_array, a_percent_array,  g_percent_array]`

## Example

```
"per_base_seq_content": {
            "x_category": ["1","2","3","4","5","6","7","8","9","10-14","15-19","20-24","25-29",...],
            "g_counts": [
            	115472,
                80383,
                76362,
                83320,
                81674,
                80780,
                61621,
                59846,
                61399,
                60254,
                67033,
                64991,
                65134,
                ...
            ],
            "c_counts": [
            	63849,
                69603,
                80458,
                75224,
                60497,
                51506,
                49210,
                57837,
                55266,
                54779,
                59960,
                60946,
                61423,
                ...
            ],
            "a_counts": [
            	38933,
                29960,
                35063,
                41477,
                50424,
                59984,
                60756,
                62327,
                58272,
                72997,
                65287,
                64525,
                64301,
                ...
            ],
            "t_counts": [
            	31663,
                69593,
                58117,
                49979,
                57405,
                57728,
                78413,
                69990,
                75063,
                61970,
                57718,
                59538,
                59142,
                ...
            ],
            "percentages": [
                [
                    12.669406242872633,
                    27.88862662750111,
                    23.2468,
                    19.991600000000002,
                    22.962,
                    23.09138473107785,
                    31.365199999999998,
                    27.996,
                    30.0252,
                    23.834518135229015,
                    22.720000000000002,
                    23.19664,
                    22.7348,
                    ...
                ],
                [
                    25.54808196321179,
                    27.892634017127584,
                    32.1832,
                    30.0896,
                    24.198800000000002,
                    20.602564820518566,
                    19.683999999999997,
                    23.1348,
                    22.1064,
                    23.88315821305314,
                    25.505119999999998,
                    25.90864,
                    27.118,
                    ...
                ],
                [
                    15.578372019510478,
                    12.006139320907753,
                    14.025199999999998,
                    16.5908,
                    20.1696,
                    23.993791950335602,
                    24.3024,
                    24.9308,
                    23.308799999999998,
                    26.43276229241967,
                    24.285680000000003,
                    24.44032,
                    25.01016,
                    ...
                ],
                [
                    46.2041397744051,
                    32.21260003446355,
                    30.5448,
                    33.328,
                    32.669599999999996,
                    32.31225849806798,
                    24.648400000000002,
                    23.9384,
                    24.5596,
                    25.849561359298175,
                    27.489200000000004,
                    26.4544,
                    25.13704,
                    ...
                ]
            ]
  },
```

