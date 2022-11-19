# [ Per Sequence GC Content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html)

## Summary

This module measures the GC content across the whole length of each sequence in a file and compares it to a modelled normal distribution of GC content.

In a normal random library you would expect to see a roughly normal distribution of GC content where the central peak corresponds to the overall GC content of the underlying genome. Since we don't know the the GC content of the genome the modal GC content is calculated from the observed data and used to build a reference distribution.

An unusually shaped distribution could indicate a contaminated library or some other kinds of biased subset. A normal distribution which is shifted indicates some systematic bias which is independent of base position. If there is a systematic bias which creates a shifted normal distribution then this won't be flagged as an error by the module since it doesn't know what your genome's GC content should be.

+ **X-category**：It's a set of mean %GC. 
+ **Y-gc distribution**：For each %GC in x-category，we count the number of corresponding sentences.
+ **Y-theoretic distribution**：  The theoretic number of corresponding sentences for each %GC in x-category.
+ **Deviation percent**：Records the sum of the deviations for all corresponding positions of y-gc distribution and y-theoretic distribution.

## Example

```
"per_seq_gc_content": {
            "x_category": [
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                ...
            ],
            "y_gc_distribution": [
                0,
                0,
                0.5,
                0.5,
                0,
                0,
                0,
                1,
                1,
                3.5,
                3.5,
                0,
                3,
                8,
               ...
            ],
            "y_theo_distribution": [
                11.301763218607428,
                13.897288489846439,
                17.029620470668664,
                20.79557306799733,
                25.306252897428912,
                30.68851136823189,
                37.08641280899335,
                44.66268842218363,
                53.60013975983366,
                64.10294913058424,
                76.39784805677265,
                90.73508881215378,
                107.3891584351843,
               ...
            ],
            "deviation_percent": 42.145681622216166
        },
```

