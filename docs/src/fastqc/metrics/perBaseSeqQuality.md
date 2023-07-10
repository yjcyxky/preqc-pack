# [Per Base Sequence Quality](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html)

## Summary

The Per Base Sequence Quality module  shows an overview of the range of quality values across all bases at each position in the FastQ file.
It divides all positions into several groups  which contain a certain digit of position or the range of positions. 
+ **X-labels**：It divides all positions into several groups  which contain a certain digit of position or the range of positions. 

- **Mean**：For each position or range of position in x-labels, It calculates the average quality scores  for all sequences.
- **Median**：For each position or range of position in x-labels, It calculates the median quality scores  for all sequences.
- **Lower Quartile**：For each position or range of position in x-labels, It calculates the  lower quartile（25%） quality scores  for all sequences.
- **Upper Quartile**：For each position or range of position in x-labels, It calculates the higher quartile（75%） quality scores  for all sequences.
- **Lowest**：For each position or range of position in x-labels, It calculates the lowest quality scores  for all sequences.
- **Highest**：For each position or range of position in x-labels, It calculates the highest quality scores  for all sequences.

## Example

```
"per_base_seq_quality": {
            "xlabels": ["1","2","3","4","5","6","7","8","9","10-14","15-19","20-24","25-29","30-34","35-39",...],
            "mean": [36.307356,36.290732,36.46844,36.51792,36.521944,36.526816,36.473336,
36.543864,36.515792,36.563705600000006,36.560803199999995,36.5081312,36.4498688,36.42268,36.3846096, ...],
            "median": [37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,...],
            "lower_quartile": [37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,... ],
            "upper_quartile": [37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,... ],
            "lowest": [37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,... ],
            "highest": [37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,... ],
        },
```

