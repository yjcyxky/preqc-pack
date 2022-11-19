# [ Per Base N Content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html)

## Summary

If a sequencer is unable to make a base call with sufficient confidence then it will normally substitute an N rather than a conventional base] call.

It's not unusual to see a very low proportion of Ns appearing in a sequence, especially nearer the end of a sequence. However, if this proportion rises above a few percent it suggests that the analysis pipeline was unable to interpret the data well enough to make valid base calls.

+ **X-category**：It divides all base positions into several groups  which contains a certain digit of position or the range of positions. 
+ **N counts**：Records the count of base 'N' for each position in x-category.
+ **Not n counts**：  Records the sum of count of base 'A', base 'T', base 'C' and base 'G' for each position in x-category.
+ **Percentages**：The proportion of each base position in x-category for which the base 'N' has been called.

## Example

```
"per_base_n_content": {
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
            "n_counts": [
                83,
                461,
                0,
                0,
                0,
                2,
                0,
                0,
                0,
                0,
                2,
                0,
                ...
            ],
            "not_n_counts": [
                249917,
                249539,
                250000,
                250000,
                250000,
                249998,
                250000,
                250000,
                250000,
                250000,
                249998,
                250000,
                ...
            ],
            "percentages": [
                0.0332,
                0.18439999999999998,
                0,
                0,
                0,
                0.0007999999999999999,
                0,
                0,
                0,
                0.00015999999999999999,
                0,
                0,
                ...
            ] 
        },
			
```

