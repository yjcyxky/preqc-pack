# [Per Sequence Quality Score](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html)

## Summary

The Per Sequence Quality Score report describes  allows you to see if a subset of your sequences have universally low quality values. It is often the case that a subset of sequences will have universally poor quality, often because they are poorly imaged (on the edge of the field of view etc), however these should represent only a small percentage of the total sequences.

If a significant proportion of the sequences in a run have overall low quality then this could indicate some kind of systematic problem - possibly with just part of the run (for example one end of a flowcell).

Results from this module will not be displayed if your input is a BAM/SAM file in which quality scores have not been recorded.

+ **X- category quality**：It's a set of distinct quality scores in all sequences processed, and all scores are sorted in an ascending order.
+ **Y-category count**：It records the quality count in all sequences processed for each quality score in x-category quality.
+ **Most frequent score**：The quality score that occurs the most

## Example

```
"per_seq_quality_score": {
            "x_category_quality": [16,17,18,19,20,21,22,23,24,25,...],
            "y_category_count": [182111,5,4,9,37,61,126,342,479,485,...],
            "most_frequent_score": 16
        },
```

