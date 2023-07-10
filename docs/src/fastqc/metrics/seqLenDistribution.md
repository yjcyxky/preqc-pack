#  [ Sequence Length Distribution](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html)

## Summary

Some high throughput sequencers generate sequence fragments of uniform length, but others can contain reads of wildly varying lengths. Even within uniform length libraries some pipelines will trim sequences to remove poor quality base calls from the end.

This module shows the distribution of fragment sizes in the file which was analysed.

In many cases this will produce a simple graph showing a peak only at one size, but for variable length FastQ files this will show the relative amounts of each different size of sequence fragment.

+ **X-category**：A set of length.
+ **Graph counts**：The number of sentences for each length in x-category.

## Example

```
"seq_len_distribution": {
            "x_categories": [
                "149",
                "150",
                "151"
            ],
            "graph_counts": [
                0,
                250000,
                0
            ]
        },
```

