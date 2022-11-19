#  [ Adapter Content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html)

## Summary

The Kmer Content module will do a generic analysis of all of the Kmers in your library to find those which do not have even coverage through the length of your reads. This can find a number of different sources of bias in the library which can include the presence of read-through adapter sequences building up on the end of your sequences.

You can however find that the presence of any overrepresented sequences in your library (such as adapter dimers) will cause the Kmer plot to be dominated by the Kmers these sequences contain, and that it's not always easy to see if there are other biases present in which you might be interested.

One obvious class of sequences which you might want to analyse are adapter sequences. It is useful to know if your library contains a significant amount of adapter in order to be able to assess whether you need to adapter trim or not. Although the Kmer analysis can theoretically spot this kind of contamination it isn't always clear. This module therefore does a specific search for a set of separately defined Kmers and will give you a view of the total proportion of your library which contain these Kmers. A results trace will always be generated for all of the sequences present in the adapter config file so you can see the adapter content of your library, even if it's low.

Once a sequence has been seen in a read it is counted as being present right through to the end of the read so the percentages you see will only increase as the read length goes on

The module also calculates an expected overall loss of sequence were the library to be deduplicated.

+ **Labels**：A set of adpaters from the file of adpater.
+ **X-labels**：It divides all positions into several groups  which contain a certain digit of position or the range of positions. 
+ **Enrichments**： It ranges in the order of adapters in `Labels` and shows a cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.

## Example

```
"adapter_content": {
            "labels": [
                "Illumina Universal Adapter",
                "Illumina Small RNA 3' Adapter",
                "Illumina Small RNA 5' Adapter",
                "Nextera Transposase Sequence",
                "SOLID Small RNA Adapter"
            ],
            "x_labels": [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10-11",
                "12-13",
                "14-15",
                "16-17",
                ...
            ],
            "enrichments": [
                [
                    0.002,
                    0.002,
                    0.002,
                    0.002,
                    0.002,
                    0.002,
                    0.002,
                    0.002,
                    0.002,
                    0.0021999999999999997,
                    0.0024,
                    0.0024,
                    0.0024,
                    ...
                ],
                [
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    ...
                ],
                [
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    ...
                ],
                [
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    ...
                ],
                [
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    ...
                ]
            ]
        },
```

