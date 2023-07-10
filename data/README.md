## How to prepare the pattern file

```bash
cd data
python patterns2json.py
```

After the patterns2json.py, tow files would be generated, patterns.bson and patterns.json.

## Patterns Format

```
         Kmer                index    ref_or_alt(0 = ref, 1 = alt)
"TCCTTGTCATATGTTTTTCTG": [     0,         0       ]
```
