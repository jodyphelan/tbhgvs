# tbhgvs

tbhgvs contains functions for converting mutations to from genome coordinates to hgvs format.

## Install
You can install the package using the following commands:
```
git clone https://github.com/jodyphelan/tbhgvs.git
cd tbhgvs
python setup.py install
```

The ```genome2hgvs``` function calls takes a list of dictionaries in the form of:
```
[{"pos":7585,"ref":"G","alt":"C"}]
```
and returns a dictionary with the form of:

 ```
 {'nucleotide_change': '7585G>C', 'type': 'missense', 'gene_name': 'gyrA', 'gene_id': 'Rv0006', 'hgvs': 'p.Ser95Thr'}
 ```
 If multiple mutations need to be annotated together (e.g. for mutations in the same codon), multiple dictionaries can be supplied in the list.

```
genome2hgvs([
    {'pos': '7585', 'ref': 'G', 'alt': 'C'},
    {'pos': '7586', 'ref': 'C', 'alt': 'G'}
])
```
Have a look at `example_script.py` to see how to define a command-line script to perform the conversion.
