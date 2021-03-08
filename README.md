# tbhgvs

tbhgvs contains functions for converting mutations to from genome coordinates to hgvs format.

## Install
You can install the package using the following commands:
```
git clone https://github.com/jodyphelan/tbhgvs.git
cd tbhgvs
python setup.py install
```

The next thing to do is to download the fasta and gff files (this only needs to be done once!):
```
tbhgvs-download-files.py
```

Next we can create the reference database class and use it to convert.
``` python
ref_db = tbhgvs.reference_db()
ref_db.hgvs2genome("p.Ser450Leu","Rv0667")

ref_db.genome2hgvs([{'pos': '761155', 'ref': 'C', 'alt': 'T'}, {'pos': '761156', 'ref': 'G', 'alt': 'A'}])
```
