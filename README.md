# tbhgvs

tbhgvs contains functions for converting mutations to from genome coordinates to hgvs format.

## Install
You can install the package using the following commands:
```
git clone https://github.com/jodyphelan/tbhgvs.git
cd tbhgvs
python setup.py install
```

## Usage

The fisrt thing to do is to download the fasta and gff files (this only needs to be done once!):
```
tbhgvs-download-files.py
```

You can make your own script or you can use you can use the `tbhgvs-covert.py` script as a commands-line tool.
If you would like to make your own script you can use this example as a template.

``` python
ref_db = tbhgvs.reference_db()
ref_db.hgvs2genome("p.Ser450Leu","Rv0667")

ref_db.genome2hgvs([{'pos': '761155', 'ref': 'C', 'alt': 'T'}, {'pos': '761156', 'ref': 'G', 'alt': 'A'}])
```
