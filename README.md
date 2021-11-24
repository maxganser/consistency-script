# consistency-script


## Preliminaries
This script requires DeSignate, a tool that detects molecular signature characters for taxon diagnoses. To use DeSignate, clone its repository to the root directory of this repository:

```
git clone https://github.com/DatabaseGroup/DeSignate
```

## Usage
To execute the script use the following command:

```
python consistency.py --alignments path/alignment_01.fasta path/alignment_02.fasta path/alignment_03.fasta --query_group path/query_group.txt --reference_group path/reference_group.txt
```
