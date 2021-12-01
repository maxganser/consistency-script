# Consistency script 
## Tool description
In a first step, the tool utilizes [DeSignate](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3498-6) to detect signature characters for a selected query group in a reference alignment and alternative alignments comprising identical sequences. Secondly, consensus signature characters congruently detected in all alignments are identified.

For more details and an example application, please read our manuscript (*submitted*) @ [MPE](https://www.journals.elsevier.com/molecular-phylogenetics-and-evolution) 

## Usage
### Requirements
This script requires DeSignate ([Hütter et al. 2020](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3498-6)), a tool that detects molecular signature characters for taxon diagnoses. To use DeSignate, clone its repository to the root directory of this repository:
```
git clone https://github.com/DatabaseGroup/DeSignate
```
### Input files
1. Alignment files in *fasta* format (**PLEASE NOTE** that sequence labels must be identical in all files. Otherwhise, the script might currently not detect potential consensus signature characters even if they are present)
   - Example:
   ```
   >Sequence-1-label
   -TTGGCTGTCACAGTGTC-
   >Sequence-2-label
   --TGGTACTGACAGTGT--
   ...
   ```
2. Two separate files with comma separated labels of sequences comprising the query and reference group (e.g., in *txt* or *csv* format)
   - Example:
   ```
   Sequence-1-label, Sequence-2-label, ...
   ```
### Output files
Work in progress:
- *consensus-sigchars.csv*: Alignment positions of consensus signature characters + DeSignate results (character states, signature type, entropy values)  
- *non-consensus-positions.csv*: Reference alignment positions of non-consensus signature characters
- *designate-results.csv*: Complete DeSignate results of reference alignment for the selected query and reference groups 

### Commands
To execute the script use the following command:
```
python consistency.py --alignments path/alignment_01.fasta path/alignment_02.fasta path/alignment_03.fasta --query_group path/query_group.txt --reference_group path/reference_group.txt
```
