# Consistency script 
## Tool description
In a first step, the tool utilizes [DeSignate](https://designate.dbresearch.uni-salzburg.at/home/) to detect signature characters for a selected query group in a reference alignment and alternative alignments comprising identical sequences. Secondly, consensus signature characters congruently detected in all alignments are identified.

For more details and an example application, please read our manuscript @ [MPE](https://doi.org/10.1016/j.ympev.2022.107433) 

## Usage
### Requirements
This script requires DeSignate ([Hütter et al. 2020](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3498-6)), a tool that detects molecular signature characters for taxon diagnoses. To use DeSignate, clone its repository to the root directory of this repository:
```
git clone https://github.com/DatabaseGroup/DeSignate
```
### Input files
1. Alignment files in *fasta* format
   - Example:
   ```
   >Sequence-1-label
   -TTGGCTGTCACAGTGTC-
   >Sequence-2-label
   --TGGTACTGACAGTGT--
   ...
   ```
2. Two separate files with comma separated sequence labels comprising the query and reference group (e.g., in *txt* or *csv* format)
   - Example:
   ```
   Sequence-1-label, Sequence-2-label, ...
   ```
   **PLEASE NOTE:** Sequence labels must be identical in the alignments and also exactly match those in the query and reference group files.
   Otherwhise, the program terminates with an error message stating the missing/wrong sequence labels.

### Output files
- **consensus-sigchars.csv** : Alignment positions of consensus signature characters + DeSignate results (character states, signature type, entropy values)  
- **non-consensus-positions.csv** : Reference alignment positions of non-consensus signature characters
- **designate-results.csv** : Complete DeSignate results of reference alignment for the selected query and reference groups 

### Commands
To execute the script use the following command:
```
python consistency.py --alignments path/alignment_01.fasta path/alignment_02.fasta path/alignment_03.fasta --query_group path/query_group.txt --reference_group path/reference_group.txt
```
```
List of commands:
--alignments : Paths to alignment files. The first file represents the reference alignment, subsequent files represent alternative alignments.
--query_group : Path to query group file.
--reference_group: Path to reference group file.
--k_window : Two position analysis, default = 1 for one position analysis.
--consider_gaps : Include gaps as a character state, default = True.
```
