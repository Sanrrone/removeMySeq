# removeMySeq
removeMySeq allow you to remove a subsequence from your sequence, it means giving a protein sequence (or nucl seq), the tool will search for the corresponding match on your fasta target and then remove the match. It will automatically recognize between DNA and AA for your searches in both sides.

## How it works

removeMySeq receive one sequence as input and a multiple sequences as target (multifasta), it uses Blast+ suite to match the exact position of your sequences, then the sequence is removed from the target sequence merging the remaments (by default: '').

## Requisites

* Python >= v3.7
* BioPython module == v1.7 (v1.8 or v1.6 doesn't work)
* Blast+ >= 2.10.0+

## Usage

#### minimal use

`removems -q myseq.fasta -s myfasta.fasta`

* *-q* input file, must be a single sequence in fasta format
* *-s* subject file, a multifasta file

#### complete use

`removems -q myseq.fasta -s myfasta.fasta -g aaaa -i 85 -a 85 -l 50 -v -t`

* **-q** input file, must be a single sequence in fasta format.
* **-s** subject file, a multifasta file.
* **-o** output file name (default 'clean.fasta'.
* **-g** glue, string to merge the sequences after the removed ones (default: '').
* **-i** identity, identity value for query subject that must match with subject sequences (default: 85).
* **-a** aligment length, minimal aligment length percent for query (default: 85). 
* **-l** minimum length, the minimum length to keep the sequence after removing the match part (default: 0).
* **-v** inverse the match to keep the query sequence or equivalent in subject sequences (default: False).
* **-t** translate the output sequence if they are DNA (default: False)

## Output

##### query sequence
>\>query seq
AAAAAAAA

##### target sequences

>\>MH708667.1
CATTATCGTAAAAAAAATGGACTGATTGTAGTTGTGTTGGTGGTGCT

\>MH708666.1
CATTATCGTTGGACTGATTGTAGTTGTGTAAAAAAAATGGTGGTGCT

\>MH708665.1
CATTATCGTTGGACTGATTGTAGTAAAAAAAATGTGTTGGTGGTGCT

#### Output
>\>MH708667.1
CATTATCGTTGGACTGATTGTAGTTGTGTTGGTGGTGCT

\>MH708666.1
CATTATCGTTGGACTGATTGTAGTTGTGTTGGTGGTGCT

\>MH708665.1
CATTATCGTTGGACTGATTGTAGTTGTGTTGGTGGTGCT

## Notes
the query sequence must have only one sequence, if you provide a multifasta as input only first sequence will work.
