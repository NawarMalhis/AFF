# Annotated Fasta Format (AFF)

The Annotated FASTA Format (AFF), with a .af file extension, is a text-based format for protein amino acid sequences that includes multiple annotation tags and a statistical summary header. This AFF Python library provides tools for reading, writing, and processing .af files, as well as reading residue-level scores, generating ROC and precision-recall plots, and calculating metrics such as AUC, average precision score (APS), and success rate.

AFF files are divided into two sections, the header section and the data section. See the AFF sample file below.

### AFF sample file
```bash
# Data Name: MC2 CV1 Dataset
#
# Optional top comments
#
# Sequences: 141
#
# Format:
#   >accession|Fasta=Fasta_ID|UniProt=UniProt_ID|PDB=PDB_ID|PHI-base=PHI-base_ID|DisProt=DisProt_ID|IntAct=IntAct_ID|ELM=ELM_ID|STRING=STRING_ID|IDEAL=IDEAL_ID|SignaLink=SignaLink_ID|Pfam=Pfam_ID|AlphaFoldDB=AlphaFoldDB_ID|DrugBank=DrugBank_ID|BioGRID=BioGRID_ID|MINT=MINT_ID|DIP=DIP_ID|KEGG=KEGG_ID|BindingDB=BindingDB_ID
#   Amino acid sequence
#   MoRFTest High-quality MoRF annotations used for testing
#   MoRFTrain Low-quality MoRF annotations used for training
#
# ID Counts:
#   ID All# Unique#
#   Fasta 141 141
#   UniProt 140 140
#   PDB 124 123
#   PHI-base 2 2
#   DisProt 56 56
#   IntAct 102 102
#   ELM 40 40
#   STRING 114 114
#   IDEAL 54 54
#   SignaLink 65 65
#   Pfam 135 108
#   AlphaFoldDB 125 125
#   DrugBank 27 25
#   BioGRID 87 87
#   MINT 82 82
#   DIP 81 81
#   KEGG 122 122
#   BindingDB 34 34
#
# Tags Counts:
#   tag Seq# Seg# '0' '1' '-'
#   MoRFTest 113 135 28,581 8,181 48,392
#   MoRFTrain 141 172 28,812 8,902 47,440
#
# Optional bottom comments
#
>DP00958|Fasta=DP00958|UniProt=B3LPU2;C7GS09;G2WGJ2;A0A0L8VMH3;N1P1V7;H0GIH2;P46984;A6ZQF1;C8ZB35;A0A6C1DTD4|AlphaFoldDB=B3LPU2|Pfam=PF08738|PDB=4WXA;4WX8|BioGRID=33578|DIP=DIP-1474N|IntAct=P46984|MINT=P46984|STRING=4932.YJL184W|KEGG=sce:YJL184W|DisProt=DP00958
MKLPVAQYSAPDGVEKSFAPIRDDPRYMTTEGRTTGPSDHVLNAGQIDRDKPSEPERTKDGSQLTYLGQLRTQLTGLQDDINEFLTGRMELAKNKKKAGADEKRIQEEINQLLDGGDGDEDAV
0--------------------000000000000000000000000000000000000000000011111111111111111111111111111111111111111111111110000000000
011111111111111111111000000000000000000000000000000000000000000011111111111111111111111111111111111111111111111110000000000
>Q8YFU1|Fasta=Q8YFU1|UniProt=A0A0H3ANC4;Q2YMM2;A0A0E1X347;A0AB36PRZ8;C0RHL4;A0A7U8K8P4;A0AAE9LAS2;C7LGK0;A0A0F6APL3;A0AAW7BBF0|AlphaFoldDB=A0A0H3ANC4|KEGG=bov:BOV_0511|Pfam=PF00576|PDB=4Q14|STRING=359391.BAB1_0532
MGKLSTHVLDTAHGTPAAAMRVELYRIAASGTPELLKRVVTNLDGRTDAPLLSGDEMRTGIYELQFHVAEYFEGRGAELAHEPFLDLIPIRFGIADEDGNYHVPLLVSPWSYSTYRGS
1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
... more data
```
## The header section
The header section outlines the format of the data annotation and provides some statistical information about these annotations. Lines in the header start with the ‘#’ character, indicating non-data lines.
The header is divided into seven parts.
1.	Data name
2.	Optional top comments
3.	Sequence count.
4.	**Data annotation format:** See below.
5.	**ID counts:** See below.
6.	**Tag counts:** See below.
7.	Optional bottom comments

The **AFF (Data annotation format)** defines the structure for representing a single annotated protein sequence and consists of three main components:
- **Header Line:** Starts with a '>' character, followed by the sequence accession and one or more identifier phrases often referencing specialized databases separated by '|'. Each identifier phrase consists of the identifier source, the '=' character, and a list of IDs separated by ';'. The first identifier phrase with the source 'FASTA' refers to the identifier used by the original fasta file, which is usually the accession.
- **Amino Acid Sequence:** Provided on the line immediately following the header. It uses one-letter amino acid codes (e.g., A, R, N, D) and is limited to a single line.
- **Annotation Line(s):** One or more lines follow the amino acid sequence, each corresponding to a specific feature or **TAG**. These lines must be the same length as the amino acid sequence, with each character aligned to a residue. Standard annotation characters include '0' (feature absent), '1' (feature present), and '-' (unknown). Additional user-defined characters can also be used to represent other annotations.

The **ID counts:** count all identifier phrases used. In the above AFF file example, 135 out of the 141 total sequences have Pfam IDs, with 108 unique IDs. Note that it only considers the first ID of each identifier phrase when computing unique IDs.

The **Tag counts:** count the three standard annotation characters '0' (feature absent), '1' (feature present), and '-' (unknown) for each annotation TAG.

## The data section

More details soon.



