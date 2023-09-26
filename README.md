# heuristicount
Supply some known barcodes (sgRNAs), and discover the rest using `heuristicount.py`.

```
❯ python heuristicount.py design_key.fasta ABA1-10_S149_L004_R1_001.fastq.gz ABA1-10_S149_L004_R2_001.fastq.gz > barcode_results.tsv
```
![image](https://github.com/ryandward/barcoder/assets/6970996/b7be2c0f-8154-41bd-9c5a-b1a8cc4b383a)
```
❯ head barcode_results.tsv | tabulate   
--------------------  ----
ATTGCGTTACCAGCATAGAT  4010
AAAGTGTGATTCAATAAACA  5280
ATTTACGTGGATCTACGCCT  5871
AGCGGTTGCATTTTAATATC  6780
CGCGCATCTTTGTTAATGTG  3734
TTAACGGCATCTACACTATC  3333
GGCATCCAGAAAGAAGCTTC  3625
AAGTGCTAGAGGACCTAGTA  3073
ATGTTATGTCATCTTGAAGG  2832
CTTTAATGCGATCGTAGGGC  1886
--------------------  ----
```

# targets
Find genomic targets for your barcodes using `targets.py`
```
❯ awk '{print ">" $1; gsub(/\*$/, "", $1); print $1}' barcode_results.tsv > heuristic_key.fasta
❯ python targets.py heuristic_key.fasta GCF_009759685.1.gb 1 > heuristic_targets.tsv
```
![image](https://github.com/ryandward/barcoder/assets/6970996/be39e064-cc05-43a7-aded-ef2210ce6232)
```
❯ head heuristic_targets.tsv | tabulate -1                   
name                  spacer                pam      len  locus_tag    chr         target                  mismatches  diff    coords              offset  sp_dir    tar_dir
--------------------  --------------------  -----  -----  -----------  ----------  --------------------  ------------  ------  ----------------  --------  --------  ---------
ATTGCGTTACCAGCATAGAT  ATTGCGTTACCAGCATAGAT  TGG       20  GO593_02110  CP046654.1  ATTGCGTTACCAGCATAGAT             0  -       431230..431250         378  F         R
AAAGTGTGATTCAATAAACA  AAAGTGTGATTCAATAAACA  AGG       20  GO593_09085  CP046654.1  AAAGTGTGATTCAATAAAtA             1  a2G     1914661..1914681       201  R         F
ATTTACGTGGATCTACGCCT  ATTTACGTGGATCTACGCCT  AGG       20  GO593_07940  CP046654.1  ATTTACGTGGATCTACGCCa             1  t1A     1661276..1661296       143  R         F
AGCGGTTGCATTTTAATATC  AGCGGTTGCATTTTAATATC  AGG       20  GO593_01430  CP046654.1  AGCGGTTGCATTTTAATgTC             1  g18A    284082..284102         170  F         R
CGCGCATCTTTGTTAATGTG  CGCGCATCTTTGTTAATGTG  CGG       20  GO593_04205  CP046654.1  CGCGCATCTTTGTTAAcGTG             1  c17T    883479..883499         165  F         R
TTAACGGCATCTACACTATC  TTAACGGCATCTACACTATC  TGG       20  GO593_11615  CP046654.1  TTAACGGCATCTtCACTATC             1  t13A    2446470..2446490       279  F         R
GGCATCCAGAAAGAAGCTTC  GGCATCCAGAAAGAAGCTTC  TGG       20  GO593_13095  CP046654.1  GGCATCCAGAAAGAAGCTTC             0  -       2736511..2736531       111  F         R
AAGTGCTAGAGGACCTAGTA  AAGTGCTAGAGGACCTAGTA  AGG       20  GO593_17785  CP046654.1  AAGTGCTAGAGcACCTAGTA             1  c12G    3724024..3724044       154  F         R
AAAGTGTGATTCAATAAATA  AAAGTGTGATTCAATAAATA  AGG       20  GO593_09085  CP046654.1  AAAGTGTGATTCAATAAATA             0  -       1914661..1914681       201  R         F

```
