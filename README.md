# heuristicount
Supply some known barcodes (sgRNAs), and discover the rest using `heuristicount.py`. Heuristicount.py is unique that it uses pattern matching to discover new barcodes in your data, while also retaining strict rules to prevent overcounting.

```vim
❯ python heuristicount.py design_key.fasta ABA1-10_S149_L004_R1_001.fastq.gz ABA1-10_S149_L004_R2_001.fastq.gz > barcode_results.tsv
[09:44:44] Initializing heuristic barcode counting                                                                                                                                                                                                          heuristicount.py:179
           Reading barcodes...                                                                                                                                                                                                                              heuristicount.py:184
           Reading results from FASTQ...                                                                                                                                                                                                                    heuristicount.py:191
[09:45:18] Finding orientation of paired-end reads...                                                                                                                                                                                                       heuristicount.py:200
           Finding forward coordinates...                                                                                                                                                                                                                   heuristicount.py:223
           Finding reverse coordinates                                                                                                                                                                                                                      heuristicount.py:231
[09:45:19] Finding forward junctions...                                                                                                                                                                                                                     heuristicount.py:242
           Finding reverse junctions...                                                                                                                                                                                                                     heuristicount.py:249
           Mapping junctions to barcodes...                                                                                                                                                                                                                 heuristicount.py:256
           Counting barcodes at junctions...                                                                                                                                                                                                                heuristicount.py:264
[09:45:41] Preparing results...                                                                                                                                                                                                                             heuristicount.py:269
           Heuristic filtering and identifying new barcodes...                                                                                                                                                                                              heuristicount.py:277
                                                                                                                                                                                                                                                            heuristicount.py:384
                           heuristicount.py                             Summary                                                                                                                                                                                                 
            ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━                                                                                                                                                                                                
                             Input & Config                                                                                                                                                                                                                                     
                                   Barcodes                    design_key.fasta                                                                                                                                                                                                 
                              Forward Reads   ABA1-10_S149_L004_R2_001.fastq.gz                                                                                                                                                                                                 
                              Reverse Reads   ABA1-10_S149_L004_R1_001.fastq.gz                                                                                                                                                                                                 
                                    Threads                                  12                                                                                                                                                                                                 
                           Operating System                               Linux                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                
                                 Heuristics                                                                                                                                                                                                                                     
                             Barcode Length                                  20                                                                                                                                                                                                 
                             Forward Offset                                 100                                                                                                                                                                                                 
                             Reverse Offset                                  82                                                                                                                                                                                                 
                           Forward Junction                         TAGT...GTTT                                                                                                                                                                                                 
                           Reverse Junction                         AAAC...ACTA                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                
                    Barcode Alignment Stats                                                                                                                                                                                                                                     
                          Barcodes Declared                                6680                                                                                                                                                                                                 
                  Documented Barcodes Found                                6666                                                                                                                                                                                                 
                Undocumented Barcodes Found                               32696                                                                                                                                                                                                 
                                Total Reads                            19275050                                                                                                                                                                                                 
                   Documented Barcode Reads                            17267584                                                                                                                                                                                                 
                 Undocumented Barcode Reads                              109911                                                                                                                                                                                                 
                        Documented Fraction                              0.8959                                                                                                                                                                                                 
                      Undocumented Fraction                              0.0057                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                
                  Top 5 Documented Barcodes                                                                                                                                                                                                                                     
                       AGGCGACTCCACAATCACTG                               18003                                                                                                                                                                                                 
                       AATTGATCAATCAGATTGGT                               17574                                                                                                                                                                                                 
                       AGCAGATCGACTTCATATTT                               17465                                                                                                                                                                                                 
                       CCATCACTACCATCCATAAT                               14902                                                                                                                                                                                                 
                       CAGCTCGGCTCAAATCATTC                               13845                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                
                Top 5 Undocumented Barcodes                                                                                                                                                                                                                                     
                      ACTCCACAAACTGAAAGATA*                                1535                                                                                                                                                                                                 
                      ATCCTCATTGGTCGAAGTTT*                                 709                                                                                                                                                                                                 
                      AGGCGACTCAACAATCACTG*                                 642                                                                                                                                                                                                 
                      AGCTTTAATTTAATCACCAC*                                 606                                                                                                                                                                                                 
                      AGGCGGCTCCACAATCACTG*                                 602                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                
                           Finished at 2023-09-26 09:45:41.374100
```

# targets
Find genomic targets for your barcodes using `targets.py`. Targets.py circularizes genomes and quickly finds all valid targets using Bowtie, providing a report to help design or interpret experiments. 
```vim
❯ python targets.py heuristic_key.fasta GCF_009759685.1.gb 1 > heuristic_targets.tsv
[09:46:35] Initializing barcode target seeker                                                                                                                                                                                                                     targets.py:298
           Circularizing genome...                                                                                                                                                                                                                                targets.py:317
           Annotating regions to identify...                                                                                                                                                                                                                      targets.py:328
           Aligning annotations to genome...                                                                                                                                                                                                                      targets.py:334
[09:46:38] Generating genome maps...                                                                                                                                                                                                                              targets.py:337
[09:46:43] Finding matches...                                                                                                                                                                                                                                     targets.py:341
[09:46:44] Cleaning up...                                                                                                                                                                                                                                         targets.py:353
                                                                                                                                                                                                                                                                  targets.py:460
                                 targets.py                   Summary                                                                                                                                                                                                           
            ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━                                                                                                                                                                                                          
                             Input & Config                                                                                                                                                                                                                                     
                                   Barcodes       heuristic_key.fasta                                                                                                                                                                                                           
                        Genbank Genome File        GCF_009759685.1.gb                                                                                                                                                                                                           
                       Number of Mismatches                         1                                                                                                                                                                                                           
                                    Threads                        12                                                                                                                                                                                                           
                           Operating System                     Linux                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                
                                 Heuristics                                                                                                                                                                                                                                     
                                   Organism   Acinetobacter baumannii                                                                                                                                                                                                           
                                Chromosomes                         2                                                                                                                                                                                                           
                                Total Genes                      3805                                                                                                                                                                                                           
                          Overlapping Genes                       724                                                                                                                                                                                                           
                      Ambiguous Coordinates                      6511                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                
                      Barcode Mapping Stats                                                                                                                                                                                                                                     
                             Genes Targeted                       408                                                                                                                                                                                                           
                         Targeting Barcodes                     12468                                                                                                                                                                                                           
                         0 Mismatch Spacers                      1615                                                                                                                                                                                                           
                         1 Mismatch Spacers                     10846                                                                                                                                                                                                           
                     Non-targeting Barcodes                     26921                                                                                                                                                                                                           
                        Intergenic Barcodes                         0                                                                                                                                                                                                           
                         Ambiguous Barcodes                         7                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                
                      Finished at 2023-09-26 09:46:44.121173   

```
