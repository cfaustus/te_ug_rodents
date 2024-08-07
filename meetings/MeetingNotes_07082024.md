Meeting Metaviromics 7 August 2024
================
laura
2024-08-07

## Meeting Notes

Attendees: Sophie Gryseels, Magdalini Bletsa, Philippe Selhorst, Jayna
Raghwani, Christina Faust, Laura Bergner

Sophie - challenging museum specimens and human samples (FFPE)  
Wildlife trafficking, wild meat sampling  
TE as screening tool

Magda - metagenomics, hepaciviruses  
FFPE - HIV using hybrid capture  
Pathogens carried by migratory birds  
Ancient DNA lab, hybrid capture for Hep B in ancient human samples

Philippe - metagenomic protocols for Nanopore

Less good experience with Twist  
Using depletion instead of enrichment  
Mastomys, positive selection worked less well than depletion, but don’t
know about what was in the panel or what was done by the company  
SISPA generates short pieces so haven’t been using adaptive sampling

Ordering probes family by family  
CoV, ParamyxoV, Flavis  
Tested on known positive samples  
Worked well with Illumina seq  
As you combine probes how does that affect efficiency

RNA viruses mainly of interest  
Balancing depth vs breadth  
38 virus families, most RNA  
Drod vir database & VIRION  
Include viruses from metagenomic surveys as well as some from
databases  
460 virus genomes

They are doing ~400 genomes per family  
Flavis started with 8000, after clustering had 800 genomes used to
design probes  
Similar for CoV and ParamyxoV

Not sure about combining DNA and RNA viruses  
Better to use nuclease treatment or sensitivity is lost  
Need to split workflow  
Could use DNA transcripts in a RNA workflow but won’t get genomes

They are using myBaits, RNA probes  
50-60K probes per family  
Smaller tiling  
Have not combined them  
150K probes in the combined set, unsure how that will affect results

Less genomes more families vs more families less genomes?

Current panels cluster at 90% identity  
Can capture different strains of the same virus (in HIV)

Have only done known positive samples so far  
Real bio controls would be better than something synthetic  
Spike in controls for non rodent pathogen on panel  
They’ve done spike in control in just pure water  
We are wanting also to quantify variation between samples  
Importance of negative controls for monitoring cross sample
contamination

They have also found more viruses with shrews  
They also found a CoV but having difficulty assembling full genome  
They were sequencing from organs but have gotten nearly full genomes
from shrew faecal samples in Glasgow  
Combine data from CoVs across teams showing geographical spread

Timeline Glasgow: doing some more metagenomic sequencing  
Start to develop panel in ~December  
First run test in new year  
Test out on samples from Oxford  
Eventually do sequencing in country

Problem with targeting what is already known in rodents, might miss out
on surprising/new  
Per family - mammalian all included  
Main question is around host switching so don’t want to exclude anything

Could use ancestral reconstruction to look for viruses that don’t exist
in contemporary dataset, include those in bait design  
Doing wider phylo analysis and ancestral reconstruction to guess at
diversity

To make it universal, extra RNA and do RT with overhang instead of
random primer. Use template switching RT step then amplify from those
overhung sites (similar to S9N)  
Then could design downstream steps including barcodes to be more
platform agnostic  
Potential to combine S9N, hybrid barcodes and twist hyb?  
Could potentially be more sensitive

Crispr based negative selection to get rid of host reads  
Can use with template switching  
Can be used to get rid of host RNA  
Developed for Mastomys and humans  
What is the most abundant RNA  
Wanting to make a central repository for Crispr panels so this data can
be shared  
Going to make an Aedes/Culex panel  
Could even combined with target enrichment  
Downside is that you need one probe or Cas9 module for each molecule  
Originally used in transcriptomics on cells to deplete rRNA  
More difficult in DNA samples  
Expensive for an initial order but includes a lot  
Order two primers, anneal and elongate  
One ds 200 bp molecule per target

Benchmarking - on the virus level similar but on the depletion level it
works better  
Compared to commercial rRNA kits  
More sensitive compared to generic hosts included in commercial panels  
5 euros per sample  
Not yet published

Get away from commercial assays to a more flexible and modular
approach  
Discovery funding - apply for something together in the future?
