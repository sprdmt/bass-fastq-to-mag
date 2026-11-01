# Outputs

```
results/
├── qc/                          trimmed reads, cutadapt reports, FastQC (raw and trimmed)
├── reports/multiqc_reads.html   read QC summary
├── assembly/<group>/
│   ├── contigs.fa               MEGAHIT co-assembly (original MEGAHIT contig names)
│   ├── final_k.txt              k of the last MEGAHIT iteration (used for the graph)
│   ├── assembly_graph.gfa       assembly graph used by GraphBin2
│   ├── assembly_graph_stats.tsv segments and links of the graph
│   └── quast/report.tsv         assembly statistics
├── genes/<group>/               Prodigal proteins (.faa), genes (.fna), GFF
├── mapping/<group>/
│   ├── bam/<sample>.bam(.bai)   reads of every sample mapped to its group's co-assembly
│   ├── coverm_contig_mean.tsv   mean coverage per contig and sample
│   └── jgi_depth.txt            MetaBAT2-style depth table
├── binning/<group>/<binner>/
│   ├── run/                     raw output of MetaBAT2 / GenomeFace / SemiBin2 / TaxVAMB
│   └── contig2bin.tsv           contig-to-bin table
├── graphbin2/<group>/<binner>/contig2bin.tsv   graph-refined bins
├── binette/<group>/
│   ├── final_bins/              consensus bins
│   └── final_bins_quality_reports.tsv   CheckM2 completeness/contamination
├── anvio/<group>/               optional manual-refinement workspace (docs/manual_refinement.md)
│   ├── contigs.db, merged/PROFILE.db   anvi'o databases with collections 'Binette' and 'Refine'
│   ├── triage.tsv               keep / refine / drop decision for every Binette bin
│   ├── HOW_TO_REFINE.md         bins to review with their anvi-refine commands
│   └── candidates/              kept + refined bins (with mags.source: anvio)
├── mags/
│   ├── candidates/              all bins entering dereplication (<group>__<bin>.fa)
│   ├── drep/                    dRep output (clusters, representatives)
│   ├── fasta/<project>_MAG_<n>.fna   final MAGs
│   ├── mags.tsv                 names, provenance, CheckM2 values, dRep cluster
│   ├── gtdbtk/                  GTDB-Tk classification
│   ├── rna/                     barrnap and tRNAscan-SE results per MAG
│   └── mag_summary.tsv          one row per MAG: quality, MIMAG class, taxonomy, coverage
└── ena/                         ENA submission files (only with the ena_* targets)
logs/        one log file per job
benchmarks/  run time and peak memory per job (used to calibrate the resource model)
```

`mag_summary.tsv` is the main result: completeness and contamination (CheckM2), quality class,
GTDB taxonomy, genome size, contigs, N50, mean coverage per sample and the dRep cluster each MAG
represents. The quality class (`mimag_quality`) is one of

| class | definition |
|---|---|
| high-quality | MIMAG: > 90 % completeness, < 5 % contamination, 23S, 16S and 5S rRNA genes and >= 18 tRNAs |
| near-complete | >= 90 % completeness and <= 5 % contamination, without the rRNA/tRNA criteria (Almeida et al. 2019) |
| medium-quality | MIMAG: >= 50 % completeness, < 10 % contamination |
| low-quality | everything else: a MAG kept at exactly 10.00 % contamination (the `mags:` limits are inclusive, MIMAG's 10 % is not) or any MAG once the thresholds in `mags:` are lowered |

`mags.tsv` carries the same class from the CheckM2 values alone (`quality_checkm2`, without
high-quality, which needs the rRNA/tRNA genes) and the column `origin`: the Binette input set(s)
a bin was taken from unchanged, `binette` for a bin that Binette built from several sets, or
`anvio_refined` for a bin refined by hand.
