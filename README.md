# Analysis of Endogenous and Domesticated Viral Elements in Diptera

This repository contains all scripts and pipelines used to analyze endogenous viral elements (EVEs) and domesticated viral sequences in parasitoid (Tachinidae) and free-living Diptera. The goal is to explore the mechanisms of integration, conservation, and domestication of these viral sequences.

---

## Data Preparation

### Objective
Assemble genomes, compute quality metrics, and infer phylogeny to establish the foundation for the analysis.

### Scripts and Steps

1. **Genome Assembly**
   - **Script:** `Snakemake_Assembly`
   - **Input:** Illumina reads
   - **Output:** Assembled genomes
   - **Description:** Assembly was performed using Megahit (v1.2.9), Heterozygous regions were homogenized using Redundans, integrating tools such as Bwa-Mem2, Sspace, and GapCloser.

2. **Assembly Quality Metrics**
   - **Pipeline:** `Snakemake_BUSCO_Quast`
   - **Additional Scripts:**
     - `BUSCO_summary.sh` and `Quast_summary.sh`: Summarize results for each species.
     - `Mean_Quast_stats.py` and `Mean_BUSCO_stats.py`: Compute average metrics.
     - `Group_quast_busco.py`: Combine QUAST and BUSCO results.
   - **Tools:**
     - BUSCO (v5.4.5) with the `--mode genome` option and `--lineage diptera_odb10` dataset.
     - QUAST (v5.0.2).

3. **Phylogenetic Inference**
   - **Scripts:**
     - `BUSCO_phylogeny_choose_genes.sh`: Selects 300 conserved BUSCO genes.
     - `Snakemake_BUSCO_phylogeny`: Performs alignment, trimming, and tree construction.
     - `catfasta2phyml.pl`: Used by `Snakemake_BUSCO_phylogeny` to create partitions.
   - **Tools:**
     - Clustal Omega (v1.2.3) for alignment.
     - TrimAL (v1.4.rev15) with `-automated1` for automated trimming.
     - IQ-TREE (v2.0.3) with options:
       - `-m MFP`: Model selection.
       - `-alrt 2500` and `-bb 2500`: Bootstrap support.
       - `-bnni`: Optimize branch length.

---

## EVE Analysis Pipeline

### Objective
Identify, validate, cluster, and analyze viral elements integrated into the genomes of the studied species.

### Scripts and Steps

#### Homology Search

1. **Identification of Candidate Loci**
   - **Pipeline:** `Snakemake_viral_homology`
   - **Scripts Used:**
     - `Make_change_strand_mmseqs.py`: Adjusts sequence strand orientation.
     - `Overlapping_sequences_BUSCO_Viral_loci2.R`: Handles overlapping loci.
     - `Extract_and_translate_loci.py`: Extracts and translates candidate loci.
   - **Tools:** MMseqs2 (commit 5daca424b162cc5fdf0b9cd151aebed86975cbf6) with options:
     - `--evalue 0.01`: E-value threshold.
     - `--sensitivity 7.5`: Sensitivity setting for homology search.
     - `--min-seq-id 0.3`: Minimum sequence identity of 30%.

#### Taxonomic Assignment

2. **Validation of Candidate Loci**
   - **Scripts:**
     - `Create_filter_taxdb.sh`: Builds a custom taxonomic database.
     - `Snakemake_Tax_assignation_TopHit`: Performs taxonomic assignment.
   - **Tools:** MMseqs2 `easy-taxonomy` module with options:
     - `--lca-mode 4`: Top-hit based assignment.
     - `--evalue 1e-5`: E-value threshold.
   - **Description:** Retains only sequences confidently assigned to viral taxa.

#### Clustering

3. **Clustering of EVEs and Viral Sequences**
   - **Scripts:**
     - `Get_fasta_from_taxid.py`: Extracts Schizophora sequences.
     - `Remove_species_from_database.py`: Filters sequences from studied species.
     - `Clustering_without_fna.sh`: Groups sequences by similarity.
     - `Filter_clusters.py`: Filters clusters containing both viral and candidate sequences.
     - `Filter_clusters_with_90per_schizophora.py`: Removes clusters with >90% Schizophora sequences.
     - `Format_headers_keep_accession.py`: Formats sequence headers.
     - `Create_faa_fna_per_cluster.py`: Generates FASTA files per cluster.
     - `Candidates_vs_Refseq_homology_search.sh`: Matches candidates to RefSeq viral proteins.
     - `Get_viral_protein_info.sh` and `Get_viral_family_info.py`: Retrieve detailed viral protein and family information.
   - **Tools:** MMseqs2 with options:
     - `--cluster-mode 1`: Connected component clustering.
     - `--cov-mode 0`: Bi-directional coverage.
     - `--min-cov 0.3`: Minimum overlap of 30%.

#### Phylogenetic Analysis

4. **Cluster Alignment and Phylogeny**
   - **Pipeline:**
     - `Snakemake_align_clusters`: Aligns sequences.
     - `Snakemake_phylogenyML_clusters`: Constructs phylogenetic trees.
   - **Tools:** Clustal Omega, TrimAL (`-automated1`), and IQ-TREE (`-m MFP`, `-alrt 1000`).

#### Genomic Environment Analysis

5. **Validation of Genomic Integration**
   - **Scripts:**
     - `Snakemake_TE_search`: Searches for transposable elements using RepeatPeps.
     - `Snakemake_mapping_depth`: Maps reads to genomes and calculates scaffold depth.
     - `Snakemake_BUSCO_TE`: Counts BUSCOs and TEs per scaffold, merging overlaps.
       - **Sub-scripts:**
         - `TE_positions.py`: Identifies TE positions.
         - `TE_merge_overlaps.sh`: Merges overlapping TEs using Bedtools.
         - `TE_count_per_scaffold.py`: Counts TEs per scaffold.
     - `Regroup_all_info_candidates.py`: Consolidates candidate information.
     - `Find_Endogenized_candidates.py`: Validates EVEs as endogenous using depth, BUSCO, and TE criteria.
     - `Remove_TE.sh`: Removes loci overlapping TEs and filters non-endogenized sequences.
   - **Tools:**
     - Bwa-mem2 for mapping.
     - Samtools for depth calculation.
     - RepeatPeps for TE identification.

#### Domestication Analysis

6. **Domestication Signature Analysis**
   - **Scripts:**
     - `Get_clusters_and_monophyletic_groups.R`: Identifies clusters and monophyletic groups to analyze.
     - `Get_clusters_min_3seq.py`: Filters clusters with at least three sequences.
     - `Snakemake_codon_ali_clusters`: Generates codon alignments using MACSE.
     - `Trimming_Codon_ali_manually.sh`: Manually trims codon alignments for quality.
     - `Snakemake_dNdS_analysis`: Performs dN/dS ratio analysis.
       - **Sub-script:** `dNdS_analysis.py`: Calculates dN/dS ratios and applies statistical tests.
   - **Tools:**
     - MACSE for codon alignment.
     - TrimAL (`-automated1`, `-backtrans`) for trimming.
     - PAML for dN/dS ratio estimation.


#### Synteny Analysis

7. **Synteny Analysis of Domesticated Sequences**
   - **Pipeline:** `Snakemake_Synteny_domesticated`
   - **Scripts:**
     - `filter_insect_sequences.py`: Extracts insect-specific sequences from cluster FASTA files based on headers.
     - `extract_scaffolds.py`: Identifies and extracts scaffolds containing the sequences of interest.
     - `predict_genes.py`: Uses Metaeuk to predict gene models from extracted scaffolds. Generates output in GFF and FASTA formats.
     - `compare_synteny.py`: Identifies upstream and downstream genes and compares synteny across species.
   - **Tools:**
     - Metaeuk for gene prediction (`--min-length 33`, `--compressed 1`).
   - **Description:** This step determines whether domesticated sequences share conserved genomic contexts across species, indicating shared integration sites and evolutionary relationships.

---

## Statistical Tests

### Objective
Analyze the impact of lifestyle and data quality on the detection and domestication of EVEs.

### Scripts and Steps

1. **Lifestyle Effect**
   - **Script:** `Statistical_test_Diptera_Tachinidae.R`
   - **Description:** Phylogenetic generalized linear mixed model (PGLMM) to test whether parasitoids have more EVEs than free-living species.

2. **Data Quality Impact**
   - **Script:** `Statistical_test_genome_quality.R`
   - **Description:** Compares Horizon and DTOL genomes to evaluate the effect of data quality on EVE detection and validation.

---
