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
     - QUAST (v5.0.2)

3. **Phylogenetic Inference**
   - **Scripts:**
     - `BUSCO_phylogeny_choose_genes.sh`: Selects 300 conserved BUSCO genes.
     - `Snakemake_BUSCO_phylogeny`: Performs alignment, trimming, and tree construction.
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
     - `Format_headers_keep_accession.py`: Formats sequence headers.
     - `Create_faa_fna_per_cluster.py`: Generates FASTA files per cluster.
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
     - `Find_Endogenized_candidates.py`: Validates EVEs using sequencing depth, BUSCO presence, and TE overlap.
     - `Regroup_all_info_candidates.py`: Consolidates candidate information.
     - `Remove_TE.sh`: Removes loci overlapping TEs.
   - **Tools:**
     - Bwa-mem2 for mapping.
     - Samtools coverage for depth calculation.
     - RepeatPeps for TE identification.

#### Domestication Analysis

6. **Domestication Signature Analysis**
   - **Scripts:**
     - `Snakemake_dNdS_analysis`: Analyzes dN/dS ratios to detect purifying selection.
     - `Get_clusters_min_3seq.py`: Filters clusters with at least three sequences.
   - **Tools:**
     - MACSE for codon alignment.
     - TrimAL (`-automated1`) for trimming.
     - PAML for dN/dS ratio estimation.

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
