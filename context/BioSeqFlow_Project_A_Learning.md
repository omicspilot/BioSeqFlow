# BioSeqFlow — NGS Toolchain Learning Plan
## Project A: Learn the Tools Before Building the Product

> **Purpose**: This is the learning phase before building the BioSeqFlow GUI. The goal is to understand every NGS tool deeply enough to wrap it reliably in a production application. You will build a working Nextflow CLI pipeline as the learning deliverable — not the final product, but the foundation of it.
>
> **Output**: A working Nextflow pipeline (CLI only) covering RNA-seq and variant calling end to end. Published on GitHub under OmicsPilot as `bioseqflow-core`. This becomes the backend engine that the BioSeqFlow GUI wraps in Project B.
>
> **Philosophy**: Every tool you learn here will become a subprocess call in the FastAPI server. Understanding what the tool does, what it outputs, and what can go wrong is what makes the GUI reliable. Do not skip the biology — the GUI is useless if you cannot interpret what it shows.

---

## The Two Tracks You Are Building

**Track A — RNA-seq (gene expression analysis)**
```
Raw FASTQ
  → QC:           FastQC + MultiQC
  → Trimming:     fastp
  → Alignment:    STAR (splice-aware)
  → Quantify:     featureCounts
  → Statistics:   pyDESeq2
  → Output:       count matrix, DE results, volcano plot
```

**Track B — DNA Variant Calling**
```
Raw FASTQ
  → QC:           FastQC + MultiQC
  → Trimming:     fastp
  → Alignment:    BWA-MEM2
  → BAM cleanup:  samtools (sort, markdup, index)
  → Calling:      GATK HaplotypeCaller
  → Annotation:   SnpEff
  → Output:       annotated VCF
```

Both tracks share QC and trimming. The Nextflow pipeline accepts a sample sheet CSV and processes both tracks in parallel.

---

## Biological Foundation (Read Before Touching Any Tool)

Understanding the biology behind each step is what separates someone who runs tools from someone who can build reliable software around them.

### What NGS actually produces
A sequencer fragments DNA or RNA into millions of short pieces (reads, typically 75-150bp) and records the base sequence + a quality score per base. The result is a FASTQ file. The analytical challenge: given millions of anonymous fragments, what do they tell us about gene expression or genetic variation?

### Why RNA-seq needs a splice-aware aligner
RNA is reverse-transcribed into cDNA before sequencing. But mRNA has introns removed — a single 150bp read can span an exon-exon junction that does not exist in the reference genome DNA. A standard DNA aligner (BWA) would fail to map that read. STAR was built specifically to handle spliced reads.

### Why PCR duplicates must be removed
During library preparation, DNA is amplified by PCR. The same original molecule can produce 5, 10, or 50 identical copies. If you count those copies as independent evidence for a variant, you inflate allele frequencies and create false positives. `samtools markdup` flags these so GATK ignores them.

### Why normalisation is required for differential expression
Two samples can have 20M and 80M reads respectively due to technical variation in sequencing depth, not biological differences. Raw counts are therefore not comparable. DESeq2 estimates size factors per sample to correct for this before testing for differential expression.

### Why variant annotation matters
A VCF record like `chr17  43044295  .  A  T` is meaningless without biological context. SnpEff tells you: this is in BRCA1, it causes an amino acid change from Lysine to Methionine, it is predicted HIGH impact (likely damaging). That context is what makes a variant clinically interpretable.

### Key file formats — know these cold
| Format | Contents | Produced by | Consumed by |
|---|---|---|---|
| FASTQ | Raw reads + Phred quality scores | Sequencer | Aligners |
| FASTA | Reference genome sequences | NCBI/Ensembl | Aligners (index) |
| SAM/BAM | Aligned reads with mapping info | STAR, BWA-MEM2 | samtools, GATK |
| VCF/BCF | Genetic variants | GATK | SnpEff, bcftools |
| GTF/GFF | Gene/transcript annotations | Ensembl | featureCounts, STAR |
| BED | Genomic intervals (regions) | Various | bedtools |
| h5ad | AnnData single-cell objects | Scanpy | CellMetPro |

---

## Environment Setup (Do This First)

```bash
# Install Miniforge if not already done
# https://github.com/conda-forge/miniforge

# Create dedicated environment
mamba create -n bioseqflow-tools -c bioconda -c conda-forge \
  python=3.11 \
  fastqc multiqc fastp \
  star bwa-mem2 samtools \
  subread \
  gatk4 \
  snpeff \
  nextflow \
  pydeseq2 pandas matplotlib seaborn \
  bcftools bedtools

mamba activate bioseqflow-tools
```

**Test datasets** — use these throughout to keep runtimes short:
- nf-core test datasets: https://github.com/nf-core/test-datasets
- NCBI SRA small RNA-seq: any GSE with paired-end reads <1GB
- GATK test data: https://console.cloud.google.com/storage/browser/gatk-test-data

**Reference files** — download chr22 only for learning (much faster):
- Human genome GRCh38 chr22 FASTA: Ensembl FTP
- Matching GTF (same Ensembl release): Ensembl FTP
- Note: genome and annotation versions MUST match

---

## Week 1 — Quality Control

### Biological question
Is my sequencing data good enough to proceed, and what preprocessing does it need?

### Tools

**FastQC**
Generates per-read quality metrics. Key outputs to understand:
- Per-base sequence quality: Q20 = 1 error in 100 bases, Q30 = 1 in 1000. Bases below Q20 at the 3' end are normal and will be trimmed.
- Adapter content: sequencing reads through the insert into the adapter sequence. Must be removed.
- GC content distribution: should match the expected distribution for your organism (~50% for human). Deviation suggests contamination or bias.
- Overrepresented sequences: if one sequence appears in >1% of reads it is flagged. Often rRNA contamination or adapter dimers.
- Docs: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

**MultiQC**
Aggregates FastQC reports across samples into a single interactive HTML. Essential for cohort-level QC — spot outlier samples instantly.
- Docs: https://multiqc.info/docs/

**fastp**
All-in-one preprocessing: auto-detects adapters, trims low-quality bases, filters short reads, deduplicates, and outputs a JSON/HTML QC report. Faster and more feature-rich than Trimmomatic.
- Key parameters: `--qualified_quality_phred 20`, `--length_required 50`, `--detect_adapter_for_pe`
- Docs: https://github.com/OpenGene/fastp

### Tasks
- [ ] Run FastQC on 2 raw FASTQ files (paired-end)
- [ ] Open the HTML report and identify: worst quality region, adapter presence, GC distribution
- [ ] Run fastp to trim and filter
- [ ] Re-run FastQC on trimmed files
- [ ] Run MultiQC on both sets of FastQC reports
- [ ] Write notes: what changed after trimming, what would you flag as problematic?

### What to know after this week
- What a good vs bad FastQC report looks like
- What adapters are and why they cause alignment failures
- Why Q30 is the industry standard threshold
- How to interpret MultiQC's summary table for a cohort

---

## Week 2 — RNA-seq Alignment & Quantification

### Biological question
Which genes are expressed and at what level?

### Tools

**STAR**
Splice-aware RNA-seq aligner. Two-step: build index (once, slow), align reads (fast).
- Key concept: splice junctions. STAR discovers novel junctions from the data and maps reads that span them.
- Key output stats (in Log.final.out): % uniquely mapped (aim >75% for human), % multi-mapped, % too short
- Key parameters: `--outSAMtype BAM SortedByCoordinate`, `--quantMode GeneCounts`
- Docs: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

**samtools**
SAM/BAM manipulation. The Unix toolkit of bioinformatics.
- `samtools sort` — sort BAM by coordinate (required for indexing and most tools)
- `samtools index` — create .bai index for random access
- `samtools flagstat` — alignment summary statistics
- `samtools view -c -F 4` — count mapped reads
- SAM flag decoder: https://broadinstitute.github.io/picard/explain-flags.html
- Docs: http://www.htslib.org/doc/samtools.html

**featureCounts (Subread)**
Counts how many reads map to each gene in the GTF annotation. Produces the count matrix that goes into DESeq2.
- Key parameter: `--isPairedEnd`, `-s 0/1/2` (strandness — must match library prep protocol)
- Output: tab-separated count matrix (genes × samples)
- Docs: https://subread.sourceforge.net/featureCounts.html

### Tasks
- [ ] Build STAR genome index for chr22 (`--runMode genomeGenerate`)
- [ ] Align trimmed reads from Week 1
- [ ] Check `Log.final.out` — what % uniquely mapped?
- [ ] Sort and index the BAM with samtools
- [ ] Run `samtools flagstat` — understand each line of output
- [ ] Run featureCounts on the BAM
- [ ] Open the count matrix — how many genes have zero counts?

### What to know after this week
- Why RNA-seq needs a splice-aware aligner (not BWA)
- What SAM flags encode and how to read them
- What a count matrix is and why it is the input to differential expression
- What % uniquely mapped is acceptable and what causes low mapping rates

---

## Week 3 — Differential Expression Analysis

### Biological question
Which genes change their expression between two conditions, and by how much?

### Tools

**pyDESeq2**
Python implementation of DESeq2, the gold standard for RNA-seq differential expression. Core statistical framework:
1. Estimate size factors (normalise for sequencing depth)
2. Estimate dispersion (model gene-level variance)
3. Fit negative binomial model per gene
4. Wald test for differential expression
5. Benjamini-Hochberg FDR correction

Key output columns: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`
- Docs: https://pydeseq2.readthedocs.io/en/latest/

**Statistical concepts to understand**
- **Negative binomial model**: RNA-seq counts are overdispersed (variance > mean), so Poisson is not appropriate. The NB model accounts for this.
- **Size factors**: per-sample scaling factors that correct for sequencing depth. A sample with 2x more reads gets size factor ~2.
- **Dispersion**: gene-level variance estimate. Highly expressed genes tend to have lower dispersion.
- **log2 fold change**: log2FC=1 means 2x expression, log2FC=2 means 4x, log2FC=-1 means 0.5x. Log scale because fold changes are multiplicative, not additive.
- **padj (FDR)**: adjusted p-value after correcting for multiple testing across ~20,000 genes. Standard cutoff: padj < 0.05.
- **Volcano plot**: log2FC on x-axis, -log10(padj) on y-axis. Upper right = upregulated significant. Upper left = downregulated significant.

**StatQuest resource**: https://www.youtube.com/watch?v=UFB993xufUU — watch before coding

### Tasks
- [ ] Load count matrix from Week 2 into pyDESeq2
- [ ] Define experimental design (two conditions from metadata)
- [ ] Run DESeq2 workflow: fit, test, results
- [ ] Filter: padj < 0.05 AND |log2FC| > 1
- [ ] Plot volcano: colour significant genes red
- [ ] Plot heatmap of top 50 DE genes (rows = genes, cols = samples)
- [ ] Run PCA on normalised counts — do samples cluster by condition?

### What to know after this week
- Why raw counts cannot be compared directly between samples
- What FDR correction is and why a p-value of 0.05 is meaningless without it
- How to read a volcano plot and identify biologically interesting genes
- What log2FC = 1 means biologically

---

## Week 4 — DNA Alignment & Variant Calling

### Biological question
What genetic variants does this sample carry relative to the reference genome?

### Tools

**BWA-MEM2**
DNA read aligner. Unlike STAR, does not need to handle splicing — maps DNA reads to genomic DNA.
- Same two-step process: index reference, then align
- Output: SAM file (convert to BAM with samtools)
- Docs: https://github.com/bwa-mem2/bwa-mem2

**samtools markdup**
Marks PCR duplicate reads so variant callers can ignore them.
- Pipeline: `samtools sort → samtools markdup → samtools index`
- Check % duplicates in the stats output — >30% is high for WGS, up to 50% is acceptable for targeted panels

**GATK HaplotypeCaller**
Standard germline variant caller. Performs local de novo assembly around candidate variant sites rather than just pileup counting — more accurate, especially for indels.
- Key outputs: GVCF (per-sample, used for joint genotyping) or VCF (final variants)
- Key VCF fields to understand: DP (depth), GQ (genotype quality), GT (genotype: 0/0=homref, 0/1=het, 1/1=homalt), AF (allele frequency)
- Docs: https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller

**SnpEff**
Annotates variants with functional consequences.
- Impact categories: HIGH (frameshift, stop gained — likely damaging), MODERATE (missense — may be damaging), LOW (synonymous — likely benign), MODIFIER (intronic, intergenic — unknown impact)
- Integrates with ClinVar, dbSNP
- Docs: https://pcingola.github.io/SnpEff/

**StatQuest resource**: https://www.youtube.com/watch?v=0OWkTzTOfKk

### Tasks
- [ ] Index reference genome with BWA-MEM2
- [ ] Align DNA reads to reference
- [ ] Sort, mark duplicates, index BAM with samtools
- [ ] Check % duplicates in markdup stats
- [ ] Run GATK HaplotypeCaller in GVCF mode
- [ ] Genotype GVCF to produce final VCF
- [ ] Annotate VCF with SnpEff
- [ ] Count: how many SNPs vs indels? How many HIGH impact variants?
- [ ] Open a few HIGH impact variants and look them up in ClinVar

### What to know after this week
- What PCR duplicates are and why they must be marked (not removed)
- The difference between a SNP and an indel
- How to read a VCF record: CHROM, POS, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample columns
- What GATK's GVCF mode is and why joint genotyping exists
- How to interpret SnpEff impact categories

---

## Week 5 — Extended Tools

### Tools

**Salmon**
Pseudo-alignment based quantification. Much faster than STAR + featureCounts because it maps k-mers to a transcript index rather than aligning to the full genome.
- When to prefer Salmon: speed is priority, BAM file not needed downstream
- When to prefer STAR + featureCounts: you need the BAM for splicing analysis, fusion detection, or visual inspection in IGV
- Docs: https://salmon.readthedocs.io/en/latest/

**bcftools**
VCF/BCF manipulation complement to samtools.
- `bcftools stats` — VCF quality statistics, Ti/Tv ratio (should be ~2.0-2.1 for WGS)
- `bcftools filter` — filter variants by QUAL, DP, GQ
- `bcftools annotate` — add INFO fields from external databases
- `bcftools merge` — merge VCF files from multiple samples
- Docs: http://samtools.github.io/bcftools/bcftools.html

**bedtools**
Genomic interval operations. Think of it as set algebra for genomic coordinates.
- `bedtools intersect` — find overlapping regions between two BED/VCF files
- `bedtools coverage` — calculate per-base depth over target regions
- Essential for targeted sequencing analysis (WES, panels)
- Docs: https://bedtools.readthedocs.io/en/latest/

### Tasks
- [ ] Add Salmon as alternative Track A quantification option
- [ ] Compare Salmon vs featureCounts counts on same dataset — how correlated are they?
- [ ] Run `bcftools stats` on the VCF from Week 4, check Ti/Tv ratio
- [ ] Use `bcftools filter` to keep only PASS variants with DP>10 and GQ>20
- [ ] Use `bedtools intersect` to find variants overlapping a gene list BED file

---

## Week 6 — Nextflow Pipeline Integration

### Goal
Wrap everything from weeks 1-5 into a proper Nextflow pipeline. This is `bioseqflow-core` v0.1 — the engine that BioSeqFlow GUI will wrap.

### Study these first
- [nf-core/rnaseq](https://github.com/nf-core/rnaseq) — reference RNA-seq pipeline implementation
- [nf-core/sarek](https://github.com/nf-core/sarek) — reference variant calling pipeline
- [nf-core/modules](https://github.com/nf-core/modules) — reusable module definitions (import rather than rewrite)

### Pipeline structure
```
bioseqflow-core/
├── main.nf                          # Entry point, workflow selection
├── nextflow.config                  # Profiles: local, docker, slurm, aws
├── modules/
│   ├── fastqc/main.nf
│   ├── fastp/main.nf
│   ├── star/align/main.nf
│   ├── star/genomegenerate/main.nf
│   ├── bwa_mem2/align/main.nf
│   ├── samtools/sort/main.nf
│   ├── samtools/markdup/main.nf
│   ├── samtools/index/main.nf
│   ├── featurecounts/main.nf
│   ├── gatk/haplotypecaller/main.nf
│   ├── snpeff/main.nf
│   └── pydeseq2/main.nf
├── workflows/
│   ├── rnaseq.nf                    # Track A
│   └── variantcalling.nf            # Track B
├── subworkflows/
│   └── qc_trim.nf                   # Shared QC+trimming
├── conf/
│   ├── params.config                # Default parameters
│   └── modules.config               # Tool-specific settings
├── assets/
│   └── sample_sheet_template.csv
├── bin/
│   └── generate_report.py
├── .github/
│   └── workflows/ci.yml
├── CHANGELOG.md
├── LICENSE
└── README.md
```

### Key Nextflow concepts to apply
- **Channels**: sample sheet read as channel, branched to RNA-seq and variant calling by sample type column
- **Processes**: each tool is one process, defined inputs/outputs/container
- **Subworkflows**: QC+trimming extracted as shared subworkflow used by both tracks
- **Profiles**: `local` (no containers), `docker` (Docker), `singularity` (HPC), `aws` (ECS)
- **Resume**: `-resume` skips completed processes — critical for long pipeline runs
- **nf-core modules**: `nf-core modules install fastqc` rather than writing from scratch

### Tasks
- [ ] Study nf-core/rnaseq structure for 2 hours before writing any code
- [ ] Install nf-core tools: `pip install nf-core`
- [ ] Import existing nf-core modules for FastQC, fastp, STAR, samtools, featureCounts
- [ ] Write custom modules for pyDESeq2 and report generation
- [ ] Build shared `qc_trim` subworkflow
- [ ] Build `rnaseq.nf` workflow wiring QC → STAR → featureCounts → pyDESeq2
- [ ] Build `variantcalling.nf` workflow wiring QC → BWA → GATK → SnpEff
- [ ] Write `main.nf` entry point with workflow selection from sample sheet
- [ ] Write `nextflow.config` with local and docker profiles
- [ ] Test RNA-seq track end to end with test dataset
- [ ] Test variant calling track end to end with test dataset
- [ ] Set up GitHub Actions CI: validate nextflow syntax, run with test profile
- [ ] Write README with quickstart, sample sheet format, output structure
- [ ] Tag and publish as v0.1.0

---

## Success Criteria

By the end of this project you should be able to:

- [ ] Run a complete RNA-seq pipeline from raw FASTQ to differential expression results
- [ ] Run a complete variant calling pipeline from raw FASTQ to annotated VCF
- [ ] Explain every tool decision and connect it to the biology
- [ ] Read and interpret a FastQC report, SAM flagstat output, and VCF record
- [ ] Describe what happens inside STAR and GATK at a conceptual level
- [ ] Explain the difference between alignment-based and pseudo-alignment quantification and when to use each
- [ ] Have `bioseqflow-core` v0.1.0 published on GitHub under OmicsPilot
- [ ] Be ready to wrap every tool as a subprocess call in the BioSeqFlow FastAPI server

---

## How This Feeds Project B

Every tool you learn here becomes:

| Tool | In BioSeqFlow GUI |
|---|---|
| FastQC | Subprocess call → parse JSON → render QC cards |
| fastp | Subprocess call → parse JSON → before/after stats |
| STAR | Long-running subprocess → stream stdout → WebSocket progress |
| samtools | Subprocess calls → parse flagstat → alignment stats panel |
| featureCounts | Subprocess call → parse output → count matrix preview |
| pyDESeq2 | Python library call → results dataframe → volcano plot |
| GATK | Long-running subprocess → stream log → WebSocket progress |
| SnpEff | Subprocess call → parse annotated VCF → variant table |

Understanding what each tool prints, what it produces, and what can go wrong is what makes the FastAPI server reliable. That is the real purpose of this learning phase.
