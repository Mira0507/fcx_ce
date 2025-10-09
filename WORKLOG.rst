fcx_ce
======

@Mira0507

- prep conda env
    - requirements: ``requirements.txt``
    - env: ``env``
    - archived: ``env.archived.yaml``

- prep input files
    - anndata manually copied
        - from ``../../trujilloae/thalamus_atlas/Subsetted_cell_types/NEU/ExNeu1_FOR_DEG.h5ad``
          to ``input/thalamus_excitatory/ExNeu1_FOR_DEG.h5ad``
        - from ``../../trujilloae/thalamus_atlas/Subsetted_cell_types/NEU/ExNeu2_FOR_DEG.h5ad``
          to ``input/thalamus_excitatory/ExNeu1_FOR_DEG.h5ad``
    - bam copied by running ``workflow/thalamus_excitatory/prep_input.py``
        - from:
            - ``../../../data/biogen/CELLRANGER``
            - ``../../../data/marsan_2023_thalamus/CELLRANGER``
            - ``../../../data/syn52383413/CELLRANGER``
        - to: 

        .. code-block:: bash

            $ tree input/thalamus_excitatory/bam/
            input/thalamus_excitatory/bam/
            ├── 3075-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 3075-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 3514-T_Control-Biogen_possorted_genome_bam.bam
            ├── 3514-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 3549-T_Control-Biogen_possorted_genome_bam.bam
            ├── 3549-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 3676-T_Control-Biogen_possorted_genome_bam.bam
            ├── 3676-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 3733-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 3733-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 3862-T_Control-Biogen_possorted_genome_bam.bam
            ├── 3862-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 3937-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 3937-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 4571-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 4571-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 4984-T_Control-Biogen_possorted_genome_bam.bam
            ├── 4984-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 5166-T_Control-Biogen_possorted_genome_bam.bam
            ├── 5166-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 5189-T_Control-Biogen_possorted_genome_bam.bam
            ├── 5189-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 6025-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 6025-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 6203-T_Control-Biogen_possorted_genome_bam.bam
            ├── 6203-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 6243-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 6243-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 6283-T_Control-Biogen_possorted_genome_bam.bam
            ├── 6283-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 6400-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 6400-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 6433-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 6433-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 6562-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 6562-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 6636-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 6636-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 6658-T_Control-Biogen_possorted_genome_bam.bam
            ├── 6658-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 6863-T_Control-Biogen_possorted_genome_bam.bam
            ├── 6863-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 7024-T_Control-Biogen_possorted_genome_bam.bam
            ├── 7024-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 7197-T_Control-Biogen_possorted_genome_bam.bam
            ├── 7197-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 7228-T_Control-Biogen_possorted_genome_bam.bam
            ├── 7228-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 7624-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 7624-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 8047-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 8047-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 8130-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 8130-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 8251-T_Control-Biogen_possorted_genome_bam.bam
            ├── 8251-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── 8699-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 8699-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 8813-T_FTD-Biogen_possorted_genome_bam.bam
            ├── 8813-T_FTD-Biogen_possorted_genome_bam.bam.bai
            ├── 8866-T_Control-Biogen_possorted_genome_bam.bam
            ├── 8866-T_Control-Biogen_possorted_genome_bam.bam.bai
            ├── Control10_Control-Marsan_possorted_genome_bam.bam
            ├── Control10_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control11_Control-Marsan_possorted_genome_bam.bam
            ├── Control11_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control1_Control-Marsan_possorted_genome_bam.bam
            ├── Control1_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control2_Control-Marsan_possorted_genome_bam.bam
            ├── Control2_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control3_Control-Marsan_possorted_genome_bam.bam
            ├── Control3_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control4_Control-Marsan_possorted_genome_bam.bam
            ├── Control4_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control5_Control-Marsan_possorted_genome_bam.bam
            ├── Control5_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control6_Control-Marsan_possorted_genome_bam.bam
            ├── Control6_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control7_Control-Marsan_possorted_genome_bam.bam
            ├── Control7_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control8_Control-Marsan_possorted_genome_bam.bam
            ├── Control8_Control-Marsan_possorted_genome_bam.bam.bai
            ├── Control9_Control-Marsan_possorted_genome_bam.bam
            ├── Control9_Control-Marsan_possorted_genome_bam.bam.bai
            ├── D19-12358_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12358_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12359_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12359_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12360_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12360_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12361_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12361_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12362_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12362_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12363_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12363_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12365_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12365_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12366_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12366_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12368_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12368_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12369_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12369_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12370_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12370_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12371_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12371_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12372_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12372_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12373_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12373_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12374_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12374_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12375_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12375_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12376_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12376_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12377_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12377_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12378_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12378_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12379_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12379_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12380_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12380_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12381_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12381_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12382_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12382_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12383_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12383_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12384_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12384_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12385_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12385_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12386_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12386_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12387_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12387_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12388_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12388_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12389_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12389_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12390_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12390_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12391_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12391_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12392_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12392_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12393_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12393_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12394_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12394_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12395_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12395_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12396_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12396_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12397_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12397_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12398_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12398_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12399_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12399_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12400_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12400_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12401_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12401_Control-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12402_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12402_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12403_AD-Mathys_possorted_genome_bam.bam
            ├── D19-12403_AD-Mathys_possorted_genome_bam.bam.bai
            ├── D19-12404_Control-Mathys_possorted_genome_bam.bam
            ├── D19-12404_Control-Mathys_possorted_genome_bam.bam.bai
            ├── FTLD-GRN10_FTD-Marsan_possorted_genome_bam.bam
            ├── FTLD-GRN10_FTD-Marsan_possorted_genome_bam.bam.bai
            ├── FTLD-GRN1_FTD-Marsan_possorted_genome_bam.bam
            ├── FTLD-GRN1_FTD-Marsan_possorted_genome_bam.bam.bai
            ├── FTLD-GRN2_FTD-Marsan_possorted_genome_bam.bam
            ├── FTLD-GRN2_FTD-Marsan_possorted_genome_bam.bam.bai
            ├── FTLD-GRN3_FTD-Marsan_possorted_genome_bam.bam
            ├── FTLD-GRN3_FTD-Marsan_possorted_genome_bam.bam.bai
            ├── FTLD-GRN4_FTD-Marsan_possorted_genome_bam.bam
            ├── FTLD-GRN4_FTD-Marsan_possorted_genome_bam.bam.bai
            ├── FTLD-GRN5_FTD-Marsan_possorted_genome_bam.bam
            ├── FTLD-GRN5_FTD-Marsan_possorted_genome_bam.bam.bai
            ├── FTLD-GRN6_FTD-Marsan_possorted_genome_bam.bam
            ├── FTLD-GRN6_FTD-Marsan_possorted_genome_bam.bam.bai
            ├── FTLD-GRN7_FTD-Marsan_possorted_genome_bam.bam
            ├── FTLD-GRN7_FTD-Marsan_possorted_genome_bam.bam.bai
            ├── FTLD-GRN8_FTD-Marsan_possorted_genome_bam.bam
            ├── FTLD-GRN8_FTD-Marsan_possorted_genome_bam.bam.bai
            ├── FTLD-GRN9_FTD-Marsan_possorted_genome_bam.bam
            └── FTLD-GRN9_FTD-Marsan_possorted_genome_bam.bam.bai

            0 directories, 194 files

    - metadata: ``input/thalamus_excitatory/combined_thalamus_metadata.csv``

- read aggregation
    - by group: 12 groups (for sashimi plot in IGV)
        - group 1: Marsan, control, ExNeu1
        - group 2: Marsan, control, ExNeu2
        - group 3: Marsan, FTD, ExNeu1
        - group 4: Marsan, FTD, ExNeu2

        - group 5: Mathys, control, ExNeu1
        - group 6: Mathys, control, ExNeu2
        - group 7: Mathys, AD, ExNeu1
        - group 8: Mathys, AD, ExNeu2

        - group 9: Biogen, control, ExNeu1
        - group 10: Biogen, control, ExNeu2
        - group 11: Biogen, FTD, ExNeu1
        - group 12: Biogen, FTD, ExNeu2
    - by sample: 97 samples (for potential differential testing)
    - by barcode: ? cells (TBD)

- genes of interest: STMN2 and UNC13A


2025-09-24
----------

@Mira0507

- ``workflow/thalamus_excitatory/Snakefile`` in progress
    - barcode table prepared by concatenating two metadata tables
        - input:
            - ``input/thalamus_excitatory/ExNeu1_FOR_DEG.h5ad``
            - ``input/thalamus_excitatory/ExNeu1_FOR_DEG.h5ad``
        - output:
            - ``workflow/thalamus_excitatory/Snakefile/results/barcodes.tsv``
    - notes
        - ['D19-12386', 'D19-12393', 'Control10', 'FTLD-GRN1'] found 
          in none of the AnnData objects
        - the output barcode table consists of a total of 27,719 rows 


2025-09-26
----------

@Mira0507

- recopied input bam files
    - conda env: ``env``
    - updated script: ``workflow/thalamus_excitatory/prep_input.py``
    - note: file names changed from 
      ``<sampleid>_<study-disease>_possorted_genome_bam.bam`` to
      ``<sampleid>_possorted_genome_bam.bam``

- update ``workflow/thalamus_excitatory/Snakefile``
    - conda env: ``env``
    - notes
        - paths to input bam files added to the metadata table 
          (``workflow/thalamus_excitatory/results/barcodes.tsv``)
        - rules ``create_header`` and ``prep_sam`` added
        - filter barcodes based on *CR* (raw) or *CB* (corrected)? 

        .. code-block:: python

            # Leave or remove `-1` suffix from the barcodes extracted from AnnData
            barcode_list_cr = [f"CR:Z:{barcode.strip('-1')}" for barcode in set(df.barcode)]
            barcode_list_cb = [f"CB:Z:{barcode.strip()}" for barcode in set(df.barcode)]
            # Join barcodes into a single string for each option
            cr_joined = "|".join(barcode_list_cr)
            cb_joined = "|".join(barcode_list_cb)

            b_dic = {'cr': cr_joined, 'cb': cb_joined}
            for key, value in c_dic.items():

                # Create filtered sam files
                cmd = [
                    "samtools view ",
                    bam,
                    " | grep -E '",
                    value,
                    "'", 
                    f" >> temp_{key}.sam"
                ]
                cmd = "".join(cmd)
                subprocess.run(cmd, shell=True)

                # Sort
                cmd = [
                    f"grep -v '^@' temp_{key}.sam | ",
                    f"sort > {key}_sorted.txt"
                ]
                subprocess.run(cmd, shell=True)


            # Filter reads unique to each sam file
            cmd1 = "comm -23 cr_sorted.txt cb_sorted.txt > unique_to_cr.txt"
            cmd2 = "comm -13 cr_sorted.txt cb_sorted.txt > unique_to_cb.txt"
            subprocess.run(cmd1, chell=True)
            subprocess.run(cmd2, chell=True)

            # $ ll | grep unique_to
            # -rw-rw---- 1 sohnm CARD_MPU    0 Sep 26 21:22 unique_to_cr.txt
            # -rw-rw---- 1 sohnm CARD_MPU 7.0K Sep 26 21:22 unique_to_cb.txt

        - more reads were filtered by *CB*. this is because *CR* only includes
          exact matches and any reads with a hamming distance of 1 or larget
          are not included. in contrast, *CB* includes reads after seq-error
          correction performed by cellranger in addition to exact matches.
          this results in more reads filtered using the *CB* field.



2025-09-29
----------

@Mira0507

- update ``workflow/thalamus_excitatory/Snakefile``
    - conda env: ``env``
    - notes
        - runtime changed from 30min to 4hrs for the rule ``prep_sam``
        - rules ``prep_sam`` and ``create_sample_bam`` ran on every 20 samples 
          to make sure these rules run error-free and temporary files 
          ranging from 20G to 50G don't exceed the quota
        - rules ``create_sample_bam`` added
        - rules ``create_group_bam`` added

- reference transcriptome for all cellranger runs:
  ``/fdb/cellranger/refdata-cellranger-2024-A/refdata-gex-GRCh38-2024-A``


2025-09-30
----------

@Mira0507

- update ``workflow/thalamus_excitatory/Snakefile``
    - conda env: ``env``
    - notes
        - add celltype-specific filtering 
          (e.g. ``"{wildcards.sample}_{wildcards.celltype}.sam"``)
        - rerun every 30 samples per Snakemake job submission
        - STMN2 nor UNC13A undetected in 3733-T_ExNeu1. This issue results in an error
          when generating the sam file. Therefore, header is added to the sam file 
          at the ``prep_sam`` rule, as shown below:

        .. code-block:: bash

            samtools view ../../input/thalamus_excitatory/bam/3733-T_possorted_genome_bam.bam |
                grep -E "CB:Z:ATTCCATTCAGGGTAG-1" > 3733-T_ExNeu1_temp.sam || true
            cat results/header.sam > 3733-T_ExNeu1.sam
            cat 3733-T_ExNeu1_temp.sam | grep -E "GN:Z:STMN2|GN:Z:UNC13A" >> 3733-T_ExNeu1.sam || true
            rm 3733-T_ExNeu1_temp.sam

            # NOTE: || true ensures that Snakemake doesn't raise an error even if the `grep`
            # command returns nothing


2025-10-01
----------

@Mira0507

- update ``workflow/thalamus_excitatory/Snakefile``
    - conda env: ``env``
    - notes
        - function ``list_group_sams`` updated
        - rule ``create_group_celltype_bam`` updated


2025-10-03
----------

@Mira0507

- update ``workflow/thalamus_excitatory/Snakefile``
    - conda env: ``env``
    - notes
        - rule ``create_group_celltype_bam`` updated
        - subset of output bam files inspected using IGV

- update ``README.md``


2025-10-06
----------

@Mira0507

- run ``workflow/thalamus_excitatory/Snakefile``
  on all input files
    - conda env: ``env``
    - dryrun:

    .. code-block:: bash

        Job stats:
        job                           count
        --------------------------  -------
        all                               1
        create_group_celltype_bam        12
        create_header                     1
        create_sample_celltype_bam      186
        prep_sam                        186
        total                           386

    - all input samples

    .. code-block:: bash

        # Print all unique sample identifiers
        $ cut -f 1 workflow/thalamus_excitatory/results/barcodes.tsv | sort | uniq | head -n 3
        3075-T
        3514-T
        3549-T

        # Print the number of unique identifiers. This returns the column name. 
        $ cut -f 1 workflow/thalamus_excitatory/results/barcodes.tsv | sort | uniq | wc -l
        94  # A total of 93 input samples analyzed except for column name...

- add rule ``extract_junctions`` to ``workflow/thalamus_excitatory/Snakefile``
    - this rule runs the ``regtools junctions extract`` command that extracts 
      splicing junctions
    - references
        - https://regtools.readthedocs.io/en/latest/commands/junctions-extract/
        - https://www.nature.com/articles/s41467-023-37266-6

- update ``README.md``


2025-10-07
----------

@Mira0507

- install ``regtools`` in ``env`` (w ``--freeze-installed`` parameter)

- splicing junctions extracted on all samples 
    - script: ``workflow/thalamus_excitatory/Snakefile``
    - conda env: ``env``
    - input: output of the ``create_sample_celltype_bam`` rule
    - output: 
      ``workflow/thalamus_excitatory/results/bed/sample/<sample>_<celltype>_regtools.bed``

    .. code-block:: bash

        $ head -n 4 workflow/thalamus_excitatory/results/bed/sample/D19-12358_ExNeu2_regtools.bed
        chr19   17601563        17602859        JUNC00000001    2       ?       17601563        17602859        255,0,0 2       35,56   0,1240
        chr19   17609229        17612687        JUNC00000002    2       ?       17609229        17612687        255,0,0 2       20,76   0,3382
        chr19   17619641        17619962        JUNC00000003    1       ?       17619641        17619962        255,0,0 2       46,45   0,276
        chr19   17620638        17621837        JUNC00000004    2       ?       17620638        17621837        255,0,0 2       84,6    0,1193

    - notes
        - this step is conducted by the ``extract_junctions`` rule
        - some input bam files have no reads. such input files return
          empty output bed.
        - output bed columns
            - chrom: chromosome 
            - chromStart: The starting position of the junction-anchor. This includes 
              the maximum overhang for the junction on the left. For the exact junction 
              start add blockSizes[0].
            - chromEnd: The ending position of the junction-anchor. This includes 
              the maximum overhang for the juncion on the left. For the exact 
              junction end subtract blockSizes[1].
            - name: The name of the junctions, the junctions are just numbered JUNC1 to 
              JUNCn.
            - score: The number of reads supporting the junction.
            - strand: Defines the strand - either '+' or '-'. This is calculated 
              using the XS tag in the BAM file.
            - thickStart: Same as chromStart
            - thinkEnd: Same as chromEnd
            - itemRgb: RGB value ("255,0,0" by default)
            - blockCount: The number of blocks, 2 by default.
            - blockSize: A comma-separated list of the block sizes. The number of items 
              in this list should correspond to blockCount.
            - blockStart: A comma-separated list of block starts. All of the blockStart
              positions should be calculated relative to chromStart. The number of items
              in this list should correspond to blockCount.

- prep IGV session files per group in ``../shipped/igv``

    .. code-block:: bash

        $ tree
        .
        ├── all.xml
        ├── biogen.xml
        ├── exneu1.xml
        ├── exneu2.xml
        ├── marsan.xml
        ├── mathys.xml
        └── results
            └── bam
                └── group
                    ├── AD-Mathys_ExNeu1_sorted.bam
                    ├── AD-Mathys_ExNeu1_sorted.bam.bai
                    ├── AD-Mathys_ExNeu2_sorted.bam
                    ├── AD-Mathys_ExNeu2_sorted.bam.bai
                    ├── Control-Biogen_ExNeu1_sorted.bam
                    ├── Control-Biogen_ExNeu1_sorted.bam.bai
                    ├── Control-Biogen_ExNeu2_sorted.bam
                    ├── Control-Biogen_ExNeu2_sorted.bam.bai
                    ├── Control-Marsan_ExNeu1_sorted.bam
                    ├── Control-Marsan_ExNeu1_sorted.bam.bai
                    ├── Control-Marsan_ExNeu2_sorted.bam
                    ├── Control-Marsan_ExNeu2_sorted.bam.bai
                    ├── Control-Mathys_ExNeu1_sorted.bam
                    ├── Control-Mathys_ExNeu1_sorted.bam.bai
                    ├── Control-Mathys_ExNeu2_sorted.bam
                    ├── Control-Mathys_ExNeu2_sorted.bam.bai
                    ├── FTD-Biogen_ExNeu1_sorted.bam
                    ├── FTD-Biogen_ExNeu1_sorted.bam.bai
                    ├── FTD-Biogen_ExNeu2_sorted.bam
                    ├── FTD-Biogen_ExNeu2_sorted.bam.bai
                    ├── FTD-Marsan_ExNeu1_sorted.bam
                    ├── FTD-Marsan_ExNeu1_sorted.bam.bai
                    ├── FTD-Marsan_ExNeu2_sorted.bam
                    └── FTD-Marsan_ExNeu2_sorted.bam.bai

    - note
        - IGV session files were prepared locally and copied to HPC
        - Relative paths were used to navigate input file paths in each IGV session file

        .. code-block:: bash

            # Relative paths point to input files
            $ head -n 8 marsan.xml
            <?xml version="1.0" encoding="UTF-8" standalone="no"?>
            <Session genome="hg38" locus="chr8:79598936-79685953" nextAutoscaleGroup="2" version="8">
                <Resources>
                    <Resource path="results/bam/group/Control-Marsan_ExNeu2_sorted.bam" type="bam"/>
                    <Resource path="results/bam/group/FTD-Marsan_ExNeu1_sorted.bam" type="bam"/>
                    <Resource path="results/bam/group/Control-Marsan_ExNeu1_sorted.bam" type="bam"/>
                    <Resource path="results/bam/group/FTD-Marsan_ExNeu2_sorted.bam" type="bam"/>
                </Resources>

        - files
            - ``all.xml``: all groups and celltypes
            - ``exneu1.xml``: ExNeu1 in all groups
            - ``exneu2.xml``: ExNeu2 in all groups
            - ``biogen.xml``: samples from biogen, both celltypes
            - ``marsan.xml``: samples from marsan, both celltypes
            - ``mathys.xml``: samples from mathys, both celltypes


2025-10-08
----------

@Mira0507

- install leafcutter
    - reference: https://davidaknowles.github.io/leafcutter/articles/Installation.html
    - conda env
        - recipe: ``lc_requirements.txt``
        - env: ``lc_env``
        - installing the R package: 

        .. code-block:: bash

            # Ensure to have ``lc_env`` activated
            $ mamba install -c davidaknowles r-leafcutter --freeze-installed``

        - env exported: ``lc_env.archived.yaml``

- add rule ``prep_juncfiles`` to ``workflow/thalamus_excitatory/Snakefile``


