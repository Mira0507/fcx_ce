import re

logfile = "Snakefile.log"

results = []

with open(logfile) as f:
    in_rule = False
    bam_count = 0
    outfile = None

    for line in f:
        if line.startswith("rule aggr_bams_group_celltype:"):
            in_rule = True
            bam_count = 0
            continue

        if in_rule:

            if line.lstrip().startswith("input:"):
                bam_count += len(re.findall(r"\.bam\b", line))

            # input list ends before output section
            if line.lstrip().startswith("output:"):
                outfile = line.split(",", 1)[0].strip()
                results.append((outfile, bam_count))
                in_rule = False

for wc, count in results:
    print(f"{wc}\t{count}")
