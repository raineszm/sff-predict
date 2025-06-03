from snakemake.script import snakemake

with open(snakemake.output.headlines, "w") as f:
    for month_headlines in snakemake.input.months:
        with open(month_headlines, "r") as f:
            for line in f:
                f.write(line)
