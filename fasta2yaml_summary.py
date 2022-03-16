#!/usr/bin/env python

import sys
import yaml

import fasta_utils

fasta_desc_lines = list()

with sys.stdin as fh:
    for line in fh:
        if fasta_utils.is_description_line(line):
            fasta_desc_lines.append(line)

summary = fasta_utils.summarise_species_protein_data(fasta_desc_lines)
yaml_text = yaml.dump(summary, explicit_start=True, default_flow_style=False)

with sys.stdout as fh:
    fh.write(yaml_text)
