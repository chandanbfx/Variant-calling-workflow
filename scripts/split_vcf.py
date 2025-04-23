#!/usr/bin/env python
"""
Enhanced VCF splitting script with:
- Better error handling
- Parallel processing
- Gzip support
"""

import gzip
import sys
import os
from pathlib import Path
import subprocess

def split_vcf(input_vcf, reference, output_dir):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    try:
        # Get sample names
        samples = subprocess.check_output(
            ["bcftools", "query", "-l", input_vcf],
            stderr=subprocess.PIPE, text=True
        ).strip().split('\n')
        
        # Process each sample in parallel
        procs = []
        for sample in samples:
            out_vcf = f"{output_dir}/{sample}.vcf.gz"
            cmd = [
                "gatk", "SelectVariants",
                "-R", reference,
                "-V", input_vcf,
                "-sn", sample,
                "-O", out_vcf
            ]
            procs.append(subprocess.Popen(cmd))
        
        # Wait for all to complete
        for p in procs:
            p.wait()
            
    except Exception as e:
        sys.stderr.write(f"ERROR: {str(e)}\n")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: split_vcf.py <input.vcf> <reference.fa> <output_dir>")
    split_vcf(sys.argv[1], sys.argv[2], sys.argv[3])
