import numpy as np
import pandas as pd
import tempfile
import os
import gc
import pysam
 
# Load BAM file
bamfile = pysam.AlignmentFile("./SRR_final_sorted.bam", "rb")
 
# Calculate read depth
def calculate_read_depth(bamfile, output_file):
    read_depth = []
    for pileupcolumn in bamfile.pileup():
        read_depth.append([bamfile.get_reference_name(pileupcolumn.reference_id), pileupcolumn.reference_pos, pileupcolumn.n])
    read_depth_df = pd.DataFrame(read_depth, columns=["chrom", "position", "depth"])
    read_depth_df.to_csv(output_file, index=False)
    del read_depth, read_depth_df
    gc.collect()
 
# Detect split reads
def detect_split_reads(bamfile, output_file):
    split_reads = []
    for read in bamfile.fetch():
        if read.has_tag('SA'):  # SA tag indicates supplementary alignments
            sa_tag = read.get_tag('SA')
            sa_parts = sa_tag.split(',')
            if len(sa_parts) > 3:
                end = int(sa_parts[1])
                svlen = abs(end - read.reference_start)
                if svlen > 50:
                    split_reads.append({
                        "chrom": bamfile.get_reference_name(read.reference_id),
                        "position": read.reference_start,
                        "end": end,
                        "svlen": svlen,
                        "sv_type": "INV"
                    })
    split_reads_df = pd.DataFrame(split_reads)
    split_reads_df.to_csv(output_file, index=False)
    del split_reads, split_reads_df
    gc.collect()
 
# Detect discordant reads
def detect_discordant_reads(bamfile, output_file):
    discordant_reads = []
    for read in bamfile.fetch():
        if read.is_paired and not read.is_proper_pair:
            end = read.next_reference_start
            svlen = abs(end - read.reference_start)
            if svlen > 50:
                discordant_reads.append({
                    "chrom": bamfile.get_reference_name(read.reference_id),
                    "position": read.reference_start,
                    "end": end,
                    "svlen": svlen,
                    "sv_type": "TRA"
                })
    discordant_reads_df = pd.DataFrame(discordant_reads)
    discordant_reads_df.to_csv(output_file, index=False)
    del discordant_reads, discordant_reads_df
    gc.collect()
 
# Filter structural variants by quality
def filter_svs(input_file, output_file, quality_threshold=50):
    svs_df = pd.read_csv(input_file)
    filtered_svs_df = svs_df[svs_df['svlen'] > quality_threshold]
    filtered_svs_df.to_csv(output_file, index=False)
    del svs_df, filtered_svs_df
    gc.collect()
 
def detect_deletions(read_depth_file, output_file, threshold=0.5, min_length=50):
    read_depth_df = pd.read_csv(read_depth_file)
    deletions = []
    mean_depth = read_depth_df['depth'].mean()
    current_region = None
 
    for index, row in read_depth_df.iterrows():
        if row['depth'] < mean_depth * threshold:
            if current_region is None:
                current_region = {
                    "chrom": row['chrom'],
                    "start": row['position'],
                    "end": row['position'],
                    "depth": row['depth']
                }
            else:
                current_region['end'] = row['position']
                current_region['depth'] += row['depth']
        else:
            if current_region is not None:
                svlen = current_region['end'] - current_region['start'] + 1
                if svlen > min_length:
                    deletions.append({
                        "chrom": current_region['chrom'],
                        "position": current_region['start'],
                        "sv_type": "DEL",
                        "end": current_region['end'],
                        "svlen": -svlen
                    })
                current_region = None
 
    if current_region is not None:
        svlen = current_region['end'] - current_region['start'] + 1
        if svlen > min_length:
            deletions.append({
                "chrom": current_region['chrom'],
                "position": current_region['start'],
                "sv_type": "DEL",
                "end": current_region['end'],
                "svlen": -svlen
            })
    
    deletions_df = pd.DataFrame(deletions)
    deletions_df.to_csv(output_file, index=False)
    del read_depth_df, deletions, deletions_df
    gc.collect()
 
def detect_duplications(read_depth_file, output_file, threshold=1.5, min_length=50):
    read_depth_df = pd.read_csv(read_depth_file)
    duplications = []
    mean_depth = read_depth_df['depth'].mean()
    current_region = None
 
    for index, row in read_depth_df.iterrows():
        if row['depth'] > mean_depth * threshold:
            if current_region is None:
                current_region = {
                    "chrom": row['chrom'],
                    "start": row['position'],
                    "end": row['position'],
                    "depth": row['depth']
                }
            else:
                current_region['end'] = row['position']
                current_region['depth'] += row['depth']
        else:
            if current_region is not None:
                svlen = current_region['end'] - current_region['start'] + 1
                if svlen > min_length:
                    duplications.append({
                        "chrom": current_region['chrom'],
                        "position": current_region['start'],
                        "sv_type": "DUP",
                        "end": current_region['end'],
                        "svlen": svlen
                    })
                current_region = None
 
    if current_region is not None:
        svlen = current_region['end'] - current_region['start'] + 1
        if svlen > min_length:
            duplications.append({
                "chrom": current_region['chrom'],
                "position": current_region['start'],
                "sv_type": "DUP",
                "end": current_region['end'],
                "svlen": svlen
            })
 
    duplications_df = pd.DataFrame(duplications)
    duplications_df.to_csv(output_file, index=False)
    del read_depth_df, duplications, duplications_df
    gc.collect()
 
def create_vcf(input_files, output_file):
    with open(output_file, 'w') as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=CustomSVDetector\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for input_file in input_files:
            svs_df = pd.read_csv(input_file)
            for _, sv in svs_df.iterrows():
                svlen = abs(sv['svlen'])
                vcf.write(f"{sv['chrom']}\t{sv['position']}\t.\tN\t<{sv['sv_type']}:{svlen}>\t.\tPASS\tSVTYPE={sv['sv_type']};END={sv['end']};SVLEN={sv['svlen']}\n")
            del svs_df
            gc.collect()
 
# Temporary file paths
read_depth_file = tempfile.mktemp()
split_reads_file = tempfile.mktemp()
discordant_reads_file = tempfile.mktemp()
filtered_split_reads_file = tempfile.mktemp()
filtered_discordant_reads_file = tempfile.mktemp()
deletions_file = tempfile.mktemp()
duplications_file = tempfile.mktemp()
 
# Calculate read depth and save to temp file
calculate_read_depth(bamfile, read_depth_file)
 
# Detect split reads and save to temp file
detect_split_reads(bamfile, split_reads_file)
 
# Detect discordant reads and save to temp file
detect_discordant_reads(bamfile, discordant_reads_file)
 
# Filter split reads and save to temp file
filter_svs(split_reads_file, filtered_split_reads_file)
 
# Filter discordant reads and save to temp file
filter_svs(discordant_reads_file, filtered_discordant_reads_file)
 
# Detect deletions and save to temp file
detect_deletions(read_depth_file, deletions_file)
 
# Detect duplications and save to temp file
detect_duplications(read_depth_file, duplications_file)
 
# Create the VCF file
create_vcf([filtered_split_reads_file, filtered_discordant_reads_file, deletions_file, duplications_file], "structural_variants.vcf")
 
# Clean up temporary files
os.remove(read_depth_file)
os.remove(split_reads_file)
os.remove(discordant_reads_file)
os.remove(filtered_split_reads_file)
os.remove(filtered_discordant_reads_file)
os.remove(deletions_file)
os.remove(duplications_file)
 
print("Structural variant detection completed. VCF file generated.")