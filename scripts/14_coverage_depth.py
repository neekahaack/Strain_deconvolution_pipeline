import subprocess
from pathlib import Path

bowtie_simulation_dir = Path("/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/bowtie_alignment_simulation")
bowtie_raw_dir = Path("/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/bowtie_alignment_raw")
output_file = Path("/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/depth_of_coverage.csv")

# Write header
with open(output_file, "w") as out:
    out.write("sample,avg_depth,coverage_breadth\n")


def compute_depth_and_breadth(bam_path):
    """
    Computes:
    - average depth (mean coverage)
    - coverage breadth = fraction of positions with depth > 0
    """

    # Average depth command
    avg_cmd = (
        f"samtools depth -a {bam_path} | "
        "awk '{{sum+=$3}} END {{if (NR>0) print sum/NR; else print 0}}'"
    )
    avg_res = subprocess.run(avg_cmd, shell=True, capture_output=True, text=True)
    avg_depth = float(avg_res.stdout.strip())

    # Coverage breadth command (≥1x coverage)
    breadth_cmd = (
        f"samtools depth -a {bam_path} | "
        "awk '{if ($3>0) covered+=1; total+=1} END {if (total>0) print covered/total; else print 0}'"
    )
    breadth_res = subprocess.run(breadth_cmd, shell=True, capture_output=True, text=True)
    breadth = float(breadth_res.stdout.strip())

    return avg_depth, breadth


# Process simulated BAMs
for bam in bowtie_simulation_dir.glob("*.bam"):
    avg, breadth = compute_depth_and_breadth(bam)
    with open(output_file, "a") as out:
        out.write(f"{bam.stem}_sim,{avg:.3f},{breadth:.4f}\n")


# Process raw BAMs
for bam in bowtie_raw_dir.glob("*.bam"):
    avg, breadth = compute_depth_and_breadth(bam)
    with open(output_file, "a") as out:
        out.write(f"{bam.stem}_raw,{avg:.3f},{breadth:.4f}\n")

print("Done.")
