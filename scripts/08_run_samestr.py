import subprocess
from pathlib import Path
import shutil


def run_samestr_per_metagenome(threads=8):

    # --------------------------------------------------------
    # PATHS
    # --------------------------------------------------------
    bowtie_base = Path(
        "/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/bowtie_alignment"
    )

    samestr_db = Path(
        "/g/typas/Personal_Folders/Neeka/Model/data/09_samestr/samestr_db"
    )

    output_base = Path(
        "/g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/samestr_results"
    )
    output_base.mkdir(exist_ok=True, parents=True)

    # sanity checks
    if not samestr_db.exists():
        raise RuntimeError(f"SameStr marker DB not found at: {samestr_db}")

    # --------------------------------------------------------
    # LOOP THROUGH METAGENOMES
    # --------------------------------------------------------
    for meta_dir in sorted(bowtie_base.iterdir()):
        if not meta_dir.is_dir():
            continue

        metagenome_name = meta_dir.name
        print(f"\n=== Processing SameStr for metagenome: {metagenome_name} ===")

        # paths
        convert_dir = output_base / metagenome_name / "convert"
        merge_dir = output_base / "merge"
        filter_dir = output_base / "filter"
        compare_dir = output_base / "compare"
        
        

        convert_dir.mkdir(parents=True, exist_ok=True)
        merge_dir.mkdir(parents=True, exist_ok=True)
        filter_dir.mkdir(parents=True, exist_ok=True)
        compare_dir.mkdir(parents=True, exist_ok=True)
        
        

        # ---------- CLEAN OLD RESULTS ----------
        for f in convert_dir.iterdir():
            if f.is_file() or f.is_symlink():
                f.unlink()
            elif f.is_dir():
                shutil.rmtree(f)

        for f in compare_dir.iterdir():
            if f.is_file() or f.is_symlink():
                f.unlink()
            elif f.is_dir():
                shutil.rmtree(f)

        for f in merge_dir.iterdir():
            if f.is_file() or f.is_symlink():
                f.unlink()
            elif f.is_dir():
                shutil.rmtree(f)
        # ----------------------------------------

        # --------------------------------------------------------
        # FILTER OUT ALL INVALID FILE TYPES
        # --------------------------------------------------------
        invalid_endings = (
            ".sam", ".sam.gz", ".sam.bz2",
            ".bowtie2", "metaphlan.bowtie2.bz2",
            ".mapout"
        )

        def is_invalid(f):
            return any(str(f).endswith(end) for end in invalid_endings)

        candidate_files = [
            f for f in meta_dir.iterdir()
            if not is_invalid(f)
        ]

        # --------------------------------------------------------
        # FIND THE ONE VALID BAM FILE
        # --------------------------------------------------------
        bam_files = [
            f for f in candidate_files
            if f.name == f"{metagenome_name}.bam"
        ]

        if len(bam_files) == 0:
            print(f"No BAM found for {metagenome_name} (after filtering).")
            continue
        if len(bam_files) > 1:
            print(f"Multiple BAMs found for {metagenome_name}:")
            for b in bam_files:
                print("   •", b)
            print("Aborting this sample.")
            continue

        bam_file = bam_files[0]
        print(f"Found BAM: {bam_file}")

        # --------------------------------------------------------
        # SAMESTR CONVERT
        # --------------------------------------------------------
        print("Running samestr convert...")

        convert_cmd = [
            "samestr", "convert",
            "--input-files", str(bam_file),
            "--marker-dir", str(samestr_db),
            "--nprocs", str(threads),
            "--min-vcov", "1",
            "--output-dir", str(convert_dir)
        ]

        print("COMMAND:", " ".join(convert_cmd))
        subprocess.run(convert_cmd, check=True)
        return()


        # --------------------------------------------------------
        # SAMESTR MERGE  
        # --------------------------------------------------------
        print("Running samestr merge...")

        # Collect converted npz files for this sample
        npz_files_convert = sorted(list((output_base).glob("**/*.npz")))
        npz_files_convert = [f for f in npz_files_convert if 'convert' in str(f)]
        if not npz_files_convert:
            print(f"merge failed — no .npz output from convert for {metagenome_name}.")
            continue

        merge_cmd = [
            "samestr", "merge",
            "--input-files",
        ] + [str(f) for f in npz_files_convert] + [
            "--marker-dir", str(samestr_db),
            "--nprocs", str(threads),
            "--output-dir", str(merge_dir),
        ]

        print("COMMAND:", " ".join(merge_cmd))
        subprocess.run(merge_cmd, check=True)


        # --------------------------------------------------------
        # SAMESTR FILTER
        # --------------------------------------------------------
        print("Running samestr filter...")

        # Filter expects .npz + .names.txt from MERGE
        npz_files_merge = sorted(list(merge_dir.rglob("*.npz")))
        name_files_merge = sorted(list(merge_dir.rglob("*.names.txt")))

        if not npz_files_merge:
            print(f"filter failed — no .npz output from merge for {metagenome_name}.")
            continue

        if not name_files_merge:
            print(f"filter failed — no .names.txt output from merge for {metagenome_name}.")
            continue

        filter_cmd = [
            "samestr", "filter",
            "--keep-mono",
            "--input-files",
        ] + [str(f) for f in npz_files_merge] + [
            "--input-names",
        ] + [str(n) for n in name_files_merge] + [
            "--marker-dir", str(samestr_db),
            "--nprocs", str(threads),
            "--output-dir", str(filter_dir),
        ]

        print("COMMAND:", " ".join(filter_cmd))
        subprocess.run(filter_cmd, check=True)


        # --------------------------------------------------------
        # SAMESTR COMPARE
        # --------------------------------------------------------
        print("Running samestr compare...")

        # Compare expects .npz + .names.txt from FILTER
        npz_files_filter = sorted(list(filter_dir.rglob("*.npz")))
        name_files_filter = sorted(list(filter_dir.rglob("*.names.txt")))

        if not npz_files_filter:
            print(f"compare failed — no .npz output from filter for {metagenome_name}.")
            continue

        if not name_files_filter:
            print(f"compare failed — no .names.txt output from filter for {metagenome_name}.")
            continue

        compare_cmd = [
            "samestr", "compare",
            "--input-files",
        ] + [str(f) for f in npz_files_filter] + [
            "--input-names",
        ] + [str(n) for n in name_files_filter] + [
            "--marker-dir", str(samestr_db),
            "--nprocs", str(threads),
            "--output-dir", str(compare_dir)
        ]

        print("COMMAND:", " ".join(compare_cmd))
        breakpoint()
        subprocess.run(compare_cmd, check=True)

        print(f"Done: SameStr results written to: {compare_dir}")


if __name__ == "__main__":
    run_samestr_per_metagenome(threads=8)


#FASTQ → Bowtie2 (SameStr markers) → BAM → samestr convert

# what i did to get databases
# wget https://cmprod1.cibio.unitn.it/databases/MetaPhlAn/mpa_vJun23_CHOCOPhlAnSGB_202307.tar
# unzipped it in the db_MetaPhlAn folder
# tar -xvf mpa_vJun23_CHOCOPhlAnSGB_202307.tar
# then concatenated SGB and VSG bz2 files into one file each (as per MetaPhlAn instructions)
# > mpa_vJun23_CHOCOPhlAnSGB_202307.fna.bz2

# next: building the SameStr database using the MetaPhlAn markers

#samestr db \
#--markers-info db_MetaPhlAn/mpa_vJun23_CHOCOPhlAnSGB_202307.pkl \
#--markers-fasta db_MetaPhlAn/mpa_vJun23_CHOCOPhlAnSGB_202307.fna.bz2 \
#--db-version db_MetaPhlAn/mpa_latest \
#--output-dir samestr_db/

# this created:
#   db_markers/                 (per-clade marker FASTAs)
#   db_taxonomy.tsv             (taxonomic mapping)
#   db_clades.json              (clade definitions)
#   db_manifest.json            (database metadata)
#   samestr_markers.fna         (all markers concatenated)      #actually built this with function in align script
#   command_log.json            (SameStr build log)

# 2) build Bowtie2 index for SameStr marker database
#cd /g/.../09_samestr/samestr_db

#bowtie2-build --large \
  #samestr_markers.fna \
  #samestr_marker_index

# this generated:
#   samestr_marker_index.1.bt2l
#   samestr_marker_index.2.bt2l
#   samestr_marker_index.3.bt2l
#   samestr_marker_index.4.bt2l
#   samestr_marker_index.rev.1.bt2l
#   samestr_marker_index.rev.2.bt2l

# (completed the samestr_db directory and alignment index)
