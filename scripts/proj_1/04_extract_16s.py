# #need to isolate 16S sequence from whole sequence

# #in gff file there are 9 columns
# #first is contigs, there is also start and stop of gene (base #)
# #can do it this way: contig of 16S rRNA is determined
#     #contig is assembled from short reads and corresponds to continuous portion of gene
# #contig is used to find corresponding sequence in fasta doc
#     #.fnn has contigs and corresponding nucleotide sequence
# #so fasta file has to be loaded in
#     #can search whole file for rRNA name then take sequence after the name and until the next > symbol
# #then make new file that contains taxonomy info and 16s sequence

import os
from pathlib import Path

# seq = "TAACAGCTGGAAACGGCTGCTAATACCGCATAAGCGCACAGAATCGCATGATTCGGTGTGAAAAGCTCCGGCAGTATAGGATGGTCCCGCGTCTGATTAGCTGGTTGGCGGGGTAACGGCCCACCAAGGCGACGATCAGTAGCCGGCTTGAGAGAGTGGACGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAAAAAGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGTGCGTAGGTGGCATGGTAAGTCAGAAGTGAAAGCCCGGGGCTTAACCCCGGGACTGCTTTTGAAACTGTCATGCTGGAGTGCAGGAGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGACTGTCACTGACACTGATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATACTAGGTGTCGGGGCCGTAGAGGCTTCGGTGCCGCAGCAAACGCAGTAAGTATTCCACCTGGGGAGTACGTTCGCAAGAATGAAACTCAAAGGAATTGACGGGGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCTGGTCTTGACATCCCAATGACCGAACCTTAACCGGTTTTTTCTTTCGAGACATTGGAGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTATCTTTAGTAGCCAGCATTACGGATGGGCACTCTAGAGAGACTGCCAGGGATAACCTGGAGGAAGGTGGGGACGACGTCAAATCATCATGCCCCTTATGGCCAGGGCTACACACGTGCTACAATGGCGTAAACAAAGGGAAGCGAAGTCGTGAGGCGAAGCAAATCCCAGAAATAACGTCTCAGTTCGGATTGTAGTCTGCAACTCGACTACATGAAGCTGGAATCGCTAGTAATCGTGAATCAGAATGTCACGGTGAATACGTTCCCGGGTCTTGTACACACCGCCCGTCACACCATGGGAGTCAGTAACGCCCGAAGTCAGTGACCCAACCTTATAGGAGGGAGCTGCCGAAGGTGGGACCGATAACTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTT"

# len(seq)

# for kw in keyword:
#     pos = text.count(kw) 
#     if pos != -1:
#         snippet = text[pos+len(kw):pos+len(kw)+30]
#         print("kw found")
#         print(kw)
#         print("Next 30 chars:", snippet)
#     else:
#         print("not found")

# for files in os.listdir(dirname):
#     if files.endswith(ext):
#         print(files)  # printing file name of desired extension
#     else:
#         continue


# def process_folders(folder: Path):
#   for item in folder.iterdir():
#     if item.is_dir():
#       # recursion for subfolder
#       process_folder(item)
#     elif item.suffix == ".txt":
#       process_file(item)

# base_indir = "/g/typas/Personal_Folders/Neeka/Model/data/test_one_genome/PROKKA_09122025"
# outdir = "/g/typas/Personal_Folders/Neeka/Model/data/test_one_genome/16S"
#base = "*.ffn" #somehow need to go through indir and into subfolders and find .ffn files
#file = "/g/typas/Personal_Folders/Neeka/Model/data/one_genome/PROKKA_09122025/PROKKA_09122025.ffn"

# base = os.path.basename(base_indir)
# print(base)

#opens file in read mode, then reads content into variable, lines, as a list of strings. 
# with open ("/g/typas/Personal_Folders/Neeka/Model/data/test_one_genome/PROKKA_09122025/PROKKA_09122025.ffn", "r") as f:
#     lines = f.readlines()

# #looking for this keyword
# keyword = "16S ribosomal RNA"

# #loops through index and element (text found at that index)
# with open("/g/typas/Personal_Folders/Neeka/Model/data/test_one_genome/{base}.ffn", "w") as f:
#     for i, line in enumerate(lines):
#         #if kw and element match, print the element
#         if keyword in line:
#             #print(line)
#             print(line.strip())     #this strips away \n and any spaces
#             #sequence is set to nothing to start
#             sequence = ""
#             for seq in lines[i+1:]:           #lines(i+1:) is slice of lines after header(i=0), seq is one of thse slices
#                 #escapes loop when next contig begins
#                 if seq.startswith(">"):
#                     break
#                 sequence += seq.strip()
#             f.write(line + sequence)



#define a function
#required (positional) arguments go first
#optional arguments with defaults after
def extract_16S(
    input_fasta_file,
    output_base_path = '/g/typas/Personal_Folders/Neeka/Model/data/04_16s_sequences/',
    keyword = "16S ribosomal RNA"
):
    #opens file in read mode, then reads content into variable, lines, as a list of strings. 
    with open (input_fasta_file, "r") as f:
        lines = f.readlines()
    
    #os.path.basename takes file path and drops all upstream path info
    #split then splits at all ., and then we take all but the last slice (ffn) with the -1
    #then joing all the parts (besides last slice) that were separated at the . by slice
    base = '.'.join(os.path.basename(input_fasta_file).split('.')[:-1])

    #
    with open(f"{output_base_path}{base}__16S_genes.ffn", "w") as f:
        num_seq_found = 0
        for i, line in enumerate(lines):
            #if kw and element match, print the element
            if keyword in line:
                num_seq_found += 1
                #print(line)
                #print(line.strip())     #this strips away \n and any spaces
                #sequence is set to nothing to start
                sequence = ""
                for seq in lines[i+1:]:           #lines(i+1:) is slice of lines after header(i=0), seq is one of thse slices
                    #escapes loop when next contig begins
                    if seq.startswith(">"):
                        break
                    sequence += seq.strip()
                line_without_newline = line.rstrip()
                line_without_newline = line_without_newline.split(' ')[0]
                line_without_newline = line_without_newline + '_____' + base + '\n'
                f.write(line_without_newline + sequence + '\n')
        print(f'Found {num_seq_found} hits!')

# extract_16S(
#     input_fasta_file = "/g/typas/Personal_Folders/Neeka/Model/data/test_one_genome/PROKKA_09122025/PROKKA_09122025.ffn",
#     output_base_path = '/g/typas/Personal_Folders/Neeka/Model/data/test_one_genome/'    
# )
print(__name__)
breakpoint()         #after this will show pdb and that is debugger prompt, type c for continue
if __name__ == '__main__':
    start_dir = "/g/typas/Personal_Folders/Neeka/Model/data/03_annotation/"
    for dirpath, dirnames, filenames in os.walk(start_dir):
        for filename in filenames:
            if filename.endswith('.ffn') and not '16S' in filename:
                full_path = os.path.join(dirpath, filename)
                # Perform your operation here
                print(f"Operating on: {full_path}")
                extract_16S(
                    input_fasta_file = full_path
                )

