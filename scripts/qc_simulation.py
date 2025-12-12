from Bio import SeqIO
import matplotlib.pyplot as plt


def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return 100.0 * (g + c) / len(seq)

def gc_distribution(fq_file, n=10000):
    gc_vals = []
    for i, record in enumerate(SeqIO.parse(fq_file, "fastq")):
        gc_vals.append(gc_content(str(record.seq)))
        if i >= n:  # limit to 10k reads for speed
            break
    return gc_vals

sim_gc = gc_distribution("simulated_reads/genome1_1.fq")
real_gc = gc_distribution("real_reads/genome1_1.fastq")

plt.hist(sim_gc, bins=40, alpha=0.5, label="Simulated")
plt.hist(real_gc, bins=40, alpha=0.5, label="Real")
plt.xlabel("GC content (%)")
plt.ylabel("Read count")
plt.legend()
plt.show()