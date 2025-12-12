align_sequences = __import__('05_align_sequences')

if __name__ == '__main__':
    align_sequences.plot_distribution(
         align_sequences.sequence_file,
         align_sequences.histogram_file)