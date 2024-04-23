FASTAFileFix formats the DNA strings of a .fasta file into records of DNA 80 characters or less.

Currently the input file is a hard coded name "final_set_of_exons.fasta".
          the output file is a hard coded name "final_set_of_exons_formatted.fasta"

HybPiper list fasta files with the content formatted as:
- Header records  : <= 25 characters
- DNA strings     : records of <= 80 characters

For each concatenated DNA string from the inputfile two steps are performed prior to writing a properly formatted gene into the outfile file:
  1. The DNA string is scrubbed to remove all line break characters (end of line, line feeds, set)
  2. The DNA string is checked to make sure the only valid letters used are: "ACGTURYSWKMBDHVN "

Upon completion of the script, the output file can them be used with the hybpiper check_targetfile utility.
