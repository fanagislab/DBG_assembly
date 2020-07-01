1. Program I:  clean_lowqual
   
This program trims the low quality parts in a reads, and output a largest block in which the average error rate is lower than a given cutoff. The illumina sequencing machine produces reads with average error rate of 1%, however, the error bases are not distributed evenly. By set a cutoff (-e) of 0.1% for the average error rate, we can filter out most of (90%) of the sequencing errors. The program also has a function to filter the extreme short reads by using a cutoff (-r); The program runs in a multiple thread mode (-t). 

The input file should be fastq format, and there are two result files: one is the high-quality reads file, and the other one is a statistics file.

More help information can be obtained by typing "clean_lowqual" or "clean_lowqual -h" in the command line.


2. Program II: clean_adapter

Besides quality checking, another important task in QC is to detect the contaminations, often including the adapters introduced in the libray construction process. This program is used to detect and filter the adapter contaminations, however, it can also be used to filter other types of contaminations, by changing the default contamination file (-a) with user-specified file. The alignment algorithm between reads and adapters used ungapped dynamic programming local alignment, using a score matrix [match:1; mismatch:-2], and reports only the best hit, the aligning score cutoff can be set by a parameter (-s). The program runs in a multiple thread mode (-t). 

The input file should be fastq format, and there are two result files: one is the clean reads file, and the other one is a statistics file.

More help information can be obtained by typing "clean_adapter" or "clean_adapter -h" in the command line.

3. Furthermore

Both these two programs are designed to trim a read, and return a clean or high-quality part as the result;  However, for the most stringency, you may want to filtered out a read if it is trimmed for any reason. This can be achieved by just setting the minimum read length (-r) parameter to the original read length.

Both these two programs do not require much computer memory, the speed depends on the amount of contaminations and error-bases in the sequencing reads.








