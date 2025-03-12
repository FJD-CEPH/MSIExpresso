MSIExpresso
A tool for RNASeq analysis of transcriptional events related to microsatellite instability and MSI status determination.
Requirements:
- Perl
- BWA
- samtools
- bcftools
- MSIExpresso reference database available here: https://fjd-ceph.box.com/shared/static/0v3p4wjx5h0usdrej7jeod6o42ivqon9.zip
Running MSIExpresso:
<path to>/MSIExpresso.pl <configuration file> [updateStatus]
The configuration file contains all the information needed to run one sample analysis.
Parameters are in the format key/value as follow : <parameter name>=<parameter value>
One parameter per line.
List of all parameters available:
Global parameters:
database dir: the location of the fasta database provided with MSIExpresso. (EG: <your path>/MSIExpresso/db).
output dir: the location for sample analysis output files (must be created before running MSIExpresso).
output files root: base name for each output file (e.g: the sample name)
path to bwa: path to bwa executable (version 0.7.15 or higher recommended) if not in the default path
bwa threads: number of threads for running bwa path to samtools: path to samtools executable (version 1.3.1 or higher recommended) if not in the default path
Analysis parameters:
For detailed analysis parameters see MSIExpressoDocumentation.txt on GitHub
