# CSGanalysis
In this repository I provide perl code that I wrote for the analysis of orchid bee chemosensory gene (CSG) repertoires for a manuscript currently under review at BMC Evolutionary Biology [Brand P, Ramirez S, Leese F, Tollrian R & Eltz T - Rapid evolution of chemosensory receptor genes in a pair of sibling species of orchid bees (Apidae: Euglossini)]. The scripts were designed to run on a Mac and will probably run on other Unix-like OSs. PC is not supported but could be easily implemented.

#Prerequisites to Run the Pipeline

MySQL 5.6 (or higher) - https://dev.mysql.com/downloads/

blast+ suite          - http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

DBI (perl module)     - http://search.cpan.org/~timb/DBI/DBI.pm

A query database containing e.g. chemosensory genes of your choice

# Function of the Scripts
or_suchev1.1.pl 

This is the main script that coordinates the steps from saving your transcriptome dataset in MySQL to screen for chemosensory genes of your choice, saving them in a specific table structure (the script does NOT create these tables - they have to be present before running. You can install them running the create-table command I provide in the CreateTables.txt file). While running or_suchev1.1.pl makes use of multiple programs provided by the blast+ suite (see Get it to Work) as well as the orf-finder.pl script provided in this repository.


orf-finder.pl

This script is designed to find all possible open reading frames (ORFs) of a batch fasta file. It is able to detect ORFs on both strands. Amino acids detected are restricted to the standard set and can be extended easily. I set the script to only output ORFs of at least 150 amino acids. You can change that by setting $lens to your desired length.

NOTE: The script loads all input sequences into memory. Don't run the scripts on fasta files that are too big for your RAM.
The script is set to run blast with 4 threads with a minimum evalue of 10*e^-6. You can adjust these values in the code ($eval and $threads).


CreateTables.txt

Create-table syntax to create tables used by the or_suchev1.1.pl script. 

# Get it to Work
The scripts are designed to interact with a MySQL database. In order to get the scripts to run make sure to include your MySQL database name, password and socket in or_suchev1.1.pl I hardcoded a specific folder structure and database structure that have to be adjusted in or_suchev1.1.pl. All scripts should be saved together in a single folder for best performance. 


If you have any questions, comments or the need to make fun of my terrible hacking, feel free to write me an email to pbrand@ucdavis.edu

Philipp
