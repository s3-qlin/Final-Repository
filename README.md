# Final-Repository
# Using BLAST KEY to Find Homologs
Initializes a BLAST protein database using the file allprotein.fas as input. The -dbtype prot flag specifies that the database contains protein sequences.
````
makeblastdb -in allprotein.fas -dbtype prot
````
Downloads the protein sequence with accession NP_000916.2 (Pyruvate Dehydrogenase E1 beta subunit) in FASTA format. The -F fasta flag ensures the output is in FASTA format.
````
ncbi-acc-download -F fasta -m protein "NP_000916.2"
````
Executes a BLASTP search to compare the query sequence (NP_000916.2.fa) against the database allprotein.fas. The -outfmt 0 flag specifies the default BLAST output format, while -max_hsps 1 ensures only the best alignment for each subject sequence is reported.
````
blastp -db ../allprotein.fas -query NP_000916.2.fa -outfmt 0 -max_hsps 1 -out pyruvate-dehydrogenase-complex.blastp.typical.out
````
Similar to the previous command but produces a tabular output format (-outfmt 6) with specific columns for sequence ID, percentage identity, alignment length, mismatches, gaps, e-value, bitscore, etc. This format is easier to parse programmatically.
````
blastp -db ../allprotein.fas -query NP_000916.2.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out pyruvate-dehydrogenase-complex.blastp.detail.out
````
Filters the BLAST results to retain only hits with an e-value less than 1e-30. The filtered sequence IDs are saved to a new file.
````
awk '{if ($6< 1e-30)print $1 }' pyruvate-dehydrogenase-complex.blastp.detail.out > pyruvate-dehydrogenase-complex.blastp.detail.filtered.out
````
Counts the number of lines in the filtered results file, which corresponds to the total number of sequences passing the e-value threshold.
````
wc -l pyruvate-dehydrogenase-complex.blastp.detail.filtered.out
````
Extracts and counts homologs for each species by matching patterns in sequence IDs, sorting them, and finding unique counts.
````
grep -o -E "^[A-Z]\.[a-z]+" pyruvate-dehydrogenase-complex.blastp.detail.filtered.out | sort | uniq -c
````
# Sequence Alignment of Gene Family
Extracts sequences matching the BLAST output file and removes any sequences containing "carpio," saving the filtered sequences into pyruvate-dehydrogenase-complex.homologs.fas.
````
seqkit grep --pattern-file ~/lab03-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.fas
````
Aligns all the homolog sequences to generate a global alignment file saved as pyruvate-dehydrogenase-complex.homologs.al.fas.
````
muscle -align ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.fas -output ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas
````
Opens the alignment in ALV for detailed visualization.
````
alv -kli ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas | less -RS
````
Displays the alignment, highlighting positions where the most common amino acid is present in at least 50% of sequences.
````
alv -kli --majority ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas | less -RS
````
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas
````
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas
````
Determines the total number of columns (alignment width) in the multiple sequence alignment.
````
alignbuddy -al ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas
````
Removes all columns containing gaps and calculates the remaining alignment width.
````
alignbuddy -trm all ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas | alignbuddy -al
````
Removes all invariant (conserved) columns and calculates the remaining alignment width.
````
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas | alignbuddy -al
````
Computes the average percent identity among all aligned sequences using T-Coffee.
````
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas -output sim
````
Computes the average percent identity among all aligned sequences using AlignBuddy and AWK for processing.
````
alignbuddy -pi ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF; i++){ sum+=$i; num++ } } END { print(100*sum/num) } '
````
# Gene Family Phylogeny Using IQ-TREE
Modifies the alignment file to replace spaces with underscores and removes sequences containing the duplicate label dupelabel. Outputs a cleaned alignment file in the Lab 5 directory.
````
sed 's/ /_/g'  ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.fas
````
Constructs a maximum likelihood phylogenetic tree using IQ-TREE. The -bb 1000 option computes ultrafast bootstrap values with 1000 replicates, and -nt 2 sets the number of threads for faster computation.
````
iqtree -s ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.fas -bb 1000 -nt 2 
````
Outputs the phylogenetic tree in ASCII format for quick inspection.
````
nw_display ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.fas.treefile
````
Generates an unrooted phylogenetic tree graphic as a PDF using an R script. Adjusts text label size (0.4) and truncates labels to 15 characters.
````
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R  ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.fas.treefile ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.fas.treefile.pdf 0.4 15
````
Roots the tree at the midpoint of the longest branch using gotree and saves the output in a new file.
````
gotree reroot midpoint -i ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile
````
Outputs the midpoint-rooted tree in ASCII format for quick inspection.
````
nw_order -c n ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile | nw_display -
````
Creates an SVG graphic of the midpoint-rooted tree with specified styling options.
````
nw_order -c n ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s >  ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile.svg -
````
Converts the SVG tree graphic into a PDF format for easier sharing and review.
````
convert  ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile.svg  ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile.pdf
````
Displays the tree as a cladogram (equal branch lengths) in SVG format.
````
nw_order -c n ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.midCl.treefile.svg -
````
Converts the cladogram SVG into a PDF for easier accessibility and inspection.
````
convert ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.midCl.treefile.pdf
````
# Reconcilation of Gene and Species Tree
Performs reconciliation of the gene tree with the species tree using Notung. The command specifies:
<br> -s: Species tree input file.
<br> -g: Gene tree input file.
<br> --reconcile: Perform reconciliation.
<br> --speciestag prefix: Specifies that species names are in the prefix of gene names.
<br> --savepng: Saves the reconciled tree as a PNG image.
<br> --events: Saves inferred duplication and loss events.
<br> --outputdir: Specifies the output directory.
````
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/pyruvate-dehydrogenase-complex/
````
Displays the species tree in ASCII format for verification.
````
nw_display ~/lab05-$MYGIT/species.tre
````
Extracts and displays internal nodes of the species tree from the Notung reconciliation output.
````
grep NOTUNG-SPECIES-TREE ~/lab06-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile.rec.ntg | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -
````
Converts the Notung reconciliation output into RecPhyloXML format, suitable for visualization and further analysis.
````
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile.rec.ntg --include.species
````
Generates a graphical representation of the reconciled gene tree within the species tree. Key options:
<br> -f: Input file in RecPhyloXML format.
<br> -o: Output file in SVG format.
````
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile.rec.ntg.xml -o ~/lab06-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile.rec.svg
````
Converts the SVG file to a high-resolution PDF for better viewing and sharing.
````
convert -density 150 ~/lab06-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile.rec.pdf
````
#Protein Domain Prediction
Removes asterisks (*), representing stop codons, from the input sequence file. The modified file is saved in the Lab 8 directory.
````
sed 's/*//' ~/lab04-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.fas > ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.fas
````
Runs RPS-BLAST to scan the protein sequences against the Pfam database. Outputs results in tabular format with six columns:
<br> Query sequence ID
<br> Query sequence length
<br> Start position of the match
<br> End position of the match
<br> E-value of the match
<br> Title of the matched domain
````
rpsblast -query ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
````
Copies the gene tree file created in Lab 5 to the Lab 8 directory for integration with protein domain annotations.
````
cp ~/lab05-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex
````
Uses an R script to produce a PDF combining the phylogenetic tree with predicted protein domains. The tree is displayed on the left, and the corresponding domains are aligned on the right, with a legend for domain descriptions.
````
Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.homologsf.al.mid.treefile ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.rps-blast.out ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.tree.rps.pdf
````
Formats and prints the RPS-BLAST output in a human-readable format, with tab-delimited columns and a compact view for easy examination.
````
mlr --inidx --ifs "\t" --opprint cat ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.rps-blast.out | tail -n +2 | less -S
````
Identifies and counts the number of domain annotations for each protein sequence in the dataset.
````
cut -f 1 ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.rps-blast.out | sort | uniq -c
````
Extracts the domain titles and counts their occurrences to determine the most frequently detected domains in the dataset.
````
cut -f 6 ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.rps-blast.out | sort | uniq -c
````
Computes the length of each annotated domain and sorts them in descending order to identify the longest domains in the dataset.
````
awk '{a=$4-$3;print $1,"\t",a;}' ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.rps-blast.out | sort -k2nr
````
Retrieves the query sequence IDs along with their associated E-values to assess the statistical significance of the domain predictions.
````
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/pyruvate-dehydrogenase-complex/pyruvate-dehydrogenase-complex.rps-blast.out
````
