# PhosphoKin
Authors: Panagiota-Aggeliki Galliou (ag.gal.work@gmail.com), Kleio-Maria Verrou (kleioverrou@yahoo.com)  

Tool for searching kinase binding sites and potential phosphorylated residues in a protein. This tool also categorizes the phosphorylated residues according to their location relatively to the proteins' known active sites (inside of, close to within six residues proximity and outside of active sites). It also links kinases to active sites according to the phosphorylated residue they are responsible for and shows their activity inside of, close to and outside of active sites. It finds the possibly responsible kinases for the experimentally observed phosphorylated residues in a protein according to the predicted phosphorylated residues for the protein.

# Reference

Under reviewing.

# Requirements

PhosphoKin is a tool written in Python 3 and, therefore, it can be run in most operating systems.

# How to run

## Needed files

PhosphoKin uses as input four .txt files; one containing the protein sequence, one containing the protein’s active sites, one containing the experimentally assigned phosphorylated residues in the protein and one containing the kinases along with their recognition motifs. 

### 1.) The Active Sites file:
The file must have the following format: Start_of_active_site_1 - End_of_active_site_1, Start_of_active_site_2 - End_of_active_site_2, etc.

e.g.

>13-22, 150-162, 1147-1458
### 2.) The Experimentally Observed Phosphorylations file:
The file must have the following format: phosphorylated_residue_1, phosphorylated_residue_2, phosphorylated_residue_3, etc.

e.g.

>S4, S156, T445
### 3.) The Motifs file:
The file must have the following format: Name_of_kinase_1[space]protein type: motif_1-motif_2-motif_3 Name_of_kinase_2[space]protein type: motif_1-motif_2-etc.

e.g.

>H1K Kinase:[pS/pT]P[R/K]-[pS/pT]PX[R/K]-[R/K][pS/pT]P
>
>ATM Kinase:pSQ-[P/L/I/M]X[L/I/E/D]pSQ-LpSQE
   
### 4.) The Protein Sequence file:
The file must have a faste format (https://en.wikipedia.org/wiki/FASTA_format)

e.g.

>`>`sp|P11047|LAMC1_HUMAN Laminin subunit gamma-1 OS=Homo sapiens OX=9606 GN=LAMC1 PE=1 SV=3 MRGSHRAAPALRPRGRLWPVLAVLAAAAAAGCAQAAMDECTDEGGRPQRCMPEFVNAAFN VTVVATNTCGTPPEEYCVQTGVTGVTKSCHLCDAGQPHLQHGAAFLTDYNNQADTTWWQS QTMLAGVQYPSSINLTLHLGKAFDITYVRLKFHTSRPESFAIYKRTREDGPWIPYQYYSG SCENTYSKANRGFIRTGGDEQQALCTDEFSDISPLTGGNVAFSTLEGRPSAYNFDNSPVL
# In Terminal

`python3 PhosphoKin.py`

PhosphoKin is quite interactive and will ask from the user during the run each file. Remember to provide the exact filename (with its extension) that corresponds to each file. Furthermore, if the files are not in the same directory were the PhosphoKin is saved, insert the whole path of the files name. 

# The output

As output the tool produces four .txt files; one for the prediction of phosphorylation sites and phosphorylated residues in the protein, one for the identification of possibly responsible kinases for the experimentally observed phosphorylated residues in the protein, one for the categorization of phosphorylated residues according to their location relatively to active sites and one for the association of kinases with active sites. For the analysis of the output, please read the reference.

# Got a Question?

Please do not hesitate to ask any question via email: Panagiota-Angeliki Galliou (ag.gal.work@gmail.com).
