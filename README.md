# PRotein COllision Detection Algorithm (PROCOD)


Identification of steric crashes in protein interfaces

To identify atomic clashes in protein interfaces, I have developed a method called PROCOD based on GEPOL.
GEPOL approximates the molecular surface of proteins by triangular tessellation of the atoms,
i.e. the surface of atoms is represented by a set of triangles. On a first step, it identify
the interface of the complex by filtering those surface regions that do not change when comparing
bound and unbound conformations. 

Subsequently, the triangles enclosed in the interface(s) of the complex are extracted, and these
are used to identify steric crashes.For each pair of triangles from different protein chains and
a distance smaller than 0.3 Angstroms we compute: 

(i) the normal vector; 
(ii) a position vector calculated as the difference between the centers of the two triangles; and 
(iii) the dot product between the normal vector,the position unit vector.

By comparing the sign of the dot product, it is possible to identify
whether the two given atoms are indeed overlapping, i.e. steric crash). 
This method is been implemented in SPServer (http://sbi.upf.edu/spserver/), a web server to assess 
the quality of protein folds and protein-protein interactions using knowledge-based potentials.


References : 
Aguirre-Plans, J., Meseguer, A., Molina-Fernandez, R. et al. SPServer: split-statistical potentials for the analysis of protein structures and protein–protein interactions. BMC Bioinformatics 22, 4 (2021). https://doi.org/10.1186/s12859-020-03770-5
