# Protein collision detection


Identification of steric crashes in protein interfaces

To identify atomic clashes in protein interfaces, I have developed a method based on GEPOL.
GEPOL approximates the molecular surface of proteins by triangular tessellation of the atoms, i.e. the
surface of atoms is represented by a set of triangles. On a first step, we identify the interface of
the complex by filtering those surface regions that do not change when comparing bound and
unbound conformations. Subsequently, the triangles enclosed in the interface(s) of the complex
are extracted, and these are used to identify steric crashes.For each pair of triangles from 
different protein chains and a distance smaller than 0.3 Angstroms we compute: 
(i) the normal vector; 
(ii) a position vector calculated as the difference
between the centers of the two triangles; and 
(iii) the dot product between the normal vector 
the position unit vector. By comparing the sign of the dot product, it is possible to identify
whether the two given atoms are indeed overlapping, i.e. steric crash (see Supplementary
Figure S1). This information is provided in the summary table under the column ‘crashes’ with a
yes/no value. If the value is yes, then a color gradient is shown together with a clickable pop-up
window identifying the residues and atoms as well as the penetration area accounting for the
steric hindrances.
