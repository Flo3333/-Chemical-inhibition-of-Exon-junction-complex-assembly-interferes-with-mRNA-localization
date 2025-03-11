Centrosome Project
==================
*Updated on 2023/10/15*

This project aims at quantifing the correlations between centrosomal and specific RNA positions in different treatment conditions.
We are working with a set of images with the given information : 


Files format
------------

filenames example : 

1. **r04c03-Y14_FLAP XYf01-Cy3-sk1fk1fl1.tiff**  
 *r04c03* : position on plate in microscope  
 *Y14* : cell line  
 *FLAP XY* : RNA targeted by fish microscopy  
 *f01* : field of view number  
 *Cy3* : channel  
 *sk1fk1fl1* : microscope parameter - unimportant  
 *tiff* : file format extension  

2. **r04c04-Y14_DMSO_DYNC1H1f01-Cy3-sk1fk1fl1.tiff**  
 *r04c04* : position on plate in microscope  
 *Y14* : cell line  
 *DMSO* : treatment
 *DYNC1H1* : RNA targeted by fish microscopy  
 *f01* : field of view number  
 *Cy3* : channel  
 *sk1fk1fl1* : microscope parameter - unimportant  
 *tiff* : file format extension  

Data sets
---------

**EJC_HeLacentrin**

3 channels : 
1. **DAPI** --> Nucleus
2. **Cy3** --> RNA single molecules + auto_fluorescence for cytoplasm
3. **EGFP** --> centrosome

4 treatments :
1. DMSO
2. CB16h
3. inh216h
4. untreated

--> Quantification of RNA localisation relative to centrosome with different treatment
--> Quantification of RNA clusters formation (gene of interest **DYNC1H1**)

**EJC_HEKY14_withoutIF**

2 channels : 
1. **DAPI** --> Nucleus
2. **Cy3** --> RNA single molecules + auto_fluorescence for cytoplasm

7 treatments :
1. DMSO
2. Inh2 pendant 16h
3. dTAG pendant 24h
4. dTAG pendant 16h
5. dTAG pendant 12h
6. dTAG pendant 8h
7. untreated