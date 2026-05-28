# allorhizobium-salinity-adaptation


Code for "Stress adaptation of free-living microbes generates novel benefits to plant hosts". Plant and microbial phenotypic data are included here, as well as ancestor assembly and processed SNP calls. Raw sequencing data can be found at  https://www.ncbi.nlm.nih.gov/bioproject/PRJNA135652. Direct questions to Kevin Ricks, at kevin.ricks@utoronto.ca

The **code** directory houses the scripts to process sequence data and recreate analyses. These analyses have been subset into folders for analyzing the phenotypic and genotypic data. 

The folder for the phenotypic analysis, "phenotypic_analysis" includes building  models around plant growth and microbial growth dependent on microbial strain history. Data to recreate these analyses can be found in the **data** folder, subfolder plant_microbe_phenotypes

The genomics folder, "genomics", includes analyses for: building and annotating an assembly of the ancestral Allorhizobium strain, calling SNPs from the evolved strains in reference to this ancestor, and characterizing the distribution of these SNPs by the historic selective treatments. These analyses are assumed to be conducted in the above described order. 

Recreating these analyses from the beginning requires the raw sequencing files, found in the associated NCBI bioproject. However, we've included the output from the assembly, annotation, and SNP calling in the **data** folder, from which subsequent analyses can be conducted.
