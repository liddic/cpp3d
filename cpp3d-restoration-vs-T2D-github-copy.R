#########################
#
# cpp3d: Compound processing potential, at the scale of individual compounds, with 3-d chemical/bioenergetic mapping visualisation (O:C, H:C, N:C)
# Craig Liddicoat | Flinders University, South Australia 
#
#########################

# record library and version info
.libPaths() # "/Library/Frameworks/R.framework/Versions/4.2/Resources/library"

R.Version()
# "R version 4.2.2 (2022-10-31)"
citation()
# R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL
# https://www.R-project.org/.

library(readxl); packageVersion("readxl") # '1.4.1'
library(plyr); packageVersion("plyr") # '1.8.8'
library(dplyr); packageVersion("dplyr") # '1.0.10'
library(vegan);packageVersion("vegan") # '2.6.4'
library(phyloseq); packageVersion("phyloseq") # '1.42.0'
library(ggplot2); packageVersion("ggplot2") # '3.4.0'
library(grid); packageVersion("grid") #  '4.2.2'
library(reshape2); packageVersion("reshape2") # '1.4.4'
library(tidyr); packageVersion("tidyr") # '1.2.1'
library(corrr); packageVersion("corrr") # '0.4.4'
library(ggforce); packageVersion("ggforce") # '0.4.1'
library(ggrepel); packageVersion("ggrepel") # '0.9.2'
library(stringdist); packageVersion("stringdist") # ‘0.9.10’
library(stringr); packageVersion("stringr") # ‘1.5.0’
library(doParallel); packageVersion("doParallel") # '1.0.17'
library(RColorBrewer); packageVersion("RColorBrewer") # '1.1.3'
library(ggpp); packageVersion("ggpp") # ‘0.5.0’ # https://cran.r-project.org/web/packages/ggpp/vignettes/grammar-extensions.html
library(corrplot)                  ;packageVersion("corrplot") #  '0.92'
library(caret)                     ;packageVersion("caret") # '6.0.93'
library(MASS)                     ;packageVersion("MASS") # ‘7.3.58.1’
library(ggsignif); packageVersion("ggsignif") # '0.6.4'
library(moments)                  ;packageVersion("moments") # ‘0.14.1’
library(ANCOMBC); packageVersion("ANCOMBC") # ‘2.0.1’
library(grDevices); packageVersion("grDevices") #  '4.2.2'
library(ggbiplot); packageVersion("ggbiplot") #  ‘0.55’
library(viridis); packageVersion("viridis") #  ‘0.6.2’
library(FSA); packageVersion("FSA") # '0.9.3'
library(rcompanion); packageVersion("rcompanion") # '2.4.18'
library(fields); packageVersion("fields") # ‘14.1’
library(car); packageVersion("car") # ‘3.1.1’
library(multcompView); packageVersion("multcompView") # ‘0.1.8’
library(gtools); packageVersion("gtools") # ‘3.9.4’
library(igraph); packageVersion("igraph") #  '1.4.2'
library(pheatmap); packageVersion("pheatmap") # '1.0.12'
library(colorspace); packageVersion("colorspace") # ‘2.1.0’
library(qvalue); packageVersion("qvalue") #  ‘2.30.0’
library(glmnet); packageVersion("glmnet") # ‘4.1.7’
library(randomForest); packageVersion("randomForest") # ‘4.7.1.1’
library(epiR); packageVersion("epiR") # ‘2.0.53’

#library(zCompositions); packageVersion("zCompositions") # '1.4.0.1'
#library(propr); packageVersion("propr") # '4.2.6'
#library(hexbin); packageVersion("hexbin") # '1.28.2'
#library(raster); packageVersion("raster") # '3.6.11' also loads package 'sp'
#detach("package:raster", unload=TRUE)
#unloadNamespace('raster')


#########################
## CAUTION!! save.image("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/cpp3d-indiv-resto-vs-t2d-WORKSPACE-v3.RData")
##      load("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/cpp3d-indiv-resto-vs-t2d-WORKSPACE-v3.RData")
#########################

workdir <- "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"
setwd(workdir)
getwd()


par.default <- par()



#### ModelSEED lookup tables
#-------------------------

seed_db_dir <- "/Users/lidd0026/WORKSPACE/DATA/ModelSEEDDB/select_files"


## Unique_ModelSEED_Reaction_ECs

rxn_ECs.lut <- read_excel(path = paste0(seed_db_dir,"/","Unique_ModelSEED_Reaction_ECs.xlsx") , range = "A1:C30834")
rxn_ECs.lut <- as.data.frame(rxn_ECs.lut)
head(rxn_ECs.lut)
names(rxn_ECs.lut) # "ModelSEED ID" "External ID"  "Source"
names(rxn_ECs.lut) <- c("rxn_id", "EC_id",  "Source")

dim(rxn_ECs.lut) # 30833     3
length(unique(rxn_ECs.lut$rxn_id)) # 25771
length(unique(rxn_ECs.lut$EC_id)) # 7354


## Unique_ModelSEED_Reaction_Pathways

# rxn_pathways.lut <- read_excel(path = paste0(seed_db_dir,"/","Unique_ModelSEED_Reaction_Pathways.xlsx") , range = "A1:C121445")
# rxn_pathways.lut <- as.data.frame(rxn_pathways.lut)
rxn_pathways.lut <- read.csv(file = paste0(seed_db_dir,"/","Unique_ModelSEED_Reaction_Pathways.txt"), header = TRUE, sep = "\t")
class(rxn_pathways.lut) # "data.frame"

head(rxn_pathways.lut)
names(rxn_pathways.lut) # "ModelSEED.ID" "External.ID"  "Source"
names(rxn_pathways.lut) <- c("rxn_id", "External_rxn_name",  "Source")

dim(rxn_pathways.lut) # 121444      3
length(unique(rxn_pathways.lut$rxn_id)) # 21734
length(unique(rxn_pathways.lut$External_rxn_name)) # 3819


## Model SEED Subsystems

subsys.lut <- read_excel(path = paste0(seed_db_dir,"/","ModelSEED_Subsystems.xlsx") , range = "A1:E9821")
subsys.lut <- as.data.frame(subsys.lut)
dim(subsys.lut) # 9820    5

head(subsys.lut)
tail(subsys.lut)
names(subsys.lut) # "Class"     "Sub-class" "Name"      "Role"      "Reaction" 
names(subsys.lut) <- c("Class",    "Subclass", "Name",  "Role",     "Reaction")

subsys.lut[sample(1:dim(subsys.lut)[1], 10,replace = FALSE), ]
#                                     Class                                    Subclass                                                Name                                                                           Role Reaction
# 2984          Clustering-based subsystems                                           -                              CBSS-196620.1.peg.2477                                      Maltose O-acetyltransferase (EC 2.3.1.79) rxn01133
# 6834 Fatty Acids, Lipids, and Isoprenoids                                 Fatty acids                       Fatty_Acid_Biosynthesis_FASII          (3R)-hydroxymyristoyl-[acyl carrier protein] dehydratase (EC 4.2.1.-) rxn05427
# 7057 Fatty Acids, Lipids, and Isoprenoids                                 Fatty acids                   Unsaturated_Fatty_Acid_Metabolism                  3-oxoacyl-[acyl-carrier-protein] synthase, KASI (EC 2.3.1.41) rxn05350
# 9384                      Stress Response                              Osmotic stress Choline_and_Betaine_Uptake_and_Betaine_Biosynthesis                   Glycine betaine ABC transport system, permease protein OpuAB rxn05181
# 5904              Experimental Subsystems                                           -                              Transporters_In_Models                 Oligopeptide transport ATP-binding protein oppF (TC 3.A.1.5.1) rxn05539
# 5267              Experimental Subsystems                                           -              Sugar_catabolome_in_Shewanella_species                                              Alpha-galactosidase (EC 3.2.1.22) rxn00818
# 2626                Cell Wall and Capsule   Capsular and extracellular polysacchrides                                 Alginate_metabolism                                      Acetoin (diacetyl) reductase (EC 1.1.1.5) rxn01685
# 8538                   Protein Metabolism                        Protein biosynthesis                       Translation_factors_bacterial                                       Methionine aminopeptidase (EC 3.4.11.18) rxn12635
# 1853                        Carbohydrates                                Fermentation                 Acetyl-CoA_fermentation_to_Butyrate                     Acyl-CoA dehydrogenase, short-chain specific (EC 1.3.99.2) rxn10012
# 776           Amino Acids and Derivatives Lysine, threonine, methionine, and cysteine                     Lysine_Biosynthesis_DAP_Pathway 2,3,4,5-tetrahydropyridine-2,6-dicarboxylate N-acetyltransferase (EC 2.3.1.89) rxn03030
## Correspondences:

# Superfocus (imported 'tab' object): Subsystem Level 1 = subsys.lut$Class
#                                     Subsystem Level 2 = subsys.lut$Subclass
#    (skip the subsys.lut$Name field)
#                                     Function = subsys.lut$Role  # replace " " with "_" to match Superfocus output

dim(subsys.lut) # 9820    5
length(unique(subsys.lut$Class )) # 29
length(unique(subsys.lut$Subclass )) # 143
length(unique(subsys.lut$Name )) # 666
length(unique(subsys.lut$Role )) # 1885
length(unique(subsys.lut$Reaction )) # 1986


## Reactions

rxns.lut <- read_excel(path = paste0(seed_db_dir,"/","reactions.xlsx") , range = "A1:V43775")
rxns.lut <- as.data.frame(rxns.lut)
dim(rxns.lut) # 43774    22
names(rxns.lut)
# [1] "id"                "abbreviation"      "name"              "code"              "stoichiometry"     "is_transport"      "equation"          "definition"       
# [9] "reversibility"     "direction"         "abstract_reaction" "pathways"          "aliases"           "ec_numbers"        "deltag"            "deltagerr"        
# [17] "compound_ids"      "status"            "is_obsolete"       "linked_reaction"   "notes"             "source" 

head(rxns.lut)

length(unique(rxns.lut$id)) # 43774
length(unique(rxns.lut$name)) # 27668
length(unique(rxns.lut$code)) # 36212
length(unique(rxns.lut$equation)) # 37310
length(unique(rxns.lut$aliases)) # 32718    # SEARCH aliases ???????
length(unique(rxns.lut$ec_numbers)) # 7607


## Compounds

compounds.lut <- read_excel(path = paste0(seed_db_dir,"/","compounds.xlsx") , range = "A1:T33993")
warnings()
compounds.lut <- as.data.frame(compounds.lut)
dim(compounds.lut) # 33992    20
names(compounds.lut)
# [1] "id"                "abbreviation"      "name"              "formula"           "mass"              "source"            "inchikey"          "charge"           
# [9] "is_core"           "is_obsolete"       "linked_compound"   "is_cofactor"       "deltag"            "deltagerr"         "pka"               "pkb"              
# [17] "abstract_compound" "comprised_of"      "aliases"           "smiles"

head(compounds.lut)

length(unique(compounds.lut$id)) # 33992
length(unique(compounds.lut$formula)) # 16763


#-------------------------


#### Atomic ratio lookup function? O, C, H, N & map all ModelSEED compounds
#    Include multi-elements & element counts per 'beyond vK' paper
#    Include all critical nutrients for microbes
#-------------------------

# from earlier
## Compounds

compounds.lut <- read_excel(path = paste0(seed_db_dir,"/","compounds.xlsx") , range = "A1:T33993") # excludes column 'U'
# warnings()
compounds.lut <- as.data.frame(compounds.lut)
dim(compounds.lut) # 33992    20
names(compounds.lut)
# [1] "id"                "abbreviation"      "name"              "formula"           "mass"              "source"            "inchikey"          "charge"           
# [9] "is_core"           "is_obsolete"       "linked_compound"   "is_cofactor"       "deltag"            "deltagerr"         "pka"               "pkb"              
# [17] "abstract_compound" "comprised_of"      "aliases"           "smiles"

head(compounds.lut)
head(compounds.lut[ ,1:7])
#         id abbreviation  name       formula mass           source                    inchikey
# 1 cpd00001          h2o   H2O           H2O   18 Primary Database XLYOFNOQVPJJNP-UHFFFAOYSA-N
# 2 cpd00002          atp   ATP C10H13N5O13P3  504 Primary Database ZKHQWZAMYRWXGA-KQYNXXCUSA-K
# 3 cpd00003          nad   NAD C21H26N7O14P2  662 Primary Database BAWFJGJZGIEFAR-NNYOXOHSSA-M
# 4 cpd00004         nadh  NADH C21H27N7O14P2  663 Primary Database BOPGDPNILDQYTO-NNYOXOHSSA-L
# 5 cpd00005        nadph NADPH C21H26N7O17P3  742 Primary Database ACFIXJIJDZMPPO-NNYOXOHSSA-J
# 6 cpd00006         nadp  NADP C21H25N7O17P3  741 Primary Database XJLXINKUBYWONI-NNYOXOHSSA-K

length(unique(compounds.lut$id)) # 33992
length(unique(compounds.lut$formula)) # 16763
length(unique(compounds.lut$name)) # 33844


## get all combinations of 'O?' 'C?' 'H?' - to isolate oxygen, carbon, hydrogen containing compounds

# https://pmc.ncbi.nlm.nih.gov/articles/PMC4100946/
#   Elemental Economy: microbial strategies for optimizing growth in the face of nutrient limitation
# 
# Merchant, S.S. & Helmann, J.D. Chapter 2 - Elemental Economy: Microbial Strategies for Optimizing Growth in the Face of Nutrient Limitation. In: Poole, R.K. (ed). Advances in Microbial Physiology, vol. 60. Academic Press, 2012, pp 91-210.
# 
# Table 1 - include 'Required for all cells' and 'Required for most cells'
# 
# Required for All Cells
# C	basis of all organic molecules
# H	H2O, organic molecules
# N	organic molecules, esp. proteins and nucleic acids
# O	H2O, organic molecules
# # Already done
# 
# P	nucleic acids, NTPs, metabolites, phospholipids
# S	proteins, glutathione and LMW thiols, biotin, lipoic acid, thiamin
# Mg	major cation; cofactor for phosphotransferase reactions
# Zn	enzyme cofactor, protein folding
# 
# Required for Most Cells
# K	major cation, common in cells
# Ca	major cation, required by many eukaryotes
# Mn	enzyme cofactor, ribonucleotide reductase, SOD, PS II
# Fe	heme, iron-sulfur cluster, non-heme enzymes
# Co	enzyme cofactor, B12-dependent enzymes
# Cu	enzyme cofactor, electron carrier, respiration, SOD
# Mo	FeMoCo cofactor (nitrogenase), Mo cofactor enzymes



temp <- compounds.lut

sel <- grep(pattern = "Na", x = temp$formula)

sel <- grep(pattern = "P", x = temp$formula)
sel <- grep(pattern = "Si", x = temp$formula)
sel <- grep(pattern = "Mg", x = temp$formula)
sel <- grep(pattern = "Zn", x = temp$formula)
sel <- grep(pattern = "K", x = temp$formula)
sel <- grep(pattern = "Ca", x = temp$formula)
sel <- grep(pattern = "Mn", x = temp$formula)
sel <- grep(pattern = "Fe", x = temp$formula)
sel <- grep(pattern = "Co", x = temp$formula)
sel <- grep(pattern = "Cu", x = temp$formula)
sel <- grep(pattern = "Mo", x = temp$formula)


temp$formula[sel[1:100]]


#.          1.      2.      3.       4.       5.       6.        7.        8.      9.     10.          11             12.        13.   14.       15.               16
form <- c("NaCl", "H2O", "CaOH", "CaOCl2Rh","HgCl","C5H8O7PR", "MoOH2", "CoCH3", "OsR", "ScO2", "C10H12N4O12P2", "C15H27N5O5", "HgR", "Zn",
          
          "C25H26N9NaO8S2", "C12H7Cl2NNaO2", "C21H26N7O14P2", "H9AlFeMgO15Si4", "C55H75N4O6Zn", "C16H16K3O6PS", "C23H30CaN3Na3O11",
          
          "C4H6MnN2R2S4", "C33H30FeN4O4", "C55H83CoN14O9R2", "C7H4ClCuNO3S", "C20H22MoN10O15P2S2", "")

atom <- "O" # avoid: Os, Og
atom <- "C" # avoid: Cs, Ca, Ce, Cr, Co, Cm, Cu, Cd, Cn, Cf, Cl, 
atom <- "H" # avoid: Hf, Hs, Hg, Ho, He
atom <- "N" # avoid: Na, Nd, Ne, Np, Ni, Nh, Nb, No

atom <- "P" # avoid: Pd, Pt, Pb, Po
atom <- "S" # avoid: Si, Sc, Se, Sr, Sn, Sb
atom <- "Mg" #
atom <- "Zn" # 

atom <- "K" # avoid Kr 
atom <- "Ca" # 
atom <- "Mn" # 
atom <- "Fe" # 
atom <- "Co" # 
atom <- "Cu" # 
atom <- "Mo" # 





atomic_no <- function(form, atom) {
  ##form=form
  ##atom="O"
  ##atom="N"
  coefs <- list()
  for (c in 1:length(form)) {
    #c<-1
    pos <- str_locate(string = form[c], pattern = atom)[1]
    if (is.na(pos)) { # not present
      ##print("zero")
      coef<-0
    } else if (pos==nchar(form[c])) { # at the end, i.e. coef = 1.    # A 
      ##print("A")
      coef<-1
    } else {
      # check for off-targets ... otherwise decide coef = 1, 2, etc?       # B
      ##print("B")
      checkstring <- substring(text = form[c], first = pos, last = pos+1)
      
      if (atom == "O" & checkstring %in% c("Os", "Og")) {   # C
        ##print("C")
        coef<-0
      } else if (atom == "C" & checkstring %in% c("Cs", "Ca", "Ce", "Cr", "Co", "Cm", "Cu", "Cd", "Cn", "Cf", "Cl")) { # D
        ##print("D")
        coef<-0
      } else if (atom == "H" & checkstring %in% c("Hf", "Hs", "Hg", "Ho", "He")) { # E 
        ##print("E")
        coef<-0
        
      } else if (atom == "N" & checkstring %in% c("Na", "Nd", "Ne", "Np", "Ni", "Nh", "Nb", "No")) { # F 
        ##print("F")
        coef<-0
        
      } else if (atom == "P" & checkstring %in% c("Pd", "Pt", "Pb", "Po")) { # G 
        ##print("G")
        coef<-0
        
      } else if (atom == "S" & checkstring %in% c("Si", "Sc", "Se", "Sr", "Sn", "Sb")) { # H 
        ##print("H")
        coef<-0
        
      } else if (atom == "Mg" & checkstring %in% c("")) { # I 
        ##print("I")
        coef<-0
        
      } else if (atom == "Zn" & checkstring %in% c("")) { # J 
        ##print("J")
        coef<-0
        
      } else if (atom == "K" & checkstring %in% c("Kr")) { # K 
        ##print("K")
        coef<-0
        
      } else if (atom == "Ca" & checkstring %in% c("")) { # L 
        ##print("L")
        coef<-0
        
      }else if (atom == "Mn" & checkstring %in% c("")) { # M 
        ##print("M")
        coef<-0
        
      }else if (atom == "Fe" & checkstring %in% c("")) { # N 
        ##print("N")
        coef<-0
        
      }else if (atom == "Co" & checkstring %in% c("")) { # O 
        ##print("O")
        coef<-0
        
      }else if (atom == "Cu" & checkstring %in% c("")) { # P 
        ##print("P")
        coef<-0
        
      }else if (atom == "Mo" & checkstring %in% c("")) { # Q 
        ##print("Q")
        coef<-0
        
      } else { # R
        ##print("R")
        coef<-1
        coef.try <- coef
        pos.end <- pos+1
        
        while (is.numeric(coef.try) & !is.na(coef.try) & pos.end <= nchar(form[c])+1 ) {
          # substring allows last index to be too long & need to allow 'pos.end <= nchar(form[c])+1' for atomic numbers at end of formula
          coef <- coef.try
          coef.try <- as.numeric( substring(text = form[c], first = pos+1, last = pos.end) )
          pos.end <- pos.end+1
          ##print("G")
        } # END while loop
      } # END else
      
    } # END else
    
    coefs[[c]] <- coef
    
  } # END c loop
  return( unlist(coefs) )
}

form
# [1] "NaCl"               "H2O"                "CaOH"               "CaOCl2Rh"           "HgCl"               "C5H8O7PR"           "MoOH2"             
# [8] "CoCH3"              "OsR"                "ScO2"               "C10H12N4O12P2"      "C15H27N5O5"         "HgR"                "Zn"                
# [15] "C25H26N9NaO8S2"     "C12H7Cl2NNaO2"      "C21H26N7O14P2"      "H9AlFeMgO15Si4"     "C55H75N4O6Zn"       "C16H16K3O6PS"       "C23H30CaN3Na3O11"  
# [22] "C4H6MnN2R2S4"       "C33H30FeN4O4"       "C55H83CoN14O9R2"    "C7H4ClCuNO3S"       "C20H22MoN10O15P2S2" "" 

atomic_no(form, "O")
# [1]  0  1  1  1  0  7  1  0  0  2 12  5  0  0  8  2 14 15  6  6 11  0  4  9  3 15  0

atomic_no(form="NaCl",atom="C")
# [1] 0
atomic_no(form=form,atom="C")
# [1]  0  0  0  0  0  5  0  0  0  0 10 15  0  0 25 12 21  0 55 16 23  4 33 55  7 20  0
atomic_no(form=form,atom="H")
# 0  2  1  0  0  8  2  3  0  0 12 27  0  0 26  7 26  9 75 16 30  6 30 83  4 22  0
atomic_no(form=form,atom="N")
# 0 0 0 0 0 0 0 0 0 0 4 5 0 0 9 1
atomic_no(form=form,atom="P")
# [1] 0 0 0 0 0 1 0 0 0 0 2 0 0 0 0 0 2 0 0 1 0 0 0 0 0 2 0
atomic_no(form=form,atom="S")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 4 0 0 1 2 0
atomic_no(form=form,atom="Mg")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
atomic_no(form=form,atom="Zn")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0
atomic_no(form=form,atom="K")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0
atomic_no(form=form,atom="Ca")
# [1] 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
atomic_no(form=form,atom="Mn")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
atomic_no(form=form,atom="Fe")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0
atomic_no(form=form,atom="Co")
# [1] 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
atomic_no(form=form,atom="Cu")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
atomic_no(form=form,atom="Mo")
# [1] 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0


df.comp <- data.frame(id=compounds.lut$id, 
                      abbrev=compounds.lut$abbreviation,
                      name=compounds.lut$name,
                      form=compounds.lut$formula,
                      OC_ratio=NA,
                      HC_ratio=NA,
                      NC_ratio=NA,
                      
                      PC_ratio=NA,
                      
                      # new
                      NP_ratio = NA,
                      O_count = NA,
                      N_count = NA,
                      P_count = NA,
                      S_count = NA,
                      mass = NA,
                      
                      SC_ratio=NA,
                      MgC_ratio=NA,
                      ZnC_ratio=NA,
                      
                      KC_ratio=NA,
                      CaC_ratio=NA,
                      MnC_ratio=NA,
                      FeC_ratio=NA,
                      CoC_ratio=NA,
                      CuC_ratio=NA,
                      MoC_ratio=NA
                      
)



for (i in 1:length(df.comp$id)) {
  #i<-1
  formx.char <- df.comp$form[i]
  
  # O:C ratio (replace Inf with NA when C not present)
  df.comp$OC_ratio[i] <- atomic_no(form = formx.char, atom = "O")/atomic_no(form = formx.char, atom = "C")
  df.comp$OC_ratio[i][is.infinite(df.comp$OC_ratio[i])] <- NA
  
  # H:C ratio (replace Inf with NA when C not present)
  df.comp$HC_ratio[i] <- atomic_no(form = formx.char, atom = "H")/atomic_no(form = formx.char, atom = "C")
  df.comp$HC_ratio[i][is.infinite(df.comp$HC_ratio[i])] <- NA
  
  # N:C ratio (replace Inf with NA when C not present)
  df.comp$NC_ratio[i] <- atomic_no(form = formx.char, atom = "N")/atomic_no(form = formx.char, atom = "C")
  df.comp$NC_ratio[i][is.infinite(df.comp$NC_ratio[i])] <- NA
  
  # PC_ratio
  df.comp$PC_ratio[i] <- atomic_no(form = formx.char, atom = "P")/atomic_no(form = formx.char, atom = "C")
  df.comp$PC_ratio[i][is.infinite(df.comp$PC_ratio[i])] <- NA
  
  # new
  # NP_ratio
  df.comp$NP_ratio[i] <- atomic_no(form = formx.char, atom = "N")/atomic_no(form = formx.char, atom = "P")
  df.comp$NP_ratio[i][is.infinite(df.comp$NP_ratio[i])] <- NA
  
  # O_count
  df.comp$O_count[i] <- atomic_no(form = formx.char, atom = "O")
  
  # N_count
  df.comp$N_count[i] <- atomic_no(form = formx.char, atom = "N")
  
  # P_count
  df.comp$P_count[i] <- atomic_no(form = formx.char, atom = "P")
  
  # S_count
  df.comp$S_count[i] <- atomic_no(form = formx.char, atom = "S")
  
  # mass
  df.comp$mass[i] <- compounds.lut$mass[i]
  
  # SC_ratio
  df.comp$SC_ratio[i] <- atomic_no(form = formx.char, atom = "S")/atomic_no(form = formx.char, atom = "C")
  df.comp$SC_ratio[i][is.infinite(df.comp$SC_ratio[i])] <- NA
  
  # MgC_ratio
  df.comp$MgC_ratio[i] <- atomic_no(form = formx.char, atom = "Mg")/atomic_no(form = formx.char, atom = "C")
  df.comp$MgC_ratio[i][is.infinite(df.comp$MgC_ratio[i])] <- NA
  
  # ZnC_ratio
  df.comp$ZnC_ratio[i] <- atomic_no(form = formx.char, atom = "Zn")/atomic_no(form = formx.char, atom = "C")
  df.comp$ZnC_ratio[i][is.infinite(df.comp$ZnC_ratio[i])] <- NA
  
  # KC_ratio
  df.comp$KC_ratio[i] <- atomic_no(form = formx.char, atom = "K")/atomic_no(form = formx.char, atom = "C")
  df.comp$KC_ratio[i][is.infinite(df.comp$KC_ratio[i])] <- NA
  
  # CaC_ratio
  df.comp$CaC_ratio[i] <- atomic_no(form = formx.char, atom = "Ca")/atomic_no(form = formx.char, atom = "C")
  df.comp$CaC_ratio[i][is.infinite(df.comp$CaC_ratio[i])] <- NA
  
  # MnC_ratio
  df.comp$MnC_ratio[i] <- atomic_no(form = formx.char, atom = "Mn")/atomic_no(form = formx.char, atom = "C")
  df.comp$MnC_ratio[i][is.infinite(df.comp$MnC_ratio[i])] <- NA
  
  # FeC_ratio
  df.comp$FeC_ratio[i] <- atomic_no(form = formx.char, atom = "Fe")/atomic_no(form = formx.char, atom = "C")
  df.comp$FeC_ratio[i][is.infinite(df.comp$FeC_ratio[i])] <- NA
  
  # CoC_ratio
  df.comp$CoC_ratio[i] <- atomic_no(form = formx.char, atom = "Co")/atomic_no(form = formx.char, atom = "C")
  df.comp$CoC_ratio[i][is.infinite(df.comp$CoC_ratio[i])] <- NA
  
  # CuC_ratio
  df.comp$CuC_ratio[i] <- atomic_no(form = formx.char, atom = "Cu")/atomic_no(form = formx.char, atom = "C")
  df.comp$CuC_ratio[i][is.infinite(df.comp$CuC_ratio[i])] <- NA
  
  # MoC_ratio
  df.comp$MoC_ratio[i] <- atomic_no(form = formx.char, atom = "Mo")/atomic_no(form = formx.char, atom = "C")
  df.comp$MoC_ratio[i][is.infinite(df.comp$MoC_ratio[i])] <- NA
  
  print(paste0("completed ",i))
  
}


saveRDS(object = df.comp, file = "df.comp__multi-element-v2.RDS")



## Add additional Compound classes from Rivas-Ubach 2018

df.comp2 <- df.comp
head(df.comp2)
# id abbrev  name          form  OC_ratio HC_ratio  NC_ratio  PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 1 cpd00001    h2o   H2O           H2O        NA       NA       NaN       NaN      NaN       1       0       0       0   18      NaN       NaN       NaN      NaN       NaN
# 2 cpd00002    atp   ATP C10H13N5O13P3 1.3000000 1.300000 0.5000000 0.3000000 1.666667      13       5       3       0  504        0         0         0        0         0
# 3 cpd00003    nad   NAD C21H26N7O14P2 0.6666667 1.238095 0.3333333 0.0952381 3.500000      14       7       2       0  662        0         0         0        0         0
# 4 cpd00004   nadh  NADH C21H27N7O14P2 0.6666667 1.285714 0.3333333 0.0952381 3.500000      14       7       2       0  663        0         0         0        0         0
# 5 cpd00005  nadph NADPH C21H26N7O17P3 0.8095238 1.238095 0.3333333 0.1428571 2.333333      17       7       3       0  742        0         0         0        0         0
# 6 cpd00006   nadp  NADP C21H25N7O17P3 0.8095238 1.190476 0.3333333 0.1428571 2.333333      17       7       3       0  741        0         0         0        0         0
# MnC_ratio FeC_ratio CoC_ratio CuC_ratio MoC_ratio
# 1       NaN       NaN       NaN       NaN       NaN
# 2         0         0         0         0         0
# 3         0         0         0         0         0
# 4         0         0         0         0         0
# 5         0         0         0         0         0
# 6         0         0         0         0         0


names(df.comp2)
# [1] "id"        "abbrev"    "name"      "form"      "OC_ratio"  "HC_ratio"  "NC_ratio"  "PC_ratio"  "NP_ratio"  "O_count"   "N_count"   "P_count"   "S_count"   "mass"     
# [15] "SC_ratio"  "MgC_ratio" "ZnC_ratio" "KC_ratio"  "CaC_ratio" "MnC_ratio" "FeC_ratio" "CoC_ratio" "CuC_ratio" "MoC_ratio"

# all start as unspecified
df.comp2$class <- "Unspecified"

# define classes using "OC_ratio"  "HC_ratio"  "NC_ratio"  "PC_ratio"  "NP_ratio"  "O_count"   "N_count"   "P_count"   "S_count"   "mass"
# based on Rivas-Ubach et al 2018 beyond vK paper

# Lipid: qty 1727
sel <- which(df.comp2$OC_ratio <= 0.6 &
               df.comp2$HC_ratio >= 1.32 & 
               df.comp2$NC_ratio <= 0.126 & 
               df.comp2$PC_ratio < 0.35 & 
               df.comp2$NP_ratio <= 5 )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Lipid"

# Protein (constraints 1): qty 4597
sel <- which(df.comp2$OC_ratio > 0.12 & df.comp2$OC_ratio <= 0.6 &
               df.comp2$HC_ratio > 0.9 & df.comp2$HC_ratio < 2.5 &
               df.comp2$NC_ratio >= 0.126 & df.comp2$NC_ratio <= 0.7 &
               df.comp2$PC_ratio < 0.17 &
               
               df.comp2$N_count >= 1 )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Protein"

# Protein (constraints 2): qty 1199
sel <- which(df.comp2$OC_ratio > 0.6 & df.comp2$OC_ratio <= 1 &
               df.comp2$HC_ratio > 1.2 & df.comp2$HC_ratio < 2.5 &
               df.comp2$NC_ratio > 0.2 & df.comp2$NC_ratio <= 0.7 &
               df.comp2$PC_ratio < 0.17 &
               
               df.comp2$N_count >= 1 )
df.comp2$class[sel] <- "Protein"

# A-Sugar: qty 317
sel <- which(df.comp2$OC_ratio >= 0.61 &
               df.comp2$HC_ratio >= 1.45 & 
               df.comp2$NC_ratio > 0.07 & df.comp2$NC_ratio <= 0.2 & 
               df.comp2$PC_ratio < 0.3 &
               df.comp2$NP_ratio <= 2 &
               df.comp2$O_count >= 3 &
               df.comp2$N_count >= 1
               )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Amino sugar"

# Carbohydrate: qty 1214
sel <- which(df.comp2$OC_ratio >= 0.8 &
               df.comp2$HC_ratio >= 1.65 & df.comp2$HC_ratio < 2.7 & 
               
               df.comp2$N_count == 0
)  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Carbohydrate"

# Nucleotide: qty 169 
sel <- which(df.comp2$OC_ratio >= 0.5 & df.comp2$OC_ratio < 1.7 &
               df.comp2$HC_ratio > 1 & df.comp2$HC_ratio < 1.8 &
               df.comp2$NC_ratio >= 0.2 & df.comp2$NC_ratio <= 0.5 &
               df.comp2$PC_ratio >= 0.1 & df.comp2$PC_ratio < 0.35 &
               df.comp2$NP_ratio > 0.6 & df.comp2$NP_ratio <= 5 &

               df.comp2$N_count >= 2 &
               df.comp2$P_count >= 1 &
               df.comp2$S_count == 0 &
               df.comp2$mass > 305 & df.comp2$mass < 523 
               )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Nucleotide"

# Phytochemical (oxy-aromatic): qty 114
sel <- which(df.comp2$OC_ratio <= 1.15 &
               df.comp2$HC_ratio < 1.32 & 
               df.comp2$NC_ratio < 0.126 & 
               df.comp2$PC_ratio <= 0.2 & 
               df.comp2$NP_ratio <= 3 )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Phytochemical"

#df.comp2[sel, ]

dim(df.comp2) # 33992    25
length(which(df.comp2$class=="Unspecified")) # 24710

table(df.comp2$class)
# Amino sugar  Carbohydrate         Lipid    Nucleotide Phytochemical       Protein   Unspecified 
#         317          1214          1727           169           114          5741         24710 

saveRDS(object = df.comp2, file = "df.comp__multi-element-v2b.RDS")


shapes.class <- c("Amino sugar" = 2 ,
                  "Carbohydrate" = 0 ,
                  "Lipid" = 3,
                  "Nucleotide" = 5,
                  "Phytochemical" = 8  ,
                  "Protein" = 4,
                  "Unspecified" = 1 )

#-------------------------


#### Sun & Badgley - post-mining - build reaction search in parallel - get_reactions & compounds
#### Using SUPER-FOCUS functional potential profiles are from Bioenergetic mapping / CPPv1 paper
#### Individual compound level !!
#-------------------------

phy <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")

## convert each row in functional tax_table to "mean van Krevelen distance to health-associated compounds"

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 30125     4


get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
  
  rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
  rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
  
  ###testing
  # #for (i in 1:dim(df.tax)[1]) {
  # for (i in 1:1000) {
  #i<-1
  #i<-22889
  #i<-22894
  #i<-31768
  #i<-9406
  #i<-9422
  #i<-8
  #1<-13
  #i<-7
  
  sub1 <- df.tax$subsys_L1[i]
  sub2 <- df.tax$subsys_L2[i]
  sub3 <- df.tax$subsys_L3[i]
  
  fxn.temp <- df.tax$fxn[i]
  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  
  # store results corresponding to each Superfocus row
  fxn.list <- list()
  fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
  
  # check for multiple functions/reactions?
  flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
  flag2 <- grepl(pattern = "_@_", x = fxn.temp)
  if (!any(flag1,flag2)==TRUE) {
    # no multiples
    fxns <- fxn.temp
  } else if (flag1==TRUE) {
    fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
  } else {
    fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
  }
  # remove underscores
  ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
  
  # process each fxn & store attributes
  #df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, min_adist_modelSEED=NA, min_amatch_modelSEED=NA, rxns=NA, tot_mean_OC_x=NA, tot_mean_HC_y=NA , tot_mean_NC_z=NA )
  
  df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA) #, compound_name=NA, compound_form=NA, compound_molar_prop=NA, OC_x=NA, HC_y=NA , NC_z=NA )
  
  # # do round brackets interfere with search? - YES
  # lookfor <- "option 4 (this one)"
  # lookuplist <- c("option 1", "option 2", "option 3 (this one)", "option 4 (this one)")
  # grep(pattern = lookfor, x = lookuplist)
  
  # Identify '/' separators with no '_'  ??
  
  for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
    #f<-1
    #f<-2
    f.in <- fxns[f]
    
    # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
    full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
    full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
    
    ## data cleaning
    
    # trim off '_#' and '_##' tags
    trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
    if (!is.na(trim_nchar) & length(trim_nchar)==1) {
      f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
    }
    
    # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
    f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
    
    #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
    #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
    
    full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
    full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
    
    sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
    
    ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
    if (length(sel.rx)>=1) {
      df.fxns$matching_method[f] <- "Exact hierachy match"
      df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
      
    } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
      # search by EC id if present
      
      f.in <- fxns[f] # this goes back to string with brackets for EC
      ## LOOK FOR MULTIPLE ECs ??????????
      # 22889
      # 22894
      # 31768
      
      how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
      
      ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
      #class(ECs)
      ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
      ECs.collapse <- paste0(ECs, collapse = "|")
      
      sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
      
      if (length(how_many_ECs)==0 | length(ECs)==0) {
        # there was a glitch, database typo, or some error in identifying the EC number
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      } else if (length(sel.rx)>=1) {
        # combined EC hits identified
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
        # treat EC hits individually
        sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
        # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
        sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else {
        # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      # END EC matching
      
      
    } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
      ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
      # 1. 'reactions' table by name: rxns.lut$name
      # 2. 'reactions' table by aliases: rxns.lut$aliases
      # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
      # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
      
      if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
        # 1a - exact match - rxns.lut$name
        sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
        #rxns.lut$name[sel.rx]
        df.fxns$matching_method[f] <- "Matched Reactions name"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
        # 2a - exact match - rxns.lut$aliases
        sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
        #rxns.lut$aliases[sel.rx]
        #rxns.lut$name[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Reactions aliases"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
        # 3a - exact match - subsys.lut$Role
        sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
        #subsys.lut$Role[sel.rx]
        #subsys.lut$Reaction[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Subsytem role"
        df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
        # 4a - exact match - rxn_pathways.lut$External_rxn_name
        sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
        
        df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
        df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
        
        
      } else {
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      
      ## DON'T RUN PARTIAL MATCHING AT THIS STAGE
      
      
    } # END function - reaction search
    
    #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
    
    print(paste0("completed fxn ", f))
    
    
    ## now investigate these reactions ...
    # Reactions lookup table: 
    # - "equation": Definition of reaction expressed using compound IDs and after protonation
    # Compounds lookup table:
    # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
    #df.fxns
    
    
    #if (df.fxns$matching_method == "No match found") {
    if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
      
      df.Rxns <- NA
      df.Compounds <- NA
      
    } else { # reaction(s) were identified
      
      # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
      f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
      
      df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
                            rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA ) #, OC_ratios=NA, HC_ratios=NA, NC_ratios=NA, coefwtmean_OC_x=NA, coefwtmean_HC_y=NA, coefwtmean_NC_z=NA)
      
      #df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns, compound_name=NA, compound_form=NA, compound_molar_prop=NA, OC_x=NA, HC_y=NA , NC_z=NA )
      
      for (r in 1:dim(df.Rxns)[1]) {
        #r<-1
        #this_rxn <- "rxn00004"
        this_rxn <- df.Rxns$rxn_id[r]
        sel <- which(rxns.lut$id == this_rxn)
        ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
        ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
        ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
        
        # extract compound info
        
        #df.Rxns$rxn_eqn[r]
        #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
        #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
        
        ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
        # 5 23 43 61
        # 6 25 46 65 83
        
        ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
        # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
        
        if (length(compds)>=1) {
          
          df.Rxns$compds[r] <- paste0(compds, collapse = ";")
          
          ## get compound coefficients?
          start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
          end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
          ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
          
          df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
          
          # get formulas of compounds
          
          formx <-filter(compounds.lut, id %in% compds )
          row.names(formx) <- formx$id
          ( formx.char <- formx[compds, ]$formula )
          # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
          # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
          # "C7H7O7" "H2O"    "C7H5O6"
          df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
          
          ( compd.names <- formx[compds, ]$name )
          # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
          
          
          # # O:C ratio (replace Inf with NA when C not present)
          # OC_ratio <- atomic_no(form = formx.char, atom = "O")/atomic_no(form = formx.char, atom = "C")
          # OC_ratio[is.infinite(OC_ratio)] <- NA
          # 
          # # H:C ratio (replace Inf with NA when C not present)
          # HC_ratio <- atomic_no(form = formx.char, atom = "H")/atomic_no(form = formx.char, atom = "C")
          # HC_ratio[is.infinite(HC_ratio)] <- NA
          # 
          # # N:C ratio (replace Inf with NA when C not present)
          # NC_ratio <- atomic_no(form = formx.char, atom = "N")/atomic_no(form = formx.char, atom = "C")
          # NC_ratio[is.infinite(NC_ratio)] <- NA
          
          temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
                                     cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
                                     #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
                                     )
          
        } else {
          # No specified reaction equation or chemical formula info
          df.Rxns$compds[r] <- NA
          df.Rxns$compd_coef[r] <- NA
          df.Rxns$chem_formx[r] <- NA
          
          temp.df.Compounds <- NA
         
        }
        
        if (r==1) { df.Compounds <- temp.df.Compounds }
        
        if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
        
        # clean up - if there are additional reactions?
        temp.df.Compounds <- NA
        
      } # END loop for r - rxn_id's per f/f.in
      
    } # END else loop when reactions identified
    
    # store results corresponding to each sub-reaction of each Superfocus row
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
    
  
  } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
  
  
  #return(fxn.list)
  
  saveRDS(object = fxn.list, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/fxn-list-",fxn.superfocus.rowlabel,".rds") ) # use readRDS()
  
  #print(paste0("COMPLETED ROW ",i," OF SUPERFOCUS FUNCTIONAL TAXA  # # # # # # # # # # # # # # # # # # # # #"))
  
} # END function to be run in parallel for each superfocus row


# # # # # # # # # # # # # # # # # #


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

#foreach(i=1:100 , .packages=c('stringr', 'dplyr')) %dopar%
foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
  get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )

stopCluster(cl)
time.finish <- Sys.time()



time.start #  
#[1] "2024-12-17 20:12:22 ACDT"
time.finish # 
#[1] "2024-12-17 22:46:46 ACDT"




## assemble results


modelSEED_rxn_result_dir <- "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv"


dim(df.tax)
# 30125     4




# read first output
i<-1
#temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-fxn_",i,".rds"))
fxn.superfocus.rowlabel <- row.names(df.tax)[i]
temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))

print( length(temp) )
print( names(temp) )
# "fxn_1"
print( paste0( class(temp[[1]])," ",names(temp)," of length ", length(temp[[1]])," named ",paste0(names(temp[[1]]), collapse = " & ") ))
# "list fxn_1 of length 3 named fxns & rxns & compounds"
class ( temp[[1]][["fxns"]] ) # data.frame
class ( temp[[1]][["rxns"]] ) # list
class ( temp[[1]][["compounds"]] ) # list

length( temp[[1]][["compounds"]] ) # 1
names( temp[[1]][["compounds"]] ) # NULL
class( temp[[1]][["compounds"]][[1]] ) # "logical"
is.na( temp[[1]][["compounds"]][[1]] )

i<-2
fxn.superfocus.rowlabel <- row.names(df.tax)[i]
temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))

length(temp) # 1
names(temp) # "fxn_2"
print( paste0( class(temp[[1]])," ",names(temp)," of length ", length(temp[[1]])," named ",paste0(names(temp[[1]]), collapse = " & ") ))
# "list fxn_2 of length 3 named fxns & rxns & compounds"
class ( temp[[1]][["fxns"]] ) # data.frame
class ( temp[[1]][["rxns"]] ) # list
class ( temp[[1]][["compounds"]] ) # list

length( temp[[1]][["compounds"]] ) # 1
names( temp[[1]][["compounds"]] ) # NULL
class( temp[[1]][["compounds"]][[1]] ) # "data.frame"

df.out <- temp[[1]][["compounds"]][[1]]

names(df.out) #
# [1] "superfocus_fxn" "f"              "f__in"          "rxn_id"         "cpd_id"         "cpd_name"       "cpd_form"       "cpd_molar_prop" 



# total results written to disk - 30125
dim(df.tax) # 30125     4
num_results_files <- dim(df.tax)[1]


# assemble all compound data outputs
# start with blank row

df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA #, 
                     #OC_x=NA, HC_y=NA, NC_z=NA
                     )

for (i in 1:num_results_files) {
  #i<-1
  #i<-2
  #i<-13
  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
  
  f_no <- length( temp[[1]][["compounds"]] )
  
  for (f in 1:f_no) {
    #f<-2
    # only add non-NA results
    if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
      
      df.temp <- temp[[1]][["compounds"]][[f]]
      ok <- complete.cases(df.temp)
      df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
      df.out <- rbind(df.out,df.temp)
    }
  }
 
  print(paste0("added df ",i," of ",num_results_files ))
   
}


str(df.out)
# 'data.frame':	1154548 obs. of  8 variables:

saveRDS(object = df.out, file = "df.out--get_rxns_and_compounds_indiv--sunbad-resto.RDS")
df.out <- readRDS(file = "df.out--get_rxns_and_compounds_indiv--sunbad-resto.RDS")

# remove NA first row
head(df.out)
#   superfocus_fxn  f                                                                         f__in   rxn_id   cpd_id
# 1            <NA> NA                                                                          <NA>     <NA>     <NA>
# 2           fxn_2  1                                                   2-methylaconitate isomerase rxn25278 cpd25681
# 3           fxn_2  1                                                   2-methylaconitate isomerase rxn25278 cpd02597
# 11          fxn_3  1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd24620
# 31          fxn_3  1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd25681
# 14          fxn_4  1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd01501
#                                        cpd_name cpd_form cpd_molar_prop      OC_x      HC_y NC_z
# 1                                          <NA>     <NA>             NA        NA        NA   NA
# 2                      2-methyl-trans-aconitate   C7H5O6              1 0.8571429 0.7142857    0
# 3                         cis-2-Methylaconitate   C7H5O6              1 0.8571429 0.7142857    0
# 11 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7              1 1.0000000 1.0000000    0
# 31                     2-methyl-trans-aconitate   C7H5O6              1 0.8571429 0.7142857    0
# 14                              2-Methylcitrate   C7H7O7              1 1.0000000 1.0000000    0

df.out <- df.out[-1, ]


# check for different cpd_molar_prop ??
hist(df.out$cpd_molar_prop)

dim(df.out) # 1154547       8


# normalise molar_prop to cpd_relabun so total of 1 per superfocus function !!

df.out$cpd_molar_prop_norm <- NA

length(unique(df.out$superfocus_fxn)) # 16518

phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 30125 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 30125 taxa by 4 taxonomic ranks ]

100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy)) # 54.83154 % of functions represented
100*(16518/30125) # 54.83154

fxns_found <- unique(df.out$superfocus_fxn)

for (k in 1:length(fxns_found)) {
  #k<-1
  this_fxn <- fxns_found[k]
  sel <- which(df.out$superfocus_fxn == this_fxn)
  
  sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
  # calculate 
  
  df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
  
  print(paste0("completed ",k))
  
}

sum(df.out$cpd_molar_prop_norm) # 16518


sample_sums(phy)
# 20C 30B 30A UMA 10B  5A 10C 20A 20B  5B UMB 10A  5C UMC 30C 
# 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 

dim(df.out) # 1154547       9




getwd() # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"

saveRDS(object = df.out, file = "df.out--tidy-compounds_indiv--cpp3d-sunbad-resto.RDS")

df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-sunbad-resto.RDS" )

#-------------------------


#### Sun & Badgley - post-mining
#    Collate CPP (compound rel abun %)
#    Test trending with restoration
#-------------------------

this_study <- "-sunbad-resto-"
header <- "cpp3d-cpdtrend"
phy <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")
df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-sunbad-resto.RDS" )
dim(df.out) # 1154547       9

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 30125     4


phy@sam_data
# Sample Data:        [15 samples by 8 sample variables]:
#   sample_name mgrast_id metagenome_id metagenome_name investigation_type       seq_meth file_name age
# 20C         20C mgl422999  mgm4679658.3       20C.fastq         metagenome pyrosequencing 20C.fastq  22
# 30B         30B mgl423005  mgm4679659.3       30B.fastq         metagenome pyrosequencing 30B.fastq  31
# 30A         30A mgl423002  mgm4679660.3       30A.fastq         metagenome pyrosequencing 30A.fastq  31
# UMA         UMA mgl423011  mgm4679661.3       UMA.fastq         metagenome pyrosequencing UMA.fastq  UM
# 10B         10B mgl422987  mgm4679662.3       10B.fastq         metagenome pyrosequencing 10B.fastq  12
# 5A           5A mgl422975  mgm4679663.3        5A.fastq         metagenome pyrosequencing  5A.fastq   6
# 10C         10C mgl422990  mgm4679664.3       10C.fastq         metagenome pyrosequencing 10C.fastq  12
# 20A         20A mgl422993  mgm4679665.3       20A.fastq         metagenome pyrosequencing 20A.fastq  22
# 20B         20B mgl422996  mgm4679666.3       20B.fastq         metagenome pyrosequencing 20B.fastq  22
# 5B           5B mgl422978  mgm4679667.3        5B.fastq         metagenome pyrosequencing  5B.fastq   6
# UMB         UMB mgl423014  mgm4679668.3       UMB.fastq         metagenome pyrosequencing UMB.fastq  UM
# 10A         10A mgl422984  mgm4679669.3       10A.fastq         metagenome pyrosequencing 10A.fastq  12
# 5C           5C mgl422981  mgm4679670.3        5C.fastq         metagenome pyrosequencing  5C.fastq   6
# UMC         UMC mgl423017  mgm4679671.3       UMC.fastq         metagenome pyrosequencing UMC.fastq  UM
# 30C         30C mgl423008  mgm4679672.3       30C.fastq         metagenome pyrosequencing 30C.fastq  31


# NOTE ADDITIONAL AVERAGE PLOT (n = 3) LEVEL DATA, available from Avera et al 2015
# New Forests (2015) 46:683–702
# DOI 10.1007/s11056-015-9502-8


sample_names(phy)
identical( sample_names(phy), colnames( as.matrix( phy@otu_table)) ) # TRUE

df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
dim(df.OTU) # 30125    15
df.OTU[1:5, 1:8]
# 20C          30B          30A          UMA          10B           5A          10C          20A
# fxn_1 5.233125e-04 5.972645e-05 1.551835e-04 1.195234e-04 2.501741e-04 5.509386e-04 3.930307e-04 4.019769e-04
# fxn_2 2.012740e-05 0.000000e+00 0.000000e+00 0.000000e+00 2.204177e-05 1.756616e-05 0.000000e+00 3.122151e-05
# fxn_3 1.363632e-04 2.526888e-05 4.877197e-05 2.988084e-05 1.146172e-04 1.250163e-04 9.707384e-05 1.951344e-05
# fxn_4 3.975162e-05 1.148586e-05 6.650723e-05 8.964252e-05 2.204177e-05 7.984617e-05 8.286791e-05 5.854033e-05
# fxn_5 3.371340e-03 2.216770e-03 2.194738e-03 1.942255e-03 2.182136e-03 2.626939e-03 2.379493e-03 2.566018e-03

mean( df.OTU[ , "20C"] ) # 0.003319502
sum( df.OTU[ , "20C"] ) # 100

sample_sums(phy) # all values of 100


1154547*15 # 17318205 = up to 17,318,205 rows


# HERE !!

## create df.cpd_vk_coords_long

# loop through each sample

# add grouping variables

# for each function, assign relative abundance across selected compounds

## later collapse to unique combos of "OC_x__HC_y__NC_z", add rel_abun, collate unique cpd_form, collate rxn_id, collate cpd_id 


get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
  #i<-1
  #phy_in = phy
  #dat.cpd = df.out
  
  this_samp <- sample_names(phy_in)[i]
  df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
  
  dat.cpd$sample <- this_samp
  
  dat.cpd$cpd_rel_abun_norm <- NA
  
  fxns_all <- row.names(df.OTU)
  
  for (k in 1:length(fxns_all)) {
    #k<-1
    this_fxn <- fxns_all[k]
    sel <- which(dat.cpd$superfocus_fxn == this_fxn)
    
    if (length(sel)>=1) {
      dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
      
    }
  } # END rel abun values for all relevant functions added
  
  saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  
} # END


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)

stopCluster(cl)
time.finish <- Sys.time()

# output 1
i<-1
this_samp <- sample_names(phy)[i]
#saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
dat <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
head(dat)

for ( i in 2:length(sample_names(phy)) ) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  #saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  temp <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  
  dat <- rbind(dat, temp)
  
  print(paste0("completed ",i))
}


saveRDS(object = dat, file = "dat.cpd-long-all-samps-cpp3d-sunbad-resto.rds" )
dat <- readRDS("dat.cpd-long-all-samps-cpp3d-sunbad-resto.rds")

rm(temp)

str(dat)
# 'data.frame':	17318205 obs. of  11 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id             : chr  "rxn25278" "rxn25278" "rxn25279" "rxn25279" ...
# $ cpd_id             : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ cpd_name           : chr  "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" ...
# $ cpd_form           : chr  "C7H5O6" "C7H5O6" "C7H7O7" "H2O" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.333 0.333 0.333 ...
# $ sample             : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun_norm  : num  1.01e-05 1.01e-05 4.55e-05 4.55e-05 4.55e-05 ...

#dat$combos <- paste0(dat$OC_x,"__",dat$HC_y,"__",dat$NC_z)

sum(dat$cpd_rel_abun_norm) # 1005.557
sum(dat$cpd_rel_abun_norm)/nsamples(phy) # 67.03712 = average 67% functional relative abundance per sample

length(which(is.na(dat$cpd_rel_abun_norm))==TRUE) # 0
length(which( dat$cpd_rel_abun_norm > 0) == TRUE) # 11960630
length(which( dat$cpd_rel_abun_norm == 0) == TRUE) # 5357575

# so this step does collect some zero relative abundances from the 'otu-table' 

names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"               "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm" "sample"              "cpd_rel_abun_norm" 

## create dat.cpd.distil
## later collapse to unique combos of "OC_x__HC_y__NC_z", add rel_abun, collate unique cpd_form, collate rxn_id, collate cpd_id 

length(unique(dat$cpd_id)) # 8370



## Collate compounds within each sample 


unique_cpd <- unique(dat$cpd_id)
samp_names <- sample_names(phy)




collate_compounds <- function(dat.cpd, unique_cpd, samp) {
  #i<-1
  #samp = samp_names[i]
  #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
  
  this_samp <- samp
  
  cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, #OC_x=NA, HC_y=NA, NC_z=NA, 
                         cpd_rel_abun=NA)
  
  for (c in 1:length(unique_cpd)) {
    #c<-1
    this_cpd <- unique_cpd[c]
    sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
    
    # cpd_data$OC_x[c] <- unique(dat.cpd$OC_x[sel.cpd])
    # cpd_data$HC_y[c] <- unique(dat.cpd$HC_y[sel.cpd])
    # cpd_data$NC_z[c] <- unique(dat.cpd$NC_z[sel.cpd])
    
    # # now collate rel abun for this sample only
    # sel <- which(dat.cpd$sample == this_samp)
    # dat.cpd <- dat.cpd[sel, ]
    
    # # again select in this sample-specific dataset
    # sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
    
    if (length(sel.cpd) >=1) {
      cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
    }
    
  } # END all compounds
  
  saveRDS(object = cpd_data, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
  
} # END


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])

stopCluster(cl)
time.finish <- Sys.time()



# output 1
i<-1
this_samp <- sample_names(phy)[i]
dat.cpd.collate <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
head(dat.cpd.collate)

for ( i in 2:length(sample_names(phy)) ) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  temp <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
  
  dat.cpd.collate <- rbind(dat.cpd.collate, temp)
  
  print(paste0("completed ",i))
}


str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  3 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...

sum(dat.cpd.collate$cpd_rel_abun) # 1005.557.   929.1978
sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample)) # 67.03712.   61.94652

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-sunbad-resto.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-sunbad-resto.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000013 0.000165 0.008009 0.001455 8.107427
# PREVIOUSLY, BASED ON C,O,H,N
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000020 0.000210 0.008008 0.002060 4.046109 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf -4.9017 -3.7832    -Inf -2.8372  0.9089
# PREVIOUSLY, BASED ON C,O,H,N
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -4.693  -3.678    -Inf  -2.686   0.607


# log10 abun
dat.cpd.collate$log10_abun <- dat.cpd.collate$cpd_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(dat.cpd.collate$log10_abun == 0) # qty 8897 7945
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(dat.cpd.collate$log10_abun[ -subsel.zero ])
  dat.cpd.collate$log10_abun[ subsel.zero ] <- zero_replace
}
dat.cpd.collate$log10_abun <- log10(dat.cpd.collate$log10_abun)

hist(dat.cpd.collate$log10_abun); summary( dat.cpd.collate$log10_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.9875 -4.9017 -3.7832 -4.1549 -2.8372  0.9089


# make group variable from sample name

dat.cpd.collate$group <- NA


for (i in 1:length(sample_names(phy))) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group) # "22" "31" "UM" "12" "6"
dat.cpd.collate$group <- factor(dat.cpd.collate$group, levels = c("6", "12", "22", "31", "UM"), ordered = TRUE)

dat.cpd.collate$group_label <- factor(dat.cpd.collate$group, 
                                      levels = c("6","12", "22", "31", "UM"),
                                      labels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"),ordered = TRUE)



levels(dat.cpd.collate$group) # "6"  "12" "22" "31" "UM"

dat.cpd.collate$ord_group <- NA
sel <- which(dat.cpd.collate$group == "6") # qty Compound: 25110   23208  Indiv 31500 ; Fxn-MEan 12240
dat.cpd.collate$ord_group[sel] <- 1
sel <- which(dat.cpd.collate$group == "12") # qty 25110
dat.cpd.collate$ord_group[sel] <- 2
sel <- which(dat.cpd.collate$group == "22") # qty 25110
dat.cpd.collate$ord_group[sel] <- 3
sel <- which(dat.cpd.collate$group == "31") # qty 25110
dat.cpd.collate$ord_group[sel] <- 4
sel <- which(dat.cpd.collate$group == "UM") # qty 25110
dat.cpd.collate$ord_group[sel] <- 5


saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")



str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...
# $ log10_abun  : num  -3.362 -1.643 -3.341 0.711 -1.804 ...
# $ group       : Ord.factor w/ 5 levels "6"<"12"<"22"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ group_label : Ord.factor w/ 5 levels "6 yr"<"12 yr"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ ord_group   : num  3 3 3 3 3 3 3 3 3 3 ...


length( unique(dat.cpd.collate$cpd_id) ) # 8370
8370*15 # 125550


dat.test <- data.frame(cpd = unique(dat.cpd.collate$cpd_id), data_for_this_cpd=NA , p_val = NA, kendall_tau = NA, trend_with_age = NA )


for (i in 1:dim(dat.test)[1]) {
  #i<-1
  this_cpd <- dat.test$cpd[i]
  sel <- which(dat.cpd.collate$cpd_id == this_cpd)
  # Kendall Tau correlation
  
  x = dat.cpd.collate$ord_group[sel]
  y = dat.cpd.collate$log10_abun[sel]
  df = as.data.frame(cbind(x,y))
  sel.ok = which(complete.cases(df)==TRUE)
  df = df[sel.ok, ]
  if (sd(df$y)==0 & length(unique(df$y))==1 & unique(df$y)[1]==min(df$y)) {  
    # disqualified: if they were all zeros before zero-replacement
    dat.test$data_for_this_cpd[i] <- "disqualified"
    
  } else if ( length(which(df$y > min(df$y))) < 0.25*length(sel) ) { # i.e. 3.75 ; 3 or less
    # low data: if non-zero cases do not comprise at least a replicate / a quarter of available scenario data
    dat.test$data_for_this_cpd[i] <- "low data"
    
  } else {
    ktcor<- cor.test(x = df$x, y = df$y, method = "kendall")
    dat.test$p_val[i] <- ktcor$p.value
    dat.test$kendall_tau[i] <- ktcor$statistic
    
  }
  
  # x = df.heat3d$ord_group[sel]
  # y = df.heat3d$log10_abun[sel]
  # df = as.data.frame(cbind(x,y))
  # ok = which(complete.cases(df)==TRUE)
  # if (length(ok)>1) { df = df[ok, ]}
  # plot(df$x, df$y)
  
  if (!(is.na(dat.test$p_val[i])|is.na(dat.test$kendall_tau[i]))) {
    if (dat.test$p_val[i] <= 0.05 & dat.test$kendall_tau[i] > 0) { dat.test$trend_with_age[i] <- "Increasing" }
    if (dat.test$p_val[i] <= 0.05 & dat.test$kendall_tau[i] < 0) { dat.test$trend_with_age[i] <- "Decreasing" }
  }
  print(paste0("Completed ",i))
}

sel.disq <- which(dat.test$data_for_this_cpd == "disqualified") # empty
sel.low <- which(dat.test$data_for_this_cpd == "low data") # 381

length(unique(dat.cpd.collate$cpd_id)) # 8370
length(unique(dat.cpd.collate$cpd_id)) - (length(sel.disq) + length(sel.low) ) # 7989

sel.nonNA <- which(!is.na(dat.test$p_val)) # 7989 applicable tests

sel.sig <- which(dat.test$p_val <= 0.05) # 2958

# only keep applicable tests; the remainder are likely due to zero replacement:
dat.test <- dat.test[ sel.nonNA, ]


summary(dat.test$p_val)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000133 0.0113463 0.1879858 0.2947939 0.5434239 1.0000000

hist(dat.test$p_val)


# ## False Discovery Rate correction (Benjamini & Hochberg 1995)
# http://www.statisticshowto.com/benjamini-hochberg-procedure/

# what is m?

m <- dim(dat.test)[1]  # 7989 applicable compounds tested
alpha <- 0.05

p_values <- dat.test[ order(dat.test$p_val, decreasing = FALSE) ,  "p_val" ]

summary(p_values)


plot(x=1:m, y=p_values, xlab="k", ylab="P(k)")

# criteria values for Benjamini & Hochberg test
test_values <- rep(NA, times=length(p_values))

for (i in 1:m) { test_values[i] <- (i/m)*alpha }
test_values

test_results <- rep(NA, times=length(p_values))

## calculate test results
for (i in 1:m) {
  if ( p_values[i] < test_values[i] ) { test_results[i] <- "yes" }
}
test_results
# get index of largest ranked p-value with "yes" result (i.e. p-value is smaller than test criteria)

idx <- rev(which(!is.na(test_results)))[1] # 2122.  compounds: 1953 ; cf  313

dat.test$sigBH <- NA
if (length(1:idx) >0) {
  dat.test[ order(dat.test$p_val,decreasing = FALSE)[1:idx] , "sigBH" ] <- "sig"
}

length(1:idx) # 2122.  compounds 1953 ;  heatmap cell 313


p_values[idx] # 0.01311466;  0.01219892 ;  0.01984315


plot(x=1:m, y=p_values, xlab="k (index of ranked P-values)", ylab="P-value(k)" ) #, xlim=c(0,100), ylim=c(0,0.05))
# title("(a) 16S", adj=0)
abline(a=0, b=(alpha/m), col="red" )
text(x = 4000, y = 0.08, adj = 0,labels = "slope = alpha/m", col = "red")

points(x=c(1:m)[1:idx], y=p_values[1:idx], col="purple" )
text(x = 1000, y = 0.1, adj=0.5, labels = "P(k) < \n(k/m)*alpha", col = "purple")


dev.print(tiff, filename = paste0(workdir,"/plots/","3d-indiv-compound-Benjamini-Hochberg-significant-p-values-",this_study,header,".tiff"),
          width = 14, height = 14, units = "cm", res=600, compression = "lzw",type="cairo" )


# extract sig results

#sel.sig <- which(dat.test$sigBH == "sig") # 60
sel.sig <- which(dat.test$p_val <= 0.05) # 2958.  compounds: 2742  ; heatmap cells 381

dat.test.sig <- dat.test[sel.sig, ]

dat.test.sig$minuslog10_p_val <- -log10(dat.test.sig$p_val)

plot(x = dat.test.sig$kendall_tau , y =dat.test.sig$minuslog10_p_val , xlab="Kendall Tau", ylab="-log10(P-value)")

sel.sigBH <- which(dat.test.sig$sigBH == "sig")
points(x=dat.test.sig$kendall_tau[sel.sigBH], y=dat.test.sig$minuslog10_p_val[sel.sigBH], col="purple" )

dev.print(tiff, filename = paste0(workdir,"/plots/","3d-indiv-compound-VolcanoPlot-Benjamini-Hochberg-significant-p-values-",this_study,header,".tiff"),
          width = 12, height = 14, units = "cm", res=600, compression = "lzw",type="cairo" )


# join compound info??

length(which(dat.test.sig$sigBH=="sig")) # 2122


dat.test.sig$cpd_names <- NA
dat.test.sig$cpd_forms <- NA

dat.test.sig$OC_x <- NA
dat.test.sig$HC_y <- NA
dat.test.sig$NC_z <- NA

dat.test.sig$mass <- NA
dat.test.sig$class <- NA

for (i in 1:dim(dat.test.sig)[1]) {
  #i<-1
  this_cpd <- dat.test.sig$cpd[i]
  
  sel.cpd <- which(df.comp2$id == this_cpd)
  
  dat.test.sig$cpd_names[i] <- df.comp2$name[sel.cpd]
  dat.test.sig$cpd_forms[i] <- df.comp2$form[sel.cpd]
  
  dat.test.sig$OC_x[i] <- df.comp2$OC_ratio[sel.cpd]
  dat.test.sig$HC_y[i] <- df.comp2$HC_ratio[sel.cpd]
  dat.test.sig$NC_z[i] <- df.comp2$NC_ratio[sel.cpd]
  
  dat.test.sig$mass[i] <- df.comp2$mass[sel.cpd]
  dat.test.sig$class[i] <- df.comp2$class[sel.cpd]
  
  print(paste0("completed ",i))
}

#write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.rds")
saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.rds")
# updated save later



## plot as Increasing or Decreasing?? in vK space


names(dat.test.sig)
# [1] "cpd"               "data_for_this_cpd" "p_val"             "kendall_tau"       "trend_with_age"    "sigBH"             "minuslog10_p_val" 
# [8] "cpd_names"         "cpd_forms"         "OC_x"              "HC_y"              "NC_z"


p <- ggplot(data = filter(dat.test.sig, sigBH == "sig")) +
  coord_equal()+
  ggtitle("Microbiota compound processing potential\n- Post-mining restoration case study")+
  xlim(0,3.4)+ ylim(0,4.1)+
  geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with age\nin functional\ncapacity (%)\nallocated to\ncompounds"))+
  
  geom_mark_rect(data= vkgrouprect, aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ 
  
  annotate(geom="text", x= 0+0.01, y= 2.3+0.02, label = "Lipid", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.2+0.01, y= 2.2+0.02, label = "Protein", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  #annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373", lineheight = 0.8) + # top-left
  annotate(geom="text", x= 0.7+0.01, y= 2.4+0.02, label = "Carbohydrate", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0+0.01, y= 0.5-0.02, label = "Condensed\naromatics", hjust=0, vjust=1, size = 2 , col="#737373", lineheight = 0.8) + # bottom-left
  annotate(geom="text", x= 0.25+0.01, y= 0.75-0.02, label = "Lignin", hjust=0, vjust=1, size = 2 , col="#737373" ) + # bottom-left
  annotate(geom="text", x= 0.67+0.01, y= 0.53-0.02, label = "Tannin", hjust=0, vjust=1, size = 2 , col="#737373") + # bottom-left
  
  theme_bw()+
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.9)) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","Compounds-indiv-vKSpace-Trend-with-Age-",this_study,header,".tiff"), width = 14, height = 14, units = "cm", res=350, compression="lzw",type="cairo")
# Removed 147 rows containing missing values or values outside the scale range (`geom_point()`). 


hist(dat.test.sig$NC_z)


min(dat.test.sig$NC_z[ dat.test.sig$NC_z > 0 ], na.rm = TRUE) # 0.01449275


dim(dat.test.sig) # 2958   15
sel.ok <- which(!is.na(dat.test.sig$NC_z) ) # qty 2757

dat.test.sig$z_layer <- NA

subsel <- which(dat.test.sig$NC_z[sel.ok] == 0) # 1184
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C = 0"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0 & dat.test.sig$NC_z[sel.ok] <= 0.2 ) # 845
# Rivas-Ubach et al 2018 Table 1 highlights N/C breaks at 0.126, 0.2, 0.5, 0.7
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0 to 0.2"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0.2) # 728
max(dat.test.sig$NC_z[sel.ok]) # 2
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 2"

unique(dat.test.sig$z_layer[sel.ok]) # "N:C = 0"       "N:C >0 to 0.2" "N:C >0.2 to 2"

dat.test.sig$z_layer <- factor(dat.test.sig$z_layer, levels = c("N:C = 0",
                                                                "N:C >0 to 0.2",
                                                                "N:C >0.2 to 2"), ordered = TRUE)


write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.rds")
saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.rds")
#dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.rds")

dim(dat.test.sig) #  2958   15
head(dat.test.sig)
dim(dat.test.sig[ which(dat.test.sig$sigBH == "sig"), ]) # 2122   15
sel <- which(dat.test.sig$sigBH == "sig" & !is.na(dat.test.sig$OC_x) ) # 1975
sel <- which(dat.test.sig$sigBH == "sig" & dat.test.sig$trend_with_age == "Decreasing") # 1365
sel <- which(dat.test.sig$sigBH == "sig" & dat.test.sig$trend_with_age == "Decreasing" & !is.na(dat.test.sig$OC_x) ) # 1255
sel <- which(dat.test.sig$sigBH == "sig" & dat.test.sig$trend_with_age == "Increasing") # 757
sel <- which(dat.test.sig$sigBH == "sig" & dat.test.sig$trend_with_age == "Increasing" & !is.na(dat.test.sig$OC_x) ) # 720

# build vk layer for each z-layer facet AND age group
for (i in 1:length(levels(dat.test.sig$z_layer))) {
  #i<-1
  temp <- vkgrouprect
  temp$z_layer <- levels(dat.test.sig$z_layer)[i]
  if (i==1) { keep <- temp }
  if (i>1) { keep <- rbind(keep, temp)}
  print(paste0("completed ",i))
}
vkgrouprect.facets <- keep
rm(keep)

str( vkgrouprect.facets )
vkgrouprect.facets$z_layer <- factor( vkgrouprect.facets$z_layer )



p <- ggplot(data = filter(dat.test.sig[sel.ok, ], sigBH == "sig")) +
  coord_equal()+
  ggtitle("Compound processing potential of microbiota - Post-mining restoration case study")+
  #xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  xlab("O:C ratio")+ ylab("H:C ratio")+
  guides(color = guide_legend(title = "Trend with age in functional capacity\n(%) allocated to compounds"))+
  
  facet_wrap(facets = vars(z_layer))+
  
  #geom_mark_rect(data= vkgrouprect, aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ 
  
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  annotate(geom="text", x= 0+0.01, y= 2.3+0.02, label = "Lipid", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.2+0.01, y= 2.2+0.02, label = "Protein", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  #annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373", lineheight = 0.8) + # top-left
  annotate(geom="text", x= 0.7+0.01, y= 2.4+0.02, label = "Carbohydrate", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0+0.01, y= 0.5-0.02, label = "Condensed\naromatics", hjust=0, vjust=1, size = 2 , col="#737373", lineheight = 0.8) + # bottom-left
  annotate(geom="text", x= 0.25+0.01, y= 0.75-0.02, label = "Lignin", hjust=0, vjust=1, size = 2 , col="#737373" ) + # bottom-left
  annotate(geom="text", x= 0.67+0.01, y= 0.53-0.02, label = "Tannin", hjust=0, vjust=1, size = 2 , col="#737373") + # bottom-left
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.9)) ,
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Age-",this_study,header,"-v2.tiff"), width = 20, height = 12, units = "cm", res=350, compression="lzw",type="cairo")


subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0.2 to 2" & dat.test.sig$sigBH[sel.ok] == "sig" &
                  dat.test.sig$OC_x[sel.ok] > 0.52 & dat.test.sig$OC_x[sel.ok] <= 1.05 &
                  dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.2 )
dat.test.sig[sel.ok[subsel] , ] # Proteins!
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Nucleotide     Protein Unspecified 
# 9         122           8
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" &
                  dat.test.sig$OC_x[sel.ok] > 0.25 & dat.test.sig$OC_x[sel.ok] <= 1.05 &
                  #dat.test.sig$OC_x[sel.ok] > 0.52 & dat.test.sig$OC_x[sel.ok] <= 1.05 &
                  dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.2 )
#dat.test.sig[sel.ok[subsel] , ] # 
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Amino sugar       Lipid     Protein Unspecified 
# 19         116          71          94
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" &
                  dat.test.sig$OC_x[sel.ok] > 0.35 & dat.test.sig$OC_x[sel.ok] <= 0.61 &
                  #dat.test.sig$OC_x[sel.ok] > 0.61 & dat.test.sig$OC_x[sel.ok] <= 1.05 & # 0.52  1.05
                  dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.2 )
#dat.test.sig[sel.ok[subsel] , ] # 
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid     Protein Unspecified 
# 46          11          14

## Check dominant compound classes in each vK zone?
## "N:C = 0"
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C = 0" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0 & dat.test.sig$OC_x[sel.ok] <= 0.6 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid Unspecified - ok Lipids!
# 58         251 

subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C = 0" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.7 & dat.test.sig$OC_x[sel.ok] <= 1.25 & dat.test.sig$HC_y[sel.ok] > 1.5 & dat.test.sig$OC_x[sel.ok] <= 2.4 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Carbohydrate  Unspecified - ok Carbohydrates
# 92           19 

# "N:C >0 to 0.2"
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0 & dat.test.sig$OC_x[sel.ok] <= 0.5 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid     Protein Unspecified 
# 143          55          66 
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0 & dat.test.sig$OC_x[sel.ok] <= 0.6 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid     Protein Unspecified - Lipids with Protein overlap
# 147          85          74
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.12 & dat.test.sig$OC_x[sel.ok] <= 0.6 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid     Protein Unspecified 
# 147          85          58 
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.12 & dat.test.sig$OC_x[sel.ok] <= 1.05 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Amino sugar       Lipid     Protein Unspecified 
#          19         147          85         121

# "N:C >0.2 to 2"
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0.2 to 2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.3 & dat.test.sig$OC_x[sel.ok] <= 1.05 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Nucleotide     Protein Unspecified 
# 9         191          17 
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0.2 to 2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.12 & dat.test.sig$OC_x[sel.ok] <= 1.05 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Nucleotide     Protein Unspecified 
# 9         220          20



## Use adjusted compound classes (adapted from Wu 2018, D'Andrilli, Rivas-Ubach 2018, and Minor et al 2015)

vkgrouprect.facets2 <- read.table(file = "cpp3d-compound-classes.tsv", header = TRUE, sep = "\t" )
vkgrouprect.facets2.labels <- read.table(file = "cpp3d-compound-classes-labels.tsv", header = TRUE, sep = "\t" )


p <- ggplot(data = filter(dat.test.sig[sel.ok, ], sigBH == "sig")) +
  coord_equal()+
  
  #ggtitle("Compound processing potential of microbiota - Post-mining restoration case study")+
  
  ##xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age, shape = class), size = 1, alpha = 0.3 ) + # 
  #scale_shape_manual(values = shapes.class, name = "Compound type")+
  
  xlab("O:C ratio")+ ylab("H:C ratio")+
  guides(color = guide_legend(title = "Trend with age in functional capacity\n(%) allocated to compounds"))+
  
  facet_wrap(facets = vars(z_layer))+
  
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T,  size = 2.5 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.5 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.5 , col="#737373" , lineheight = 0.8)+
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.9)) ,
    
    legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"),
    title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
  )

p
#dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Age-",this_study,header,"-v2c.tiff"), width = 20, height = 12, units = "cm", res=350, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Age-Sunbad-resto-NoGGtitle-v2c.tiff"), width = 20, height = 12, units = "cm", res=600, compression="lzw",type="cairo")





### also visualise all compounds for a single sample

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

# isolate a single sample
df1 <- dat.cpd.collate
unique(df1$sample)
#  "20C" "30B" "30A" "UMA" "10B" "5A"  "10C" "20A" "20B" "5B"  "UMB" "10A" "5C"  "UMC" "30C"
sel <- which(df1$sample == "UMA")
df1 <- df1[sel, ]
names(df1)
# [1] "cpd_id"       "sample"       "cpd_rel_abun" "log10_abun"   "group"        "group_label"  "ord_group"    "cpd_names"    "cpd_forms"   
# [10] "OC_x"         "HC_y"         "NC_z" 

df1$cpd_names <- NA
df1$cpd_forms <- NA
df1$OC_x <- NA
df1$HC_y <- NA
df1$NC_z <- NA

for (i in 1:dim(df1)[1]) {
  #i<-1
  this_cpd <- df1$cpd_id[i]
  
  sel.cpd <- which(df.comp2$id == this_cpd)
  
  df1$cpd_names[i] <- df.comp2$name[sel.cpd]
  df1$cpd_forms[i] <- df.comp2$form[sel.cpd]
  
  df1$OC_x[i] <- df.comp2$OC_ratio[sel.cpd]
  df1$HC_y[i] <- df.comp2$HC_ratio[sel.cpd]
  df1$NC_z[i] <- df.comp2$NC_ratio[sel.cpd]
  
  print(paste0("completed ",i))
}

hist(df1$NC_z)


min(df1$NC_z[ df1$NC_z > 0 ], na.rm = TRUE) # 0.004201681

dim(df1) # 8370 12
dim(df1[df1$cpd_rel_abun > 0, ]) # 7826   13

sel.ok <- which(!is.na(df1$NC_z) ) # qty 7736

df1$z_layer <- NA

subsel <- which(df1$NC_z[sel.ok] == 0) # 3661
df1$z_layer[sel.ok [subsel]] <- "N:C = 0"

subsel <- which(df1$NC_z[sel.ok] > 0 & df1$NC_z[sel.ok] <= 0.2 ) # 2290
# Rivas-Ubach et al 2018 Table 1 highlights N/C breaks at 0.126, 0.2, 0.5, 0.7
df1$z_layer[sel.ok [subsel]] <- "N:C >0 to 0.2"

subsel <- which(df1$NC_z[sel.ok] > 0.2) # 1785
max(df1$NC_z[sel.ok]) # 3
#df1$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 2"
df1$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 3"

unique(df1$z_layer[sel.ok]) # "N:C = 0"       "N:C >0 to 0.2" "N:C >0.2 to 3"

df1$z_layer <- factor(df1$z_layer, levels = c("N:C = 0",
                                                                "N:C >0 to 0.2",
                                                                "N:C >0.2 to 3"), ordered = TRUE)


# compound class guides
vkgrouprect.facets3 <- vkgrouprect.facets2
vkgrouprect.facets3.labels <- vkgrouprect.facets2.labels

sel <- which(vkgrouprect.facets3$z_layer == "N:C >0.2 to 2")
vkgrouprect.facets3$z_layer[sel] <- "N:C >0.2 to 3"

sel <- which(vkgrouprect.facets3.labels$z_layer == "N:C >0.2 to 2")
vkgrouprect.facets3.labels$z_layer[sel] <- "N:C >0.2 to 3"



p <- ggplot(data = df1[sel.ok, ]) +
  coord_equal()+
  
  #ggtitle("Compound processing potential of microbiota - Post-mining restoration case study")+
  
  ##xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age, shape = class), size = 1, alpha = 0.3 ) + # 
  #scale_shape_manual(values = shapes.class, name = "Compound type")+
  
  xlab("O:C ratio")+ ylab("H:C ratio")+
  
  #guides(color = guide_legend(title = "Log10 functional capacity (%)\nallocated to compounds"))+
  
  scale_color_continuous(type = "viridis", direction = -1) + 
  #guides(color = guide_colorbar(direction = "horizontal", barheight = 0.65, barwidth = 3.5, vjust = 0, hjust = 0.5, title = "Log10 functional capacity (%)\nallocated to compounds"))+
  guides(color = guide_colorbar(direction = "horizontal", barheight = 0.65, barwidth = 3.5, vjust = 0, hjust = 0.5, title = "Log10 functional rel abun (%)\nallocated to compounds"))+   
  
  facet_wrap(facets = vars(z_layer))+
  
  geom_mark_rect(data= filter(vkgrouprect.facets3, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets3, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets3, z_layer == "N:C >0.2 to 3" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  geom_point(aes(x = OC_x, y = HC_y, color = log10_abun), size = 1, alpha = 0.3 ) + # 
  
  geom_text(data = filter(vkgrouprect.facets3.labels, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T,  size = 2.5 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets3.labels, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.5 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets3.labels, z_layer == "N:C >0.2 to 3" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.5 , col="#737373" , lineheight = 0.8)+
  
  theme_bw()+
  theme(
    #legend.position = "right",
    
    legend.position = c(0.5, 0.05), # as fraction of plot dimensions
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(0.7)),
    
    #legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.background = element_rect(fill = "transparent"),
    
    #legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.6)) ,
    
    legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"),
    #title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
  )

p
#dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Age-",this_study,header,"-v2c.tiff"), width = 20, height = 12, units = "cm", res=350, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-1Sample-UMA-Sunbad-resto-NoGGtitle-v2c.tiff"), width = 20, height = 12, units = "cm", res=450, compression="lzw",type="cairo")




#-------------------------

#### Sun & Badgley - post-mining - CPP as phyloseq object
#    PCoA using CPP vs Functions
#-------------------------

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

data_in <- dat.cpd.collate

length( unique(data_in$cpd_id) ) # 8370
8370*15 # 125550


### get data into phyloseq object ...

head(data_in)
#     cpd_id sample cpd_rel_abun log10_abun group group_label ord_group
# 1 cpd25681    20C 0.0004345842 -3.3619261    22       22 yr         3
# 2 cpd02597    20C 0.0227382574 -1.6432428    22       22 yr         3
# 3 cpd24620    20C 0.0004564208 -3.3406346    22       22 yr         3
# 4 cpd00001    20C 5.1354106008  0.7105752    22       22 yr         3
# 5 cpd01501    20C 0.0157048650 -1.8039658    22       22 yr         3
# 6 cpd00851    20C 0.0121620567 -1.9149930    22       22 yr         3

df.wide <- dcast(data_in, formula = sample + group ~ cpd_id , value.var = "cpd_rel_abun" )

df.wide[1:5, 1:10]
#   sample group cpd00001 cpd00002 cpd00003 cpd00004  cpd00005  cpd00006  cpd00007 cpd00008
# 1    10A    12 5.126634 2.287831 1.058127 1.031673 0.7565038 0.7592450 0.3678576 1.393837
# 2    10B    12 5.109487 2.333880 1.058379 1.032674 0.7348372 0.7370765 0.3472805 1.392368
# 3    10C    12 5.103210 2.319224 1.070352 1.044885 0.7498325 0.7520871 0.3398280 1.399435
# 4    20A    22 5.138510 2.323322 1.069859 1.043509 0.7560155 0.7585376 0.3598548 1.399336
# 5    20B    22 5.173074 2.355075 1.064497 1.037668 0.7483448 0.7511583 0.3623488 1.419401

unique(paste0(df.wide$sample,"--",df.wide$group))
# [1] "10A--12" "10B--12" "10C--12" "20A--22" "20B--22" "20C--22" "30A--31" "30B--31" "30C--31" "5A--6"  
# [11] "5B--6"   "5C--6"   "UMA--UM" "UMB--UM" "UMC--UM"

# save group variable
samp <- df.wide[ ,1:2]
row.names(samp) <- samp$sample

# transpose
df.wide <- t(df.wide[ ,-2]) # minus 'group' column

head(df.wide)
# [,1]        [,2]        [,3]        [,4]        [,5]        [,6]        [,7]       
# sample   "10A"       "10B"       "10C"       "20A"       "20B"       "20C"       "30A"      
# cpd00001 "5.126634"  "5.109487"  "5.103210"  "5.138510"  "5.173074"  "5.135411"  "5.036718" 
# cpd00002 "2.287831"  "2.333880"  "2.319224"  "2.323322"  "2.355075"  "2.318166"  "2.309780" 
# cpd00003 "1.058127"  "1.058379"  "1.070352"  "1.069859"  "1.064497"  "1.060200"  "1.067475" 
# cpd00004 "1.031673"  "1.032674"  "1.044885"  "1.043509"  "1.037668"  "1.033088"  "1.041951" 
# cpd00005 "0.7565038" "0.7348372" "0.7498325" "0.7560155" "0.7483448" "0.7564190" "0.7424559"
# [,8]        [,9]        [,10]       [,11]       [,12]       [,13]       [,14]      
# sample   "30B"       "30C"       "5A"        "5B"        "5C"        "UMA"       "UMB"      
# cpd00001 "4.977622"  "5.054107"  "5.136126"  "5.128531"  "5.155902"  "4.988787"  "5.045806" 
# cpd00002 "2.272223"  "2.332804"  "2.335298"  "2.345242"  "2.334887"  "2.298596"  "2.329824" 
# cpd00003 "1.093609"  "1.077584"  "1.054323"  "1.056899"  "1.075898"  "1.098387"  "1.096881" 
# cpd00004 "1.068711"  "1.051745"  "1.026151"  "1.028856"  "1.047329"  "1.072772"  "1.071776" 
# cpd00005 "0.7597370" "0.7483150" "0.7482304" "0.7517447" "0.7519912" "0.7931363" "0.7843985"
# [,15]      
# sample   "UMC"      
# cpd00001 "5.064984" 
# cpd00002 "2.339306" 
# cpd00003 "1.089745" 
# cpd00004 "1.065067" 
# cpd00005 "0.7746962"

samp_names <- df.wide[1, ]
tax_names <- row.names(df.wide[-1, ])
head(tax_names) # "cpd00001" "cpd00002" "cpd00003" "cpd00004" "cpd00005" "cpd00006"
otu.df <- df.wide[-1, ] # remove sample labels in 1st row
# this is necessary to create numeric matrix

colnames(otu.df) <- samp_names

# convert OTU table to matrix
class(otu.df) # "matrix" "array" 
#otu.df <- as.matrix(otu.df)

# convert to numeric matrix
# https://stackoverflow.com/questions/20791877/convert-character-matrix-into-numeric-matrix
otu.df <- apply(otu.df, 2, as.numeric)

rownames(otu.df) # NULL
dim(otu.df) #  8370   15
rownames(otu.df) <- tax_names

## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
OTU <- otu_table(otu.df, taxa_are_rows = TRUE)



# # convert Taxonomy table to matrix  

tax <- data.frame(cpd_id = tax_names)
row.names(tax) <- tax_names

tax <- as.matrix(tax)

identical( row.names(otu.df), row.names(tax) ) # TRUE


## Create 'taxonomyTable'
#  tax_table - Works on any character matrix.
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
TAX <- tax_table(tax)


## Create a phyloseq object, merging OTU & TAX tables
phy.cpp = phyloseq(OTU, TAX)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8370 taxa and 15 samples ]
# tax_table()   Taxonomy Table:    [ 8370 taxa by 1 taxonomic ranks ]


sample_names(phy.cpp)
# "10A" "10B" "10C" "20A" "20B" "20C" "30A" "30B" "30C" "5A"  "5B"  "5C"  "UMA" "UMB" "UMC"

identical(sample_names(phy.cpp), samp$sample) # TRUE


# row.names need to match sample_names() from phyloseq object
row.names(samp) <- samp$sample




### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

SAMP <- sample_data(samp)


### Combine SAMPDATA into phyloseq object
phy.cpp <- merge_phyloseq(phy.cpp, SAMP)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8370 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 8370 taxa by 1 taxonomic ranks ]

phy.cpp@sam_data
# Sample Data:        [15 samples by 2 sample variables]:
#   sample group
# 10A    10A    12
# 10B    10B    12
# 10C    10C    12
# 20A    20A    22
# 20B    20B    22
# 20C    20C    22
# 30A    30A    31
# 30B    30B    31
# 30C    30C    31
# 5A      5A     6
# 5B      5B     6
# 5C      5C     6
# UMA    UMA    UM
# UMB    UMB    UM
# UMC    UMC    UM




phy_in <- phy.cpp

phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8370 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 8370 taxa by 1 taxonomic ranks ]

min(taxa_sums(phy_in)) # 8.136846e-08
sum(sample_sums(phy_in)) # 1005.557
sample_sums(phy_in)
# 10A      10B      10C      20A      20B      20C      30A      30B      30C       5A       5B 
# 66.37627 67.28058 67.43371 67.14898 67.44508 67.38783 67.00958 66.68865 67.16063 67.19161 67.09898 
# 5C      UMA      UMB      UMC 
# 67.64787 66.59001 66.44038 66.65667

summary( sample_sums(phy_in) )
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 66.38   66.67   67.15   67.04   67.33   67.65 

sd( sample_sums(phy_in) )
# 0.3958827

max(taxa_sums(phy_in)) # 120.7394


# don't rarefy - already in form of relative abundance %


table(phy_in@sam_data$group)
# 6 12 22 31 UM 
# 3  3  3  3  3



## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
#ord <- ordinate(r1.ps, "NMDS", "bray")
ord <- ordinate(phy_in, "PCoA", "bray")


ord

unique(phy_in@sam_data$group)
# [1] 12 22 31 6  UM
# Levels: 6 < 12 < 22 < 31 < UM


p <- plot_ordination(phy_in, ord, type="samples", color="group")
p

p$labels$x # "Axis.1   [61.5%]"
x_lab <- "PCo1 (61.5%)"

p$labels$y # "Axis.2   [17.9%]"
y_lab <- "PCo2 (17.9%)"

#temp <- r1.ps
p_df <- p$data


cols.group <- c("6" = "#9e0142",
                "12" = "#d53e4f",
                "22" = "#fdae61",
                "31" = "#abdda4",
                "UM" = "#3288bd")



p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = group))+
  theme_bw()+
  geom_point()+
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Reveg\nage (yr)") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "Restoration: CPP\n\n", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-sunbad-resto-v2.tiff"), width = 9, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# "sample" "group" 

# Adonis test
set.seed(123)
adonis2(bray ~ group , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group, data = sampledf)
#          Df  SumOfSqs      R2      F Pr(>F)    
# group     4 0.0044641 0.75881 7.8651  0.001 ***
# Residual 10 0.0014190 0.24119                  
# Total    14 0.0058831 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



beta <- betadisper(bray, sampledf$group)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     4 4.4289e-05 1.1072e-05 0.3636    999  0.842
# Residuals 10 3.0451e-04 3.0452e-05





## Functions

phy_in <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")

phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 30125 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 30125 taxa by 4 taxonomic ranks ]

min(taxa_sums(phy_in)) # 2.521458e-06
sum(sample_sums(phy_in)) # 1500
sample_sums(phy_in)
# 20C 30B 30A UMA 10B  5A 10C 20A 20B  5B UMB 10A  5C UMC 30C 
# 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100

summary( sample_sums(phy_in) )

max(taxa_sums(phy_in)) # 9.941053




# do not rarefy as already normalized using rel abun (%)


table(phy_in@sam_data$age)
# 12 22 31  6 UM 
# 3  3  3  3  3 

phy_in@sam_data$age <- factor(phy_in@sam_data$age,
                              levels = c("6", "12", "22", "31", "UM"),
                              ordered = TRUE)




## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
#ord <- ordinate(r1.ps, "NMDS", "bray")
ord <- ordinate(phy_in, "PCoA", "bray")


ord

unique(phy_in@sam_data$age)
# [1] 12 22 31 6  UM
# Levels: 6 < 12 < 22 < 31 < UM


p <- plot_ordination(phy_in, ord, type="samples", color="age")
p

p$labels$x # "Axis.1   [49.6%]"
x_lab <- "PCo1 (49.6%)"

p$labels$y # "Axis.2   [18.8%]"
y_lab <- "PCo2 (18.8%)"

#temp <- r1.ps
p_df <- p$data


cols.group <- c("6" = "#9e0142",
                "12" = "#d53e4f",
                "22" = "#fdae61",
                "31" = "#abdda4",
                "UM" = "#3288bd")



p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = age))+ # NOTE change from group to age
  theme_bw()+
  geom_point()+
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Reveg\nage (yr)") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "Restoration: Fxns", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Functions-sunbad-resto-v2.tiff"), width = 9, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# [1] "sample_name"        "mgrast_id"          "metagenome_id"      "metagenome_name"    "investigation_type" "seq_meth"          
# [7] "file_name"          "age" 

# Adonis test
set.seed(123)
adonis2(bray ~ age , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ age, data = sampledf)
#          Df SumOfSqs      R2      F Pr(>F)   
# age       4 0.033804 0.67039 5.0846  0.002 **
# Residual 10 0.016621 0.32961                 
# Total    14 0.050425 1.00000                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



beta <- betadisper(bray, sampledf$age)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     4 0.00032954 8.2385e-05 0.3447    999  0.838
# Residuals 10 0.00238990 2.3899e-04 



#-------------------------


#### Sun & Badgely - post-mining - Lasso to model carbon sequestration rates ?
#-------------------------

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...
# $ log10_abun  : num  -3.362 -1.643 -3.341 0.711 -1.804 ...
# $ group       : Ord.factor w/ 5 levels "6"<"12"<"22"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ group_label : Ord.factor w/ 5 levels "6 yr"<"12 yr"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ ord_group   : num  3 3 3 3 3 3 3 3 3 3 ...


# later join to df.comp for compound information

data_in <- dat.cpd.collate[which(dat.cpd.collate$cpd_rel_abun > 0), c("sample", "group", "cpd_id", "cpd_rel_abun")]


unique(data_in$group)
# [1] 22 31 UM 12 6 
# Levels: 6 < 12 < 22 < 31 < UM

data_in$o_hor_rate <- NA  # Mg/ha/yr C sequestration in O Horizon
data_in$soil_c_rate <- NA # Mg/ha/yr C sequestration in mineral soil

# data from Avera et al 2015, table 3
sel <- which(data_in$group == "6") # n = 23187
data_in$o_hor_rate[sel] <- 0.91
data_in$soil_c_rate[sel] <- 0.49

sel <- which(data_in$group == "12") # n = 23236
data_in$o_hor_rate[sel] <- 0.53
data_in$soil_c_rate[sel] <- 0.70

sel <- which(data_in$group == "22") # n = 23274
data_in$o_hor_rate[sel] <- 0.43
data_in$soil_c_rate[sel] <- 0.58

sel <- which(data_in$group == "31") # n = 23236
data_in$o_hor_rate[sel] <- 0.23
data_in$soil_c_rate[sel] <- 0.33


sel <- which(data_in$group == "UM") # n = 23720
# assume these have reached steady state, and are no longer building organic matter?
data_in$o_hor_rate[sel] <- 0
data_in$soil_c_rate[sel] <- 0


head(data_in)
#   sample group   cpd_id cpd_rel_abun o_hor_rate soil_c_rate
# 1    20C    22 cpd25681 0.0004345842       0.43        0.58
# 2    20C    22 cpd02597 0.0227382574       0.43        0.58
# 3    20C    22 cpd24620 0.0004564208       0.43        0.58
# 4    20C    22 cpd00001 5.1354106008       0.43        0.58
# 5    20C    22 cpd01501 0.0157048650       0.43        0.58
# 6    20C    22 cpd00851 0.0121620567       0.43        0.58

# Convert to wide format

#df.wide <- dcast(data_in, formula = sample ~ cpd_id , value.var = "cpd_rel_abun" )
df.wide <- dcast(data_in, formula = sample + group + o_hor_rate + soil_c_rate  ~ cpd_id , value.var = "cpd_rel_abun" )

dim(df.wide) # 15 8374

head(df.wide[ ,c(1:10)])
#   sample group o_hor_rate soil_c_rate cpd00001 cpd00002 cpd00003 cpd00004  cpd00005  cpd00006
# 1    10A    12       0.53        0.70 5.126634 2.287831 1.058127 1.031673 0.7565038 0.7592450
# 2    10B    12       0.53        0.70 5.109487 2.333880 1.058379 1.032674 0.7348372 0.7370765
# 3    10C    12       0.53        0.70 5.103210 2.319224 1.070352 1.044885 0.7498325 0.7520871
# 4    20A    22       0.43        0.58 5.138510 2.323322 1.069859 1.043509 0.7560155 0.7585376
# 5    20B    22       0.43        0.58 5.173074 2.355075 1.064497 1.037668 0.7483448 0.7511583
# 6    20C    22       0.43        0.58 5.135411 2.318166 1.060200 1.033088 0.7564190 0.7590545

sel.cpd.names <- grep(pattern = "cpd", x = names(df.wide)) # n = 8370
cpd.names <- names(df.wide)[sel.cpd.names]

# are there any compounds with missing data?
# https://stackoverflow.com/questions/20364450/find-names-of-columns-which-contain-missing-values
missing1 <- names(which(colSums(is.na(df.wide)) > 0)) # n = 1259
missing2 <- colnames(df.wide)[ apply(df.wide, 2, anyNA) ]
identical(missing1, missing2) # TRUE
rm(missing1, missing2)
missing <- colnames(df.wide)[ apply(df.wide, 2, anyNA) ]
df.wide[ , missing[1]]

cpd.names <- cpd.names[-which(cpd.names %in% missing)] # n = 7111

# Lasso model C-sequestration rates using cpd_id ??

library(glmnet); packageVersion("glmnet") # ‘4.1.7’
library(epiR); packageVersion("epiR") # ‘2.0.53’


# cv.glmnet(x = as.matrix( dat[ subsel.tr[[k]] , vars.mod[vars.used] ]) ,
#           y = dat[ subsel.tr[[k]] , 1 ] , standardize=FALSE, 
#           nfolds=length(subsel.tr[[k]]),  family="gaussian")

set.seed(123)
cvFit <- cv.glmnet(x = as.matrix( df.wide[  , cpd.names ]) ,
                   y = df.wide[ , "soil_c_rate" ] , standardize=TRUE, 
                   nfolds = 3,  family="gaussian")

plot(cvFit)
coef(cvFit,s="lambda.1se")

coef.df<-as.data.frame(as.matrix( coef(cvFit,s="lambda.1se")  )) 
head(coef.df)
names(coef.df)<- "coefft"
coef.df$cpd_id<-rownames(coef.df)

sel <- which( abs(coef.df$coefft) > 0 )

coef.df[sel, ]
# coefft      cpd_id
# (Intercept)  -0.4836442 (Intercept)
# cpd00642     27.2640989    cpd00642
# cpd02993     58.7123346    cpd02993
# cpd23762    -35.9091387    cpd23762


vars.act.used <-  coef.df$cpd_id[!coef.df$coefft == 0][-1] # ignore intercept
vars.act.used
# "cpd00642" "cpd02993" "cpd23762"



x <- df.wide$soil_c_rate
y <- as.data.frame(predict(object = cvFit, newx = as.matrix(df.wide[  , cpd.names ]), type = "response", s="lambda.1se"))
y<- as.numeric(y[,1])

plot(x,y)

cor.test(x, y)
# Pearson's product-moment correlation
# 
# data:  x and y
# t = 15.091, df = 13, p-value = 1.281e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9174916 0.9910881
# sample estimates:
#   cor 
# 0.9726242 

# Lin's concordance correlation coefficient
ccc <- epi.ccc(x, y, ci = "asymptotic", conf.level = 0.95)
ccc$rho.c[,1] # 0.4769526




lasso_results <- coef.df[ abs(coef.df$coefft)>0 , ]

lasso_results$Name <- NA
lasso_results$Formula <- NA

lasso_results$Compartment <- "Mineral soil"

for (i in 1:dim(lasso_results)[1]) {
  #i<-3 
  this_cpd <- lasso_results$cpd_id[i]
  
  sel.cpd <- which(df.comp$id == this_cpd)
  
  if (length(sel.cpd) > 0) {
    lasso_results$Name[i] <- df.comp$name[sel.cpd]
    lasso_results$Formula[i] <- df.comp$form[sel.cpd]
    
  }
  print(paste0("completed ",i))
}

lasso_results
# coefft      cpd_id                                           Name        Formula  Compartment
# (Intercept)  -0.4836442 (Intercept)                                           <NA>           <NA> Mineral soil
# cpd00642     27.2640989    cpd00642                                   L-Rhamnulose        C6H12O5 Mineral soil
# cpd02993     58.7123346    cpd02993                                    Lipid IV(A) C68H126N2O23P2 Mineral soil
# cpd23762    -35.9091387    cpd23762 2,4-diacetamido-2,4,6-trideoxy-D-mannopyranose     C10H18N2O5 Mineral soil

# Prediction formula: C-sequestration (Mg/ha/yr) 
# = 27.26*[CPP:L-Rhamnulose] + 58.71*[CPP:Lipid IV(A)]
# - -35.91*[CPP:2,4-diacetamido-2,4,6-trideoxy-D-mannopyranose] -0.48

saveRDS(lasso_results, file = "lasso_results-Mineral-Soil-C-sequestration-model-Sunbad.RDS")
saveRDS(cvFit, file = "cvFit-lasso_results-Mineral-Soil-C-sequestration-model-Sunbad.RDS")
cvFit <- readRDS( "cvFit-lasso_results-Mineral-Soil-C-sequestration-model-Sunbad.RDS" )


# Plot Soil c sequestration results 
cor_test<- cor.test(x = x, y = y)
cor_test
# Pearson's product-moment correlation
# 
# data:  x and y
# t = 15.091, df = 13, p-value = 1.281e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9174916 0.9910881
# sample estimates:
#       cor 
# 0.9726242 

pval <- ifelse(test = cor_test$p.value < 0.001, yes = paste0("P < 0.001"), no = paste0("P = ",round(cor_test$p.value,3)) )
#test_result <- paste0("Pearsons correlation\nCor. = ",round(cor_test$estimate,3),";\n",pval) #
test_result <- paste0("Pearsons Cor. = ",round(cor_test$estimate,3),", ",pval,"\nLin's CCC = ",round(ccc$rho.c[,1],3),"\nTRAINING ONLY\n\nLASSO") #

p <- ggplot(data = data.frame(Observed = x, Predicted = y) , aes(x=Observed, y = Predicted))+ # log10abun
  ggtitle( "Mineral soil - C sequestration (Mg/ha/yr)" )+
  geom_point()+
  geom_point(shape = 1)+
  geom_smooth(method="lm")+
  #geom_smooth(method="loess")+
  theme_bw()+
  xlab("Observed")+ ylab("Predicted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  #scale_x_continuous(labels= names(age_vec))+
  #annotate(geom="text_npc", npcx = "middle", npcy = "top", label = this_var, size = 3 )+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 )+
  annotate(geom="text_npc", npcx = "middle", npcy = "bottom", 
           label = "Prediction formula: C-sequestration (Mg/ha/yr)\n= 27.26*[CPP:L-Rhamnulose] + 58.71*[CPP:Lipid IV(A)]\n-35.91*[CPP:2,4-diacetamido-2,4,6-trideoxy-D-mannopyranose] -0.48", size = 2, color = "darkgrey" )+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lasso-Mineral-Soil-Predicted-vs-Observed-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## model O Horizon C-sequestration

df.wide$o_hor_rate

set.seed(123)
cvFit <- cv.glmnet(x = as.matrix( df.wide[  , cpd.names ]) ,
                   y = df.wide[ , "o_hor_rate" ] , standardize=TRUE, 
                   nfolds = 3,  family="gaussian")

plot(cvFit)
coef(cvFit,s="lambda.1se")

coef.df<-as.data.frame(as.matrix( coef(cvFit,s="lambda.1se")  )) 
head(coef.df)
names(coef.df)<- "coefft"
coef.df$cpd_id<-rownames(coef.df)

sel <- which( abs(coef.df$coefft) > 0 )

coef.df[sel, ]
# coefft      cpd_id
# (Intercept)  3.779121 (Intercept)
# cpd00011    -2.752136    cpd00011
# cpd00015    -4.456938    cpd00015
# cpd00262    -4.781530    cpd00262
# cpd00982    -2.541185    cpd00982
# cpd03443    -5.770835    cpd03443


coef.df$cpd_id[!coef.df$coefft == 0]
# "(Intercept)" "cpd00011"    "cpd00015"    "cpd00262"    "cpd00982"    "cpd03443"   
vars.act.used <-  coef.df$cpd_id[!coef.df$coefft == 0][-1] # ignore intercept
vars.act.used
# "cpd00011" "cpd00015" "cpd00262" "cpd00982" "cpd03443"



x <- df.wide$o_hor_rate
y <- as.data.frame(predict(object = cvFit, newx = as.matrix(df.wide[  , cpd.names ]), type = "response", s="lambda.1se"))
y<- as.numeric(y[,1])

plot(x,y)

cor.test(x, y)
# Pearson's product-moment correlation
# 
# data:  x and y
# t = 20.217, df = 13, p-value = 3.316e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9526086 0.9949634
# sample estimates:
#       cor 
# 0.9844659

# Lin's concordance correlation coefficient
ccc <- epi.ccc(x, y, ci = "asymptotic", conf.level = 0.95)
ccc$rho.c[,1] # 0.7993679


lasso_results <- coef.df[ abs(coef.df$coefft)>0 , ]

lasso_results$Name <- NA
lasso_results$Formula <- NA

lasso_results$Compartment <- "O Horizon"

for (i in 1:dim(lasso_results)[1]) {
  #i<-3 
  this_cpd <- lasso_results$cpd_id[i]
  
  sel.cpd <- which(df.comp$id == this_cpd)
  
  if (length(sel.cpd) > 0) {
    lasso_results$Name[i] <- df.comp$name[sel.cpd]
    lasso_results$Formula[i] <- df.comp$form[sel.cpd]
    
  }
  print(paste0("completed ",i))
}

lasso_results
# coefft      cpd_id                           Name        Formula Compartment
# (Intercept)  3.779121 (Intercept)                           <NA>           <NA>   O Horizon
# cpd00011    -2.752136    cpd00011                            CO2            CO2   O Horizon
# cpd00015    -4.456938    cpd00015                            FAD  C27H31N9O15P2   O Horizon
# cpd00262    -4.781530    cpd00262                     Oxalyl-CoA C23H31N7O19P3S   O Horizon
# cpd00982    -2.541185    cpd00982                          FADH2  C27H33N9O15P2   O Horizon
# cpd03443    -5.770835    cpd03443 3-Octaprenyl-4-hydroxybenzoate       C47H69O3   O Horizon

# Prediction formula: C-sequestration (Mg/ha/yr) 
# 3.78 -2.75*[CPP:CO2] -4.46*[CPP:FAD] -4.78*[CPP:Oxalyl-CoA]
# -2.54*[CPP:FADH2] -5.77*[3-Octaprenyl-4-hydroxybenzoate]


saveRDS(lasso_results, file = "lasso_results-O-Horizon-C-sequestration-model-Sunbad.RDS")
saveRDS(cvFit, file = "cvFit-lasso_results-O-Horizon-C-sequestration-model-Sunbad.RDS")
cvFit <- readRDS("cvFit-lasso_results-O-Horizon-C-sequestration-model-Sunbad.RDS")

# Plot Soil c sequestration results 
cor_test<- cor.test(x = x, y = y)
cor_test
# Pearson's product-moment correlation
# 
# data:  x and y
# t = 20.217, df = 13, p-value = 3.316e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9526086 0.9949634
# sample estimates:
#   cor 
# 0.9844659 

pval <- ifelse(test = cor_test$p.value < 0.001, yes = paste0("P < 0.001"), no = paste0("P = ",round(cor_test$p.value,3)) )
#test_result <- paste0("Pearsons correlation\nCor. = ",round(cor_test$estimate,3),";\n",pval) #
test_result <- paste0("Pearsons Cor. = ",round(cor_test$estimate,3),", ",pval,"\nLin's CCC = ",round(ccc$rho.c[,1],3),"\nTRAINING ONLY") #


p <- ggplot(data = data.frame(Observed = x, Predicted = y) , aes(x=Observed, y = Predicted))+ # log10abun
  ggtitle( "O Horizon - C sequestration (Mg/ha/yr)" )+
  geom_point()+
  geom_point(shape = 1)+
  geom_smooth(method="lm")+
  #geom_smooth(method="loess")+
  theme_bw()+
  xlab("Observed")+ ylab("Predicted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  #scale_x_continuous(labels= names(age_vec))+
  #annotate(geom="text_npc", npcx = "middle", npcy = "top", label = this_var, size = 3 )+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 )+
  annotate(geom="text_npc", npcx = "middle", npcy = "bottom", 
           label = "Prediction formula: C-sequestration (Mg/ha/yr)\n3.78 -2.75*[CPP:CO2] -4.46*[CPP:FAD] -4.78*[CPP:Oxalyl-CoA]\n-2.54*[CPP:FADH2] -5.77*[3-Octaprenyl-4-hydroxybenzoate]", color = "darkgrey", size = 2 )+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(c)", x = unit(0.04, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lasso-O-Horizon-Predicted-vs-Observed-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## plot Reveg age vs C-sequestration

age_vec
# 6 12 22 31 UM 
# 1  2  3  4  5 

cseq <- data.frame(Age = rep(age_vec, times = 2), Type = rep(c("Mineral soil","O Horizon", "C sequestration (Mg/ha/yr)", "C pool (Mg/ha)"), each = 5), 
                   value = c(0.49,0.7,0.58,0.33,0,
                             0.91,0.53,0.43,0.23,0
                             
                             ))


p <- ggplot(data = cseq , aes(x=Age, y = value, color = Type))+
  ggtitle( "Observed C sequestration (Mg/ha/yr)" )+
  geom_point(shape = 2)+
  geom_line(aes(group = Type), linetype = "dashed")+
  theme_bw()+
  xlab("Reveg age (years)")+ ylab("C sequestration (Mg/ha/yr)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  scale_color_manual(values = c("Mineral soil" = "#a6761d" , "O Horizon" = "#7fc97f" ))+
  #annotate(geom="text_npc", npcx = "middle", npcy = "top", label = this_var, size = 3 )+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 )+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    legend.position = c(0.25,0.2),
    axis.title = element_text(size = rel(0.9))
  )
p

#grid.text(label = "(a)", x = unit(0.04, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Observed-C-sequestration-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


# Mineral soil pool and sequestration
# https://r-graph-gallery.com/line-chart-dual-Y-axis-ggplot2.html

cseq2 <- data.frame(Age = age_vec, 
                    `C sequestration\n(Mg/ha/yr)` = c(0.49,0.7,0.58,0.33,0),
                    `C pool (Mg/ha)` = c(2.45,7.67,12.24,9.76,36.32)
                    )
names(cseq2) <- c("Age","C sequestration\n(Mg/ha/yr)","C pool (Mg/ha)")
36/0.7 # 51.42857
coeff <- 50

SequestrationColor = "#a6761d"
PoolColor = "#006d2c" # "#7fc97f"

p <- ggplot(data = cseq2 , aes(x=Age))+
  ggtitle( "Soil C sequestration and C pool" )+
  geom_point(aes(y = `C sequestration\n(Mg/ha/yr)`), shape = 2, color=SequestrationColor)+
  geom_point(aes(y = `C pool (Mg/ha)` / coeff), shape = 2, color=PoolColor)+
  geom_line(aes(y = `C sequestration\n(Mg/ha/yr)`), linetype = "dashed", color=SequestrationColor)+
  geom_line(aes(y = `C pool (Mg/ha)` / coeff), linetype = "dashed", color=PoolColor)+
  theme_bw()+
  xlab("Reveg age (years)")+ #ylab("C sequestration (Mg/ha/yr)")+
  
  scale_y_continuous(
    # Features of the first axis
    name = "C sequestration (Mg/ha/yr)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="C pool (Mg/ha)")
  )+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  #scale_color_manual(values = c(`C sequestration\n(Mg/ha/yr)` = "#a6761d" , `C pool (Mg/ha)` = "#7fc97f" ))+
  
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    
    axis.title.y = element_text(color = SequestrationColor),
    axis.title.y.right = element_text(color = PoolColor),
    
    #legend.position = c(0.25,0.2),
    axis.title = element_text(size = rel(0.9))
  )
p


#grid.text(label = "(a)", x = unit(0.04, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Observed-C-sequestration-C-pool-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------


#### Sun & Badgely - post-mining - Random Forest to model carbon sequestration rates ?
#-------------------------

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...
# $ log10_abun  : num  -3.362 -1.643 -3.341 0.711 -1.804 ...
# $ group       : Ord.factor w/ 5 levels "6"<"12"<"22"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ group_label : Ord.factor w/ 5 levels "6 yr"<"12 yr"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ ord_group   : num  3 3 3 3 3 3 3 3 3 3 ...


# later join to df.comp for compound information

data_in <- dat.cpd.collate[which(dat.cpd.collate$cpd_rel_abun > 0), c("sample", "group", "cpd_id", "cpd_rel_abun")]


unique(data_in$group)
# [1] 22 31 UM 12 6 
# Levels: 6 < 12 < 22 < 31 < UM

data_in$o_hor_rate <- NA  # Mg/ha/yr C sequestration in O Horizon
data_in$soil_c_rate <- NA # Mg/ha/yr C sequestration in mineral soil

# data from Avera et al 2015, table 3
sel <- which(data_in$group == "6") # n = 23187
data_in$o_hor_rate[sel] <- 0.91
data_in$soil_c_rate[sel] <- 0.49

sel <- which(data_in$group == "12") # n = 23236
data_in$o_hor_rate[sel] <- 0.53
data_in$soil_c_rate[sel] <- 0.70

sel <- which(data_in$group == "22") # n = 23274
data_in$o_hor_rate[sel] <- 0.43
data_in$soil_c_rate[sel] <- 0.58

sel <- which(data_in$group == "31") # n = 23236
data_in$o_hor_rate[sel] <- 0.23
data_in$soil_c_rate[sel] <- 0.33


sel <- which(data_in$group == "UM") # n = 23720
# assume these have reached steady state, and are no longer building organic matter?
data_in$o_hor_rate[sel] <- 0
data_in$soil_c_rate[sel] <- 0


head(data_in)
#   sample group   cpd_id cpd_rel_abun o_hor_rate soil_c_rate
# 1    20C    22 cpd25681 0.0004345842       0.43        0.58
# 2    20C    22 cpd02597 0.0227382574       0.43        0.58
# 3    20C    22 cpd24620 0.0004564208       0.43        0.58
# 4    20C    22 cpd00001 5.1354106008       0.43        0.58
# 5    20C    22 cpd01501 0.0157048650       0.43        0.58
# 6    20C    22 cpd00851 0.0121620567       0.43        0.58

# Convert to wide format

#df.wide <- dcast(data_in, formula = sample ~ cpd_id , value.var = "cpd_rel_abun" )
df.wide <- dcast(data_in, formula = sample + group + o_hor_rate + soil_c_rate  ~ cpd_id , value.var = "cpd_rel_abun" )

dim(df.wide) # 15 8374

head(df.wide[ ,c(1:10)])
#   sample group o_hor_rate soil_c_rate cpd00001 cpd00002 cpd00003 cpd00004  cpd00005  cpd00006
# 1    10A    12       0.53        0.70 5.126634 2.287831 1.058127 1.031673 0.7565038 0.7592450
# 2    10B    12       0.53        0.70 5.109487 2.333880 1.058379 1.032674 0.7348372 0.7370765
# 3    10C    12       0.53        0.70 5.103210 2.319224 1.070352 1.044885 0.7498325 0.7520871
# 4    20A    22       0.43        0.58 5.138510 2.323322 1.069859 1.043509 0.7560155 0.7585376
# 5    20B    22       0.43        0.58 5.173074 2.355075 1.064497 1.037668 0.7483448 0.7511583
# 6    20C    22       0.43        0.58 5.135411 2.318166 1.060200 1.033088 0.7564190 0.7590545

sel.cpd.names <- grep(pattern = "cpd", x = names(df.wide)) # n = 8370
cpd.names <- names(df.wide)[sel.cpd.names]

# are there any compounds with missing data?
# https://stackoverflow.com/questions/20364450/find-names-of-columns-which-contain-missing-values
missing1 <- names(which(colSums(is.na(df.wide)) > 0)) # n = 1259
missing2 <- colnames(df.wide)[ apply(df.wide, 2, anyNA) ]
identical(missing1, missing2) # TRUE
rm(missing1, missing2)
missing <- colnames(df.wide)[ apply(df.wide, 2, anyNA) ]
df.wide[ , missing[1]]

cpd.names <- cpd.names[-which(cpd.names %in% missing)] # n = 7111

# Lasso model C-sequestration rates using cpd_id ??

library(randomForest); packageVersion("randomForest") # ‘4.7.1.1’
library(epiR); packageVersion("epiR") # ‘2.0.53’


set.seed(123)
rfFit <- randomForest(soil_c_rate ~ . , data = df.wide[  , c("soil_c_rate",cpd.names) ], ntree = 100, mtry = 5)
# Warning message:
#   In randomForest.default(m, y, ...) :
#   The response has five or fewer unique values.  Are you sure you want to do regression?


rfFit
# Call:
#   randomForest(formula = soil_c_rate ~ ., data = df.wide[, c("soil_c_rate",      cpd.names)], ntree = 100, mtry = 5) 
# Type of random forest: regression
# Number of trees: 100
# No. of variables tried at each split: 5
# 
# Mean of squared residuals: 0.03234982
# % Var explained: 44.87


rfFit$confusion

plot(rfFit)

varImpPlot(rfFit, main = "Variable importance plot - CPP predictors\nof C-sequestration in Mineral Soil")

grid.text(label = "(a)", x = unit(0.1, "npc") , y = unit(0.94,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","RF-varImpPlot-Mineral-Soil-sunbad-resto-v2.tiff"), width = 14, height = 14, units = "cm", res=450, compression="lzw",type="cairo")



x <- df.wide$soil_c_rate
y <- as.data.frame(predict(rfFit, newdata = df.wide[  , c("soil_c_rate",cpd.names) ] ))
y<- as.numeric(y[,1])

plot(x,y)

cor.test(x, y)
# data:  x and y
# t = 17.548, df = 13, p-value = 1.964e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9378905 0.9933545
# sample estimates:
#   cor 
# 0.9795379

# Lin's concordance correlation coefficient
ccc <- epi.ccc(x, y, ci = "asymptotic", conf.level = 0.95)
ccc$rho.c[,1] # 0.9290218









# Plot Soil c sequestration results 
cor_test<- cor.test(x = x, y = y)
cor_test
# Pearson's product-moment correlation
# 
# data:  x and y
# t = 17.548, df = 13, p-value = 1.964e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9378905 0.9933545
# sample estimates:
#   cor 
# 0.9795379 

pval <- ifelse(test = cor_test$p.value < 0.001, yes = paste0("P < 0.001"), no = paste0("P = ",round(cor_test$p.value,3)) )
#test_result <- paste0("Pearsons correlation\nCor. = ",round(cor_test$estimate,3),";\n",pval) #
test_result <- paste0("Pearsons Cor. = ",round(cor_test$estimate,3),", ",pval,"\nLin's CCC = ",round(ccc$rho.c[,1],3),"\nTRAINING ONLY\n\nRANDOM FOREST") #

p <- ggplot(data = data.frame(Observed = x, Predicted = y) , aes(x=Observed, y = Predicted))+ # log10abun
  ggtitle( "Mineral soil - C sequestration (Mg/ha/yr)" )+
  geom_point()+
  geom_point(shape = 1)+
  geom_smooth(method="lm")+
  #geom_smooth(method="loess")+
  theme_bw()+
  xlab("Observed")+ ylab("Predicted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  #scale_x_continuous(labels= names(age_vec))+
  #annotate(geom="text_npc", npcx = "middle", npcy = "top", label = this_var, size = 3 )+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 )+
  # annotate(geom="text_npc", npcx = "middle", npcy = "bottom", 
  #          label = "Prediction formula: C-sequestration (Mg/ha/yr)\n= 27.26*[CPP:L-Rhamnulose] + 58.71*[CPP:Lipid IV(A)]\n-35.91*[CPP:2,4-diacetamido-2,4,6-trideoxy-D-mannopyranose] -0.48", size = 2, color = "darkgrey" )+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","RF-Mineral-Soil-Predicted-vs-Observed-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## model O Horizon C-sequestration

df.wide$o_hor_rate

set.seed(123)
cvFit <- cv.glmnet(x = as.matrix( df.wide[  , cpd.names ]) ,
                   y = df.wide[ , "o_hor_rate" ] , standardize=TRUE, 
                   nfolds = 3,  family="gaussian")

plot(cvFit)
coef(cvFit,s="lambda.1se")

coef.df<-as.data.frame(as.matrix( coef(cvFit,s="lambda.1se")  )) 
head(coef.df)
names(coef.df)<- "coefft"
coef.df$cpd_id<-rownames(coef.df)

sel <- which( abs(coef.df$coefft) > 0 )

coef.df[sel, ]
# coefft      cpd_id
# (Intercept)  3.779121 (Intercept)
# cpd00011    -2.752136    cpd00011
# cpd00015    -4.456938    cpd00015
# cpd00262    -4.781530    cpd00262
# cpd00982    -2.541185    cpd00982
# cpd03443    -5.770835    cpd03443


coef.df$cpd_id[!coef.df$coefft == 0]
# "(Intercept)" "cpd00011"    "cpd00015"    "cpd00262"    "cpd00982"    "cpd03443"   
vars.act.used <-  coef.df$cpd_id[!coef.df$coefft == 0][-1] # ignore intercept
vars.act.used
# "cpd00011" "cpd00015" "cpd00262" "cpd00982" "cpd03443"



x <- df.wide$o_hor_rate
y <- as.data.frame(predict(object = cvFit, newx = as.matrix(df.wide[  , cpd.names ]), type = "response", s="lambda.1se"))
y<- as.numeric(y[,1])

plot(x,y)

cor.test(x, y)
# Pearson's product-moment correlation
# 
# data:  x and y
# t = 20.217, df = 13, p-value = 3.316e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9526086 0.9949634
# sample estimates:
#       cor 
# 0.9844659

# Lin's concordance correlation coefficient
ccc <- epi.ccc(x, y, ci = "asymptotic", conf.level = 0.95)
ccc$rho.c[,1] # 0.7993679


lasso_results <- coef.df[ abs(coef.df$coefft)>0 , ]

lasso_results$Name <- NA
lasso_results$Formula <- NA

lasso_results$Compartment <- "O Horizon"

for (i in 1:dim(lasso_results)[1]) {
  #i<-3 
  this_cpd <- lasso_results$cpd_id[i]
  
  sel.cpd <- which(df.comp$id == this_cpd)
  
  if (length(sel.cpd) > 0) {
    lasso_results$Name[i] <- df.comp$name[sel.cpd]
    lasso_results$Formula[i] <- df.comp$form[sel.cpd]
    
  }
  print(paste0("completed ",i))
}

lasso_results
# coefft      cpd_id                           Name        Formula Compartment
# (Intercept)  3.779121 (Intercept)                           <NA>           <NA>   O Horizon
# cpd00011    -2.752136    cpd00011                            CO2            CO2   O Horizon
# cpd00015    -4.456938    cpd00015                            FAD  C27H31N9O15P2   O Horizon
# cpd00262    -4.781530    cpd00262                     Oxalyl-CoA C23H31N7O19P3S   O Horizon
# cpd00982    -2.541185    cpd00982                          FADH2  C27H33N9O15P2   O Horizon
# cpd03443    -5.770835    cpd03443 3-Octaprenyl-4-hydroxybenzoate       C47H69O3   O Horizon

# Prediction formula: C-sequestration (Mg/ha/yr) 
# 3.78 -2.75*[CPP:CO2] -4.46*[CPP:FAD] -4.78*[CPP:Oxalyl-CoA]
# -2.54*[CPP:FADH2] -5.77*[3-Octaprenyl-4-hydroxybenzoate]


saveRDS(lasso_results, file = "lasso_results-O-Horizon-C-sequestration-model-Sunbad.RDS")
saveRDS(cvFit, file = "cvFit-lasso_results-O-Horizon-C-sequestration-model-Sunbad.RDS")
cvFit <- readRDS("cvFit-lasso_results-O-Horizon-C-sequestration-model-Sunbad.RDS")

# Plot Soil c sequestration results 
cor_test<- cor.test(x = x, y = y)
cor_test
# Pearson's product-moment correlation
# 
# data:  x and y
# t = 20.217, df = 13, p-value = 3.316e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9526086 0.9949634
# sample estimates:
#   cor 
# 0.9844659 

pval <- ifelse(test = cor_test$p.value < 0.001, yes = paste0("P < 0.001"), no = paste0("P = ",round(cor_test$p.value,3)) )
#test_result <- paste0("Pearsons correlation\nCor. = ",round(cor_test$estimate,3),";\n",pval) #
test_result <- paste0("Pearsons Cor. = ",round(cor_test$estimate,3),", ",pval,"\nLin's CCC = ",round(ccc$rho.c[,1],3),"\nTRAINING ONLY") #


p <- ggplot(data = data.frame(Observed = x, Predicted = y) , aes(x=Observed, y = Predicted))+ # log10abun
  ggtitle( "O Horizon - C sequestration (Mg/ha/yr)" )+
  geom_point()+
  geom_point(shape = 1)+
  geom_smooth(method="lm")+
  #geom_smooth(method="loess")+
  theme_bw()+
  xlab("Observed")+ ylab("Predicted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  #scale_x_continuous(labels= names(age_vec))+
  #annotate(geom="text_npc", npcx = "middle", npcy = "top", label = this_var, size = 3 )+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 )+
  annotate(geom="text_npc", npcx = "middle", npcy = "bottom", 
           label = "Prediction formula: C-sequestration (Mg/ha/yr)\n3.78 -2.75*[CPP:CO2] -4.46*[CPP:FAD] -4.78*[CPP:Oxalyl-CoA]\n-2.54*[CPP:FADH2] -5.77*[3-Octaprenyl-4-hydroxybenzoate]", color = "darkgrey", size = 2 )+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(c)", x = unit(0.04, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lasso-O-Horizon-Predicted-vs-Observed-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## plot Reveg age vs C-sequestration

age_vec
# 6 12 22 31 UM 
# 1  2  3  4  5 

cseq <- data.frame(Age = rep(age_vec, times = 2), Type = rep(c("Mineral soil","O Horizon"), each = 5), 
                   value = c(0.49,0.7,0.58,0.33,0,
                             0.91,0.53,0.43,0.23,0))


p <- ggplot(data = cseq , aes(x=Age, y = value, color = Type))+
  ggtitle( "Observed C sequestration (Mg/ha/yr)" )+
  geom_point(shape = 2)+
  geom_line(aes(group = Type), linetype = "dashed")+
  theme_bw()+
  xlab("Reveg age (years)")+ ylab("C sequestration (Mg/ha/yr)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  #annotate(geom="text_npc", npcx = "middle", npcy = "top", label = this_var, size = 3 )+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 )+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    legend.position = c(0.25,0.2),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(a)", x = unit(0.04, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Observed-C-sequestration-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------


##########################
##########################
##########################
##########################


#### Forslund-SWE-T2D - build reaction search in parallel - get_reactions & compounds
#### Individual compound level !!
#-------------------------

phy <- readRDS("phy-phyloseq-object-Forslund-SWE-T2D.RDS")

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 19099    4


get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
  
  rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
  rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
  
  ###testing
  # #for (i in 1:dim(df.tax)[1]) {
  # for (i in 1:1000) {
  #i<-1
  #i<-22889
  #i<-22894
  #i<-31768
  #i<-9406
  #i<-9422
  #i<-8
  #1<-13
  #i<-7
  
  sub1 <- df.tax$subsys_L1[i]
  sub2 <- df.tax$subsys_L2[i]
  sub3 <- df.tax$subsys_L3[i]
  
  fxn.temp <- df.tax$fxn[i]
  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  
  # store results corresponding to each Superfocus row
  fxn.list <- list()
  fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
  
  # check for multiple functions/reactions?
  flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
  flag2 <- grepl(pattern = "_@_", x = fxn.temp)
  if (!any(flag1,flag2)==TRUE) {
    # no multiples
    fxns <- fxn.temp
  } else if (flag1==TRUE) {
    fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
  } else {
    fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
  }
  # remove underscores
  ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
  
  # process each fxn & store attributes
  #df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, min_adist_modelSEED=NA, min_amatch_modelSEED=NA, rxns=NA, tot_mean_OC_x=NA, tot_mean_HC_y=NA , tot_mean_NC_z=NA )
  
  df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA) #, compound_name=NA, compound_form=NA, compound_molar_prop=NA, OC_x=NA, HC_y=NA , NC_z=NA )
  
  # # do round brackets interfere with search? - YES
  # lookfor <- "option 4 (this one)"
  # lookuplist <- c("option 1", "option 2", "option 3 (this one)", "option 4 (this one)")
  # grep(pattern = lookfor, x = lookuplist)
  
  # Identify '/' separators with no '_'  ??
  
  for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
    #f<-1
    #f<-2
    f.in <- fxns[f]
    
    # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
    full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
    full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
    
    ## data cleaning
    
    # trim off '_#' and '_##' tags
    trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
    if (!is.na(trim_nchar) & length(trim_nchar)==1) {
      f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
    }
    
    # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
    f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
    
    #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
    #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
    
    full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
    full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
    
    sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
    
    ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
    if (length(sel.rx)>=1) {
      df.fxns$matching_method[f] <- "Exact hierachy match"
      df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
      
    } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
      # search by EC id if present
      
      f.in <- fxns[f] # this goes back to string with brackets for EC
      ## LOOK FOR MULTIPLE ECs ??????????
      # 22889
      # 22894
      # 31768
      
      how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
      
      ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
      #class(ECs)
      ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
      ECs.collapse <- paste0(ECs, collapse = "|")
      
      sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
      
      if (length(how_many_ECs)==0 | length(ECs)==0) {
        # there was a glitch, database typo, or some error in identifying the EC number
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      } else if (length(sel.rx)>=1) {
        # combined EC hits identified
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
        # treat EC hits individually
        sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
        # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
        sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else {
        # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      # END EC matching
      
      
    } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
      ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
      # 1. 'reactions' table by name: rxns.lut$name
      # 2. 'reactions' table by aliases: rxns.lut$aliases
      # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
      # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
      
      if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
        # 1a - exact match - rxns.lut$name
        sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
        #rxns.lut$name[sel.rx]
        df.fxns$matching_method[f] <- "Matched Reactions name"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
        # 2a - exact match - rxns.lut$aliases
        sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
        #rxns.lut$aliases[sel.rx]
        #rxns.lut$name[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Reactions aliases"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
        # 3a - exact match - subsys.lut$Role
        sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
        #subsys.lut$Role[sel.rx]
        #subsys.lut$Reaction[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Subsytem role"
        df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
        # 4a - exact match - rxn_pathways.lut$External_rxn_name
        sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
        
        df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
        df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
        
        
      } else {
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      
      ## DON'T RUN PARTIAL MATCHING AT THIS STAGE
      
      
    } # END function - reaction search
    
    #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
    
    print(paste0("completed fxn ", f))
    
    
    ## now investigate these reactions ...
    # Reactions lookup table: 
    # - "equation": Definition of reaction expressed using compound IDs and after protonation
    # Compounds lookup table:
    # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
    #df.fxns
    
    
    #if (df.fxns$matching_method == "No match found") {
    if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
      
      df.Rxns <- NA
      df.Compounds <- NA
      
    } else { # reaction(s) were identified
      
      # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
      f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
      
      df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
                            rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA ) #, OC_ratios=NA, HC_ratios=NA, NC_ratios=NA, coefwtmean_OC_x=NA, coefwtmean_HC_y=NA, coefwtmean_NC_z=NA)
      
      #df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns, compound_name=NA, compound_form=NA, compound_molar_prop=NA, OC_x=NA, HC_y=NA , NC_z=NA )
      
      for (r in 1:dim(df.Rxns)[1]) {
        #r<-1
        #this_rxn <- "rxn00004"
        this_rxn <- df.Rxns$rxn_id[r]
        sel <- which(rxns.lut$id == this_rxn)
        ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
        ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
        ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
        
        # extract compound info
        
        #df.Rxns$rxn_eqn[r]
        #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
        #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
        
        ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
        # 5 23 43 61
        # 6 25 46 65 83
        
        ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
        # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
        
        if (length(compds)>=1) {
          
          df.Rxns$compds[r] <- paste0(compds, collapse = ";")
          
          ## get compound coefficients?
          start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
          end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
          ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
          
          df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
          
          # get formulas of compounds
          
          formx <-filter(compounds.lut, id %in% compds )
          row.names(formx) <- formx$id
          ( formx.char <- formx[compds, ]$formula )
          # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
          # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
          # "C7H7O7" "H2O"    "C7H5O6"
          df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
          
          ( compd.names <- formx[compds, ]$name )
          # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
          
          
          # # O:C ratio (replace Inf with NA when C not present)
          # OC_ratio <- atomic_no(form = formx.char, atom = "O")/atomic_no(form = formx.char, atom = "C")
          # OC_ratio[is.infinite(OC_ratio)] <- NA
          # 
          # # H:C ratio (replace Inf with NA when C not present)
          # HC_ratio <- atomic_no(form = formx.char, atom = "H")/atomic_no(form = formx.char, atom = "C")
          # HC_ratio[is.infinite(HC_ratio)] <- NA
          # 
          # # N:C ratio (replace Inf with NA when C not present)
          # NC_ratio <- atomic_no(form = formx.char, atom = "N")/atomic_no(form = formx.char, atom = "C")
          # NC_ratio[is.infinite(NC_ratio)] <- NA
          
          temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
                                          cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
                                          #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
          )
          
        } else {
          # No specified reaction equation or chemical formula info
          df.Rxns$compds[r] <- NA
          df.Rxns$compd_coef[r] <- NA
          df.Rxns$chem_formx[r] <- NA
          
          temp.df.Compounds <- NA
          
        }
        
        if (r==1) { df.Compounds <- temp.df.Compounds }
        
        if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
        
        # clean up - if there are additional reactions?
        temp.df.Compounds <- NA
        
      } # END loop for r - rxn_id's per f/f.in
      
    } # END else loop when reactions identified
    
    # store results corresponding to each sub-reaction of each Superfocus row
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
    
    
  } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
  
  
  #return(fxn.list)
  
  saveRDS(object = fxn.list, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv/fxn-list-",fxn.superfocus.rowlabel,".rds") ) # use readRDS()
  
  #print(paste0("COMPLETED ROW ",i," OF SUPERFOCUS FUNCTIONAL TAXA  # # # # # # # # # # # # # # # # # # # # #"))
  
} # END function to be run in parallel for each superfocus row

# HERE

# # # # # # # # # # # # # # # # # #


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

#foreach(i=1:100 , .packages=c('stringr', 'dplyr')) %dopar%
foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
  get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )

stopCluster(cl)
time.finish <- Sys.time()



time.start #  
#[1] "2024-12-17 20:12:22 ACDT"
time.finish # 
#[1] "2024-12-17 22:46:46 ACDT"


## assemble results


modelSEED_rxn_result_dir <- "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv"


dim(df.tax)
# 19099     4


# read first output
i<-1
#temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-fxn_",i,".rds"))
fxn.superfocus.rowlabel <- row.names(df.tax)[i]
temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))

print( length(temp) )
print( names(temp) )
# "fxn_1"
print( paste0( class(temp[[1]])," ",names(temp)," of length ", length(temp[[1]])," named ",paste0(names(temp[[1]]), collapse = " & ") ))
# "list fxn_1 of length 3 named fxns & rxns & compounds"
class ( temp[[1]][["fxns"]] ) # data.frame
class ( temp[[1]][["rxns"]] ) # list
class ( temp[[1]][["compounds"]] ) # list

length( temp[[1]][["compounds"]] ) # 1
names( temp[[1]][["compounds"]] ) # NULL
class( temp[[1]][["compounds"]][[1]] ) # "logical"
is.na( temp[[1]][["compounds"]][[1]] )

i<-2
fxn.superfocus.rowlabel <- row.names(df.tax)[i]
temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))

length(temp) # 1
names(temp) # "fxn_2"
print( paste0( class(temp[[1]])," ",names(temp)," of length ", length(temp[[1]])," named ",paste0(names(temp[[1]]), collapse = " & ") ))
# "list fxn_2 of length 3 named fxns & rxns & compounds"
class ( temp[[1]][["fxns"]] ) # data.frame
class ( temp[[1]][["rxns"]] ) # list
class ( temp[[1]][["compounds"]] ) # list

length( temp[[1]][["compounds"]] ) # 1
names( temp[[1]][["compounds"]] ) # NULL
class( temp[[1]][["compounds"]][[1]] ) # "data.frame"

df.out <- temp[[1]][["compounds"]][[1]]

names(df.out) #
# [1] "superfocus_fxn" "f"              "f__in"          "rxn_id"         "cpd_id"         "cpd_name"       "cpd_form"       "cpd_molar_prop" 



# total results written to disk - 30125
dim(df.tax) # 19099     4
num_results_files <- dim(df.tax)[1]
#num_results_files <- 100


# assemble all compound data outputs
# start with blank row

df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA #, 
                     #OC_x=NA, HC_y=NA, NC_z=NA
)

for (i in 1:num_results_files) {
  #i<-1
  #i<-2
  #i<-13
  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
  
  f_no <- length( temp[[1]][["compounds"]] )
  
  for (f in 1:f_no) {
    #f<-2
    # only add non-NA results
    if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
      
      df.temp <- temp[[1]][["compounds"]][[f]]
      ok <- complete.cases(df.temp)
      df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
      df.out <- rbind(df.out,df.temp)
    }
  }
  
  print(paste0("added df ",i," of ",num_results_files ))
  
}


str(df.out)
# 'data.frame':	545807 obs. of  8 variables:

saveRDS(object = df.out, file = "df.out--get_rxns_and_compounds_indiv--Forslund-SWE-T2D.RDS")
df.out <- readRDS(file = "df.out--get_rxns_and_compounds_indiv--Forslund-SWE-T2D.RDS")

# remove NA first row
head(df.out)
# superfocus_fxn  f                                                                         f__in   rxn_id   cpd_id
# 1           <NA> NA                                                                          <NA>     <NA>     <NA>
#   2          fxn_2  1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd24620
# 3          fxn_2  1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd00001
# 4          fxn_2  1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd25681
# 5          fxn_3  1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd01501
# 6          fxn_3  1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd00001
# cpd_name cpd_form cpd_molar_prop
# 1                                         <NA>     <NA>             NA
# 2 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7              1
# 3                                          H2O      H2O              1
# 4                     2-methyl-trans-aconitate   C7H5O6              1
# 5                              2-Methylcitrate   C7H7O7              1
# 6                                          H2O      H2O              1

df.out <- df.out[-1, ]


# check for different cpd_molar_prop ??
hist(df.out$cpd_molar_prop)

dim(df.out) # 545806      8


# normalise molar_prop to cpd_relabun so total of 1 per superfocus function !!

df.out$cpd_molar_prop_norm <- NA

length(unique(df.out$superfocus_fxn)) # 10576

phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 19099 taxa and 145 samples ]
# sample_data() Sample Data:       [ 145 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 19099 taxa by 4 taxonomic ranks ]

100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy)) # 55.37463 % of functions represented
100*(10576/19099) # 55.37463

fxns_found <- unique(df.out$superfocus_fxn)

for (k in 1:length(fxns_found)) {
  #k<-1
  this_fxn <- fxns_found[k]
  sel <- which(df.out$superfocus_fxn == this_fxn)
  
  sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
  # calculate 
  
  df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
  
  print(paste0("completed ",k))
  
}

sum(df.out$cpd_molar_prop_norm) # 10576


sample_sums(phy)
# all 100

dim(df.out) # 545806      9




getwd() # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"

saveRDS(object = df.out, file = "df.out--tidy-compounds_indiv--cpp3d-Forslund-SWE-T2D.RDS")

df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-Forslund-SWE-T2D.RDS" )


#-------------------------


#### Forslund-SWE-T2D - get cpd rel abun per sample
#### Individual compound level !!
#-------------------------

this_study <- "-Forslund-SWE-T2D-"
phy <- readRDS("phy-phyloseq-object-Forslund-SWE-T2D.RDS")

df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-Forslund-SWE-T2D.RDS" )
dim(df.out) # 545806      9

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 19099     4

phy@sam_data


samples.char.T2D <- sample_names(phy)
saveRDS(object = samples.char.T2D , file = "samples.char.T2D.RDS")
samples.char <- readRDS("samples.char.T2D.RDS")

df.OTU <- as.data.frame( phy@otu_table )
saveRDS(object = df.OTU , file = "df.OTU.T2D.RDS")
df.OTU <- readRDS("df.OTU.T2D.RDS")


table(phy@sam_data$Status)
# ND CTRL T2D metformin- T2D metformin+ 
#   92             33             20

sample_names(phy)
#  [1] "ERR260132" "ERR260133" "ERR260134" "ERR260135"  etc
nsamples(phy) # 145


# SRA Runs Metadata: 
# "Metadata-Forslund-ERP002469-SWE-samples--SraRunTable"  

meta <- read_excel(path= "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/data/Forslund-2015-Type2-diabetes/Metadata-Forslund-ERP002469-SWE-samples--SraRunTable.xlsx",
                   sheet=1, range="A1:AI146", col_names = TRUE)
meta <- as.data.frame(meta)
str(meta)
# 'data.frame':	145 obs. of  35 variables:

identical(sample_names(phy), meta$Run) # TRUE


df.samp <- as.data.frame(phy@sam_data)

identical(df.samp$Run, meta$Run) # TRUE

names(meta)
# [1] "Run"                    "Assay Type"             "AvgSpotLen"             "Bases"                 
# [5] "BioProject"             "BioSample"              "Bytes"                  "Center Name"           
# [9] "Consent"                "DATASTORE filetype"     "DATASTORE provider"     "DATASTORE region"      
# [13] "ENA-FIRST-PUBLIC (run)" "ENA-FIRST-PUBLIC"       "ENA-LAST-UPDATE (run)"  "ENA-LAST-UPDATE"       
# [17] "Experiment"             "External_Id"            "INSDC_center_alias"     "INSDC_center_name"     
# [21] "INSDC_first_public"     "INSDC_last_update"      "INSDC_status"           "Instrument"            
# [25] "Library Name"           "LibraryLayout"          "LibrarySelection"       "LibrarySource"         
# [29] "Organism"               "Platform"               "ReleaseDate"            "Sample Name"           
# [33] "Sample_Name"            "SRA Study"              "Submitter_Id"


# but additional data from Karlsson et al 2013 includes impaired glucose tolerance (IGT; n = 49) or normal glucose tolerance (NGT; n = 43)

# this lookup has IGT / NGT classification
t2d.class <- read_excel(path= "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/data/Forslund-2015-Type2-diabetes/Karlsson et al 2013--41586_2013_BFnature12198_MOESM507_ESM.xlsx",
                        sheet="Supplementary Table 3", range="A2:AD147", col_names = TRUE)
t2d.class <- as.data.frame(t2d.class)
str(t2d.class)
# 'data.frame':	145 obs. of  30 variables:
# $ Sample ID                                                                             : num  51 53 54 58 59 60 77 80 88 92 ...
# $ Age (years)                                                                           : num  69.1 70.3 69.9 70.2 69.4 ...
# $ Classification                                                                        : chr  "IGT" "NGT" "IGT" "NGT" ...
# etc ...

# this lookup has no of reads / bases
t2d.nbases <- read_excel(path= "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/data/Forslund-2015-Type2-diabetes/Karlsson et al 2013--41586_2013_BFnature12198_MOESM507_ESM.xlsx",
                         sheet="Supplementary Table 4", range="A2:C147", col_names = TRUE)
t2d.nbases <- as.data.frame(t2d.nbases)
str(t2d.nbases)
# 'data.frame':	145 obs. of  3 variables:
# $ Sample ID                 : num  51 53 54 58 59 60 77 80 88 92 ...
# $ Number of reads           : num  10653221 13971868 13561331 13213300 12460015 ...
# $ Total number of bases (bp): num  2.15e+09 2.82e+09 2.74e+09 2.67e+09 2.52e+09 ...

identical(t2d.class$`Sample ID`,t2d.nbases$`Sample ID`) # TRUE


length(unique(meta$Bases)) # 145

unique_samps <- unique(df.samp$Run) # qty 145

temp <- df.samp

df.samp$group_new <- NA
df.samp$age <- NA

for (i in 1:dim(df.samp)[1]) {
  #i<-1
  this_run <- unique_samps[i]
  sel.bases.row <- which(t2d.nbases$`Total number of bases (bp)` == df.samp$Bases[i])
  
  df.samp$group_new[i] <- paste0(df.samp$Status[i],"__", t2d.class$Classification[sel.bases.row])
  
  df.samp$age[i] <- t2d.class$`Age (years)`[sel.bases.row]
  
  print(paste0("completed ", i))
  
}


unique(df.samp$group_new) # "ND CTRL__IGT"        "T2D metformin-__T2D" "ND CTRL__NGT"        "T2D metformin+__T2D"


# T2D (n = 53), impaired glucose tolerance (IGT; n = 49) or normal glucose tolerance (NGT; n = 43)

df.samp$group_new <- gsub(pattern = "ND CTRL__NGT", replacement = "Normal", x = df.samp$group_new)
df.samp$group_new <- gsub(pattern = "ND CTRL__IGT", replacement = "IGT", x = df.samp$group_new)
df.samp$group_new <- gsub(pattern = "T2D metformin-__T2D", replacement = "T2D met neg", x = df.samp$group_new)
df.samp$group_new <- gsub(pattern = "T2D metformin+__T2D", replacement = "T2D met pos", x = df.samp$group_new)
# did not work?
sel <- which(df.samp$group_new == "T2D metformin+__T2D")
df.samp$group_new[sel] <- "T2D met pos"

unique(df.samp$group_new)
# "IGT"         "T2D met neg" "Normal"      "T2D met pos"

df.samp$group_new <- factor(df.samp$group_new, levels = c("T2D met neg", "T2D met pos", "IGT", "Normal"), ordered=TRUE)

identical( phy@sam_data$Run , df.samp$Run ) # TRUE


saveRDS(object = df.samp, file = "df.samp.with-t2dclass-age-Forslund-SWE-T2D.RDS")


sample_names(phy)
identical( sample_names(phy), colnames( as.matrix( phy@otu_table)) ) # TRUE

df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
dim(df.OTU) # 19099   145
df.OTU[1:5, 1:8]
# ERR260132   ERR260133   ERR260134  ERR260135    ERR260136   ERR260137   ERR260138   ERR260139
# fxn_1 0.00000000 0.000000000 0.000000000 0.00000000 0.0000000000 0.000000000 0.000000000 0.000000000
# fxn_2 0.00000000 0.000000000 0.000000000 0.00000000 0.0000000000 0.000000000 0.000000000 0.000000000
# fxn_3 0.00000000 0.000000000 0.000000000 0.00000000 0.0000000000 0.000000000 0.000000000 0.000000000
# fxn_4 0.00000000 0.000000000 0.000000000 0.00000000 0.0000000000 0.000000000 0.000000000 0.000000000
# fxn_5 0.01602096 0.004141166 0.004530646 0.00113941 0.0007990692 0.001384355 0.001487265 0.009191508

sample_sums(phy) # all values of 100



## create df.cpd_vk_coords_long

# loop through each sample

# add grouping variables

# for each function, assign relative abundance across selected compounds

## later collapse to unique combos of "OC_x__HC_y__NC_z", add rel_abun, collate unique cpd_form, collate rxn_id, collate cpd_id 


get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
  #i<-1
  #phy_in = phy
  #dat.cpd = df.out
  
  this_samp <- sample_names(phy_in)[i]
  df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
  
  dat.cpd$sample <- this_samp
  
  dat.cpd$cpd_rel_abun_norm <- NA
  
  fxns_all <- row.names(df.OTU)
  
  for (k in 1:length(fxns_all)) {
    #k<-1
    this_fxn <- fxns_all[k]
    sel <- which(dat.cpd$superfocus_fxn == this_fxn)
    
    if (length(sel)>=1) {
      dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
      
    }
  } # END rel abun values for all relevant functions added
  
  saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  
} # END


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)

stopCluster(cl)
time.finish <- Sys.time()

# output 1
i<-1
this_samp <- sample_names(phy)[i]
#saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
dat <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
head(dat)

for ( i in 2:length(sample_names(phy)) ) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  #saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  temp <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  
  dat <- rbind(dat, temp)
  
  print(paste0("completed ",i))
}


saveRDS(object = dat, file = "dat.cpd-long-all-samps-cpp3d-Forslund-SWE-T2D.rds" )
dat <- readRDS("dat.cpd-long-all-samps-cpp3d-Forslund-SWE-T2D.rds")

rm(temp)

str(dat)
# 'data.frame':	79141870 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_2" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn25279" "rxn25279" "rxn25279" "rxn03060" ...
# $ cpd_id             : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ cpd_name           : chr  "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "2-methyl-trans-aconitate" "2-Methylcitrate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.3333 0.3333 0.3333 0.0556 0.0556 ...
# $ sample             : chr  "ERR260132" "ERR260132" "ERR260132" "ERR260132" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...

#dat$combos <- paste0(dat$OC_x,"__",dat$HC_y,"__",dat$NC_z)

sum(dat$cpd_rel_abun_norm) # 10255.36
sum(dat$cpd_rel_abun_norm)/nsamples(phy) # 70.72665 = average 70.7% functional relative abundance per sample

length(which(is.na(dat$cpd_rel_abun_norm))==TRUE) # 0
length(which( dat$cpd_rel_abun_norm > 0) == TRUE) #  27263110
length(which( dat$cpd_rel_abun_norm == 0) == TRUE) # 51878760

# so this step does collect some zero relative abundances from the 'otu-table' 

names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"               "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm" "sample"              "cpd_rel_abun_norm" 

## create dat.cpd.distil
## later collapse to unique combos of "OC_x__HC_y__NC_z", add rel_abun, collate unique cpd_form, collate rxn_id, collate cpd_id 

length(unique(dat$cpd_id)) # 7261


## Collate compounds within each sample 


unique_cpd <- unique(dat$cpd_id)
samp_names <- sample_names(phy)



collate_compounds <- function(dat.cpd, unique_cpd, samp) {
  #i<-1
  #samp = samp_names[i]
  #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
  
  this_samp <- samp
  
  cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, #OC_x=NA, HC_y=NA, NC_z=NA, 
                         cpd_rel_abun=NA)
  
  for (c in 1:length(unique_cpd)) {
    #c<-1
    this_cpd <- unique_cpd[c]
    sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
    
    # cpd_data$OC_x[c] <- unique(dat.cpd$OC_x[sel.cpd])
    # cpd_data$HC_y[c] <- unique(dat.cpd$HC_y[sel.cpd])
    # cpd_data$NC_z[c] <- unique(dat.cpd$NC_z[sel.cpd])
    
    # # now collate rel abun for this sample only
    # sel <- which(dat.cpd$sample == this_samp)
    # dat.cpd <- dat.cpd[sel, ]
    
    # # again select in this sample-specific dataset
    # sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
    
    if (length(sel.cpd) >=1) {
      cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
    }
    
  } # END all compounds
  
  saveRDS(object = cpd_data, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
  
} # END


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])

stopCluster(cl)
time.finish <- Sys.time()




# output 1
i<-1
this_samp <- sample_names(phy)[i]
dat.cpd.collate <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
head(dat.cpd.collate)

for ( i in 2:length(sample_names(phy)) ) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  temp <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Forslund-SWE-T2D-R-working-files-indiv/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
  
  dat.cpd.collate <- rbind(dat.cpd.collate, temp)
  
  print(paste0("completed ",i))
}


str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  3 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260132" "ERR260132" "ERR260132" "ERR260132" ...
# $ cpd_rel_abun: num  0 5.647932 0 0 0.000186 ...

sum(dat.cpd.collate$cpd_rel_abun) # 10255.36
sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample)) # 70.72665

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-Forslund-SWE-T2D.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-Forslund-SWE-T2D.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000134 0.009741 0.001372 7.397862

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf -3.8740    -Inf -2.8626  0.8691


# log10 abun
dat.cpd.collate$log10_abun <- dat.cpd.collate$cpd_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(dat.cpd.collate$log10_abun == 0) # qty 284945
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(dat.cpd.collate$log10_abun[ -subsel.zero ])
  dat.cpd.collate$log10_abun[ subsel.zero ] <- zero_replace
}
dat.cpd.collate$log10_abun <- log10(dat.cpd.collate$log10_abun)

hist(dat.cpd.collate$log10_abun); summary( dat.cpd.collate$log10_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.4820 -8.4820 -3.8740 -4.8087 -2.8626  0.8691


# dat.sel <- dat.cpd.collate[ which(dat.cpd.collate$cpd_rel_abun > 0), ]
# n <- length(unique(dat.sel$sample))
# sum(dat.sel$abun) # 982.155
# sum(dat.sel$abun)/n # 65.477
# dat.sel$abun <- dat.sel$abun/n


# make group variable from sample name

dat.cpd.collate$group <- NA


# from above

identical( phy@sam_data$Run , df.samp$Run ) # TRUE
identical( sample_names(phy), df.samp$Run ) # TRUE
unique(df.samp$group_new)
# [1] IGT         T2D met neg Normal      T2D met pos
# Levels: T2D met neg < T2D met pos < IGT < Normal

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( df.samp$Run )) {
  #i<-1
  #this_samp <- sample_names(phy)[i]
  this_samp <- df.samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  #dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  dat.cpd.collate$group[sel] <- as.character( df.samp$group_new[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group) # "IGT"         "T2D met neg" "Normal"      "T2D met pos"
dat.cpd.collate$group <- factor(dat.cpd.collate$group, levels = c("Normal", "IGT", "T2D met pos", "T2D met neg"), ordered = TRUE)

dat.cpd.collate$group_label <- factor(dat.cpd.collate$group, 
                                      levels = c("Normal", "IGT", "T2D met pos", "T2D met neg"),
                                      labels = c("Normal", "IGT", "T2D met+", "T2D met-"),ordered = TRUE)



levels(dat.cpd.collate$group) # "Normal"      "IGT"         "T2D met pos" "T2D met neg"

dat.cpd.collate$ord_group <- factor(dat.cpd.collate$group, 
                                    levels = c("Normal", "IGT", "T2D met pos", "T2D met neg"),
                                    labels = c("1", "2", "3", "4"),ordered = TRUE)
dat.cpd.collate$ord_group <- as.integer(dat.cpd.collate$ord_group)
unique(dat.cpd.collate$ord_group) # 2 4 1 3

head(dat.cpd.collate)


saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")



str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260132" "ERR260132" "ERR260132" "ERR260132" ...
# $ cpd_rel_abun: num  0 5.647932 0 0 0.000186 ...
# $ log10_abun  : num  -8.482 0.752 -8.482 -8.482 -3.73 ...
# $ group       : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ group_label : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ ord_group   : int  2 2 2 2 2 2 2 2 2 2 ...


length( unique(dat.cpd.collate$cpd_id) ) # 7261
7261*145 # 1052845

#-------------------------


#### Forslund-SWE-T2D - CPP as phyloseq object
#    PCoA using CPP vs Functions
#-------------------------

data_in <- dat.cpd.collate.T2DNORM
str(data_in)
# 'data.frame':	551836 obs. of  7 variables:

length( unique(data_in$cpd_id) ) # 7261
length( unique(data_in$cpd_id[data_in$cpd_rel_abun > 0]) ) # 7031
length( unique(data_in$sample) ) # 76
7261*76 # 551836


### get data into phyloseq object ...

head(data_in)
#         cpd_id    sample cpd_rel_abun log10_abun       group group_label ord_group
# 50828 cpd24620 ERR260139 0.0000000000 -8.4820250 T2D met neg    T2D met-         4
# 50829 cpd00001 ERR260139 4.9744050062  0.6967411 T2D met neg    T2D met-         4
# 50830 cpd25681 ERR260139 0.0000000000 -8.4820250 T2D met neg    T2D met-         4
# 50831 cpd01501 ERR260139 0.0000000000 -8.4820250 T2D met neg    T2D met-         4
# 50832 cpd02597 ERR260139 0.0001838302 -3.7355832 T2D met neg    T2D met-         4
# 50833 cpd00851 ERR260139 0.0012068230 -2.9183564 T2D met neg    T2D met-         4

df.wide <- dcast(data_in, formula = sample + group_label ~ cpd_id , value.var = "cpd_rel_abun" )

df.wide[1:5, 1:10]
# sample group_label cpd00001 cpd00002  cpd00003  cpd00004  cpd00005  cpd00006   cpd00007 cpd00008
# 1 ERR260139    T2D met- 4.974405 3.235450 0.5465927 0.4750047 0.6423355 0.6418188 0.08380036 1.855040
# 2 ERR260140    T2D met- 5.843524 2.466131 0.5055806 0.4590337 0.5003057 0.5005527 0.06382374 1.418439
# 3 ERR260144    T2D met- 5.060686 3.206847 0.5376955 0.4692535 0.5938320 0.5946577 0.07292624 1.790827
# 4 ERR260147      Normal 4.784980 2.081487 0.6373313 0.6027531 0.4733389 0.4775388 0.14749815 1.298642
# 5 ERR260151    T2D met- 5.011624 2.989416 0.4911650 0.4350064 0.5745852 0.5751227 0.06705570 1.619055

unique(paste0(df.wide$sample,"--",df.wide$group_label))
#  [1] "ERR260139--T2D met-" "ERR260140--T2D met-" "ERR260144--T2D met-" "ERR260147--Normal"   "ERR260151--T2D met-"
# [6] "ERR260152--T2D met-" "ERR260153--Normal"   "ERR260159--T2D met-" "ERR260161--T2D met-" "ERR260162--T2D met-"
# [11] "ERR260163--Normal"   "ERR260165--T2D met-" "ERR260166--T2D met-" "ERR260167--T2D met-" "ERR260169--T2D met-"
# [16] "ERR260170--Normal"   "ERR260171--Normal"   "ERR260173--T2D met-" "ERR260174--T2D met-" "ERR260175--Normal"  
# [21] "ERR260179--T2D met-" "ERR260180--Normal"   "ERR260181--T2D met-" "ERR260185--T2D met-" "ERR260186--T2D met-"
# [26] "ERR260188--T2D met-" "ERR260189--T2D met-" "ERR260190--T2D met-" "ERR260193--Normal"   "ERR260198--T2D met-"
# [31] "ERR260199--T2D met-" "ERR260201--T2D met-" "ERR260203--T2D met-" "ERR260204--Normal"   "ERR260205--Normal"  
# [36] "ERR260206--T2D met-" "ERR260207--T2D met-" "ERR260209--Normal"   "ERR260210--T2D met-" "ERR260215--Normal"  
# [41] "ERR260216--Normal"   "ERR260217--Normal"   "ERR260218--Normal"   "ERR260221--Normal"   "ERR260223--Normal"  
# [46] "ERR260224--Normal"   "ERR260225--Normal"   "ERR260226--Normal"   "ERR260227--Normal"   "ERR260230--Normal"  
# [51] "ERR260231--Normal"   "ERR260234--Normal"   "ERR260241--T2D met-" "ERR260242--Normal"   "ERR260243--Normal"  
# [56] "ERR260244--Normal"   "ERR260246--Normal"   "ERR260250--Normal"   "ERR260251--Normal"   "ERR260252--Normal"  
# [61] "ERR260253--Normal"   "ERR260255--Normal"   "ERR260256--Normal"   "ERR260258--Normal"   "ERR260259--Normal"  
# [66] "ERR260260--Normal"   "ERR260263--Normal"   "ERR260264--Normal"   "ERR260265--Normal"   "ERR260266--Normal"  
# [71] "ERR260267--Normal"   "ERR260268--Normal"   "ERR260271--T2D met-" "ERR260273--T2D met-" "ERR260276--T2D met-"
# [76] "ERR275252--T2D met-"

# save group variable
samp <- df.wide[ ,1:2]
row.names(samp) <- samp$sample

# transpose
df.wide <- t(df.wide[ ,-2]) # minus 'group' column

head(df.wide)
# [,1]        [,2]        [,3]        [,4]        [,5]        [,6]        [,7]        [,8]        [,9]        [,10]      
# sample   "ERR260139" "ERR260140" "ERR260144" "ERR260147" "ERR260151" "ERR260152" "ERR260153" "ERR260159" "ERR260161" "ERR260162"
# cpd00001 "4.974405"  "5.843524"  "5.060686"  "4.784980"  "5.011624"  "4.989076"  "5.322874"  "5.874802"  "5.137622"  "5.105185" 
# cpd00002 "3.235450"  "2.466131"  "3.206847"  "2.081487"  "2.989416"  "3.108336"  "2.811890"  "2.553363"  "3.058027"  "3.032850" 
# cpd00003 "0.5465927" "0.5055806" "0.5376955" "0.6373313" "0.4911650" "0.5464072" "0.5231293" "0.4834641" "0.5765393" "0.5268389"
# cpd00004 "0.4750047" "0.4590337" "0.4692535" "0.6027531" "0.4350064" "0.4873457" "0.4703010" "0.4325929" "0.5190591" "0.4700231"
# cpd00005 "0.6423355" "0.5003057" "0.5938320" "0.4733389" "0.5745852" "0.6171960" "0.5335373" "0.4984751" "0.6426151" "0.5759009"
# [,11]       [,12]       [,13]       [,14]       [,15]       [,16]       [,17]       [,18]       [,19]       [,20]      
# sample   "ERR260163" "ERR260165" "ERR260166" "ERR260167" "ERR260169" "ERR260170" "ERR260171" "ERR260173" "ERR260174" "ERR260175"
# cpd00001 "5.143965"  "5.428956"  "5.919074"  "6.528776"  "5.591585"  "5.935547"  "5.948920"  "5.418080"  "5.389667"  "5.347925" 
# cpd00002 "3.048017"  "2.895433"  "2.637628"  "2.247076"  "2.692172"  "2.471387"  "2.643477"  "3.106761"  "2.956955"  "3.145791" 
# cpd00003 "0.5626358" "0.5560412" "0.5022639" "0.5077893" "0.4878529" "0.5509693" "0.5484417" "0.5302777" "0.4921912" "0.5256009"
# cpd00004 "0.5068864" "0.5001304" "0.4598232" "0.4622794" "0.4313041" "0.4996299" "0.5103907" "0.4783177" "0.4446042" "0.4581342"
# cpd00005 "0.5969206" "0.6009239" "0.5328602" "0.4533267" "0.5665398" "0.5085005" "0.5499079" "0.6126746" "0.5977709" "0.6134421"

samp_names <- df.wide[1, ]
tax_names <- row.names(df.wide[-1, ])
head(tax_names) # "cpd00001" "cpd00002" "cpd00003" "cpd00004" "cpd00005" "cpd00006"
otu.df <- df.wide[-1, ] # remove sample labels in 1st row
# this is necessary to create numeric matrix

colnames(otu.df) <- samp_names

# convert OTU table to matrix
class(otu.df) # "matrix" "array" 
#otu.df <- as.matrix(otu.df)

# convert to numeric matrix
# https://stackoverflow.com/questions/20791877/convert-character-matrix-into-numeric-matrix
otu.df <- apply(otu.df, 2, as.numeric)

rownames(otu.df) # NULL
dim(otu.df) #  7261   76
rownames(otu.df) <- tax_names

## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
OTU <- otu_table(otu.df, taxa_are_rows = TRUE)


# # convert Taxonomy table to matrix  

tax <- data.frame(cpd_id = tax_names)
row.names(tax) <- tax_names

tax <- as.matrix(tax)

identical( row.names(otu.df), row.names(tax) ) # TRUE


## Create 'taxonomyTable'
#  tax_table - Works on any character matrix.
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
TAX <- tax_table(tax)


## Create a phyloseq object, merging OTU & TAX tables
phy.cpp = phyloseq(OTU, TAX)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7261 taxa and 76 samples ]
# tax_table()   Taxonomy Table:    [ 7261 taxa by 1 taxonomic ranks ]


sample_names(phy.cpp)
# [1] "ERR260139" "ERR260140" "ERR260144" "ERR260147" ... etc.

identical(sample_names(phy.cpp), samp$sample) # TRUE


# row.names need to match sample_names() from phyloseq object
row.names(samp) <- samp$sample


### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

SAMP <- sample_data(samp)


### Combine SAMPDATA into phyloseq object
phy.cpp <- merge_phyloseq(phy.cpp, SAMP)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7261 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 7261 taxa by 1 taxonomic ranks ]

# But see below - only 7031 taxa - after remove taxa with zero reads (fix legacy from extra classes T2D Met+ and IGT)

phy.cpp@sam_data
# Sample Data:        [76 samples by 2 sample variables]:
# sample group_label
# ERR260139 ERR260139    T2D met-
# ERR260140 ERR260140    T2D met-
# ERR260144 ERR260144    T2D met-
# ERR260147 ERR260147      Normal
# ERR260151 ERR260151    T2D met-
# ERR260152 ERR260152    T2D met-
# ERR260153 ERR260153      Normal
# ERR260159 ERR260159    T2D met-
# ERR260161 ERR260161    T2D met-
# ERR260162 ERR260162    T2D met-
# ERR260163 ERR260163      Normal
# etc. ...

T2DNorm.samps <- as.data.frame(phy.cpp@sam_data)


phy_in <- phy.cpp

phy_in
# as above 

min(taxa_sums(phy_in)) # 0
sort(taxa_sums(phy_in))[1:1000]

## NOTE THERE ARE SOME Compounds WITH zero relative abundance ... likely a carryover from IGT and Met+ classes

# prune taxa that have zero sequence reads
phy_in <- prune_taxa(taxa = taxa_sums(phy_in) > 0, x = phy_in)
phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7031 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 7031 taxa by 1 taxonomic ranks ]

saveRDS(object = phy_in, file = "phy.cpp-cleaned-Forslund-SWE-T2D-v3.RDS")


sum(sample_sums(phy_in)) # 5361.087
sample_sums(phy_in)
# ERR260139 ERR260140 ERR260144 ERR260147 ERR260151 ERR260152 ERR260153 ERR260159 ERR260161 ERR260162 ERR260163 ERR260165 ERR260166 
# 74.74261  66.15145  74.32332  59.87642  69.96858  72.87472  69.27684  67.49364  73.71519  72.42302  71.95372  71.45717  68.92647 
# ERR260167 ERR260169 ERR260170 ERR260171 ERR260173 ERR260174 ERR260175 ERR260179 ERR260180 ERR260181 ERR260185 ERR260186 ERR260188 
# 65.60642  70.39758  68.21662  69.81802  73.84764  71.39022  74.05353  73.07446  75.32151  68.92975  73.20277  72.15931  72.55811 
# ERR260189 ERR260190 ERR260193 ERR260198 ERR260199 ERR260201 ERR260203 ERR260204 ERR260205 ERR260206 ERR260207 ERR260209 ERR260210 
# 72.21470  70.59005  70.96755  72.15018  71.46823  71.79915  67.38239  68.91987  73.16533  72.18530  73.94387  72.74923  72.22642 
# ERR260215 ERR260216 ERR260217 ERR260218 ERR260221 ERR260223 ERR260224 ERR260225 ERR260226 ERR260227 ERR260230 ERR260231 ERR260234 
# 73.35239  66.25152  70.92966  67.75426  65.44894  71.10523  70.26225  70.61357  73.37748  65.08487  66.38871  69.72384  64.54567 
# ERR260241 ERR260242 ERR260243 ERR260244 ERR260246 ERR260250 ERR260251 ERR260252 ERR260253 ERR260255 ERR260256 ERR260258 ERR260259 
# 71.63385  69.66454  73.15731  72.48794  68.91389  67.89545  69.15521  70.20690  72.56340  68.90731  69.17781  70.07831  71.12961 
# ERR260260 ERR260263 ERR260264 ERR260265 ERR260266 ERR260267 ERR260268 ERR260271 ERR260273 ERR260276 ERR275252 
# 71.48764  71.76574  68.93297  63.89533  72.99892  73.05355  70.20026  74.79608  71.40596  69.18567  72.03342

summary( sample_sums(phy_in) )
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 59.88   68.93   71.12   70.54   72.56   75.32

sd( sample_sums(phy_in) )
# 2.875939

max(taxa_sums(phy_in)) # 500.2104


# don't rarefy - already in form of relative abundance %


table(phy_in@sam_data$group_label)
# Normal T2D met- 
#   43       33 


## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
#ord <- ordinate(r1.ps, "NMDS", "bray")
ord <- ordinate(phy_in, "PCoA", "bray")


ord

unique(phy_in@sam_data$group_label)
# [1] T2D met- Normal  
# Levels: Normal < T2D met-


p <- plot_ordination(phy_in, ord, type="samples", color="group_label")
p

p$labels$x # "Axis.1   [40.6%]"
x_lab <- "PCo1 (40.6%)"

p$labels$y # "Axis.2   [22.9%]"
y_lab <- "PCo2 (22.9%)"

40.6 + 22.9 # 63.5


#temp <- r1.ps
p_df <- p$data


cols.group <- c("Normal" = "#7fbf7b",
                "T2D met-" = "#af8dc3")



p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = group_label))+
  theme_bw()+
  geom_point()+
  
  stat_ellipse(linetype = "dashed")+
  
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Diagnosis") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "T2D: CPP", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(d)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

#dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-Forslund-SWE-T2D-v3.tiff"), width = 9, height = 8, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-Forslund-SWE-T2D-v3.tiff"), width = 9.5, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# "sample" "group_label" 

# Adonis test
set.seed(123)
adonis2(bray ~ group_label , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group_label, data = sampledf)
#             Df SumOfSqs      R2      F Pr(>F)
# group_label  1  0.00816 0.01972 1.4884  0.174
# Residual    74  0.40554 0.98028              
# Total       75  0.41370 1.00000



beta <- betadisper(bray, sampledf$group)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.000001 1.040e-06 9e-04    999  0.979
# Residuals 74 0.088576 1.197e-03



### phy.cpp - but test only compounds that consistently trend with disturbed soils and in T2D?


phy_in <- readRDS("phy.cpp-cleaned-Forslund-SWE-T2D-v3.RDS")
phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7031 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 7031 taxa by 1 taxonomic ranks ]

# get consistently trending compounds?

dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.rds")
head(dat.test.sig)
unique(dat.test.sig$trend_group)
# [1] Decreasing in T2D (reduced exposure in quality ecosystems)     Increasing in T2D (increased exposure in disturbed ecosystems)
# [3] Decreasing in T2D (reduced exposure in disturbed ecosystems)   Increasing in T2D (increased exposure in quality ecosystems)  
# 4 Levels: Decreasing in T2D (reduced exposure in quality ecosystems) < ...
sel <- which(dat.test.sig$trend_group %in% c("Decreasing in T2D (reduced exposure in disturbed ecosystems)",
                                             "Increasing in T2D (increased exposure in disturbed ecosystems)"))
# qty 128
keep_taxa <- dat.test.sig$cpd[sel]

phy_in <- prune_taxa( taxa = keep_taxa, x = phy_in)
min(taxa_sums(phy_in)) # 7.551019e-05

table(phy_in@sam_data$group_label)
# Normal T2D met- 
#   43       33 

## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
ord <- ordinate(phy_in, "PCoA", "bray")
ord
unique(phy_in@sam_data$group_label)
# [1] T2D met- Normal  
# Levels: Normal < T2D met-

p <- plot_ordination(phy_in, ord, type="samples", color="group_label")
p

p$labels$x # "Axis.1   [58%]"
x_lab <- "PCo1 (58%)"

p$labels$y # "Axis.2   [20.8%]"
y_lab <- "PCo2 (20.8%)"

58 + 20.8 # 78.8

p_df <- p$data

cols.group <- c("Normal" = "#7fbf7b",
                "T2D met-" = "#af8dc3")

p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = group_label))+
  theme_bw()+
  geom_point()+
  
  stat_ellipse(linetype = "dashed")+
  
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Diagnosis") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D: CPP\nConsistent with\nsoil trends", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

#dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-Forslund-SWE-T2D-v3.tiff"), width = 9, height = 8, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-T2D-Soil-Consistent-trends-Forslund-SWE-T2D-v3.tiff"), width = 10, height = 8, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","a-PCoA-Cpp3d-T2D-Soil-Consistent-trends-Forslund-SWE-T2D-v3.tiff"), width = 10, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# "sample" "group_label" 

# Adonis test
set.seed(123)
adonis2(bray ~ group_label , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group_label, data = sampledf)
#             Df SumOfSqs      R2      F Pr(>F)   
# group_label  1  0.02414 0.05675 4.4521   0.01 **
# Residual    74  0.40117 0.94325                 
# Total       75  0.42531 1.00000                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



beta <- betadisper(bray, sampledf$group_label)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.000026 2.603e-05 0.0175    999  0.909
# Residuals 74 0.109889 1.485e-03 



### test all trending compounds ?

phy_in <- readRDS("phy.cpp-cleaned-Forslund-SWE-T2D-v3.RDS")
phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7031 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 7031 taxa by 1 taxonomic ranks ]

# get consistently trending compounds?

dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.rds")
head(dat.test.sig)
dim(dat.test.sig) # 276  25
unique(dat.test.sig$trend_group)
# [1] Decreasing in T2D (reduced exposure in quality ecosystems)     Increasing in T2D (increased exposure in disturbed ecosystems)
# [3] Decreasing in T2D (reduced exposure in disturbed ecosystems)   Increasing in T2D (increased exposure in quality ecosystems)  
# 4 Levels: Decreasing in T2D (reduced exposure in quality ecosystems) < ...

table(dat.test.sig$trend_group, useNA = "ifany" )
# Decreasing in T2D (reduced exposure in quality ecosystems)   Decreasing in T2D (reduced exposure in disturbed ecosystems) 
# 98                                                             70 
# Increasing in T2D (increased exposure in disturbed ecosystems)   Increasing in T2D (increased exposure in quality ecosystems) 
# 58                                                             50 

# sel <- which(dat.test.sig$trend_group %in% c("Decreasing in T2D (reduced exposure in disturbed ecosystems)",
#                                              "Increasing in T2D (increased exposure in disturbed ecosystems)"))

# these all have a relationship in t2d and with soil condition
keep_taxa <- dat.test.sig$cpd

phy_in <- prune_taxa( taxa = keep_taxa, x = phy_in)
min(taxa_sums(phy_in)) # 7.551019e-05

table(phy_in@sam_data$group_label)
# Normal T2D met- 
#   43       33 

## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
ord <- ordinate(phy_in, "PCoA", "bray")
ord
unique(phy_in@sam_data$group_label)
# [1] T2D met- Normal  
# Levels: Normal < T2D met-

p <- plot_ordination(phy_in, ord, type="samples", color="group_label")
p

p$labels$x # "Axis.1   [53.7%]"
x_lab <- "PCo1 (53.7%)"

p$labels$y # "Axis.2   [18.1%]"
y_lab <- "PCo2 (18.1%)"

53.7 + 18.1 # 71.8

p_df <- p$data

cols.group <- c("Normal" = "#7fbf7b",
                "T2D met-" = "#af8dc3")

p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = group_label))+
  theme_bw()+
  geom_point()+
  
  stat_ellipse(linetype = "dashed")+
  
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Diagnosis") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "T2D: CPP\nAny trend in\nT2D and soil", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

#dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-Forslund-SWE-T2D-v3.tiff"), width = 9, height = 8, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-T2D-Any-trends-T2D-soil-Forslund-SWE-T2D-v3.tiff"), width = 10, height = 8, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","b-PCoA-Cpp3d-T2D-Any-trends-T2D-soil-Forslund-SWE-T2D-v3.tiff"), width = 10, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# "sample" "group_label" 

# Adonis test
set.seed(123)
adonis2(bray ~ group_label , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group_label, data = sampledf)
#             Df SumOfSqs     R2      F Pr(>F)   
# group_label  1  0.03077 0.0595 4.6813  0.008 **
# Residual    74  0.48643 0.9405                 
# Total       75  0.51720 1.0000                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(bray, sampledf$group_label)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.00020 0.00020047 0.1232    999  0.756
# Residuals 74 0.12042 0.00162728




## Functions

phy_in <- readRDS("phy-phyloseq-object-Forslund-SWE-T2D.RDS")

phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 19099 taxa and 145 samples ]
# sample_data() Sample Data:       [ 145 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 19099 taxa by 4 taxonomic ranks ]

# only analyse T2D Met- and Normal

keep_samps <- unique(T2DNorm.samps$sample)

phy_in <- prune_samples(keep_samps, phy_in)

head( phy_in@sam_data )
# Sample Data:        [6 samples by 5 sample variables]:
#   Sample Country.subset         Status      Bases       Run
# ERR260139 NG-5636_334            SWE T2D metformin- 2036676514 ERR260139
# ERR260140 NG-5636_344            SWE T2D metformin- 1935856900 ERR260140
# ERR260144 NG-5636_353            SWE T2D metformin- 2483902494 ERR260144
# ERR260147 NG-5636_365            SWE        ND CTRL 2821768300 ERR260147
# ERR260151 NG-5636_378            SWE T2D metformin- 2630431274 ERR260151
# ERR260152 NG-5636_380            SWE T2D metformin- 1813559434 ERR260152

identical(row.names(phy_in@sam_data),T2DNorm.samps$sample ) # TRUE

phy_in@sam_data$group_label <- T2DNorm.samps$group_label
phy_in@sam_data


min(taxa_sums(phy_in)) # 0


# prune taxa that have zero sequence reads
phy_in <- prune_taxa(taxa = taxa_sums(phy_in) > 0, x = phy_in)
phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17962 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 17962 taxa by 4 taxonomic ranks ]


sum(sample_sums(phy_in)) # 7600
sample_sums(phy_in)
# all 100

summary( sample_sums(phy_in) )

max(taxa_sums(phy_in)) # 232.0692


# do not rarefy as already normalized using rel abun (%)


table(phy_in@sam_data$group_label)
# Normal T2D met- 
#   43       33

phy_in@sam_data$group_label
# [1] T2D met- T2D met- T2D met- Normal   T2D met- T2D met- Normal   T2D met- T2D met- T2D met- Normal   T2D met- T2D met- T2D met-
# [15] T2D met- Normal   Normal   T2D met- T2D met- Normal   T2D met- Normal   T2D met- T2D met- T2D met- T2D met- T2D met- T2D met-
# [29] Normal   T2D met- T2D met- T2D met- T2D met- Normal   Normal   T2D met- T2D met- Normal   T2D met- Normal   Normal   Normal  
# [43] Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal   T2D met- Normal   Normal   Normal  
# [57] Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal   Normal  
# [71] Normal   Normal   T2D met- T2D met- T2D met- T2D met-
# Levels: Normal < T2D met-



## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
#ord <- ordinate(r1.ps, "NMDS", "bray")
ord <- ordinate(phy_in, "PCoA", "bray")


ord

unique(phy_in@sam_data$group_label)
# [1] T2D met- Normal  
# Levels: Normal < T2D met-


p <- plot_ordination(phy_in, ord, type="samples", color="group_label")
p

p$labels$x # "Axis.1   [28.3%]"
x_lab <- "PCo1 (28.3%)"

p$labels$y # "Axis.2   [21.4%]"
y_lab <- "PCo2 (21.4%)"

#temp <- r1.ps
p_df <- p$data


cols.group <- c("Normal" = "#7fbf7b",
                "T2D met-" = "#af8dc3")



p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = group_label))+
  theme_bw()+
  geom_point()+
  
  stat_ellipse(linetype = "dashed")+
  
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Diagnosis") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "T2D: Fxns", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(c)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

#dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Functions-Forslund-SWE-T2D-v3.tiff"), width = 9, height = 8, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Functions-Forslund-SWE-T2D-v3.tiff"), width = 9.5, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# "Sample"         "Country.subset" "Status"         "Bases"          "Run"            "group_label"   

# Adonis test
set.seed(123)
# adonis2(bray ~ group_label , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group_label, data = sampledf)
#             Df SumOfSqs      R2      F Pr(>F)
# group_label  1  0.05283 0.02054 1.5516  0.136
# Residual    74  2.51948 0.97946              
# Total       75  2.57231 1.00000 



beta <- betadisper(bray, sampledf$group_label)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.00086 0.0008571 0.1555    999  0.712
# Residuals 74 0.40779 0.0055107



#-------------------------


#### Forslund-SWE-T2D - MANN-WHITNEY TESTS - CPP FOR ALL AVAILABLE COMPOUNDS
#    Test sig diff Normal vs T2D only  - keep all p <= 0.05
#-------------------------

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260132" "ERR260132" "ERR260132" "ERR260132" ...
# $ cpd_rel_abun: num  0 5.647932 0 0 0.000186 ...
# $ log10_abun  : num  -8.482 0.752 -8.482 -8.482 -3.73 ...
# $ group       : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ group_label : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ ord_group   : int  2 2 2 2 2 2 2 2 2 2 ...

# select only Normal and T2D
unique(dat.cpd.collate$group)
# [1] IGT         T2D met neg Normal      T2D met pos
# Levels: Normal < IGT < T2D met pos < T2D met neg

sel <- which(dat.cpd.collate$group %in% c("T2D met neg", "Normal"))

dat.cpd.collate.T2DNORM <- dat.cpd.collate[sel, ]

length( unique(dat.cpd.collate.T2DNORM$cpd_id) ) # 7261
length( unique(dat.cpd.collate.T2DNORM$sample) ) # 76

data_in <- dat.cpd.collate.T2DNORM

unique(data_in$group) # T2D met neg Normal     

dat.test <- data.frame(cpd = unique(dat.cpd.collate.T2DNORM$cpd_id), data_for_this_cpd=NA , p_val = NA, mean_t2d = NA, mean_Normal = NA, t_statistic = NA, estimate_log_diff = NA, perc_diff = NA, trend_with_disease = NA )

dat.test <- data.frame(cpd = unique(data_in$cpd_id), data_for_this_cpd=NA , 
                       median_t2d = NA, median_normal = NA,  
                       alt = NA,
                       p_val = NA, W_statistic = NA, hl_effect_wilcox = NA,
                       trend_with_disease = NA
)

for (i in 1:dim(dat.test)[1]) {
  #i<-1
  this_cpd <- dat.test$cpd[i]
  sel <- which(data_in$cpd_id == this_cpd)
  # prepare data
  df = data.frame(group = as.character( data_in$group[sel]), 
                  value = data_in$log10_abun[sel], 
                  value_perc = data_in$cpd_rel_abun[sel] )
  
  x <- df$value[df$group == "T2D met neg"]
  y <- df$value[df$group == "Normal"]
  
  if ( length(which(x == min(x) )) > 0.5*length(x) & length(which(y == min(y) )) > 0.5*length(y) ) {
    # low data: if at least one dataset does not have at least 50% non-zero cases
    dat.test$data_for_this_cpd[i] <- "low data"
    
  } else {
    
    # # test for homogeneity of variances
    # var.test(x,y) # e.g., j<-100: F = 2.3477, num df = 32, denom df = 42, p-value = 0.009872
    
    dat.test$median_t2d[i] <- median( x ) # log10 abun  
    dat.test$median_normal[i] <- median( y ) # log10 abun
    
    alt <- NA
    if (median( x ) > median( y ) ) {
      alt <- "greater"
    } else {
      alt <- "less"
    }
    dat.test$alt[i] <- alt
    
    # Wilcoxon-Mann-Whitney Test
    wmw.test <- wilcox.test(x, y, alternative = alt, paired = FALSE, conf.int = TRUE) # based on log10 abun
    
    dat.test$p_val[i] <- wmw.test$p.value
    dat.test$W_statistic[i] <- wmw.test$statistic
    
    # https://search.r-project.org/CRAN/refmans/DescTools/html/HodgesLehmann.html
    # https://aakinshin.net/posts/r-hodges-lehmann-problems
    
    dat.test$hl_effect_wilcox[i] <- wmw.test$estimate
    #dat.test$hl_effect_hlfxn[i] <- hl(x, y)
    
    if (!(is.na(dat.test$p_val[i])|is.na(dat.test$W_statistic[i]))) {
      if (dat.test$p_val[i] <= 0.05 & alt == "greater") { dat.test$trend_with_disease[i] <- "Increasing" }
      else if (dat.test$p_val[i] <= 0.05 & alt == "less") { dat.test$trend_with_disease[i] <- "Decreasing" }
      else { dat.test$trend_with_disease[i] <- "No trend" }
    }
  }
  print(paste0("Completed ",i))
}


sel.low <- which(dat.test$data_for_this_cpd == "low data") # 1709

length(unique(dat.cpd.collate.T2DNORM$cpd_id)) # 7261
length(unique(dat.cpd.collate.T2DNORM$cpd_id)) - length(sel.low) # 5552

sel.nonNA <- which(!is.na(dat.test$p_val)) # 5552 applicable tests


# what are most interesting compounds - that decrease in restoration (increase with land disturbance) and increase with T2D ??
filter( dat.test[ order(dat.test$W_statistic, decreasing = TRUE),  ], decResto_incT2D == 1 )

#dat.test %>% filter(incResto_incT2D == 1) %>% arrange(., p_val) %>% slice(1:10) %>% pull(W_statistic) %>% mean() 
dat.test %>% arrange(., p_val) %>% slice(1:50)
# cpd data_for_this_cpd median_t2d median_normal     alt        p_val W_statistic hl_effect_wilcox trend_with_disease
# 1  cpd27894              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 2  cpd02158              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 3  cpd20917              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 4  cpd23905              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 5  cpd28017              <NA> -2.8952325    -2.9711761 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 6  cpd28332              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 7  cpd24099              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 8  cpd01826              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 9  cpd26049              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 10 cpd27062              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 11 cpd22506              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 12 cpd23628              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 13 cpd25951              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 14 cpd25952              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 15 cpd25953              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 16 cpd09713              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 17 cpd27027              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 18 cpd29256              <NA> -3.1962625    -3.2722061 greater 0.0006807261      1011.0     7.716439e-02         Increasing
# 19 cpd23904              <NA> -3.1962625    -3.2722061 greater 0.0007074925      1010.0     7.456159e-02         Increasing
# 20 cpd02894              <NA> -1.7941598    -1.8490690 greater 0.0007074925      1010.0     5.612854e-02         Increasing
# 21 cpd21088              <NA> -3.3534294    -2.9736340    less 0.0008034725       408.0    -3.296790e-01         Decreasing
# 22 cpd02103              <NA> -8.4820250    -4.7306480    less 0.0008046816       446.0    -4.926577e-05         Decreasing
# 23 cpd02662              <NA> -8.4820250    -7.1224077    less 0.0010740946       442.0    -5.571444e-01         Decreasing
# 24 cpd03725              <NA> -2.8922259    -2.9623715 greater 0.0013354952       993.0     7.168170e-02         Increasing
# 25 cpd00465              <NA> -8.4820250    -4.1219361    less 0.0016839162       444.5    -1.071271e+00         Decreasing
# 26 cpd25874              <NA> -8.4820250    -6.3322478    less 0.0016922065       440.0    -1.376087e+00         Decreasing
# 27 cpd25883              <NA> -8.4820250    -6.3322478    less 0.0016922065       440.0    -1.376087e+00         Decreasing
# 28 cpd04326              <NA> -5.0509443    -4.2042575    less 0.0018381743       434.0    -7.890715e-01         Decreasing
# 29 cpd00397              <NA> -1.2743146    -1.3449430 greater 0.0019068864       983.0     7.018496e-02         Increasing
# 30 cpd01059              <NA> -8.4820250    -6.0783314    less 0.0020113046       444.0    -8.939130e-01         Decreasing
# 31 cpd23041              <NA> -3.3805631    -3.0416911    less 0.0021356790       436.5    -3.207468e-01         Decreasing
# 32 cpd23037              <NA> -3.3805631    -3.0416911    less 0.0021356790       436.5    -3.207468e-01         Decreasing
# 33 cpd01239              <NA> -3.7367518    -3.5740411    less 0.0022553076       438.0    -2.445223e-01         Decreasing
# 34 cpd12566              <NA> -8.4820250    -5.0879427    less 0.0022625690       470.0    -8.027574e-05         Decreasing
# 35 cpd28212              <NA> -8.4820250    -5.0879427    less 0.0022625690       470.0    -8.027574e-05         Decreasing
# 36 cpd21918              <NA> -8.4820250    -5.0879427    less 0.0022625690       470.0    -8.027574e-05         Decreasing
# 37 cpd28215              <NA> -8.4820250    -4.6487704    less 0.0023269553       455.0    -1.152886e+00         Decreasing
# 38 cpd30007              <NA> -3.7402571    -3.5748068    less 0.0023305134       439.0    -2.439347e-01         Decreasing
# 39 cpd03722              <NA> -3.7402571    -3.5748068    less 0.0023305134       439.0    -2.439347e-01         Decreasing
# 40 cpd29990              <NA> -3.7189127    -3.4346533    less 0.0024079826       440.0    -2.982317e-01         Decreasing
# 41 cpd00039              <NA> -0.9908031    -0.9741944    less 0.0024284195       443.0    -3.034832e-02         Decreasing
# 42 cpd28591              <NA> -8.4820250    -4.3784655    less 0.0024368204       472.0    -4.131022e-05         Decreasing
# 43 cpd28584              <NA> -8.4820250    -4.3784655    less 0.0024368204       472.0    -4.131022e-05         Decreasing
# 44 cpd00847              <NA> -8.4820250    -4.9524968    less 0.0024460548       466.0    -2.963855e-01         Decreasing
# 45 cpd22343              <NA> -8.4820250    -7.1581522    less 0.0026352075       472.0    -6.750047e-05         Decreasing
# 46 cpd04041              <NA> -8.4820250    -7.1581522    less 0.0028332583       474.0    -2.581399e-05         Decreasing
# 47 cpd04042              <NA> -8.4820250    -7.1581522    less 0.0028332583       474.0    -2.581399e-05         Decreasing
# 48 cpd06715              <NA> -8.4820250    -7.1581522    less 0.0028332583       474.0    -2.581399e-05         Decreasing
# 49 cpd14978              <NA> -8.4820250    -7.1581522    less 0.0028332583       474.0    -2.581399e-05         Decreasing
# 50 cpd14979              <NA> -8.4820250    -6.8571222    less 0.0028332583       474.0    -3.953920e-05         Decreasing

write.csv(x = dat.test, file = "dat.test-cpp3d-indiv-Compounds-Forslund-SWE-T2D--All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c.csv")

saveRDS(dat.test, file = "dat.test-cpp3d-indiv-Compounds-Forslund-SWE-T2D--All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c.RDS")

# only keep applicable tests; the remainder are likely due to zero replacement:
#dat.test <- dat.test[ sel.nonNA, ]

# extract sig results ... (no p-adjustment)

sel.sig <- which(dat.test$p_val <= 0.05) # 1243

dat.test.sig <- dat.test[sel.sig, ]

dat.test.sig$minuslog10_p_val <- -log10(dat.test.sig$p_val)

plot(x = dat.test.sig$W_statistic , y =dat.test.sig$minuslog10_p_val , xlab="W statistic", ylab="-log10(P-value)")

dev.print(tiff, filename = paste0(workdir,"/plots/","3d-indiv-compound-VolcanoPlot-P-values--All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c.tiff"),
          width = 12, height = 14, units = "cm", res=600, compression = "lzw",type="cairo" )


# join compound info??

dat.test.sig$cpd_names <- NA
dat.test.sig$cpd_forms <- NA

dat.test.sig$OC_x <- NA
dat.test.sig$HC_y <- NA
dat.test.sig$NC_z <- NA

dat.test.sig$mass <- NA
dat.test.sig$class <- NA

for (i in 1:dim(dat.test.sig)[1]) {
  #i<-1
  this_cpd <- dat.test.sig$cpd[i]

  sel.cpd <- which(df.comp2$id == this_cpd)

  dat.test.sig$cpd_names[i] <- df.comp2$name[sel.cpd]
  dat.test.sig$cpd_forms[i] <- df.comp2$form[sel.cpd]

  dat.test.sig$OC_x[i] <- df.comp2$OC_ratio[sel.cpd]
  dat.test.sig$HC_y[i] <- df.comp2$HC_ratio[sel.cpd]
  dat.test.sig$NC_z[i] <- df.comp2$NC_ratio[sel.cpd]

  dat.test.sig$mass[i] <- df.comp2$mass[sel.cpd]
  dat.test.sig$class[i] <- df.comp2$class[sel.cpd]

  print(paste0("completed ",i))
}

write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c.rds")


hist(dat.test.sig$NC_z); summary(dat.test.sig$NC_z)


dim(dat.test.sig) # 1243   17
sel.ok <- which(!is.na(dat.test.sig$NC_z) ) # qty 1158
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.0000  0.0000  0.1069  0.1667  1.0000      85 

dat.test.sig$z_layer <- NA

subsel <- which(dat.test.sig$NC_z[sel.ok] == 0) # 616
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C = 0"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0 & dat.test.sig$NC_z[sel.ok] <= 0.2 ) # 306
# Rivas-Ubach et al 2018 Table 1 highlights N/C breaks at 0.126, 0.2, 0.5, 0.7
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0 to 0.2"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0.2) # 236
max(dat.test.sig$NC_z[sel.ok]) # 1
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 1"

unique(dat.test.sig$z_layer[sel.ok]) # "N:C = 0"       "N:C >0 to 0.2" "N:C >0.2 to 1"

dat.test.sig$z_layer <- factor(dat.test.sig$z_layer, levels = c("N:C = 0",
                                                                "N:C >0 to 0.2",
                                                                "N:C >0.2 to 1"), ordered = TRUE)

saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c.rds")

# HERE !!

dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c.rds")


dim(dat.test.sig) #  1243   18
head(dat.test.sig)
dim(dat.test.sig[ which(dat.test.sig$p_val < 0.05), ]) # 1243   18
sel <- which(!is.na(dat.test.sig$OC_x) ) # 1158
sel <- which(dat.test.sig$trend_with_disease == "Decreasing") # 915
sel <- which(dat.test.sig$trend_with_disease == "Decreasing" & !is.na(dat.test.sig$OC_x) ) # 864
sel <- which(dat.test.sig$trend_with_disease == "Increasing") # 328
sel <- which(dat.test.sig$trend_with_disease == "Increasing" & !is.na(dat.test.sig$OC_x) ) # 294



## plot as Increasing or Decreasing?? in vK space
## Use adjusted compound classes (adapted from Wu 2018, D'Andrilli, Rivas-Ubach 2018, and Minor et al 2015)

# # Zones and labels as above
# vkgrouprect.facets2 <- read.table(file = "cpp3d-compound-classes.tsv", header = TRUE, sep = "\t" )
# vkgrouprect.facets2.labels <- read.table(file = "cpp3d-compound-classes-labels.tsv", header = TRUE, sep = "\t" )




p <- ggplot(data = dat.test.sig[sel.ok, ]) +
  coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study")+
  #xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + #
  xlab("O:C ratio")+ ylab("H:C ratio")+
  
  guides(color = guide_legend(title = "Trend with disease in functional capacity\n(%) allocated to compounds"))+
  
  facet_wrap(facets = vars(z_layer))+
  
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  #geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 )+ # 
  
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T,  size = 2 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2 , col="#737373" , lineheight = 0.8)+
  #geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2 , col="#737373" , lineheight = 0.8)+
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.9)) ,
    title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Disease----All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-large.tiff"), width = 20, height = 12, units = "cm", res=500, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Disease----All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-lim.tiff"), width = 22, height = 12, units = "cm", res=500, compression="lzw",type="cairo")

dev.print(tiff, file = paste0(workdir,"/plots/","noGgtitle-3d-Compounds-indiv-vKSpace-Trend-with-Disease----All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-lim.tiff"), width = 20, height = 12, units = "cm", res=500, compression="lzw",type="cairo")

#-------------------------


#### Tests - Forslund-SWE-T2D
#    Test sig diff Normal vs T2D only, with subset informed by soil quality trends  
#    FOCUS ONLY ON SIGNIFICANTLY DIFFERENT COMPOUNDS PREVIOUSLY SELECTED FROM ECOSYSTEM RESTORATION (P-ADJUSTED)
#    NON-PARAMETRIC (Wilcoxon-Mann-Whitney Test) and threshold for 'low data' exclusions
#-------------------------

#saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.rds")
dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.rds")
max(dat.test.sig$p_val) # 0.04552315

str(dat.test.sig)
# 'data.frame':	2958 obs. of  14 variables:
# $ cpd              : chr  "cpd25681" "cpd24620" "cpd00001" "cpd00851" ...
# $ data_for_this_cpd: chr  NA NA NA NA ...
# $ p_val            : num  0.000575 0.015075 0.015075 0.006249 0.006249 ...
# $ kendall_tau      : num  -3.44 -2.43 -2.43 -2.73 -2.73 ...
# $ trend_with_age   : chr  "Decreasing" "Decreasing" "Decreasing" "Decreasing" ...
# $ sigBH            : chr  "sig" NA NA "sig" ...
# $ minuslog10_p_val : num  3.24 1.82 1.82 2.2 2.2 ...
# $ cpd_names        : chr  "2-methyl-trans-aconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "trans-4-Hydroxy-L-proline" ...
# $ cpd_forms        : chr  "C7H5O6" "C7H7O7" "H2O" "C5H9NO3" ...
# $ OC_x             : num  0.857 1 NA 0.6 0.6 ...
# $ HC_y             : num  0.714 1 NA 1.8 1.8 ...
# $ NC_z             : num  0 0 NaN 0.2 0.2 ...
# $ z_layer          : Ord.factor w/ 3 levels "N:C = 0"<"N:C >0 to 0.2"<..: 1 1 NA 2 2 3 3 3 3 NA ...
# $ mass             : num  188 206 18 130 130 89 175 147 147 1 ...

sel.sigBH <- which( dat.test.sig$sigBH == "sig") # 2122

resto.sig.cpds <- unique(dat.test.sig$cpd[sel.sigBH]) # qty 2122

dat.test.sig.resto <- dat.test.sig

saveRDS(resto.sig.cpds, file = "resto.sig.cpds-sigBH-from-sunbad-resto-v2.rds")


#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260132" "ERR260132" "ERR260132" "ERR260132" ...
# $ cpd_rel_abun: num  0 5.647932 0 0 0.000186 ...
# $ log10_abun  : num  -8.482 0.752 -8.482 -8.482 -3.73 ...
# $ group       : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ group_label : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ ord_group   : int  2 2 2 2 2 2 2 2 2 2 ...

length(unique(dat.cpd.collate$cpd_id)) # 7261
length(unique(dat.cpd.collate$group)) # 4
length(unique(dat.cpd.collate$sample)) # 145
7261*145 # 1052845
table(dat.cpd.collate$group)/7261
# Normal         IGT T2D met pos T2D met neg 
#     43          49          20          33 

# select only Normal and T2D
unique(dat.cpd.collate$group)
# [1] IGT         T2D met neg Normal      T2D met pos
# Levels: Normal < IGT < T2D met pos < T2D met neg

sel <- which(dat.cpd.collate$group %in% c("T2D met neg", "Normal")) # 
dat.cpd.collate.T2DNORM <- dat.cpd.collate[sel, ]

length( unique(dat.cpd.collate.T2DNORM$cpd_id) ) # 7261 
table(dat.cpd.collate.T2DNORM$group)/7261
# Normal         IGT T2D met pos T2D met neg 
# 43           0           0          33 

length( unique(dat.cpd.collate.T2DNORM$sample) ) # 76

sel <- which(dat.cpd.collate.T2DNORM$cpd_id %in% resto.sig.cpds)
dat.cpd.collate.T2DNORM <- dat.cpd.collate.T2DNORM[sel, ]

length( unique(dat.cpd.collate.T2DNORM$cpd_id) ) # 2052
table(dat.cpd.collate.T2DNORM$group)/2052
# Normal         IGT T2D met pos T2D met neg 
#     43           0           0          33

length( unique(dat.cpd.collate.T2DNORM$sample[ dat.cpd.collate.T2DNORM$group == "T2D met neg" ]) ) # 33
length( unique(dat.cpd.collate.T2DNORM$sample[ dat.cpd.collate.T2DNORM$group == "Normal" ]) ) # 43

data_in <- dat.cpd.collate.T2DNORM

dat.test <- data.frame(cpd = unique(data_in$cpd_id), data_for_this_cpd=NA , 
                       median_t2d = NA, median_normal = NA,  
                       alt = NA,
                       p_val = NA, W_statistic = NA, hl_effect_wilcox = NA, hl_effect_hlfxn = NA,
                       trend_with_disease = NA, 
                       
                       incResto_incT2D = NA, decResto_decT2D = NA, # same direction
                       incResto_decT2D = NA, decResto_incT2D = NA, # opposing direction
                       incResto_notT2D = NA, decResto_notT2D = NA,  # no trend in T2D
                       notResto_incT2D = NA, notResto_decT2D = NA, notResto_notT2D = NA, # not significant in resto
                       anyResto_incT2D = NA, anyResto_decT2D = NA, anyResto_notT2D = NA # any/all ecosystem scenarios
)

# define fxn: Hodges-Lehmann estimator of effect size - https://aakinshin.net/posts/r-hodges-lehmann-problems/
hl <- function(x, y = NULL) {
  if (is.null(y)) {
    walsh <- outer(x, x, "+") / 2
    median(walsh[lower.tri(walsh, diag = TRUE)])
  } else {
    median(outer(x, y, "-"))
  }
}



for (i in 1:dim(dat.test)[1]) {
  #i<-100 i<-1
  this_cpd <- dat.test$cpd[i]
  sel <- which(data_in$cpd_id == this_cpd)
  # prepare data
  df = data.frame(group = as.character( data_in$group[sel]), 
                  value = data_in$log10_abun[sel], 
                  value_perc = data_in$cpd_rel_abun[sel] )
  
  x <- df$value[df$group == "T2D met neg"]
  y <- df$value[df$group == "Normal"]
  
  if ( length(which(x == min(x) )) > 0.5*length(x) & length(which(y == min(y) )) > 0.5*length(y) ) {
    # low data: if at least one dataset does not have at least 50% non-zero cases
    dat.test$data_for_this_cpd[i] <- "low data"
    
  } else {
    
    # # test for homogeneity of variances
    # var.test(x,y) # e.g., j<-100: F = 2.3477, num df = 32, denom df = 42, p-value = 0.009872
    
    dat.test$median_t2d[i] <- median( x ) # log10 abun  
    dat.test$median_normal[i] <- median( y ) # log10 abun
    
    alt <- NA
    if (median( x ) > median( y ) ) {
      alt <- "greater"
    } else {
      alt <- "less"
    }
    dat.test$alt[i] <- alt
    
    # Wilcoxon-Mann-Whitney Test
    wmw.test <- wilcox.test(x, y, alternative = alt, paired = FALSE, conf.int = TRUE) # based on log10 abun
    
    dat.test$p_val[i] <- wmw.test$p.value
    dat.test$W_statistic[i] <- wmw.test$statistic
    
    # https://search.r-project.org/CRAN/refmans/DescTools/html/HodgesLehmann.html
    # https://aakinshin.net/posts/r-hodges-lehmann-problems
    
    dat.test$hl_effect_wilcox[i] <- wmw.test$estimate
    dat.test$hl_effect_hlfxn[i] <- hl(x, y)
    
    if (!(is.na(dat.test$p_val[i])|is.na(dat.test$W_statistic[i]))) {
      if (dat.test$p_val[i] <= 0.05 & alt == "greater") { dat.test$trend_with_disease[i] <- "Increasing" }
      else if (dat.test$p_val[i] <= 0.05 & alt == "less") { dat.test$trend_with_disease[i] <- "Decreasing" }
      else { dat.test$trend_with_disease[i] <- "No trend" }
    }
    
    # detect if same or opposing pattern to restoration?
    
    # trend with restoration?
    sel.resto <- which( dat.resto$cpd == this_cpd ) # which( dat.resto$cpd == "blahblah" )
    #dat.test.sig.resto[sel.resto, ]
    
    if (length(sel.resto) == 0) {
      
      # not significant in restoration
      if ( dat.test$trend_with_disease[i] == "Increasing") { dat.test$notResto_incT2D[i] <- 1 }
      if ( dat.test$trend_with_disease[i] == "Decreasing") { dat.test$notResto_decT2D[i] <- 1 }
      if ( dat.test$trend_with_disease[i] == "No trend") { dat.test$notResto_notT2D[i] <- 1 }
      
    } else {
      
      # same direction
      if ( dat.resto$trend_with_age[sel.resto] == "Increasing" & dat.test$trend_with_disease[i] == "Increasing") { dat.test$incResto_incT2D[i] <- 1 }
      if ( dat.resto$trend_with_age[sel.resto] == "Decreasing" & dat.test$trend_with_disease[i] == "Decreasing") { dat.test$decResto_decT2D[i] <- 1 }
      
      # opposing direction
      if ( dat.resto$trend_with_age[sel.resto] == "Increasing" & dat.test$trend_with_disease[i] == "Decreasing") { dat.test$incResto_decT2D[i] <- 1 }
      if ( dat.resto$trend_with_age[sel.resto] == "Decreasing" & dat.test$trend_with_disease[i] == "Increasing") { dat.test$decResto_incT2D[i] <- 1 }
      
      # no trend with T2D
      if ( dat.resto$trend_with_age[sel.resto] == "Increasing" & dat.test$trend_with_disease[i] == "No trend") { dat.test$incResto_notT2D[i] <- 1 }
      if ( dat.resto$trend_with_age[sel.resto] == "Decreasing" & dat.test$trend_with_disease[i] == "No trend") { dat.test$decResto_notT2D[i] <- 1 }
      
    }
    
    # now capture any/all ecosystem scenarios
    if ( dat.test$trend_with_disease[i] == "Increasing") { dat.test$anyResto_incT2D[i] <- 1 }
    if ( dat.test$trend_with_disease[i] == "Decreasing") { dat.test$anyResto_decT2D[i] <- 1 }
    if ( dat.test$trend_with_disease[i] == "No trend") { dat.test$anyResto_notT2D[i] <- 1 }
  
  }
  print(paste0("Completed i ",i))
  
} # END i for each compound



identical(dat.test$hl_effect_wilcox, dat.test$hl_effect_hlfxn)
plot(dat.test$hl_effect_wilcox, dat.test$hl_effect_hlfxn)
summary(dat.test$hl_effect_wilcox)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -0.8590 -0.0184 -0.0001 -0.0153  0.0222  0.1422     398 
summary(dat.test$hl_effect_hlfxn)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -0.8590 -0.0184  0.0000 -0.0152  0.0222  0.1422     398 

plot(dat.test$hl_effect_hlfxn, dat.test$p_val)
plot(dat.test$W_statistic, dat.test$p_val)

  


sel.low <- which(dat.test$data_for_this_cpd == "low data") # 398

length(unique(dat.cpd.collate.T2DNORM$cpd_id)) # 2052
#length(unique(dat.cpd.collate.T2DNORM$cpd_id)) - (length(sel.disq) + length(sel.low) ) # 1601
length(unique(dat.cpd.collate.T2DNORM$cpd_id)) - length(sel.low) # 1654

sel.nonNA <- which(!is.na(dat.test$p_val)) # 1654 applicable tests




# what are most interesting compounds - that decrease in restoration (increase with land disturbance) and increase with T2D ??
filter( dat.test[ order(dat.test$W_statistic, decreasing = TRUE),  ], decResto_incT2D == 1 )
#         cpd data_for_this_cpd   median_t2d median_normal     alt       p_val W_statistic hl_effect_wilcox hl_effect_hlfxn trend_with_disease incResto_incT2D decResto_decT2D
# 1  cpd02113              <NA> -1.594876896   -1.63381389 greater 0.005481583       951.0       0.04517931      0.04517931         Increasing              NA              NA
# 2  cpd22159              <NA> -2.096515544   -2.15593158 greater 0.005481583       951.0       0.05017562      0.05017562         Increasing              NA              NA
# 3  cpd03198              <NA> -1.143692673   -1.19725370 greater 0.006200811       947.0       0.05422913      0.05422913         Increasing              NA              NA
# 4  cpd33296              <NA> -2.030155050   -2.06565540 greater 0.006393020       946.0       0.03518134      0.03518134         Increasing              NA              NA
# 5  cpd35273              <NA> -2.030155050   -2.06565540 greater 0.006393020       946.0       0.03518134      0.03518134         Increasing              NA              NA
# 6  cpd01982              <NA> -1.275343989   -1.30053583 greater 0.008370273       937.0       0.02897927      0.02897927         Increasing              NA              NA
# 7  cpd01311              <NA> -2.002251707   -2.07290279 greater 0.008619683       936.0       0.05873564      0.05873564         Increasing              NA              NA
# 8  cpd00224              <NA> -0.872160949   -0.90094195 greater 0.011169028       927.0       0.04407164      0.04407164         Increasing              NA              NA
# 9  cpd00082              <NA> -0.815764168   -0.84588333 greater 0.012847012       922.0       0.05325226      0.05325226         Increasing              NA              NA
# 10 cpd00033              <NA> -0.749923397   -0.76530330 greater 0.013207264       921.0       0.01839419      0.01839419         Increasing              NA              NA
# 11 cpd00023              <NA> -0.003681729   -0.04444413 greater 0.015554592       915.0       0.02496598      0.02496598         Increasing              NA              NA
# 12 cpd01693              <NA> -1.202284109   -1.24849363 greater 0.015554592       915.0       0.06639789      0.06639789         Increasing              NA              NA
# 13 cpd00069              <NA> -1.113884830   -1.15144646 greater 0.015978303       914.0       0.02643750      0.02643750         Increasing              NA              NA
# 14 cpd02394              <NA> -1.580316876   -1.62033253 greater 0.018730844       908.0       0.03361335      0.03361335         Increasing              NA              NA
# 15 cpd29078              <NA> -3.221567930   -3.27280746 greater 0.020778970       904.0       0.05228592      0.05228592         Increasing              NA              NA
# 16 cpd12672              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 17 cpd12673              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 18 cpd12674              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 19 cpd12697              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 20 cpd27089              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 21 cpd27073              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 22 cpd27088              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 23 cpd29537              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 24 cpd29536              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 25 cpd29314              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 26 cpd29315              <NA> -2.442486889   -2.48602625 greater 0.021319283       903.0       0.04235738      0.04235738         Increasing              NA              NA
# 27 cpd00259              <NA> -1.583655418   -1.63772883 greater 0.021871276       902.0       0.04871078      0.04871078         Increasing              NA              NA
# 28 cpd00024              <NA> -0.326522472   -0.36300871 greater 0.022435134       901.0       0.02181594      0.02181594         Increasing              NA              NA
# 29 cpd00282              <NA> -1.210866794   -1.29050987 greater 0.022435134       901.0       0.03795462      0.03795462         Increasing              NA              NA
# 30 cpd00288              <NA> -0.774948336   -0.81174799 greater 0.023011044       900.0       0.02914844      0.02914844         Increasing              NA              NA
# 31 cpd00242              <NA> -1.059282213   -1.07304815 greater 0.025438984       896.0       0.02108114      0.02108114         Increasing              NA              NA
# 32 cpd15378              <NA> -2.176917237   -2.18854997 greater 0.026730211       894.0       0.05826438      0.05826438         Increasing              NA              NA
# 33 cpd02678              <NA> -1.343401382   -1.36763131 greater 0.026730211       894.0       0.02788244      0.02788244         Increasing              NA              NA
# 34 cpd09846              <NA> -3.155894728   -3.23496612 greater 0.027395815       893.0       0.06703895      0.06703895         Increasing              NA              NA
# 35 cpd09847              <NA> -2.854864732   -2.93393612 greater 0.027395815       893.0       0.06703895      0.06703895         Increasing              NA              NA
# 36 cpd09848              <NA> -2.854864732   -2.93393612 greater 0.027395815       893.0       0.06703895      0.06703895         Increasing              NA              NA
# 37 cpd09849              <NA> -3.155894728   -3.23496612 greater 0.027395815       893.0       0.06703895      0.06703895         Increasing              NA              NA
# 38 cpd15891              <NA> -3.155894728   -3.23496612 greater 0.027395815       893.0       0.06703895      0.06703895         Increasing              NA              NA
# 39 cpd28546              <NA> -2.854864732   -2.93393612 greater 0.027395815       893.0       0.06703895      0.06703895         Increasing              NA              NA
# 40 cpd28547              <NA> -2.854864732   -2.93393612 greater 0.027395815       893.0       0.06703895      0.06703895         Increasing              NA              NA
# 41 cpd00424              <NA> -1.593607732   -1.64121133 greater 0.028075005       892.0       0.03910868      0.03910868         Increasing              NA              NA
# 42 cpd12543              <NA> -1.632646402   -1.66880349 greater 0.028767979       891.0       0.04271343      0.04271343         Increasing              NA              NA
# 43 cpd12848              <NA> -1.627038351   -1.66420678 greater 0.029474932       890.0       0.04259876      0.04259876         Increasing              NA              NA
# 44 cpd00247              <NA> -1.332698079   -1.39180468 greater 0.030931571       888.0       0.03827772      0.03827772         Increasing              NA              NA
# 45 cpd00076              <NA> -1.182489505   -1.25902186 greater 0.037357168       880.0       0.05947599      0.05947599         Increasing              NA              NA
# 46 cpd11597              <NA> -2.484312496   -2.55269248 greater 0.038709442       878.5       0.07087141      0.07085951         Increasing              NA              NA
# 47 cpd02654              <NA> -1.831004539   -1.88680131 greater 0.040028258       877.0       0.04972258      0.04972258         Increasing              NA              NA
# 48 cpd24592              <NA> -3.241517166   -3.28920284 greater 0.040951934       876.0       0.04116872      0.04116872         Increasing              NA              NA
# 49 cpd23593              <NA> -2.673105520   -2.71277855 greater 0.040951934       876.0       0.03933646      0.03933646         Increasing              NA              NA
# 50 cpd22517              <NA> -2.639821026   -2.67653652 greater 0.043825735       873.0       0.04033845      0.04033845         Increasing              NA              NA
# 51 cpd02820              <NA> -1.539829824   -1.58358312 greater 0.044818599       872.0       0.04280860      0.04280860         Increasing              NA              NA
# 52 cpd03200              <NA> -1.578503413   -1.63624654 greater 0.045829264       871.0       0.04789511      0.04789511         Increasing              NA              NA
# 53 cpd00773              <NA> -1.530542798   -1.57256837 greater 0.047904802       869.0       0.03120599      0.03120599         Increasing              NA              NA
# 54 cpd00193              <NA> -1.413597330   -1.46003043 greater 0.047904802       869.0       0.04088073      0.04088073         Increasing              NA              NA
# 55 cpd00522              <NA> -1.105296529   -1.14682308 greater 0.047904802       869.0       0.03525662      0.03525662         Increasing              NA              NA

write.csv(x = filter( dat.test[ order(dat.test$W_statistic, decreasing = TRUE),  ], decResto_incT2D == 1 ), file = "decResto_incT2D--dat.test-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.csv")

saveRDS(dat.test, file = "dat.test-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.RDS")

# only keep applicable tests; the remainder are likely due to zero replacement:
#dat.test <- dat.test[ sel.nonNA, ]

sel.sig <- which(dat.test$p_val <= 0.05) # 276

summary(dat.test$p_val)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0008  0.0922  0.2080  0.2431  0.3538  0.8831     398

hist(dat.test$p_val)

names(dat.test)


sum(dat.test[ ,c("incResto_incT2D", 
                 "decResto_decT2D", 
                 "incResto_decT2D", 
                 "decResto_incT2D", 
                 "incResto_notT2D", 
                 "decResto_notT2D", 
                 "notResto_incT2D", 
                 "notResto_decT2D")], na.rm = TRUE) # 1654

sum(dat.test[ ,c("incResto_incT2D", 
                 "decResto_decT2D", 
                 "incResto_decT2D", 
                 "decResto_incT2D", 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 "notResto_incT2D", 
                 "notResto_decT2D")], na.rm = TRUE) # 276

sum(dat.test[ ,c(#"incResto_incT2D", 
                 #"decResto_decT2D", 
                 #"incResto_decT2D", 
                 #"decResto_incT2D", 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 "notResto_incT2D", 
                 "notResto_decT2D"
                 )], na.rm = TRUE) # 0

sum(dat.test[ ,c("incResto_incT2D", 
                 "decResto_decT2D", 
                 "incResto_decT2D", 
                 "decResto_incT2D" #, 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 #"notResto_incT2D", 
                 #"notResto_decT2D"
                 )], na.rm = TRUE) # 276

sum(dat.test[ ,c("incResto_incT2D", 
                 #"decResto_decT2D", 
                 #"incResto_decT2D", 
                 "decResto_incT2D" #, 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 #"notResto_incT2D", 
                 #"notResto_decT2D"
)], na.rm = TRUE) # 108

sum(dat.test[ ,c(#"incResto_incT2D", 
                 "decResto_decT2D", 
                 "incResto_decT2D" #, 
                 #"decResto_incT2D" #, 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 #"notResto_incT2D", 
                 #"notResto_decT2D"
)], na.rm = TRUE) # 268




# extract sig results

#sel.sig <- which(dat.test$sigBH == "sig") # 60
sel.sig <- which(dat.test$p_val <= 0.05) # 276

dat.test.sig <- dat.test[sel.sig, ]

dat.test.sig$minuslog10_p_val <- -log10(dat.test.sig$p_val)

plot(x = dat.test.sig$diff_mean_perc , y =dat.test.sig$minuslog10_p_val , xlab="Difference (functional %)", ylab="-log10(P-value)")

plot(x = dat.test.sig$W_statistic , y =dat.test.sig$minuslog10_p_val , xlab="W statistic", ylab="-log10(P-value)")


dev.print(tiff, filename = paste0(workdir,"/plots/","3d-indiv-compound-VolcanoPlot-P-values--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"),
          width = 12, height = 14, units = "cm", res=600, compression = "lzw",type="cairo" )


# join compound info??

#length(which(dat.test.sig$sigBH=="sig")) # 2122 1953


dat.test.sig$cpd_names <- NA
dat.test.sig$cpd_forms <- NA

dat.test.sig$OC_x <- NA
dat.test.sig$HC_y <- NA
dat.test.sig$NC_z <- NA

dat.test.sig$mass <- NA

dat.test.sig$trend_group <- NA


for (i in 1:dim(dat.test.sig)[1]) {
  #i<-1
  this_cpd <- dat.test.sig$cpd[i]
  
  sel.cpd <- which(df.comp$id == this_cpd)
  
  # dat.test.sig$cpd_names[i] <- paste0( unique( dat.cpd.distil$cpd_names[sel.cpd] ), collapse = "||")
  # dat.test.sig$cpd_forms[i] <- paste0( unique( dat.cpd.distil$cpd_forms[sel.cpd] ), collapse = "||")
  
  dat.test.sig$cpd_names[i] <- df.comp$name[sel.cpd]
  dat.test.sig$cpd_forms[i] <- df.comp$form[sel.cpd]
  
  dat.test.sig$OC_x[i] <- df.comp$OC_ratio[sel.cpd]
  dat.test.sig$HC_y[i] <- df.comp$HC_ratio[sel.cpd]
  dat.test.sig$NC_z[i] <- df.comp$NC_ratio[sel.cpd]
  
  # sel.lut <- which(compounds.lut$id == this_cpd)
  # dat.test.sig$mass[i] <- compounds.lut$mass[sel.lut]
  
  dat.test.sig$mass[i] <- df.comp$mass[sel.cpd]
  
  if (!is.na(dat.test.sig$incResto_incT2D[i]) & dat.test.sig$incResto_incT2D[i] == 1) {dat.test.sig$trend_group[i] <- "Increasing in T2D (increased exposure in quality ecosystems)" } # (associated with ecosystem quality)
  if (!is.na(dat.test.sig$decResto_decT2D[i]) & dat.test.sig$decResto_decT2D[i] == 1) {dat.test.sig$trend_group[i] <- "Decreasing in T2D (reduced exposure in quality ecosystems)" } # (associated with ecosystem disturbance)
  if (!is.na(dat.test.sig$incResto_decT2D[i]) & dat.test.sig$incResto_decT2D[i] == 1) {dat.test.sig$trend_group[i] <- "Decreasing in T2D (reduced exposure in disturbed ecosystems)" } # (associated with ecosystem quality)
  if (!is.na(dat.test.sig$decResto_incT2D[i]) & dat.test.sig$decResto_incT2D[i] == 1) {dat.test.sig$trend_group[i] <- "Increasing in T2D (increased exposure in disturbed ecosystems)" } # (associated with ecosystem disturbance)
  
  print(paste0("completed ",i))
}

write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.rds")

#dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.rds")




## plot as Increasing or Decreasing?? in vK space


names(dat.test.sig)
# [1] "cpd"                "data_for_this_cpd"  "p_val"              "median_t2d"         "median_Normal"      "W_statistic"        "diff_median_perc"  
# [8] "diff_mean_perc"     "alt"                "trend_with_disease" "incResto_incT2D"    "decResto_decT2D"    "incResto_decT2D"    "decResto_incT2D"   
# [15] "incResto_notT2D"    "decResto_notT2D"    "minuslog10_p_val"   "cpd_names"          "cpd_forms"          "OC_x"               "HC_y"              
# [22] "NC_z"               "mass"               "trend_group" 

unique( dat.test.sig$trend_group )
# [1] "Decreasing in T2D (reduced exposure in quality ecosystems)"   "Increasing in T2D (increased exposure in disturbed ecosystems)"
# [3] "Decreasing in T2D (reduced exposure in disturbed ecosystems)"   "Increasing in T2D (increased exposure in quality ecosystems)" 

## superceded
# [1] "Decreasing in T2D (associated with ecosystem disturbance)" "Increasing in T2D (associated with ecosystem disturbance)"
# [3] "Decreasing in T2D (associated with ecosystem quality)"     "Increasing in T2D (associated with ecosystem quality)" 


dat.test.sig$trend_group <- factor(dat.test.sig$trend_group,
                                   levels = c("Decreasing in T2D (reduced exposure in quality ecosystems)",
                                              "Decreasing in T2D (reduced exposure in disturbed ecosystems)",
                                              "Increasing in T2D (increased exposure in disturbed ecosystems)",
                                              "Increasing in T2D (increased exposure in quality ecosystems)" ),
                                   ordered = TRUE)

col.trend_group <- c("Decreasing in T2D (reduced exposure in quality ecosystems)" = "#fdb462", 
                     "Decreasing in T2D (reduced exposure in disturbed ecosystems)" = "#e31a1c",
                     "Increasing in T2D (increased exposure in disturbed ecosystems)" = "#33a02c" ,
                     "Increasing in T2D (increased exposure in quality ecosystems)" = "#1f78b4" )

p <- #ggplot(data = filter(dat.test.sig, sigBH == "sig")) +
  ggplot(data = dat.test.sig) +
  coord_equal()+
  ggtitle("Microbiota compound processing potential\n- Type 2 Diabetes case study\n(Restoration-significant compounds)")+
  #xlim(0,3.4)+ ylim(0,4.1)+
  
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional\ncapacity (%)\nallocated to\ncompounds", nrow = 4))+
  
  geom_mark_rect(data= vkgrouprect, aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ 
  
  annotate(geom="text", x= 0+0.01, y= 2.3+0.02, label = "Lipid", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.2+0.01, y= 2.2+0.02, label = "Protein", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  #annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373", lineheight = 0.8) + # top-left
  annotate(geom="text", x= 0.7+0.01, y= 2.4+0.02, label = "Carbohydrate", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0+0.01, y= 0.5-0.02, label = "Condensed\naromatics", hjust=0, vjust=1, size = 2 , col="#737373", lineheight = 0.8) + # bottom-left
  annotate(geom="text", x= 0.25+0.01, y= 0.75-0.02, label = "Lignin", hjust=0, vjust=1, size = 2 , col="#737373" ) + # bottom-left
  annotate(geom="text", x= 0.67+0.01, y= 0.53-0.02, label = "Tannin", hjust=0, vjust=1, size = 2 , col="#737373") + # bottom-left
  
  theme_bw()+
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM-Wilcox",this_study,header,".tiff"), width = 14, height = 14, units = "cm", res=350, compression="lzw",type="cairo")
# Removed 147 rows containing missing values or values outside the scale range (`geom_point()`). 


hist(dat.test.sig$NC_z)


min(dat.test.sig$NC_z[ dat.test.sig$NC_z > 0 ], na.rm = TRUE) # 0.01449275


dim(dat.test.sig) # 276 24
sel.ok <- which(!is.na(dat.test.sig$NC_z) ) # qty 260

dat.test.sig$z_layer <- NA

subsel <- which(dat.test.sig$NC_z[sel.ok] == 0) # 91
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C = 0"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0 & dat.test.sig$NC_z[sel.ok] <= 0.2 ) # 100
# Rivas-Ubach et al 2018 Table 1 highlights N/C breaks at 0.126, 0.2, 0.5, 0.7
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0 to 0.2"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0.2) # 69
max(dat.test.sig$NC_z[sel.ok]) # 1
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 1"

unique(dat.test.sig$z_layer[sel.ok]) # "N:C = 0"       "N:C >0 to 0.2" "N:C >0.2 to 1"

dat.test.sig$z_layer <- factor(dat.test.sig$z_layer, levels = c("N:C = 0",
                                                                "N:C >0 to 0.2",
                                                                "N:C >0.2 to 1"), ordered = TRUE)

saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.rds")


dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.rds")

dim(dat.test.sig) #  276  25
head(dat.test.sig)
dim(dat.test.sig[ which(dat.test.sig$p_val < 0.05), ]) # 276  25
sel <- which(!is.na(dat.test.sig$OC_x) ) # 260
unique(dat.test.sig$trend_group)
# [1] Decreasing in T2D (reduced exposure in quality ecosystems)    
# [2] Increasing in T2D (increased exposure in disturbed ecosystems)
# [3] Decreasing in T2D (reduced exposure in disturbed ecosystems)  
# [4] Increasing in T2D (increased exposure in quality ecosystems)  
dat.test.sig$trend_group <- factor( dat.test.sig$trend_group,
                                    levels = c(
                                      "Decreasing in T2D (reduced exposure in disturbed ecosystems)",
                                      "Increasing in T2D (increased exposure in disturbed ecosystems)",
                                      "Decreasing in T2D (reduced exposure in quality ecosystems)",
                                      "Increasing in T2D (increased exposure in quality ecosystems)"),
                                    ordered = TRUE)

sel <- which(dat.test.sig$trend_group == "Decreasing in T2D (reduced exposure in disturbed ecosystems)") # 70
sel <- which(dat.test.sig$trend_group == "Decreasing in T2D (reduced exposure in disturbed ecosystems)" & !is.na(dat.test.sig$OC_x) ) # 66
sel <- which(dat.test.sig$trend_group == "Increasing in T2D (increased exposure in disturbed ecosystems)") # 58
sel <- which(dat.test.sig$trend_group == "Increasing in T2D (increased exposure in disturbed ecosystems)" & !is.na(dat.test.sig$OC_x) ) # 52

sel <- which(dat.test.sig$trend_group == "Decreasing in T2D (reduced exposure in quality ecosystems)") # 98
sel <- which(dat.test.sig$trend_group == "Decreasing in T2D (reduced exposure in quality ecosystems)" & !is.na(dat.test.sig$OC_x) ) # 93
sel <- which(dat.test.sig$trend_group == "Increasing in T2D (increased exposure in quality ecosystems)") # 50
sel <- which(dat.test.sig$trend_group == "Increasing in T2D (increased exposure in quality ecosystems)" & !is.na(dat.test.sig$OC_x) ) # 49




# compound class zones & labels - loaded earlier
# vkgrouprect.facets2 <- read.table(file = "cpp3d-compound-classes.tsv", header = TRUE, sep = "\t" )
# vkgrouprect.facets2.labels <- read.table(file = "cpp3d-compound-classes-labels.tsv", header = TRUE, sep = "\t" )



p <- ggplot(data = dat.test.sig[sel.ok, ]) +
  coord_equal()+
  
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(considering only compounds with significant CPP variation in ecosystem restoration)")+
  #ggtitle("Compounds with trending CPP in T2D and soil-ecosystem quality)")+
  
  #xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  
  xlab("O:C ratio")+ ylab("H:C ratio")+
  
  facet_wrap(facets = vars(z_layer))+
  
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  #geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 )+ # 
  #guides(color = guide_legend(title = "Trend with disease in functional capacity\n(%) allocated to compounds")) +
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  #guides(color = guide_legend(title = "Trend with disease\nin functional capacity(%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", ncol = 1 )) +
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T,  size = 2 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2 , col="#737373" , lineheight = 0.8)+
  #geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2 , col="#737373" , lineheight = 0.8)+
  
  theme_bw()+
  theme(
    
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
    
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-large.tiff"), width = 20, height = 12, units = "cm", res=350, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-lim.tiff"), width = 20, height = 12, units = "cm", res=350, compression="lzw",type="cairo")

dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-TrendGroups-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-large.tiff"), width = 20, height = 13, units = "cm", res=350, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-TrendGroups-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-lim.tiff"), width = 20, height = 13, units = "cm", res=350, compression="lzw",type="cairo")

dev.print(tiff, file = paste0(workdir,"/plots/","noGgtitle-3d-Compounds-indiv-vKSpace-TrendGroups-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-lim.tiff"), width = 20, height = 13, units = "cm", res=600, compression="lzw",type="cairo")


## ZOOM-IN
# ZOOM1

p <- ggplot(data = filter( dat.test.sig[sel.ok, ], z_layer == "N:C >0 to 0.2" ) ) +
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.2, 0.6)+ ylim(1.5,2.0)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  

  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM1-3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")


# ZOOM2

p <- ggplot(data = filter( dat.test.sig[sel.ok, ], z_layer == "N:C = 0" ) ) +
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.8, 1.1)+ ylim(1.6,2.1)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
 
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM2-3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------


#### ZOOM plots and tests on example cpp3d regions?
#    Zoom1 - decreasing proteins
#    Zoom2 - increasing carbs
#    Etc - 4 groups
#-------------------------

#saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.rds")

dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox.rds")

sel.ok <- which(!is.na(dat.test.sig$NC_z) ) # qty 260


# build vk layer for each z-layer facet AND age group
for (i in 1:length(levels(dat.test.sig$z_layer))) {
  #i<-1
  temp <- vkgrouprect
  temp$z_layer <- levels(dat.test.sig$z_layer)[i]
  if (i==1) { keep <- temp }
  if (i>1) { keep <- rbind(keep, temp)}
  print(paste0("completed ",i))
}
vkgrouprect.facets <- keep
rm(keep)

str( vkgrouprect.facets )
vkgrouprect.facets$z_layer <- factor( vkgrouprect.facets$z_layer )


vkgrouprect.facets


p <- ggplot(data = dat.test.sig[sel.ok, ]) +
  coord_equal()+
  ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  #xlim(0,3.4)+ ylim(0,4.1)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  facet_wrap(facets = vars(z_layer))+
  
  #geom_mark_rect(data= vkgrouprect, aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ 
  
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  annotate(geom="text", x= 0+0.01, y= 2.3+0.02, label = "Lipid", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.2+0.01, y= 2.2+0.02, label = "Protein", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  #annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373", lineheight = 0.8) + # top-left
  annotate(geom="text", x= 0.7+0.01, y= 2.4+0.02, label = "Carbohydrate", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0+0.01, y= 0.5-0.02, label = "Condensed\naromatics", hjust=0, vjust=1, size = 2 , col="#737373", lineheight = 0.8) + # bottom-left
  annotate(geom="text", x= 0.25+0.01, y= 0.75-0.02, label = "Lignin", hjust=0, vjust=1, size = 2 , col="#737373" ) + # bottom-left
  annotate(geom="text", x= 0.67+0.01, y= 0.53-0.02, label = "Tannin", hjust=0, vjust=1, size = 2 , col="#737373") + # bottom-left
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"), width = 20, height = 12, units = "cm", res=350, compression="lzw",type="cairo")




## ZOOM-IN - Proteins? N:C >0 to 0.2
# ZOOM1

p <- ggplot(data = filter( dat.test.sig[sel.ok[-subsel.diab.prot2], ], z_layer == "N:C >0 to 0.2" ) ) + # sel.ok[subsel.diab.prot1]
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.2, 0.6)+ ylim(1.5,2.0)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 2, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  annotate(geom = "label_npc", npcx = 0.02, npcy = 0.97, label = "N:C >0 to 0.2", size = 3.5)+ # geom = "text_npc"

theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM1-3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")


# check "Proteins" line1 & line2 ??

# examples
sel.diab.prot1 <- which(dat.test.sig$cpd %in% c("cpd11511","cpd11528","cpd11524","cpd11549"))
dat.test.sig[sel.diab.prot1, ]
# cpd data_for_this_cpd      p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc  alt trend_with_disease
# 1647 cpd11511              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing
# 1651 cpd11524              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing
# 1653 cpd11528              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing
# 1661 cpd11549              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing
# incResto_incT2D decResto_decT2D incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                cpd_names
# 1647              NA              NA               1              NA              NA              NA         1.688342 10-methyl-dodecanoyl-ACP
# 1651              NA              NA               1              NA              NA              NA         1.688342    5-methyl-hexanoyl-ACP
# 1653              NA              NA               1              NA              NA              NA         1.688342    7-methyl-octanoyl-ACP
# 1661              NA              NA               1              NA              NA              NA         1.688342   4-methyl-pentanoyl-ACP
# cpd_forms      OC_x     HC_y       NC_z mass                                                  trend_group       z_layer
# 1647 C24H45N2O8PRS 0.3333333 1.875000 0.08333333  553 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1651 C18H33N2O8PRS 0.4444444 1.833333 0.11111111  469 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1653 C20H37N2O8PRS 0.4000000 1.850000 0.10000000  497 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1661 C17H31N2O8PRS 0.4705882 1.823529 0.11764706  455 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
dat.test.sig$HC_y[sel.diab.prot1]/dat.test.sig$OC_x[sel.diab.prot1]
(1.833333 - 1.875000)/(0.4444444 - 0.3333333) # -0.375

# estimate of line: y = -0.375x + 2

sel.diab.prot1 <- which(dat.test.sig$HC_y == (-0.375*dat.test.sig$OC_x + 2) & dat.test.sig$z_layer == "N:C >0 to 0.2" & dat.test.sig$trend_with_disease == "Decreasing")
# qty 10
# use rounding to 4 d.p
subsel.diab.prot1 <- which(round(dat.test.sig$HC_y[sel.ok], 4) == round(-0.375*dat.test.sig$OC_x[sel.ok] + 2, 4) & dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$trend_with_disease[sel.ok] == "Decreasing")
length(subsel.diab.prot1) # 12

dat.test.sig[sel.ok[subsel.diab.prot1], ]
# cpd data_for_this_cpd      p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc  alt trend_with_disease incResto_incT2D decResto_decT2D
# 1641 cpd11499              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1643 cpd11503              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1645 cpd11507              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1647 cpd11511              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1651 cpd11524              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1653 cpd11528              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1655 cpd11532              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1657 cpd11536              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1661 cpd11549              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1663 cpd11553              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1665 cpd11557              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1667 cpd11561              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                cpd_names     cpd_forms      OC_x     HC_y       NC_z mass
# 1641               1              NA              NA              NA         1.688342    4-methyl-hexanoyl-ACP C18H33N2O8PRS 0.4444444 1.833333 0.11111111  469
# 1643               1              NA              NA              NA         1.688342    6-methyl-octanoyl-ACP C20H37N2O8PRS 0.4000000 1.850000 0.10000000  497
# 1645               1              NA              NA              NA         1.688342    8-methyl-decanoyl-ACP C22H41N2O8PRS 0.3636364 1.863636 0.09090909  525
# 1647               1              NA              NA              NA         1.688342 10-methyl-dodecanoyl-ACP C24H45N2O8PRS 0.3333333 1.875000 0.08333333  553
# 1651               1              NA              NA              NA         1.688342    5-methyl-hexanoyl-ACP C18H33N2O8PRS 0.4444444 1.833333 0.11111111  469
# 1653               1              NA              NA              NA         1.688342    7-methyl-octanoyl-ACP C20H37N2O8PRS 0.4000000 1.850000 0.10000000  497
# 1655               1              NA              NA              NA         1.688342    9-methyl-decanoyl-ACP C22H41N2O8PRS 0.3636364 1.863636 0.09090909  525
# 1657               1              NA              NA              NA         1.688342 11-methyl-dodecanoyl-ACP C24H45N2O8PRS 0.3333333 1.875000 0.08333333  553
# 1661               1              NA              NA              NA         1.688342   4-methyl-pentanoyl-ACP C17H31N2O8PRS 0.4705882 1.823529 0.11764706  455
# 1663               1              NA              NA              NA         1.688342   6-methyl-heptanoyl-ACP C19H35N2O8PRS 0.4210526 1.842105 0.10526316  483
# 1665               1              NA              NA              NA         1.688342    8-methyl-nonanoyl-ACP C21H39N2O8PRS 0.3809524 1.857143 0.09523810  511
# 1667               1              NA              NA              NA         1.688342 10-methyl-undecanoyl-ACP C23H43N2O8PRS 0.3478261 1.869565 0.08695652  539
# trend_group       z_layer
# 1641 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1643 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1645 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1647 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1651 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1653 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1655 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1657 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1661 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1663 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1665 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1667 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2

sel.diab.prot2 <- which(dat.test.sig$cpd %in% c("cpd11572","cpd11531","cpd11473","cpd11548"))
dat.test.sig[sel.diab.prot2, ]
# cpd data_for_this_cpd      p_val kendall_tau trend_with_age minuslog10_p_val                       cpd_names
# 5724 cpd11473              <NA> 0.04351584   -2.018725     Decreasing         1.361353             (2E)-Hexenoyl-[acp]
# 5742 cpd11531              <NA> 0.04351584   -2.018725     Decreasing         1.361353  9-methyl-trans-dec-2-enoyl-ACP
# 5750 cpd11548              <NA> 0.04351584   -2.018725     Decreasing         1.361353 4-methyl-trans-pent-2-enoyl-ACP
# 5762 cpd11572              <NA> 0.04351584   -2.018725     Decreasing         1.361353       trans-Octodec-2-enoyl-ACP
# cpd_forms      OC_x     HC_y       NC_z mass       z_layer
# 5724 C17H29N2O8PRS 0.4705882 1.705882 0.11764706  453 N:C >0 to 0.2
# 5742 C22H39N2O8PRS 0.3636364 1.772727 0.09090909  523 N:C >0 to 0.2
# 5750 C17H29N2O8PRS 0.4705882 1.705882 0.11764706  453 N:C >0 to 0.2
# 5762 C29H53N2O8PRS 0.2758621 1.827586 0.06896552  621 N:C >0 to 0.2

# estimate of line: y = -0.625x + 2.0
# round to 4 d.p.
sel.diab.prot2 <- which(round(dat.test.sig$HC_y,4) == round(-0.625*dat.test.sig$OC_x + 2,4) & dat.test.sig$z_layer == "N:C >0 to 0.2" & dat.test.sig$trend_with_age == "Decreasing")
dat.test.sig[sel.diab.prot2, ]
# cpd data_for_this_cpd      p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc  alt trend_with_disease incResto_incT2D decResto_decT2D
# 1639 cpd11473              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1654 cpd11531              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1660 cpd11548              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1670 cpd11572              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                       cpd_names     cpd_forms      OC_x     HC_y       NC_z mass
# 1639               1              NA              NA              NA         1.688342             (2E)-Hexenoyl-[acp] C17H29N2O8PRS 0.4705882 1.705882 0.11764706  453
# 1654               1              NA              NA              NA         1.688342  9-methyl-trans-dec-2-enoyl-ACP C22H39N2O8PRS 0.3636364 1.772727 0.09090909  523
# 1660               1              NA              NA              NA         1.688342 4-methyl-trans-pent-2-enoyl-ACP C17H29N2O8PRS 0.4705882 1.705882 0.11764706  453
# 1670               1              NA              NA              NA         1.688342       trans-Octodec-2-enoyl-ACP C29H53N2O8PRS 0.2758621 1.827586 0.06896552  621
# trend_group       z_layer
# 1639 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1654 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1660 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1670 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2


# estimate of line 2: y = -0.625x + 2.0

subsel.diab.prot2 <- which(round(dat.test.sig$HC_y[sel.ok], 4) == round(-0.625*dat.test.sig$OC_x[sel.ok] + 2, 4) & dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$trend_with_disease[sel.ok] == "Decreasing")
length(subsel.diab.prot2) # 23


dat.test.sig[ sel.ok[unique(c(subsel.diab.prot1, subsel.diab.prot2))], ]
# cpd data_for_this_cpd      p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc  alt trend_with_disease incResto_incT2D decResto_decT2D
# 1641 cpd11499              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1643 cpd11503              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1645 cpd11507              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1647 cpd11511              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1651 cpd11524              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1653 cpd11528              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1655 cpd11532              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1657 cpd11536              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1661 cpd11549              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1663 cpd11553              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1665 cpd11557              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1667 cpd11561              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1636 cpd11475              <NA> 0.01851005  -4.252574     -4.068426         510    -2.952102e-05  -4.215527e-05 less         Decreasing              NA              NA
# 1637 cpd11465              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1638 cpd11469              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1639 cpd11473              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1640 cpd11498              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1642 cpd11502              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1644 cpd11506              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1646 cpd11510              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1648 cpd11514              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1649 cpd11518              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1650 cpd11523              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1652 cpd11527              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1654 cpd11531              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1656 cpd11535              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1658 cpd11539              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1659 cpd11543              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1660 cpd11548              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1662 cpd11552              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1664 cpd11556              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1666 cpd11560              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1668 cpd11564              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1669 cpd11568              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1670 cpd11572              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
#      incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                             cpd_names     cpd_forms      OC_x     HC_y       NC_z mass
# 1641               1              NA              NA              NA         1.688342                 4-methyl-hexanoyl-ACP C18H33N2O8PRS 0.4444444 1.833333 0.11111111  469
# 1643               1              NA              NA              NA         1.688342                 6-methyl-octanoyl-ACP C20H37N2O8PRS 0.4000000 1.850000 0.10000000  497
# 1645               1              NA              NA              NA         1.688342                 8-methyl-decanoyl-ACP C22H41N2O8PRS 0.3636364 1.863636 0.09090909  525
# 1647               1              NA              NA              NA         1.688342              10-methyl-dodecanoyl-ACP C24H45N2O8PRS 0.3333333 1.875000 0.08333333  553
# 1651               1              NA              NA              NA         1.688342                 5-methyl-hexanoyl-ACP C18H33N2O8PRS 0.4444444 1.833333 0.11111111  469
# 1653               1              NA              NA              NA         1.688342                 7-methyl-octanoyl-ACP C20H37N2O8PRS 0.4000000 1.850000 0.10000000  497
# 1655               1              NA              NA              NA         1.688342                 9-methyl-decanoyl-ACP C22H41N2O8PRS 0.3636364 1.863636 0.09090909  525
# 1657               1              NA              NA              NA         1.688342              11-methyl-dodecanoyl-ACP C24H45N2O8PRS 0.3333333 1.875000 0.08333333  553
# 1661               1              NA              NA              NA         1.688342                4-methyl-pentanoyl-ACP C17H31N2O8PRS 0.4705882 1.823529 0.11764706  455
# 1663               1              NA              NA              NA         1.688342                6-methyl-heptanoyl-ACP C19H35N2O8PRS 0.4210526 1.842105 0.10526316  483
# 1665               1              NA              NA              NA         1.688342                 8-methyl-nonanoyl-ACP C21H39N2O8PRS 0.3809524 1.857143 0.09523810  511
# 1667               1              NA              NA              NA         1.688342              10-methyl-undecanoyl-ACP C23H43N2O8PRS 0.3478261 1.869565 0.08695652  539
# 1636               1              NA              NA              NA         1.732592                   (2E)-Decenoyl-[acp] C21H37N2O8PRS 0.3809524 1.761905 0.09523810  509
# 1637               1              NA              NA              NA         1.688342    But-2-enoyl-[acyl-carrier protein] C15H25N2O8PRS 0.5333333 1.666667 0.13333333  425
# 1638               1              NA              NA              NA         1.688342                 (2E)-Dodecenoyl-[acp] C23H41N2O8PRS 0.3478261 1.782609 0.08695652  537
# 1639               1              NA              NA              NA         1.688342                   (2E)-Hexenoyl-[acp] C17H29N2O8PRS 0.4705882 1.705882 0.11764706  453
# 1640               1              NA              NA              NA         1.688342        4-methyl-trans-hex-2-enoyl-ACP C18H31N2O8PRS 0.4444444 1.722222 0.11111111  467
# 1642               1              NA              NA              NA         1.688342        6-methyl-trans-oct-2-enoyl-ACP C20H35N2O8PRS 0.4000000 1.750000 0.10000000  495
# 1644               1              NA              NA              NA         1.688342        8-methyl-trans-dec-2-enoyl-ACP C22H39N2O8PRS 0.3636364 1.772727 0.09090909  523
# 1646               1              NA              NA              NA         1.688342     10-methyl-trans-dodec-2-enoyl-ACP C24H43N2O8PRS 0.3333333 1.791667 0.08333333  551
# 1648               1              NA              NA              NA         1.688342 12-methyl-trans-tetra-dec-2-enoyl-ACP C26H47N2O8PRS 0.3076923 1.807692 0.07692308  579
# 1649               1              NA              NA              NA         1.688342  14-methyl-trans-hexa-dec-2-enoyl-ACP C28H51N2O8PRS 0.2857143 1.821429 0.07142857  607
# 1650               1              NA              NA              NA         1.688342        5-methyl-trans-hex-2-enoyl-ACP C18H31N2O8PRS 0.4444444 1.722222 0.11111111  467
# 1652               1              NA              NA              NA         1.688342        7-methyl-trans-oct-2-enoyl-ACP C20H35N2O8PRS 0.4000000 1.750000 0.10000000  495
# 1654               1              NA              NA              NA         1.688342        9-methyl-trans-dec-2-enoyl-ACP C22H39N2O8PRS 0.3636364 1.772727 0.09090909  523
# 1656               1              NA              NA              NA         1.688342     11-methyl-trans-dodec-2-enoyl-ACP C24H43N2O8PRS 0.3333333 1.791667 0.08333333  551
# 1658               1              NA              NA              NA         1.688342 13-methyl-trans-tetra-dec-2-enoyl-ACP C26H47N2O8PRS 0.3076923 1.807692 0.07692308  579
# 1659               1              NA              NA              NA         1.688342  15-methyl-trans-hexa-dec-2-enoyl-ACP C28H51N2O8PRS 0.2857143 1.821429 0.07142857  607
# 1660               1              NA              NA              NA         1.688342       4-methyl-trans-pent-2-enoyl-ACP C17H29N2O8PRS 0.4705882 1.705882 0.11764706  453
# 1662               1              NA              NA              NA         1.688342       6-methyl-trans-hept-2-enoyl-ACP C19H33N2O8PRS 0.4210526 1.736842 0.10526316  481
# 1664               1              NA              NA              NA         1.688342        8-methyl-trans-non-2-enoyl-ACP C21H37N2O8PRS 0.3809524 1.761905 0.09523810  509
# 1666               1              NA              NA              NA         1.688342     10-methyl-trans-undec-2-enoyl-ACP C23H41N2O8PRS 0.3478261 1.782609 0.08695652  537
# 1668               1              NA              NA              NA         1.688342    12-methyl-trans-tridec-2-enoyl-ACP C25H45N2O8PRS 0.3200000 1.800000 0.08000000  565
# 1669               1              NA              NA              NA         1.688342  14-methyl-trans-pentadec-2-enoyl-ACP C27H49N2O8PRS 0.2962963 1.814815 0.07407407  593
# 1670               1              NA              NA              NA         1.688342             trans-Octodec-2-enoyl-ACP C29H53N2O8PRS 0.2758621 1.827586 0.06896552  621
#                                                       trend_group       z_layer
# 1641 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1643 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1645 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1647 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1651 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1653 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1655 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1657 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1661 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1663 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1665 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1667 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1636 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1637 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1638 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1639 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1640 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1642 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1644 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1646 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1648 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1649 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1650 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1652 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1654 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1656 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1658 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1659 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1660 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1662 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1664 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1666 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1668 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1669 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1670 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2

temp <- dat.test.sig[ sel.ok[unique(c(subsel.diab.prot1, subsel.diab.prot2))], ]
temp$cpd_names <- gsub(pattern = ",", replacement = ";", x = temp$cpd_names)

#write.csv(x = dat.test.sig[ sel.ok[unique(c(subsel.diab.prot1, subsel.diab.prot2))], ], file = "Zoom1-Forslund-SWE-T2D-DECREASING-protein.csv", quote = FALSE, row.names = FALSE )
write.csv(x = temp, file = "Zoom1-Forslund-SWE-T2D-DECREASING-protein.csv", quote = FALSE, row.names = FALSE )

t2d.zoom1.decreasing.proteins <- dat.test.sig$cpd[ sel.ok[unique(c(subsel.diab.prot1, subsel.diab.prot2)) ]]

length(t2d.zoom1.decreasing.proteins) # 35


## Zoom1 - a) plot in restoration ; b) plot in T2D


# a) in restoration

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate
sel <- which(df$cpd_id %in% t2d.zoom1.decreasing.proteins) # 525
df <- df[sel, ]

str(df)
# 'data.frame':	525 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)

  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))

  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = 3.6459, p-value = 0.0002665
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.7406561 

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)


str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  xlab("Reveg age (years)")+ ylab("CPP of indicator group 1a - \U03A3 rel abun (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    )

p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-a-Indicator-group1a-PROTEINS-Trend-with-Age.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")




# b) in T2D

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 145
length(unique(df$cpd_id)) # 7261
145*7261 # 1052845

sel <- which(df$group %in% c("Normal","T2D met neg"))
df <- df[sel, ]
length(unique(df$sample)) # 76
sel <- which(df$cpd_id %in% t2d.zoom1.decreasing.proteins)
df <- df[sel, ]
length(unique(df$cpd_id)) # 35
35*76 # 2660

str(df)
# 'data.frame':	2660 obs. of  7 variables:


res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  4 variables:

unique(res$group) # "T2D met neg" "Normal"
res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # n = 33
y <- res$sum_rel_abun[ res$group == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP of indicator group 1a - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-b-Indicator-group1a-PROTEINS-Trend-with-T2D.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")


# # # # # # 

# ZOOM1
# Proteins group 1b

p <- ggplot(data = filter( dat.test.sig[sel.ok[-subsel.diab.prot1b], ], z_layer == "N:C >0 to 0.2" ) ) + # sel.ok[subsel.diab.prot1]. [-subsel.diab.prot2]
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.2, 0.6)+ ylim(1.5,2.0)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 2, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  annotate(geom = "label_npc", npcx = 0.02, npcy = 0.97, label = "N:C >0 to 0.2", size = 3.5)+ # geom = "text_npc"
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
#dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM1-3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")

# group 1b - increasing prteins which are a side offshoot of decreasing proteins

subsel.diab.prot1b <- which(dat.test.sig$HC_y[sel.ok] > 1.5 & dat.test.sig$HC_y[sel.ok] < 1.85 & 
                              dat.test.sig$OC_x[sel.ok] > 0.2 & dat.test.sig$OC_x[sel.ok] < 0.4 &
                              dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" &
                              dat.test.sig$trend_with_disease[sel.ok] == "Increasing")
length(subsel.diab.prot1b) # 14

dat.test.sig[sel.ok[subsel.diab.prot1b], ]
#           cpd data_for_this_cpd       p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc     alt trend_with_disease incResto_incT2D decResto_decT2D incResto_decT2D
# 537  cpd27645              <NA> 0.037357168  -2.204061     -2.244497         880     5.557184e-04   7.075248e-04 greater         Increasing               1              NA              NA
# 538  cpd27644              <NA> 0.032446511  -2.209225     -2.247546         886     5.216812e-04   7.229196e-04 greater         Increasing               1              NA              NA
# 539  cpd23307              <NA> 0.032446511  -2.209225     -2.247546         886     5.216812e-04   7.229196e-04 greater         Increasing               1              NA              NA
# 540  cpd23318              <NA> 0.032446511  -2.209225     -2.247546         886     5.216812e-04   7.229196e-04 greater         Increasing               1              NA              NA
# 541  cpd23346              <NA> 0.032446511  -2.510255     -2.548576         886     2.608406e-04   3.614598e-04 greater         Increasing               1              NA              NA
# 681  cpd01337              <NA> 0.048970075  -2.251675     -2.297176         868     5.572016e-04   4.817884e-04 greater         Increasing               1              NA              NA
# 1269 cpd35273              <NA> 0.006393020  -2.030155     -2.065655         946     7.322578e-04   6.106556e-04 greater         Increasing              NA              NA              NA
# 1509 cpd01311              <NA> 0.008619683  -2.002252     -2.072903         936     1.493606e-03   1.127907e-03 greater         Increasing              NA              NA              NA
# 1538 cpd23593              <NA> 0.040951934  -2.673106     -2.712779         876     1.853190e-04   1.449014e-04 greater         Increasing              NA              NA              NA
# 1539 cpd29078              <NA> 0.020778970  -3.221568     -3.272807         904     6.681670e-05   5.659618e-05 greater         Increasing              NA              NA              NA
# 1795 cpd24377              <NA> 0.030196064  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05 greater         Increasing               1              NA              NA
# 1797 cpd24558              <NA> 0.030196064  -3.204333     -3.317360         889     1.431453e-04   1.103263e-04 greater         Increasing               1              NA              NA
# 1798 cpd24559              <NA> 0.030196064  -3.204333     -3.317360         889     1.431453e-04   1.103263e-04 greater         Increasing               1              NA              NA
# 1800 cpd31095              <NA> 0.030196064  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05 greater         Increasing               1              NA              NA
# decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                                                                                                         cpd_names
# 537               NA              NA              NA         1.427626                                                                               NAcMur-Peptide-NAcGlc-Undecaprenols
# 538               NA              NA              NA         1.488832                                                                              NAcMur-4Peptide-NAcGlc-Undecaprenols
# 539               NA              NA              NA         1.488832 N-acetylglucosamine--N-acetylmuramoyl-(tetrapeptide) diphospho-undecaprenol dimer with D,D cross-link (S. aureus)
# 540               NA              NA              NA         1.488832                                            a peptidoglycan with D,D cross-links (meso-diaminopimelate containing)
# 541               NA              NA              NA         1.488832                                                                  a peptidoglycan with D,D cross-link (E. faecium)
# 681               NA              NA              NA         1.310069                                                                                                           Romicil
# 1269               1              NA              NA         2.194294                                                                                      alpha-Kdo8N-(2->6)-lipid IVA
# 1509               1              NA              NA         2.064509                                                                                                      Dethiobiotin
# 1538               1              NA              NA         1.387726                                                                                       CDP-1,2-dipalmitoylglycerol
# 1539               1              NA              NA         1.682376                                                                                    CDP-1-18:1(9Z)-2-16:0-glycerol
# 1795              NA              NA              NA         1.520050                                                                                                     GlcNAc-PP-C50
# 1797              NA              NA              NA         1.520050                                                                                                          CPD-5165
# 1798              NA              NA              NA         1.520050                                                                                                          CPD-5166
# 1800              NA              NA              NA         1.520050                   alpha-L-Rhamnopyranosyl-(1->3)-N-acetyl-alpha-D-glucosaminyl-diphospho-trans,octacis-decaprenol
# cpd_forms      OC_x     HC_y       NC_z         mass                                                    trend_group       z_layer
# 537    C94H152N8O26P2R 0.2765957 1.617021 0.08510638 10000000.000   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 538    C91H147N7O25P2R 0.2747253 1.615385 0.07692308 10000000.000   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 539  C303H493N55O104P4 0.3432343 1.627063 0.18151815     5552.000   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 540   C267H421N31O96P4 0.3595506 1.576779 0.11610487     4610.000   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 541   C279H457N43O92P4 0.3297491 1.637993 0.15412186     4868.000   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 681         C35H62NO12 0.3428571 1.771429 0.02857143      688.000   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 1269    C76H139N3O29P2 0.3815789 1.828947 0.03947368     1620.884 Increasing in T2D (increased exposure in disturbed ecosystems) N:C >0 to 0.2
# 1509        C10H17N2O3 0.3000000 1.700000 0.20000000      213.000 Increasing in T2D (increased exposure in disturbed ecosystems) N:C >0 to 0.2
# 1538     C44H79N3O15P2 0.3409091 1.795455 0.06818182      953.000 Increasing in T2D (increased exposure in disturbed ecosystems) N:C >0 to 0.2
# 1539     C46H81N3O15P2 0.3260870 1.760870 0.06521739      978.104 Increasing in T2D (increased exposure in disturbed ecosystems) N:C >0 to 0.2
# 1795      C58H95NO12P2 0.2068966 1.637931 0.01724138     1061.000   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 1797   C132H218N2O47P2 0.3560606 1.651515 0.01515152     2646.000   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 1798   C138H228N2O52P2 0.3768116 1.652174 0.01449275     2808.000   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 1800     C64H105NO16P2 0.2500000 1.640625 0.01562500     1207.707   Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2

temp <- dat.test.sig[ sel.ok[subsel.diab.prot1b], ] 
temp$cpd_names <- gsub(pattern = ",", replacement = ";", x = temp$cpd_names)

#write.csv(x = dat.test.sig[ sel.ok[subsel.diab.prot1b], ], file = "Zoom1b-Forslund-SWE-T2D-Increasing-protein.csv", quote = FALSE, row.names = FALSE )
write.csv(x = temp, file = "Zoom1b-Forslund-SWE-T2D-Increasing-protein.csv", quote = FALSE, row.names = FALSE )

t2d.zoom1b.increasing.proteins <- dat.test.sig$cpd[ sel.ok[subsel.diab.prot1b] ]

length(t2d.zoom1b.increasing.proteins) # 14


## Zoom1b - c) plot in restoration ; d) plot in T2D


# c) in restoration

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate

sel <- which(df$cpd_id %in% t2d.zoom1b.increasing.proteins) # 210
df <- df[sel, ]

str(df)
# 'data.frame':	210 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = 2.228, p-value = 0.02588
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.4526232 

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)


str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  xlab("Reveg age (years)")+ ylab("CPP of indicator group 1b - \U03A3 rel abun (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(c)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-c-Indicator-group1b-PROTEINS-Trend-with-Age.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")




# b) in T2D

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 145
length(unique(df$cpd_id)) # 7261
145*7261 # 1052845

sel <- which(df$group %in% c("Normal","T2D met neg"))
df <- df[sel, ]
length(unique(df$sample)) # 76
sel <- which(df$cpd_id %in% t2d.zoom1b.increasing.proteins)
df <- df[sel, ]
length(unique(df$cpd_id)) # 14
14*76 # 2660

str(df)
# 'data.frame':	1064 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  4 variables:

unique(res$group) # "T2D met neg" "Normal"
res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # n = 33
y <- res$sum_rel_abun[ res$group == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %   "less"

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP of indicator group 1b - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(d)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-d-Indicator-group1b-PROTEINS-Trend-with-T2D.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")


# # # # # #


# ZOOM1
# Proteins group 1c

p <- ggplot(data = filter( dat.test.sig[sel.ok[subsel.diab.prot1c], ], z_layer == "N:C >0 to 0.2" ) ) + # [-subsel.diab.prot1c]
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.2, 0.6)+ ylim(1.5,2.0)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 2, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  annotate(geom = "label_npc", npcx = 0.02, npcy = 0.97, label = "N:C >0 to 0.2", size = 3.5)+ # geom = "text_npc"
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
#dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM1-3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")

# group 1c - increasing prteins which are a side offshoot of decreasing proteins

subsel.diab.prot1c <- which(dat.test.sig$HC_y[sel.ok] > 1.49 & dat.test.sig$HC_y[sel.ok] < 1.6 & 
                              dat.test.sig$OC_x[sel.ok] > 0.5 & dat.test.sig$OC_x[sel.ok] < 0.6 &
                              dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" &
                              dat.test.sig$trend_with_disease[sel.ok] == "Increasing")
length(subsel.diab.prot1c) # 5

dat.test.sig[sel.ok[subsel.diab.prot1c], ]
# cpd data_for_this_cpd      p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc
# 533  cpd15506              <NA> 0.03244651  -2.510255     -2.548576         886     2.608406e-04   3.614598e-04
# 534  cpd15515              <NA> 0.03244651  -2.510255     -2.548576         886     2.608406e-04   3.614598e-04
# 535  cpd15508              <NA> 0.03244651  -2.510255     -2.548576         886     2.608406e-04   3.614598e-04
# 536  cpd15510              <NA> 0.03244651  -2.209225     -2.247546         886     5.216812e-04   7.229196e-04
# 1792 cpd24627              <NA> 0.03019606  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05
# alt trend_with_disease incResto_incT2D decResto_decT2D incResto_decT2D decResto_incT2D incResto_notT2D
# 533  greater         Increasing               1              NA              NA              NA              NA
# 534  greater         Increasing               1              NA              NA              NA              NA
# 535  greater         Increasing               1              NA              NA              NA              NA
# 536  greater         Increasing               1              NA              NA              NA              NA
# 1792 greater         Increasing               1              NA              NA              NA              NA
# decResto_notT2D minuslog10_p_val
# 533               NA         1.488832
# 534               NA         1.488832
# 535               NA         1.488832
# 536               NA         1.488832
# 1792              NA         1.520050
# cpd_names
# 533                                                           two disacharide linked murein units, tetrapeptide corsslinked tetrapeptide (A2pm->D-ala) (middle of chain)
# 534    three disacharide linked murein units (pentapeptide crosslinked tetrapeptide (A2pm->D-ala) tetrapeptide corsslinked tetrapeptide (A2pm->D-ala)) (middle of chain)
# 535  three disacharide linked murein units (tetrapeptide crosslinked tetrapeptide (A2pm->D-ala) & tetrapeptide corsslinked tetrapeptide (A2pm->D-ala)) (middle of chain)
# 536                                                                   two linked disacharide pentapeptide and tetrapeptide murein units (uncrosslinked, middle of chain)
# 1792                                                                                                                                                            CPD-6242
#            cpd_forms      OC_x     HC_y      NC_z     mass
# 533    C74H112N14O39 0.5270270 1.513514 0.1891892     1820
# 534   C114H172N22O59 0.5175439 1.508772 0.1929825     2792
# 535   C111H167N21O58 0.5225225 1.504505 0.1891892     2721
# 536  C77H119N15O42R2 0.5454545 1.545455 0.1948052     1909
# 1792    C13H20N2O7R2 0.5384615 1.538462 0.1538462 10000000
# trend_group       z_layer
# 533  Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 534  Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 535  Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 536  Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2
# 1792 Increasing in T2D (increased exposure in quality ecosystems) N:C >0 to 0.2


temp <- dat.test.sig[sel.ok[subsel.diab.prot1c], ] 
temp$cpd_names <- gsub(pattern = ",", replacement = ";", x = temp$cpd_names)

#write.csv(x = dat.test.sig[ sel.ok[subsel.diab.prot1c], ], file = "Zoom1c-Forslund-SWE-T2D-Increasing-protein.csv", quote = FALSE, row.names = FALSE )
write.csv(x = temp, file = "Zoom1c-Forslund-SWE-T2D-Increasing-protein.csv", quote = FALSE, row.names = FALSE )

t2d.zoom1c.increasing.proteins <- dat.test.sig$cpd[ sel.ok[subsel.diab.prot1c] ]

length(t2d.zoom1c.increasing.proteins) # 5


## Zoom1b - c) plot in restoration ; d) plot in T2D


# e) in restoration

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate

sel <- which(df$cpd_id %in% t2d.zoom1c.increasing.proteins) # 210
df <- df[sel, ]

str(df)
# 'data.frame':	75 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = 2.8357, p-value = 0.004573
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.5760658

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)


str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  xlab("Reveg age (years)")+ ylab("CPP of indicator group 1c - \U03A3 rel abun (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(e)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-e-Indicator-group1c-PROTEINS-Trend-with-Age.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")




# f) in T2D

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 145
length(unique(df$cpd_id)) # 7261
145*7261 # 1052845

sel <- which(df$group %in% c("Normal","T2D met neg"))
df <- df[sel, ]
length(unique(df$sample)) # 76
sel <- which(df$cpd_id %in% t2d.zoom1c.increasing.proteins)
df <- df[sel, ]
length(unique(df$cpd_id)) # 5
5*76 # 380

str(df)
# 'data.frame':	380 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  4 variables:

unique(res$group) # "T2D met neg" "Normal"
res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # n = 33
y <- res$sum_rel_abun[ res$group == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %   "less"

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP of indicator group 1c - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(f)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-f-Indicator-group1c-PROTEINS-Trend-with-T2D.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")





# # # # # #
# # # # # #

# ZOOM2

p <- ggplot(data = filter( dat.test.sig[sel.ok[subsel.diab.carbline], ], z_layer == "N:C = 0" ) ) +
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.825, 1.05)+ ylim(1.65,2.1)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 2, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  annotate(geom = "label_npc", npcx = 0.02, npcy = 0.97, label = "N:C = 0", size = 3.5)+ # geom = "text_npc"
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM2-3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")



sel.diab.carbs <- which(dat.test.sig$cpd %in% c("cpd00082","cpd00076","cpd03200","cpd01399"))

dat.test.sig[sel.diab.carbs, ]
# cpd data_for_this_cpd       p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc     alt trend_with_disease incResto_incT2D decResto_decT2D
# 498 cpd03200              <NA> 0.045829264 -1.5785034    -1.6362465         871      0.003285949    0.002303270 greater         Increasing              NA              NA
# 646 cpd00076              <NA> 0.037357168 -1.1824895    -1.2590219         880      0.010613702    0.008435426 greater         Increasing              NA              NA
# 703 cpd00082              <NA> 0.012847012 -0.8157642    -0.8458833         922      0.010240515    0.024767669 greater         Increasing              NA              NA
# 954 cpd01399              <NA> 0.008875502 -1.2448601    -1.3119454         935      0.008144642    0.007196424 greater         Increasing               1              NA
# incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val     cpd_names cpd_forms      OC_x     HC_y NC_z mass
# 498              NA               1              NA              NA         1.338857 Manninotriose C18H32O16 0.8888889 1.777778    0  504
# 646              NA               1              NA              NA         1.427626       Sucrose C12H22O11 0.9166667 1.833333    0  342
# 703              NA               1              NA              NA         1.891198    D-Fructose   C6H12O6 1.0000000 2.000000    0  180
# 954              NA              NA              NA              NA         2.051807 Maltotetraose C24H42O21 0.8750000 1.750000    0  666
# trend_group z_layer
# 498 Increasing in T2D (increased exposure in disturbed ecosystems) N:C = 0
# 646 Increasing in T2D (increased exposure in disturbed ecosystems) N:C = 0
# 703 Increasing in T2D (increased exposure in disturbed ecosystems) N:C = 0
# 954   Increasing in T2D (increased exposure in quality ecosystems) N:C = 0


# estimate of line: y = 2x
# (round to 4 d.p.?)
subsel.diab.carbline <- which( round(dat.test.sig$HC_y[sel.ok],4) == round(2*dat.test.sig$OC_x[sel.ok],4) & dat.test.sig$OC_x[sel.ok] > 0.85 & dat.test.sig$OC_x[sel.ok] <= 1 & dat.test.sig$HC_y[sel.ok] > 1.7 & dat.test.sig$HC_y[sel.ok] <= 2 & dat.test.sig$NC_z[sel.ok] == 0 & dat.test.sig$trend_with_disease[sel.ok] == "Increasing")
length(subsel.diab.carbline) # 8

dat.test.sig[ sel.ok[subsel.diab.carbline], ]
#          cpd data_for_this_cpd       p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc     alt trend_with_disease incResto_incT2D decResto_decT2D
# 109 cpd00224              <NA> 0.011169028 -0.8721609    -0.9009420         927     0.0086069568   0.0174619005 greater         Increasing              NA              NA
# 497 cpd03198              <NA> 0.006200811 -1.1436927    -1.1972537         947     0.0083342513   0.0097083138 greater         Increasing              NA              NA
# 498 cpd03200              <NA> 0.045829264 -1.5785034    -1.6362465         871     0.0032859493   0.0023032704 greater         Increasing              NA              NA
# 638 cpd00259              <NA> 0.021871276 -1.5836554    -1.6377288         902     0.0030534290   0.0026951598 greater         Increasing              NA              NA
# 646 cpd00076              <NA> 0.037357168 -1.1824895    -1.2590219         880     0.0106137020   0.0084354263 greater         Increasing              NA              NA
# 699 cpd01074              <NA> 0.003507291 -2.4035918    -2.4494292         965     0.0003954819   0.0004409342 greater         Increasing               1              NA
# 703 cpd00082              <NA> 0.012847012 -0.8157642    -0.8458833         922     0.0102405149   0.0247676690 greater         Increasing              NA              NA
# 954 cpd01399              <NA> 0.008875502 -1.2448601    -1.3119454         935     0.0081446419   0.0071964243 greater         Increasing               1              NA
#     incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val     cpd_names cpd_forms      OC_x     HC_y NC_z mass
# 109              NA               1              NA              NA         1.951985   L-Arabinose   C5H10O5 1.0000000 2.000000    0  150
# 497              NA               1              NA              NA         2.207552     Melibiose C12H22O11 0.9166667 1.833333    0  342
# 498              NA               1              NA              NA         1.338857 Manninotriose C18H32O16 0.8888889 1.777778    0  504
# 638              NA               1              NA              NA         1.660126    D-Lyxulose   C5H10O5 1.0000000 2.000000    0  150
# 646              NA               1              NA              NA         1.427626       Sucrose C12H22O11 0.9166667 1.833333    0  342
# 699              NA              NA              NA              NA         2.455028      Nigerose C12H22O11 0.9166667 1.833333    0  342
# 703              NA               1              NA              NA         1.891198    D-Fructose   C6H12O6 1.0000000 2.000000    0  180
# 954              NA              NA              NA              NA         2.051807 Maltotetraose C24H42O21 0.8750000 1.750000    0  666
# trend_group z_layer
# 109 Increasing in T2D (increased exposure in disturbed ecosystems) N:C = 0
# 497 Increasing in T2D (increased exposure in disturbed ecosystems) N:C = 0
# 498 Increasing in T2D (increased exposure in disturbed ecosystems) N:C = 0
# 638 Increasing in T2D (increased exposure in disturbed ecosystems) N:C = 0
# 646 Increasing in T2D (increased exposure in disturbed ecosystems) N:C = 0
# 699   Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 703 Increasing in T2D (increased exposure in disturbed ecosystems) N:C = 0
# 954   Increasing in T2D (increased exposure in quality ecosystems) N:C = 0

write.csv(x = dat.test.sig[ sel.ok[subsel.diab.carbline], ], file = "Zoom2-Forslund-SWE-T2D-Increasing-Carbs.csv", quote = FALSE, row.names = FALSE )

t2d.zoom2.increasing.carbs <- dat.test.sig$cpd[ sel.ok[subsel.diab.carbline] ]

sort( dat.test.sig[ sel.ok[subsel.diab.carbline], "cpd_names" ])
# "D-Fructose"    "D-Lyxulose"    "L-Arabinose"   "Maltotetraose" "Manninotriose" "Melibiose"     "Nigerose"      "Sucrose"    

# adjusted so ONLY THOSE in trend group Increasing in T2D (increased exposure in disturbed ecosystems)
t2d.zoom2B.increasing.carbs <- c("cpd00224", "cpd03198", "cpd03200", "cpd00259", "cpd00076", "cpd00082")

dat.test.sig$cpd_names[ which(dat.test.sig$cpd %in% t2d.zoom2B.increasing.carbs) ]
# "L-Arabinose"   "Melibiose"     "Manninotriose" "D-Lyxulose"    "Sucrose"       "D-Fructose"  


## Zoom1 - a) plot in restoration ; b) plot in T2D


# a) in restoration

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate

#sel <- which(df$cpd_id %in% t2d.zoom2.increasing.carbs) # 120
sel <- which(df$cpd_id %in% t2d.zoom2B.increasing.carbs) # 90
df <- df[sel, ]

str(df)
# 'data.frame':	120 obs. of  7 variables:
# 'data.frame':	90 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = -4.1522, p-value = 3.292e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.843525 

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

#length(t2d.zoom2.increasing.carbs) # 8
length(t2d.zoom2B.increasing.carbs) # 6

str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  xlab("Reveg age (years)")+ ylab("CPP of indicator group 2 - \U03A3 rel abun (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom2-a-Indicator-group2-CARBS-Trend-with-Age.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")



# b) in T2D

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 145
length(unique(df$cpd_id)) # 7261
145*7261 # 1052845

sel <- which(df$group %in% c("Normal","T2D met neg"))
df <- df[sel, ]
length(unique(df$sample)) # 76
#sel <- which(df$cpd_id %in% t2d.zoom2.increasing.carbs)
sel <- which(df$cpd_id %in% t2d.zoom2B.increasing.carbs)
df <- df[sel, ]
length(unique(df$cpd_id)) # 6 8

str(df)
# 'data.frame':	608 obs. of  7 variables:
# 'data.frame':	456 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  4 variables:

unique(res$group) # "T2D met neg" "Normal"
res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # n = 33
y <- res$sum_rel_abun[ res$group == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun % #  "less"

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP of indicator group 2 - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom2-b-Indicator-group2-CARBS-Trend-with-T2D.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")


# # # # # #
# # # # # #


# ZOOM3

p <- ggplot(data = filter( dat.test.sig[sel.ok[-subsel.diab.arom2], ], z_layer == "N:C = 0" ) ) + # [-subsel.diab.arom1] [subsel.diab.carbline] subsel.diab.lignin
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0, 0.75)+ ylim(0.4,1.3)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 2, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  annotate(geom = "label_npc", npcx = 0.94, npcy = 0.05, label = "N:C = 0", size = 3.5)+ # geom = "text_npc"
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM3-3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")


subsel.diab.lignin <- which( dat.test.sig$OC_x[sel.ok] > 0.3 & dat.test.sig$OC_x[sel.ok] < 0.7 & 
                               dat.test.sig$HC_y[sel.ok] > 0.45 & dat.test.sig$HC_y[sel.ok] < 1.35  & 
                               dat.test.sig$NC_z[sel.ok] == 0 & dat.test.sig$trend_with_disease[sel.ok] == "Increasing")
length(subsel.diab.lignin) # 17

dat.test.sig[ sel.ok[subsel.diab.lignin], ]
# cpd data_for_this_cpd      p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc     alt trend_with_disease incResto_incT2D decResto_decT2D
# 568  cpd15017              <NA> 0.03735717  -3.404144     -3.500101         880     7.817171e-05   5.279877e-05 greater         Increasing               1              NA
# 891  cpd00373              <NA> 0.03402135  -3.505067     -3.618390         884     7.178544e-05   5.351324e-05 greater         Increasing               1              NA
# 1779 cpd09265              <NA> 0.03019606  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05 greater         Increasing               1              NA
# 1780 cpd03340              <NA> 0.03019606  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05 greater         Increasing               1              NA
# 1781 cpd09286              <NA> 0.03019606  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05 greater         Increasing               1              NA
# 1782 cpd02450              <NA> 0.03019606  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05 greater         Increasing               1              NA
# 1783 cpd11241              <NA> 0.03019606  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05 greater         Increasing               1              NA
# 1784 cpd04125              <NA> 0.03019606  -3.028242     -3.141269         889     2.147179e-04   1.654894e-04 greater         Increasing               1              NA
# 1785 cpd12896              <NA> 0.03019606  -3.204333     -3.317360         889     1.431453e-04   1.103263e-04 greater         Increasing               1              NA
# 1786 cpd15020              <NA> 0.03019606  -3.204333     -3.317360         889     1.431453e-04   1.103263e-04 greater         Increasing               1              NA
# 1787 cpd15019              <NA> 0.03019606  -3.028242     -3.141269         889     2.147179e-04   1.654894e-04 greater         Increasing               1              NA
# 1793 cpd27078              <NA> 0.03019606  -3.204333     -3.317360         889     1.431453e-04   1.103263e-04 greater         Increasing               1              NA
# 1794 cpd27084              <NA> 0.03019606  -3.204333     -3.317360         889     1.431453e-04   1.103263e-04 greater         Increasing               1              NA
# 1796 cpd24928              <NA> 0.03019606  -3.204333     -3.317360         889     1.431453e-04   1.103263e-04 greater         Increasing               1              NA
# 1799 cpd25021              <NA> 0.03019606  -3.204333     -3.317360         889     1.431453e-04   1.103263e-04 greater         Increasing               1              NA
# 1801 cpd00865              <NA> 0.03019606  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05 greater         Increasing               1              NA
# 1802 cpd31120              <NA> 0.03019606  -3.505363     -3.618390         889     7.157263e-05   5.516313e-05 greater         Increasing               1              NA
# incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                                                          cpd_names   cpd_forms
# 568               NA              NA              NA              NA         1.427626                     Delphinidin 3-O-(6-O-malonyl-beta-D-glucoside)   C24H21O15
# 891               NA              NA              NA              NA         1.468248                                                           Sinapate    C11H11O5
# 1779              NA              NA              NA              NA         1.520050                                               Sophoraflavonoloside   C27H29O16
# 1780              NA              NA              NA              NA         1.520050                                                       Isotrifoliin   C21H19O12
# 1781              NA              NA              NA              NA         1.520050                                                              QUOSP   C27H29O17
# 1782              NA              NA              NA              NA         1.520050                                           Flavonol 3-O-D-glucoside    C21H20O8
# 1783              NA              NA              NA              NA         1.520050               Flavonol 3-O-beta-D-glucosyl-(1->2)-beta-D-glucoside   C27H30O13
# 1784              NA              NA              NA              NA         1.520050                                                    cis-p-Coumarate      C9H7O3
# 1785              NA              NA              NA              NA         1.520050                               4'-O-beta-D-Glucosyl-cis-p-coumarate    C15H17O8
# 1786              NA              NA              NA              NA         1.520050 Delphinidin 3-O-(6''-O-malonyl)-beta-glucoside-3'-O-beta-glucoside   C30H31O20
# 1787              NA              NA              NA              NA         1.520050                                                        Ternatin C5   C36H41O25
# 1793              NA              NA              NA              NA         1.520050                                        Flavonol-3-O-B-D-Glucosides  C21H11O8R9
# 1794              NA              NA              NA              NA         1.520050                                       Flavonol-glucosyl-glucosides C27H21O13R9
# 1796              NA              NA              NA              NA         1.520050                     delphinidin 3-O-(6''-O-malonyl)-beta-glucoside   C24H21O15
# 1799              NA              NA              NA              NA         1.520050                    kaempferol 3-O-beta-D-glucosyl-(1->2)-glucoside   C27H29O16
# 1801              NA              NA              NA              NA         1.520050                                        1-O-Sinapoyl-beta-D-glucose   C17H22O10
# 1802              NA              NA              NA              NA         1.520050                                       1-O-Vanilloyl-beta-D-glucose    C14H18O9
# OC_x      HC_y NC_z         mass                                                  trend_group z_layer
# 568  0.6250000 0.8750000    0 5.470000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 891  0.4545455 1.0000000    0 2.230000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1779 0.5925926 1.0740741    0 6.100000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1780 0.5714286 0.9047619    0 4.640000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1781 0.6296296 1.0740741    0 6.260000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1782 0.3809524 0.9523810    0 4.000000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1783 0.4814815 1.1111111    0 5.620000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1784 0.3333333 0.7777778    0 1.630000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1785 0.5333333 1.1333333    0 3.250000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1786 0.6666667 1.0333333    0 7.100000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1787 0.6944444 1.1388889    0 8.720000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1793 0.3809524 0.5238095    0 1.000000e+07 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1794 0.4814815 0.7777778    0 1.000000e+07 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1796 0.6250000 0.8750000    0 5.500000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1799 0.5925926 1.0740741    0 6.100000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1801 0.5882353 1.2941176    0 3.860000e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 1802 0.6428571 1.2857143    0 3.300951e+02 Increasing in T2D (increased exposure in quality ecosystems) N:C = 0

temp <- dat.test.sig[ sel.ok[subsel.diab.lignin], ]
temp$cpd_names <- gsub(pattern = ",", replacement = ";", x = temp$cpd_names)

#write.csv(x = dat.test.sig[ sel.ok[subsel.diab.lignin], ], file = "Zoom3a-Forslund-SWE-T2D-Increasing-Lignin.csv", quote = FALSE, row.names = FALSE )
write.csv(x = temp, file = "Zoom3a-Forslund-SWE-T2D-Increasing-Lignin.csv", quote = FALSE, row.names = FALSE )

t2d.zoom3a.increasing.lignin <- dat.test.sig$cpd[ sel.ok[subsel.diab.lignin] ]

sort( dat.test.sig[ sel.ok[subsel.diab.lignin], "cpd_names" ])
# [1] "1-O-Sinapoyl-beta-D-glucose"                                        "1-O-Vanilloyl-beta-D-glucose"                                      
# [3] "4'-O-beta-D-Glucosyl-cis-p-coumarate"                               "cis-p-Coumarate"                                                   
# [5] "Delphinidin 3-O-(6-O-malonyl-beta-D-glucoside)"                     "delphinidin 3-O-(6''-O-malonyl)-beta-glucoside"                    
# [7] "Delphinidin 3-O-(6''-O-malonyl)-beta-glucoside-3'-O-beta-glucoside" "Flavonol 3-O-beta-D-glucosyl-(1->2)-beta-D-glucoside"              
# [9] "Flavonol 3-O-D-glucoside"                                           "Flavonol-3-O-B-D-Glucosides"                                       
# [11] "Flavonol-glucosyl-glucosides"                                       "Isotrifoliin"                                                      
# [13] "kaempferol 3-O-beta-D-glucosyl-(1->2)-glucoside"                    "QUOSP"                                                             
# [15] "Sinapate"                                                           "Sophoraflavonoloside"                                              
# [17] "Ternatin C5"


## Zoom1 - a) plot in restoration ; b) plot in T2D


# a) in restoration

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate

sel <- which(df$cpd_id %in% t2d.zoom3a.increasing.lignin) # 255
df <- df[sel, ]

str(df)
# 'data.frame':	255 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = 3.1395, p-value = 0.001692
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.6377872 

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

length(t2d.zoom3a.increasing.lignin ) # 17

str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  xlab("Reveg age (years)")+ ylab("CPP of indicator group 3a - \U03A3 rel abun (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom3-a-Indicator-group3a-Lignin-Trend-with-Age.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")



# b) in T2D

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 145
length(unique(df$cpd_id)) # 7261
145*7261 # 1052845

sel <- which(df$group %in% c("Normal","T2D met neg"))
df <- df[sel, ]
length(unique(df$sample)) # 76
sel <- which(df$cpd_id %in% t2d.zoom3a.increasing.lignin)
df <- df[sel, ]
length(unique(df$cpd_id)) # 17

str(df)
# 'data.frame':	1292 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  4 variables:

unique(res$group) # "T2D met neg" "Normal"
res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # n = 33
y <- res$sum_rel_abun[ res$group == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun % #  "less"

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP of indicator group 3a - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom3-b-Indicator-group3a-Lignin-Trend-with-T2D.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")



# # # # # #

# Group 3b

subsel.diab.arom1 <- which( dat.test.sig$OC_x[sel.ok] > 0.1 & dat.test.sig$OC_x[sel.ok] < 0.4 & 
                               dat.test.sig$HC_y[sel.ok] > 0.4 & dat.test.sig$HC_y[sel.ok] < 1.3  & 
                               dat.test.sig$NC_z[sel.ok] == 0 & #dat.test.sig$trend_with_disease[sel.ok] == "Decreasing")
                               dat.test.sig$trend_group[sel.ok] == "Decreasing in T2D (reduced exposure in disturbed ecosystems)")
length(subsel.diab.arom1) # 10

dat.test.sig[ sel.ok[subsel.diab.arom1], ]
# cpd data_for_this_cpd       p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc  alt trend_with_disease incResto_incT2D decResto_decT2D
# 117  cpd00430              <NA> 0.036499442  -2.252382     -2.194067       538.0    -8.036997e-04  -6.499769e-04 less         Decreasing              NA              NA
# 409  cpd03336              <NA> 0.033892782  -8.482025     -6.616385       542.0    -2.418881e-07  -2.798300e-06 less         Decreasing              NA              NA
# 436  cpd22764              <NA> 0.008839917  -8.482025     -6.917415       493.0    -1.209441e-07   1.082224e-05 less         Decreasing              NA              NA
# 941  cpd27851              <NA> 0.026871296  -8.482025     -4.919994       540.0    -1.202282e-05  -1.576338e-04 less         Decreasing              NA              NA
# 942  cpd27852              <NA> 0.026871296  -8.482025     -4.618964       540.0    -2.404564e-05  -3.152676e-04 less         Decreasing              NA              NA
# 945  cpd01554              <NA> 0.027655868  -8.482025     -5.548132       537.0    -2.830530e-06  -3.467139e-05 less         Decreasing              NA              NA
# 946  cpd01722              <NA> 0.027655868  -8.482025     -5.548132       537.0    -2.830530e-06  -3.467139e-05 less         Decreasing              NA              NA
# 947  cpd06913              <NA> 0.026871296  -8.482025     -5.618964       540.0    -2.404564e-06  -3.152676e-05 less         Decreasing              NA              NA
# 948  cpd34699              <NA> 0.026871296  -8.482025     -5.618964       540.0    -2.404564e-06  -3.152676e-05 less         Decreasing              NA              NA
# 1753 cpd00435              <NA> 0.005059234  -8.482025     -7.041254       492.5    -9.093812e-08  -2.617036e-06 less         Decreasing              NA              NA
# incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                       cpd_names cpd_forms      OC_x      HC_y NC_z        mass
# 117                1              NA              NA              NA         1.437714                            PACT    C8H8O2 0.2500000 1.0000000    0 1.36000e+02
# 409                1              NA              NA              NA         1.469893                3-Chlorocatechol  C6H5ClO2 0.3333333 0.8333333    0 1.44000e+02
# 436                1              NA              NA              NA         2.053552          3-chlorobenzyl alcohol   C7H7ClO 0.1428571 1.0000000    0 1.42000e+02
# 941                1              NA              NA              NA         1.570711                 Phenolic-Donors    C6H5OR 0.1666667 0.8333333    0 1.00000e+07
# 942                1              NA              NA              NA         1.570711 Phenoxyl-rad-of-phenolic-donors    C6H4OR 0.1666667 0.6666667    0 1.00000e+07
# 945                1              NA              NA              NA         1.558213                 Sinapyl alcohol  C11H14O4 0.3636364 1.2727273    0 2.10000e+02
# 946                1              NA              NA              NA         1.558213              p-Coumaryl alcohol   C9H10O2 0.2222222 1.1111111    0 1.50000e+02
# 947                1              NA              NA              NA         1.570711                       Baicalein   C15H9O5 0.3333333 0.6000000    0 2.69000e+02
# 948                1              NA              NA              NA         1.570711            6,7-dehydrobaicalein   C15H7O5 0.3333333 0.4666667    0 2.68225e+02
# 1753               1              NA              NA              NA         2.295915                  Phenylcarbinol     C7H8O 0.1428571 1.1428571    0 1.08000e+02
# trend_group z_layer
# 117  Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 409  Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 436  Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 941  Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 942  Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 945  Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 946  Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 947  Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 948  Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 1753 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0

write.csv(x = dat.test.sig[ sel.ok[subsel.diab.arom1], ], file = "Zoom3b-Forslund-SWE-T2D-Decreasing-Arom1.csv", quote = FALSE, row.names = FALSE )

t2d.zoom3b.decreasing.arom1 <- dat.test.sig$cpd[ sel.ok[subsel.diab.arom1] ]

sort( dat.test.sig[ sel.ok[subsel.diab.arom1], "cpd_names" ])
# [1] "3-chlorobenzyl alcohol"          "3-Chlorocatechol"                "6,7-dehydrobaicalein"            "Baicalein"                      
# [5] "p-Coumaryl alcohol"              "PACT"                            "Phenolic-Donors"                 "Phenoxyl-rad-of-phenolic-donors"
# [9] "Phenylcarbinol"                  "Sinapyl alcohol" 


## Zoom1 - a) plot in restoration ; b) plot in T2D


# a) in restoration

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate

sel <- which(df$cpd_id %in% t2d.zoom3b.decreasing.arom1) # 150
df <- df[sel, ]

str(df)
# 'data.frame':	150 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = 3.9497, p-value = 7.825e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.8023774 

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

length(t2d.zoom3b.decreasing.arom1 ) # 10

str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  xlab("Reveg age (years)")+ ylab("CPP of indicator group 3b - \U03A3 rel abun (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(c)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom3-c-Indicator-group3b-Arom1-Trend-with-Age.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")



# b) in T2D

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 145
length(unique(df$cpd_id)) # 7261
145*7261 # 1052845

sel <- which(df$group %in% c("Normal","T2D met neg"))
df <- df[sel, ]
length(unique(df$sample)) # 76
sel <- which(df$cpd_id %in% t2d.zoom3b.decreasing.arom1)
df <- df[sel, ]
length(unique(df$cpd_id)) # 10

str(df)
# 'data.frame':	760 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  4 variables:

unique(res$group) # "T2D met neg" "Normal"
res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # n = 33
y <- res$sum_rel_abun[ res$group == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun % #  "less" "greater"

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP of indicator group 3b - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(d)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom3-d-Indicator-group3b-Arom1-Trend-with-T2D.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")



# # # # # #

# Group 3c


subsel.diab.arom2 <- which( dat.test.sig$OC_x[sel.ok] >= 0 & dat.test.sig$OC_x[sel.ok] < 0.4 & 
                              dat.test.sig$HC_y[sel.ok] > 0.4 & dat.test.sig$HC_y[sel.ok] < 1.3  & 
                              dat.test.sig$NC_z[sel.ok] == 0 &
                              dat.test.sig$trend_with_disease[sel.ok] == "Decreasing")
                              #dat.test.sig$trend_group[sel.ok] == "Decreasing in T2D (reduced exposure in disturbed ecosystems)")
length(subsel.diab.arom2) # 24

dat.test.sig[ sel.ok[subsel.diab.arom2], ]
# cpd data_for_this_cpd       p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc  alt trend_with_disease incResto_incT2D decResto_decT2D
# 117  cpd00430              <NA> 0.036499442  -2.252382     -2.194067       538.0    -8.036997e-04  -6.499769e-04 less         Decreasing              NA              NA
# 371  cpd03480              <NA> 0.014026894  -8.482025     -5.905827       509.0    -1.242148e-06  -3.672569e-05 less         Decreasing              NA               1
# 409  cpd03336              <NA> 0.033892782  -8.482025     -6.616385       542.0    -2.418881e-07  -2.798300e-06 less         Decreasing              NA              NA
# 410  cpd02022              <NA> 0.008809138  -8.482025     -7.117507       502.5    -7.629450e-08  -1.796933e-06 less         Decreasing              NA               1
# 413  cpd19545              <NA> 0.008809138  -8.482025     -7.117507       502.5    -7.629450e-08  -1.796933e-06 less         Decreasing              NA               1
# 436  cpd22764              <NA> 0.008839917  -8.482025     -6.917415       493.0    -1.209441e-07   1.082224e-05 less         Decreasing              NA              NA
# 941  cpd27851              <NA> 0.026871296  -8.482025     -4.919994       540.0    -1.202282e-05  -1.576338e-04 less         Decreasing              NA              NA
# 942  cpd27852              <NA> 0.026871296  -8.482025     -4.618964       540.0    -2.404564e-05  -3.152676e-04 less         Decreasing              NA              NA
# 945  cpd01554              <NA> 0.027655868  -8.482025     -5.548132       537.0    -2.830530e-06  -3.467139e-05 less         Decreasing              NA              NA
# 946  cpd01722              <NA> 0.027655868  -8.482025     -5.548132       537.0    -2.830530e-06  -3.467139e-05 less         Decreasing              NA              NA
# 947  cpd06913              <NA> 0.026871296  -8.482025     -5.618964       540.0    -2.404564e-06  -3.152676e-05 less         Decreasing              NA              NA
# 948  cpd34699              <NA> 0.026871296  -8.482025     -5.618964       540.0    -2.404564e-06  -3.152676e-05 less         Decreasing              NA              NA
# 1710 cpd04041              <NA> 0.002833258  -8.482025     -7.158152       474.0    -6.947809e-08  -1.451314e-06 less         Decreasing              NA               1
# 1711 cpd04042              <NA> 0.002833258  -8.482025     -7.158152       474.0    -6.947809e-08  -1.451314e-06 less         Decreasing              NA               1
# 1713 cpd14978              <NA> 0.002833258  -8.482025     -7.158152       474.0    -6.947809e-08  -1.451314e-06 less         Decreasing              NA               1
# 1714 cpd14979              <NA> 0.002833258  -8.482025     -6.857122       474.0    -1.389562e-07  -2.676475e-06 less         Decreasing              NA               1
# 1717 cpd19522              <NA> 0.002833258  -8.482025     -7.158152       474.0    -6.947809e-08  -1.451314e-06 less         Decreasing              NA               1
# 1718 cpd19536              <NA> 0.002833258  -8.482025     -7.158152       474.0    -6.947809e-08  -1.451314e-06 less         Decreasing              NA               1
# 1719 cpd19551              <NA> 0.002833258  -8.482025     -7.158152       474.0    -6.947809e-08  -1.451314e-06 less         Decreasing              NA               1
# 1720 cpd19546              <NA> 0.002833258  -8.482025     -6.857122       474.0    -1.389562e-07  -2.902629e-06 less         Decreasing              NA               1
# 1721 cpd19550              <NA> 0.002833258  -8.482025     -7.158152       474.0    -6.947809e-08  -1.451314e-06 less         Decreasing              NA               1
# 1726 cpd03479              <NA> 0.002833258  -8.482025     -6.681031       474.0    -2.084343e-07  -4.353943e-06 less         Decreasing              NA               1
# 1727 cpd01533              <NA> 0.002833258  -8.482025     -6.681031       474.0    -2.084343e-07  -4.353943e-06 less         Decreasing              NA               1
# 1753 cpd00435              <NA> 0.005059234  -8.482025     -7.041254       492.5    -9.093812e-08  -2.617036e-06 less         Decreasing              NA              NA
# incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                                     cpd_names cpd_forms      OC_x      HC_y NC_z
# 117                1              NA              NA              NA         1.437714                                          PACT    C8H8O2 0.2500000 1.0000000    0
# 371               NA              NA              NA              NA         1.853038                        2-Hydroxyphenylacetate    C8H7O3 0.3750000 0.8750000    0
# 409                1              NA              NA              NA         1.469893                              3-Chlorocatechol  C6H5ClO2 0.3333333 0.8333333    0
# 410               NA              NA              NA              NA         2.055067                         Phenanthrene-3,4-diol  C14H10O2 0.1428571 0.7142857    0
# 413               NA              NA              NA              NA         2.055067                   4,5-Dihydroxybenzo[a]pyrene  C20H12O2 0.1000000 0.6000000    0
# 436                1              NA              NA              NA         2.053552                        3-chlorobenzyl alcohol   C7H7ClO 0.1428571 1.0000000    0
# 941                1              NA              NA              NA         1.570711                               Phenolic-Donors    C6H5OR 0.1666667 0.8333333    0
# 942                1              NA              NA              NA         1.570711               Phenoxyl-rad-of-phenolic-donors    C6H4OR 0.1666667 0.6666667    0
# 945                1              NA              NA              NA         1.558213                               Sinapyl alcohol  C11H14O4 0.3636364 1.2727273    0
# 946                1              NA              NA              NA         1.558213                            p-Coumaryl alcohol   C9H10O2 0.2222222 1.1111111    0
# 947                1              NA              NA              NA         1.570711                                     Baicalein   C15H9O5 0.3333333 0.6000000    0
# 948                1              NA              NA              NA         1.570711                          6,7-dehydrobaicalein   C15H7O5 0.3333333 0.4666667    0
# 1710              NA              NA              NA              NA         2.547714                                          DDMU  C14H9Cl3 0.0000000 0.6428571    0
# 1711              NA              NA              NA              NA         2.547714                                          DDMS C14H11Cl3 0.0000000 0.7857143    0
# 1713              NA              NA              NA              NA         2.547714           1-Hydro-1,1a-dihydroxy-9-fluorenone  C13H10O3 0.2307692 0.7692308    0
# 1714              NA              NA              NA              NA         2.547714              2,3-Dihydroxy-2'-carboxybiphenyl   C13H9O4 0.3076923 0.6923077    0
# 1717              NA              NA              NA              NA         2.547714 cis-3,4-Phenanthrenedihydrodiol-4-carboxylate  C15H11O4 0.2666667 0.7333333    0
# 1718              NA              NA              NA              NA         2.547714            Benzo[a]pyrene-cis-4,5-dihydrodiol  C20H14O2 0.1000000 0.7000000    0
# 1719              NA              NA              NA              NA         2.547714        Benzo[a]pyrene-trans-11,12-dihydrodiol  C20H14O2 0.1000000 0.7000000    0
# 1720              NA              NA              NA              NA         2.547714                 11,12-Dihydroxybenzo[a]pyrene  C20H12O2 0.1000000 0.6000000    0
# 1721              NA              NA              NA              NA         2.547714          Benzo[a]pyrene-cis-11,12-dihydrodiol  C20H14O2 0.1000000 0.7000000    0
# 1726              NA              NA              NA              NA         2.547714                                        Rattex    C9H6O2 0.2222222 0.6666667    0
# 1727              NA              NA              NA              NA         2.547714                               Dihydrocoumarin    C9H8O2 0.2222222 0.8888889    0
# 1753               1              NA              NA              NA         2.295915                                Phenylcarbinol     C7H8O 0.1428571 1.1428571    0
#             mass                                                  trend_group z_layer
# 117  1.36000e+02 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 371  1.51000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 409  1.44000e+02 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 410  2.10000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 413  2.84000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 436  1.42000e+02 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 941  1.00000e+07 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 942  1.00000e+07 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 945  2.10000e+02 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 946  1.50000e+02 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 947  2.69000e+02 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 948  2.68225e+02 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0
# 1710 2.82000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1711 2.84000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1713 2.14000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1714 2.29000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1717 2.56000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1718 2.86000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1719 2.86000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1720 2.84000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1721 2.86000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1726 1.46000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1727 1.48000e+02   Decreasing in T2D (reduced exposure in quality ecosystems) N:C = 0
# 1753 1.08000e+02 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C = 0

temp <- dat.test.sig[ sel.ok[subsel.diab.arom2], ]
temp$cpd_names <- gsub(pattern = ",", replacement = ";", x = temp$cpd_names)

#write.csv(x = dat.test.sig[ sel.ok[subsel.diab.arom2], ], file = "Zoom3c-Forslund-SWE-T2D-Decreasing-Arom2.csv", quote = FALSE, row.names = FALSE )
write.csv(x = temp, file = "Zoom3c-Forslund-SWE-T2D-Decreasing-Arom2.csv", quote = FALSE, row.names = FALSE )

t2d.zoom3c.decreasing.arom2 <- dat.test.sig$cpd[ sel.ok[subsel.diab.arom2] ]

sort( dat.test.sig[ sel.ok[subsel.diab.arom2], "cpd_names" ])
# [1] "1-Hydro-1,1a-dihydroxy-9-fluorenone"           "11,12-Dihydroxybenzo[a]pyrene"                 "2-Hydroxyphenylacetate"                       
# [4] "2,3-Dihydroxy-2'-carboxybiphenyl"              "3-chlorobenzyl alcohol"                        "3-Chlorocatechol"                             
# [7] "4,5-Dihydroxybenzo[a]pyrene"                   "6,7-dehydrobaicalein"                          "Baicalein"                                    
# [10] "Benzo[a]pyrene-cis-11,12-dihydrodiol"          "Benzo[a]pyrene-cis-4,5-dihydrodiol"            "Benzo[a]pyrene-trans-11,12-dihydrodiol"       
# [13] "cis-3,4-Phenanthrenedihydrodiol-4-carboxylate" "DDMS"                                          "DDMU"                                         
# [16] "Dihydrocoumarin"                               "p-Coumaryl alcohol"                            "PACT"                                         
# [19] "Phenanthrene-3,4-diol"                         "Phenolic-Donors"                               "Phenoxyl-rad-of-phenolic-donors"              
# [22] "Phenylcarbinol"                                "Rattex"                                        "Sinapyl alcohol" 


## Zoom1 - a) plot in restoration ; b) plot in T2D


# a) in restoration

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate

sel <- which(df$cpd_id %in% t2d.zoom3c.decreasing.arom2) # 150
df <- df[sel, ]

str(df)
# 'data.frame':	360 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = 3.342, p-value = 0.0008317
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.6789347

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

length(t2d.zoom3c.decreasing.arom2 ) # 24

str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  xlab("Reveg age (years)")+ ylab("CPP of indicator group 3c - \U03A3 rel abun (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(e)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom3-e-Indicator-group3c-Arom2-Trend-with-Age.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")



# b) in T2D

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 145
length(unique(df$cpd_id)) # 7261
145*7261 # 1052845

sel <- which(df$group %in% c("Normal","T2D met neg"))
df <- df[sel, ]
length(unique(df$sample)) # 76
sel <- which(df$cpd_id %in% t2d.zoom3c.decreasing.arom2)
df <- df[sel, ]
length(unique(df$cpd_id)) # 24

str(df)
# 'data.frame':	1824 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  4 variables:

unique(res$group) # "T2D met neg" "Normal"
res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # n = 33
y <- res$sum_rel_abun[ res$group == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun % #  "less" "greater"

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP of indicator group 3c - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(f)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom3-f-Indicator-group3c-Arom2-Trend-with-T2D.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")



# # # # # #
# # # # # #

# Group 4

## ZOOM-IN - Proteins? N:C >0 to 0.2
# ZOOM4

p <- ggplot(data = filter( dat.test.sig[sel.ok[-subsel.diab.lip], ], z_layer == "N:C >0 to 0.2" ) ) + # [subsel.diab.lip]
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.05, 0.25)+ ylim(1.5,2.5)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 2, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  annotate(geom = "label_npc", npcx = 0.02, npcy = 0.97, label = "N:C >0 to 0.2", size = 3.5)+ # geom = "text_npc"
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM4-3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-VS-NORM--Wilcox",this_study,header,".tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")






subsel.diab.lip <- which( dat.test.sig$OC_x[sel.ok] > 0.1 & dat.test.sig$OC_x[sel.ok] < 0.2 &
                            dat.test.sig$HC_y[sel.ok] > 1.7 & dat.test.sig$HC_y[sel.ok] < 2.25 &
                            dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & 
                            dat.test.sig$trend_with_disease[sel.ok] == "Decreasing")
length(subsel.diab.lip) # 7


dat.test.sig[ sel.ok[subsel.diab.lip], ]
# cpd data_for_this_cpd       p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc  alt trend_with_disease incResto_incT2D decResto_decT2D
# 91  cpd01879              <NA> 0.037357168  -3.108548     -2.954499       539.0    -3.316072e-04  -2.222514e-04 less         Decreasing              NA               1
# 94  cpd00449              <NA> 0.036143812  -2.571833     -2.508260       537.5    -4.225042e-04  -7.345017e-04 less         Decreasing              NA               1
# 95  cpd00213              <NA> 0.032941270  -2.582182     -2.482603       533.5    -6.744343e-04  -1.283255e-03 less         Decreasing              NA               1
# 96  cpd19172              <NA> 0.047823070  -3.675007     -3.632468       550.0    -2.174899e-05  -6.194412e-05 less         Decreasing              NA               1
# 98  cpd27418              <NA> 0.038231077  -2.587774     -2.562983       540.0    -1.517695e-04  -1.572807e-04 less         Decreasing              NA               1
# 99  cpd27764              <NA> 0.047823070  -3.107432     -3.031200       550.0    -1.498286e-04  -2.357112e-04 less         Decreasing              NA               1
# 423 cpd25067              <NA> 0.009401624  -5.663850     -5.405817       485.0    -1.759653e-06  -1.041943e-05 less         Decreasing              NA               1
# incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                             cpd_names      cpd_forms      OC_x     HC_y       NC_z
# 91               NA              NA              NA              NA         1.427626                  3-Dehydrosphinganine      C18H38NO2 0.1111111 2.111111 0.05555556
# 94               NA              NA              NA              NA         1.441966                      Dihydrolipoamide      C8H17NOS2 0.1250000 2.125000 0.12500000
# 95               NA              NA              NA              NA         1.482260                             Lipoamide      C8H15NOS2 0.1250000 1.875000 0.12500000
# 96               NA              NA              NA              NA         1.320363       Enzyme N6-(dihydrolipoyl)lysine     C8H16NORS2 0.1250000 2.000000 0.12500000
# 98               NA              NA              NA              NA         1.417583        Lipoyl-Protein-N6-lipoyllysine C14H25N2O2R2S2 0.1428571 1.785714 0.14285714
# 99               NA              NA              NA              NA         1.320363 Oxo-glutarate-dehydrogenase-DH-lipoyl C15H29N2O2R2S2 0.1333333 1.933333 0.13333333
# 423              NA              NA              NA              NA         2.026797     1-16:0-2-18:2-phosphatidylcholine     C42H80NO8P 0.1904762 1.904762 0.02380952
# mass                                                trend_group       z_layer
# 91       300 Decreasing in T2D (reduced exposure in quality ecosystems) N:C >0 to 0.2
# 94       207 Decreasing in T2D (reduced exposure in quality ecosystems) N:C >0 to 0.2
# 95       205 Decreasing in T2D (reduced exposure in quality ecosystems) N:C >0 to 0.2
# 96  10000000 Decreasing in T2D (reduced exposure in quality ecosystems) N:C >0 to 0.2
# 98  10000000 Decreasing in T2D (reduced exposure in quality ecosystems) N:C >0 to 0.2
# 99  10000000 Decreasing in T2D (reduced exposure in quality ecosystems) N:C >0 to 0.2
# 423      757 Decreasing in T2D (reduced exposure in quality ecosystems) N:C >0 to 0.2

write.csv(x = dat.test.sig[ sel.ok[subsel.diab.lip], ], file = "Zoom4-Forslund-SWE-T2D-DECREASING-lip.csv", quote = FALSE, row.names = FALSE )

t2d.zoom4.decreasing.lip <- dat.test.sig$cpd[ sel.ok[subsel.diab.lip]]

length(t2d.zoom4.decreasing.lip ) # 7


## Zoom1 - a) plot in restoration ; b) plot in T2D


# a) in restoration

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate
sel <- which(df$cpd_id %in% t2d.zoom4.decreasing.lip) # 105
df <- df[sel, ]

str(df)
# 'data.frame':	105 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = -3.5446, p-value = 0.0003932
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#        tau 
# -0.7200823 

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)


str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  xlab("Reveg age (years)")+ ylab("CPP of indicator group 4 - \U03A3 rel abun (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom4-Indicator-group4-Lip-Trend-with-Age.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")




# b) in T2D

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 145
length(unique(df$cpd_id)) # 7261
145*7261 # 1052845

sel <- which(df$group %in% c("Normal","T2D met neg"))
df <- df[sel, ]
length(unique(df$sample)) # 76
sel <- which(df$cpd_id %in% t2d.zoom4.decreasing.lip)
df <- df[sel, ]
length(unique(df$cpd_id)) # 7
7*76 # 532

str(df)
# 'data.frame':	532 obs. of  7 variables:


res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  4 variables:

unique(res$group) # "T2D met neg" "Normal"
res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # n = 33
y <- res$sum_rel_abun[ res$group == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP of indicator group 4 - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom4-b-Indicator-group4-Lip-Trend-with-T2D.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------


#### Cross-check compounds with functions ??
#    in example of decreasing ACP proteins
#-------------------------

#saveRDS(object = df.out, file = "df.out--tidy-compounds_indiv--cpp3d-Forslund-SWE-T2D.RDS")
df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-Forslund-SWE-T2D.RDS" )

str(df.out)
# 'data.frame':	545806 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_2" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn25279" "rxn25279" "rxn25279" "rxn03060" ...
# $ cpd_id             : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ cpd_name           : chr  "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "2-methyl-trans-aconitate" "2-Methylcitrate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.3333 0.3333 0.3333 0.0556 0.0556 ...

head(df.out)
# superfocus_fxn f                                                                         f__in   rxn_id   cpd_id                                     cpd_name cpd_form
# 2          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 3          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd00001                                          H2O      H2O
# 4          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 5          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 6          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd00001                                          H2O      H2O
# 7          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# cpd_molar_prop cpd_molar_prop_norm
# 2              1          0.33333333
# 3              1          0.33333333
# 4              1          0.33333333
# 5              1          0.05555556
# 6              1          0.05555556
# 7              1          0.05555556



# look for these compounds across all T2D functions?

t2d.zoom1.decreasing.proteins
# [1] "cpd11499" "cpd11503" "cpd11507" "cpd11511" "cpd11524" "cpd11528" "cpd11532" "cpd11536" "cpd11549" "cpd11553" "cpd11557" "cpd11561" "cpd11475" "cpd11465" "cpd11469"
# [16] "cpd11473" "cpd11498" "cpd11502" "cpd11506" "cpd11510" "cpd11514" "cpd11518" "cpd11523" "cpd11527" "cpd11531" "cpd11535" "cpd11539" "cpd11543" "cpd11548" "cpd11552"
# [31] "cpd11556" "cpd11560" "cpd11564" "cpd11568" "cpd11572"



#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

df.abun <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")
str(df.abun)
# 'data.frame':	1052845 obs. of  7 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260132" "ERR260132" "ERR260132" "ERR260132" ...
# $ cpd_rel_abun: num  0 5.647932 0 0 0.000186 ...
# $ log10_abun  : num  -8.482 0.752 -8.482 -8.482 -3.73 ...
# $ group       : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ group_label : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ ord_group   : int  2 2 2 2 2 2 2 2 2 2 ...

unique(df.abun$group)
# [1] IGT         T2D met neg Normal      T2D met pos
# Levels: Normal < IGT < T2D met pos < T2D met neg

# only consider Normal and T2D met neg
sel <- which(df.abun$group %in% c("Normal", "T2D met neg"))
df.abun <- df.abun[sel, ]

# now consider only compounds of interest
sel <- which(df.abun$cpd_id %in% t2d.zoom1.decreasing.proteins)
df.abun <- df.abun[sel, ]
str(df.abun)
# 'data.frame':	2660 obs. of  7 variables:
# $ cpd_id      : chr  "cpd11475" "cpd11465" "cpd11469" "cpd11473" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  3.54e-06 3.54e-06 3.54e-06 3.54e-06 3.54e-06 ...
# $ log10_abun  : num  -5.45 -5.45 -5.45 -5.45 -5.45 ...
# $ group       : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 4 4 4 4 4 4 4 4 4 4 ...
# $ group_label : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 4 4 4 4 4 4 4 4 4 4 ...
# $ ord_group   : int  4 4 4 4 4 4 4 4 4 4 ...

table(df.abun$group)
# Normal         IGT T2D met pos T2D met neg 
# 1505           0           0        1155 

p <- ggplot(df.abun,aes(x = cpd_rel_abun, fill = group_label))+
  theme_bw()+
  guides(fill = guide_legend(title = "Diagnosis"))+
  geom_histogram(alpha = 0.6)+
  xlab("\U03A3 CPP (%)")+ ylab("Count")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p
dev.print(tiff, filename = paste0(workdir,"/plots/","Histogram-Decreasing-proteins-T2D-vs-Normal.tiff"),
          width = 10, height = 8, units = "cm", res=600, compression = "lzw",type="cairo" )


# subset function - compound info to compounds of interest
head( df.out )
# superfocus_fxn f                                                                         f__in   rxn_id   cpd_id                                     cpd_name cpd_form
# 2          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 3          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd00001                                          H2O      H2O
# 4          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 5          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 6          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd00001                                          H2O      H2O
# 7          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# cpd_molar_prop cpd_molar_prop_norm
# 2              1          0.33333333
# 3              1          0.33333333
# 4              1          0.33333333
# 5              1          0.05555556
# 6              1          0.05555556
# 7              1          0.05555556

sel <- which(df.out$cpd_id %in% t2d.zoom1.decreasing.proteins) # 72
df.out <- df.out[sel, ]


unique(df.out$superfocus_fxn)
# "fxn_9916"  "fxn_9950"  "fxn_9951"  "fxn_10029"



unique(df.out$f) # 2 1

unique(df.out$f__in)
# [1] "Trans-2-decenoyl-[acyl-carrier-protein] isomerase (EC 5.3.3.14)" 
# "Enoyl-.acyl-carrier-protein. reductase .NADH. .EC 1.3.1.9."     
# [3] "Enoyl-.acyl-carrier-protein. reductase .NADPH. .EC 1.3.1.10." 

unique(df.out$rxn_id)
# [1] "rxn07455" "rxn05322" "rxn05324" "rxn05326" "rxn05327" "rxn05433" "rxn05434" "rxn05435" "rxn05436" "rxn05437" "rxn05438" "rxn05439" "rxn05440" "rxn05441" "rxn05442"
# [16] "rxn05443" "rxn05444" "rxn05445" "rxn05446" "rxn05447" "rxn05448" "rxn05449" "rxn05450" "rxn05464" "rxn05353" "rxn05355" "rxn05356" "rxn05357" "rxn05362" "rxn05366"
# [31] "rxn05370" "rxn05374" "rxn05378" "rxn05382" "rxn05387" "rxn05391" "rxn05395" "rxn05399" "rxn05403" "rxn05407" "rxn05412" "rxn05416" "rxn05420" "rxn05424" "rxn05428"
# [46] "rxn05432" "rxn05463"


phy <- readRDS("phy-phyloseq-object-Forslund-SWE-T2D.RDS")

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 19099    4

sel <- which(row.names(df.tax) %in% c("fxn_9916",  "fxn_9950",  "fxn_9951",  "fxn_10029"))
df.tax[sel, ]
#                                      subsys_L1   subsys_L2                         subsys_L3
# fxn_9916  Fatty Acids, Lipids, and Isoprenoids Fatty acids     Fatty Acid Biosynthesis FASII
# fxn_9950  Fatty Acids, Lipids, and Isoprenoids Fatty acids     Fatty Acid Biosynthesis FASII
# fxn_9951  Fatty Acids, Lipids, and Isoprenoids Fatty acids     Fatty Acid Biosynthesis FASII
# fxn_10029 Fatty Acids, Lipids, and Isoprenoids Fatty acids Unsaturated Fatty Acid Metabolism
# fxn
# fxn_9916  3-hydroxyacyl-[acyl-carrier-protein]_dehydratase,_FabA_form_(EC_4.2.1.59)_@_Trans-2-decenoyl-[acyl-carrier-protein]_isomerase_(EC_5.3.3.14)
# fxn_9950                                                                                   Enoyl-[acyl-carrier-protein]_reductase_[NADH]_(EC_1.3.1.9)
# fxn_9951                                                                                 Enoyl-[acyl-carrier-protein]_reductase_[NADPH]_(EC_1.3.1.10)
# fxn_10029 3-hydroxyacyl-[acyl-carrier-protein]_dehydratase,_FabA_form_(EC_4.2.1.59)_@_Trans-2-decenoyl-[acyl-carrier-protein]_isomerase_(EC_5.3.3.14)


# linked to these enzymes: 
# EC 5.3.3.14 
# EC 1.3.1.9
# EC 1.3.1.10


# EC 5.3.3.14 - https://www.genome.jp/dbget-bin/www_bget?ec:5.3.3.14
# While the enzyme from Escherichia coli is highly specific for the 10-carbon enoyl-ACP, the enzyme from 
# Streptococcus pneumoniae can also use the 12-carbon enoyl-ACP as substrate in vitro but not 14- or 16-carbon enoyl-ACPs [3]. 
# ACP can be replaced by either CoA or N-acetylcysteamine thioesters. The cis-3-enoyl product is required to form 
# unsaturated fatty acids, such as palmitoleic acid and cis-vaccenic acid, in dissociated (or type II) fatty-acid biosynthesis.

# Orthology
# K18474 - trans-2-decenoyl-[acyl-carrier protein] isomerase
# - https://www.genome.jp/entry/K18474
# Reference	
# PMID:12237320
# Authors	
# Marrakchi H, Choi KH, Rock CO
# Title	
# A new mechanism for anaerobic unsaturated fatty acid formation in Streptococcus pneumoniae.
# Journal	
# J Biol Chem 277:44809-16 (2002)
# DOI:10.1074/jbc.M208920200

# K01716 - 3-hydroxyacyl-[acyl-carrier protein] dehydratase / trans-2-decenoyl-[acyl-carrier protein] isomerase
# - https://www.genome.jp/entry/K01716
# Reference	
# PMID:2832401
# Authors	
# Cronan JE Jr, Li WB, Coleman R, Narasimhan M, de Mendoza D, Schwab JM
# Title	
# Derived amino acid sequence and identification of active site residues of Escherichia coli beta-hydroxydecanoyl thioester dehydrase.
# Journal	
# J Biol Chem 263:4641-6 (1988)
# DOI:10.1016/S0021-9258(18)68830-1


# EC 1.3.1.9 - https://www.genome.jp/dbget-bin/www_bget?enzyme+1.3.1.9
# The enzyme catalyses an essential step in fatty acid biosynthesis, the reduction of the 2,3-double bond 
# in enoyl-acyl-[acyl-carrier-protein] derivatives of the elongating fatty acid moiety. The enzyme from the 
# bacterium Escherichia coli accepts substrates with carbon chain length from 4 to 18 [3]. 
# The FAS-I enzyme from the bacterium Mycobacterium tuberculosis prefers substrates with carbon chain length from 12 to 24 carbons.


# EC 1.3.1.10 - https://www.genome.jp/dbget-bin/www_bget?ec:1.3.1.10
# One of the activities of EC 2.3.1.86, fatty-acyl-CoA synthase system, an enzyme found in yeasts (Ascomycota and Basidiomycota). 
# Catalyses the reduction of enoyl-acyl-[acyl-carrier protein] derivatives of carbon chain length from 4 to 16. 
# The yeast enzyme is Si-specific with respect to NADP+. cf. EC 1.3.1.39, enoyl-[acyl-carrier-protein] reductase (NADPH, Re-specific) and 
# EC 1.3.1.104, enoyl-[acyl-carrier-protein] reductase (NADPH), which describes enzymes whose stereo-specificity towards NADPH is not known. See also EC 1.3.1.9, enoyl-[acyl-carrier-protein] reductase (NADH).



### Collate compounds of interest and append related functions

df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-Forslund-SWE-T2D.RDS" )

str(df.out)
# 'data.frame':	545806 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_2" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn25279" "rxn25279" "rxn25279" "rxn03060" ...
# $ cpd_id             : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ cpd_name           : chr  "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "2-methyl-trans-aconitate" "2-Methylcitrate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.3333 0.3333 0.3333 0.0556 0.0556 ...

phy <- readRDS("phy-phyloseq-object-Forslund-SWE-T2D.RDS")
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 19099    4


zoom_results <- list()

## Zoom 1a
temp <- read.csv(file = "Zoom1-Forslund-SWE-T2D-DECREASING-protein.csv", header = TRUE)
temp$zoom_group <- "1a"
temp$zoom_group_description <- "Decreasing proteins"

dim(temp) # 35 27

names(temp)
# [1] "cpd"                    "data_for_this_cpd"      "p_val"                  "median_t2d"             "median_Normal"          "W_statistic"            "diff_median_perc"      
# [8] "diff_mean_perc"         "alt"                    "trend_with_disease"     "incResto_incT2D"        "decResto_decT2D"        "incResto_decT2D"        "decResto_incT2D"       
# [15] "incResto_notT2D"        "decResto_notT2D"        "minuslog10_p_val"       "cpd_names"              "cpd_forms"              "OC_x"                   "HC_y"                  
# [22] "NC_z"                   "mass"                   "trend_group"            "z_layer"                "zoom_group"             "zoom_group_description"

# add functions?
temp$fxns <- NA
temp$superfocus_fxns <- NA
temp$subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  print(paste0("completed ",i))
}
temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
zoom_results[["1a. Decreasing proteins"]] <- temp
write.csv(x = temp, file = "Zoom1-Forslund-SWE-T2D-DECREASING-protein--with-functions.csv", quote = FALSE, row.names = FALSE )




## Zoom2-Forslund-SWE-T2D-Increasing-Carbs.csv
temp <- read.csv(file = "Zoom2-Forslund-SWE-T2D-Increasing-Carbs.csv", header = TRUE)
temp$zoom_group <- "2"
temp$zoom_group_description <- "Increasing carbs"

dim(temp) # 8 27

# add functions?
temp$fxns <- NA
temp$superfocus_fxns <- NA
temp$subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  print(paste0("completed ",i))
}

temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
zoom_results[["2. Increasing carbs"]] <- temp
write.csv(x = temp, file = "Zoom2-Forslund-SWE-T2D-Increasing-Carbs--with-functions.csv", quote = FALSE, row.names = FALSE )



## Zoom1b-Forslund-SWE-T2D-Increasing-protein.csv
temp <- read.csv(file = "Zoom1b-Forslund-SWE-T2D-Increasing-protein.csv", header = TRUE)
#temp <- read.table(file = "Zoom1b-Forslund-SWE-T2D-Increasing-protein.csv", header = TRUE, sep = ",")
temp$zoom_group <- "1b"
temp$zoom_group_description <- "Increasing proteins"

dim(temp) # 14 27

# add functions?
temp$fxns <- NA
temp$superfocus_fxns <- NA
temp$subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  print(paste0("completed ",i))
}
temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
zoom_results[["1b. Increasing proteins"]] <- temp
write.csv(x = temp, file = "Zoom1b-Forslund-SWE-T2D-Increasing-protein--with-functions.csv", quote = FALSE, row.names = FALSE )



## Zoom1c-Forslund-SWE-T2D-Increasing-protein.csv
temp <- read.csv(file = "Zoom1c-Forslund-SWE-T2D-Increasing-protein.csv", header = TRUE)
temp$zoom_group <- "1c"
temp$zoom_group_description <- "Increasing proteins"

dim(temp) # 5 27

# add functions?
temp$fxns <- NA
temp$superfocus_fxns <- NA
temp$subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  print(paste0("completed ",i))
}
temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
zoom_results[["1c. Increasing proteins"]] <- temp
write.csv(x = temp, file = "Zoom1c-Forslund-SWE-T2D-Increasing-protein--with-functions.csv", quote = FALSE, row.names = FALSE )




## Zoom3a-Forslund-SWE-T2D-Increasing-Lignin.csv
temp <- read.csv(file = "Zoom3a-Forslund-SWE-T2D-Increasing-Lignin.csv", header = TRUE)
temp$zoom_group <- "3a"
temp$zoom_group_description <- "Increasing lignins"

dim(temp) # 17 27

# add functions?
temp$fxns <- NA
temp$superfocus_fxns <- NA
temp$subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  print(paste0("completed ",i))
}
temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
zoom_results[["3a. Increasing lignins"]] <- temp
write.csv(x = temp, file = "Zoom3a-Forslund-SWE-T2D-Increasing-Lignin--with-functions.csv", quote = FALSE, row.names = FALSE )




## Zoom3c-Forslund-SWE-T2D-Decreasing-Arom2.csv
temp <- read.csv(file = "Zoom3c-Forslund-SWE-T2D-Decreasing-Arom2.csv", header = TRUE)
temp$zoom_group <- "3c"
temp$zoom_group_description <- "Decreasing aromatics"

dim(temp) # 24 27

# add functions?
temp$fxns <- NA
temp$superfocus_fxns <- NA
temp$subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  print(paste0("completed ",i))
}
temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
zoom_results[["3c. Decreasing aromatics"]] <- temp
write.csv(x = temp, file = "Zoom3c-Forslund-SWE-T2D-Decreasing-Arom2--with-functions.csv", quote = FALSE, row.names = FALSE )




## Zoom4-Forslund-SWE-T2D-DECREASING-lip

temp <- read.csv(file = "Zoom4-Forslund-SWE-T2D-DECREASING-lip.csv", header = TRUE)
temp$zoom_group <- "4"
temp$zoom_group_description <- "Decreasing lipids"

dim(temp) # 7 27

# add functions?
temp$fxns <- NA
temp$superfocus_fxns <- NA
temp$subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  print(paste0("completed ",i))
}
temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
zoom_results[["4. Decreasing lipids"]] <- temp
write.csv(x = temp, file = "Zoom4-Forslund-SWE-T2D-DECREASING-lip--with-functions.csv", quote = FALSE, row.names = FALSE )





#-------------------------


##########################
##########################
##########################
##########################


#### Test for selected compounds and hormones in restoration and t2d?
#-------------------------

lookfor <- as.data.frame( read_excel(path = "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/compounds-to-look-for.xlsx", sheet = 1, range = "A1:D44", col_names = TRUE) )
class(lookfor) # data.frame

lookfor$cpd_name
# [1] "cholecystokinin"                               "secretin"                                      "somatostatin"                                 
# [4] "motilin"                                       "ghrelin"                                       "glucagon-like peptide 1"                      
# [7] "glucose-dependent insulinotropic peptide"      "insulin-like peptide 5"                        "peptide YY"                                   
# [10] "gastrin"                                       "serotonin"                                     "neurotensin"                                  
# [13] "growth differentiation factor 15"              "fibroblast growth factor 19"                   "guanylin"                                     
# [16] "uroguanylin"                                   "oxyntomodulin"                                 "acetate"                                      
# [19] "propionate"                                    "butyrate"                                      "ammonia"                                      
# [22] "indole"                                        "trimethylamine N-oxide"                        "trimethylamine"                               
# [25] "p-cresyl sulfate"                              "indoxyl sulfate"                               "hydrogen sulfide"                             
# [28] "methane"                                       "menaquinone"                                   "cellulose"                                    
# [31] "xylan"                                         "pectin"                                        "amylopectin"                                  
# [34] "amylose"                                       "Carbohydrate response element binding protein" "glucosamine"                                  
# [37] "muramic acid"                                  "galactose"                                     "mannose"                                      
# [40] "arabinose"                                     "xylose"                                        "Carbon dioxide"                               
# [43] "adenosine triphosphate" 

## compound-level rel abundances in restoration and T2D

dat.cpd.res <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")
str(dat.cpd.res)
# 'data.frame':	125550 obs. of  7 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...
# $ log10_abun  : num  -3.362 -1.643 -3.341 0.711 -1.804 ...
# $ group       : Ord.factor w/ 5 levels "6"<"12"<"22"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ group_label : Ord.factor w/ 5 levels "6 yr"<"12 yr"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ ord_group   : num  3 3 3 3 3 3 3 3 3 3 ...

dat.cpd.t2d <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")
str(dat.cpd.t2d)
# 'data.frame':	1052845 obs. of  7 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260132" "ERR260132" "ERR260132" "ERR260132" ...
# $ cpd_rel_abun: num  0 5.647932 0 0 0.000186 ...
# $ log10_abun  : num  -8.482 0.752 -8.482 -8.482 -3.73 ...
# $ group       : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ group_label : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ ord_group   : int  2 2 2 2 2 2 2 2 2 2 ...

age_vec
# 6 12 22 31 UM 
# 1  2  3  4  5 

# select only Normal and T2D
unique(dat.cpd.t2d$group)
# [1] IGT         T2D met neg Normal      T2D met pos
# Levels: Normal < IGT < T2D met pos < T2D met neg
sel <- which(dat.cpd.t2d$group %in% c("T2D met neg", "Normal"))
dat.cpd.t2d <- dat.cpd.t2d[sel, ]

unique(dat.cpd.t2d$group_label)
# [1] T2D met- Normal  
# Levels: Normal < IGT < T2D met+ < T2D met-

# change order of factor levels
dat.cpd.t2d$group_label <- factor(dat.cpd.t2d$group_label, levels = c("T2D met-", "Normal"),
                                  ordered = TRUE)


length( unique(dat.cpd.t2d$cpd_id) ) # 7261
length( unique(dat.cpd.t2d$sample) ) # 76



# compound info is in here
names(df.comp)
# [1] "id"        "abbrev"    "name"      "form"      "OC_ratio"  "HC_ratio"  "NC_ratio"  "PC_ratio"  "NP_ratio"  "O_count"   "N_count"   "P_count"   "S_count"   "mass"     
# [15] "SC_ratio"  "MgC_ratio" "ZnC_ratio" "KC_ratio"  "CaC_ratio" "MnC_ratio" "FeC_ratio" "CoC_ratio" "CuC_ratio" "MoC_ratio"
head(df.comp$abbrev) # "h2o"   "atp"   "nad"   "nadh"  "nadph" "nadp"
head(df.comp$name) # "H2O"   "ATP"   "NAD"   "NADH"  "NADPH" "NADP" 
head(df.comp$form) # "H2O"           "C10H13N5O13P3" "C21H26N7O14P2" "C21H27N7O14P2" "C21H26N7O17P3" "C21H25N7O17P3"


## "Carbon dioxide"                               
sel.cpd <- which(df.comp$name == "CO2")
this_var <- "CO2"

df.comp[sel.cpd, ]
#          id abbrev name
# 11 cpd00011    co2  CO2

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00011")
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 4.1522, p-value = 3.292e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#      tau 
# 0.843525
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00011")
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "adenosine triphosphate" 

sel.cpd <- which(df.comp$name == "ATP")
this_var <- "ATP"

df.comp[sel.cpd, ]
#         id abbrev name          form
# 2 cpd00002    atp  ATP C10H13N5O13P3

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00002")
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -1.3166, p-value = 0.188
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#        tau 
# -0.2674591
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00002")
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "cholecystokinin" 
sel.cpd <- which(df.comp$name == "cholecystokinin")
sel.cpd <- grep(pattern = "*holecystokinin", x = df.comp$name)
sel.cpd <- grep(pattern = "CCK|cck", x = df.comp$abbrev)
df.comp[sel.cpd, ]
#             id abbrev   name
# 13205 cpd14835  CCK-8  CCK-8
# 13224 cpd14854 CCK-33 CCK-33
this_var <- "CCK-8"
this_var <- "CCK-33"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14835") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd14854") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14835") # empty
sel <- which(dat.cpd.t2d$cpd_id == "cpd14854")



## "secretin" (KEGG Compound C13523: Secretin; Vitrum)

sel.cpd <- which(df.comp$name == "secretin") # empty
sel.cpd <- grep(pattern = "*ecretin|*itrum", x = df.comp$name)
#sel.cpd <- grep(pattern = "CCK|cck", x = df.comp$abbrev)
df.comp[sel.cpd, ]
#             id abbrev   name             form
# 12827 cpd14242 Vitrum Vitrum C128H219N44O38R2
this_var <- "Secretin (Vitrum)"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14242") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14242") # empty



## "somatostatin"                                 

sel.cpd <- which(df.comp$name == "Somatostatin") # empty
sel.cpd <- grep(pattern = "*omatostatin", x = df.comp$name)
df.comp[sel.cpd, ]
#             id         abbrev           name
# 13113 cpd14743 Somatostatin-2 Somatostatin-2
# 13114 cpd14744 Somatostatin-1 Somatostatin-1
this_var <- "Somatostatin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14743") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd14744") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14743") # empty
sel <- which(dat.cpd.t2d$cpd_id == "cpd14744") # empty




## "motilin" 
sel.cpd <- which(df.comp$name == "Motilin") # ok
#sel.cpd <- grep(pattern = "*otilin", x = df.comp$name)
df.comp[sel.cpd, ]
#             id  abbrev    name
# 13138 cpd14768 Motilin Motilin
this_var <- "Motilin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14768") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14768") # empty


## "ghrelin"

sel.cpd <- which(df.comp$name == "Ghrelin") # ok
#sel.cpd <- grep(pattern = "*hrelin", x = df.comp$name)
df.comp[sel.cpd, ]
#             id  abbrev    name
# 13117 cpd14747 Ghrelin Ghrelin
this_var <- "Ghrelin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14747") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14747") # empty



## "glucagon-like peptide 1" ( or GLP-1;GLP1)                
sel.cpd <- which(df.comp$name == "GLP-1") # ok
sel.cpd <- grep(pattern = "GLP|*lucagon-like peptide 1", x = df.comp$name)
df.comp[sel.cpd, ]
#             id                  abbrev                    name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio
# 13139 cpd14769 Glucagon-like peptide 1 Glucagon-like peptide 1
this_var <- "GLP-1"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14769") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14769") # empty



## "glucose-dependent insulinotropic peptide"
sel.cpd <- which(df.comp$name == "Glucose-dependent insulinotropic peptide") # empty
sel.cpd <- grep(pattern = "GIP|*lucose-dependent insulinotropic peptide", x = df.comp$name) # ok - GIP
df.comp[sel.cpd, ]
#             id abbrev name
# 13012 cpd14636    GIP  GIP
this_var <- "GIP"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14636") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14636") # empty



## "insulin-like peptide 5"
sel.cpd <- which(df.comp$name == "Glucose-dependent insulinotropic peptide") # empty
sel.cpd <- grep(pattern = "INSl5|*nsulin-like peptide 5", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#             id                 abbrev                   name
# 13268 cpd14898 Insulin-like peptide 5 Insulin-like peptide 5
this_var <- "INSL5"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14898") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14898") # empty


## "peptide YY"                                   
sel.cpd <- which(df.comp$name == "Peptide YY") # ok
sel.cpd <- grep(pattern = "PYY|*eptide YY", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#             id     abbrev       name
# 13209 cpd14839 Peptide YY Peptide YY
this_var <- "PYY"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14839") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14839") # empty


## "gastrin"  
sel.cpd <- which(df.comp$name == "Gastrin") # ok
#sel.cpd <- grep(pattern = "|*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#             id  abbrev    name
# 13204 cpd14834 Gastrin Gastrin
this_var <- "Gastrin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14834") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14834") # empty


## "serotonin" 
sel.cpd <- which(df.comp$name == "Serotonin") # ok
#sel.cpd <- grep(pattern = "*erotonin", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#           id    abbrev      name      form
# 568 cpd00579 Serotonin Serotonin C10H13N2O
this_var <- "Serotonin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00579") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = -1.8229, p-value = 0.06831
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# -0.370328
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00579") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "neurotensin"                                  
sel.cpd <- which(df.comp$name == "Neurotensin") # ok
#sel.cpd <- grep(pattern = "|*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#             id      abbrev        name
# 11780 cpd11960 Neurotensin Neurotensin
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11960") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11960") # empty


## "growth differentiation factor 15" 
sel.cpd <- which(df.comp$name == "Growth differentiation factor 15") # empty
sel.cpd <- grep(pattern = "GDF15|GDF-15", x = df.comp$name) # empty
df.comp[sel.cpd, ]
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "") # empty


## "fibroblast growth factor 19"
sel.cpd <- which(df.comp$name == "Fibroblast growth factor 19") # empty
sel.cpd <- grep(pattern = "FGF|*ibroblast growth", x = df.comp$name) # empty
df.comp[sel.cpd, ]
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "") # empty


## "guanylin"                                     
sel.cpd <- which(df.comp$name == "Guanylin") # ok
#sel.cpd <- grep(pattern = "|*", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#             id   abbrev     name
# 13099 cpd14729 Guanylin Guanylin
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14729") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14729") # empty


## "uroguanylin"
sel.cpd <- which(df.comp$name == "Uroguanylin") # ok
#sel.cpd <- grep(pattern = "|*roguanylin", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#             id      abbrev        name
# 13098 cpd14728 Uroguanylin Uroguanylin
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14728") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14728") # empty


## "oxyntomodulin"
sel.cpd <- which(df.comp$name == "Oxyntomodulin") # ok
#sel.cpd <- grep(pattern = "|*", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#             id        abbrev          name
# 16751 cpd19466 Oxyntomodulin Oxyntomodulin
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd19466") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd19466") # empty


## "acetate"
sel.cpd <- which(df.comp$name == "Acetate") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#          id abbrev    name   form
# 29 cpd00029     ac Acetate C2H3O2
this_var <- "Acetate"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00029") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 3.5446, p-value = 0.0003932
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.7200823
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00029") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "propionate"
sel.cpd <- which(df.comp$name == "Propionate") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#           id abbrev       name   form
# 140 cpd00141    ppa Propionate C3H5O2
this_var <- "Propionate"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00141") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 2.1268, p-value = 0.03344
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.4320494 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00141") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "butyrate"
sel.cpd <- which(df.comp$name == "Butyrate") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#           id abbrev     name   form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 210 cpd00211    but Butyrate C4H7O2
this_var <- "Butyrate"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00211") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 3.2408, p-value = 0.001192
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#      tau 
# 0.658361
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00211") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "ammonia"
#sel.cpd <- which(df.comp$name == "Ammonia") # empty
#sel.cpd <- grep(pattern = "NH3|*mmonia", x = df.comp$name) # ok
sel.cpd <- which(df.comp$name == "NH3") # ok
df.comp[sel.cpd, ]
#          id abbrev name
# 13 cpd00013    nh4  NH3
this_var <- "Ammonia"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00013") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = -1.7217, p-value = 0.08513
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#        tau 
# -0.3497543
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)


# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by ord_group
# Kruskal-Wallis chi-squared = 8.9, df = 4, p-value = 0.06365


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00013") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "indole" (KEGG C00463: 2,3-Benzopyrrole)
#sel.cpd <- which(df.comp$name == "indole") # empty
#sel.cpd <- which(df.comp$name == "benzopyrrole") # empty
sel.cpd <- which(df.comp$abbrev == "indole" & df.comp$form == "C8H7N") # ok
#sel.cpd <- grep(pattern = "Indole", x = df.comp$abbrev) # empty
#sel.cpd <- grep(pattern = "*ndole", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#           id abbrev  name  form 
# 356 cpd00359 indole indol C8H7N

this_var <- "Indole"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00359") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 0.10127, p-value = 0.9193
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.02057378
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by ord_group
# Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00359") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "trimethylamine N-oxide"
#sel.cpd <- which(df.comp$name == "TMAO") # empty
sel.cpd <- which(df.comp$form == "C3H9NO") # ok
#sel.cpd <- grep(pattern = "*rimethylamine", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#           id abbrev     name   form
# 797 cpd00811   tmao (CH3)3NO C3H9NO

this_var <- "TMAO"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00811") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 2.8357, p-value = 0.004573
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.5760658
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00811") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "trimethylamine"                               
sel.cpd <- which(df.comp$abbrev == "tma") # ok
#sel.cpd <- which(df.comp$form == "C3H9NO") # ok
#sel.cpd <- grep(pattern = "*rimethylamine", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#           id abbrev    name   form 
# 437 cpd00441    tma (CH3)3N C3H10N

this_var <- "TMA"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00441") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 2.7344, p-value = 0.006249
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.5554921
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00441") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "p-cresyl sulfate" - instead use precursor: p-Cresol
sel.cpd <- which(df.comp$abbrev == "p-Cresol") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*resol", x = df.comp$abbrev) # empty
df.comp[sel.cpd, ]
# id   abbrev     name  form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass
# 1022 cpd01042 p-Cresol p-Cresol C7H8O

this_var <- "p-Cresol"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd01042") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 0.20255, p-value = 0.8395
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.04114756
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd01042") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "hydrogen sulfide"                             
sel.cpd <- which(df.comp$name == "H2S") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$abbrev) # empty
df.comp[sel.cpd, ]
#           id abbrev name
# 237 cpd00239    h2s  H2S

this_var <- "H2S"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00239") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.0382, p-value = 0.00238
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.6172134
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00239") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "methane"
sel.cpd <- which(df.comp$name == "Methane") # ok
#sel.cpd <- which(df.comp$form == "CH4") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$abbrev) # empty
df.comp[sel.cpd, ]
#            id abbrev    name form
# 1004 cpd01024  metha Methane  CH4

this_var <- "CH4"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd01024") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 2.3293, p-value = 0.01984
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.4731969 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd01024") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "menaquinone"
sel.cpd <- which(df.comp$abbrev == "Menaquinones") # ok
#sel.cpd <- which(df.comp$form == "CH4") # ok
#sel.cpd <- grep(pattern = "*enaquinones", x = df.comp$abbrev) # multiple
#sel.cpd <- grep(pattern = "C16H16O2*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id       abbrev         name      form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count
# 24771 cpd27501 Menaquinones Menaquinones C16H15O2R

this_var <- "Menaquinones"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd27501") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 1.9242, p-value = 0.05433
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.3909018
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by ord_group
# Kruskal-Wallis chi-squared = 7.6333, df = 4, p-value = 0.106


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd27501") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "cellulose"                                    
sel.cpd <- which(df.comp$name == "Cellulose") # ok
#sel.cpd <- which(df.comp$form == "CH4") # ok
sel.cpd <- grep(pattern = "*ellulose", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "C16H16O2*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id    abbrev      name      form
# 11571 cpd11746 Cellulose Cellulose C6H10O5R2

this_var <- "Cellulose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11746") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -0.50637, p-value = 0.6126
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.1028689 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11746") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "xylan"
sel.cpd <- which(df.comp$name == "Xylan") # ok
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "xylan", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
# id                                    abbrev                                      name        form  OC_ratio  HC_ratio
# 10924 cpd11070 3beta-Hydroxylanostane-7,11-dione acetate 3beta-Hydroxylanostane-7,11-dione acetate    C32H52O4 0.1250000 1.6250000
# 11789 cpd11970                              Arabinoxylan                              Arabinoxylan  C10H16O8R2 0.8000000 1.6000000
# 12068 cpd12254                     Glucuronoarabinoxylan                     Glucuronoarabinoxylan        null       NaN       NaN
# 12217 cpd12408                            Glucuronoxylan                            Glucuronoxylan C21H31O18R2 0.8571429 1.4761905
# 12354 cpd12550                  4-O-methylglucuronoxylan                  4-O-methylglucuronoxylan C22H33O18R2 0.8181818 1.5000000
# 17230 cpd19945                     11-Deoxylandomycinone                     11-Deoxylandomycinone    C19H14O5 0.2631579 0.7368421
# 19578 cpd22295                               Acetylxylan                               Acetylxylan C27H39O21R3 0.7777778 1.4444444
# 19670 cpd22387                              Arabinoxylan                              Arabinoxylan C35H56O29R2 0.8285714 1.6000000
# 24444 cpd27172    Glucuronoarabinoxylan-Oligosaccharides    Glucuronoarabinoxylan-Oligosaccharides C27H41O23R2 0.8518519 1.5185185
# 24445 cpd27173                    Glucuronoarabinoxylans                    Glucuronoarabinoxylans C52H81O43R2 0.8269231 1.5576923
# 24447 cpd27175           Glucuronoxylan-Oligosaccharides           Glucuronoxylan-Oligosaccharides C22H33O19R2 0.8636364 1.5000000
# 24448 cpd27176                           Glucuronoxylans                           Glucuronoxylans C37H57O31R2 0.8378378 1.5405405
# 25918 cpd28650                      (1, 4-beta-D-xylan)n                      (1, 4-beta-D-xylan)n 

sel.cpd <- which(df.comp$name == "Arabinoxylan") # ok
df.comp[sel.cpd, ]
#             id       abbrev         name        form
# 11789 cpd11970 Arabinoxylan Arabinoxylan  C10H16O8R2 - NOT PRESENT IN DATA
# 19670 cpd22387 Arabinoxylan Arabinoxylan C35H56O29R2


this_var <- "Arabinoxylan"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd22387") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -2.228, p-value = 0.02588
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.4526232
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd22387") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## "Hemicellulose"  
# Hemicelluloses include xyloglucans, xylans, mannans and glucomannans
sel.cpd <- which(df.comp$name == "Hemicellulose") # ok
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ellulose", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id        abbrev          name form
# 27115 cpd29869 Hemicellulose Hemicellulose null
this_var <- "Hemicellulose"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd29869") # empty - Hemicellulose
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd29869") # empty


sel.cpd <- which(df.comp$name == "Xyloglucan") # ok
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*yloglucan", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id     abbrev       name        form  
# 11577 cpd11752 Xyloglucan Xyloglucan        null  - Not in data
# 25623 cpd28355 Xyloglucan Xyloglucan C39H61O33R5
this_var <- "Xyloglucan"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd28355") # 15

sel <- which(dat.cpd.res$cpd_id == "cpd11752")



x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 1.6204, p-value = 0.1052
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.3291805
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd28355") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")







## "pectin"
sel.cpd <- which(df.comp$name == "Pectin" & df.comp$form == "C26H34O24R2") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name        form
# 11434 cpd11601 pectin Pectin C26H34O24R2
this_var <- "Pectin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11601") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -1.9242, p-value = 0.05433
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.3909018 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11601") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "amylopectin"                                  
sel.cpd <- which(df.comp$name == "Amylopectin") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#           id      abbrev        name      form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count
# 262 cpd00265 Amylopectin Amylopectin C30H52O26
this_var <- "Amylopectin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00265") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 1.6204, p-value = 0.1052
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.3291805 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00265") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "amylose"
sel.cpd <- which(df.comp$name == "Amylose") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev    name      form
# 11560 cpd11735 14glun Amylose C6H10O5R2
this_var <- "Amylose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11735") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 1.2153, p-value = 0.2243
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.2468854
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11735") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "glucosamine" - instead use Chitin because glucosamine derives from Chitin
# D-Glucosamine is made naturally in the form of glucosamine-6-phosphate
#sel.cpd <- which(df.comp$name == "Glucosamine") # ok - but not in samples
#sel.cpd <- grep(pattern = "D-glucosamine-6-phosphate", x = df.comp$name) # D-Glucosamine is made naturally in the form of glucosamine-6-phosphate
sel.cpd <- grep(pattern = "*D-glucosamine", x = df.comp$name)
sel.cpd <- which(df.comp$name == "alpha-D-glucosamine 6-phosphate") 
# "alpha-D-glucosamine" "alpha-D-glucosamine 6-phosphate"
#sel.cpd <- grep(pattern = "*mino-2-deoxy-glucose", x = df.comp$name)
#sel.cpd <- which(df.comp$name == "2-amino-2-deoxy-glucose") # ok - but not in samples
#sel.cpd <- which(df.comp$form == "C6H13NO5") # ok - but not in samples

unique(df.comp$name[sel.cpd])

df.comp[sel.cpd, ]
#             id                          abbrev                            name      form 
# 21178 cpd23898 alpha-D-glucosamine 6-phosphate alpha-D-glucosamine 6-phosphate C6H13NO8P


this_var <- "D-glucosamine 6-phosphate"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd23898") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -2.5318, p-value = 0.01135
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.5143445
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd23898") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "chitin"
sel.cpd <- which(df.comp$name == "Chitin") #
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*lucosamine", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name       form 
# 11508 cpd11683 chitin Chitin C8H13NO5R2
this_var <- "Chitin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11683") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.8484, p-value = 0.0001189
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.7818036 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11683") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "muramic acid" - occurs naturally as N-acetylmuramic acid in peptidoglycan, whose primary function is a structural component of many typical bacterial cell walls
sel.cpd <- which(df.comp$name == "Muramic acid") # ok but none in samples
sel.cpd <- which(df.comp$name == "N-Acetylmuramic acid") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*cetylmuramic", x = df.comp$name) # multiple - use this N-acetylmuramic acid-1-phosphate
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id                           abbrev                                                  name        form
# 13716 cpd15396                            anhgm N-Acetyl-D-glucosamine(anhydrous)N-Acetylmuramic acid C19H29N2O12
# 23540 cpd26263 N-acetylmuramic acid-1-phosphate                      N-acetylmuramic acid-1-phosphate C11H17NO11P
this_var <- "N-acetylmuramic acid-1-phosphate"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd26263") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 0.60764, p-value = 0.5434
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.1234427 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd26263") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "galactose"
sel.cpd <- which(df.comp$name == "Galactose") # 
sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*", x = df.comp$name) # multiple - use this N-acetylmuramic acid-1-phosphate
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#           id  abbrev      name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 108  cpd00108     gal Galactose C6H12O6        1        2        0        0      NaN       6       0       0       0  180        0         0         0        0         0
# 1090 cpd01112 cbs_337 Galactose C6H12O6
this_var <- "Galactose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00108") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.0382, p-value = 0.00238
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.6172134
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00108") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "mannose"                                      

sel.cpd <- which(df.comp$name == "D-Mannose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*annose$", x = df.comp$name) # multiple - use this D-mannose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#           id abbrev      name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 138 cpd00138    man D-Mannose C6H12O6
this_var <- "D-Mannose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00138") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -2.9369, p-value = 0.003315
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.5966396
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00138") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "arabinose"

sel.cpd <- which(df.comp$name == "L-Arabinose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*rabinose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#           id abbrev        name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio
# 223 cpd00224 arab-L L-Arabinose C5H10O5
this_var <- "L-Arabinose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00224") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -4.1522, p-value = 3.292e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.843525
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00224") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "xylose"

sel.cpd <- which(df.comp$name == "Xylose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*Xylose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id             abbrev               name        form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count
# 153   cpd00154              xyl-D             Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 1047  cpd01068           L-Xylose           L-Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 1397  cpd01422            cll_392      beta-D-Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 1460  cpd01487     alpha-D-Xylose     alpha-D-Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 24104 cpd26831           D-Xylose           D-Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 25629 cpd28361 Xyloglucans-Xylose Xyloglucans-Xylose C39H64O33R2 0.8461538 1.641026        0        0      NaN      33       0
# 26016 cpd28748         (+)-Xylose         (+)-Xylose     C5H10O5 
this_var <- "Xylose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00154") # Xylose
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.7471, p-value = 0.0001789
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.7612299 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00154") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## glucose
sel.cpd <- which(df.comp$name == "D-Glucose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*lucose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
# id    abbrev      name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 27    cpd00027     glc-D D-Glucose C6H12O6        1        2        0        0      NaN       6       0       0       0  180        0         0         0        0         0
# 24094 cpd26821 D-Glucose D-Glucose C6H12O6
this_var <- "D-Glucose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00027") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.9497, p-value = 7.825e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.8023774
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00027") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Fructose
sel.cpd <- which(df.comp$name == "D-Fructose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ructose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#          id abbrev       name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 82 cpd00082    fru D-Fructose C6H12O6
this_var <- "D-Fructose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00082") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.7471, p-value = 0.0001789
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.7612299 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00082") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Sucrose
sel.cpd <- which(df.comp$name == "Sucrose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ucrose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#          id abbrev    name      form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 76 cpd00076   sucr Sucrose C12H22O11
this_var <- "Sucrose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00076") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.8484, p-value = 0.0001189
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.7818036
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00076") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## leptin
sel.cpd <- which(df.comp$name == "Leptin") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*eptin$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count  mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 16743 cpd19458 Leptin Leptin
this_var <- "Leptin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd19458") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd19458") # empty




## lignin

sel.cpd <- which(df.comp$name == "Lignin") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ucrose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 12548 cpd12745 Lignin Lignin null
this_var <- "Lignin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd12745") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 3.9497, p-value = 7.825e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.8023774
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd12745") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## palmitoleic acid (monounsaturated fatty acid)
# Where the balance of the balance of palmitic acid (saturated fatty acid) and palmitoleic acid (monounsaturated fatty acid) can affect proinflammatory immune response [Tsai et al 2021 - https://www.plefa.com/article/S0952-3278(21)00033-8/fulltext ]

sel.cpd <- which(df.comp$name == "Palmitoleic acid") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*almitoleic", x = df.comp$name) # 
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#            id           abbrev             name     form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio
# 5179 cpd05274 Palmitoleic acid Palmitoleic acid C16H29O2
this_var <- "Palmitoleic acid"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd05274") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -2.3293, p-value = 0.01984
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.4731969 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd05274") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## palmitic acid (saturated fatty acid), in KEGG: Hexadecanoic acid / Hexadecanoate / Hexadecylic acid / Palmitate / Cetylic acid

sel.cpd <- which(df.comp$name == "Palmitate") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*exadecanoate$", x = df.comp$name) # 
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#           id abbrev      name     form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 213 cpd00214   hdca Palmitate C16H31O2
this_var <- "Palmitate (Palmitic acid)"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00214") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -0.70892, p-value = 0.4784
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.1440165 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00214") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")





## vaccenic acid
#Plasma and tissue levels of palmitoleic acid and vaccenic acid are largely produced through de novo lipogenesis and regulated by numerous hormones including insulin [7]. Critically, de novo lipogenic dysregulation results in aberrant fatty acid levels [7], which, in turn, may predict future clinical outcomes. - https://www.sciencedirect.com/science/article/abs/pii/S1262363619301776?via%3Dihub 

sel.cpd <- which(df.comp$name == "Vaccenic acid") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*$", x = df.comp$name) # 
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#            id        abbrev          name     form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio
# 5184 cpd05279 Vaccenic acid Vaccenic acid C18H33O2
this_var <- "Vaccenic acid"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd05279") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd05279") # empty




## Glycogen (animal starch) - Glycogen is a multibranched polysaccharide of glucose that serves as a form of energy storage in animals,[2] fungi, and bacteria

sel.cpd <- which(df.comp$name == "Glycogen") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*$", x = df.comp$name) # 
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#           id   abbrev     name      form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 154 cpd00155 glycogen Glycogen C24H42O21
this_var <- "Glycogen"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00155") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 2.9369, p-value = 0.003315
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.5966396
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00155") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Amylum (plant starch)

sel.cpd <- which(df.comp$name == "Starch") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*$", x = df.comp$name) # 
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name        form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count    mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 11482 cpd11657 starch Starch C12H20O10R2 0.8333333 1.666667        0        0      NaN      10       0       0       0 9.9e+02        0         0         0        0         0
# 25462 cpd28193 Starch Starch        null
this_var <- "Starch"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11657") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 0.91147, p-value = 0.3621
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.185164
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11657") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




### fructan

sel.cpd <- which(df.comp$name == "Fructan") # empty
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ructan", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#            id               abbrev                 name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count
# 24363 cpd27091             Fructans             Fructans null      NaN      NaN      NaN      NaN      NaN       0       0       0
# 26031 cpd28763 (2,1-beta-D-fructan) (2,1-beta-D-fructan) null      NaN      NaN      NaN      NaN      NaN       0       0       0
# 26589 cpd29322   2,6-beta-D-fructan   2,6-beta-D-fructan null
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd27091") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd28763") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd29322") # empty




### inulin

                       
sel.cpd <- which(df.comp$name == "Inulin") # ok
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ructosyl", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name       form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass
# 11435 cpd11602 Inulin Inulin  C24H42O21 0.8750000 1.750000        0        0      NaN      21       0       0       0  342
# 24584 cpd27312 Inulin Inulin C18H31O15R
this_var <- "Inulin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11602") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd27312") # empty




### galactan

                            
sel.cpd <- which(df.comp$name == "Galactan") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id   abbrev     name        form
# 12578 cpd12777 Galactan Galactan C12H20O11R2
this_var <- "Galactan"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd12777") # empty




## TESTING T2D INDEX ??

t2d.zoom1.decreasing.proteins
# [1] "cpd11499" "cpd11503" "cpd11507" "cpd11511" "cpd11524" "cpd11528" "cpd11532" "cpd11536" "cpd11549" "cpd11553" "cpd11557" "cpd11561" "cpd11475"
# [14] "cpd11465" "cpd11469" "cpd11473" "cpd11498" "cpd11502" "cpd11506" "cpd11510" "cpd11514" "cpd11518" "cpd11523" "cpd11527" "cpd11531" "cpd11535"
# [27] "cpd11539" "cpd11543" "cpd11548" "cpd11552" "cpd11556" "cpd11560" "cpd11564" "cpd11568" "cpd11572"
t2d.zoom3b.decreasing.arom1
# "cpd00430" "cpd03336" "cpd22764" "cpd27851" "cpd27852" "cpd01554" "cpd01722" "cpd06913" "cpd34699" "cpd00435"
t2d.zoom2B.increasing.carbs
# "cpd00224" "cpd03198" "cpd03200" "cpd00259" "cpd00076" "cpd00082"


## sugars - D-Fructose: "cpd00082", Sucrose : "cpd00076" , L-Arabinose: "cpd00224"
#sugars <- c("cpd00082", "cpd00076", "cpd00224")
sugars <- t2d.zoom2B.increasing.carbs

## starches - starch: "cpd11657", amylopectin: "cpd00265", amylose: "cpd11735"
starches <- c("cpd11657", "cpd00265", "cpd11735")
## ammonia "cpd00013"
ammonia <- "cpd00013"
## acps - acyl carrier proteins - list: t2d.zoom1.decreasing.proteins
acps <- t2d.zoom1.decreasing.proteins
## fiber - lignin "cpd12745" + pectin "cpd11601"
fiber <- c("cpd12745", "cpd11601")
#fiber <- "cpd12745"

## tmao "cpd00811"
tmao <- "cpd00811"

aroms <- t2d.zoom3b.decreasing.arom1


temp <- dat.cpd.t2d
head(temp)
# cpd_id    sample cpd_rel_abun log10_abun       group group_label ord_group
# 50828 cpd24620 ERR260139 0.0000000000 -8.4820250 T2D met neg    T2D met-         4
# 50829 cpd00001 ERR260139 4.9744050062  0.6967411 T2D met neg    T2D met-         4
# 50830 cpd25681 ERR260139 0.0000000000 -8.4820250 T2D met neg    T2D met-         4
# 50831 cpd01501 ERR260139 0.0000000000 -8.4820250 T2D met neg    T2D met-         4
# 50832 cpd02597 ERR260139 0.0001838302 -3.7355832 T2D met neg    T2D met-         4
# 50833 cpd00851 ERR260139 0.0012068230 -2.9183564 T2D met neg    T2D met-         4

unique_samps <- unique(temp$sample)

df.calc <- data.frame(sample = unique_samps, group_label=NA, sugar=NA, starch=NA, nh3=NA, acp=NA, fib=NA, tmao=NA, arom = NA)
for (i in 1:length(unique_samps)) {
  #i<-22
  this_samp <- unique_samps[i]
  sel <- which(temp$sample == this_samp)
  df.calc$group_label[i] <- as.character( unique(temp$group_label[sel]) )
  
  sel.sugar <- which(temp$cpd_id %in% sugars & temp$sample == this_samp)
  df.calc$sugar[i] <- sum( temp$cpd_rel_abun[sel.sugar], na.rm = TRUE )
  
  sel.starch <- which(temp$cpd_id %in% starches & temp$sample == this_samp)
  df.calc$starch[i] <- sum( temp$cpd_rel_abun[sel.starch], na.rm = TRUE )
  
  sel.nh3 <- which(temp$cpd_id == ammonia & temp$sample == this_samp)
  df.calc$nh3[i] <- sum( temp$cpd_rel_abun[sel.nh3], na.rm = TRUE )
  
  sel.acp <- which(temp$cpd_id %in% acps & temp$sample == this_samp)
  df.calc$acp[i] <- sum( temp$cpd_rel_abun[sel.acp], na.rm = TRUE )
  
  sel.fib <- which(temp$cpd_id %in% fiber & temp$sample == this_samp)
  #sel.fib <- which(temp$cpd_id == fiber & temp$sample == this_samp)
  df.calc$fib[i] <- sum( temp$cpd_rel_abun[sel.fib], na.rm = TRUE )
  
  sel.tmao <- which(temp$cpd_id == tmao & temp$sample == this_samp)
  df.calc$tmao[i] <- sum( temp$cpd_rel_abun[sel.tmao], na.rm = TRUE )
  
  sel.arom <- which(temp$cpd_id %in% aroms & temp$sample == this_samp)
  df.calc$arom[i] <- sum( temp$cpd_rel_abun[sel.arom], na.rm = TRUE )
  
  print(paste0("completed ",i))
  
}


## zero replace and scale variables between 1-100

df.calc.scaled <- df.calc
vars <- names(df.calc.scaled)[-c(1,2)]

for (i in 1:length(vars)) {
  #i<-1
  values <- df.calc.scaled[ ,vars[i]]
  scaled_values <- scales::rescale(values, to = c(1, 100))
  df.calc.scaled[ ,vars[i]] <- scaled_values
  print(paste0("completed ",i))
}


df.calc.scaled$index1 <- ( df.calc.scaled$sugar * df.calc.scaled$starch  ) / ( df.calc.scaled$acp ) # 0.009
df.calc.scaled$index1 <- ( df.calc.scaled$sugar * df.calc.scaled$starch * df.calc.scaled$nh3 ) / ( df.calc.scaled$acp ) # P = 0.003
df.calc.scaled$index1 <- ( df.calc.scaled$sugar * df.calc.scaled$starch * df.calc.scaled$nh3 ) / ( df.calc.scaled$acp * df.calc.scaled$fib  ) # P = 0.005

df.calc.scaled$index1 <- ( df.calc.scaled$sugar * df.calc.scaled$starch * df.calc.scaled$nh3 ) / ( df.calc.scaled$acp * df.calc.scaled$tmao ) # P = 0.004
df.calc.scaled$index1 <- ( df.calc.scaled$sugar * df.calc.scaled$starch * df.calc.scaled$nh3 ) / ( df.calc.scaled$acp * df.calc.scaled$arom ) # W = 997, P = 0.001

df.calc.scaled$index1 <- ( df.calc.scaled$sugar * df.calc.scaled$starch * df.calc.scaled$nh3 ) / ( df.calc.scaled$acp * df.calc.scaled$arom * df.calc.scaled$fib ) # W = 992, P = 0.001
df.calc.scaled$index1 <- ( df.calc.scaled$sugar * df.calc.scaled$starch * df.calc.scaled$nh3 ) / ( df.calc.scaled$acp * df.calc.scaled$arom * df.calc.scaled$tmao ) # W = 984, P = 0.002

this_var <- "Energy harvest / Energy management"
df.calc.scaled$index1 <- ( df.calc.scaled$sugar * df.calc.scaled$starch * df.calc.scaled$nh3 ) / ( df.calc.scaled$acp * df.calc.scaled$arom * df.calc.scaled$tmao * df.calc.scaled$fib ) # W = 1001, P = <0.001
# > wmw.test$p.value
# [1] 0.0009950782

this_var <- "Energy harvest"
df.calc.scaled$index1 <- ( df.calc.scaled$sugar * df.calc.scaled$starch * df.calc.scaled$nh3 ) # W = 990 P = 0.001
# > wmw.test$p.value
# [1] 0.001488098

this_var <- "1 / Energy management"
df.calc.scaled$index1 <- 1 / ( df.calc.scaled$acp * df.calc.scaled$arom * df.calc.scaled$tmao * df.calc.scaled$fib ) # W = 970, P = 0.003
# > wmw.test$p.value
# [1] 0.002973089

df.calc.scaled$log10_index1 <- log10( df.calc.scaled$index1 )


df.calc.scaled$group_label <- factor(df.calc.scaled$group_label, levels = c("T2D met-", "Normal"))

# T2D
subsel.t2d <- which(df.calc.scaled$group_label == "T2D met-")
subsel.norm <- which(df.calc.scaled$group_label == "Normal")

df.calc.scaled
df.calc.scaled[subsel.t2d ,]
df.calc.scaled[subsel.norm ,]

#this_var <- "Index1"
# x <- df.calc.scaled$index1[ subsel.t2d ] # n = 33
# y <- df.calc.scaled$index1[ subsel.norm ] # n = 43
x <- df.calc.scaled$log10_index1[ subsel.t2d ] # n = 33
y <- df.calc.scaled$log10_index1[ subsel.norm ] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = df.calc.scaled , aes(x=group_label, y = log10_index1))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  
  #xlab("Diagnosis")+ ylab("Log10 ratio of scaled CPP values")+
  xlab("Diagnosis")+ ylab("Log10 multiplied scaled CPP values")+
  #xlab("Diagnosis")+ ylab("Log10 1/multiplied scaled CPP values")+
  
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
grid.text(label = "(c)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

dev.print(tiff, file = paste0(workdir,"/plots/","T2D-Index-a-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 10, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","T2D-Index-b-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 10, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","T2D-Index-c-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v2.tiff"), width = 8, height = 10, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------


#### Sample summary stats
#-------------------------

## Restoration

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...
# $ log10_abun  : num  -3.362 -1.643 -3.341 0.711 -1.804 ...
# $ group       : Ord.factor w/ 5 levels "6"<"12"<"22"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ group_label : Ord.factor w/ 5 levels "6 yr"<"12 yr"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ ord_group   : num  3 3 3 3 3 3 3 3 3 3 ...

length( unique(dat.cpd.collate$cpd_id) ) # 8370
length( unique(dat.cpd.collate$cpd_id[ dat.cpd.collate$cpd_rel_abun > 0] ) ) # 8370
length( unique(dat.cpd.collate$sample) ) # 15

data_in <- dat.cpd.collate

head(data_in)
# cpd_id sample cpd_rel_abun log10_abun group group_label ord_group
# 1 cpd25681    20C 0.0004345842 -3.3619261    22       22 yr         3
# 2 cpd02597    20C 0.0227382574 -1.6432428    22       22 yr         3
# 3 cpd24620    20C 0.0004564208 -3.3406346    22       22 yr         3
# 4 cpd00001    20C 5.1354106008  0.7105752    22       22 yr         3
# 5 cpd01501    20C 0.0157048650 -1.8039658    22       22 yr         3
# 6 cpd00851    20C 0.0121620567 -1.9149930    22       22 yr         3

dim(data_in) # 125550      7
8370*15 # 125550

unique_samps <- unique(data_in$sample)

no_compounds <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  print(paste0("completed ",i))
}

mean(no_compounds) # 7776.867
sd(no_compounds) # 92.37878

phy <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 30125 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 30125 taxa by 4 taxonomic ranks ]

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 15
NonZeroFxns

mean(NonZeroFxns) # 20327.93
sd(NonZeroFxns) # 767.4159
#rm(NonZeroFxns)

# check
unique_samps <- sample_names(phy)

unique_samps[1] # "20C"
a <- prune_samples(unique_samps[1], phy)
min(taxa_sums(a)) # 0
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a

unique_samps[15] # "30C"
a <- prune_samples(unique_samps[15], phy)
min(taxa_sums(a)) # 0
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a






## T2D

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-Forslund-SWE-T2D.rds")

str(dat.cpd.collate)
# 'data.frame':	1052845 obs. of  7 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260132" "ERR260132" "ERR260132" "ERR260132" ...
# $ cpd_rel_abun: num  0 5.647932 0 0 0.000186 ...
# $ log10_abun  : num  -8.482 0.752 -8.482 -8.482 -3.73 ...
# $ group       : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ group_label : Ord.factor w/ 4 levels "Normal"<"IGT"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ ord_group   : int  2 2 2 2 2 2 2 2 2 2 ...

# select only Normal and T2D
unique(dat.cpd.collate$group)
# [1] IGT         T2D met neg Normal      T2D met pos
# Levels: Normal < IGT < T2D met pos < T2D met neg

sel <- which(dat.cpd.collate$group %in% c("T2D met neg", "Normal"))

dat.cpd.collate.T2DNORM <- dat.cpd.collate[sel, ]

length( unique(dat.cpd.collate.T2DNORM$cpd_id) ) # 7261
length( unique(dat.cpd.collate.T2DNORM$cpd_id[ dat.cpd.collate.T2DNORM$cpd_rel_abun > 0] ) ) # 7031
length( unique(dat.cpd.collate.T2DNORM$sample) ) # 76

data_in <- dat.cpd.collate.T2DNORM

head(data_in)
#         cpd_id    sample cpd_rel_abun log10_abun       group group_label ord_group
# 50828 cpd24620 ERR260139 0.0000000000 -8.4820250 T2D met neg    T2D met-         4
# 50829 cpd00001 ERR260139 4.9744050062  0.6967411 T2D met neg    T2D met-         4
# 50830 cpd25681 ERR260139 0.0000000000 -8.4820250 T2D met neg    T2D met-         4
# 50831 cpd01501 ERR260139 0.0000000000 -8.4820250 T2D met neg    T2D met-         4
# 50832 cpd02597 ERR260139 0.0001838302 -3.7355832 T2D met neg    T2D met-         4
# 50833 cpd00851 ERR260139 0.0012068230 -2.9183564 T2D met neg    T2D met-         4

dim(data_in) # 551836      7
7261*76 # 551836

unique_samps <- unique(data_in$sample)

no_compounds <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  print(paste0("completed ",i))
}

mean(no_compounds) # 5292.368
sd(no_compounds) # 583.9855

phy <- readRDS("phy-phyloseq-object-Forslund-SWE-T2D.RDS")
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 19099 taxa and 145 samples ]
# sample_data() Sample Data:       [ 145 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 19099 taxa by 4 taxonomic ranks ]

phy <- prune_samples(unique_samps, phy)
min(taxa_sums(phy)) # o
phy

# prune taxa that have zero sequence reads
phy <- prune_taxa(taxa = taxa_sums(phy) > 0, x = phy)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17962 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 17962 taxa by 4 taxonomic ranks ]

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 76
NonZeroFxns

mean(NonZeroFxns) # 7861.039
sd(NonZeroFxns) # 1698.778
#rm(NonZeroFxns)

# check
unique_samps[1] # "ERR260139"
a <- prune_samples(unique_samps[1], phy)
min(taxa_sums(a)) # o
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a

unique_samps[76] # "ERR275252"
a <- prune_samples(unique_samps[76], phy)
min(taxa_sums(a)) # o
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a

# check ages of subjects?
df.samp <- readRDS("df.samp.with-t2dclass-age-Forslund-SWE-T2D.RDS")
head(df.samp)

sel <- which(row.names(df.samp) %in% unique_samps)
summary(df.samp$age[sel])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 68.96   69.98   70.34   70.44   71.10   71.84
table(df.samp$group_new[sel])
# T2D met neg T2D met pos         IGT      Normal 
# 33           0           0          43



#-------------------------


