library(tidyverse)
#install.packages('reticulate')
library(reticulate)
library(rcdk)
library(tidyverse)
#use_python("/usr/local/bin/python")
use_python("C:/Users/HelenSepman/AppData/Local/r-miniconda/envs/r-reticulate/python.exe")
#use_python("C:/Users/karpa/AppData/Local/Microsoft/WindowsApps/python.exe")


# reticulate::py_install("mordred")
# reticulate::py_install("rdkit")
# reticulate::py_install("pandas", pip = TRUE)


#---- PaDEL ----
PaDEL_original = function(standards) {
  SMILES_list = standards %>% 
    select(SMILES) %>% 
    unique() %>%
    na.omit() %>%
    mutate(Name = row_number())
  standards = standards %>%
    left_join(SMILES_list) 
  
  write_delim(SMILES_list %>% select(SMILES) %>% unique(),
              "SMILES.smi",
              col_names = FALSE)
  
  command = "java -jar PaDEL-Descriptor/PaDEL-Descriptor.jar -dir" #file name where they will be calculated
  command_final = paste(command, "SMILES.smi", "-file", "descs_calc.csv", "-2d", sep =" ") #makes text for command prompt
  javaOutput = system(command_final, intern = TRUE) #goes into commant prompt
  #PaDEL saves the descriptors to the local folder
  descs = read_delim("descs_calc.csv",
                     delim = ",",
                     col_names = TRUE)
  
  descs = descs %>%
    mutate_all(~replace(., str_detect(., "Infinity"), as.numeric(0))) %>%
    group_by(Name) %>%
    mutate(Name = str_split(Name, pattern = "_")[[1]][length(str_split(Name, pattern = "_")[[1]])]) %>%
    ungroup() %>%
    mutate(Name = as.integer(Name)) 
    
  cols <- names(descs)
  descs[cols] <- lapply(descs[cols], as.numeric)
  
  
  descs = descs %>%
    left_join(standards) %>%
    select(colnames(standards), everything()) %>%
    select(-Name)
  
  # write_delim(descs,
  #             "data/descs_calc.csv",
  #             delim = ",")
  
  return(descs)
}


#---- Mordred ----

source_python("C:/Users/HelenSepman/OneDrive - Kruvelab/computational/IE mudeli script ja failid/students/Helen/Mordred_python.py")
np = import("numpy", convert = FALSE)
pandas = import("pandas", convert = FALSE)
rdkit = import("rdkit", conver = FALSE)
mordred = import("mordred", convert = FALSE)

Mordred_calculation_R = function(SMILES) {
  descriptors = py_to_r(Mordred_calculation(SMILES))
  return(descriptors)
}

unlisting <- function(listed) {
  if (is.numeric(listed[[1]])){
    new <- listed[[1]]
  } else{
    new <- NA
  }
  return(new)
}

Mordred_descriptors = function(SMILES_list) {
  descriptors = tibble()
  for (SMILES in SMILES_list$SMILES) {
    descriptor = Mordred_calculation_R(SMILES)
    descriptor_this = tibble(value = descriptor$`_values`)
    if (dim(descriptors)[1] == 0) {
      descriptors = descriptor_this
    } else {
      descriptors = descriptors %>%
        bind_cols(descriptor_this)
    }
  }
  names_Mordred = read_delim("Mordred_names.txt",
                             delim = ",",
                             col_names = FALSE) %>%
    rename(desc_name = X1)
  descriptors = descriptors %>%
    bind_cols(names_Mordred) 
  descriptors = descriptors%>%
    gather(key = "compound", value = "value", -dim(descriptors)[2]) %>%
    spread(desc_name, value) %>%
    bind_cols(tibble("SMILES" = SMILES_list$SMILES)) %>%
    select(-compound) %>%
    select(SMILES, everything())
  
  descriptors = descriptors%>%
    group_by(SMILES) %>%
    mutate_if(is.list, funs(unlisting(.)))%>%
    ungroup()
  return(descriptors)
}

#SMILES_list = tibble("SMILES" = c('c1ccccc1N', 'c1ccccc1Cl', 'c1ccccc1Cl'))
#descriptors = Mordred_descriptors(SMILES_list)


#---- Morgan-2 ----
source_python("C:/Users/HelenSepman/OneDrive - Kruvelab/computational/IE mudeli script ja failid/students/Helen/Morgan2_python.py")

Morgan2_calculation_R = function(SMILES) {#, radius = 2, nBits = 1024) {
  descriptors = Morgan2_calculation(SMILES)#, radius, nBits))
  return(descriptors)
}

Morgan2_descriptors = function(SMILES_list) {
  descriptors = tibble()
  for (SMILES in SMILES_list$SMILES) {
    descriptor = Morgan2_calculation_R(SMILES)
    descriptor_this = tibble(descriptor)
    if (dim(descriptors)[1] == 0) {
      descriptors = descriptor_this
    } else {
      descriptors = descriptors %>%
        bind_cols(descriptor_this)
    }
  }
  names_Morgan2 = tibble(desc_name = 1:1024) %>%
    mutate(desc_name = str_c("Morgan2_", desc_name))
  descriptors = descriptors %>%
    bind_cols(names_Morgan2) 
  descriptors = descriptors%>%
    gather(key = "compound", value = "value", -dim(descriptors)[2]) %>%
    spread(desc_name, value) %>%
    bind_cols(tibble("SMILES" = SMILES_list$SMILES)) %>%
    select(-compound) %>%
    select(SMILES, everything())
  return(descriptors)
}


#---- Map4 ----

# with Linux-based Virtual Machine

#---- SIRIUS fingerprints ----


# ---- RCDK fingerprints ----

# Cannot calculate first 18 descriptors as SMILES do not have 3D info
rcdk_fingerprints <- function(SMILES_list) {
  dn <- get.desc.names(type = "all")
  mols <- parse.smiles(SMILES_list$SMILES)
  allDescs <- eval.desc(mols, dn)
  allDescs <- allDescs %>%
    rownames_to_column(var = "SMILES")
  return(allDescs)
}


#---- pretraining data cleaning ----
cleaning_descriptors = function(descs,
                         nearZeroVar_freqCut = 80/20,
                         highlyCorrelated_cutoff = 0.75) {
 
   # Removing columns with missing values
  descs = descs %>%
    select_if(~ sum(is.na(.))< 10,) %>%
    drop_na()
  
   # new_descs <- descs
   # new_descs[, sapply(new_descs, class) != "environment"]
  
  smiles = descs %>%
    select(SMILES)
  
  descs = descs %>%
    dplyr::select(-SMILES)
  
  # Checking that any of the categorical values would not have more than 80% of existing/non-existing values
  descs = descs %>%  
    select(-c(nearZeroVar(descs, 
                          freqCut = nearZeroVar_freqCut)))
  
  # Removing columns that are correlated to each other
  correlationMatrix <- cor(descs, use = "complete.obs")
  highlyCorrelated <- findCorrelation(correlationMatrix, 
                                      cutoff = highlyCorrelated_cutoff)
  
  descs <- descs %>%
    dplyr::select(-highlyCorrelated) %>%
    bind_cols(smiles)
  
  return(descs)
}


Fingerprint_calc <- function(compoundslist){
  library(rcdklibs)
  library(rcdk)
  
  compoundslist_onlySMILES <- compoundslist %>%
    select(SMILES) %>%
    unique()
  
  # substructure fingerprints
  row <- c(1:307)
  list <- as.data.frame(row)
  
  col_names_substr <- paste("Un", c(55:361), sep = "")
  
  substr_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  error_comp <- tibble()
  
  for (n in 1:length(compoundslist_onlySMILES$SMILES)){
    #col <- compoundslist$id[[n]]
    SMILES <- compoundslist_onlySMILES$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      substr_fingerprints <- get.fingerprint(mol2,
                                             type = "substructure")
      
      table <- as.data.frame(substr_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      suppressMessages(datarow <- list %>%
                         left_join(table) %>%
                         select(-row))
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      substr_table <- substr_table %>%
        add_row(newrow)
    }
    
    else{
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>%
        rbind(wrong_smile)
      print(SMILES)
    }
  }
  
  substr_table[is.na(substr_table)] <- 0
  
  substr_table <- substr_table %>%
    filter(rowname != "row")
  
  names(substr_table) <- c("SMILES", col_names_substr)
  
  #MACCS fingerprints
  row <- c(1:166)
  
  list <- as.data.frame(row)
  
  col_names_maccs <- paste("Un", c(362:527), sep = "")
  
  maccs_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist_onlySMILES$SMILES)){
    #col <- compoundslist$id[[n]]
    SMILES <- compoundslist_onlySMILES$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      maccs_fingerprints <- get.fingerprint(mol2,
                                            type = "maccs")
      
      table <- as.data.frame(maccs_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      suppressMessages(datarow <- list %>%
                         left_join(table) %>%
                         select(-row))
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      maccs_table <- maccs_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  maccs_table[is.na(maccs_table)] <- 0
  
  maccs_table <- maccs_table %>%
    filter(rowname != "row")
  
  names(maccs_table) <- c("SMILES", col_names_maccs)
  
  
  #pubchem fingerprints
  row <- c(1:881) 
  
  list <- as.data.frame(row)
  
  col_names_pubchem <- paste("Un", c(528:1408), sep = "")
  
  pubchem_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist_onlySMILES$SMILES)){
    #col <- compoundslist$id[[n]]
    SMILES <- compoundslist_onlySMILES$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      pubchem_fingerprints <- get.fingerprint(mol2,
                                              type = "pubchem")
      
      table <- as.data.frame(pubchem_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      suppressMessages(datarow <- list %>%
                         left_join(table) %>%
                         select(-row))
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      pubchem_table <- pubchem_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  pubchem_table[is.na(pubchem_table)] <- 0
  
  pubchem_table <- pubchem_table %>%
    filter(rowname != "row")
  
  names(pubchem_table) <- c("SMILES", col_names_pubchem)
  
  #KlekRoth fingerprints
  
  row <- c(as.numeric(1:4860))
  
  list <- as.data.frame(row)
  
  col_names_KlekotaRoth <- paste("Un", c(1409:6268), sep = "")
  
  
  KlekotaRoth_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist_onlySMILES$SMILES)){
    #col <- compoundslist$id[[n]]
    SMILES <- compoundslist_onlySMILES$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      KlekotaRoth_fingerprints <- get.fingerprint(mol2,
                                                  type = "kr")
      
      table <- as.data.frame(KlekotaRoth_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      suppressMessages(datarow <- list %>%
                         left_join(table) %>%
                         select(-row))
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      KlekotaRoth_table <- KlekotaRoth_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  KlekotaRoth_table[is.na(KlekotaRoth_table)] <- 0
  
  KlekotaRoth_table <- KlekotaRoth_table %>%
    filter(rowname != "row")
  
  names(KlekotaRoth_table) <- c("SMILES", col_names_KlekotaRoth)
  
  
  #custommadeSMARTS
  row <- c(1:202)
  list <- as.data.frame(row)
  
  smartslist <- c("C1(CCC3C2CCCC3)CCCCC12","CC(C)CC(C(~O)[O,N])N","CC(C(C(~O)[O,N])N)O",
                  "[OH0]=[CH0](O)[CH1][CH2][OH1]","[OH0]=[CH0](N)[CH1](N)[CH2][OH1]","[OH0]=[CH0](N)[CH1]N",
                  "[CH2]([CH1]([CH2][OH0][PH0](~O)(~O)~O)~O)~O","[CH2]([CH1]([CH1]([CH1]([CH1]1~[#7])~O)~O)[OH0]1)[OH0][PH0](~O)(~O)~O",
                  "[#6]~C([CH0]([NH1][CH1](~[#6])[CH0](~[!#1])~O)~[OH0D1])~[!#1]",
                  "[CH1](~[!#1])[CH1]([CH1]([CH1]([CH1]([CH2]~O)~O)~O)~O)[NH1][CH0](~[CH3D1])~[OH0D1]","c(:c:c:n1[CH2][CH2]~C):c1",
                  "[CH2](c(:c:n:c1):n1)[CH1]([CH0](~[!#1])~O)~N","[CH2]([CH2]~[!#1])[NH0](~[#6])~[#6]","[#6]~C([NH1][CH1]([CH2]c(:c:n:c1:c:c:c:c2):c12)[CH0](~[!#1])~O)~[OH0D1]",
                  "[#6]~S[CH2][CH2][CH1]([CH0](~[!#1])~O)~N","[#6]~C([NH1][CH1]([CH2][CH2][CH0](~[!#1])~O)[CH0](~[!#1])~O)~[OH0D1]",
                  "[#6]~Cc(:c:c:c(:c1~[!#1])~[!#1]):c1","[#6]~C([CH0](~[!#1])~[!#1])[NH0](~[#6])[CH0](~[!#1])~O","c(:c(~[!#1]):n:c(:n1)~[!#1])(:c1~[!#1])~[!#1]",
                  "[#6]~C[CH0](~[!#1])[OH0][CH1](~[#6])~[#6]","[#6]~C[CH2][CH2][CH2][CH0]([NH1]~[#6])~[OH0D1]",
                  "[#6]~C[OH0][CH1]([CH1]([CH1]([CH1]([CH2]1)~O)~O)~O)[OH0]1","[CH1]([CH0](~C)(~[#6])~[!#1])(~C)[OH0][CH1](~C)~[!#1]",
                  "[#6]~C(~[!#1])[CH2][CH1]([OH0]~[#6])[OH0]~[#6]","[#6]~C([CH0](~[!#1])~O)[NH1][CH0]([CH2]~[!#1])~[OH0D1]","[#6]~n(:c:c:n1):c1",
                  "[#6]~Oc(:c:c:c:c1~[!#1]):c1~[!#1]","[#6]~C(~C)[OH0][PH0](~O)(~O)~O","c(:c:c:c(:c1)[NH0](~[!#1])~[!#1]):c1",
                  "[#6]~C([CH1](~[#6])~[!#1])[NH1][CH0](~[#6])~[OH0D1]","[CH2](~[#6])[NH1][CH0](~[!#1])~[!#1]","[#6]~O[CH1]([CH0](~C)(~[#6])~[!#1])~[#6]",
                  "[#6]~O[CH1]([CH2]~[!#1])[CH2]~[!#1]","[#6]~C(~[!#1])[NH0](~[#6])~[#6]","[CH2]([CH1]([CH1]([CH1]([CH1]1~[!#1])~O)~O)[OH0]1)~O",
                  "[#6]~Cc(:c:c:c(:c1)~[!#1]):c1~[!#1]","[#6]~Cc(:c:c:c:c1~[!#1]):c1~[!#1]","[#6]~C([CH1]([CH1]([CH1]([CH1]1[OH0]c(:c:c:c:c2):c2)~O)~O)~O)[OH0]1",
                  "[CH2]([CH0](~[!#1])~[!#1])[OH0]~[#6]","[#6]~C(~[#6])[CH0](~[!#1])[OH0]~[#6]","c(:c:c:c(:c1)[CH0](~[!#1])(~[!#1])~[!#1]):c1",
                  "[#6]~C[OH0]c(:c:c:c(:c1)~[!#1]):c1","[#6]~C=[CH1][CH2][CH2][CH2][CH0](~[!#1])~O","[#6]~Oc(:c:c(:c:c1)[CH2]~[!#1]):c1~[!#1]",
                  "[#6]~N[CH1]([CH2]~[#6])[CH0](~[!#1])~O","[#6]~C[CH2][CH0](~[!#1])[OH0][CH1](~[#6])~[#6]","[#6]~C[CH2][OH0][CH0](~[!#1])~[#6]",
                  "[#6]~C([CH2][OH0][PH0](~O)(~O)~O)~[!#1]","[#6]~C([CH0](~[!#1])~O)[NH1][CH0]([CH1]([CH1](~[#6])~[!#1])~[!#1])~[OH0D1]",
                  "c(:c:c(:c([CH0]~[!#1]):c1~[!#1])~[!#1]):c1","[#6]~Cc(:c:c:c(:c1[CH2][CH1]=[CH0](~[#6])~[#6])~[#8]):c1","c(:c:c(:c:c1)~[!#1]):c1[CH1]~[!#1]",
                  "[#6]~C[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH0](~[OH0D1])[OH0][CH1](~[#6])~[#6]","[#6]~O[CH0](~[!#1])[CH2]~[!#1]",
                  "[#6]~O[CH1]([CH1]([CH1]([CH1]([CH1]1[CH2]~O)~[OH1D1])~[OH1D1])~O)[OH0]1","[#6]~N[CH2][CH2][CH2][CH1]~[!#1]",
                  "c(:c:c(:c:c1)~[!#1]):c1[CH2]~[!#1]","[#6]~C([CH1]([CH0](~C)~[!#1])~[!#1])[OH0]~[#6]","[CH2]([CH2][NH1][CH0](~[!#1])~[OH0D1])[CH0](~[!#1])~O",
                  "[#6]~C([CH2][CH2][CH1]([CH0](~[#6])(~[#6])[CH1]([CH2][CH2]1)~O)[CH0]12~[#6])[CH1]2[CH2]~[#6]","[#6]~Cc(:c:c(:c:c1)~[!#1]):c1~[!#1]",
                  "[CH2]([CH0]([CH2][OH0][CH1]1~[OH0D2])([CH1]1~[OH1D1])~[OH1D1])~[OH1D1]","[#6]~Oc(:c:c:c(:c1)~[!#1]):c1~[!#1]",
                  "[#6]~N[CH2][CH1]([CH0](~[!#1])~[!#1])~[#6]","[CH2]([CH2][CH1]([CH0]1(~[CH3D1])~C)~C)[CH1]1~[!#1]","[#6]~O[CH0](~[!#1])[CH1](~[#6])~[!#1]",
                  "[#6]~C[OH0][PH0](~O)(~[!#1])[OH0][CH2]~[#6]","[#6]~c(:c(:c:c(:c1)[OH0]~[#6])~[!#1]):c1~[!#1]","[#6]~C([CH0](~[!#1])~[!#1])[OH0][CH0](~[!#1])~[#6]",
                  "[#6]~C[CH2][NH0](~[#6])~[!#1]","[CH1](~[!#1])[CH1]([CH1]([CH1]([CH1]([CH2]~[!#1])~O)~O)~O)~[!#1]","[#6]~C[CH1]([CH2]~[#6])[OH0]~[#6]",
                  "[#6]~C([CH2]~[#6])[NH0](~[#6])~[#6]","[#6]~C[CH1](~[#6])[CH0](~[#6])[CH0](~[!#1])[OH0]~[#6]","[#6]~C[OH0][CH1](~[#6])[OH0][CH2]~[#6]",
                  "[#6]~C[CH2][CH0](~O)[OH0][CH2][CH1]~[#6]","[#6]~c(:c:c(:c(:c1)~[!#1])~[!#1]):c1~[!#1]","[#6]~c(:c:c:c:c1):c1[OH0]~[#6]",
                  "[#6]~c(:c:c(:c:c1~[!#1])~[!#1]):c1","[#6]~C[CH1]=[CH1][CH2][CH2][CH0](~[!#1])~O","[#6]~C(~[#6])(~[#6])[OH0]~[#6]",
                  "[#6]~N[CH0](~[!#1])[NH1]~[#6]","c(:c:c(:c:c1)~[!#1]):c1[CH0](~[!#1])~[!#1]",
                  "[#6]~C([CH1]([CH1]([CH1]([CH1]1~[OH0D2])~O)[OH0]~[#6])~O)[OH0]1","[#6]~C[OH0][CH1](~[#6])[OH0][CH1](~[#6])[CH0](~[!#1])~[!#1]",
                  "[#6]~c(:c:c:c(:c1[OH0]~[#6])~[!#1]):c1","[#6]~c(:c:c:c:c1):c1[OH0]c(:c:c:c:c2):c2","[#6]~c(:c:c:c:c1[CH2]~[!#1]):c1",
                  "c(:c(:c:c(:c1[CH1]~[!#1])~[!#1])~[!#1]):c1~[!#1]","[#6]~C([CH0]([NH1][CH1]([CH2]c(:c:c:c:c1):c1)[CH0](~[!#1])~O)~[OH0D1])~[!#1]",
                  "c(:c:c:c(:c1)[CH1](~[#6])~[!#1]):c1","c(:c:c(:c:c1~[!#1])[CH0]~[!#1]):c1","[#6]~C(~[!#1])[CH1](~[#6])[OH0][CH0](~[!#1])~[#6]",
                  "[#6]~Sc(:c:c:c:c1):c1","[CH2]([CH2][OH0][PH0](~O)~O)~[!#1]","[#6]~C[OH0][CH1]([CH0](~C)(~[#6])~[!#1])~[!#1]",
                  "[#6]~C(c(:c:c:c:c1~[!#1]):c1)~[!#1]","c(:c:c(:c:c1~[!#1])[CH2]~[!#1]):c1","[#6]~C[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH0](~O)[OH0]~[#6]",
                  "[#6]~C[CH2][CH1]([CH0](~[!#1])~O)[NH1][CH0](~[#6])~[OH0D1]","[#6]~c(:c:c:c(:c1~[!#1])~[!#1]):c1~[!#1]","[#6]~C[OH0][CH0](~[!#1])[CH1]~[!#1]",
                  "c(:c:c(:c(:c1)[NH0](~[!#1])~[!#1])~[!#1]):c1~[!#1]","[CH2]([CH1]([CH0](~[!#1])~O)[NH1][CH0]([CH1]([CH2]~[#6])~[#7])~[OH0D1])~[!#1]",
                  "[#6]~c(:c:c:c:n1~[!#1]):c1","c(:c([CH2]~[!#1]):c(:c:c1~[!#1])~[!#1]):c1~[!#1]","[#6]~Nc(:c:c:c:c1):c1",
                  "[#6]~C[CH0](~[#6])([CH0](~[!#1])[OH0]~[#6])~[!#1]","[#6]~O[CH0](~[!#1])[CH2][CH2]~[#6]",
                  "c(:c:c:c(:c1[CH1]~[!#1])~[!#1]):c1","[#6]~C(=[CH1]~[!#1])[CH0](~[!#1])[OH0]~[#6]",
                  "[#6]~c(:c(:c:c:c1~[!#1])~[!#1]):c1[OH0]~[#6]","[#6]~C([CH1]([CH1]([CH1]([CH1]~[!#1])~[!#1])~[!#1])[OH0]~[#6])~[!#1]",
                  "[CH2]([CH1]([CH0]([NH1]~[#6])~[OH0D1])~[!#1])~[#6]","[#6]~C[NH1][CH2][CH2]~[!#1]","[#6]~C[CH1]([CH1](~[#6])[CH0](~[!#1])[OH0]~[#6])~[!#1]",
                  "[#6]~C([NH1][CH1]([CH2]~[!#1])[CH0](~[!#1])~O)~O","[#6]~c(:c:c(:c:c1)[OH0]~[#6]):c1~[!#1]",
                  "[#6]~C(c(:c:c:c(:c1)~[!#1]):c1)~[!#1]","[#6]~C([CH0](~[!#1])~O)[NH1][CH0]([CH1]([CH2]c(:c:c:c:c1):c1)~[#7])~[OH0D1]",
                  "[#6]~c(:c:c:c(:c1)[CH2]~[!#1]):c1~[!#1]","[#6]~c(:c:c:c:c1):c1[CH2]~[!#1]","c(:c:c:c(:c1)[OH0][CH1]~[!#1]):c1",
                  "c(:c(:c:c(:c1~[!#1])~[!#1])[CH0]~[!#1]):c1~[!#1]","[#6]~C(~[#6])[CH0](~[!#1])[OH0][CH1](~[#6])~[#6]","[#6]~C([CH1]([CH1]([CH1]([CH1]1~O)[OH0]~[#6])~O)~O)[OH0]1",
                  "[#6]~Oc(:c:c(:c:c1~[!#1])~[!#1]):c1","[#6]~C[OH0][CH0](~[!#1])[CH2][CH2]~[#6]","[#6]~C[OH0][CH1](~[#6])[OH0][CH0](~[!#1])~[#6]",
                  "[#6]~C([CH0]([NH1][CH1]([CH2][CH0](~[!#1])~O)[CH0](~[!#1])~O)~O)~N","[#6]~C(~[#6])[NH0](~C)~[!#1]","[CH1][NH1][CH0](=[OH0])[CH1]([CH2]c)[NH2]",
                  "[CH3][CH2][CH1]=[CH1][CH2][CH1]=[CH1][CH2][CH1]=[CH1][CH2][CH1]","c[OH0][CH1]([CH1])c","[CH1][OH0]c",
                  "[CH3][CH1]([CH1])[OH1]","[CH3][CH2][CH2][CH2][CH2][CH1]","[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH1]","[CH1][CH2]c",
                  "[CH1][CH1]=[CH1][CH1]","[CH0]([CH3])[CH2][CH2][CH0]","[CH1][CH0](=[OH0])[NH1][CH1]([CH2]c)[CH0](=[OH0])[OH1]","[CH1]n",
                  "[CH3][OH0][CH1]","[CH1][CH2][CH1]([CH1])[OH1]","[CH1][CH1]([CH1]([CH2][OH1])[OH1])[OH1]","c[CH3]",
                  "[CH1][CH2][CH2]c","[CH1][CH2][CH2][CH2][CH0](=[OH0])[OH1]","[CH1][CH2][CH0](=[OH0])[NH2]",
                  "[CH1][CH2][CH0][CH3]","[CH0][CH0](=[OH0])[OH1]","[CH0]([CH1])[CH0](=[OH0])[OH0][CH3]","[CH1][CH2][CH2][CH0][OH1]",
                  "[CH0]([CH2][OH1])[OH0][CH1]","[CH1][CH2][CH1]","[CH1][CH2][CH2][CH1]","[CH1]c","[CH1][CH2][OH0][CH1]","[CH1][CH1]([CH2][CH1])[CH1]",
                  "[CH2]([CH2][CH0](=[OH0])[OH1])[CH1]","c[CH0][OH1]","[CH3][OH0][CH0](=[OH0])c","[CH1][CH0](=[OH0])[OH0][CH1]","[CH0][CH2][OH1]",
                  "[NH0][CH0](=[OH0])[CH1]","[CH1][CH2][CH0](=[OH0])[OH1]","[CH0]([CH1])c","[CH1][CH2][CH2][CH2][CH2][NH2]","[CH0][CH0](=[OH0])[OH0][CH1]",
                  "c[CH0]","[CH1][CH2][CH1]=[CH0]","[CH1][CH1]([CH3])[OH0][CH1]","c[CH0](=[OH0])[OH0][CH1]","c[CH0](=[OH0])c","[CH1][CH2][CH2][CH0]",
                  "[CH0][CH1]","[CH1][NH0]","[CH1][CH2][CH0](=[OH0])c","[CH0][CH0]","[CH1][CH2][CH0]","c[CH0](=[OH0])[CH1]","c[CH1]=[CH1][CH0]",
                  "[NH0]c","[CH1][CH2][OH0][PH0](=[OH0])([OH1])[OH1]","c[OH0][CH0]","[CH0][CH0][CH3]","[CH3][CH1]([CH3])[CH1]","c[CH1]=[OH0]",
                  "[CH3][CH2][CH2][CH2][CH2][CH1]=[CH1][CH2][CH1]","c[OH0][CH2][CH0]","[CH0]([CH1])[OH1]","c[CH1]=[CH1][CH0](=[OH0])[OH0][CH1]",
                  "[CH3][CH0](=[OH0])[NH1][CH1]","[CH0](c)c","[CH0][NH2]","[CH0]([CH0])[OH1]","[CH0][NH1]c","[CH0]([CH3])[OH0][CH1]",
                  "[CH3][OH0][CH0][CH0]","c[CH0](=[OH0])[OH1]","[CH1][OH0][CH0](=[OH0])[CH3]")
  
  col_names_custom <- c("Un8179","Un8181","Un8182","Un8183","Un8184","Un8185","Un8187","Un8190","Un8191","Un8192","Un8193","Un8194",
                        "Un8197","Un8198","Un8199","Un8200","Un8202","Un8203","Un8205","Un8206","Un8207","Un8210","Un8213","Un8215",
                        "Un8217","Un8218","Un8219","Un8220","Un8221","Un8222","Un8224","Un8225","Un8226","Un8228","Un8229","Un8230",
                        "Un8232","Un8233","Un8234","Un8235","Un8236","Un8237","Un8239","Un8240","Un8241","Un8242","Un8243","Un8244",
                        "Un8245","Un8246","Un8247","Un8248","Un8250","Un8251","Un8252","Un8253","Un8254","Un8255","Un8256","Un8258",
                        "Un8259","Un8260","Un8261","Un8262","Un8263","Un8264","Un8266","Un8267","Un8268","Un8269","Un8270","Un8271",
                        "Un8272","Un8273","Un8274","Un8275","Un8276","Un8277","Un8280","Un8282","Un8283","Un8284","Un8285","Un8286",
                        "Un8287","Un8290","Un8292","Un8294","Un8295","Un8296","Un8297","Un8299","Un8303","Un8304","Un8308","Un8309",
                        "Un8310","Un8311","Un8313","Un8314","Un8315","Un8316","Un8318","Un8319","Un8320","Un8321","Un8323","Un8324",
                        "Un8325","Un8327","Un8329","Un8332","Un8333","Un8336","Un8337","Un8338","Un8339","Un8340","Un8341","Un8342",
                        "Un8345","Un8347","Un8348","Un8349","Un8351","Un8353","Un8354","Un8355","Un8356","Un8359","Un8360","Un8363",
                        "Un8367","Un8374","Un8376","Un8378","Un8379","Un8380","Un8381","Un8382","Un8383","Un8384","Un8385","Un8386",
                        "Un8387","Un8388","Un8390","Un8394","Un8395","Un8396","Un8397","Un8399","Un8400","Un8401","Un8402","Un8403",
                        "Un8405","Un8406","Un8407","Un8408","Un8409","Un8410","Un8411","Un8412","Un8413","Un8414","Un8415","Un8416",
                        "Un8417","Un8420","Un8421","Un8422","Un8423","Un8424","Un8426","Un8427","Un8428","Un8431","Un8432","Un8433",
                        "Un8434","Un8436","Un8437","Un8438","Un8439","Un8440","Un8441","Un8443","Un8444","Un8445","Un8446","Un8447",
                        "Un8448","Un8450","Un8451","Un8453","Un8454","Un8455","Un8456","Un8457","Un8460","Un8461")
  
  
  custom_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist_onlySMILES$SMILES)){
    #col <- compoundslist$id[[n]]
    SMILES <- compoundslist_onlySMILES$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      custom_fingerprints <- get.fingerprint(mol2,
                                             type = "substructure",
                                             substructure.pattern = smartslist)
      
      table <- as.data.frame(custom_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      suppressMessages(datarow <- list %>%
                         left_join(table) %>%
                         select(-row))
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      custom_table <- custom_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  custom_table[is.na(custom_table)] <- 0
  
  custom_table <- custom_table %>%
    filter(rowname != "row")
  
  names(custom_table) <- c("SMILES", col_names_custom)
  
  
  #ringsystems
  # list <- read_delim("ringsystem_fingerprints_names.txt", delim = ",") %>%
  #   mutate(row = as.numeric(rowRing)) %>%
  #   select(row)
  
  row <- c(1:113)
  list <- as.data.frame(row)
  
  smartslistring <- c("[CH1]([CH1][OH0][CH1]-1~[#7])([CH1]-1~[#8])~[#8]","c(:c(:c(:c:c:1):o:c:c:2):c:2):c:1","C(C(C(C(C-1~[#8])~[#8])~[#8])~[#8])(C-1~[#8])~[#8]",
                      "c(:c:o:c(:c:c(:c:c:1~[#8])~[#8]):c:1:2):c:2~[#8]","c(:c(:c:c:c:1~[#17])~[#8]):c:1~[#8]","c(:c(:n:c(:n:1)~[#8])~[#7])(:c:1~[#8])~[#7]",
                      "c(:c(:c:c:c:1)~[#8]):c:1","C(CC(OC=1)~[#8])C=1","c(:c:c:c:c:1:c:c:n:2):c:1:2","c(:c(:c:c:c:1)~[#16]):c:1",
                      "[CH1]([CH1][CH1]([CH2]-1)~[#8])[CH1]-1~[#8]","c(:c:o:c:c:1):c:1~[#8]","c(:c(:n:c(:n:1)~[#7])~[#7]):c:1~[#8]",
                      "C(CCC(C-1CCC-2)C-2)(C-1CC=C-3C(C(CC-4)CCC-5)C-5)C-3-4",  "c(:c:c(:c:n:n:1):c:1:c:2):c:2","c(:c(:c:c:c:1CCCO-2)~[#8]):c:1-2",
                      "C(CCOC-1~[#8])C-1","c(:c:c:c:c:1):c:1~[#7]","c(:c:c:c(:c:1C(C-2)~[#8])N-2):c:1","c(:c(:c(:c:c:1)~[#8])~[#8])(:c:1)~[#8]",
                      "C(CCC-1)(C-1)~[#8]","c(:c(:c:c(:c:1)~[#8])~[#8])(:c:1~[#8])~[#8]","c(:n:c:c:1):n:1","c(:n:c:n:1):n:1",
                      "C(C(C(OC-1)~[#8])~[#7])(C-1~[#8])~[#8]","[CH1]([CH2][CH2][CH0][CH1]-1)[CH1]-1","c(:c(:c:c:c:1:c:c:c(:o:2)~[#8])~[#8]):c:1:2",
                      "c(:n:c:n:c:1~[#7])(:c:1)~[#7]","n(:c:c:c:c:1):c:1","C(COC-1)C-1","C(CCCC-1CC(CC-2)~[#8])C-1-2","C(C(COC-1)~[#8])(C-1~[#8])~[#8]",
                      "c(:c:o:c(:c:c(:c:c:1~[#8])~[#8]):c:1:2)(:c:2)~[#8]","[CH1]([CH0][CH2][CH2][CH2]-1)[CH0]-1","c(:c(:c:c(:c:1[CH0]([CH1][CH1][OH0]-2)~[#8])~[#8])~[#8]):c:1-2",
                      "C(CCC(C-1CCC-2CCC-3)C-2-3)C-1","c(:c:c:n:c:1~[#8]):n:1","n(:c:c:c(:c:1:c:c:c:c:2)~[#8]):c:1:2",   "[CH1]([CH1][CH1]([CH2]-1)~[#8])[CH0]-1~[#8]",
                      "c(:c(:c:c(:c:1)~[#8])~[#8]):c:1~[#8]","c(:c:c:c(CCCC-1):c-1:2):c:2","C(CCCC-1~[#8])(C-1)~[#8]",
                      "[CH1]([CH2][CH1][OH0][CH1]-1)[CH1]-1","c(:c:c(:c(:c:c:c:o:1):c:1:2)~[#8])(:c:2)~[#8]",
                      "[CH0]([CH1]([CH1]([CH1]=[CH1][OH0]-1)[CH1]-2~[#8])[CH1]-1~[#8])([CH1]-2-3)[OH0]-3",
                      "C(CCC(C-1CC(C-2~[#8])~[#8])C-2)(C-1CC=C-3C(C(CC-4)CCC-5)C-5)C-3-4","[CH0]([CH2][CH2][CH1][CH0]-1)[CH1]-1",
                      "C(CC(C-1)~[#8])C-1~[#8]","c(:c:c:c(:c:1)~[#9])(:c:1)~[#7]","C(COC-1)(C-1)~[#8]","c(:n:c:n:c:1~[#7])(:c:1)~[#8]",    
                      "c(:c:c:n:c:1~[#7]):n:1","c(:c(:c(:c(:c:1)~[#8])[CH0]([CH2][CH1]-2)~[#8])[OH0]-2):c:1~[#8]","C(CCCC-1~[#8])C-1",
                      "n(:c:c:c:1):c:1","c(:c:c(:n:c:1)~[#8]):n:1","c(:c:o:c(:c:c:c:c:1):c:1:2)(:c:2~[#8])~[#8]",
                      "[CH1]([CH2][CH2][CH0][CH1]-1)[CH2]-1","[CH0]([CH1][CH2][CH2]-1)[CH0]-1~[#8]","c(:n:c:n:c:1)(:c:1)~[#7]",
                      "C(C(COc:1:c:c(:c:c:2)~[#8])~[#8])c:1:2","c(:c:c(:c(:c:1:c:c:c:2~[#8]):c:2)~[#8]):o:1","c(:c(:c:c:c:1~[#8])~[#8])(:c:1)~[#8]",
                      "C(C=CCC-1CCCC-2)C-1-2","C(CCC(C-1CC(C-2CCC-3)~[#8])C-2-3)C-1","c(:c:c(:c(:c:1:c:c:c:2~[#8]):c:2~[#8])~[#8]):o:1",
                      "C(CCCC-1CCCC-2)C-1-2","[CH0]([CH2][CH2][CH1]([CH0]-1)~[#8])([CH1]-1[CH2][CH2][CH0]-2)[CH1]-2","C(CCC-1CCCC-2)C-1-2",
                      "C(CC(CC-1~[#8])~[#8])(C-1~[#8])~[#8]","C(CC(C-1)~[#8])(O-1)~[#7]","C(CC(CC-1~[#8])~[#8])O-1","c(:c:c:c:c:1OCCC-2~[#8])(:c:1-2)~[#8]",
                      "C(CCC-1C(C(CC-2)C(C(C-3)CC(C-4)~[#8])C-4)C-3)C-1-2",
                      "c(:n:c(:n:c:1)~[#8])(:c:1)~[#7]","C(C(C(C-1)~[#8])~[#8])O-1",
                      "[CH0]([CH1]([CH2][CH2][CH0]([CH1]-1[CH2][CH1]=[CH0]-2[CH1]([CH0]([CH2][CH2]-3)[CH2][CH2][CH0]-4)[CH2]-4)[CH0]-2-3)[CH0]-1[CH2][CH2]-5)[CH1]-5~[#8]",
                      "[CH1]([CH0]-1)[OH0]-1","C(=COCC-1CCC-2)C-1-2","n(:c:n:c:1:c(:n:c:n:2)~[#7]):c:1:2","C(CCc(:c:c:c:c:1):c:1-2)O-2",
                      "c(:c:o:c(:c:1:c:c:c:2~[#8]):c:2)(:c:1~[#8])~[#8]","[CH0](=[CH1][OH0][CH1][CH1]-1)[CH1]-1","c(:c:c(:n:c:1~[#7])~[#7]):n:1",
                      "c(:c:c(:c:c:1)~[#8]):n:1","C(C=CC-1~[#8])O-1","c(:n:n:c:1):c:1","c(:c(:c(:c(:c(:c:c(:c:1)~[#8])~[#8]):c:1:2)~[#8])~[#8]):o:2",
                      "[CH1]([CH1][CH1][CH1]-1)[OH0]-1","[CH0]([CH1][CH0][CH2][CH2]-1)[CH0]-1","c(:c:c:c:c:1:n:c:n:2):c:1:2","n(:c:n:c:1:c(:n:c:n:2)~[#8]):c:1:2",
                      "c(:n:c:n:c:1:n:c:n:2):c:1:2","c(:n:c:c:c:1~[#8])(:n:1)~[#7]",  "[CH0]([CH2][CH2][CH0][CH1]-1)[CH2]-1","[CH0](=[CH1][CH0]([OH0]-1)~[#8])[CH2]-1","c(:c:n:c(:n:1)~[#8]):c:1~[#8]",
                      "C(COc(:c:c:c:c:1):c:1-2)(C-2~[#8])~[#8]","c(:c(:c:c(:c:c:c(:c:1)~[#8]):c:1:2)~[#8]):o:2","[CH2]([CH2][CH1][CH0][CH1]-1)[CH1]-1",
                      "c(:c(:c:c(:c:1:o:c:c:2):c:2~[#8])~[#8])(:c:1)~[#8]","C(COc(:c:c:c:c:1):c:1-2)C-2~[#8]","C(CC=CO-1)C-1",
                      "[CH0]([CH2][CH2][CH1]([CH0]-1)~[#8])[CH1]-1","C(COC-1~[#8])C-1","C(COC-1~[#8])(C-1~[#8])~[#8]","n(:c:c:n:1):c:1~[#16]",
                      "[CH1]([CH1][CH1][CH1][CH2]-1)[OH0]-1","c(:n:c:n:c:1:c:c:c:c:2)(:c:1:2)~[#8]","[CH1]([CH1][CH1]([CH2][CH1]-1~[#8])~[#8])[OH0]-1",
                      "c1ccccc1-c","[cD2]1ccc[cD2]c1-c(c)c",
                      "[C](-,=[C]-,=1)-,=[C](-,=[C]-,=2)-[C](-,=[C]-,=[C]-,=[C]-,=2-,=[O])-,=[C](-,=[C]-,=3)-,=[C]-,=1-,=[C](-,=[C]-,=4)-,=[C](-,=[C]-,=3)-,=[C]-,=[C]-,=4")
  
  
  col_names_ring <- c("Un8468","Un8480","Un8483","Un8484","Un8486","Un8489","Un8491","Un8494","Un8495","Un8496","Un8497","Un8498",
                      "Un8501","Un8502","Un8505","Un8516","Un8521","Un8524","Un8529","Un8539","Un8540","Un8541","Un8542","Un8546",
                      "Un8547","Un8553","Un8556","Un8558","Un8559","Un8563","Un8566","Un8568","Un8575","Un8581","Un8582","Un8584",
                      "Un8585","Un8587","Un8588","Un8589","Un8606","Un8607","Un8610","Un8611","Un8618","Un8619","Un8620","Un8632",
                      "Un8635","Un8637","Un8638","Un8641","Un8647","Un8648","Un8651","Un8655","Un8663","Un8664","Un8670","Un8673",
                      "Un8677","Un8678","Un8682","Un8686","Un8691","Un8692","Un8698","Un8706","Un8710","Un8723","Un8724","Un8730",
                      "Un8736","Un8740","Un8743","Un8749","Un8750","Un8751","Un8758","Un8760","Un8767","Un8770","Un8772","Un8777",
                      "Un8779","Un8782","Un8783","Un8788","Un8789","Un8797","Un8802","Un8808","Un8814","Un8820","Un8834","Un8843",
                      "Un8849","Un8852","Un8854","Un8858","Un8860","Un8864","Un8871","Un8874","Un8881","Un8888","Un8889","Un8901",
                      "Un8903","Un8908","Un8922","Un8923","Un8924")
  
  ring_table <- as.data.frame(t(list)) %>%
    rownames_to_column()
  
  for (n in 1:length(compoundslist_onlySMILES$SMILES)){
    #col <- compoundslist$id[[n]]
    SMILES <- compoundslist_onlySMILES$SMILES[[n]]
    
    mol2 <- parse.smiles(SMILES)[[1]]
    
    if(!is.null(mol2)){
      
      ring_fingerprints <- get.fingerprint(mol2,
                                           type = "substructure",
                                           substructure.pattern = smartslistring)
      
      table <- as.data.frame(ring_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)
      
      table <- table %>%
        mutate(row = as.numeric(row))
      
      suppressMessages(datarow <- list %>%
                         left_join(table) %>%
                         select(-row))
      
      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()
      
      ring_table <- ring_table %>%
        add_row(newrow)
    }
    
    else{
      print(SMILES)
    }
  }
  
  ring_table[is.na(ring_table)] <- 0
  
  ring_table <- ring_table %>%
    filter(rowname != "row")
  
  names(ring_table) <- c("SMILES", col_names_ring)
  
  
  
  #Combining all FP together
  suppressMessages(allFP_table <- substr_table %>%
                     left_join(maccs_table) %>%
                     left_join(pubchem_table) %>%
                     left_join(KlekotaRoth_table) %>%
                     left_join(custom_table) %>%
                     left_join(ring_table))
  
  suppressMessages(final_table <- compoundslist %>%
                     left_join(allFP_table))
  
  return(final_table)
}

