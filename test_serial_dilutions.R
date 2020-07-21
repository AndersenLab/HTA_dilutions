df <- data.frame(
    drug = "pzq",
    drugstock = 330,
    diluent = "DMSO",
    wellconc = "0,250,500,1000,1500,2000,3000",
    version = "v2",
    drugplates = 9,
    controlplates = 0,
    wells = 12,
    intconc = 330
)

df <- data.frame(
    drug = "emodepside",
    drugstock = 10,
    diluent = "DMSO",
    wellconc = "0,0.0098,0.0196,0.0391,0.0781,0.1563,0.3125",
    version = "v2",
    drugplates = 6,
    controlplates = 0,
    wells = 12,
    intconc = 10
)

      
# change int conc
dose <- df %>%
    tidyr::separate_rows(wellconc, sep = ",") %>%
    dplyr::mutate(version = as.character(version),
                  workingconc = as.numeric(ifelse(version == "v2", wellconc, wellconc*3)),
                  wellconc = as.numeric(wellconc),
                  total_volume = ifelse(version == "v2", 
                                        drugplates*wells*50*1.1,
                                        drugplates*wells*25*1.1),
                  lysate = 0.99*total_volume,
                  control_drug_ul = 0,
                  control_dil_ul = 0.01*total_volume) %>%
    dplyr::arrange(desc(wellconc)) %>%
    # change for serial dilutions, start conc and volume
    dplyr::mutate(intconc = dplyr::lag(workingconc)/1000,
                  intconc = ifelse(is.na(intconc), drugstock, intconc),
                  # total_volume = 0.01*total_volume) %>%
                  total_volume = ifelse(wellconc != max(wellconc), 0.01*total_volume, total_volume)) %>%
    dplyr::mutate(drug_ul = workingconc*total_volume/(intconc*1000),
                  drug_dilute = ceiling((drug_ul * intconc) / drugstock),
                  dil_dilute = (drug_dilute * drugstock/intconc) - drug_dilute,
                  dil_ul = control_dil_ul - drug_ul) %>%
    dplyr::arrange(wellconc) %>%
    dplyr::mutate(total_volume = min(total_volume))

for(i in 2:nrow(dose)) {
    j <- i-1
    if this is the highest concentration, change the volume to add lysate
    if(dose$wellconc[i] == max(dose$wellconc)) {
        new_volume <- (dose$control_dil_ul[i] + dose$drug_ul[j])*100
    } else {
        new_volume <- dose$control_dil_ul[i] + dose$drug_ul[j]
    }

    # recalculate with new volume
    dose$drug_ul[i] = dose$workingconc[i]*new_volume/(dose$intconc[i]*1000)
    dose$drug_dilute[i] = ceiling((dose$drug_ul[i] * dose$intconc[i]) / dose$drugstock[i])
    dose$dil_dilute[i] = (dose$drug_dilute[i] * dose$drugstock[i]/dose$intconc[i]) - dose$drug_dilute[i]
    
    if(dose$wellconc[i] == max(dose$wellconc)) {
        dose$dil_ul[i] = new_volume/100 - dose$drug_ul[i]
    } else {
        dose$dil_ul[i] = new_volume - dose$drug_ul[i]
    }
}



    # update total volume to make sure you have enough for dilutions
    dplyr::mutate(dilution_volume = ifelse(wellconc == max(wellconc), dplyr::lead(drug_ul) + control_dil_ul, total_volume + dplyr::lead(drug_ul)))




    dplyr::mutate(drug_ul = workingconc*total_volume/(intconc*1000),
                  drug_dilute = ceiling((drug_ul * intconc) / drugstock),
                  dil_dilute = (drug_dilute * drugstock/intconc) - drug_dilute,
                  dil_ul = control_dil_ul - drug_ul)


# no changes?
drugdilute <- dose %>%
    dplyr::group_by(drug) %>%
    dplyr::mutate(dilution = paste0("1:", drugstock / intconc),
                  drug_need = drug_ul / (drugstock / intconc),
                  total_drug = sum(drug_need),
                  dilution = ifelse(dilution == "1:1", "Do not dilute.", dilution)) %>%
    dplyr::select(drug, diluent, `drugstock (mM)` = drugstock, `intermediate_conc (mM)` = intconc, dilution_factor = dilution, total_drug) %>%
    dplyr::distinct(drug, .keep_all = T) %>%
    dplyr::select(Drug = drug, Diluent = diluent, `Stock Conc. (mM)` = `drugstock (mM)`, `Working Conc. (mM)` = `intermediate_conc (mM)`,
                  Dilution = dilution_factor, `Total Drug (uL)` = total_drug)

# distinct by drug only - but keep any dilutions if there are

# CHANGE THIS
drugplate <- dose %>%
    dplyr::mutate(intconc = ifelse(wellconc == max(wellconc), "Stock concentration", intconc * 1000),
                  intconc = ifelse(wellconc == 0, "NA", intconc)) %>%
    # dplyr::mutate(start_conc_uM = paste0("1:", drugstock/intconc),
                  # start_conc_uM = ifelse(start_conc_uM == "1:1", "Stock Concentration", start_conc_uM)) %>%
    dplyr::select(condition = drug, diluent, start_conc_uM = intconc, concentration_uM = wellconc, plates = drugplates, lysate, dil_ul, drug_ul, mix_to_well_uL = control_dil_ul)  %>%
    dplyr::arrange(desc(concentration_uM)) %>%
    dplyr::mutate(concentration_uM = as.character(concentration_uM) )

# split by dilution factor
# allplates <- split(drugplate, drugplate$dilution_factor)
allplates <- drugplate %>%
    dplyr::mutate(drug_ul = round(drug_ul, digits = 2),
                  dil_ul = round(dil_ul, digits = 2),
                  lysate = round(lysate, digits = 2)) %>%
    dplyr::select(Condition = condition, Diluent = diluent, `Start Concentration (uM)` = start_conc_uM, `Concentration (uM)` = concentration_uM, Plates = plates, `Food (uL)` = lysate,
                  `Diluent (uL)` = dil_ul, `Drug (uL)` = drug_ul, `Drug + Diluent to food (uL)` = mix_to_well_uL)
