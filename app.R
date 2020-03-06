
library(shiny)
library(shinythemes)
library(tidyverse)
library(gsheet)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
    theme = shinythemes::shinytheme('yeti'),
    
    # Application title
    titlePanel("HTA Dilutions App"),
    
    # sidebar layout
    sidebarLayout(
        sidebarPanel(
            
            h2("Setup"),
            
            # V2 or V3 assay
            shiny::radioButtons(inputId = "version", label = "HTA Version:", choiceNames = c("v2 (96h)", "v3 (48h)"), choiceValues = c("v2", "v3"), inline = T),
            
            # bacteria or lysate
            shiny::radioButtons(inputId = "food", label = "Food type:", choiceNames = c("Lysate", "HB101 (OD 100)"), choiceValues = c("Lysate", "HB101 (OD 100)"), inline = T),
            
            # dose response? Default is false
            shiny::checkboxInput(inputId = "dose", label = "Dose response?", value = F),
            
            # minimum pipette volume
            shiny::numericInput(inputId = "min", label = "Minimum volume to pipette:", value = 0.5, min = 0, max = NA),
            
            # choose number of drugs
            shiny::numericInput(inputId = "drugs", label = "Number of drugs:", value = 1, min = 1, max = NA),
            
            # go button
            shiny::actionButton(inputId = "go", label = "Setup Complete")
            
            
        ),
        
        mainPanel(
            
            # create tabs
            tabsetPanel(id = "dilution_tabs",
                        type = "tabs",
                        tabPanel("setup",
                                 tagList(
                                     # Add link to google sheet with drug information
                                     h6(em("Not sure what concentration your drug is at? Check out this", 
                                           a("Google Sheet", href = "https://docs.google.com/spreadsheets/d/1nS3_Ra2XfguACxCdp2Ni4vUDm2dwsvP_qwumqhzF0tE/edit#gid=211787800"))),
                                     
                                     # based on number of drugs, allow input information for drug, condition, control, dose, plates, etc.
                                     shiny::uiOutput("setup"),
                                     
                                     # add action button
                                     shiny::actionButton(inputId = "drug_done", label = "Calculate")
                                
                                 )),
                        tabPanel("dilutions",
                                 tagList(
                                     # lysate setup
                                     shiny::uiOutput("lysate_setup"),
                                     
                                     # drug dilutions
                                     h2("Drug Dilutions"),
                                     shiny::tableOutput("drugdilutions"),
                                     
                                     # Add link to google sheet with drug information
                                     h6(em("Not sure how many drug tubes do you need? Check out this", 
                                           a("Google Sheet", href = "https://docs.google.com/spreadsheets/d/1nS3_Ra2XfguACxCdp2Ni4vUDm2dwsvP_qwumqhzF0tE/edit#gid=211787800"))),
                                     br(),
                                     
                                     # dilutions dataframe
                                     h2("Plate Dilutions"),
                                     conditionalPanel(condition = "input.version == 'v2'", 
                                                      h5("Add 50 uL food + drug to wells.")),
                                     conditionalPanel(condition = "input.version == 'v3'", 
                                                      h5("Add 25 uL food + drug to wells with 50 uL worms.")),
                                     shiny::uiOutput("dilutions"),
                                     br(),
                                     
                                     # download button
                                     shiny::downloadButton("download", "Download Dilutions")
                                 )))

        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    # read in drug inventory
    tubes <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1nS3_Ra2XfguACxCdp2Ni4vUDm2dwsvP_qwumqhzF0tE/edit#gid=211787800")
    names(tubes) <- tubes[1,]
    tubes <- tubes[2:nrow(tubes),]
    
    # when go button is pushed, do this
    shiny::observeEvent(input$go, {
        setup_plates()
    })
    
    # setup plates
    setup_plates <- reactive({
        
        # how many drugs?
        numdrugs <- input$drugs
        
        # allow user to input controls
        output$setup <- shiny::renderUI({
            
            # update when button is pushed
            input$go
            
            if(input$dose) {
                lapply(1:numdrugs, function(i) {
                    inputPanel(
                        tagList(
                            h2(glue::glue("Drug {i}")),
                            
                            # drug name
                            shiny::textInput(inputId = glue::glue("drug_name{i}"), label = "Drug:", value = NA),
                            
                            # drug stock concentration
                            shiny::numericInput(inputId = glue::glue("drug_stock{i}"), label = "Drug stock concentration (mM):", min = 0, value = NA),
                            
                            # diluent
                            shiny::radioButtons(inputId = glue::glue("control_name{i}"), label = "Control:", choices = c("DMSO", "Water", "Other")),
                            
                            # drug concentration
                            shiny::textInput(inputId = glue::glue("drug_dose{i}"), label = "Drug concentrations (uM): (e.g. 0,5,10,15,20)"),
                            
                            # number of plates
                            shiny::numericInput(inputId = glue::glue("drugplates{i}"), label = "Number of drug plates:", min = 1, value = 1),
                            
                            # number of wells per condition per plate
                            shiny::numericInput(inputId = glue::glue("wells{i}"), label = "Number of wells per condition/plate:", min = 1, max = 96, value = 16)
                        )
                    )
                })
            } else {
                lapply(1:numdrugs, function(i) {
                    inputPanel(
                        tagList(
                            h2(glue::glue("Drug {i}")),
                            
                            # drug name
                            shiny::textInput(inputId = glue::glue("drug_name{i}"), label = "Drug:", value = NA),
                            
                            # drug stock concentration
                            shiny::numericInput(inputId = glue::glue("drug_stock{i}"), label = "Drug stock concentration (mM):", min = 0, value = NA),
                            
                            # diluent
                            shiny::radioButtons(inputId = glue::glue("control_name{i}"), label = "Control:", choices = c("DMSO", "Water", "Other")),
                            
                            # drug concentration
                            shiny::numericInput(inputId = glue::glue("drug_dose{i}"), label = "Drug concentration (uM):", min = 0, value = NA),
                            
                            # number of plates
                            shiny::numericInput(inputId = glue::glue("drugplates{i}"), label = "Number of drug plates:", min = 1, value = 1),
                            
                            # number of control plates
                            shiny::numericInput(inputId = glue::glue("controlplates{i}"), label = "Number of control plates:", min = 0, value = 1),
                            
                            # number of wells per condition per plate
                            shiny::numericInput(inputId = glue::glue("wells{i}"), label = "Number of wells per condition/plate:", min = 1, max = 96, value = 96)
                        )
                    )
                })
            }
        })
    })
    
    
    # when go button is pushed, do this
    shiny::observeEvent(input$drug_done, {
        dosetable <- calculate_dose()
        
        if(input$dose) {
            lapply(1:input$drugs, function(d) {
                output[[paste0('T', d)]] <- shiny::renderTable({
                    dosetable[[1]] %>%
                        dplyr::mutate(Condition = as.character(Condition)) %>%
                        dplyr::filter(Condition == input[[glue::glue("drug_name{d}")]])
                })
            })
            
            output$dilutions <- renderUI({
                tagList(lapply(1:input$drugs, function(i) {
                    shiny::tableOutput(paste0('T', i))
                }))
            })
        } else {
            # output dataframe
            output$dilutions <- shiny::renderTable({
                dosetable[[1]]
            })
        }
        
        # also switch tabs
        shiny::updateTabsetPanel(session, "dilution_tabs", selected = "dilutions")

    })
    
    # calculate dose of drugs
    calculate_dose <- reactive({
        
        # how many drugs?
        numdrugs <- input$drugs
        
        if(input$dose) {
            # create dataframe
            alldose <- NULL
            for(i in 1:numdrugs) {
                df <- data.frame(
                    drug = input[[glue::glue("drug_name{i}")]],
                    drugstock = input[[glue::glue("drug_stock{i}")]],
                    diluent = input[[glue::glue("control_name{i}")]],
                    wellconc = input[[glue::glue("drug_dose{i}")]],
                    version = input$version,
                    drugplates = input[[glue::glue("drugplates{i}")]],
                    controlplates = 0,
                    wells = input[[glue::glue("wells{i}")]],
                    intconc = input[[glue::glue("drug_stock{i}")]]
                )
                
                # check tenfold dilutions
                min <- input$min
                counts <- NULL
                dilutions <- c(1, 10, 100, 1000, 10000)
                for(i in dilutions) {
                    dose <- df %>%
                        tidyr::separate_rows(wellconc, sep = ",") %>%
                        dplyr::mutate(version = as.character(version),
                                      wellconc = as.numeric(wellconc),
                                      workingconc = ifelse(version == "v2", wellconc, wellconc*3),
                                      total_volume = ifelse(version == "v2", 
                                                            drugplates*wells*50*1.1,
                                                            drugplates*wells*25*1.1),
                                      intconc = intconc / i,
                                      lysate = 0.99*total_volume,
                                      drug_ul = workingconc*total_volume/(intconc*1000),
                                      drug_dilute = ceiling((drug_ul * intconc) / drugstock),
                                      dil_dilute = (drug_dilute * drugstock/intconc) - drug_dilute,
                                      dil_ul = (0.01*total_volume) - drug_ul,
                                      control_drug_ul = 0,
                                      control_dil_ul = 0.01*total_volume) %>%
                        dplyr::filter(drug_ul > min)
                    counts <- append(counts, nrow(dose))
                }
                
                # determine which intermediate concentration
                dil <- dilutions[which(counts == max(counts))][1]
                
                # use dilution
                dose <- df %>%
                    tidyr::separate_rows(wellconc, sep = ",") %>%
                    dplyr::mutate(version = as.character(version),
                                  wellconc = as.numeric(wellconc),
                                  workingconc = ifelse(version == "v2", wellconc, wellconc*3),
                                  total_volume = ifelse(version == "v2", 
                                                        drugplates*wells*50*1.1,
                                                        drugplates*wells*25*1.1),
                                  intconc = intconc / dil,
                                  lysate = 0.99*total_volume,
                                  drug_ul = workingconc*total_volume/(intconc*1000),
                                  drug_dilute = ceiling((drug_ul * intconc) / drugstock),
                                  dil_dilute = (drug_dilute * drugstock/intconc) - drug_dilute,
                                  dil_ul = (0.01*total_volume) - drug_ul,
                                  control_drug_ul = 0,
                                  control_dil_ul = 0.01*total_volume)
                
                # check to see if the higher concentrations exceed 1% control
                toohigh <- dose %>%
                    dplyr::filter(dil_ul < 0)
                
                justright <- dose %>%
                    dplyr::filter(dil_ul > 0)
                
                if(nrow(toohigh) > 0) {
                    counts <- NULL
                    dilutions <- c(1, 10, 100, 1000, 10000)
                    for(i in dilutions) {
                        high <- toohigh %>%
                            dplyr::select(names(df))%>%
                            dplyr::mutate(intconc = drugstock) %>%
                            dplyr::mutate(version = as.character(version),
                                          wellconc = as.numeric(wellconc),
                                          workingconc = ifelse(version == "v2", wellconc, wellconc*3),
                                          total_volume = ifelse(version == "v2", 
                                                                drugplates*wells*50*1.1,
                                                                drugplates*wells*25*1.1),
                                          intconc = intconc / i,
                                          lysate = 0.99*total_volume,
                                          drug_ul = workingconc*total_volume/(intconc*1000),
                                          drug_dilute = ceiling((drug_ul * intconc) / drugstock),
                                          dil_dilute = (drug_dilute * drugstock/intconc) - drug_dilute,
                                          dil_ul = (0.01*total_volume) - drug_ul,
                                          control_drug_ul = 0,
                                          control_dil_ul = 0.01*total_volume) %>%
                            dplyr::filter(drug_ul > min)
                        counts <- append(counts, nrow(high))
                    }
                    
                    # determine which intermediate concentration
                    dil <- dilutions[which(counts == max(counts))][1]
                    
                    # use dilution for high rows
                    high <- toohigh %>%
                        dplyr::select(names(df))%>%
                        dplyr::mutate(intconc = drugstock) %>%
                        dplyr::mutate(version = as.character(version),
                                      wellconc = as.numeric(wellconc),
                                      workingconc = ifelse(version == "v2", wellconc, wellconc*3),
                                      total_volume = ifelse(version == "v2", 
                                                            drugplates*wells*50*1.1,
                                                            drugplates*wells*25*1.1),
                                      intconc = intconc / dil,
                                      lysate = 0.99*total_volume,
                                      drug_ul = workingconc*total_volume/(intconc*1000),
                                      drug_dilute = ceiling((drug_ul * intconc) / drugstock),
                                      dil_dilute = (drug_dilute * drugstock/intconc) - drug_dilute,
                                      dil_ul = (0.01*total_volume) - drug_ul,
                                      control_drug_ul = 0,
                                      control_dil_ul = 0.01*total_volume)
                    
                    dose <- rbind(justright, high)
                }
                counts <- NULL
                toohigh <- NULL
                alldose <- rbind(alldose, dose)
            }
            
            dose <- alldose
            
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
            
            
            drugplate <- dose %>%
                dplyr::mutate(dilution_factor = paste0("1:", drugstock/intconc),
                              dilution_factor = ifelse(dilution_factor == "1:1", "Stock Concentration", dilution_factor)) %>%
                dplyr::select(condition = drug, diluent, dilution_factor, concentration_uM = wellconc, plates = drugplates, lysate, dil_ul, drug_ul)  %>%
                dplyr::arrange(concentration_uM) %>%
                dplyr::mutate(concentration_uM = as.character(concentration_uM))
            
            # split by dilution factor
            # allplates <- split(drugplate, drugplate$dilution_factor)
            allplates <- drugplate %>%
                dplyr::mutate(drug_ul = round(drug_ul, digits = 2),
                              dil_ul = round(dil_ul, digits = 2),
                              lysate = round(lysate, digits = 2)) %>%
                dplyr::select(Condition = condition, Diluent = diluent, `Dilution Factor` = dilution_factor, `Concentration (uM)` = concentration_uM, Plates = plates, `Food (uL)` = lysate,
                              `Diluent (uL)` = dil_ul, `Drug (uL)` = drug_ul)
            
            # calculate lysate needs
            plates <- ceiling(sum(dose$wells) / 96)
            
        } else {
            # create dataframe
            df <- NULL
            for(i in 1:numdrugs) {
                df2 <- data.frame(
                    drug = input[[glue::glue("drug_name{i}")]],
                    drugstock = input[[glue::glue("drug_stock{i}")]],
                    diluent = input[[glue::glue("control_name{i}")]],
                    wellconc = input[[glue::glue("drug_dose{i}")]],
                    version = input$version,
                    drugplates = input[[glue::glue("drugplates{i}")]],
                    controlplates = input[[glue::glue("controlplates{i}")]],
                    wells = input[[glue::glue("wells{i}")]],
                    intconc = input[[glue::glue("drug_stock{i}")]]
                )
                df <- rbind(df, df2)
            }
            
            dose <- df %>%
                dplyr::mutate(version = as.character(version),
                              workingconc = ifelse(version == "v2", wellconc, wellconc*3),
                              total_volume = ifelse(version == "v2", 
                                                    drugplates*wells*50*1.1,
                                                    drugplates*wells*25*1.1),
                              lysate = 0.99*total_volume,
                              drug_ul = workingconc*total_volume/(intconc*1000),
                              drug_dilute = ceiling((drug_ul * intconc) / drugstock),
                              dil_dilute = (drug_dilute * drugstock/intconc) - drug_dilute,
                              dil_ul = (0.01*total_volume) - drug_ul,
                              control_drug_ul = 0,
                              control_dil_ul = 0.01*total_volume)
            
            # Check to make sure dilution is enough, should not pipette less than the minimum amount
            min <- input$min
            tooSmall <- dose %>%
                dplyr::filter(drug_ul < min, drug_ul > 0)
            dose2 <- dose
            
            while(nrow(tooSmall) > 0) {
                dose2 <- dose2 %>%
                    dplyr::select(names(df))
                oldconcs <- dose2 %>%
                    dplyr::filter(!(drug %in% unique(tooSmall$drug)))
                concs2 <- dose2 %>%
                    dplyr::filter(drug %in% unique(tooSmall$drug)) %>%
                    dplyr::mutate(intconc = intconc / 10)
                newconcs <- rbind(oldconcs, concs2)
                
                dose2 <- newconcs %>%
                    dplyr::mutate(version = as.character(version),
                                  workingconc = ifelse(version == "v2", wellconc, wellconc*3),
                                  total_volume = ifelse(version == "v2", 
                                                        drugplates*wells*50*1.1,
                                                        drugplates*wells*25*1.1),
                                  lysate = 0.99*total_volume,
                                  drug_ul = workingconc*total_volume/(intconc*1000),
                                  dil_ul = (0.01*total_volume) - drug_ul,
                                  drug_dilute = ceiling((drug_ul * intconc) / drugstock),
                                  dil_dilute = (drug_dilute * drugstock/intconc) - drug_dilute,
                                  control_drug_ul = 0,
                                  control_dil_ul = 0.01*total_volume)
                
                tooSmall <- dose2 %>%
                    dplyr::filter(drug_ul < min, drug_ul > 0)
            }
            
            dose <- dose2
            
            drugdilute <- dose %>%
                dplyr::group_by(drug) %>%
                dplyr::mutate(dilution = paste0("1:", drugstock / intconc),
                              drug_need = drug_ul / (drugstock / intconc),
                              total_drug = sum(drug_need),
                              dilution = ifelse(dilution == "1:1", "Do not dilute.", dilution)) %>%
                dplyr::select(drug, diluent, `drugstock (mM)` = drugstock, `intermediate_conc (mM)` = intconc, dilution_factor = dilution, total_drug) %>%
                dplyr::distinct() %>%
                dplyr::select(Drug = drug, Diluent = diluent, `Stock Conc. (mM)` = `drugstock (mM)`, `Working Conc. (mM)` = `intermediate_conc (mM)`,
                              Dilution = dilution_factor, `Total Drug (uL)` = total_drug)
            
            drugplate <- dose %>%
                dplyr::select(condition = drug, diluent, concentration_uM = wellconc, plates = drugplates, lysate, dil_ul, drug_ul) %>%
                dplyr::mutate(concentration_uM = as.character(concentration_uM))
            
            controlplate <- dose %>%
                dplyr::select(condition = diluent, diluent, concentration_uM = wellconc, plates = controlplates, lysate, dil_ul = control_dil_ul, drug_ul = control_drug_ul, wells) %>%
                dplyr::mutate(diluent = "None",
                              concentration_uM = "1%",
                              total_volume = ifelse(input$version == "v2", 
                                                    plates*wells*50*1.1,
                                                    plates*wells*25*1.1),
                              lysate = 0.99*total_volume,
                              drug_ul = 0,
                              dil_ul = 0.01*total_volume) %>%
                dplyr::filter(plates > 0) %>%
                dplyr::select(condition, diluent, concentration_uM , plates, lysate, dil_ul, drug_ul)
            
            allplates <- drugplate %>%
                dplyr::full_join(controlplate) %>%
                dplyr::mutate(drug_ul = round(drug_ul, digits = 2),
                              dil_ul = round(dil_ul, digits = 2),
                              lysate = round(lysate, digits = 2)) %>%
                dplyr::select(Condition = condition, Diluent = diluent, `Concentration (uM)` = concentration_uM, Plates = plates, `Food (uL)` = lysate,
                              `Diluent (uL)` = dil_ul, `Drug (uL)` = drug_ul)
            
            # calculate lysate needs
            plates <- sum(dose$drugplates) + sum(dose$controlplates)
        }
        
        # lysate or bacteria?
        if(input$food == "Lysate") {
            # lysate
            if(input$version == "v2") {
                # make 10 mg/mL for 96 hours
                #If lysate is exact, add another tube
                lysate <- sum(allplates$`Food (uL)`)
                if((lysate %% 10000) == 0) {
                    lysate <- lysate + 10000
                }
                
                lysate_mix <- (paste0("Total lysate mix needed: ", round(lysate / 1000, digits = 2), " mL (10 mg/mL), actually make: ", round(ceiling(lysate / 10000)*10, digits = 2), " mL"))
                lysate_tubes <- (paste0("Number of lysate tubes needed: ", ceiling(lysate / 10000)))
                lysate <- round(ceiling(lysate / 10000) *10000, digits = 2)
                lysate_K <- (paste0("Add ",  round(lysate / 10000 * 9, digits = 2), " mL of K medium to lysate."))
                kan <- (paste0("Add ", round(50*lysate/80000, digits = 2), " uL of Kanamycin to lysate mix."))
                
                lysate_dilutions <- c(lysate_mix, lysate_tubes, lysate_K, kan)
            } else if(input$version == "v3") {
                # make 15 mg/mL which will be diluted 1:3 in the well
                #If lysate is exact, add another tube
                lysate <- sum(allplates$`Food (uL)`)
                if((lysate %% 6666) == 0) {
                    lysate <- lysate + 6666
                }
                
                lysate_mix <- (paste0("Total lysate mix needed: ", round(lysate / 1000, digits = 2), " mL (15 mg/mL), actually make: ", round(ceiling(lysate / 6666)*6.666, digits = 2), " mL"))
                lysate_tubes <- (paste0("Number of lysate tubes needed: ", round(ceiling(lysate / 6666), digits = 2)))
                lysate <- round(ceiling(lysate / 6666) *6666, digits = 2)
                lysate_K <- (paste0("Add ",  round(lysate/1000 - ceiling(lysate / 6666), digits = 2), " mL of K medium to lysate."))
                kan <- (paste0("Add ", round(3*50*lysate/80000, digits = 2), " uL of Kanamycin to lysate mix."))
                
                lysate_dilutions <- c(lysate_mix, lysate_tubes, lysate_K, kan)
            }
        } else {
            # bacteria -- OD 100 to OD 10
            if(input$version == "v2") {
                # make OD 20 for 96 hours
                #If lysate is exact, add another tube
                lysate <- sum(allplates$`Food (uL)`)
                if((lysate %% 5000) == 0) {
                    lysate <- lysate + 5000
                }
                
                lysate_mix <- (paste0("Total HB101 needed: ", round(lysate / 1000, digits = 2), " mL (OD 20), actually make: ", round(ceiling(lysate / 5000)*5, digits = 2), " mL"))
                lysate_tubes <- (paste0("Number of HB101 tubes needed: ", ceiling(lysate / 5000)))
                lysate <- round(ceiling(lysate / 5000) *5000, digits = 2)
                lysate_K <- (paste0("Add ",  round(lysate / 5000 * 5, digits = 2) - (lysate / 5000), " mL of K medium to HB101"))
                kan <- (paste0("Add ", round(50*lysate/80000, digits = 2), " uL of Kanamycin to HB101."))
                
                lysate_dilutions <- c(lysate_mix, lysate_tubes, lysate_K, kan)
            } else if(input$version == "v3") {
                # make OD 30 which will be diluted 1:3 in the well
                #If bacteria is exact, add another tube
                lysate <- sum(allplates$`Food (uL)`)
                if((lysate %% 3333) == 0) {
                    lysate <- lysate + 3333
                }
                
                lysate_mix <- (paste0("Total HB101 needed: ", round(lysate / 1000, digits = 2), " mL (OD 30), actually make: ", round(ceiling(lysate / 3333)*3.3333, digits = 2), " mL"))
                lysate_tubes <- (paste0("Number of HB101 tubes needed: ", round(ceiling(lysate / 3333), digits = 2)))
                lysate <- round(ceiling(lysate / 3333) *3333, digits = 2)
                lysate_K <- (paste0("Add ",  round(lysate/1000 - ceiling(lysate / 3333), digits = 2), " mL of K medium to HB101"))
                kan <- (paste0("Add ", round(3*50*lysate/80000, digits = 2), " uL of Kanamycin to HB101."))
                
                lysate_dilutions <- c(lysate_mix, lysate_tubes, lysate_K, kan)
            }
        }

        
        return(list(allplates, lysate_dilutions, drugdilute))

        
    })
    
    # lysate setup
    output$lysate_setup <- shiny::renderUI({
        
        dosetable <- calculate_dose()
        
        # drug dilutions
        output$drugdilutions <- shiny::renderTable({
            dosetable[[3]]
        })

        lysate_instructions <- dosetable[[2]]

        tagList(
            h2("Food Setup"),
            h6(glue::glue("Version: {input$version}")),
            h6(glue::glue("Food: {input$food}")),
            h6(lysate_instructions[1]),
            h6(lysate_instructions[2]),
            h6(lysate_instructions[3]),
            h6(lysate_instructions[4]),
            br()
            )
        
        
    })
    
    # download report
    output$download <- downloadHandler(
        # For PDF output, change this to "report.pdf"
        filename = "report.pdf",
        content = function(file) {
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            tempReport <- file.path(tempdir(), "report.Rmd")
            file.copy("report.Rmd", tempReport, overwrite = TRUE)

            # load inputs into global environment for markdown            
            inputs <- names(input)
            for(i in 1:length(inputs)) {
                assign(inputs[i], input[[inputs[i]]])
            }
            
            
            # Knit the document, passing in the `params` list, and eval it in a
            # child of the global environment (this isolates the code in the document
            # from the code in this app).
            rmarkdown::render(tempReport, output_file = file)
            
        }
    )
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)

