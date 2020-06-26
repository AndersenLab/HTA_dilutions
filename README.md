# High-throughput assay dilutions

An R shiny web app was developed to calculate drug dilutions for the high-throughput drug-response assays and can be found [here](https://katiesevans9.shinyapps.io/HTA_V3_dilution/). 

### Setup
The user selects:
- type of assay (V2 or V3)
- type of food (lysate or live bacteria)
- whether or not they will be performing a dose response with several drug doses or a standard assay with just one drug dose
- minimum volume to pipette 
- number of drugs in the assay

Once these settings are complete, the user will press the “Setup Complete” button and a new window will appear to the right under the “setup” tab. 

(insert image here)

In the setup tab, the user will input the information for the assay including:
- drug name 
- drug stock concentration 
- drug final concentration(s) 
- name of diluent (control) 
- number of plates/wells to be run in the assay

(insert image here)

When complete, the user pushes the “Calculate” button and the application automatically calculates the drug dilutions based on pre-set specifications and transfers the user to the “dilution” tab. Here, the application details the necessary dilutions for 1) food, 2) drug stock and 3) drug plate dilutions. A printable PDF can be generated by clicking the “Download Dilutions” button at the bottom.

(insert image here)
