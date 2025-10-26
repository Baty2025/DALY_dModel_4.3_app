# **DALY Calculation app for neurocysticercosis**

This Shiny application calculates **Disability-Adjusted Life Years (DALYs)** using a **transitional disease model** for neurocysticercosis (NCC).

### App Interface

![App interface](https://github.com/Baty2025/DALY_dModel_4.3_app/blob/main/WWW/exmpl_interface.PNG)

### Computational Disease Model

The underlying concepts and disease progression structure are represented in the following computational model: ![Model](https://github.com/Baty2025/DALY_dModel_4.3_app/blob/main/WWW/model_4.3.jpg)

You can create similar model [here](https://baty2025.shinyapps.io/Code/)

# **How to use this app**

### **1. Clone the repository**

```         
 git clone https://github.com/Baty2025/DALY_dModel_4.3_app.git 
```

**Alternatively,**

-   Download all files as ZIP file from here: <https://github.com/Baty2025/DALY_dModel_4.3_app>
-   Extract the ZIP file and keep **all files in the same folder**

### 2. Set up the environment

-   Open R project file: `DALY_dModel_4.3_app.Rproj`

-   Before running the app, make sure `renv` is installed, then restore the project environment:

```         
install.packages("renv") 
renv::restore()
```

This will automatically install all the required packages used in the app.

### 3. Launch app

-   Run the following command in your R Console:

```         
shiny::runApp("app.R")
```

This will launch the Shiny app in your default web browser.

## Calculate DALY

-   The app comes pre-loaded with **example data** representing the Indian context..

-   Click the **“Calculate Burden”** button to generate DALY estimates and probabilistic sensitivity analysis results. Example outputs:

    ![Results](https://github.com/Baty2025/DALY_dModel_4.3_app/blob/main/WWW/results.PNG)

-   Additional analyses can be viewed in subsequent tabs or panes.

-   You can modify any input values from the sidebar and recalculate the burden interactively.

# Why this app?

Typically, NCC-associated DALYs are calculated using an attributional model. *dModel 4.3* is a hypothetical framework designed to work when highly granular data on NCC prevalence in the population are available. However, such detailed estimates remain out of reach because there is currently no rapid, non-imaging test to measure NCC prevalence in the general population. Brain imaging—the only reliable diagnostic tool—cannot ethically be applied to asymptomatic individuals.

Given this limitation, this app allows users to first calculate DALYs using the data they currently have. Then, they can adjust different neurological parameters to explore how these changes impact the overall DALY estimates. This interactive approach makes it easier to visualize and understand how NCC-associated DALY estimates vary under different assumptions.

## **Citation**

If you use this app in your research, please cite:

> Mohammad Shah Jalal, University of Montreal (2025). *DALY_dModel_4.3_app: A Shiny Application for Estimating the Burden of Neurocysticercosis.*\
> GitHub repository: <https://github.com/Baty2025/DALY_dModel_4.3_app>

## **Contact**

For questions, suggestions, or bug reports, please contact:\
📧 **shah.jalal.baty\@gmail.com**
