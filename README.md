# **DALY Calculation app for neurocysticercosis**

This Shiny App is to calculate disability-adjusted life years (DALYs) using transitional disease model of neurocysticercosis. The interface of the app is as follows:

![App interface](https://github.com/Baty2025/DALY_dModel_4.3_app/blob/main/WWW/exmpl_interface.PNG)

The concepts behind this calculation are as shown in following computational disease model: ![Model](https://github.com/Baty2025/DALY_dModel_4.3_app/blob/main/WWW/model_4.3.jpg)

# **How to use this app**

## **Clone the repository**

```         
 git clone https://github.com/Baty2025/DALY_dModel_4.3_app.git 
```

**Alternatively,**

-   Download all files as ZIP file from here: <https://github.com/Baty2025/DALY_dModel_4.3_app>
-   Un-zip and Keep all files within **same** folder.

## Setup

-   Open R project: `DALY_dModel_4.3_app.Rproj`

-   Before running the app, make sure you have `renv` installed:

```         
install.packages("renv") 
renv::restore()
```

This will install all required packages for the app.

## Launch app

-   Run following code from **Console**

```         
shiny::runApp("app.R")
```

This will bring you to the app's interface.

## Calculate DALY

-   You will see all of required data have already been auto populated. These data are collected for Indian context and given in the app as an example data set.
-   If you click, "Calculate Burden" button, the app will give you the summary results and probabilistic sensitivity analysis in the main pane. You will see some more analysis in subsequent panes.
-   You can change any values from left side bar and recalculate the burden.

# Why this app?

Typically, NCC-associated DALYs are calculated using an attributional model. *dModel 4.3* is a hypothetical framework designed to work when highly granular data on NCC prevalence in the population are available. However, such detailed estimates remain out of reach because there is currently no rapid, non-imaging test to measure NCC prevalence in the general population. Brain imaging—the only reliable diagnostic tool—cannot ethically be applied to asymptomatic individuals.

Given this limitation, this app allows users to first calculate DALYs using the data they currently have. Then, they can adjust different neurological parameters to explore how these changes impact the overall DALY estimates. This interactive approach makes it easier to visualize and understand how NCC-associated DALY estimates vary under different assumptions.
