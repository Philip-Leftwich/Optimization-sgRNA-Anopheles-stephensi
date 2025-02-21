# Model tables ====

source("scripts/cutting_analysis.R")
source("scripts/homing_mosaic_analysis.R")

sjPlot::tab_model(homing.model,
                  cut.model.1, 
                  mosaic.model.1,
                  pred.labels=c(
                    "Intercept (Line 1626)",
                    "Cas9 parent (Male)",
                    "Cas9 grandparent (Male)",
                    "Line 1757",
                    "Line 1758",
                    "Line 1759",
                    "Cas9 parent:grandparent",
                    "Parent:1757",
                    "Parent:1758",
                    "Parent:1759",
                    "Grandparent:1757",
                    "Grandparent:1758",
                    "Grandparent: 1759"
                    ),
                  dv.labels=c("Homing against WT",
                              "Homing against KO",
                              "Cutting against KO"),
                  collapse.ci=TRUE,
                  transform=NULL,
                  p.style="stars")

