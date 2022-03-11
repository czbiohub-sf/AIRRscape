# AIRRscape: an interactive web tool for exploring B-cell receptor repertoires and antibody responses  

## Dependencies 
Have RStudio downloaded: https://www.rstudio.com/products/rstudio/download/

## Running AIRRscape
To run AIRRscape, clone the repo and open the app.R file in your RStudio, then click "Run App". As a Shiny app, it can run as a window of RStudio, or as a tab in a web browser (recommended).

Basic instructions are shown on the right hand side of the window. When you first start AIRRscape you will see a set of heatmaps - to choose between datasets shown select the dataset from the options in the top selectable list. Datasets of multiple repertoires can be visualized either as separate heatmap panels, or as a single combined heatmap (labeled 'combined').

The heatmaps split antibody repertoires into bins based on their germline V-gene and J-gene assignments (x-axis), and their CDR3 lengths (y-axis). The second selectable list on the right lets you select the fill color for the bins: either by average somatic hypermutation (SHM), maximum SHM, or by percent of total (i.e. what percent of antibodies in that panel are found in that bin).

You can interactively explore the data in the heatmap: you can 1) hover over a bin to get a popup window showing some basic stats, you can 2) create a bounding box to get a table of all the antibodies within the box, and finally you can 3) click on a single bin to get a table of antibodies in just that bin. The bounding box table and clicked table will appear in two different spots below the heatmaps. From either table, you can sort the table by any of the columns. You can also use the search bar to search for an antibody of interest, or even to limit the table to a single attribute (i.e. one germline).

From the table of antibodies after clicking on a selected bin, you can download all or selected antibodies in the chosen bin, download the distance matrix of all antibodies, or create topologies of selected antibodies based on their CDR3 amino acid motifs. There are multiple options for making topologies. The most straightforward is to manually choose between 3 and ~500 antibodies from the table, and choose either the 'NJ' or 'parsimony' options. A toplogy with your selected antibodies will appear below. Another option is to select just a single antibody of interest, and the choose one of the last four topology options. These will search the entire table, and show you the topology of all CDR# motifs in that bin that are within the chosen identity threshold (based on aa distance). Four thresholds are available ranging from 50% - 100% identity.

The toplogies default to a particular height and width that may not be optimal given the number of antibodies, or the size of your screen. You can adjust the height and width using the sliders just above the topology. The topology tips include the antibody names but also the CDR3 motifs (in Courier New font so as to be aligned), and the V-gene assignments. Lastly, the button below the topology that will save a screenshot of the entire page including heatmaps, tables, and topology.

Note that after exploring one table for some time, it is possible to unwittingly have multiple antibodies selected but not in view - this will affect the topology-making options. A simple solution is to click on another bin, and then click back on the bin of interest thus refreshing the table. Also note that the topologies require the selected antibodies to all have the same CDR3 length - otherwise the calculation will fail with an error:
```{r}
Warning: Error in <-: length of 'dimnames' [1] not equal to array extent
```

## Datasets and code for manuscript
The repo is split into the AIRRscape app and loaded datasets (shinyapp folder), and partially processed datasets & code used to process these loaded datasets for the manuscript (paper_assets folder). There are two R scripts in the paper_assets folder: the airrscape_preprocessing.R script includes the initial R code for combining datasets from repositories in AIRR format (datasets here include all columns). The airrscape_processing.R script includes code that will combine partially processed datasets to the loaded sets directly used in AIRRscape. Importantly, this script also includes a custom function 'shinyprocess' that does this trimming and calculation of AIRRscape-specific columns for any AIRR-formatted repertoire dataset (in a .tab format).

## Tips
We recommend using AIRRscape as a tab on a web browser, as wide as possible. By removing reads with identical CDR3 motifs & germline assignments during processing, each loaded repertoire contains up to 200,000 sequences. The heatmaps of these repertoires do not take more than a few seconds to load, but the final option of all datasets loaded will take longer.
When making topologies there is an upper limit for making these on the fly. We recommend limiting these to no more than 500 sequences. When finding the most closely related CDR3 motifs to an antibody sequence of interest, note that the topology will be limited to 500 sequences. Similarly, when viewing the largest combined datasets, note the number of sequences in your bin of interest. On a typical laptop running RStudio searching a bin of 1000 sequences for the most closely related CDR3 motifs will take about 1 minute. A search of the largest bins in the largest combined datasets will take even longer.

## Citation
To cite AIRRscape in publications, use:
tktk

(modified using RStudio)
