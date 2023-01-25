![](airrscape_logo_and_name2.png)
![](fig2b_forcomms.png)

## Dependencies & R Session Info
**AIRRscape is available as a web tool at https://airrscape.czbiohub.org.**
To run AIRRscape locally, R & RStudio must be installed. AIRRscape was developed using R v4.2.0 with RStudio v2022.02.2.
Attached base R packages:	_stats, graphics, grDevices, utils, datasets, methods, base_

Other attached packages:
	_shinyscreenshot_0.2.0, bslib_0.3.1, shinycssloaders_1.0.0, phangorn_2.8.1, ape_5.6-2, seqinr_4.2-16, DT_0.23, forcats_0.5.1, stringr_1.4.0, dplyr_1.0.9, purrr_0.3.4, readr_2.1.2, tidyr_1.2.0, tibble_3.1.7, tidyverse_1.3.1, alakazam_1.2.0, ggplot2_3.3.6, shiny_1.7.1_

## Running AIRRscape
**Launch AIRRscape at https://airrscape.czbiohub.org.** Basic instructions are shown in the left sidebar. When first launching, AIRRscape may take 10-20 seconds to load.

To run AIRRscape locally, clone this repo, open the ```app.R``` file in RStudio, and then click "Run App". AIRRscape can run in a window of RStudio or as a tab in a web browser (recommended).

### _AIRRscape_ tab
Upon starting AIRRscape, the main panel has two possible tabs: _AIRRscape_ and _Import Data_, with the _AIRRscape_ tab initially showing a heatmap from a pre-loaded dataset ("SARS-CoV2 mAbs - heavy chains & light chains"). Choose between datasets by selecting from the options in the top dropdown menu of the sidebar. Datasets of multiple repertoires can be visualized either as separate heatmap panels (labeled "IgH"), or as a single combined heatmap (labeled "IgH combined").
To view custom datasets (.tab or .tsv files in AIRR-C format), first use the _Import Data_ tab to upload datasets, convert+combine them, and then download the processed files. Next, upload these files and select "Custom datasets - IgH or TR" or "Custom datasets - IgH or TR combined" from the dropdown menu.

The heatmaps split BCR or TCR repertoires into bins based on either their germline V-gene family + J-gene assignments or individual V-gene assignments (x-axis), chosen via the second dropdown menu in the sidebar, and their CDR3 lengths (y-axis). Select the fill color for the bins using the third dropdown menu in the sidebar: color by 1) average somatic hypermutation (SHM), 2) maximum SHM, or 3) percent of total (i.e. what percent of sequences in the panel are found in that bin).

To interactively explore the data using the heatmap: 1) hover over a bin to get a popup window showing some basic stats, 2) create a bounding box to get a (upper) table of all the sequences within the box, or 3) click on a single bin to get a (lower) table of its sequences. Note that these tables will appear in two different spots below the heatmap, and they can be sorted by any of the columns. From either table, use the search bar to search for an sequence of interest or to limit the table to a single attribute (i.e. one germline V-gene). Note that the upper table displaying sequences within the bounding box can only be used for searches of those sequences. For further functionality, use the lower table after clicking on a single bin.

Available options from the lower table of sequences are: 1) download all or 2) selected sequences in the chosen bin, 3) download the distance matrix of all sequences, or 4) create topologies of selected sequences based on their CDR3 amino acid motifs. There are multiple options for making topologies. The most straightforward is to manually choose between 3 and ~500 sequences from the table, and choose either the 'NJ' or 'parsimony' options. A topology with the selected sequences will appear below. Another option is to select a single BCR or TCR sequence, and then choose one of the last four topology options. These will search the entire table and display the topology of all CDR3 motifs in that bin within the chosen identity threshold (based on amino acid distance). Four thresholds are available, ranging from 50% - 100% identity.

The default height and width of topologies may not be optimal for viewing, depending on the number of sequences or the screen size. Use the sliders to adjust the height and width. The topology tips include the sequence names, as well as the CDR3 motifs and V-gene assignments. Lastly, the button below the topology labeled "Take a screenshot" will save an image of the entire page including the heatmap, table, and topology.

![](airrscape_video_45seconds.gif)

### _Import Data_ tab
To process (convert & combine) datasets, first upload each dataset separately (maximum 6). As long as they follow AIRR-C standards (.tab or .tsv) and include some required columns, they will be automatically converted for viewing in AIRRscape. Metadata is not required; simply enter the name of each dataset, which will be used as labels in each faceted dataset. Then, click the Step 1 combine button to begin processing the datasets. Note the processing occurs in the background, so next just click the Step 2 & 3 download buttons separately to get the processed files (wait for step 2 to finish before pressing the Step 3 button). For larger datasets with >500,000 sequences it might take 1-2 minutes before the download begins. Finally, import these two files just below in the Import Data section for viewing in the _AIRRscape_ tab. The two bottom upload panes can be used to view any converted files processed in the current or a previous session.

Required dataset columns:
 _v_call, j_call, v_identity, junction_aa_

 ![](AIRRscape_importing_customdatasets_95sec.gif)

## Datasets and code for manuscript
This repository is split into sections. The 'shinyapp' folder contains the AIRRscape app and loaded datasets. The 'paper_assets' folder contains partially processed datasets used in the AIRRscape publication and R scripts used for their processing. The ```airrscape_preprocessing.R``` script contains code used to combine datasets from repositories in AIRR-C format (datasets here include all columns). The ```airrscape_processing.R``` script contains code used to convert & combine partially processed datasets into the files directly used in AIRRscape. This script also includes a custom function ```AIRRscapeprocess``` that trims datasets and creates AIRRscape-specific columns for any AIRR-C formatted repertoire dataset. Note that the  _Import Data_ tab includes this function to convert user-uploaded datasets.

## Tips
We recommend using AIRRscape as a tab on a web browser (and as wide as possible). The sizes of processed repertoires are reduced by removing reads with identical CDR3 motifs & germline assignments, and each pre-loaded repertoire contains up to 200,000 sequences. The heatmaps of these repertoires take only few seconds to load, but the largest dataset "SARS-CoV2 HIV & Dengue datasets - IgH combined" will take longer to load.

For topology building, there is an upper limit for calculating these in real-time. We recommend limiting these to no more than 500 sequences. When finding the most closely related CDR3 motifs to a sequence of interest, note that the topology will be limited to 500 sequences. Similarly, when viewing the largest combined datasets, take note of the number of sequences in the bin of interest. On a typical laptop running RStudio searching a bin of 1000 sequences for the most closely related CDR3 motifs will take about 1 minute, and about 40-60 seconds using the web portal. Note that searches of larger bins will take considerably longer - searching ~7,500 sequences will take 12-15 minutes and may not finish before timing out on the web portal, which will disconnect after 15 minutes of inactivity.

Note that after exploring one table for some time, it is possible to unwittingly have multiple sequences selected but not in view - this will affect the topology-making options. A simple solution is to click on another bin, and then click back on the bin of interest thus refreshing the table. Also note that the topologies require the selected sequences to all have the same CDR3 length - otherwise the calculation will fail with an error: ```Warning: Error in <-: length of 'dimnames' [1] not equal to array extent```. This error is likely to result after clicking on a bin while viewing repertoires displayed across multiple panels - in that case, the table will sometimes include sequences with multiple CDR3 lengths. In this scenario, finding the most closely related CDR3 motifs to a sequence of interest will fail, though manually creating NJ or parsimony topologies of sequences with the same CDR3 length will work.

Other warnings may occur if there are only 1 or 2 closely related CDR3 motifs: ```Warning: Error in [[: subscript out of bounds``` & ```Warning: Error in nj: cannot build an NJ tree with less than 3 observations```. If viewing selected sequences from the (upper) table made via a bounding box, the topology functions will not work and instead show an error: ```Error: argument is of length zero```. If processing very large datasets in the _Import Data_ tab using the web portal, the error ```Failed - server problem``` may appear after clicking the step 2 button. In that case click on step 3, wait a similar amount of time and then see if that file downloads. Additional alternating clicks of step 2 and step 3 buttons (after 60-90 seconds to allow for processing) often succeeds in such cases.

## Citation
If you use AIRRscape, please cite our publication:

**AIRRscape: an interactive tool for exploring B-cell receptor repertoires and antibody responses**

Eric Waltari, Saba Nafees, Krista M. McCutcheon, Joan Wong, John E. Pak

PLOS Computational Biology 2022; doi: https://doi.org/10.1371/journal.pcbi.1010052
