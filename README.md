# TermPick
Term-Seq analysis of S. aureus Term-Seq data, with Jai Tree and Sylvania Wu. 
  
The latest scripts are in the folder Source/Call_peaks_v2


## Scripts
The scripts, in order of operations:

1. read_gr_files.Rmd: Read the bedgraph like files (.gr) from different replicates and experimental conditions and collate them into a single file. 
* Each .gr files have the following columns:
  + Column 1: genome sequence ID
  + Column 2: Position
  + Column 3: raw read depth count

2. pick_peak_step.Rmd: Read in the "log_cpm_wider_tbl.tab" file and perform peak picking using the "detect.spikes" function from the "peakPick" library. I have tried to reformat the "detect.spikes" from the "peakPick" R library to run multicores properly, but this failed to run properly, so can only run on one cores. It is still fast enough as long as we process them on HPC. 


Please read description of the algorithm from the peakPick vignette for more details (https://cran.rstudio.com/web/packages/peakPick/vignettes/peakPick-vignette.html). This is my interpretation of the algorithm here:

The detect.spike algorithms look at the window before and after the current position and calculate the standard deviation (SD) of the values within the window. The current position is identified as a peak if the current value is above X number of standard deviastions (e.g. include as peak if the current value is four times the SD). The algorithm is iterative in that peaks that has been identified will be excluded from the next round of analysis. Peaks are identified in subsequent iterations until the total number of peaks identified is stable (and not significantly change).



3. process_peak_pick.Rmd: Read in the files generated by pick_peak_step.Rmd. This groups the peaks together if they are within X nucleotides of each other (default window size of 5). This results in a list of group of sites which were further analysed to extract summary statistics, as described below.
  * This script perform the following calculations:
  + Perform Venn Diagram analysis by comparing the termination sites found by each experimental conditions and the overlaps among the termination sites found by multiple experimental conditions. It then counts the number of termination sites in each section of the Venn diagram. The files are in the directory Results/Sites/UpSet_Counts directory.
  + Calculate the ratio of the number of group of sites common to all experimental conditions divided by the total number of all group of sites. The files are in the directory Results/Sites/UpSet_Ratios directory.
  + Calculate the above numbers for different combinations of grouping distance window (default of 5 nucleotides) and counts per millions (CPM) cutoff value above which peaks are kept for further analysis. (Automatically covert to log2(CPM) as we now use this for analysis.)   
  
4. analyse_param_scan.Rmd: Plot the numbers calculated from the script "process_peak_pick.Rmd" to determine the best combination of parameters. The recommended grouping window size is set to 5 and log2(CPM) cutoff set to 0 as this gives the best number of peaks. 


5. Once the best peak picking window size (window_length parameter) and standard deviation (std_dev_cutoff parameter) the needs to be updated in the "process_peak_pick.Rmd" script. The flag "get_final_results" parameter also needs to be set to TRUE in the "process_peak_pick.Rmd" script and then run again. For example:
  - window_length <- 8;  # Window length for picking (on either side of the current position).
  - std_dev_cutoff <-  7; # Multiple of standard deviation cutoff for peak picking
  - cpm_cutoff_list <- c(1); # CPM threshold above which a peak is kept for further analysis. All peaks with CPM below is values are discarded.
  - distance_window_list <- c( 5); # Consecutive peaks are put together into the same group if they are within this distance apart.
  

This will then produce three output files: 
* Results/Sites/Log_Selected_Parameters/all_positions_and_cpm_win8_sd4_cpm0.tab: 
  + Column 1: position: genomic position 
  + Column 2: Sample: Sample name
  + Column 3: Is_Spike: Is it a peak, 1 = is a peak, 0 = not a peak
  + Column 4: condition (e.g. cont)	
  + Column 5: replicate: replicate number 
  + Column 6: strand: '-' for negative strand or '+' for positive strand  
  + Column 7: depth: Raw read count
  + Column 8: norm_depth: In log2(CPM)
* Results/Sites/Log_Selected_Parameters/condition_vs_position_clusters_win8_sd4_cpm0_gw5.tab
  + Column 1: group: The group ID for a group of peaks	
  + Column 2: cont:	Is there a peak in this group in the control 
 
* Results/Sites/Log_Selected_Parameters/groups_for_all_peaks_win8_sd4_cpm0_gw5.tab
 + Column 1: strand: '-' for negative strand or '+' for positive strand  
 + Column 2: position: genomic position 
 + Column 3: group: The group ID for a group of peaks
 + Example: 
 
|strand	|position	|group |
|------|-----|---------|
|-	|2779	|1 |
|-	|2782	|1 |
|-	|2784	|1 |
|-	|2785	|1 |
|-	|2822	|2 |
|-	|3279	|3 |
|-	|3280	|3 |
  
6. merge_master_table.Rmd: Merge all the result files generated from "process_peak_pick.Rmd" and from "read_gr_files.Rmd" to create the final output table.
This script will find the position with the highest read depth across all experimental conditions in each group and use this is the final position. 

  
## README for output master table:
* Column 1: group     : Group ID for the peaks. Peaks are grouped if they are within 5 nucleotides of each other.
* Column 2: position  : Genomic position of the peak
* Column 3: strand    : Indicates wether the peak is on the positive (+) or negative (-) strand.
* Columns 5 to 7: The following colums show whether this group of peaks was identified in each of the experimental condition. A number one (1) means the peak has been identified in the condition, and zero (0) otherwise.
  + cont

* Columns 8 to 23: The following columns (depth_*, norm_depth_*) provides information for each individual peak without any grouping. These peaks were identified using the detect.spikes() function from the peakPick R library. For each nucleotide location, the standard deviation of the specific position was calculated from the log CPM values of all positions within a window on either side of the specified position. A peak is called if the standard deviation is above the specified threshold. The algortihm is iterative in that it will exclude the log CPM values of identified peaks, then repeat the peak picking process repeatedly, until the current total number of detected peaks does not change significantly.
  + depth_* : Columns with the prefix depth_*. The number at the end of the column name indicates the replicate number. This column records the original read counts corresponding to that position and strand. A value of NA indicates that a peak has not been identified at the corresponding position and strand. This may be NA despite having read counts at the position and strand.
  + norm_depth* : Columns with the prefix norm_depth_*. The number at the end of the column name indicates the replicate numbe. This column records the read deapths converted to log counts per million (log CPM) for the corresponding sample. Peaks must have a CPM of at least 1. A value of NA indicates that a peak has not been identified at the corresponding position and strand. This may be NA despite having read counts at the position and strand.


