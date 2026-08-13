Nuclear membrane quantification of smFISH signal in C. elegans worms
================
Sam Zavislan-Pullaro and Erin Osborne Nishimura
2026-08-13

- [Step 1: libraries, import, and data
  munging](#step-1-libraries-import-and-data-munging)
  - [Load libraries](#load-libraries)
  - [Import data](#import-data)
  - [Create annotation columns and add
    timepoints](#create-annotation-columns-and-add-timepoints)
  - [Merge all data and apply fixed coordinate
    alignment](#merge-all-data-and-apply-fixed-coordinate-alignment)
- [Step 2: Normalize](#step-2-normalize)
  - [Color palettes](#color-palettes)
- [Step 3: Generate individual normalized line
  scans](#step-3-generate-individual-normalized-line-scans)
  - [Save the individual linescan
    plots](#save-the-individual-linescan-plots)
- [Step 4: Generate mean line scans](#step-4-generate-mean-line-scans)
  - [Save the mean line plots](#save-the-mean-line-plots)
- [Step 5: Calculate log2 fold change
  enrichment](#step-5-calculate-log2-fold-change-enrichment)
  - [Generate merged boxplot, x-axis by
    treatment](#generate-merged-boxplot-x-axis-by-treatment)
  - [Summary metrics](#summary-metrics)
  - [Statistics - Pairwise Wilcoxon](#statistics---pairwise-wilcoxon)
  - [Reps and N-values](#reps-and-n-values)
  - [Plot boxplot](#plot-boxplot)
  - [Save the box plots](#save-the-box-plots)
- [Step 6: Export plots and stats](#step-6-export-plots-and-stats)
- [Session info](#session-info)

# Step 1: libraries, import, and data munging

## Load libraries

- tidyverse
- data.table
- rstatix
- ggrepel
- gridExtra

``` r
library(tidyverse)
library(data.table)
library(rstatix)
library(ggrepel)
library(gridExtra)
```

------------------------------------------------------------------------

## Import data

Note: channel 1 = imb-2 mRNA \| channel 2 = set-3 mRNA

Strain: N2

``` r
# Load the data
# Only need to replace name of file. EX: ("../01_input/<NAME_OF_DATAFILE.txt>",header = FALSE, sep = "\t")

df_1 <- read.table("../01_input/2026-7-20_1510_datafile_for_250806A_set3_imb2_rep1.txt",
                        header = FALSE, sep = "\t")
df_2 <- read.table("../01_input/2026-7-20_1518_datafile_for_250806B_set3_imb2_rep2.txt",
                        header = FALSE, sep = "\t")
df_3 <- read.table("../01_input/2026-7-20_1525_datafile_for_250807_set3_imb2_rep3.txt",
                        header = FALSE, sep = "\t")
df_4 <- read.table("../01_input/2026-8-3_128_datafile_for_260527_N2_imb2_set3.txt",
                        header = FALSE, sep = "\t")

# Check dimensions
# In this case, all data files should have vector length = 121. This is determined by the length of the ROI box in FIJI

dim(df_1)
```

    ## [1]   9 121

``` r
#dim(df_2)
#dim(df_3)
#dim(df_4)

df_1 <- df_1[,1:120]
df_2 <- df_2[,1:120]
df_3 <- df_3[,1:120]
df_4 <- df_4[,1:120]

# The value 119 is actually determined by the dimensions, just subtract 2 from what you get in dimensions

col_times <- paste(rep(1:118), sep = "") 

# Set column names for the dataframes 

colnames(df_1) <- c("file", "channel", col_times)
colnames(df_2) <- c("file", "channel", col_times)
colnames(df_3) <- c("file", "channel", col_times)
colnames(df_4) <- c("file", "channel", col_times)

head(df_1)
```

    ##                             file channel         1          2          3
    ## 1                           file channel     0.000     0.1072     0.2144
    ## 2 250806A_N2_imb2_set3_02_R3D.dv     ch1 10216.989 10238.6611 10240.1514
    ## 3 250806A_N2_imb2_set3_02_R3D.dv     ch2  7901.378  7922.0454  7945.9189
    ## 4 250806A_N2_imb2_set3_06_R3D.dv     ch1  8530.940  8509.3584  8516.5479
    ## 5 250806A_N2_imb2_set3_06_R3D.dv     ch2  9997.788 10118.3477 10241.3223
    ## 6 250806A_N2_imb2_set3_08_R3D.dv     ch1  9188.416  9061.3389  8944.8965
    ##            4          5          6          7          8          9         10
    ## 1     0.3215     0.4287     0.5359     0.6431     0.7503     0.8574     0.9646
    ## 2 10230.8896 10240.9062 10259.0605 10302.1396 10325.9023 10350.3682 10375.5088
    ## 3  7961.6060  7965.7173  7984.8101  8014.6577  8044.3516  8051.6558  8061.3760
    ## 4  8556.0654  8591.8154  8638.3271  8697.7998  8818.8213  9034.0752  9308.5420
    ## 5 10393.2393 10566.8779 10745.1729 10891.3525 10958.1016 10952.8193 10935.8203
    ## 6  8909.6562  8896.1611  8925.8457  8950.6357  8968.6445  8992.1289  9008.0107
    ##           11        12         13         14         15         16         17
    ## 1     1.0718     1.179     1.2862     1.3933     1.5005     1.6077     1.7149
    ## 2 10367.4805 10352.462 10345.5010 10358.2686 10372.8779 10415.4414 10429.0010
    ## 3  8080.8018  8093.791  8124.0186  8140.7026  8181.9009  8206.4648  8217.6074
    ## 4  9628.7217  9936.962 10216.8232 10406.7578 10509.1309 10549.9102 10516.6982
    ## 5 10959.8467 10994.402 11047.3037 11054.4619 11042.9297 11007.3896 10994.8857
    ## 6  9020.2695  9015.352  8964.6475  8944.3057  8889.9961  8845.8193  8824.4121
    ##           18         19         20         21         22        23         24
    ## 1     1.8221     1.9292     2.0364     2.1436     2.2508     2.358     2.4651
    ## 2 10456.6475 10484.7617 10503.0225 10533.2588 10558.4111 10570.197 10594.7920
    ## 3  8241.2041  8260.2031  8262.1572  8291.5762  8321.1543  8352.180  8389.7910
    ## 4 10391.3633 10288.9424 10263.3545 10304.0684 10427.9375 10603.376 10798.2363
    ## 5 10995.2031 10984.8564 10994.5781 10941.0918 10909.2920 10879.887 10875.9531
    ## 6  8797.4268  8831.5830  8887.4756  8962.2793  9054.2383  9119.044  9183.1045
    ##           25         26         27         28        29         30         31
    ## 1     2.5723     2.6795     2.7867     2.8939     3.001     3.1082     3.2154
    ## 2 10645.5625 10671.5010 10715.8691 10792.7627 10862.693 10934.4072 11021.4863
    ## 3  8432.0967  8452.0205  8499.7129  8552.1895  8602.146  8640.9580  8688.0557
    ## 4 10945.4678 11010.7266 11033.7764 11022.4434 10970.410 10971.5488 10991.9365
    ## 5 10874.4736 10844.6484 10803.6455 10774.7686 10754.482 10747.0547 10744.1357
    ## 6  9232.3779  9165.4697  9122.3828  9066.6143  8995.291  8952.9355  8984.8174
    ##           32         33         34         35         36         37         38
    ## 1     3.3226     3.4298     3.5369     3.6441     3.7513     3.8585     3.9657
    ## 2 11148.8857 11278.9492 11395.7686 11544.3760 11709.9062 11821.3408 11946.7236
    ## 3  8749.1338  8792.0361  8832.6553  8895.5664  8960.3125  9011.9580  9064.6416
    ## 4 11004.3125 11046.2246 11096.1172 11103.5410 11080.9824 11048.7725 11042.5898
    ## 5 10751.3496 10772.5566 10776.5078 10762.2598 10760.7334 10761.8613 10754.8262
    ## 6  8983.7051  8960.9512  9022.7441  9194.2832  9371.9824  9465.6494  9460.3662
    ##           39        40         41         42         43         44         45
    ## 1     4.0728     4.180     4.2872     4.3944     4.5016     4.6087     4.7159
    ## 2 12113.4219 12272.608 12459.7432 12620.5264 12825.7607 12998.7090 13201.9365
    ## 3  9110.2695  9163.968  9217.7803  9275.0059  9345.5195  9413.9365  9480.7754
    ## 4 11039.8643 11083.689 11165.9072 11259.9346 11335.6826 11398.7627 11479.8584
    ## 5 10733.3613 10719.160 10699.4521 10716.0479 10724.7402 10746.8457 10768.3916
    ## 6  9401.9551  9381.799  9441.1133  9524.0137  9591.8867  9685.1982  9732.7480
    ##           46         47         48         49         50        51         52
    ## 1     4.8231     4.9303     5.0375     5.1446     5.2518     5.359     5.4662
    ## 2 13474.3477 13720.1396 13984.1650 14325.3574 14682.3789 15013.988 15331.2256
    ## 3  9544.0537  9620.4482  9698.4805  9798.6162  9892.4004 10010.659 10108.5400
    ## 4 11607.7490 11715.3887 11807.3086 11932.8037 12098.2842 12314.970 12563.0146
    ## 5 10806.8936 10833.0049 10859.4854 10906.3906 10951.8418 11006.058 11073.6543
    ## 6  9777.2207  9897.3740 10021.8398 10077.7939 10061.4678 10080.348 10119.6641
    ##           53         54         55         56         57         58         59
    ## 1     5.5734     5.6805     5.7877     5.8949     6.0021     6.1093     6.2164
    ## 2 15634.1221 15959.5186 16146.1406 16308.8984 16393.6973 16336.5840 16138.7529
    ## 3 10195.2891 10304.9580 10396.1748 10501.1924 10569.1729 10634.9854 10686.4736
    ## 4 12822.7285 13043.2080 13225.5996 13420.8203 13555.6602 13529.9824 13404.7686
    ## 5 11126.5430 11167.6533 11222.2148 11276.3809 11360.6787 11452.5127 11552.2021
    ## 6 10299.6592 10532.9502 10733.9287 10975.8262 11146.6250 11191.5186 11208.4805
    ##           60         61        62         63         64         65         66
    ## 1     6.3236     6.4308     6.538     6.6452     6.7523     6.8595     6.9667
    ## 2 16022.1377 15950.0352 15643.948 15063.8525 14360.0664 13627.2178 13006.5469
    ## 3 10733.4229 10786.1611 10805.247 10827.2646 10877.6357 10925.5713 10967.1162
    ## 4 13233.4629 12969.3516 12604.298 12162.1904 11621.8340 11110.3857 10752.4385
    ## 5 11647.2773 11785.4990 11935.345 12065.4854 12203.9238 12196.6309 12137.3018
    ## 6 11085.4053 10922.6357 10868.836 10953.4551 11090.2559 11106.1484 10956.2598
    ##           67         68         69         70         71         72        73
    ## 1     7.0739     7.1811     7.2882     7.3954     7.5026     7.6098     7.717
    ## 2 12587.2354 12295.6016 12031.2373 11812.3662 11617.8086 11416.7998 11244.055
    ## 3 10950.9053 10940.9258 10908.8457 10896.0107 10902.0771 10940.7959 11009.748
    ## 4 10412.9414 10178.9102  9963.5293  9746.5801  9576.5391  9440.5742  9328.665
    ## 5 12027.9463 11835.9922 11585.3750 11323.5879 11122.6689 10953.1895 10863.327
    ## 6 10583.7324 10153.0430  9799.7285  9565.0586  9368.8604  9105.1006  8775.599
    ##           74         75         76         77         78        79         80
    ## 1     7.8241     7.9313     8.0385     8.1457     8.2529     8.360     8.4672
    ## 2 11060.0850 10879.9014 10691.1631 10450.9873 10295.1846 10235.323 10190.2617
    ## 3 11126.0859 11264.8740 11414.4854 11571.0742 11711.5117 11736.849 11619.8076
    ## 4  9176.9883  8989.6426  8814.1377  8695.1982  8647.7158  8612.657  8561.7344
    ## 5 10769.6436 10688.8047 10621.5996 10531.0684 10476.2559 10454.329 10441.0957
    ## 6  8471.9658  8218.8193  8073.6655  7996.8223  7961.0156  7894.716  7771.3438
    ##           81         82         83         84         85         86         87
    ## 1     8.5744     8.6816     8.7888     8.8959     9.0031     9.1103     9.2175
    ## 2 10196.6387 10287.7852 10329.8818 10268.9717 10146.3555 10046.2178  9911.6426
    ## 3 11588.9854 11617.7627 11631.3779 11743.1338 11892.9922 12030.2217 11907.6240
    ## 4  8493.0283  8413.6133  8370.3506  8349.2344  8316.0107  8295.6582  8305.8887
    ## 5 10420.3789 10380.0967 10367.1172 10337.8711 10304.1934 10260.8223 10236.8877
    ## 6  7632.9663  7494.6646  7390.7690  7325.2275  7271.2402  7216.7485  7187.0894
    ##           88         89        90         91         92         93         94
    ## 1     9.3247     9.4318     9.539     9.6462     9.7534     9.8606     9.9677
    ## 2  9793.9180  9727.2627  9748.450  9826.5918  9909.4902  9939.4102  9923.5049
    ## 3 11613.6943 11367.4932 11230.640 11175.0742 11139.7617 11134.2559 11139.2832
    ## 4  8342.0547  8393.2832  8508.171  8615.2891  8700.4922  8714.3643  8686.7266
    ## 5 10233.7822 10231.4102 10255.409 10310.9854 10445.4932 10629.7832 10853.0986
    ## 6  7149.7002  7130.4536  7117.369  7120.3164  7129.5488  7157.4829  7200.2402
    ##           95         96         97         98         99        100       101
    ## 1    10.0749    10.1821    10.2893    10.3965    10.5036    10.6108    10.718
    ## 2  9951.8594 10037.9404 10057.1230  9966.9580  9839.5811  9711.5957  9546.867
    ## 3 11148.1465 11231.7598 11350.6465 11477.5244 11561.0322 11573.8076 11508.905
    ## 4  8620.9609  8551.4365  8507.5986  8494.1338  8464.3320  8408.5107  8382.686
    ## 5 10989.4941 10849.6445 10545.6396 10224.7393  9940.1836  9706.0098  9583.413
    ## 6  7265.0527  7392.6108  7523.2334  7574.7158  7539.8276  7456.7388  7380.733
    ##          102        103        104        105        106        107        108
    ## 1    10.8252    10.9324    11.0395    11.1467    11.2539    11.3611    11.4683
    ## 2  9423.8262  9326.8418  9280.8096  9199.0762  9088.0029  8982.5225  8882.0898
    ## 3 11328.1279 11081.1631 10996.1064 10973.5693 10959.8770 10941.9092 11040.6357
    ## 4  8374.8398  8383.9199  8357.9922  8301.5459  8223.8086  8194.5547  8199.9648
    ## 5  9482.5430  9413.3477  9346.4434  9292.1875  9246.7256  9211.5264  9176.9229
    ## 6  7322.6333  7340.1719  7408.0854  7490.5386  7517.7100  7476.2749  7428.7861
    ##          109        110        111       112        113        114        115
    ## 1    11.5754    11.6826    11.7898    11.897    12.0042    12.1113    12.2185
    ## 2  8729.3926  8624.7686  8568.8174  8522.869  8498.4365  8455.6299  8422.5938
    ## 3 11150.8818 11137.3564 10999.9111 10906.081 10798.2578 10611.4336 10437.1104
    ## 4  8202.4346  8257.7070  8309.3320  8376.038  8385.7676  8285.0869  8179.2510
    ## 5  9135.9688  9104.2129  9088.2051  9112.773  9135.0127  9126.7148  9119.6094
    ## 6  7382.4927  7367.7539  7383.3335  7395.102  7434.9209  7495.5205  7526.0327
    ##          116        117        118
    ## 1    12.3257    12.4329    12.5401
    ## 2  8441.0264  8496.3994  8595.7783
    ## 3 10347.9746 10298.0361 10245.6768
    ## 4  8076.9590  7989.9775  7903.5869
    ## 5  9042.2812  8960.5684  8913.6914
    ## 6  7569.5186  7567.1694  7541.3110

------------------------------------------------------------------------

## Create annotation columns and add timepoints

Timepoints are treated as positions along line scan

``` r
# Helper function: pivot longer and annotate 
annotate_df <- function(df) {
  df[2:nrow(df), 1:120] %>% # 120 comes from dimensions
    separate_wider_delim(file, delim = "_",
                         names = c("date", "strain", "mRNA1", "mRNA2", "embryoID", NA)) %>%
    pivot_longer(cols = `1`:`118`, names_to = "xpoint", values_to = "intensity") %>%
    mutate(unique_id = paste(date, strain, mRNA1, mRNA2, embryoID, channel, sep = "_"))
}

# Helper function: add timepoints from the header row of any loaded table
addInTimepoints <- function(tibble, ref_table) {
  timepoints <- as.vector(unlist(ref_table[1, 3:120]))
  repnumb    <- nrow(tibble)/118
  tibble::add_column(tibble, timepoints = rep(timepoints, repnumb), .after = "xpoint")
}

# Add timepoints as proxy for pixel position along line scan 
df_1_pt <- addInTimepoints(annotate_df(df_1), df_1) # Use first dataframe as reference 
df_2_pt <- addInTimepoints(annotate_df(df_2), df_1)
df_3_pt <- addInTimepoints(annotate_df(df_3), df_1)
df_4_pt <- addInTimepoints(annotate_df(df_4), df_1)

# Sanity check 
head(df_4_pt)
```

    ## # A tibble: 6 × 10
    ##   date   strain mRNA1 mRNA2 embryoID channel xpoint timepoints intensity
    ##   <chr>  <chr>  <chr> <chr> <chr>    <chr>   <chr>       <dbl>     <dbl>
    ## 1 260527 N2     imb2  set3  image04  ch1     1           0         4031.
    ## 2 260527 N2     imb2  set3  image04  ch1     2           0.107     4028.
    ## 3 260527 N2     imb2  set3  image04  ch1     3           0.214     4001.
    ## 4 260527 N2     imb2  set3  image04  ch1     4           0.322     3981.
    ## 5 260527 N2     imb2  set3  image04  ch1     5           0.429     3960.
    ## 6 260527 N2     imb2  set3  image04  ch1     6           0.536     3956.
    ## # ℹ 1 more variable: unique_id <chr>

``` r
#df_4_pt[110:120,]
```

------------------------------------------------------------------------

## Merge all data and apply fixed coordinate alignment

``` r
# Concatenate the dataframes 
all_data <- rbind(
  df_1_pt, df_2_pt, df_3_pt, df_4_pt
)

center <- ceiling(118/2) # Using dims, find center point 
center
```

    ## [1] 59

``` r
dt <- as.data.table(all_data)
dt$xpoint <- as.integer(dt$xpoint)
dt[, aligned_row := xpoint - center]

total_align_long <- dt %>%
  select(strain, aligned_row, unique_id, intensity)

cat("Total rows:", nrow(total_align_long), "\n")
```

    ## Total rows: 4248

``` r
cat("Unique embryo IDs:", n_distinct(total_align_long$unique_id), "\n")
```

    ## Unique embryo IDs: 36

``` r
# head(total_align_long)
```

------------------------------------------------------------------------

# Step 2: Normalize

Each embryo is normalized to its own mean intensity (per embryo per
channel).

``` r
norm_total <- total_align_long %>%
  separate_wider_delim(unique_id, delim = "_",
                       names = c("date", NA, "mRNA1", "mRNA2", "embryoID", "channel")) %>%
  group_by(date, embryoID, channel, mRNA1, mRNA2) %>%
  mutate(normalized_intensity = intensity / mean(intensity, na.rm = TRUE)) %>%
  ungroup()

treatment_order <- c("ch2", "ch1")
norm_total$channel <- factor(norm_total$channel, levels = treatment_order)

# ch1 = imb-2, ch2 = set-3
norm_total <- norm_total %>%
  mutate(treatment_type = case_when(
    channel == "ch1" ~ "imb-2",
    channel == "ch2" ~ "set-3",
    TRUE ~ NA_character_
  ))

table(norm_total$channel, norm_total$treatment_type)
```

    ##      
    ##       imb-2 set-3
    ##   ch2     0  2124
    ##   ch1  2124     0

------------------------------------------------------------------------

## Color palettes

``` r
treatment_colors <- c(
  "ch1" = "#2166AC",   # imb-2 — dark blue
  "ch2" = "#4D4D4D"    # set-3 — dark gray
)
type_colors <- treatment_colors

channel_labels <- c("ch1" = "imb-2 mRNA", "ch2" = "set-3 mRNA")
```

------------------------------------------------------------------------

# Step 3: Generate individual normalized line scans

All embryo-level line scans, faceted by channel and treatment.

``` r
p_linescan_individual <- ggplot(norm_total,
       aes(x = aligned_row, y = normalized_intensity,
           group = interaction(date, embryoID),
           color = channel)) +
  geom_line(alpha = 0.6, linewidth = 0.3) +
  scale_color_manual(values = treatment_colors) +
  facet_wrap(~ channel, labeller = labeller(channel = channel_labels)) +
  labs(
    title = "Individual normalized line scans — all treatments",
    x     = "Position along line scan",
    y     = "Normalized intensity",
    color = "Treatment"
  ) +
  guides(color = "none") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))
p_linescan_individual
```

![](260813_imb2_set3_analysis_SZ_EON_files/figure-gfm/linescans_individual-1.png)<!-- -->

------------------------------------------------------------------------

## Save the individual linescan plots

``` r
# Create an output directory
output_dir <- "../03_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save a filename
today <- format(Sys.Date(), "%y%m%d")
filename <- paste(output_dir, "/", today, "_linescan_individual.pdf", sep = "")

# Save the plot
pdf(filename, width = 7, height = 6)
p_linescan_individual
dev.off()
```

    ## quartz_off_screen 
    ##                 2

# Step 4: Generate mean line scans

``` r
# Calculate the mean signal of all the linscan plots and calculate the standard deviation
mean_signal <- norm_total %>%
  group_by(aligned_row, channel) %>%
  summarize(
    mean_signal = mean(normalized_intensity, na.rm = TRUE),
    sd_signal   = sd(normalized_intensity,   na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(
    ymeanmax = mean_signal + sd_signal,
    ymeanmin = mean_signal - sd_signal
  )

# Plot the means and sd as linescans and ribbons
p_linescan_mean <- ggplot(mean_signal,
       aes(x = aligned_row, y = mean_signal,
           group = channel)) +
  geom_ribbon(aes(ymin = ymeanmin, ymax = ymeanmax,
                  fill = channel), alpha = 0.2, color = NA) +
  geom_line(aes(color = channel), linewidth = 0.8) +
  scale_color_manual(values = treatment_colors, labels = channel_labels) +
  scale_fill_manual(values  = treatment_colors, labels = channel_labels) +
  labs(
    title = "Mean +- SD line scans, per channel",
    x     = "Position",
    y     = "Mean normalized intensity"
  ) +
  theme_bw(base_size = 11)

 p_linescan_mean
```

![](260813_imb2_set3_analysis_SZ_EON_files/figure-gfm/linescans_mean_overlaid-1.png)<!-- -->

------------------------------------------------------------------------

## Save the mean line plots

``` r
# Save a filename
today <- format(Sys.Date(), "%y%m%d")
filename <- paste(output_dir, "/", today, "_linescan_mean.pdf", sep = "")

# Save the plot
pdf(filename, width = 7, height = 6)
p_linescan_mean
dev.off()
```

    ## quartz_off_screen 
    ##                 2

------------------------------------------------------------------------

# Step 5: Calculate log2 fold change enrichment

Ratio of normalized intensity at position 1 versus the mean of positions
-50 and +50. Higher values indicate more enriched localization signal
(consistent with membrane localization).

``` r
peaks_and_valleys <- norm_total %>%
  filter(aligned_row %in% c(-50, 1, 50))

nest_pandv <- peaks_and_valleys %>%
  group_by(date, mRNA1, mRNA2, embryoID, channel) %>%
  nest()

my_calc2 <- function(df) {
  df$normalized_intensity[df$aligned_row == 1] /
    mean(df$normalized_intensity[df$aligned_row %in% c(-50, 50)])
}

foldChange_calc <- nest_pandv %>%
  mutate(fc_enrich = map_dbl(data, my_calc2)) %>%
  mutate(
    log2_fc        = log(fc_enrich, base = 2),
    treatment_type = case_when(
      channel == "ch1" ~ "imb-2",
      channel == "ch2" ~ "set-3",
      TRUE ~ NA_character_
    ),
    channel        = factor(channel, levels = treatment_order)
  )

summary(foldChange_calc$log2_fc)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 0.06885 0.15457 0.20954 0.29430 0.42134 0.74802

## Generate merged boxplot, x-axis by treatment

``` r
foldChange_calc <- foldChange_calc %>%
  mutate(treatment_type = factor(treatment_type, levels = c("set-3", "imb-2")))
```

------------------------------------------------------------------------

## Summary metrics

``` r
summary_stats_table <- foldChange_calc %>%
  group_by(treatment_type, channel) %>%
  summarise(
    n             = n(),
    min_log2fc    = round(min(log2_fc,    na.rm = TRUE), 3),
    median_log2fc = round(median(log2_fc, na.rm = TRUE), 3),
    mean_log2fc   = round(mean(log2_fc,   na.rm = TRUE), 3),
    max_log2fc    = round(max(log2_fc,    na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(channel, treatment_type)

print(summary_stats_table, n = Inf)
```

    ## # A tibble: 2 × 7
    ##   treatment_type channel     n min_log2fc median_log2fc mean_log2fc max_log2fc
    ##   <fct>          <fct>   <int>      <dbl>         <dbl>       <dbl>      <dbl>
    ## 1 set-3          ch2        18      0.069         0.159       0.156      0.242
    ## 2 imb-2          ch1        18      0.076         0.427       0.433      0.748

------------------------------------------------------------------------

## Statistics - Pairwise Wilcoxon

``` r
# Wilcoxon test comparing set-3 vs imb-2 enrichment
pairwise_stats <- foldChange_calc %>%
  ungroup() %>%
  wilcox_test(log2_fc ~ treatment_type) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

pairwise_stats
```

    ## # A tibble: 1 × 9
    ##   .y.     group1 group2    n1    n2 statistic          p      p.adj p.adj.signif
    ##   <chr>   <chr>  <chr>  <int> <int>     <dbl>      <dbl>      <dbl> <chr>       
    ## 1 log2_fc set-3  imb-2     18    18        27 0.00000317 0.00000317 ****

``` r
# pairwise_stats$p
```

------------------------------------------------------------------------

## Reps and N-values

How many replicates and samples are represented in the dataset?

``` r
# select the metadata
metadata <- norm_total %>% 
  select("strain", "date", "embryoID")

# pick only the unique ones
metadata <- unique(metadata)

dim(metadata)
```

    ## [1] 18  3

``` r
str(metadata)
```

    ## tibble [18 × 3] (S3: tbl_df/tbl/data.frame)
    ##  $ strain  : chr [1:18] "N2" "N2" "N2" "N2" ...
    ##  $ date    : chr [1:18] "250806A" "250806A" "250806A" "250806A" ...
    ##  $ embryoID: chr [1:18] "02" "06" "08" "09" ...

``` r
head(metadata)
```

    ## # A tibble: 6 × 3
    ##   strain date    embryoID
    ##   <chr>  <chr>   <chr>   
    ## 1 N2     250806A 02      
    ## 2 N2     250806A 06      
    ## 3 N2     250806A 08      
    ## 4 N2     250806A 09      
    ## 5 N2     250806B 01      
    ## 6 N2     250806B 04

``` r
# tabulate the number of n-values. Divide by 2 because both set-3 and erm-1 are both counted as a data point
rep_and_n <- table(metadata$date)
rep_and_n
```

    ## 
    ## 250806A 250806B  250807  260527 
    ##       4       4       2       8

## Plot boxplot

``` r
p_merged <- ggplot(foldChange_calc,
                   aes(x = treatment_type,
                       y = log2_fc,
                       fill = treatment_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = treatment_type),
              position = position_jitterdodge(jitter.width = 0.25),
              size = 1.5, alpha = 0.7) +
  scale_fill_manual(
    values = c("imb-2" = "#2166AC", "set-3" = "#4D4D4D"),
    guide  = "none"
  ) +
  scale_color_manual(
    values = c("imb-2" = "black", "set-3" = "black"),
    guide  = "none"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  labs(
    title = "log2 enrichment at position 0 vs ±50",
    x     = "mRNA",
    y     = "log2 (intensity at 0 / mean intensity at ±50)",
    fill  = NULL
  ) +
  theme_classic(base_size = 12) +
  annotate("text", x = 1.5, y = 1.25, label = paste("p = ", pairwise_stats$p, sep = "")) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1, face = "italic"),
    strip.text      = element_text(face = "bold"),
    legend.position = "none"
  )

p_merged
```

![](260813_imb2_set3_analysis_SZ_EON_files/figure-gfm/plot-1.png)<!-- -->

------------------------------------------------------------------------

## Save the box plots

``` r
# Save a filename
today <- format(Sys.Date(), "%y%m%d")
filename <- paste(output_dir, "/", today, "_boxplot_merged.pdf", sep = "")

# Save the plot
pdf(filename, width = 7, height = 6)
p_merged
dev.off()
```

    ## quartz_off_screen 
    ##                 2

------------------------------------------------------------------------

# Step 6: Export plots and stats

``` r
# Save some files as .csvs
write_csv(dt, file.path(output_dir,paste0(today, "_aligned_raw_data.csv")))
write_csv(summary_stats_table, file.path(output_dir, paste0(today, "_summary_stats.csv")))
write_csv(pairwise_stats, file.path(output_dir, paste0(today, "_pairwise_wilcoxon_stats.csv")))


# Save the rep and n nubmers
today <- format(Sys.Date(), "%y%m%d")
filename <- paste(output_dir, "/", today, "_rep_and_n.txt", sep = "")
write.table(rep_and_n, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
```

------------------------------------------------------------------------

# Session info

``` r
sessionInfo()
```

    ## R version 4.6.0 (2026-04-24)
    ## Platform: aarch64-apple-darwin23
    ## Running under: macOS Tahoe 26.5.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Denver
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3.1   ggrepel_0.9.8     rstatix_1.1.0     data.table_1.18.4
    ##  [5] lubridate_1.9.5   forcats_1.0.1     stringr_1.6.0     dplyr_1.2.1      
    ##  [9] purrr_1.2.2       readr_2.2.0       tidyr_1.3.2       tibble_3.3.1     
    ## [13] ggplot2_4.0.3     tidyverse_2.0.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.6         generics_0.1.4     stringi_1.8.9      hms_1.1.4         
    ##  [5] digest_0.6.39      magrittr_2.0.5     evaluate_1.0.5     grid_4.6.0        
    ##  [9] timechange_0.4.0   RColorBrewer_1.1-3 fastmap_1.2.0      backports_1.5.1   
    ## [13] Formula_1.2-6      scales_1.4.0       abind_1.4-8        cli_3.6.6         
    ## [17] crayon_1.5.3       rlang_1.3.0        bit64_4.8.2        withr_3.0.3       
    ## [21] yaml_2.3.12        otel_0.2.0         parallel_4.6.0     tools_4.6.0       
    ## [25] tzdb_0.5.0         broom_1.0.13       vctrs_0.7.3        R6_2.6.1          
    ## [29] lifecycle_1.0.5    bit_4.6.0          car_3.1-5          vroom_1.7.1       
    ## [33] pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6       glue_1.8.1        
    ## [37] Rcpp_1.1.2         xfun_0.60          tidyselect_1.2.1   rstudioapi_0.19.0 
    ## [41] knitr_1.51         farver_2.1.2       htmltools_0.5.9    labeling_0.4.3    
    ## [45] rmarkdown_2.31     carData_3.0-6      compiler_4.6.0     S7_0.2.2
