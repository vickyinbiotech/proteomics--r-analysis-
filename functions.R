# ============================================
# Custom R Functions
# Author: Victoria Djana
# Description: Utility functions for
# proteomics data analysis
# ============================================

# Check and install missing packages
check_packages <- function(packages) {
missing <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(missing) > 0) {
message("Installing missing packages: ", paste(missing, collapse = ", "))
install.packages(missing)
}
}

# Summary statistics per group
group_summary <- function(data, groups) {
result <- data.frame(
protein = rownames(data),
mean_ctrl = rowMeans(data[, groups == "CTRL"]),
mean_treat = rowMeans(data[, groups == "TREAT"]),
sd_ctrl = apply(data[, groups == "CTRL"], 1, sd),
sd_treat = apply(data[, groups == "TREAT"], 1, sd)
)
result$fold​​​​​​​​​​​​​​​​
