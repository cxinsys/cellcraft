#!/usr/bin/env Rscript

# Load required packages
library(ggplot2)
library(plotly)
library(jsonlite)

# Get file paths from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
pseudotime_path <- args[1]  # pseudotime.csv
cell_select_path <- args[2]  # cellSelect.csv
output_file_path <- args[3]  # pseudotime_density.json

# Load data
pseudotime <- read.csv(pseudotime_path, header = FALSE, col.names = c("Pseudotime"))
cell_select <- read.csv(cell_select_path, header = FALSE, col.names = c("CellSelect"))

# Merge the two data frames
df <- cbind(pseudotime, cell_select)

# Convert CellSelect values to categorical labels
df$Status <- factor(ifelse(df$CellSelect == 1, "Selected", "Unselected"))

# Create a multiple density plot
p <- ggplot(df, aes(x = Pseudotime, fill = Status, color = Status)) +
  geom_density(alpha = 0.5) +  # Density plot
  scale_fill_manual(values = c("Selected" = "#1f77b4", "Unselected" = "#ff7f0e")) +
  scale_color_manual(values = c("Selected" = "#1f77b4", "Unselected" = "#ff7f0e")) +
  ggtitle("Pseudotime Distribution - Density Plot") +
  labs(x = "Pseudotime", y = "Density") +
  theme_minimal()

# Convert ggplot2 object to plotly
p_plotly <- ggplotly(p)

# Safely remove 'frame' attributes
if (!is.null(p_plotly$x$data)) {
  for (i in seq_along(p_plotly$x$data)) {
    p_plotly$x$data[[i]] <- modifyList(p_plotly$x$data[[i]], list(frame = NULL))  # remove frame
  }
}

if ("frame" %in% names(p_plotly$x)) {
  p_plotly$x <- modifyList(p_plotly$x, list(frame = NULL))  # remove frame from layout
}

# Remove 'linetype' property from layout shapes
if (!is.null(p_plotly$x$layout$shapes)) {
  for (i in seq_along(p_plotly$x$layout$shapes)) {
    if ("line" %in% names(p_plotly$x$layout$shapes[[i]])) {
      p_plotly$x$layout$shapes[[i]]$line$linetype <- NULL  # remove linetype
    }
  }
}

# Force all font attributes to use "Arial"
p_plotly <- layout(p_plotly, 
                   font = list(family = "Arial"), 
                   legend = list(
                     font = list(family = "Arial"),
                     title = list(font = list(family = "Arial"))
                   ),
                   title = list(font = list(family = "Arial")),
                   xaxis = list(title = list(font = list(family = "Arial")), tickfont = list(family = "Arial")),
                   yaxis = list(title = list(font = list(family = "Arial")), tickfont = list(family = "Arial"))
)

# Clean up frame and linetype properties in the final JSON
clean_plotly_json <- function(json_content) {
  json_content <- gsub('"frame":\\s*null,?', "", json_content)  # remove frame property
  json_content <- gsub('"linetype":\\s*null,?', "", json_content)  # remove linetype property
  json_content <- gsub(',\\s*}', "}", json_content)  # remove trailing commas
  json_content <- gsub('"family":\\s*""', '"family":"Arial"', json_content)  # replace empty font family with Arial
  return(json_content)
}

# Convert plot to JSON using plotly_json()
json_data <- plotly_json(p_plotly, pretty = TRUE)

# Apply cleanup for frame and font attributes
json_data <- clean_plotly_json(json_data)

# Save as JSON file
write(json_data, output_file_path)

cat(paste("Pseudotime density plot saved to", output_file_path, "\n"))
