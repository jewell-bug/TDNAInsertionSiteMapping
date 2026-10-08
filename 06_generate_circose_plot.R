library(circlize)


# read data files

karyotype <- read.table("circos_input_files/karyotype.txt", stringsAsFactors = FALSE)
genes     <- read.table("circos_input_files/genes.txt", stringsAsFactors = FALSE)
mutants   <- read.table("circos_input_files/mutants.txt", stringsAsFactors = FALSE)

# Format Karyotype Data
karyo_df <- data.frame(
  chr   = karyotype$V3,
  start = as.numeric(karyotype$V5),
  end   = as.numeric(karyotype$V6)
)
karyo_df <- karyo_df[order(-karyo_df$end), ]

# FILTER: Keep top 37 contigs (>= 260 kb)
karyo_df <- karyo_df[karyo_df$end >= 260000, ]

# Parse Mutants & Genes
genes_df   <- genes[genes$V1 %in% karyo_df$chr, ]
genes_df   <- data.frame(chr = genes_df$V1, start = as.numeric(genes_df$V2), end = as.numeric(genes_df$V3))

mutants_df <- mutants[mutants$V1 %in% karyo_df$chr, ]
mutants_df <- data.frame(
  chr    = mutants_df$V1,
  start  = as.numeric(mutants_df$V2),
  end    = as.numeric(mutants_df$V3),
  color  = sub(".*fill_color=([^,]+).*", "\\1", mutants_df$V4),
  label  = sub(".*label=([^,]+).*", "\\1", mutants_df$V4)
)

color_map <- c("orange" = "#FF7F00", "purple" = "#984EA3", "green" = "#4DAF4A", "red" = "#E41A1C", "blue" = "#377EB8")
mutants_df$r_color <- ifelse(mutants_df$color %in% names(color_map), color_map[mutants_df$color], mutants_df$color)



pdf("Whole_Genome_5_Mutants_Final.pdf", width = 11, height = 11)

circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.01, 0.01))

# Track 0: Chromosome Axes
circos.genomicInitialize(karyo_df, plotType = c("axis"), labels.cex = 0.5)

# Clean Sector Labels (NODE_1 to NODE_37)
circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    sector.name <- CELL_META$sector.index
    display_name <- sub("(_length_.*)", "", sector.name)
    circos.text(CELL_META$xcenter, 0.5, display_name, facing = "inside", niceFacing = TRUE, cex = 0.5, font = 2)
  },
  track.height = 0.05, bg.border = NA
)

# Track 1: Gene Density / Exons
circos.genomicTrack(
  genes_df, ylim = c(0, 1), track.height = 0.08,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "#2B8CBE", border = NA, ytop = 0.9, ybottom = 0.1)
  }
)

# Track 2: Staggered Mutants  (E3, F8, G9 on NODE_2)
circos.genomicTrack(
  mutants_df, ylim = c(0, 1), track.height = 0.18,
  panel.fun = function(region, value, ...) {
    current_contig  <- CELL_META$sector.index
    current_mutants <- mutants_df[mutants_df$chr == current_contig, ]
    n_mutants       <- nrow(current_mutants)
    
    if (n_mutants > 0) {
      for (i in 1:n_mutants) {
        sub_region <- data.frame(start = current_mutants$start[i], end = current_mutants$end[i])
        
        # Vertical stacking
        y_bottom <- (i - 1) / n_mutants
        y_top    <- i / n_mutants
        
        circos.genomicRect(sub_region, col = current_mutants$r_color[i], border = "black", ytop = y_top, ybottom = y_bottom)
        
        # Vertical label staggering
        text_y <- 1.1 + ((i - 1) * 0.35)
        circos.genomicText(sub_region, labels = current_mutants$label[i], y = text_y, cex = 0.75, font = 2, col = current_mutants$r_color[i])
      }
    }
  }
)

title("Aureobasidium pullulans: T-DNA Integration Map", cex.main = 1.2)
legend("bottomright", legend = mutants_df$label, fill = mutants_df$r_color, bty = "n", cex = 0.85, title = "Mutant Lines")

dev.off()
circos.clear()

print("PDF rendered: Whole_Genome_5_Mutants_Final.pdf")




library(circlize)


karyotype <- read.table("circos_input_files/karyotype.txt", stringsAsFactors = FALSE)
genes     <- read.table("circos_input_files/genes.txt", stringsAsFactors = FALSE)
mutants   <- read.table("circos_input_files/mutants.txt", stringsAsFactors = FALSE)

# Format Karyotype Data
karyo_df <- data.frame(
  chr   = karyotype$V3,
  start = as.numeric(karyotype$V5),
  end   = as.numeric(karyotype$V6)
)

# Parse Mutants & Genes
genes_df   <- data.frame(chr = genes$V1, start = as.numeric(genes$V2), end = as.numeric(genes$V3))

mutants_df <- data.frame(
  chr    = mutants$V1,
  start  = as.numeric(mutants$V2),
  end    = as.numeric(mutants$V3),
  color  = sub(".*fill_color=([^,]+).*", "\\1", mutants$V4),
  label  = sub(".*label=([^,]+).*", "\\1", mutants$V4)
)

color_map <- c("orange" = "#FF7F00", "purple" = "#984EA3", "green" = "#4DAF4A", "red" = "#E41A1C", "blue" = "#377EB8")
mutants_df$r_color <- ifelse(mutants_df$color %in% names(color_map), color_map[mutants_df$color], mutants_df$color)


# loop through contigs

target_contigs <- unique(mutants_df$chr)

for (contig_id in target_contigs) {
  
  # Subset data for this specific contig only
  karyo_subset   <- karyo_df[karyo_df$chr == contig_id, ]
  genes_subset   <- genes_df[genes_df$chr == contig_id, ]
  mutants_subset <- mutants_df[mutants_df$chr == contig_id, ]
  
  display_name <- sub("(_length_.*)", "", contig_id)
  base_filename <- paste0("T-DNA_Contig_", display_name)
  
  
  generate_plot <- function() {
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 0, track.margin = c(0.01, 0.01))
    
    circos.genomicInitialize(karyo_subset, plotType = c("axis"), labels.cex = 0.6)
    
    # Track 0: Clean Sector Label
    circos.track(
      ylim = c(0, 1),
      panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, 0.5, display_name, facing = "inside", niceFacing = TRUE, cex = 0.9, font = 2)
      },
      track.height = 0.08, bg.border = NA
    )
    
    # Track 1: Gene Density / Exons
    if (nrow(genes_subset) > 0) {
      circos.genomicTrack(
        genes_subset, ylim = c(0, 1), track.height = 0.12,
        panel.fun = function(region, value, ...) {
          circos.genomicRect(region, value, col = "#2B8CBE", border = NA, ytop = 0.9, ybottom = 0.1)
        }
      )
    }
    
    # Track 2: Mutants Stacking
    circos.genomicTrack(
      mutants_subset, ylim = c(0, 1), track.height = 0.25,
      panel.fun = function(region, value, ...) {
        n_mutants <- nrow(mutants_subset)
        if (n_mutants > 0) {
          for (i in 1:n_mutants) {
            sub_region <- data.frame(start = mutants_subset$start[i], end = mutants_subset$end[i])
            y_bottom <- (i - 1) / n_mutants
            y_top    <- i / n_mutants
            
            circos.genomicRect(sub_region, col = mutants_subset$r_color[i], border = "black", ytop = y_top, ybottom = y_bottom)
            
            text_y <- 1.1 + ((i - 1) * 0.4)
            circos.genomicText(sub_region, labels = mutants_subset$label[i], y = text_y, cex = 0.85, font = 2, col = mutants_subset$r_color[i])
          }
        }
      }
    )
    
    title(paste0("T-DNA Insertion Contig: ", display_name), cex.main = 1.2)
    legend("bottomright", legend = mutants_subset$label, fill = mutants_subset$r_color, bty = "n", cex = 0.85, title = "Mutant Lines")
  }
  
  # Export as PDF
  pdf(paste0(base_filename, ".pdf"), width = 8, height = 8)
  generate_plot()
  dev.off()
  
  # Export as JPG 
  jpeg(paste0(base_filename, ".jpg"), width = 2400, height = 2400, res = 300)
  generate_plot()
  dev.off()
  
  circos.clear()
  print(paste0("Exported: ", base_filename, ".pdf and .jpg"))
