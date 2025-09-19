# DE + Expression Viewer

An interactive **Shiny app** for exploring differential expression (DE) results and visualizing expression patterns across samples.  
The app allows you to upload metadata, expression matrices, and DESeq2 results to dynamically generate plots and tables for your genes of interest.  

---

## Features
- **Upload support**:
  - Metadata (CSV/TXT/TSV/XLSX) with a `sample` column (IDs).
  - Expression matrix (CSV/XLSX) with a `gene` column (rows = genes; other columns = samples).
  - Differential expression (DE) results:  
    - RDS file (DESeqResults or `data.frame` with `rownames=genes`), OR  
    - CSV/XLSX table with a `gene` column (gene symbols or IDs).  

- **Visualization**:
  - Expression plots across groups defined in the metadata.  
  - Boxplots with aligned jitter points.  
  - Flexible X-axis and facet selection from any metadata column(s).  
  - Customizable main and axis titles.  

- **DE Table Viewer**:
  - Automatically maps common DESeq2 column names.  
  - Allows searching and filtering for specific genes of interest.  

---

## Example Files
- `example_metadata.csv` → metadata with sample information.  
- `example_counts.csv` → expression matrix (genes × samples).  
- `example_deseq2_results.xlsx` → differential expression results table.  

---

## Installation & Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/de-expression-viewer.git
   cd de-expression-viewer
   ```

2. Install required R packages:
   ```r
   install.packages(c("shiny", "tidyverse", "readxl", "data.table", "DT"))
   ```

3. Run the app:
   ```r
   shiny::runApp("app.R")
   ```

---

## Requirements
- R ≥ 4.0  
- Shiny ≥ 1.7    

---

## Contact
Questions, feature requests, or bug reports:  
- Email: [aysevilpektas@gmail.com](mailto:aysevilpektas@gmail.com)  
- GitHub Issues: [de-expression-viewer](https://github.com/aysevllpkts/de-expression-viewer/issues)

