"""
# **Genome Visualization Tool**

## **Overview**
This tool visualizes genomic data by combining transcript alignment, read coverage, and gene annotation information. It creates a comprehensive figure containing three panels: 
1. Transcript alignment from PSL files.
2. Read coverage histogram.
3. Gene annotation from GTF files.

---

## **Features**
1. **Transcript Alignment Visualization**:
   - Parses PSL files to extract alignment information.
   - Visualizes alignments as horizontal lines with block structures for exons.

2. **Read Coverage Histogram**:
   - Computes nucleotide-level coverage from PSL files.
   - Displays coverage data as a histogram.

3. **Gene Annotation Visualization**:
   - Parses GTF files for exon and CDS annotations.
   - Visualizes annotations as bars in a dedicated panel.

4. **Output**:
   - Saves the visualization as a high-resolution PNG file.

---

## **Dependencies**
- **Python 3.x**
- **Matplotlib** (`pip install matplotlib`)
- **NumPy** (`pip install numpy`)

---

## **Input Files**
1. **PSL File**:
   - Contains transcript alignment data (e.g., `BME163_Input_Data_6.psl`).
2. **GTF File**:
   - Contains gene annotations (e.g., `gencode.vM12.annotation.gtf`).

---

## **Usage**
1. **Input Arguments**:
   - `--outputFile` (`-o`): Name of the output PNG file. Default: `g.png`
   - `--inputFile` (`-i`): Path to the PSL input file. Default: `/path/to/your/data.psl`
   - `--outTextFile` (`-t`): Name of the output text file. Default: `Lecture4.data`

2. **Execution**:
   Run the script in a Python environment:
   ```bash
   python genome_visualization.py --outputFile output.png --inputFile input.psl --outTextFile output.txt
