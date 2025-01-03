"""
# **Gene Expression Analysis Tool**

## **Overview**
This tool analyzes gene expression data from control and treatment conditions. It performs normalization, computes statistical metrics, calculates log2 fold changes, and evaluates differential gene expression significance.

---

## **Features**
1. **Normalization**:
   - Normalizes gene expression counts by dividing each gene's expression by the total expression in each sample.

2. **Mean and Median Expression**:
   - Computes mean and median normalized expression for each gene in control and treatment groups.

3. **Log2 Fold Change**:
   - Calculates log2 fold changes to assess relative expression differences.

4. **Statistical Testing (Tier 2)**:
   - Performs a Mann-Whitney U test for each gene to calculate p-values for significant expression changes.

5. **Output**:
   - Generates a ranked, tab-delimited file of gene statistics, including fold changes and p-values.

---

## **Dependencies**
- **Python 3.x**

---

## **Input Files**
1. **Control Directory**:
   - Contains 100 CSV files for control conditions.
2. **Treatment Directory**:
   - Contains 100 CSV files for treatment conditions.

---

## **Usage**
1. **Input Arguments**:
   - `<path_to_controls>`: Path to the control directory.
   - `<path_to_treatments>`: Path to the treatment directory.

2. **Execution**:
   Run the script in a Python environment:
   ```bash
   python Firstname_Lastname_Tier1.py <path_to_controls> <path_to_treatments>
