import sys
import csv
import math

# Function to read and normalize data from multiple CSV files
def normalize_data(filepaths):
    normalized_data = {}  # Initialize a dictionary to store normalized gene data
    
    for filepath in filepaths:
        # Open each CSV file
        with open(filepath, 'r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header row
            total_expression = 0  # Variable to keep track of the total expression in a sample
            gene_data = {}  # Temporary dictionary to store the expression data for this file
            
            # Loop through each row of the CSV file, which represents a gene's expression
            for row in reader:
                gene = row[0]  # The first column is the gene name
                expression = float(row[1])  # The second column is the gene's expression value (convert to float)
                gene_data[gene] = expression  # Store the expression value for the gene
                total_expression += expression  # Add this gene's expression to the total for normalization
            
            # Normalize each gene's expression by dividing by the total expression for this sample
            for gene, expression in gene_data.items():
                normalized_expression = expression / total_expression  # Normalize by total expression
                if gene not in normalized_data:
                    normalized_data[gene] = []  # Initialize an empty list if this gene hasn't been seen before
                normalized_data[gene].append(normalized_expression)  # Append the normalized expression for this gene
    
    return normalized_data  # Return the dictionary containing normalized data for all genes

# Function to compute mean and median of normalized expression data
def compute_mean_median(data):
    summary = {}  # Dictionary to store mean and median for each gene
    for gene, expressions in data.items():
        mean_expression = sum(expressions) / len(expressions)  # Calculate the mean expression for each gene
        median_expression = sorted(expressions)[len(expressions) // 2]  # Calculate the median (middle value) of the expressions
        summary[gene] = (mean_expression, median_expression)  # Store the mean and median in a tuple for each gene
    
    return summary  # Return the summary of mean and median for each gene

# Function to compute the log2 fold change between control and treatment groups
def compute_log2_fold_change(control_summary, treatment_summary):
    log_fold_changes = {}  # Dictionary to store the log2 fold change for each gene
    for gene in control_summary.keys():
        control_mean = control_summary[gene][0]  # Get the mean expression for the gene in control group
        treatment_mean = treatment_summary[gene][0]  # Get the mean expression for the gene in treatment group
        
        # Avoid division by zero by checking if both means are greater than 0
        if control_mean > 0 and treatment_mean > 0:
            fold_change = treatment_mean / control_mean  # Compute the fold change (ratio of treatment to control)
            log2_fold_change = math.log2(fold_change)  # Calculate the log2 of the fold change
        else:
            log2_fold_change = float('-inf')  # If either mean is zero, assign negative infinity to the fold change
        
        log_fold_changes[gene] = log2_fold_change  # Store the log2 fold change for this gene
    
    return log_fold_changes  # Return the dictionary containing log2 fold changes for all genes

# Main function that orchestrates reading, processing, and outputting the gene expression data
def main(control_dir, treatment_dir):
    # Construct the file paths for control and treatment files ('experiment_#.csv')
    control_files = [f'{control_dir}/experiment_{i}.csv' for i in range(1, 101)]
    treatment_files = [f'{treatment_dir}/experiment_{i}.csv' for i in range(1, 101)]
    
    # Normalize gene expression data for both control and treatment groups
    control_normalized = normalize_data(control_files)
    treatment_normalized = normalize_data(treatment_files)
    
    # Compute the mean and median normalized expression for control and treatment groups
    control_summary = compute_mean_median(control_normalized)
    treatment_summary = compute_mean_median(treatment_normalized)
    
    # Compute log2 fold changes between control and treatment groups
    log_fold_changes = compute_log2_fold_change(control_summary, treatment_summary)
    
    # Write the results to a file in a tab-delimited format
    with open("Yashesha_Kothari_output_Tier1.txt", 'w') as output_file:
        # Write the header line
        output_file.write("#gene\t(mean normalized control expression)\t(median normalized control expression)\t(mean normalized treatment expression)\t(median normalized treatment expression)\t(logFoldChange)\n")
        
        # Write the data for each gene, sorted by log2 fold change
        for gene in sorted(log_fold_changes, key=log_fold_changes.get):
            control_mean, control_median = control_summary[gene]  # Get mean and median for control
            treatment_mean, treatment_median = treatment_summary[gene]  # Get mean and median for treatment
            log_fold_change = log_fold_changes[gene]  # Get log2 fold change
            
            # Write the gene data to the output file
            output_file.write(f"{gene}\t{control_mean:.6f}\t{control_median:.6f}\t{treatment_mean:.6f}\t{treatment_median:.6f}\t{log_fold_change:.6f}\n")

# Script entry point
if __name__ == "__main__":
    # Ensure the correct number of command-line arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python Firstname_Lastname_Tier1.py path/to/controls path/to/treatments")
        sys.exit(1)  # Exit if the number of arguments is incorrect
    
    # Get the control and treatment directory paths from command-line arguments
    control_dir = sys.argv[1]
    treatment_dir = sys.argv[2]
    
    # Call the main function to process the gene expression data
    main(control_dir, treatment_dir)
