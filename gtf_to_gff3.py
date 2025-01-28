from Bio import SeqIO

def convert_gtf_to_gff3(input_gtf, output_gff3):
    """
    Converts a GTF file to GFF3 format using Biopython.

    Args:
        input_gtf (str): Path to the input GTF file.
        output_gff3 (str): Path to the output GFF3 file.
    """
    with open(input_gtf, "r") as input_handle, open(output_gff3, "w") as output_handle:
        # Write the GFF3 header
        output_handle.write("##gff-version 3\n")
        
        # Track parent-child relationships
        feature_dict = {}  # Stores parent features and their children
        
        for line in input_handle:
            # Skip comment lines
            if line.startswith("#"):
                continue
            
            # Split the GTF line into columns
            columns = line.strip().split("\t")
            if len(columns) != 9:
                continue  # Skip malformed lines
            
            # Extract attributes
            attributes = columns[8]
            attr_dict = {}
            for attr in attributes.strip(";").split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                key, value = attr.split(" ", 1)
                value = value.strip('"')
                attr_dict[key] = value
            
            # Handle ID and Parent attributes
            if "gene_id" in attr_dict:
                attr_dict["ID"] = attr_dict["gene_id"]
                if "transcript_id" in attr_dict:
                    attr_dict["Parent"] = attr_dict["gene_id"]
                    attr_dict["ID"] = attr_dict["transcript_id"]
            
            # Rebuild the attributes column in GFF3 format
            new_attributes = ";".join([f"{key}={value}" for key, value in attr_dict.items()])
            columns[8] = new_attributes
            
            # Write the modified line to the output file
            output_handle.write("\t".join(columns) + "\n")

if __name__ == "__main__":
    # Input and output file paths
    input_gtf = "input.gtf"  # Replace with your input GTF file path
    output_gff3 = "output.gff3"  # Replace with your output GFF3 file path
    
    # Run the conversion
    convert_gtf_to_gff3(input_gtf, output_gff3)
    print(f"Conversion complete! Output saved to {output_gff3}")