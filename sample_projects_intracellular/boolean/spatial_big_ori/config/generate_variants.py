import os

def generate_xml_variants(input_file, concentrations, diffusions):
    # Get the directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Read the original XML content
    with open(input_file, 'r') as file:
        xml_content = file.read()

    # Get the base name of the file without the extension
    base_name = os.path.splitext(os.path.basename(input_file))[0]

    # Generate variants
    variant_index = 0
    for c in concentrations:
        for d in diffusions:
            # Replace placeholders with current concentration, diffusion values, and output name
            variant_content = xml_content.replace("%%C%%", str(c)).replace("%%D%%", str(d))
            output_placeholder = f"output_{variant_index}"
            variant_content = variant_content.replace("%%output%%", output_placeholder)
            
            # Define the output file name with the original name and variant index
            output_file = os.path.join(script_dir, f"{base_name}_{variant_index}.xml")
            
            # Write the modified XML content to the new file
            with open(output_file, 'w') as file:
                file.write(variant_content)
            
            print(f"Generated variant: {output_file}")
            variant_index += 1

# Example usage
input_file = "PhysiCell_settings.xml"  # Replace with your actual XML file path
concentrations = [0.0000001,0.000001,0.00001]  # List of concentration values to use
diffusions = [2000]   # List of diffusion values to use

generate_xml_variants(input_file, concentrations, diffusions)

