import argparse
import os
import pandas as pd
import plotly.express as px

def parse_arguments():
    parser = argparse.ArgumentParser(description="Create heatmaps from CSV files.")
    parser.add_argument("-i", "--input", required=True, help="Path to input directory containing CSV files")
    parser.add_argument("-o", "--output", required=True, help="Path to output directory for heatmaps")
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    for root, _, files in os.walk(args.input):
        for file in files:
            if not file.endswith('.csv'):
                continue
                
            file_path = os.path.join(root, file)
            file_name = os.path.splitext(file)[0]  # Get filename without extension
            
            try:
                # Read the CSV file
                df = pd.read_csv(file_path)
                
                # Check if the DataFrame has an 'Unnamed: 0' column to use as index
                if 'Unnamed: 0' in df.columns:
                    df.set_index('Unnamed: 0', inplace=True)
                
                # Generate appropriate heatmap based on file type
                if 'GD_' in file:
                    fig = px.imshow(
                        df,
                        color_continuous_scale=px.colors.sequential.Blues,
                        range_color=[0, 1],
                        width = 800,
                        height = 600,
                    )
                    fig.update_layout(
                        title="Genetic Distance: " + file_name,
                        xaxis_title="Samples",
                        yaxis_title="Samples"
                    )
                    output_file = os.path.join(args.output, f"{file_name}_heatmap_GD.html")
                    
                elif 'CorrectedFRD_' in file:
                    fig = px.imshow(
                        df,
                        color_continuous_scale=px.colors.sequential.Blues,
                        range_color=[0, 1],
                        width = 800,
                        height = 600,
                    )
                    fig.update_layout(
                        title="Corrected FRD: " + file_name,
                        xaxis_title="Samples",
                        yaxis_title="Samples"
                    )
                    output_file = os.path.join(args.output, f"{file_name}_heatmap_corrFRD.html")
                    
                elif 'CorrectedRD_' in file:
                    max_value = df.max().max()
                    fig = px.imshow(
                        df,
                        color_continuous_scale=px.colors.sequential.Blues,
                        range_color=[0, max_value],
                        width = 800,
                        height = 600,
                    )
                    fig.update_layout(
                        title="Corrected RD: " + file_name,
                        xaxis_title="Samples",
                        yaxis_title="Samples"
                    )
                    output_file = os.path.join(args.output, f"{file_name}_heatmap_corrRD.html")
                
                else:
                    # Skip files that don't match our patterns
                    continue
                
                # Save the heatmap
                fig.write_html(output_file)
                fig.write_image(str(output_file).replace('.html', '.png'), scale = 2)
                print(f"Created heatmap: {output_file}")
                
            except Exception as e:
                print(f"Error processing {file}: {str(e)}")

if __name__ == "__main__":
    main()