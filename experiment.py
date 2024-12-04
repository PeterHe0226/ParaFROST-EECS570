import os
import subprocess
import shutil

def run_parafrost_on_cnf_files(base_directory, results_directory):

    # Ensure the results directory exists and is empty
    if os.path.exists(results_directory):
        shutil.rmtree(results_directory)  
    os.makedirs(results_directory)  
    
    # Walk through all subdirectories
    for root, _, files in os.walk(base_directory):
        for file in files:
            if file.endswith(".cnf"):
                cnf_file_path = os.path.join(root, file)
                
                # Output file path in the results folder
                output_file_name = f"{file}.out"
                output_file_path = os.path.join(results_directory, output_file_name)
                
                try:
                    # Run the parafrost command
                    with open(output_file_path, 'w') as output_file:
                        subprocess.run(["parafrost", cnf_file_path], stdout=output_file, stderr=subprocess.STDOUT)
                    print(f"Processed: {cnf_file_path} -> {output_file_path}")
                except Exception as e:
                    print(f"Error processing {cnf_file_path}: {e}")

# Example usage
if __name__ == "__main__":
    base_directory = "/path/to/your/folder"  # Replace 
    results_directory = "/path/to/results/folder"  # Replace 
    run_parafrost_on_cnf_files(base_directory, results_directory)
