import subprocess

# Path to the .exe program
exe_path = "./k-means_openmp.exe"

# File to redirect output to
output_file = "Results.csv"

nthread = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]
npunti = ["20", "2000", "20000", "2000000", "4000000"]

# Open the file for appending
with open(output_file, "a") as file:
    for t in nthread:
        for p in npunti:
            for i in range(35):
                # Command-line arguments
                arguments = [p, "2", "3", t]

                # Launch the .exe program with arguments and redirect output to file
                subprocess.run([exe_path] + arguments, stdout=file, stderr=subprocess.PIPE)
