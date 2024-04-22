import subprocess

# Path to the .exe program
exe_path = "./k-medoids_parallelized_windows.exe"

# File to redirect output to
output_file = "results/Results_Medoids_windows_2.csv"

nthread = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"]
npunti = ["500", "1000", "2500", "5000", "10000"]
ncluster = ["3", "5", "10"]

ntest = 16*5*3
test = 1

# Open the file for appending
with open(output_file, "a") as file:
    for p in npunti:
        for c in ncluster:
            for t in nthread:
                # Command-line arguments
                arguments = [p, "2", c, t]
                # Launch the .exe program with arguments and redirect output to file
                subprocess.run([exe_path] + arguments, stdout=file, stderr=subprocess.PIPE)
                print("TEST COMPLETATI: [", test, "/", ntest,"]")
                test+=1
                
