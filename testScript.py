import subprocess

# Path to the .exe program
exe_path = "./k-means_array.exe"

# File to redirect output to
output_file = "Results_Array.csv"

nthread = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"]
npunti = ["2000", "20000", "200000", "2000000", "20000000"]

ntest = 16*5*36
test = 1

# Open the file for appending
with open(output_file, "a") as file:
    for t in nthread:
        for p in npunti:
            for i in range(36):
                # Command-line arguments
                arguments = [p, t]

                # Launch the .exe program with arguments and redirect output to file
                subprocess.run([exe_path] + arguments, stdout=file, stderr=subprocess.PIPE)
                print("TEST COMPLETATI: [", test, "/", ntest,"]")
                test+=1
                
