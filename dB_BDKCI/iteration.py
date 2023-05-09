import os


folderPath = "./matrices/"

fileNames = []
for filename in os.listdir(folderPath):
    if filename.endswith(".txt"):
        fileNames.append(filename)



command = "python3 run_bp_xor3.py -xor4c 2.4 -xor3c 1.625 -xor2c 1.0 -iterations 1000 -path matrices -matrix "


for i in range(len(fileNames)):
    os.system(command + fileNames[i])

    


