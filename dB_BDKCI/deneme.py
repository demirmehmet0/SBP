import os


folderPath = "./matrices"

fileNames = []
for filename in os.listdir(folderPath):
    if filename.endswith(".txt"):
        fileNames.append(filename)

for i in range(len(fileNames)):
    with open(folderPath + "/" + fileNames[i], 'r') as file:
        filedata = file.read()
    filedata = filedata.replace('12 12', '20 20')
    filedata = filedata.replace(']', '')
    with open(folderPath + "/" + fileNames[i], 'w') as file:
        file.write(filedata)
    
