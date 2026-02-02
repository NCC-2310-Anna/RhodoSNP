#Read the files
def ReadFile(file):
     fobj_file = open(file, "r")
     file = fobj_file.readlines()
     fobj_file.close()
     return file
    