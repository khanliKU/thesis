import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='matlabToInventor.py')
    # Adding optional argument
    parser.add_argument("-m", "--mFile", type=str, help = "Input Matlab .m file")
    parser.add_argument("-x", "--csvFile", type=str, help = "Inventor compatible .csv file")
    args = parser.parse_args()
    print(args)

    oLines = []
    
    # Using readlines()
    if args.mFile:
        input = open(args.mFile, 'r')
        Lines = input.readlines()
        
        skip = False
        count = 0
        # Strips the newline character
        for line in Lines:
            count += 1
            if "%{" in line and not skip:
                skip = True
            elif "%}" in line and skip:
                skip = False
            
            if skip:
                continue

            if "m2" in line: #Skip areas
                line = line.replace("m2","m*m")
            
            if "=" in line:
                if "e-3" in line:
                    line = line.replace("e-3","*0.001")
                tokens = line.split("=")
                name = tokens[0].strip()
                tokens = tokens[1].split("%")
                value = tokens[0].strip()
                unit = "ul"
                if len(tokens) > 1:
                    unit = tokens[1].strip()
                print("Param:{}\tValue:{}\tUnit:{}".format(name, value, unit))
                oLines.append("{},{},{}\n".format(name, value, unit))
        input.close()
                
    # writing to file
    if args.csvFile:
        output = open(args.csvFile, 'w')
        output.writelines(oLines)
        output.close()