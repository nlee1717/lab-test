import sys

counts_file = sys.argv[1]
counts_repaired_file = counts_file + "_repaired"

counts_repaired = ["bp\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous\n"]
empty_nt = "\t0\t0\t0\t0\t0\t0\t.\t.\n"

with open(counts_file, 'r') as cf:
    cf.readline()
    current_line = cf.readline()
    counter = 1
    reader = int(current_line.split()[0])

    while True:
        if counter < reader:
            counts_repaired.append(str(counter) + empty_nt)
        else:
            counts_repaired.append(current_line)
            current_line = cf.readline()
            if current_line == "":
                break
            else:
                reader = int(current_line.split()[0])
        counter += 1
cf.close()

with open(counts_repaired_file, 'w') as cpf:
    cpf.writelines(counts_repaired)
cpf.close()
