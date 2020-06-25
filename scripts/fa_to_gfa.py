import sys

if __name__ == "__main__":

    id_ctr = 0
    print("\t".join(["H", "VZ:i:2.0"]))
    with open(sys.argv[1],"r") as ifi:
        for line in ifi:
            if not line.startswith(">"):
                id_ctr += 1
                print("\t".join(["S", str(id_ctr), str(len(line.strip())), line.strip(), ""]))

