import sys

if __name__ == "__main__":
    
    node_file_name = "nodes.tsv"
    edge_file_name = "edges.tsv"
    name_to_seq = {}
    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
            line = line.strip()
            if line.startswith("H"):
                continue
            elif line.startswith("S"):
                tokens = line.split("\t")
                name_to_seq[tokens[1]] = tokens[3]
            elif line.startswith("E"):
                tokens = line.split("\t")
                print(tokens[2][0:-1], tokens[3][0:-1])

    with open(node_file_name, "w") as ofi:
        for i in name_to_seq:
            l = "\t".join([i, name_to_seq[i]]) + "\n"
            ofi.write(l)

