#include <vector>
#include <string>
#include <map>
#include <cstdint>
#include <unordered_map>
#include <zlib.h>
#include "gfakluge.hpp"
#include "kseq.h"
#include "tinyfa.hpp"


/**
 * Sketch:
 * Create a breakpoint for each variant
 * Create a breakpoint for the maxNodeSizes
 * Iterate along the breakpoints
 * set the node at each breakpoint
 * Make a node, give it an id, etc
 * maintain nodes in previous level
 * if we are past the first node, wire in edges as we go
 * along.
 * wire up snp/insertion edges first
 *
 * Lastly, wire up deletio / inversion edges and emit the graph
 *
 */
struct VCF_Variant{
    VCF_Variant(){

    };
    VCF_Variant(std::string line){

    };
    std::string seq;
    int pos;
    std::string ref;
    std::vector<std::string> alts;
    std::map<std::string, std::string> infos;
};

std::pair<int, int> variant_to_breakpoint(const VCF_Variant& var){
        int front = 0;
        int back = 0;
        if ( !var.info["SVTYPE"].empty() && var.info["SVLEN"].empty()){
            cerr << "No end info. Skipping sv." << endl;
            return make_pair(-1, -1);
        }
        if (var.info["SVTYPE"] == "DEL"){
            front = var.pos - 1;
            back = var.pos - 1 + (uint32_t) stoi(var.info["SVLEN"]);
        }
        else if (var.info["SVTYPE"] == "INS"){
            front = var.pos - 1;
            back = var.pos  - 1;
        }
        else if (var.info["SVTYPE"] == "DUP"){
            front = var.pos - 1;
            back = var.pos - 1 + (uint32_t) stoi(var.info["SVLEN"]);

        }
        else if (var.info["SVTYPE"] == "INV"){
            cerr << "inversion parsed" << endl;
            front = var.pos - 1;
            back = var.pos - 1 + (uint32_t) stoi(var.info["SVLEN"]);
        }
        else{
            front = var.pos - 1;
            back = var.pos - 1 + var.ref.length();
        }
        return std::make_pair(front, back);
};

inline void make_contig_to_breakpoints(char* vcf_file,
        map<string, vector<*VCF_Variant>>& contig_to_variants,
        map<string, map<int, vector<VCF_Variant*>>>& contig_to_breakpoints_to_variants,
        map<string, vector<int>>& contig_to_breakpoints,
        vector<*VCF_Variant>& insertions){
    
    ifstream ifi;
    ifi.open(vcf_file);

    if (ifi.good()){
       
        std::string line;
        while (std::getline(ifi, line)){
            if (line[0] != "#"){
                VCF_Variant* vv = new VCF_Variant(line);
                std::string contig = vv->seq;
                std::pair<int, int> breaks = variant_to_breakpoint(*VCF_Variant);
                if (breaks.first != -1){
                    contig_to_breakpoints[contig].push_back(breaks.first);
                    contig_to_breakpoints[contig].push_back(breaks.second);
                    contig_to_breakpoints_to_variants[contig][breaks.first].push_back(vv);
                }
                    
                if (vv.info["SVTYPE"].empty()){
                    insertions.push_back(vv);
                }
                else if (vv.info["SVTYPE"] == "INS"){
                    insertions.push_back(vv);
                }
            }
        }
    }
};

void make_breakpoints(const std::string& contig, char* fasta_file,
        const std::vector<VCF_Variant*>& variants,
        std::vector<int>& breakpoints,
        map<int, vector<VCF_Variant*>& bp_to_variants,
        const int& max_node_length){

    int numvars = variants.size();
    breakpoints.reserve(numvars * 3);
    for (int i = 0; i < numvars; ++i){
        std::pair<int, int> bps = variant_to_breakpoint(*variants[i]);
        breakpoints.push_back(bps.first);
        breakpoints.push_back(bps.second);
        bp_to_variants[bps.first].push_back(variants[i]);
    }

    int seq_len = 0; 
    TFA::tiny_faidx_t tf;
    if (checkFAIndexFileExists(fasta_file)){

    }
    else{
        
    }

    tf.getSequenceLength(tf, contig.c_str(), seq_len);
    
    for (int i = max_node_length; i < seq_len; i += max_node_length){
        breakpoints.push_back(i);
    }

}

void construct_contig_graph(char* fasta_file, vector<int>& breakpoints, map<int, vector<VCF_Variant>>& bp_to_variants, GFAKluge gg){
    
    int numbp = breakpoints.size();
    int prev = 0;
    for (int i = 0; i < numbp; ++i){

        if (i > 0){

        }
    }


};

void construct_gfa(char* fasta_file, char* vcf_file){

};

