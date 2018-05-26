#ifndef GFA_BUILDER_HPP
#define GFA_BUILDER_HPP
#include <vector>
#include <string>
#include <map>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <zlib.h>
#include "gfakluge.hpp"
#include "tinyfa.hpp"
#include "pliib.hpp"


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
        std::vector<string> splits = pliib::split(line, '\t');

    };
    std::string seq;
    int pos;
    std::string ref;
    std::vector<std::string> alts;
    std::map<std::string, std::string> info;
};

inline void set_gfa_edge_defaults(gfak::edge_elem& e){
    e.source_name = "";
    e.id = "*";
    e.sink_name = "";
    e.source_orientation_forward = true;
    e.sink_orientation_forward = true;
    e.type = 1;
    e.alignment = "0M";
    e.ends.set(0, 1);
    e.ends.set(0, 1);
    e.ends.set(2, 0);
    e.ends.set(3, 0);
};



inline std::pair<int, int> variant_to_breakpoint(const VCF_Variant& var){
    int front = 0;
    int back = 0;
    if ( !var.info.count("SVTYPE") && var.info.count("SVLEN")){
        cerr << "No end info. Skipping sv." << endl;
        return make_pair(-1, -1);
    }
    if (var.info.at("SVTYPE") == "DEL"){
        front = var.pos - 1;
        back = var.pos - 1 + (uint32_t) stoi(var.info.at("SVLEN"));
    }
    else if (var.info.at("SVTYPE") == "INS"){
        front = var.pos - 1;
        back = var.pos  - 1;
    }
    else if (var.info.at("SVTYPE") == "DUP"){
        front = var.pos - 1;
        back = var.pos - 1 + (uint32_t) stoi(var.info.at("SVLEN"));

    }
    else if (var.info.at("SVTYPE") == "INV"){
        cerr << "inversion parsed" << endl;
        front = var.pos - 1;
        back = var.pos - 1 + (uint32_t) stoi(var.info.at("SVLEN"));
    }
    else{
        front = var.pos - 1;
        back = var.pos - 1 + var.ref.length();
    }
    return std::make_pair(front, back);
};

inline void make_contig_to_breakpoints(char* vcf_file,
        std::map<string, std::vector<VCF_Variant*>>& contig_to_variants,
        std::map<string, std::map<int, std::vector<VCF_Variant*>>>& contig_to_breakpoints_to_variants,
        std::map<string, std::vector<int>>& contig_to_breakpoints,
        std::vector<VCF_Variant*>& insertions){

    std::ifstream ifi;
    ifi.open(vcf_file);

    if (ifi.good()){

        std::string line;
        while (std::getline(ifi, line)){
            if (line[0] != '#'){
                VCF_Variant* vv = new VCF_Variant(line);
                std::string contig = vv->seq;
                std::pair<int, int> breaks = variant_to_breakpoint(*vv);
                if (breaks.first != -1){
                    contig_to_breakpoints[contig].push_back(breaks.first);
                    contig_to_breakpoints[contig].push_back(breaks.second);
                    contig_to_breakpoints_to_variants[contig][breaks.first].push_back(vv);
                }

                if (vv->info.at("SVTYPE").empty()){
                    insertions.push_back(vv);
                }
                else if (vv->info.at("SVTYPE") == "INS"){
                    insertions.push_back(vv);
                }
            }
        }
    }
};

inline void make_breakpoints(const std::string& contig, char* fasta_file,
        const std::vector<VCF_Variant*>& variants,
        std::vector<int>& breakpoints,
        map<int, vector<VCF_Variant*>>& bp_to_variants,
        const int& max_node_length){

    int numvars = variants.size();
    breakpoints.reserve(numvars * 3);
    for (int i = 0; i < numvars; ++i){
        std::pair<int, int> bps = variant_to_breakpoint(*variants[i]);
        breakpoints.push_back(bps.first);
        breakpoints.push_back(bps.second);
        bp_to_variants[bps.first].push_back(variants[i]);
    }

    uint32_t seq_len = 0; 
    TFA::tiny_faidx_t tf;
    if (TFA::checkFAIndexFileExists(fasta_file)){
        TFA::parseFAIndex(fasta_file, tf);
    }
    else{
        TFA::createFAIndex(fasta_file, tf);   
    }

    getSequenceLength(tf, contig.c_str(), seq_len);

    for (int i = max_node_length; i < seq_len; i += max_node_length){
        breakpoints.push_back(i);
    }
    breakpoints.push_back(seq_len);

    std::set<int> bpset (breakpoints.begin(), breakpoints.end());

    breakpoints = std::vector<int>(bpset.begin(), bpset.end());

    std::sort(breakpoints.begin(), breakpoints.end());

}

inline void construct_contig_graph(string contig_name,
        char* contig_seq,
        uint32_t seq_len,
        vector<int>& breakpoints,
        map<int, vector<VCF_Variant*>>& bp_to_variants,
        vector<VCF_Variant*> variants,
        char* insertion_fasta,
        gfak::GFAKluge gg){

    int numbp = breakpoints.size();
    int prev = 0;

    TFA::tiny_faidx_t* insert_tf;
    if (insertion_fasta != NULL){
        insert_tf = new TFA::tiny_faidx_t();

        if (TFA::checkFAIndexFileExists(insertion_fasta)){
            TFA::parseFAIndex(insertion_fasta, *insert_tf);
        }
        else{
            TFA::createFAIndex(insertion_fasta, *insert_tf);
        }
    }


    int current_id = 0;
    int current_edge_id = 0;
    int current_pos = 0;
    std::vector<int64_t> ins_nodes;
    bool snptrip;



    std::vector<uint32_t> contig_node_ids;
    std::map<std::string, std::unordered_set<int>> path_to_nodes;
    std::map<int, std::string> node_to_path;
    std::unordered_map<int, int> bp_to_node_id;


    for (int i = 0; i < numbp; ++i){

        int ins_offset = 0;
        int bp = breakpoints[i];
        ++current_id;

        gfak::sequence_elem s;
        s.sequence.assign(contig_seq + current_pos, bp - current_pos);
        s.id = current_id;
        node_to_path[current_id] = contig_name;
        path_to_nodes[contig_name].insert(s.id);

        if (i > 0){
            if (snptrip){
                int prev_ref = contig_node_ids.size() - 2;
                while ( path_to_nodes[contig_name].find(prev_ref) == path_to_nodes[contig_name].end()){
                    gfak::edge_elem e;
                    set_gfa_edge_defaults(e);
                    e.source_name = std::to_string(contig_node_ids[prev_ref]);
                    e.sink_name = std::to_string(current_id);
                    prev_ref--;
                }
                snptrip = false;
            }

            // SNPs / insertions aren't on the reference path
            if ( path_to_nodes[contig_name].find( contig_node_ids.back() ) == path_to_nodes[contig_name].end() ){
                int prev_ref = contig_node_ids.size() - 1;
                while ( path_to_nodes[contig_name].find(prev_ref) != path_to_nodes[contig_name].end() ){
                    gfak::edge_elem e;
                    set_gfa_edge_defaults(e);
                    e.source_name = std::to_string(contig_node_ids[prev_ref]);
                    e.sink_name = std::to_string(current_id);
                    --prev_ref;
                }
                gfak::edge_elem e;
                e.source_name = std::to_string(contig_node_ids[prev_ref]);
                e.sink_name = std::to_string(current_id);
                snptrip = true;
            }
            else{
                gfak::edge_elem e;
                e.source_name = std::to_string(contig_node_ids.back());
                e.sink_name = std::to_string(current_id);
            }
        }
        contig_node_ids.push_back(current_id);
        std::vector<VCF_Variant*> bp_vars = bp_to_variants[bp];
        for (int i = 0; i < bp_vars.size(); ++i){
            VCF_Variant* bvar = bp_vars[i];
            if (bvar->info.find("SVTYPE") == bvar->info.end()){
                //int max_alt_size = 0;
                for (int altp = 0; altp < bvar->alts.size(); ++altp){
                    gfak::sequence_elem s;
                    int prev_id = current_id;
                    s.sequence.assign(bvar->alts[altp]);
                    s.id = ++current_id;
                    //max_alt_size = max(max_alt_size, bvar->alts[altp]);

                    gfak::edge_elem e;
                    set_gfa_edge_defaults(e);
                    e.source_name = prev_id;
                    e.sink_name = current_id;
                    //ins_offset = max(ins_offset, max_alt_size);
                    ins_offset = max(ins_offset, (int) bvar->alts[altp].size());
                }
            }
            else if (bvar->info.at("SVTYPE") == "INS"){
               for (int altp = 0; altp < bvar->alts[altp].size(); ++altp){
                    gfak::sequence_elem ins_node;
                    std::string seq = bvar->alts[altp];
                    if (seq[0] == '<' && insertion_fasta != NULL){
                        pliib::strip(seq, '<');
                        pliib::strip(seq, '>');
                        if (insert_tf->hasSeqID(seq.c_str())){
                            char* ins_seq;
                            getSequence(*insert_tf, seq.c_str(), ins_seq);
                        }
                    }
                    else if (pliib::canonical(seq)){
                        
                    }
                    else{
                        seq = "N";
                        continue;
                    }
                    ins_node.sequence.assign(seq);
                    
                    ins_node.id = ++current_id;

                    gfak::edge_elem e;
                    set_gfa_edge_defaults(e);
                    e.source_name = to_string(s.id);
                    e.sink_name = to_string(ins_node.id);
                    ins_offset = max(ins_offset, (int) seq.length());

               }
            }
            else{
                continue;
            }
        }
        // update current_pos here??
    }

    for (auto vvar : variants){
        if (vvar->info.find("SVTYPE") != vvar->info.end()){
            if (vvar->info.at("SVTYPE") == "DEL"){
                gfak::edge_elem e_from;
                e_from.source_name = bp_to_node_id[ vvar->pos - 1 - 1];
                e_from.sink_name = bp_to_node_id[vvar->pos -1 + stoi(vvar->info.at("SVLEN"))];
            }
            else if (vvar->info.at("SVTYPE") == "INS"){
                gfak::edge_elem e_to;
                //e_to.source_id = insertion_id_to_node_id.at( vvar-> 
                // TODO this won't work as-is
                //
            }
            else if (vvar->info.at("SVTYPE") == "DUP"){

            }
            else if (vvar->info.at("SVTYPE") == "INV"){
                gfak::edge_elem e_from;
                gfak::edge_elem e_to;
            }
        }
    }

    for (auto p : path_to_nodes){

    }


};

inline void construct_gfa(char* fasta_file, char* vcf_file, gfak::GFAKluge gg){

    // Read in vcf file and transform to variants

    // For each contig, get the relevant sequence from 
    // the FASTA file

    // 

};

#endif
