#ifndef GFA_BUILDER_HPP
#define GFA_BUILDER_HPP
#include <vector>
#include <string>
#include <map>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <set>
#include <zlib.h>
#include <ostream>
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
    std::string seq;
    int pos;
    std::string ref;
    std::vector<std::string> alts;
    std::map<std::string, std::string> info;
    std::string raw_info_line;

    VCF_Variant(){

    };
    VCF_Variant(std::string line){
        std::vector<std::string> splits = pliib::split(line, '\t');
        seq = splits[0];
        pos = stoi(splits[1]);
        ref = splits[3];
        alts = pliib::split(splits[4], ',');

        raw_info_line = splits[7];
        std::vector<std::string> isplits = pliib::split(splits[7], ';');
        for (auto i : isplits){
            std::vector<string> xsplits = pliib::split(i, '=');
            if (xsplits.size() == 1){
                this->info[xsplits[0]] = xsplits[0];
            }
            else{
                this->info[xsplits[0]] = xsplits[1];
            }
        }
    };
    std::string to_string(){
        std::stringstream st;
        st << seq << '\t' <<
            pos << '\t' << "." << '\t' << 
            ref << '\t' <<
            endl;
        return st.str();
    };
};

struct dummy_node{
    uint32_t id = 0;
    uint32_t length = 0;
    char* path = NULL;
};

struct bp_allele{
    uint32_t prev_node_length;
    uint64_t prev_node_id;
    vector<char*> alt_seqs;
    bool isRef = false;
};

inline void set_gfa_edge_defaults(gfak::edge_elem& e, uint32_t base_edge_id = 0){
    e.source_name = "";
    e.id = to_string(base_edge_id);
    e.sink_name = "";
    e.source_orientation_forward = true;
    e.sink_orientation_forward = true;
    e.type = 1;
    e.alignment = "0M";
    e.ends.set(0, 1);
    e.ends.set(1, 1);
    e.ends.set(2, 0);
    e.ends.set(3, 0);
};



inline void set_gfa_node_values(gfak::sequence_elem& s){
    s.name = to_string(s.id);
    s.length = s.sequence.length();
};



std::pair<int, int> variant_to_breakpoint(const VCF_Variant& var){
    int front = 0;
    int back = 0;
    if ( var.info.find("SVTYPE") == var.info.end() &&
        var.info.find("SVLEN") == var.info.end() && 
        var.info.find("END") == var.info.end()){
        // This is an SNV or indel;
        // We'll need a break at both ends of its REF
        // sequence, and we'll wire any edges/or nodes
        // to those breakpoints.
        front = var.pos - 1;
        back = var.pos - 1 + var.ref.length();
    }
    else if (var.info.at("SVTYPE") == "DEL"){
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
    else if (var.info.at("SVTYPE") == "TRA"){
        front = var.pos - 1;
        back = -1;
    }
    else if (var.info.at("SVTYPE") == "BND"){
        front = var.pos - 1;
        back = -1;
    }
    else{
        cerr << "No end info. Skipping variant." << endl;
        return make_pair(-1, -1);
    }
    return std::make_pair(front, back);
};

void make_contig_to_breakpoints(char* vcf_file,
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
                //cerr << vv->seq << " " << vv->pos << " " << breaks.first << " " << breaks.second << endl; 
                if (breaks.first != -1 && breaks.second != -1){
                    contig_to_breakpoints[contig].push_back(breaks.first);
                    contig_to_breakpoints[contig].push_back(breaks.second);
                    contig_to_breakpoints_to_variants[contig][breaks.first].push_back(vv);
                    contig_to_variants[contig].push_back(vv);
                    
                }

                if (vv->info.find("SVTYPE") == vv->info.end()){
                    insertions.push_back(vv);
                }
                else if ( vv->info.at("SVTYPE").empty() || 
                            vv->info.at("SVTYPE") == "SNP" ||
                            vv->info.at("SVTYPE") == "INS"){
                    insertions.push_back(vv);
                }
                else if (vv->info.at("SVTYPE") == "TRA" ||
                            vv->info.at("SVTYPE") == "BND"){
                    // Get the relevant information from the alt field.
                    // if our breaks are on the same chromosome, just tuck
                    // them in like usual.
                    cerr << "TRA/BND not yet implemented. Skipping variant." << endl;
                }
            }
        }
    }
    else{
        cerr << "ERROR [gfak build] : could not open VCF file " << vcf_file << endl;
        exit(9);
    }
};

void make_breakpoints(const std::string& contig, char* fasta_file,
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

    TFA::getSequenceLength(tf, contig.c_str(), seq_len);

    for (int i = max_node_length; i < seq_len; i += max_node_length){
        breakpoints.push_back(i);
    }
    breakpoints.push_back(seq_len);

    std::set<int> bpset (breakpoints.begin(), breakpoints.end());

    breakpoints = std::vector<int>(bpset.begin(), bpset.end());

    std::sort(breakpoints.begin(), breakpoints.end());

}

void construct_contig_graph(string contig_name,
        char* contig_seq,
        uint32_t seq_len,
        vector<int>& breakpoints,
        map<int, vector<VCF_Variant*>>& bp_to_variants,
        vector<VCF_Variant*>& variants,
        char* insertion_fasta,
        gfak::GFAKluge& gg,
        std::ostream& os,
        int& base_seq_id,
        int base_edge_id){

    char* dummy_name = new char[contig_name.size() + 1];
    memcpy(dummy_name, contig_name.c_str(), contig_name.size() + 1);

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


    int current_id = base_seq_id;
    uint32_t current_length = 0;
    int current_edge_id = base_edge_id;
    int current_pos = 0;
    std::vector<int64_t> ins_nodes;
    bool snptrip = false;




    std::vector<dummy_node*> contig_nodes;
    std::map<std::string, std::set<int>> path_to_nodes;
    std::unordered_map<int, int> bp_to_node_id;
    std::unordered_map<int, int> node_id_to_length;
    std::unordered_map<string, int> insertion_id_to_node_id;


    for (int i = 0; i < numbp; ++i){

        int ins_offset = 0;
        int bp = breakpoints[i];
        ++current_id;

        gfak::sequence_elem s;
        s.sequence.assign(contig_seq + current_pos, bp - current_pos);
        s.id = current_id;
        set_gfa_node_values(s);
        current_length = s.length;
        
        bp_to_node_id[current_pos] = current_id;
        bp_to_node_id[bp - 1] = current_id;
        node_id_to_length[current_id] = current_length;

        gg.write_element(os, s);
        os << endl;

        if (i > 0){
            if (snptrip){
                int prev_ref = 0;
                prev_ref = contig_nodes.size() - 2;
                while ( contig_nodes[prev_ref]->path == NULL){
                    // cerr << "WARN: not handling snp nodes right now" << endl;
                    // break;
                    gfak::edge_elem e;
                    set_gfa_edge_defaults(e, ++base_edge_id);
                    e.source_name = std::to_string(contig_nodes[prev_ref]->id);
                    e.sink_name = std::to_string(current_id);
                    e.source_begin = current_length;
                    e.source_end = current_length;
                    
                    gg.write_element(os, e);
                    os << endl;
                    prev_ref--;
                }
                snptrip = false;
            }

            // SNPs / insertions aren't on the reference path
            if ( contig_nodes.back()->path == NULL ){
                int prev_ref = contig_nodes.size() - 1;
                while ( (contig_nodes[prev_ref])->path == NULL ){
                    // SNP / insertion nodes are already connected properly.
                    // gfak::edge_elem e;
                    // set_gfa_edge_defaults(e, ++base_edge_id);
                    // e.source_name = std::to_string((contig_nodes[prev_ref])->id);
                    // e.sink_name = std::to_string(s.id);
                    //e.source_begin = contig_nodes[prev_ref]->length;
                    //e.source_end = contig_nodes[prev_ref]->length;

                    //gg.write_element(os, e);
                    //os << endl;
                    --prev_ref;
                }
                // Connect up to the last reference nodes.
                gfak::edge_elem e;
                set_gfa_edge_defaults(e, ++base_edge_id);
                e.source_name = std::to_string((contig_nodes[prev_ref])->id);
                e.sink_name = std::to_string(current_id);
                e.source_begin = contig_nodes[prev_ref]->length;
                e.source_end = contig_nodes[prev_ref]->length;
                snptrip = true;
                gg.write_element(os, e);
                os << endl;
            }
            else{
                gfak::edge_elem e;
                set_gfa_edge_defaults(e, ++base_edge_id);
                e.source_name = std::to_string((contig_nodes.back())->id);
                e.sink_name = std::to_string(current_id);
                e.source_begin = (contig_nodes.back())->length;
                e.source_end = (contig_nodes.back())->length;
                gg.write_element(os, e);
                os << endl;
            }
        }

        dummy_node* dn = new dummy_node();
        dn->id = s.id;
        dn->length = current_length;
        dn->path = dummy_name;

        contig_nodes.push_back(dn);
        path_to_nodes[contig_name].insert(s.id);

        std::vector<VCF_Variant*> bp_vars = bp_to_variants[bp];
        for (int i = 0; i < bp_vars.size(); ++i){
            VCF_Variant* bvar = bp_vars[i];
            if (bvar->info.find("SVTYPE") == bvar->info.end()){
                //int max_alt_size = 0;
                for (int altp = 0; altp < bvar->alts.size(); ++altp){
                    gfak::sequence_elem snp_s;
                    snp_s.sequence.assign(bvar->alts[altp]);
                    snp_s.id = ++current_id;
                    set_gfa_node_values(snp_s);
                    dummy_node* snn = new dummy_node();
                    snn->id = snp_s.id;
                    snn->length = snp_s.length;
                    snn->path = NULL;
                    contig_nodes.push_back(snn);

                    gg.write_element(os, snp_s);
                    os << endl;

                    gfak::edge_elem e;
                    set_gfa_edge_defaults(e, ++base_edge_id);
                    e.source_name = std::to_string( s.id );
                    e.sink_name = std::to_string(snp_s.id);
                    e.source_begin = s.length;
                    e.source_end = s.length;
                    //ins_offset = max(ins_offset, max_alt_size);
                    ins_offset = max(ins_offset, (int) bvar->alts[altp].size());

                    gg.write_element(os, e);
                    os << endl;

                    //exit(1);
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
                            TFA::getSequence(*insert_tf, seq.c_str(), ins_seq);
                        }
                    }
                    else if (pliib::canonical(seq)){
                        // Don't do anything, as the alt sequence was already valid DNA
                    }
                    else{
                        seq = "N";
                        continue;
                    }

                    string insert_id = "INS_" + to_string(bvar->pos) + "_" + to_string(altp);
                    ins_node.sequence.assign(seq);
                    ins_node.id = ++current_id;
                    insertion_id_to_node_id[insert_id] = ins_node.id;
                    set_gfa_node_values(ins_node);

                    gg.write_element(os, ins_node);
                    os << endl;

                    gfak::edge_elem e;
                    set_gfa_edge_defaults(e, ++base_edge_id);
                    e.source_name = to_string(s.id);
                    e.sink_name = to_string(ins_node.id);
                    ins_offset = max(ins_offset, (int) seq.length());

                    gg.write_element(os, e);
                    os << endl;

               }
            }
            else{
                continue;
            }
        }
        // update current_pos here??
        current_pos = breakpoints[i];
    }

    //for (auto k : bp_to_node_id){
    //    cout << k.first << " " << k.second << endl;
    //}
    
    for (auto vvar : variants){
        if (vvar->info.find("SVTYPE") != vvar->info.end()){
            
            if (vvar->info.at("SVTYPE") == "DEL"){
                gfak::edge_elem e_from;
                set_gfa_edge_defaults(e_from, ++base_edge_id);
                e_from.source_name = to_string(bp_to_node_id.at( vvar->pos - 1 - 1));
                e_from.sink_name = to_string(bp_to_node_id.at(vvar->pos -1 + stoi(vvar->info.at("SVLEN"))));
                e_from.source_begin = node_id_to_length.at(bp_to_node_id.at( vvar->pos - 1 - 1));
                e_from.source_end = node_id_to_length.at(bp_to_node_id.at( vvar->pos - 1 - 1));

                gg.write_element(os, e_from);
                os << endl;
            }
            else if (vvar->info.at("SVTYPE") == "INS"){
                gfak::edge_elem e_to;
                set_gfa_edge_defaults(e_to, ++base_edge_id);
                string ins_id;
                for (int altp = 0; altp < vvar->alts.size(); ++altp){
                    ins_id = "INS_" + to_string(vvar->pos) + "_" + to_string(altp);
                }
                e_to.source_name = to_string(bp_to_node_id.at(vvar->pos - 1));
                e_to.source_begin = node_id_to_length.at(bp_to_node_id.at( vvar->pos - 1));
                e_to.source_end = node_id_to_length.at(bp_to_node_id.at( vvar->pos - 1));
                e_to.sink_name = to_string(insertion_id_to_node_id.at(ins_id));
                // TODO double check this bit above, as it might be funky.

                gg.write_element(os, e_to);
                os << endl;
            }
            else if (vvar->info.at("SVTYPE") == "DUP"){
                gfak::edge_elem e_cycle;
                set_gfa_edge_defaults(e_cycle, ++base_edge_id);
                e_cycle.source_name = to_string(bp_to_node_id.at(vvar->pos - 1 - 1 + stoi(vvar->info.at("SVLEN"))));
                e_cycle.sink_name = to_string(bp_to_node_id.at(vvar->pos - 1));
                e_cycle.source_begin = node_id_to_length.at(bp_to_node_id.at(vvar->pos - 1));
                e_cycle.source_end = node_id_to_length.at(bp_to_node_id.at(vvar->pos - 1));
                gg.write_element(os, e_cycle);
                os << endl;
            }
            else if (vvar->info.at("SVTYPE") == "INV"){
                gfak::edge_elem e_from;
                gfak::edge_elem e_to;

                set_gfa_edge_defaults(e_from, ++base_edge_id);
                set_gfa_edge_defaults(e_to, ++base_edge_id);
                
                e_from.source_name = to_string(bp_to_node_id[vvar->pos - 1 - 1]);
                e_from.sink_name = to_string(bp_to_node_id[vvar->pos -1 + stoi(vvar->info.at("SVLEN")) - 1]);
                e_from.source_orientation_forward = true;
                e_from.sink_orientation_forward = false;
                e_from.source_begin = node_id_to_length.at(bp_to_node_id.at(vvar->pos - 1 - 1));
                
                e_to.source_name = to_string(bp_to_node_id[vvar->pos - 1]);
                e_to.sink_name = to_string(bp_to_node_id[ vvar->pos - 1 + stoi(vvar->info.at("SVLEN"))]);
                e_to.source_orientation_forward = false;
                e_to.sink_orientation_forward = true;
                // Set node begin / ends

                gg.write_element(os, e_from);
                os << endl;
                gg.write_element(os, e_to);
                os << endl;
            }
            else if (vvar->info.at("SVTYPE") == "TRA"){

            }
            else if (vvar->info.at("SVTYPE") == "BND"){
                
            }
        }
        else{
            // SNV / indels
        }
    }

    for (auto p : path_to_nodes){
        os << "O" << '\t' << p.first << '\t';
        vector<int> zs (p.second.begin(), p.second.end());
        vector<string> z_to_str(zs.size());
        for (int i = 0; i < zs.size(); ++i){
            z_to_str[i] = to_string(zs[i]) + "+";
        }
        string zstring = pliib::join(z_to_str, ' ');
        os << zstring << endl;
    }

    for (auto cc : contig_nodes){
        delete cc;
    }


};

void construct_gfa(char* fasta_file, char* vcf_file, char* insertion_fasta, gfak::GFAKluge gg, int max_node_length = 128){

    
    std::map<std::string, std::vector<VCF_Variant*>> contig_to_variants;
    std::map<std::string, std::map<int, std::vector<VCF_Variant*>>> contig_to_breakpoints_to_variants;
    std::vector<VCF_Variant*> insertions;
    std::map<string, std::vector<int>> contig_to_breakpoints;
    // Read in vcf file and transform to variants

    make_contig_to_breakpoints(vcf_file, 
                        contig_to_variants,
                        contig_to_breakpoints_to_variants,
                        contig_to_breakpoints,
                        insertions);

    //cerr << contig_to_variants.size() << " contigs to process" << endl;


    // For each contig, get the relevant sequence from 
    // the FASTA file
    TFA::tiny_faidx_t tf;
    if (TFA::checkFAIndexFileExists(fasta_file)){
        TFA::parseFAIndex(fasta_file, tf);
    }
    else{
        TFA::createFAIndex(fasta_file, tf);
        TFA::writeFAIndex(fasta_file, tf);
    }

    int base_seq_id = 0;
    int base_edge_id = 0;
    
    std:: cout << gg.header_string();
    for (auto contig : contig_to_variants){
        //cerr << "Processing contig: " << contig.first << endl;
        if (tf.hasSeqID(contig.first.c_str())){
            
            char* seq;
            uint32_t len;
            TFA::getSequence(tf, contig.first.c_str(), seq);
            TFA::getSequenceLength(tf, contig.first.c_str(), len);

            std::vector<VCF_Variant*> vars = contig.second;
            std::vector<int> bps;
            std::map<int, vector<VCF_Variant*>> bp_to_var;
            make_breakpoints(contig.first, fasta_file, vars, bps, bp_to_var, max_node_length);
            construct_contig_graph(contig.first, seq, len, bps,
                    bp_to_var, vars, insertion_fasta,
                    gg, std::cout,
                    base_seq_id, base_edge_id);


        }
    }
    


    /* void construct_contig_graph(string contig_name,
        char* contig_seq,
        uint32_t seq_len,
        vector<int>& breakpoints,
        map<int, vector<VCF_Variant*>>& bp_to_variants,
        vector<VCF_Variant*> variants,
        char* insertion_fasta,
        gg,
        cout); **/



};

#endif
