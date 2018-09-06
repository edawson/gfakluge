#ifndef PLIIB_D_H
#define PLIIB_D_H
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <functional>
#include <omp.h>
#include <functional>

using namespace std;

namespace pliib{
    // Char table to test for canonical bases
    static const int valid_dna[127] = {
        1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,0,1,0,1,1,1,
        0,1,1,1,1,1,1,1,1,1,
        1,1,1,0,1,1,1,1,1,1,
        1,1,1,1,1,1,0,1,0,1,
        1,1,0,1,1,1,1,1,1,1,
        1,1,1,1,1,0,1,1,1,1,
        1,1,1,1,1,1
    };

    // Reverse complement lookup table
    static char complement_array [26] = {
        84, 66, 71, 68, 69,
        70, 67, 72, 73, 74,
        75, 76, 77, 78, 79,
        80, 81, 82, 83, 65,
        85, 86, 87, 88, 89, 90
    };

    enum BREAKEND_TYPES {
        DEL,
        INS,
        INV,
        DUP,
        INTERCHROM,
        COMPLEX
    };

    // Check a string (as a char*) for non-canonical DNA bases
    inline bool canonical(const char* x, int len){
        bool trip = false;
        for (int i = 0; i < len; ++i){
            trip |= valid_dna[ static_cast<int> (x[i]) ];
        }
        return !trip;
    };

    inline bool canonical(string seq){
        const char* x = seq.c_str();
        int len = seq.length();
        return canonical(x, len);
    };


    inline void reverse_complement(const char* seq, char* ret, int len){
        for (int i = len - 1; i >=0; i--){
            ret[ len - 1 - i ] = (char) complement_array[ (int) seq[i] - 65];
        }
    }

    inline char base_complement(char c){
        return (char) complement_array[ (int) c - 65];    
    }

    /* Capitalize all characters in a string */
    /* Capitalize a C string */
    inline void to_upper(char* seq, int length){
        for (int i = 0; i < length; i++){
            char c = seq[i];
            seq[i] = ( (c - 91) > 0 ? c - 32 : c);
        }
    };

    /* Capitalize a string */
    inline string to_upper(string& seq){
        for (size_t i = 0; i < seq.length(); i++){
            char c = seq[i];
            seq[i] =  ((c - 91) > 0 ? c - 32 : c);
        }
        return seq;
    };

    inline void countChar(const char* s, char c, int& ret) {
        ret = 0;
        while(*s++ != '\0') { //Until the end of the string
            if(*s == c) {
                ++ret;
            }
        }
    };

    inline void destroy_splits(char**& splits, const int& num_splits, int*& split_sizes){
        delete [] splits;
        delete [] split_sizes;
    };

    // Modified from: https://techoverflow.net/2017/01/23/zero-copy-in-place-string-splitting-in-c/
    inline void split(char*& s, char delimiter, char**& ret, int& retsize, int*& split_sizes){
        int num_delim = 0;
        countChar(s, delimiter, num_delim);

        ret = new char*[num_delim + 1];
        retsize = num_delim + 1;
        split_sizes = new int[num_delim + 1];

        ret[0] = s;

        int i = 1;
        char* hit = s;
        while((hit = strchr(hit, delimiter)) != NULL){
            *hit = '\0';
            ++hit;

            ret[i] = hit;
            // Save the length of each string as well.
            // catch a special case: the string is the empty string
            if (ret[i - 1][0] == '\0'){
                *(split_sizes + i - 1) = 0;
            }
            else{
                *(split_sizes + (i - 1)) = strlen(ret[i - 1]); 
            }
            ++i;
        }
        // We need to handle the final string
        split_sizes[retsize - 1] = strlen(ret[retsize - 1]);
    };

    inline void split(string s, char delim, vector<string>& ret){

        int slen = s.length();
        char* s_to_split = new char[slen + 1];
        strncpy(s_to_split, s.c_str(), slen);
        s_to_split[slen] = '\0';

        char** splitret;
        int retsz;
        int* splitsz;

        split(s_to_split, delim, splitret, retsz, splitsz);
        ret.resize(retsz);

        for (int i = 0; i < retsz; ++i){
            ret[i].assign(string( splitret[i])); 
        }
        destroy_splits(splitret, retsz, splitsz);
        delete [] s_to_split;

    };

    inline vector<string> split(const string s, const char delim){
    
        vector<string> ret;
        int slen = s.length();
        char* s_to_split = new char[slen + 1];

        strncpy(s_to_split, s.c_str(), slen);
        s_to_split[slen] = '\0';

        char** splitret;
        int retsz;
        int* splitsz;


        split(s_to_split, delim, splitret, retsz, splitsz);

        ret.resize(retsz);

        for (int i = 0; i < retsz; ++i){
            ret[i].assign(string( splitret[i])); 
        }
        destroy_splits(splitret, retsz, splitsz);
        delete [] s_to_split;

        return ret;

    };
    inline vector<string> slow_split(string s, char delim){
        vector<string> ret;
        stringstream sstream(s);
        string temp;
        while(getline(sstream, temp, delim)){
            ret.push_back(temp);
        }
        return ret;

    }
    inline string join(const vector<string>& splits, char glue){
        stringstream ret;
        for (size_t i = 0; i < splits.size(); i++){
            if (i != 0){
                ret << glue;
            }
            ret << splits[i];
        }

        return ret.str();
    }

    /** Returns a string s', which is the substring of s before the
     * first appearance of 'delim'.
     * If delim is the first character, returns an empty string. **/

    inline void trim_after_char(char*&s, const int& len, char delim){
        int d_index = 0;
        while(  d_index < len && s[d_index] != delim){
            ++d_index;
        }
        char* ret = new char[d_index + 1];
        ret[d_index] = '\0';
        if (d_index > 0){
            memcpy(ret, s, d_index * sizeof(char));
        }
        delete [] s;
        s = ret;

    };

    /** Removes a character 'r' when it is seen at either the start or end of a string
     *  Returns a new string with all such occurrences of 'r' removed.
     *  s is destroyed in the process and replaced by a new string.
     *  Works like python's "strip(char)" function. **/
    inline void strip(char*& s, const int& len, char r){
        int start = 0;
        int end = len - 1;
        while (start < len && s[start] == r){
            ++start;
        }
        while(end > 0 && s[end] == r){
            --end;
        }
        int rsz = end - start + 1;
        char* ret = new char[rsz + 1];
        memcpy(ret, s + start, rsz * sizeof(char));
        // cerr << ret << endl;
        ret[rsz] = '\0';
        delete [] s;
        s = ret;
    };

    inline void strip(std::string& s, char r){
        int slen = s.length();
        char* t = new char[slen + 1];
        memcpy(t, s.c_str(), slen);
        strip(t, slen, r);
        s = std::string(t);
    };


    /** Removes a character from within a string **/
    inline void remove_char(char*& s, const int& len, char r){
        int write_index = 0;
        int read_index = 0;
        for (int i = 0; i < len, read_index < len; ++i){
            if (s[i] != r){
                s[write_index] = s[read_index];
                ++write_index;
            }
            ++read_index;
        }
        s[write_index] = '\0';
    };

    /** Removes any null or whitespace characters from within a string
     *  where whitespace is ' ' or '\t' or '\n'.
     *  Complexity: linear time in |s| **/
    inline void remove_nulls_and_whitespace(char*& s, const int& len){
        size_t write_index = 0;
        for (size_t i = 0; i < len; ++i){
            if (s[i] != '\n' && s[i] != '\t' && s[i] != '\0' && s[i] != ' '){
                s[write_index] = s[i];
                ++write_index;
            }
        }
        s[write_index] = '\0';
    };

    inline void paste(const char** splits, const int& numsplits, const int* splitlens, char*& ret){

        int sz = 0;
        for (int i = 0; i < numsplits; ++i){
            sz += splitlens[i];
        }
        ret = new char[sz + 1];
        ret[sz] = '\0';

        int r_pos = 0;
        for (int i = 0; i < numsplits; ++i){
            for (int j = 0; j < splitlens[i]; ++j){
                ret[r_pos] = splits[i][j];
                ++r_pos;
            }
        }
    };

    inline string join(vector<string> splits, string glue){
        stringstream ret;
        for (int i = 0; i < splits.size(); i++){
            if (i != 0){
                ret << glue;
            }
            ret << splits[i];
        }

        return ret.str();
    }

    inline void parse_breakend_field(const char* bend,
        const int& len,
        char*& contig,
        uint32_t& position,
        int& type,
        bool forward = false){
            int first_bracket_index = -1;
            int last_bracket_index = len;
            int colon_index = -1;

            int i = 0;
            while(i < len){
                char c = bend[i];
                if (c == ']' || c == '['){
                    if (first_bracket_index == -1){
                        first_bracket_index = i;
                    }
                    else{
                        last_bracket_index = i;
                        break;
                    }
                }
                if (c == ':'){
                    colon_index = i;
                }
                ++i;
            }
        if (contig != NULL){
            delete [] contig;
        }
        contig = new char[colon_index - first_bracket_index];
        contig[colon_index - first_bracket_index - 1] = '\0';
        memcpy(contig, bend + first_bracket_index + 1, colon_index - first_bracket_index - 1);
        
        char* pstr = new char[last_bracket_index - colon_index];
        pstr[last_bracket_index - colon_index - 1] = '\0';
        memcpy(pstr, bend + colon_index + 1, last_bracket_index - colon_index - 1);
        position = atoi(pstr);

    };


    /** 
     * Applies a lambda function <lambda> to all elements of a vector v, returning
     * a results vector the same size and type as v.
     */
    template <typename DataType, typename A>
        inline std::vector<DataType, A> p_vv_map(const std::vector<DataType, A> v, typename std::function<DataType(DataType)> lambda){
            std::vector<DataType> results(v.size());
            size_t sz = v.size();
            #pragma omp parallel for
            for (size_t i = 0; i < sz; ++i){
                results[i] = lambda(v[i]);
            }

            return results;
        }

    /**
     *  Applies a lambda function <lambda> to all elements of v, and returns those elements in v
     *  for which lambda returns true.
     */
    template<typename DataType, typename A>
        inline std::vector<DataType, A> p_vv_filter(const std::vector<DataType, A>& v, typename std::function<bool(DataType)> lambda){
            std::vector<DataType> results;
            size_t sz = v.size();
            for (size_t i = 0; i < sz; ++i){
                if (lambda(v[i])){
                    results.push_back(v[i]);
                }
            }
            return results;
        };


    /**
     * This function applies a lambda function L to all elements of
     * vector v, modifying the elements of v in place.
     */
    template<typename DataType, typename A>
        inline void p_vv_apply(std::vector<DataType, A>& v, typename std::function<DataType(DataType)> lambda){
            #pragma omp parallel for //private(i)
            for (size_t i = 0; i < v.size(); i++){
                auto r = lambda(v[i]);

                #pragma omp atomic write
                v[i] = r;
            }
        }

}
#endif
