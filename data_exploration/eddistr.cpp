//
//  eddistr.cpp
//  
//
//  g++ -std=c++14 -m64 -g -O3 -I. eddistr.cpp -o eddistr
//

#include <getopt.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <cassert>
#include <string>

#include "wfa_lm.hpp"

using namespace std;
using namespace std::chrono;
using namespace wfalm;

void print_usage() {
    cerr << "usage: eddistr [options] seqs.fa" << endl;
    cerr << endl;
    cerr << "options:" << endl;
    cerr << " -n, --num-samples     how many distance to sample [1000]" << endl;
    cerr << " -o, --output          save identities to this file" << endl;
    cerr << " -h, --help            print this message and exit" << endl;
}

vector<string> parse_fasta(const string& fasta_name) {
    
    vector<string> to_return;
    
    string line;
    
    ifstream strm(fasta_name);
    if (!strm) {
        cerr << "ERROR: failed to open " << fasta_name << endl;
        exit(1);
    }
    bool seen_header = false;
    while (strm) {
        line.clear();
        getline(strm, line);
        if (!line.empty() && line.front() == '>') {
            to_return.emplace_back();
        }
        else {
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            to_return.back().append(line);
        }
    }
    
    return to_return;
}

int main(int argc, char **argv) {
    
    cerr << "Executing command:" << endl;
    for (size_t i = 0; i < argc; ++i) {
        if (i) {
            cerr << ' ';
        }
        cerr << argv[i];
    }
    cerr << endl;
    
    int n = 1000;
    string output;
    
    int c;
    while (true){
        static struct option long_options[] =
        {
            {"num-samples", required_argument, 0, 'n'},
            {"output", required_argument, 0, 'o'},
            {0,0,0,0}
        };
        
        int option_index = 0;
        c = getopt_long (argc, argv, "o:n:h",
                         long_options, &option_index);
        if (c == -1){
            break;
        }
        
        switch(c){
            case 'n':
                n = atoi(optarg);
                break;
            case 'o':
                output = optarg;
                break;
            case 'h':
            case '?':
                print_usage();
                return 0;
            default:
                print_usage();
                return 1;
        }
    }
    
    if (optind + 1 != argc) {
        // no positional arguments
        cerr << "ERROR: expected 1 positional argument, but got " << (argc - optind) << endl;
        print_usage();
        return 1;
    }
    
    ofstream out;
    if (!output.empty()) {
        out.open(output);
        if (!out.good()) {
            cerr << "ERROR: could not write to " << output << endl;
            return 1;
        }
    }
    
    string fasta(argv[optind++]);
    
    auto sequences = parse_fasta(fasta);
    
    random_device rd;
    auto seed = rd();
    mt19937 gen(seed);
    uniform_int_distribution<size_t> distribution(0, sequences.size() - 1);
    
    auto aligner = make_linear_wfaligner(1, 1);
    
    vector<double> identities;
    
    while (identities.size() < n) {
        size_t i = distribution(gen);
        size_t j = distribution(gen);
        if (i == j) {
            continue;
        }
        
        const auto& seq1 = sequences[i];
        const auto& seq2 = sequences[j];
        
        vector<CIGAROp> cigar;
        int32_t score;
        tie(cigar, score) = aligner.wavefront_align(seq1.c_str(), seq1.size(),
                                                    seq2.c_str(), seq2.size());
        
        identities.push_back(1.0 - double(score) / max(seq1.size(), seq2.size()));
    }
    
    double k1 = 0.0;
    double k2 = 0.0;
    for (auto ident : identities) {
        if (!output.empty()) {
            out << ident << '\n';
        }
        k1 += ident;
        k2 += ident * ident;
    }
    
    double mean = k1 / identities.size();
    double var = k2 / identities.size() - mean * mean;
    
    cerr << "mean identity: " << mean << endl;
    cerr << "std dev: " << sqrt(var) << endl;
    
    return 0;
}

