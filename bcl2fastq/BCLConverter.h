//
//  BCLConverter.h
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 7/11/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#pragma once

#include "ReadStructure.h"
#include "SampleSheetEntry.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <dirent.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <chrono>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;

struct DecodedBase{
    unsigned char base;
    unsigned char qual;
    
    DecodedBase(){
        base = 0;
        qual = 0;
    };
    DecodedBase(const unsigned char b, const unsigned char q){
        base = b;
        qual = q;
    };
};

class BCLConverter{
    //parameters from RunInfo.xml
    string instrument;
    string run_id;
    string flowcell_id;
    ReadStructure** read_parts;
    unsigned char index1_pos;
    unsigned char index2_pos;
    unsigned char read1_pos;
    unsigned char read2_pos;
    
    //parameters from SampleSheet.csv
    unsigned int umi1_length;
    unsigned int umi2_length;
    unsigned int num_samples;
    SampleSheetEntry** samples;
    
    //parameters from the directory structure
    string project_path;
    unsigned char num_lanes;
    unsigned int num_cycles;
    unsigned char bcls_per_cycle;
    vector<string> cbcl_names;
    
    //demultiplexing vector
    vector<vector<unsigned int>> demuxed_indexes;
    
    //multithreading parameters
    Semaphore* sem;
    unsigned int number_of_bytes;
    
public:
    BCLConverter();
    ~BCLConverter();
    
    void setProjectPath(const string path);
    void setThreads(const unsigned int threads);
    void setNumberOfBytes(const unsigned int bytes);
    
    bitset<3>* convertIndex(const string index, const unsigned int length) const;
    unsigned int checkAndAddBCLName(const string bcl_name);
    unsigned char decodeIndexBase(const unsigned char byte) const;
    DecodedBase decodeSequenceBase(const unsigned char byte, const unsigned int* bins) const;
    
    unsigned int parseRunInfo(const string runinfo_path = "");
    unsigned int parseSampleSheet(const string ss_path = "");
    unsigned int checkBCLDirectory();
    unsigned char* readLocs(const string path, const unsigned int startpos = 0);
    unsigned char* readFilter(const string path, const unsigned int startpos = 0);
    unsigned int demuxData(const unsigned char barcode_misses);
    unsigned int demuxPart(const unsigned int startpos, const unsigned int num_reads, const unsigned char barcode_misses, Semaphore* s);
    unsigned int writeData(const unsigned int lane, const unsigned int tile, const bool pair2, const unsigned char* locs_buffer, const unsigned char* filter_buffer);
    unsigned int writePart(const unsigned int number, const string path, const bool pair2, const string header, const unsigned char* locs_buffer, const unsigned char* filter_buffer);
    unsigned int convert(const unsigned char barcode_misses = 1, const bool read_filters = true, const bool read_positions = true);
};
