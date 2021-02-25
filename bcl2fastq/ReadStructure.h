//
//  ReadStructure.h
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 8/15/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#pragma once

#include "CBCLReader.h"
#include <stdio.h>
#include <thread>

//Class for processing a read part
class ReadStructure{
    bool is_index;
    unsigned int length;
    unsigned int start;

    CBCLReader** cbcls;

public:
    ReadStructure();
    ReadStructure(const bool idx, const unsigned int l, const unsigned int s, CBCLReader** c);
    ReadStructure(const ReadStructure &rs);
    ~ReadStructure();
    
    bool getIndex() const;
    unsigned int getLength() const;
    unsigned int getStart() const;
    CBCLReader** getCBCLs() const;
    
    void clear();
    
    void readHeaders(const string project_path, const string bclname, const unsigned int lane, Semaphore* sem);
    void readCBCLs(const string lane_directory, const string bclname, const unsigned int tile, const unsigned int startpos, const unsigned int num_reads, Semaphore* sem);
};
