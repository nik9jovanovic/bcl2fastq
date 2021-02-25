//
//  SampleSheetEntry.h
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 8/17/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#pragma once

#include <string>
#include <bitset>
using namespace std;

//Container class for storing information about SampleSheet.csv entries
class SampleSheetEntry{
    bitset<3>* index1;
    //index[0] - 0=valid, 1=invalid; index[1],index[2] - bits for base
    bitset<3>* index2;
    string sample_id;
    string sample_name;

public:
    SampleSheetEntry();
    SampleSheetEntry(const string sid, const string sname, bitset<3>* i1, bitset<3>* i2);
    SampleSheetEntry(const SampleSheetEntry &s);
    ~SampleSheetEntry();
    
    bitset<3>* getIndex1() const;
    bitset<3>* getIndex2() const;
    string getSampleID() const;
    string getSampleName() const;
    
    bool compareIndexes(SampleSheetEntry* rs, unsigned int index1_length, unsigned int index2_length) const;
};
