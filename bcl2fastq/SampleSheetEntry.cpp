//
//  SampleSheetEntry.cpp
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 8/17/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#include "SampleSheetEntry.h"


SampleSheetEntry::SampleSheetEntry(){
    index1 = NULL;
    index2 = NULL;
    sample_id = "";
    sample_name = "";
};

SampleSheetEntry::SampleSheetEntry(const string sid, const string sname, bitset<3>* i1, bitset<3>* i2){
    index1 = i1;
    index2 = i2;
    sample_id = sid;
    sample_name = sname;
};

SampleSheetEntry::SampleSheetEntry(const SampleSheetEntry &s){
    index1 = s.getIndex1();
    index2 = s.getIndex2();
    sample_id = s.getSampleID();
    sample_name = s.getSampleName();
};

SampleSheetEntry::~SampleSheetEntry(){
    delete index1;
    delete index2;
};

bitset<3>* SampleSheetEntry::getIndex1() const{
    return index1;
};

bitset<3>* SampleSheetEntry::getIndex2() const{
    return index2;
};

string SampleSheetEntry::getSampleID() const{
    return sample_id;
};

string SampleSheetEntry::getSampleName() const{
    return sample_name;
};

bool SampleSheetEntry::compareIndexes(SampleSheetEntry* rs, unsigned int index1_length, unsigned int index2_length) const{
    unsigned int index1_it = 0;
    for(;index1_it < index1_length; index1_it++){
        if(((index1[index1_it][0] != rs->getIndex1()[index1_it][0]) || (index1[index1_it][1] != rs->getIndex1()[index1_it][1])) || (index1[index1_it][2] != rs->getIndex1()[index1_it][2])) break;
    };
    if(index1_it == index1_length){
        unsigned int index2_it = 0;
        for(;index2_it < index2_length; index2_it++){
            if(((index2[index2_it][0] != rs->getIndex2()[index2_it][0]) || (index2[index2_it][1] != rs->getIndex2()[index2_it][1])) || (index2[index2_it][2] != rs->getIndex2()[index2_it][2])) break;
        };
        if(index2_it == index2_length) return true;
    };
    return false;
};
