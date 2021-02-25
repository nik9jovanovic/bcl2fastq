//
//  ReadStructure.cpp
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 8/15/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#include "ReadStructure.h"

ReadStructure::ReadStructure(){
    is_index = false;
    length = 0;
    start = 0;
    cbcls = NULL;
};

ReadStructure::ReadStructure(const bool idx, const unsigned int l, const unsigned int s, CBCLReader** c){
    is_index = idx;
    length = l;
    start = s;
    cbcls = c;
};

ReadStructure::ReadStructure(const ReadStructure &rs){
    is_index = rs.getIndex();
    length = rs.getLength();
    start = rs.getStart();
    cbcls = rs.getCBCLs();
};

ReadStructure::~ReadStructure(){
    delete [] cbcls;
};

bool ReadStructure::getIndex() const{
    return is_index;
};

unsigned int ReadStructure::getLength() const{
    return length;
};

unsigned int ReadStructure::getStart() const{
    return start;
};

CBCLReader** ReadStructure::getCBCLs() const{
    return cbcls;
};

void ReadStructure::clear(){
    for(int clr = 0; clr < length; clr++) cbcls[clr]->clearBuffer();
};

void ReadStructure::readHeaders(const string project_path, const string bclname, const unsigned int lane, Semaphore* sem){
    thread* running_threads = new thread [length];
    for(int thr = start; thr < (start + length); thr++){
        sem->wait();
        running_threads[thr - start] = thread(&CBCLReader::parseHeader, cbcls[thr - start], project_path + "/Data/Intensities/BaseCalls/L00" + to_string(lane) + "/C" + to_string(thr + 1) + ".1/" + bclname, sem);
    }
    for(int thr = 0; thr < length; thr++) running_threads[thr].join();
    delete [] running_threads;
};

void ReadStructure::readCBCLs(const string lane_directory, const string bclname, const unsigned int tile, const unsigned int startpos, const unsigned int num_reads, Semaphore* sem){
    thread* running_threads = new thread [length];
    for(int thr = start; thr < (start + length); thr++){
        sem->wait();
        running_threads[thr - start] = thread(&CBCLReader::readCBCL, cbcls[thr - start], lane_directory + "/C" + to_string(thr + 1) + ".1/" + bclname, tile, startpos, num_reads, sem);
    }
    for(int thr = 0; thr < length; thr++) running_threads[thr].join();
    delete [] running_threads;
};
