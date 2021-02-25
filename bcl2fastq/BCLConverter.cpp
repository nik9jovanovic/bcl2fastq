//
//  BCLConverter.cpp
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 8/7/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#include "BCLConverter.h"

BCLConverter::BCLConverter(){
    char cCurrentPath[FILENAME_MAX];
    if(!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))) project_path = "";
    else project_path = string(cCurrentPath);
    
    read_parts = new ReadStructure* [4];
    for(int read_part = 0; read_part < 4; read_part++) read_parts[read_part] = NULL;
    
    sem = new Semaphore(thread::hardware_concurrency());
    
    cbcl_names = vector<string>();
    
    number_of_bytes = 0;
    index1_pos = 4;
    index2_pos = 4;
    read1_pos = 4;
    read2_pos = 4;
};

BCLConverter::~BCLConverter(){
    delete read_parts;
    delete sem;
    cbcl_names.clear();
    delete samples;
};

void BCLConverter::setProjectPath(const string path){
    project_path = path;
};

void BCLConverter::setThreads(const unsigned int threads){
    sem->setCount(threads);
};

void BCLConverter::setNumberOfBytes(const unsigned int bytes){
    number_of_bytes = bytes;
};

bitset<3>* BCLConverter::convertIndex(const string index, unsigned int length) const{
    bitset<3>* idx = new bitset<3> [length];
    for(int base = 0; base < length; base++){
        switch(index[base]){
            case('A'):{
                idx[base] = bitset<3> ("000");
                continue;
            };
            case('C'):{
                idx[base] = bitset<3> ("001");
                continue;
            };
            case('G'):{
                idx[base] = bitset<3> ("010");
                continue;
            };
            case('T'):{
                idx[base] = bitset<3> ("011");
                continue;
            };
            default:{
                idx[base] = bitset<3> ("100");
                continue;
            };
        };
    }
    return idx;
};

unsigned int BCLConverter::checkAndAddBCLName(const string cbcl_name){
    for(std::vector<string>::iterator cbcl = cbcl_names.begin(); cbcl != cbcl_names.end(); ++cbcl){
        if((*cbcl).compare(cbcl_name) == 0) return 1;
    }
    cbcl_names.push_back(cbcl_name);
    bcls_per_cycle++;
    return 0;
};

unsigned char BCLConverter::decodeIndexBase(const unsigned char byte) const{
    if((byte & 0b1100) == 0) return 'N';
    switch(byte & 0b11){
        case(0): return 'A';
        case(1): return 'C';
        case(2): return 'G';
        case(3): return 'T';
        default: return 'N';
    };
};

DecodedBase BCLConverter::decodeSequenceBase(const unsigned char byte, const unsigned int* bins) const{
    if((byte & 0b1100) == 0) return DecodedBase('N', 35);
    unsigned char quality = 33 + bins[(byte >> 2) & 0b11];
    switch(byte & 0b11){
        case(0): return DecodedBase('A', quality);
        case(1): return DecodedBase('C', quality);
        case(2): return DecodedBase('G', quality);
        case(3): return DecodedBase('T', quality);
        default: return DecodedBase('N', 35);
    };
};

//messy - change for something neater maybe
unsigned int BCLConverter::parseRunInfo(const string runinfo_path){
    string ri_path = runinfo_path;
    if(ri_path.empty()) ri_path = project_path + "/RunInfo.xml";
    ifstream runinfo(ri_path);
    if(runinfo.fail()){
        cerr << "RunInfo.xml not found.";
        exit(1);
    };
    bool reading_reads = false;
    unsigned int cycles_so_far = 0;
    unsigned char read_part = 0;
    for(string line; getline(runinfo, line); ){
        size_t leading_tabs = line.find_first_of('<');
        if (run_id.empty()){
            if (line.find("Run") != string::npos){
                size_t start = line.find("Number=", (leading_tabs + 1));
                if (start == string::npos) continue;
                size_t end = line.find('>', start);
                if (end == string::npos) continue;
                run_id = line.substr((start + 8), (end - 1) - (start + 8));
                continue;
            };
        };
        if(instrument.empty()){
            if (line.find("Instrument") != string::npos){
                size_t start = line.find('>', (leading_tabs + 1));
                if (start == string::npos) continue;
                size_t end = line.find('<', (start + 1));
                if (end == string::npos) continue;
                instrument = line.substr((start + 1), (end - (start + 1)));
                continue;
            };
        };
        if(flowcell_id.empty()){
            if (line.find("Flowcell") != string::npos){
                size_t start = line.find('>', (leading_tabs + 1));
                if (start == string::npos) continue;
                size_t end = line.find('<', (start + 1));
                if (end == string::npos) continue;
                flowcell_id = line.substr( (start + 1), (end - (start + 1)));
                continue;
            };
        };
        if(line.find("Read ", (leading_tabs + 1)) != string::npos){
            reading_reads = true;
            size_t num_cycles_pos = line.find("NumCycles");
            unsigned int numcycles = 0;
            if(num_cycles_pos != string::npos) numcycles = stoi(line.substr((num_cycles_pos + 11), line.find('\"', (num_cycles_pos + 11)) - (num_cycles_pos + 11)));
            size_t indexedread = line.find("IsIndexedRead");
            bool idx = (indexedread != string::npos && line[indexedread + 15] == 'Y') ? true : false;
            if(numcycles == 0) continue;
            if (read_part == 4){
                cerr << "The RunInfo.xml contains more than 4 read parts - not currently supported.";
                exit(1);
            }
            else{
                CBCLReader** tmp_readers = new CBCLReader* [numcycles];
                for (int cycle = 0; cycle < numcycles; cycle++) tmp_readers[cycle] = new CBCLReader();
                read_parts[read_part] = new ReadStructure(idx, numcycles, cycles_so_far, tmp_readers);
                cycles_so_far += numcycles;
                if(idx){
                    if(index1_pos == 4) index1_pos = read_part;
                    else index2_pos = read_part;
                }
                else{
                    if(read1_pos == 4) read1_pos = read_part;
                    else read2_pos = read_part;
                }
                read_part++;
                continue;
            };
        };
        if(reading_reads) break;
    };
    runinfo.close();
    for(int rg=0; rg < 4; rg++){
        if(read_parts[rg]) num_cycles += read_parts[rg]->getLength();
    };
    return 0;
};

unsigned int BCLConverter::parseSampleSheet(string samplesheet_path){
    string ss_path = samplesheet_path;
    if(ss_path.empty()) ss_path = project_path + "/SampleSheet.csv";
    ifstream ss(ss_path);
    if(ss.fail()){
        return 1;
    };
    bool reading_settings = false;
    bool num_samples_counting = false;
    bool reading_index_header = false;
    int index1_comma = -1;
    int index2_comma = -1;
    int sampleid_comma = -1;
    int samplename_comma = -1;
    unsigned int index1_length = (index1_pos < 4) ? read_parts[index1_pos]->getLength() : 0;
    unsigned int index2_length = (index2_pos < 4) ? read_parts[index2_pos]->getLength() : 0;
    streampos datastart = 0;
    for(string line; getline(ss, line); ){
        if(line.find("[Settings]") != string::npos){
            reading_settings = true;
            continue;
        };
        if(line.find("[Data]") != string::npos){
            reading_settings = false;
            reading_index_header = true;
            num_samples_counting = true;
            datastart = ss.tellg();
            continue;
        };
        if(num_samples_counting){
            if(reading_index_header){
                reading_index_header = false;
                continue;
            };
            if(line.find(',') != string::npos) num_samples++;
            continue;
        };
        if(reading_settings){
            stringstream ss(line);
            string token;
            while (getline(ss,token, ',')){
                if(!token.compare("Read1UMILength")){
                    getline(ss,token, ',');
                    umi1_length = stoi(token);
                    break;
                };
                if(!token.compare("Read2UMILength")){
                    getline(ss,token, ',');
                    umi2_length = stoi(token);
                    break;
                };
            };
        };
    };
    ss.clear();
    if(datastart != 0) ss.seekg(datastart);
    samples = new SampleSheetEntry* [num_samples];
    reading_index_header = true;
    bool reading_indexes = false;
    unsigned int index_num = 0;
    string token = "";
    for(string line; getline(ss, line); ){
        stringstream ss (line);
        unsigned int num_parts = 0;
        if(reading_index_header){
            while (getline(ss,token, ',')){
                num_parts++;
                if(!token.compare("index")){
                    index1_comma = num_parts;
                    continue;
                };
                if(!token.compare("index2")){
                    index2_comma = num_parts;
                    continue;
                };
                if(!token.compare("Sample_ID")){
                    sampleid_comma = num_parts;
                    continue;
                };
                if(!token.compare("Sample_Name")){
                    samplename_comma = num_parts;
                    continue;
                };
            };
            reading_index_header = false;
            reading_indexes = true;
            continue;
        };
        if(reading_indexes){
            string sampleid = "";
            string samplename = "";
            bitset<3>* index1 = NULL;
            bitset<3>* index2 = NULL;
            if(sampleid_comma == -1){
                cerr << "No Sample_ID column found in the samplesheet.";
                exit(1);
            };
            while (getline(ss,token, ',')){
                num_parts++;
                if(num_parts == sampleid_comma){
                    sampleid = token;
                    continue;
                };
                if(num_parts == samplename_comma){
                    samplename = token;
                    continue;
                };
                if(num_parts == index1_comma){
                    index1 = convertIndex(token, index1_length);
                    continue;
                };
                if(num_parts == index2_comma){
                    index2 = convertIndex(token, index2_length);
                    continue;
                };
            };
            samples[index_num] = new SampleSheetEntry(sampleid, samplename, index1, index2);
            for(int sample = 0; sample < index_num; sample++){
                if(samples[index_num]->compareIndexes(samples[sample], index1_length, index2_length)){
                    cerr << "Multiple of the same index entries found in the samplesheet.";
                    exit(1);
                };
            };
            index_num++;
        };
    };
    return 0;
};

unsigned int BCLConverter::checkBCLDirectory(){
    string basecalls_path = project_path + "/Data/Intensities/BaseCalls";
    DIR* basecalls_dir = opendir(basecalls_path.c_str());
    if (basecalls_dir == NULL){
        cerr << "Could not open " + project_path + "/Data/Intensities/BaseCalls directory.";
        exit(1);
    };
    dirent* basecalls_ent;
    while ((basecalls_ent = readdir (basecalls_dir)) != NULL) {
        if (string(basecalls_ent->d_name).find("L00") != string::npos){
            DIR* lane_dir = opendir((basecalls_path + "/" + string(basecalls_ent->d_name)).c_str());
            if(lane_dir == NULL) continue;
            num_lanes++;
            unsigned int read_cycles = 0;
            dirent* lane_ent;
            while ((lane_ent = readdir (lane_dir)) != NULL) {
                if (string(lane_ent->d_name).find("C") != string::npos){
                    read_cycles++;
                    //fetch only directories
                    DIR* cycle_dir = opendir((basecalls_path + "/" + string(basecalls_ent->d_name) + "/" + string(lane_ent->d_name)).c_str());
                    if(cycle_dir == NULL) continue;
                    dirent* cycle_ent;
                    while ((cycle_ent = readdir (cycle_dir)) != NULL) {
                        if (string(cycle_ent->d_name).find(".cbcl") != string::npos){
                            checkAndAddBCLName(string(cycle_ent->d_name));
                        };
                    };
                    closedir(cycle_dir);
                };
            };
            closedir(lane_dir);
            if(num_cycles != read_cycles) cout << "Incompatible number of cycles in lane " + string(basecalls_ent->d_name) + " with reported in RunInfo.xml.";
        };
    };
    std::sort (cbcl_names.begin(), cbcl_names.end());
    closedir (basecalls_dir);
    return 0;
};

unsigned char* BCLConverter::readLocs(const string path, const unsigned int startpos){
    FILE* locs_file = fopen(path.c_str(), "r");
    if(locs_file == NULL) return NULL;
    fseek(locs_file, 8, SEEK_SET);
    unsigned int num_clusters = fgetc(locs_file) + (fgetc(locs_file) << 8) + (fgetc(locs_file) << 16) + (fgetc(locs_file) << 24);
    
    if(startpos > num_clusters) return NULL;
    fseek(locs_file, 12 + (8 * startpos), SEEK_SET);
    
    unsigned int num_reads = 2 * number_of_bytes;
    if(num_reads == 0) num_reads = num_clusters;
    if((startpos + num_reads) > num_clusters) num_reads = num_clusters - startpos;
    
    unsigned char* locsBuffer = new unsigned char [num_reads * 8];
    size_t result = fread (locsBuffer, 1, num_reads * 8, locs_file);
    if (result != (num_reads * 8)) return NULL;
    
    return locsBuffer;
};

unsigned char* BCLConverter::readFilter(const string path, const unsigned int startpos){
    FILE* filter_file = fopen(path.c_str(), "r");
    if(filter_file == NULL) return NULL;
    fseek(filter_file, 8, SEEK_SET);
    unsigned int num_clusters = fgetc(filter_file) + (fgetc(filter_file) << 8) + (fgetc(filter_file) << 16) + (fgetc(filter_file) << 24);
    
    if(startpos > num_clusters) return NULL;
    fseek(filter_file, 12 + startpos, SEEK_SET);
    
    unsigned int num_reads = 2 * number_of_bytes;
    if(num_reads == 0) num_reads = num_clusters;
    if((startpos + num_reads) > num_clusters) num_reads = num_clusters - startpos;
    
    unsigned char* filterBuffer = new unsigned char [num_reads];
    size_t result = fread (filterBuffer, 1, num_reads, filter_file);
    if (result != num_reads) return NULL;
    
    return filterBuffer;;
};

unsigned int BCLConverter::demuxData(const unsigned char barcode_misses){
    
    unsigned int reads_per_thread = number_of_bytes / sem->getCount();
    vector<thread> running_threads;
    if(number_of_bytes < sem->getCount()){
        reads_per_thread = 2;
    };
    
    unsigned int reads_so_far = 0;
    unsigned int reads_to_process = 0;
    unsigned int threads = 0;
    Semaphore* s = new Semaphore(1);
    while(reads_so_far < (2 * number_of_bytes)){
        if(((2 * number_of_bytes) - reads_so_far) > reads_per_thread) reads_to_process = reads_per_thread;
        else reads_to_process = (2 * number_of_bytes) - reads_so_far;
        sem->wait();
        running_threads.push_back(thread(&BCLConverter::demuxPart, this, reads_so_far, reads_to_process, barcode_misses, s));
        threads++;
        reads_so_far += reads_to_process;
    };
    for(int thr = 0; thr < threads; thr++) running_threads.at(thr).join();
    return 0;
};

unsigned int BCLConverter::demuxPart(const unsigned int startpos, const unsigned int num_reads, const unsigned char barcode_misses, Semaphore* s){
    ReadStructure* index1 = (index1_pos < 4) ? read_parts[index1_pos] : NULL;
    ReadStructure* index2 = (index2_pos < 4) ? read_parts[index2_pos] : NULL;
    for(unsigned int pos = startpos; pos < startpos + num_reads; pos++){
        unsigned int sample = 0;
        for(; sample < num_samples; sample++){
            bitset<3>* idx1 = samples[sample]->getIndex1();
            CBCLReader** idx1_cbcls = index1->getCBCLs();
            unsigned int i1_misses = 0;
            unsigned int i2_misses = 0;
            unsigned int base = 0;
            unsigned char tmp = 0;
            for(; base < index1->getLength(); base++){
                if(idx1[base][2] == 1) continue;
                if(idx1_cbcls[base]->getBuffer()){
                    tmp = idx1_cbcls[base]->getBuffer()[pos / 2] >> ( (pos % 2) * 4);
                    if((tmp & 0b1100) == 0){
                        if(++i1_misses > barcode_misses) break;
                        continue;
                    }
                    if((idx1[base][0] == (tmp & 1)) && (idx1[base][1] == ((tmp >> 1) & 1))) continue;
                    if(++i1_misses > barcode_misses) break;
                };
            };
            if(i1_misses > barcode_misses) continue;
            if(base == index1->getLength() && index2){
                base = 0;
                bitset<3>* idx2 = samples[sample]->getIndex2();
                CBCLReader** idx2_cbcls = index2->getCBCLs();
                for(; base < index2->getLength(); base++){
                    if(idx2[base][2] == 1) continue;
                    if(idx2_cbcls[base]->getBuffer()){
                        tmp = idx2_cbcls[base]->getBuffer()[pos / 2] >> ((pos % 2) * 4);
                        if((tmp & 0b1100) == 0){
                            if(++i2_misses > barcode_misses) break;
                            continue;
                        }
                        if((idx2[base][0] == (tmp & 1)) && (idx2[base][1] == ((tmp >> 1) & 1))) continue;
                        if(++i2_misses > barcode_misses) break;
                    };
                    if(++i2_misses > barcode_misses) break;
                };
                if(i2_misses > barcode_misses) continue;
                if(base == index2->getLength()){
                    s->wait();
                    demuxed_indexes.at(sample + 1).push_back(pos);
                    s->notify();
                    break;
                };
            };
            if(base == index1->getLength()){
                s->wait();
                demuxed_indexes.at(sample + 1).push_back(pos);
                s->notify();
                break;
            };
        };
        if(sample == num_samples){
            s->wait();
            demuxed_indexes.at(0).push_back(pos);
            s->notify();
        };
    };
    sem->notify();
    return 0;
};

unsigned int BCLConverter::writeData(const unsigned int lane, const unsigned int tile_num, const bool pair2, const unsigned char* locs_buffer, const unsigned char* filter_buffer){
    vector<thread> running_threads;
    for(int sample = 0; sample < num_samples + 1; sample++){
        if(!demuxed_indexes.at(sample).empty()){
            std::sort (demuxed_indexes.at(sample).begin(), demuxed_indexes.at(sample).end());
            string header = "@" + instrument + ":" + run_id + ":" + flowcell_id + ":" + to_string(lane) + ":" + to_string(tile_num);
            string samplename = (sample == 0) ? "Undetermined" : (samples[sample - 1]->getSampleName().empty() ? samples[sample - 1]->getSampleID() : samples[sample - 1]->getSampleName());
            string path = project_path + "/Data/Intensities/BaseCalls/" + samplename + "_S" + to_string(sample) + "_L00" + to_string(lane) + "_R" + to_string(pair2 + 1) + "_001.fastq";
            sem->wait();
            running_threads.push_back(thread(&BCLConverter::writePart, this, sample, path, pair2, header, locs_buffer, filter_buffer));
        };
    };
    for(std::vector<thread>::iterator thread = running_threads.begin(); thread != running_threads.end(); ++thread) thread->join();
    running_threads.clear();
    return 0;
};

unsigned int BCLConverter::writePart(const unsigned int number, const string path, const bool pair2, const string header, const unsigned char* locs_buffer, const unsigned char* filter_buffer){
    ofstream fastq;

    unsigned int readpos = pair2 ? read2_pos : read1_pos;
    if(readpos == 4) return 1;
    for(std::vector<unsigned int>::iterator read = demuxed_indexes.at(number).begin(); read != demuxed_indexes.at(number).end(); ++read){
        if(filter_buffer){
            if(filter_buffer[*read] == 0){
                continue;
            };
        };
        if(!fastq.is_open()){
            fastq.open(path, std::ofstream::out | std::ofstream::app | std::ofstream::ate);
            if(!fastq.is_open()) return 1;
        };
        string index_sequence = "";
        if(index1_pos < 4){
            for(unsigned int base = 0; base < read_parts[index1_pos]->getLength(); base++){
                if(read_parts[index1_pos]->getCBCLs()[base]->getBuffer()){
                    index_sequence += decodeIndexBase( ((*read) % 2) ? read_parts[index1_pos]->getCBCLs()[base]->getBuffer()[*read / 2] >> 4 : read_parts[index1_pos]->getCBCLs()[base]->getBuffer()[*read / 2]);
                }
                else index_sequence += 'N';
            };
        };
        if(index2_pos < 4){
            index_sequence += '+';
            for(unsigned int base = 0; base < read_parts[index2_pos]->getLength(); base++){
                if(read_parts[index2_pos]->getCBCLs()[base]->getBuffer()){
                    index_sequence += decodeIndexBase( ((*read) % 2) ? read_parts[index2_pos]->getCBCLs()[base]->getBuffer()[*read / 2] >> 4 : read_parts[index2_pos]->getCBCLs()[base]->getBuffer()[*read / 2]);
                }
                else index_sequence += 'N';
            };
        };
        string sequence = "";
        string quality = "";
        for(int base = 0; base < read_parts[readpos]->getLength(); base++){
            if(read_parts[readpos]->getCBCLs()[base]->getBuffer()){
                DecodedBase tmp_base = decodeSequenceBase(((*read) % 2) ? read_parts[readpos]->getCBCLs()[base]->getBuffer()[*read / 2] >> 4 : read_parts[readpos]->getCBCLs()[base]->getBuffer()[*read / 2], read_parts[readpos]->getCBCLs()[base]->getBins());
                sequence += tmp_base.base;
                quality += tmp_base.qual;
            }
            else{
                sequence += 'N';
                quality += 35;
            };
        };
        string tmpheader = header;
        if(locs_buffer){
            float tmp = static_cast<int>(round((*((float*)(locs_buffer + ((*read) * 8))) * 10 + 1000)));
            tmpheader += ":" + to_string(static_cast<int>(tmp));
            tmp = static_cast<int>(round((*((float*)(locs_buffer + ((*read) * 8) + 4)) * 10 + 1000)));
            if(tmp != 1000){
                tmpheader += ":" + to_string(static_cast<int>(tmp));
            }
            else{
                tmpheader += ":" + to_string(*read);
            };
        };
        fastq << tmpheader << ' ' << to_string(pair2 + 1) << ":N:0:" << index_sequence << '\n' << sequence << "\n+\n" << quality << '\n';
    };
    if(fastq.is_open()) fastq.close();
    sem->notify();
    return 0;
};

unsigned int BCLConverter::convert(const unsigned char barcode_misses, const bool read_filters, const bool read_positions){
    TileInfo* tilesToProcess;
    unsigned char* locs_buffer = NULL;
    unsigned char* filter_buffer = NULL;
    auto start_time = chrono::high_resolution_clock::now();
    demuxed_indexes = vector<vector<unsigned int>>(num_samples + 1, vector<unsigned int>());
    for(unsigned int lane = 1; lane < num_lanes + 1; lane++){
        for(unsigned int bcl = 0; bcl < bcls_per_cycle; bcl++){
            read_parts[read1_pos]->readHeaders(project_path, cbcl_names[bcl], lane, sem);
            unsigned int pos = 0;
            for(;!read_parts[read1_pos]->getCBCLs()[pos]->getTiles();pos++);
            tilesToProcess = read_parts[read1_pos]->getCBCLs()[pos]->getTiles();
            unsigned int num_tiles = read_parts[read1_pos]->getCBCLs()[pos]->getNumTiles();
            unsigned int position = 0;
            for(unsigned int tile = 0; tile < num_tiles; tile++){
                auto tile_start_time = chrono::high_resolution_clock::now();
                cout << "Processing tile " << tilesToProcess[tile].tile_num << " of lane " << to_string(lane) << " started. " << std::chrono::duration_cast<std::chrono::milliseconds>(tile_start_time - start_time).count() << "\n";
                unsigned int reads_processed = 0;
                if(!number_of_bytes) number_of_bytes = tilesToProcess[tile].num_clusters / 2;
                unsigned int bytes_to_process = number_of_bytes;
                while(reads_processed < tilesToProcess[tile].num_clusters){
                    if(((tilesToProcess[tile].num_clusters - reads_processed) / 2) > bytes_to_process) bytes_to_process = (tilesToProcess[tile].num_clusters - reads_processed) / 2;
                    if(index1_pos < 4)read_parts[index1_pos]->readCBCLs(project_path + "/Data/Intensities/BaseCalls/L00" + to_string(lane), cbcl_names[bcl], tile, position, bytes_to_process, sem);
                    if(index2_pos < 4)read_parts[index2_pos]->readCBCLs(project_path + "/Data/Intensities/BaseCalls/L00" + to_string(lane), cbcl_names[bcl], tile, position, bytes_to_process, sem);
                    auto index_read_finish_time = chrono::high_resolution_clock::now();
                    cout << "Reading CBCLs for Indexes for tile " << tilesToProcess[tile].tile_num << " of lane " << to_string(lane) << " done. " << std::chrono::duration_cast<std::chrono::milliseconds>(index_read_finish_time - start_time).count() << "\n";
                    if(index1_pos == 4 && index2_pos == 4){
                        for (unsigned int adder = 0; adder < bytes_to_process * 2; adder++) demuxed_indexes.at(0).push_back(adder);
                    }
                    else demuxData(barcode_misses);
                    auto demux_finish_time = chrono::high_resolution_clock::now();
                    cout << "Demultiplexing data for tile " << tilesToProcess[tile].tile_num << " of lane " << to_string(lane) << " done. " << std::chrono::duration_cast<std::chrono::milliseconds>(demux_finish_time - start_time).count() << "\n";
                    if(read_filters) filter_buffer = this->readFilter(project_path + "/Data/Intensities/BaseCalls/L00" + to_string(lane) + "/s_" + to_string(lane) + "_" + to_string(tilesToProcess[tile].tile_num) + ".filter", position);
                    if(read_positions) locs_buffer = this->readLocs(project_path + "/Data/Intensities/s.locs", position);
                    if(read1_pos < 4){
                        read_parts[read1_pos]->readCBCLs(project_path + "/Data/Intensities/BaseCalls/L00" + to_string(lane), cbcl_names[bcl], tile, position, bytes_to_process, sem);
                        auto read1_read_finish_time = chrono::high_resolution_clock::now();
                        cout << "Reading CBCLs for Read1 for tile " << tilesToProcess[tile].tile_num << " of lane " << to_string(lane) << " done. " << std::chrono::duration_cast<std::chrono::milliseconds>(read1_read_finish_time - start_time).count() << "\n";
                        writeData(lane, tilesToProcess[tile].tile_num, false, locs_buffer, filter_buffer);
                        auto read1_writing_finish_time = chrono::high_resolution_clock::now();
                        cout << "Writing FASTQ files for Read1 for tile " << tilesToProcess[tile].tile_num << " of lane " << to_string(lane) << "done. " << std::chrono::duration_cast<std::chrono::milliseconds>(read1_writing_finish_time - start_time).count() << "\n";
                        read_parts[read1_pos]->clear();
                    };
                    if(read2_pos < 4){
                        read_parts[read2_pos]->readCBCLs(project_path + "/Data/Intensities/BaseCalls/L00" + to_string(lane), cbcl_names[bcl], tile, position, bytes_to_process, sem);
                        auto read2_read_finish_time = chrono::high_resolution_clock::now();
                        cout << "Reading CBCLs for Read2 for tile " << tilesToProcess[tile].tile_num << " of lane " << to_string(lane) << " done. " << std::chrono::duration_cast<std::chrono::milliseconds>(read2_read_finish_time - start_time).count() << "\n";
                        writeData(lane, tilesToProcess[tile].tile_num, true, locs_buffer, filter_buffer);
                        auto read2_writing_finish_time = chrono::high_resolution_clock::now();
                        cout << "Writing FASTQ files for Read2 for tile " << tilesToProcess[tile].tile_num << " of lane " << to_string(lane) << " done. " << std::chrono::duration_cast<std::chrono::milliseconds>(read2_writing_finish_time - start_time).count() << "\n";
                        read_parts[read2_pos]->clear();
                    };
                    for(unsigned int demux_iter = 0; demux_iter < num_samples + 1; demux_iter++){
                        if(!demuxed_indexes.at(demux_iter).empty()){
                            demuxed_indexes[demux_iter].clear();
                        };
                    };
                    if(index1_pos < 4) read_parts[index1_pos]->clear();
                    if(index2_pos < 4) read_parts[index2_pos]->clear();
                    position += bytes_to_process;
                    reads_processed += bytes_to_process * 2;
                };
            };
        };
    };
    auto conversion_finish_time = chrono::high_resolution_clock::now();
    cout << "Conversion done. " << std::chrono::duration_cast<std::chrono::milliseconds>(conversion_finish_time - start_time).count() << "\n";
    return 0;
};
