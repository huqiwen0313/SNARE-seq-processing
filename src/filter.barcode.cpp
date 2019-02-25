// Prefiltering reads in fastq files with low count cells - SNARE-seq format
// Qiwen Hu 2018
// Usage: ./filter.barcode.cpp fastq_1 fastq_2 fastq_3

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

using namespace std;
std::unordered_map<string, int> bcode_map;
int barcode_length = 12;

// count reads frequencies in each barcode
void barcode_stats(string file){
  string r1_file = file;
  string r1_line;
  ifstream r1;
  r1.open(r1_file.c_str());
  if(r1.is_open()){
    while(getline(r1, r1_line)){
      if(r1_line.compare(0,1,"@") == 0){
       getline(r1, r1_line);
       string barcode = r1_line.substr(0, barcode_length);
       auto pos = bcode_map.find(barcode);
       if(pos == bcode_map.end()){
          bcode_map.emplace(barcode, 1);
       }
       else{
         pos->second = pos->second + 1;
       }
      }
     }
    r1.close();
  }
}

int main(int argc, char **argv){
  // input options
  // fastq file contains barcode
  string r1_file = argv[1];
  
  // pair-end fastq files
  string r2_file = argv[2];
  string r3_file = argv[3];

  // only select cells with reads more than a cutoff
  int cell_cutoff = 8000;
  
  // count reads fall into each cell
  barcode_stats(r1_file);

 // output read counts for each barcode
  string out_file_name = r1_file + ".bacode.sum.txt";
  //string line;
  ofstream out_file(out_file_name);
  if(out_file.is_open()){
     for(auto pos = bcode_map.begin(); pos != bcode_map.end(); ++pos){
        out_file << pos->first << "\t" << pos->second << "\n";
     }
  }
  else cout << "Unable to open file";
  out_file.close();
 
 // filter fastq files with cells less than cutoff value
 string r1_line; string r2_line; string r3_line;
 string r1_out = r1_file; string r2_out = r2_file;
 string r3_out = r3_file;
 r1_out = r1_out.replace(r1_out.end()-6, r1_out.end(), "") + ".filter.fastq";
 r2_out = r2_out.replace(r2_out.end()-6, r2_out.end(), "") + ".filter.fastq";
 r3_out = r3_out.replace(r3_out.end()-6, r3_out.end(), "") + ".filter.fastq";

 ifstream r1; ifstream r2; ifstream r3;
 ofstream out_file_r1(r1_out);
 ofstream out_file_r2(r2_out);
 ofstream out_file_r3(r3_out); 

 r1.open(r1_file.c_str());
 r2.open(r2_file.c_str());
 r3.open(r3_file.c_str());

 //cout << "begin open " << r1_file << "\n"; 
 string r1_header; string r2_header; string r3_header; 
 if(r1.is_open()){
    while(getline(r1, r1_line)){
      getline(r2, r2_line); getline(r3, r3_line);
      if(r1_line.compare(0,1,"@") == 0){
       r1_header = r1_line; r2_header = r2_line; r3_header = r3_line;
       getline(r1, r1_line); getline(r2, r2_line); getline(r3, r3_line);
       
       // searching barcode harsh table
       string barcode = r1_line.substr(0, barcode_length);
       auto pos = bcode_map.find(barcode);
       if(pos->second >= cell_cutoff){
         string r1_seq = r1_line; 
         string r2_seq = r2_line;
         string r3_seq = r3_line;
         getline(r1, r1_line); getline(r2, r2_line); getline(r3, r3_line);
         getline(r1, r1_line); getline(r2, r2_line); getline(r3, r3_line);
	 string r1_qual = r1_line;
	 string r2_qual = r2_line;
	 string r3_qual = r3_line;

	 out_file_r1 << r1_header << "\n";
	 out_file_r1 << r1_seq << "\n" << "+\n";
	 out_file_r1 << r1_qual << "\n";
	 
	 out_file_r2 << r2_header << "\n";
	 out_file_r2 << r2_seq << "\n" << "+\n";
	 out_file_r2 << r2_qual << "\n";
	 
	 out_file_r3 << r3_header << "\n";
         out_file_r3 << r3_seq << "\n" << "+\n";
	 out_file_r3 << r3_qual << "\n";
        }
       }
     }
     r1.close();
     r2.close();
     r3.close();
  }
  else cout << "Unable to open file:" << r1_file << "\n"; 
  //}

  return 0;
}
	
	
		
