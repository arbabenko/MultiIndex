#include <string>
#include <vector>
#include "data_util.h"

using std::string;

struct RerankBytes {
  unsigned char bytes[8];
};

template<>
inline RerankBytes readFromFile(std::ifstream& input) {
  RerankBytes record;
  for(int q = 0; q < 8; ++q){
      input.read((char*)&(record.bytes[q]), sizeof(unsigned char));
  }
  return record;
}

int main() {
  string rerank_bytes_filename = "/sata/ResearchData/rerBytesGlobal_1000000000.dat";
  string cell_starts_filename = "/sata/ResearchData/cellStarts_1000000000.dat";
  string quantizations_filename = "/sata/ResearchData/cq_1000000000.dat";
  string index_filename = "/sata/ResearchData/indexGlobal_1000000000.dat";
  vector<RerankBytes> rerank_info;
  fillVector<RerankBytes>(rerank_bytes_filename, 1000000000, 32, &rerank_info);
  vector<int> cell_starts;
  fillVector<int>(cell_starts_filename, 16384 * 16384, 32, &cell_starts);
  vector<short> quantizations;
  fillVector<short>(quantizations_filename, 2000000000, 32, &quantizations);
  vector<int> points_written(16384 * 16384);
  vector<RerankADC8> index(1000000000);
  for(int pid = 0; pid < 1000000000; ++pid) {
      if(pid % 1000000 == 0) {
        std::cout << pid << std::endl;
      }
      int cell_id = (int)quantizations[pid] * 16384 + quantizations[pid + 1000000000];
      index[cell_starts[cell_id] + points_written[cell_id]].pid = pid;
      memcpy(&(index[cell_starts[cell_id] + points_written[cell_id]].quantizations[0]), &(rerank_info[pid].bytes[0]), 8);
      //for(int q = 0; q < 8; ++q) {
      //    index[cell_starts[cell_id] + points_written[cell_id]].quantizations[q] = rerank_info[pid].bytes[q];
      //}
      points_written[cell_id] += 1;
  }
  std::ofstream out(index_filename.c_str(), std::ios::out | std::ios::binary);
  for(int pid = 0; pid < 1000000000; ++pid) {
    out.write((char*)&(index[pid].pid), sizeof(int));
    for(int q = 0; q < 8; ++q) {
      out.write((char*)&(index[pid].quantizations[q]), sizeof(unsigned char));
    }
  }
  out.close();
  return 0;
}
