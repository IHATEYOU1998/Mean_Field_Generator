
#include "CSV_reader.h"
 

CSV_reader::CSV_reader(){}
CSV_reader::~CSV_reader(){}
double CSV_reader::CSV_find(int Row, int Col, string CSV_file_name)
  {

    std::ifstream  data(CSV_file_name);
    std::string line;
    std::vector<std::vector<std::string>> parsedCsv;
    while(std::getline(data,line))
    {
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<std::string> parsedRow;
        while(std::getline(lineStream,cell,','))
        {
            parsedRow.push_back(cell);
        }

        parsedCsv.push_back(parsedRow);
    }
    return atof(parsedCsv[Row][Col].c_str());
}
