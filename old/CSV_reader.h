#ifndef __CSV_READER_H__
#define __CSV_READER_H__

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <list>
using namespace std;

class CSV_reader{
 public:
  CSV_reader();
  ~CSV_reader();
  double CSV_find(int Row, int Col, string Spec_CSV_file); //string CSV_file_name);
};

#endif
