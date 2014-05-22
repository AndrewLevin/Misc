//this macro finds trigger menus that have a given trigger in them
//you can use a trigger name with or without a version number

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//std::string trigger_name =  "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL";
std::string trigger_name =  "HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL";

std::vector<std::string> trigger_menus;

int main()
{
  std::vector<std::string> lines;

  std::string line;
  std::fstream f("log_query_confdb_2011.dat");
  while(!f.eof()){
    f >> line;
    lines.push_back(line);
    if(line.find(trigger_name) != std::string::npos){
      int line_number = lines.size() - 2;
      while(lines[line_number].find("cdaq") == std::string::npos) line_number--;
      
      //std::cout << "lines[line_number] = " << lines[line_number] << std::endl
      //return 1;
      if(trigger_menus.size()== 0)
	trigger_menus.push_back(lines[line_number]);
      else if(lines[line_number] != trigger_menus[trigger_menus.size()-1])
	trigger_menus.push_back(lines[line_number]);
      
    }
  }

  for(int i = 0; i < trigger_menus.size(); i++)
    std::cout << "trigger_menus[" << i << "] = " << trigger_menus[i] << std::endl;



}
