/**
 * Copyright (c) 2015 Jochen Singer
 *
 * This is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or any later version.
 *
 * This software is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this software. If not, see <http://www.gnu.org/licenses/>.
 */

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>


int main(int argc, const char* argv[]){

    std::vector<std::string> fileNames;
    fileNames.push_back(argv[1]);
    fileNames.push_back(argv[2]);

    std::string leftOnly = argv[3];
    std::string rightOnly = argv[4];
    std::string both = argv[5];

    std::vector<std::string> headers;
    headers.push_back("");
    headers.push_back("");

    std::vector<std::map<std::pair<std::string, unsigned>, std::string> > vcfMap;
    vcfMap.resize(2);
    std::string currentLine;

    for (unsigned i = 0; i < fileNames.size(); ++i)
    {
        std::ifstream vcfFileIn (fileNames[i]);
        if (vcfFileIn.is_open())
        {
            while (getline(vcfFileIn, currentLine))
            {
                if (currentLine[0] != '#')
                {
                    std::size_t posChr = currentLine.find("\t", 0);
                    std::size_t posPos = currentLine.find("\t", posChr + 1);
                    vcfMap[i].insert(std::pair<std::pair<std::string, unsigned>, std::string>(std::pair<std::string, unsigned>(currentLine.substr(0, posChr),std::stoi(currentLine.substr(posChr + 1, posPos))), currentLine));
                }
                else
                    headers[i]+=currentLine + "\n";
            }
        }
        else 
        {
            std::cout << "Unable to open file";
        }
    }

    typedef std::map<std::pair<std::string, unsigned>, std::string>::iterator TIter;

    std::ofstream vcfFileOutBoth(both, std::ofstream::out | std::ofstream::app);
    vcfFileOutBoth << headers[0];

    std::ofstream vcfFileOutLeftOnly(leftOnly, std::ofstream::out | std::ofstream::app);
    vcfFileOutLeftOnly << headers[0];

    std::ofstream vcfFileOutRightOnly(rightOnly, std::ofstream::out | std::ofstream::app);

    vcfFileOutRightOnly << headers[1];


    TIter compIt;
    for (TIter it = vcfMap[0].begin(); it != vcfMap[0].end(); ++it)
    {
        compIt = vcfMap[1].find(it->first);
        if (compIt != vcfMap[1].end())
            vcfFileOutBoth << it->second + "\n";
        else
            vcfFileOutLeftOnly << it->second + "\n";
    }
    vcfFileOutBoth.close();
    vcfFileOutLeftOnly.close();

    for (TIter it = vcfMap[1].begin(); it != vcfMap[1].end(); ++it)
    {
        compIt = vcfMap[0].find(it->first);
        if (compIt == vcfMap[0].end())
            vcfFileOutRightOnly << it->second + "\n";
    }
    vcfFileOutRightOnly.close();
}

