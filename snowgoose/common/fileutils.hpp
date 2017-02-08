#ifndef __FILEUTILS_HPP__
#define __FILEUTILS_HPP__

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <sstream>

#include "sgerrcheck.hpp"

namespace snowgoose {

    class FileUtils {
    public:

        /**
         * Remove files in a folder.
         * @param dir folder name
         * @param pattern string pattern to match against the file name
         */
        static void removeFiles(const std::string &dir, const std::string &pattern) {
            std::string cmd = "cd ";
            cmd += dir;
            cmd += "; rm -f ";
            cmd += pattern;
            int rc = system(cmd.c_str());
            if (rc != 0)
                SG_ERROR_REPORT("failed to remove files");
        }

        /**
         * Reading a string from a file
         * @param fname file name 
         * @param str resulting string
         */
        static void getStringFromFile(const char* fname, std::string& str) {
            std::ifstream is(fname);
            std::ostringstream os;
            if (is.is_open()) {
                /*/
                while (!is.eof()) {
                    std::string tmp;
                    is >> tmp;
                    std::cout << "Read " << tmp << "\n";
                    str += tmp;
                }
                 */
                os << is.rdbuf();
                str += os.str();
                is.close();
            } else {
                SG_ERROR_REPORT("Unable to open the input data file\n");
            }
        }

        /**
         * Creates a file with the given name and content
         * if file is opened updates its contents
         * 
         * @param fname file name
         * @param content file content 
         */
        static void updateFileWithContent(const char* fname, const char* content) {
            std::ofstream os(fname, std::ios_base::app);
            if (os.is_open()) {
                os << content;
                os.close();
            } else {
                SG_ERROR_REPORT("Unable to open the output file\n");
            }
        }
    };
}
#endif
