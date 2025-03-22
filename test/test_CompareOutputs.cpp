#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>

bool compareFiles(const std::string& file1, const std::string& file2, double accuracy) {
    std::ifstream f1(file1);
    std::ifstream f2(file2);

    if (!f1.is_open() || !f2.is_open()) {
        std::cerr << "Error opening files." << std::endl;
        return false;
    }

    std::string line1, line2;
    // Skip headers
    std::getline(f1, line1);
    std::getline(f2, line2);

    bool fl1 = static_cast<bool>(std::getline(f1, line1));
    bool fl2 = static_cast<bool>(std::getline(f2, line2));
    while (fl1 && fl2) {
        std::istringstream iss1(line1);
        std::istringstream iss2(line2);

        double num1, num2;
        while (iss1 >> num1 && iss2 >> num2) {
					//std::cout << "Comparing " << num1 << " and " << num2 << std::endl;
            if (std::fabs(num1 - num2) > accuracy) {
								std::cout << "Difference found: " << num1 << " vs " << num2 << std::endl;
                return false;
            }
        }

        if (iss1 >> num1 || iss2 >> num2) {
						std::cerr << "Files have different number of columns." << std::endl;
            return false; // Different number of columns
        }

        fl1 = static_cast<bool>(std::getline(f1, line1));
        fl2 = static_cast<bool>(std::getline(f2, line2));
    }

    if (f1.good() != f2.good()) {
				std::cerr << "Files have different number of rows." << std::endl;
				return false; // Different number of rows
    }

    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <file1> <file2> <accuracy>" << std::endl;
        return 1;
    }

    std::string file1 = argv[1];
    std::string file2 = argv[2];
    double accuracy = 1.e-6;
		if (argc > 3)
			accuracy = std::stod(argv[3]);

    if (compareFiles(file1, file2, accuracy)) {
        //std::cout << "Files are similar within the given accuracy." << std::endl;
    } else {
        //std::cout << "Files differ beyond the given accuracy." << std::endl;
				return 1;
    }

    return 0;
}
