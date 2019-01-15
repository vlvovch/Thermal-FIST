/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/Utility.h"

#include <iostream>
#include <string>
#include <sstream>

#include "ThermalFISTConfig.h"

using namespace std;

namespace thermalfist {

  bool Disclaimer::PrintDisclaimer()
  {
    for (int i = 0; i < 79; ++i)
      cout << "#";
    cout << endl;

    cout << "#";
    for (int i = 0; i < 77; ++i)
      cout << " ";
    cout << "#";
    cout << endl;

    ostringstream sstr;

    sstr << "# "
         << "This is Thermal-FIST version "
         << ThermalFIST_VERSION_MAJOR
         << "." 
         << ThermalFIST_VERSION_MINOR;

    if (ThermalFIST_VERSION_DEVEL != 0)
      sstr << "." 
         << ThermalFIST_VERSION_DEVEL;

    string tmpstr = sstr.str();
    while (tmpstr.size() < 78)
      tmpstr += " ";
    tmpstr += "#";

    cout << tmpstr << endl;

    cout << "#";
    for (int i = 0; i < 77; ++i)
      cout << " ";
    cout << "#";
    cout << endl;

    // e-mail obfuscation (just in case)
    string email = "";
    email += 'v';
    email += char(email[0] - 7);
    email += email[0];
    email += char('a' + 2);
    email += char(email[email.size() - 1] + 5);
    email += "enk";
    email += "@";
    email += "o";
    char ch1 = email[email.size() - 2], ch2 = email[email.size() - 1];
    email[email.size() - 1] = ch1;
    email[email.size() - 2] = ch2;
    email += "fias";
    email += ".";
    email += "uni-frankfurt.de";

    sstr = ostringstream();

    sstr << "# "
      << "Copyright (c) 2019 Volodymyr Vovchenko "
      << "<" << email << ">";

    tmpstr = sstr.str();
    while (tmpstr.size() < 78)
      tmpstr += " ";
    tmpstr += "#";

    cout << tmpstr << endl;

    cout << "#";
    for (int i = 0; i < 77; ++i)
      cout << " ";
    cout << "#";
    cout << endl;

    sstr = ostringstream();

    sstr << "# "
      << "Distributed under the GNU General Public License 3.0 (GPLv3 or later)";

    tmpstr = sstr.str();
    while (tmpstr.size() < 78)
      tmpstr += " ";
    tmpstr += "#";

    cout << tmpstr << endl;

    cout << "#";
    for (int i = 0; i < 77; ++i)
      cout << " ";
    cout << "#";
    cout << endl;

    sstr = ostringstream();

    sstr << "# "
      << "The latest version is available at "
      << "https://github.com/vlvovch/Thermal-FIST";

    tmpstr = sstr.str();
    while (tmpstr.size() < 78)
      tmpstr += " ";
    tmpstr += "#";

    cout << tmpstr << endl;

    cout << "#";
    for (int i = 0; i < 77; ++i)
      cout << " ";
    cout << "#";
    cout << endl;

    for (int i = 0; i < 79; ++i)
      cout << "#";
    cout << endl;

    cout << endl;

    return Disclaimer::DisclaimerPrinted = true;
  }

  bool Disclaimer::DisclaimerPrinted = PrintDisclaimer();

} // namespace thermalfist