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

  namespace {
    string NumberToString(int num) {
      stringstream strs;
      strs << num;
      return strs.str();
    }

    string OutputString(const string &strin) {
      string ret = "";
      ret += "# ";
      ret += strin;
      while (ret.size() < 78)
        ret += " ";
      ret += "#";
      return ret;
    }
  }

  bool Disclaimer::PrintDisclaimer()
  {
    cout << string(79, '#') << endl;

    cout << "#" << string(77, ' ') << "#" << endl;

    string tmpstr = "";

    tmpstr += "This is Thermal-FIST version ";
    tmpstr += NumberToString(ThermalFIST_VERSION_MAJOR);
    tmpstr += ".";
    tmpstr += NumberToString(ThermalFIST_VERSION_MINOR);

    if (ThermalFIST_VERSION_DEVEL != 0) {
      tmpstr += ".";
      tmpstr += NumberToString(ThermalFIST_VERSION_DEVEL);
    }

    tmpstr = OutputString(tmpstr);

    cout << tmpstr << endl;

    cout << "#" << string(77, ' ') << "#" << endl;

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

    tmpstr = "Copyright (c) 2019 Volodymyr Vovchenko <" + email + ">";

    tmpstr = OutputString(tmpstr);

    cout << tmpstr << endl;

    cout << "#" << string(77, ' ') << "#" << endl;

    tmpstr = "Distributed under the GNU General Public License 3.0 (GPLv3 or later)";

    tmpstr = OutputString(tmpstr);

    cout << tmpstr << endl;

    cout << "#" << string(77, ' ') << "#" << endl;

    tmpstr = "The latest version is available at https://github.com/vlvovch/Thermal-FIST";

    tmpstr = OutputString(tmpstr);

    cout << tmpstr << endl;

    cout << "#" << string(77, ' ') << "#" << endl;

    cout << string(79, '#') << endl;

    cout << endl;

    return Disclaimer::DisclaimerPrinted = true;
  }

  bool Disclaimer::DisclaimerPrinted = PrintDisclaimer();

} // namespace thermalfist