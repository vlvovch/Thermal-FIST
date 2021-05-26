/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/Utility.h"

#include <iostream>
#include <sstream>

#include "ThermalFISTConfig.h"


// For time keeping
// Windows
#ifdef _WIN32
#include <Windows.h>
#else
#include <time.h>
#include <sys/time.h>
#endif

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
    if (Disclaimer::DisclaimerPrinted)
      return true;
    
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

    tmpstr = "Copyright (c) 2021 Volodymyr Vovchenko <" + email + ">";

    tmpstr = OutputString(tmpstr);

    cout << tmpstr << endl;

    cout << "#" << string(77, ' ') << "#" << endl;

    tmpstr = "Distributed under the GNU General Public License 3.0 (GPLv3 or later)";

    tmpstr = OutputString(tmpstr);

    cout << tmpstr << endl;

    cout << "#" << string(77, ' ') << "#" << endl;

    tmpstr = "Please cite when using this code:";
    tmpstr = OutputString(tmpstr);
    cout << tmpstr << endl;

    tmpstr = "V. Vovchenko, H. Stoecker, Comput. Phys. Commun. 244, 295 (2019)";
    tmpstr = OutputString(tmpstr);
    cout << tmpstr << endl;

    //tmpstr = "V. Vovchenko, H. Stoecker, arXiv:1901.05249 [nucl-th]";
    //tmpstr = OutputString(tmpstr);
    //cout << tmpstr << endl;

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

  long long stringToLongLong(const string &str) {
    long long ret = 0;
    int ist = 0, mn = 1;
    if (str.size() > 0 && str[0] == '-') {
      mn = -1;
      ist = 1;
    }
    for (size_t i = ist; i < str.size(); ++i) {
      if (str[i] >= '0' && str[i] <= '9') {
        ret *= 10;
        ret += static_cast<long long>(str[i] - '0');
      }
    }
    return ret * mn;
  }

  // Time keeping
  // Windows
  #ifdef _WIN32
  double get_wall_time(){
      LARGE_INTEGER time,freq;
      if (!QueryPerformanceFrequency(&freq)){
          //  Handle error
          return 0;
      }
      if (!QueryPerformanceCounter(&time)){
          //  Handle error
          return 0;
      }
      return (double)time.QuadPart / freq.QuadPart;
  }

  double get_cpu_time(){
      FILETIME a,b,c,d;
      if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
          //  Returns total user time.
          //  Can be tweaked to include kernel times as well.
          return
              (double)(d.dwLowDateTime |
              ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
      }else{
          //  Handle error
          return 0;
      }
  }

  //  Posix/Linux
  #else
  double get_wall_time(){
      struct timeval time;
      if (gettimeofday(&time,NULL)){
          //  Handle error
          return 0;
      }
      return (double)time.tv_sec + (double)time.tv_usec * .000001;
  }
  double get_cpu_time(){
      return (double)clock() / CLOCKS_PER_SEC;
  }
  #endif

} // namespace thermalfist