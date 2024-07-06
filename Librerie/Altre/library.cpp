#include <iostream>
#include "library.h"
using namespace std;
using namespace arma;

void Progress_Bar(int progress, int total, int bar_width) {

    float percentage = static_cast<float>(progress) / total;
    int pos = static_cast<int>(bar_width * percentage);

    std::string bar;
    for (int i = 0; i < bar_width; ++i) {
        if (i <= pos) bar += "◼︎";
        else bar += " ";
    }

    cout::print("|{}| {:3d} %\r", bar, int(percentage * 100.0));
    std::fflush(stdout);
}
