#include <nvector/nvector_col.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <dColVector.h>

int main() {


int length = 7;
std::vector<double> ex = {2,3,4,5,6,8,4};
ColumnVector data (length);
for (octave_f77_int_type i = 0; i < length; i++)
      data(i) = i;

N_Vector v;
SUNContext ctx;

int retval = SUNContext_Create(NULL, &ctx);
v = N_VMake_Col(length, ex.data(), ctx);

std::cout<<"built succesfully !!";

return 0;
}
