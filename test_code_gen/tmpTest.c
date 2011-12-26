#include <math.h>

int test_Atan2(const double* ind, double* dep) {
    double v_aux;
    int comparison = 0;
    double v_0_1, v_0_2, v_0_3, v_0_4, v_0_5, v_0_6, v_0_7, v_0_8, v_0_9, v_0_10, v_0_11, v_0_12, v_0_13, v_0_14, v_0_15, v_0_16, v_0_17, v_0_18, v_0_19;
    v_0_1 = ind[0];
    v_0_2 = sin(v_0_1);
    v_0_3 = cos(v_0_1);
    v_0_5 = sin(v_0_1);
    v_0_4 = cos(v_0_1);
    if (v_0_3 > ((double) 0)) v_0_6 = v_0_3;
    else if (v_0_3 < ((double) 0)) v_0_6 = -v_0_3;
    else v_0_6 = ((double) 0);
    if (v_0_5 > ((double) 0)) v_0_7 = v_0_5;
    else if (v_0_5 < ((double) 0)) v_0_7 = -v_0_5;
    else v_0_7 = ((double) 0);
    v_0_8 = (v_0_7) / v_0_6;
    v_0_10 = atan(v_0_8);
    v_0_9 = ((double) 1) + v_0_8 * v_0_8;
    v_0_11 = (v_0_6) / v_0_7;
    v_0_13 = atan(v_0_11);
    v_0_12 = ((double) 1) + v_0_11 * v_0_11;
    v_0_14 = 1.5708 - v_0_13;
    if (v_0_6 > v_0_7) v_0_15 = v_0_10;
    else v_0_15 = v_0_14;
    v_0_16 = 3.14159 - v_0_15;
    if (v_0_3 <= 0) v_0_17 = v_0_16;
    else v_0_17 = v_0_15;
    v_0_18 = 0 - v_0_17;
    if (v_0_5 <= 0) v_0_19 = v_0_18;
    else v_0_19 = v_0_17;
    dep[0] = v_0_19;
    return comparison;
}


#include <math.h>

int test_Atan2(const double* ind, double* dep) {
double v_aux;
int comparison = 0 ;
double v_0_1, v_0_2, v_0_3, v_0_4, v_0_5, v_0_6, v_0_7, v_0_8, v_0_9, v_0_10, v_0_11, v_0_12, v_0_13, v_0_14, v_0_15, v_0_16, v_0_17, v_0_18, v_0_19;
v_0_1 = ind[0];
v_0_2 = sin(v_0_1);
v_0_3 = cos(v_0_1);
v_0_5 = sin(v_0_1);
v_0_4 = cos(v_0_1);
if(v_0_3 < ((double)0)) v_0_6 = -v_0_3;
else v_0_6 = v_0_3;
if(v_0_5 < ((double)0)) v_0_7 = -v_0_5;
else v_0_7 = v_0_5;
v_0_8 = (v_0_7) / v_0_6;
v_0_10 = atan(v_0_8);
v_0_9 = ((double)1) + v_0_8 * v_0_8;
v_0_11 = (v_0_6) / v_0_7;
v_0_13 = atan(v_0_11);
v_0_12 = ((double)1) + v_0_11 * v_0_11;
v_0_14 = 1.5708 - v_0_13;
if(v_0_6 > v_0_7) v_0_15 = v_0_10;
 else v_0_15 = v_0_14;
v_0_16 = 3.14159 - v_0_15;
if(v_0_3 <= 0) v_0_17 = v_0_16;
 else v_0_17 = v_0_15;
v_0_18 = 0 - v_0_17;
if(v_0_5 <= 0) v_0_19 = v_0_18;
 else v_0_19 = v_0_17;
dep[0] = v_0_19;
return comparison;
}