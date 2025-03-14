#ifndef LINALG_FUNCS_H
#define LINALG_FUNCS_H

#include <Arduino.h>
#include <math.h>

float get_norm(const float (&arr)[3]);
void unit(float (&arr)[3]);
void cross(const float (&arr1)[3], const float (&arr2)[3], float (&arr3)[3]);

#endif