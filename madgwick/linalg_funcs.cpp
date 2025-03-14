#include "linalg_funcs.h"

float get_norm(const float (&arr)[3]) {
  float sum = 0;
  for (int i = 0; i < 3; i++) {
    sum += arr[i] * arr[i];
  }
  return sqrt(sum);
}

void unit(float (&arr)[3]) {
  float norm = get_norm(arr);
  if (norm >= 1e-12) {
    for (int i = 0; i < 3; i++) {
      arr[i] /= norm;
    }
  }
}

void cross(const float (&arr1)[3], const float (&arr2)[3], float (&arr3)[3]) {
  arr3[0] = arr1[1] * arr2[2] - arr1[2] * arr2[1];
  arr3[1] = arr2[0] * arr1[2] - arr2[2] * arr1[0];
  arr3[2] = arr1[0] * arr2[1] - arr1[1] * arr2[0];
}