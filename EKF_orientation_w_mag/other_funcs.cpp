#include "other_funcs.h"

void print_data(float (&arr)[3]) {
    for (int i = 0; i < 3; i++) {
        Serial.print(arr[i]);
        Serial.print(",");
    }
    Serial.println();
}
