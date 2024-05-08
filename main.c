#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cJSON.h"

const int PERIODIC_TABLE_SIZE = 100;

struct ChemicalElement {
    char name[50];
    int number;
    int actualNeutrons;
    double functionNeutrons;
};

void loadFile(char **data){
    FILE *file = fopen("..\\PeriodicTableJSON.json", "r");
    if (!file) {
        fprintf(stderr, "No se pudo abrir el archivo.\n");
        return;
    }

    // Read the content into data string
    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);
    *data = (char*)malloc(length + 1);
    fread(*data, 1, length, file);
    (*data)[length] = '\0';
    fclose(file);
}

void loadData(struct ChemicalElement periodicTable[]){
    char *data;
    loadFile(&data);
    cJSON *json =  cJSON_Parse(data);
    if (!json) {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr != NULL) {
            fprintf(stderr, "Error antes de: %s\n", error_ptr);
        }
        return;
    }

    cJSON *elements = cJSON_GetObjectItem(json, "elements");

    // Iterating the elements
    cJSON *element;
    int index = 0;
    cJSON_ArrayForEach(element, elements) {
        if(index >= PERIODIC_TABLE_SIZE)
            break;
        cJSON *name = cJSON_GetObjectItem(element, "name");
        cJSON *number = cJSON_GetObjectItem(element, "number");
        cJSON *atomic_mass = cJSON_GetObjectItem(element, "atomic_mass");

        // Assigning values to periodicTable
        snprintf(periodicTable[index].name, 50, "%s", name->valuestring);
        periodicTable[index].number = number->valueint;
        periodicTable[index].actualNeutrons = (int)round(atomic_mass->valuedouble - number->valueint);

        index ++;
    }

    cJSON_Delete(json);
    free(data);
}

void calculateSum(struct ChemicalElement periodicTable[], double *xSum, double *xSquaredSum, double *ySum, double *xTimesYSum){
    for(int i=0; i < PERIODIC_TABLE_SIZE; i++){
        double logX = log(periodicTable[i].number);
        double logY = (periodicTable[i].actualNeutrons != 0) ? log(periodicTable[i].actualNeutrons) : 0;
        *xSum += logX;
        *xSquaredSum += pow(logX, 2);
        *ySum += logY;
        *xTimesYSum += logX * logY;
    }
}

void solveToFunction(struct ChemicalElement periodicTable[]){
    double xSum = 0;
    double xSquaredSum = 0;
    double ySum = 0;
    double xTimesYSum = 0;
    calculateSum(periodicTable, &xSum, &xSquaredSum, &ySum, &xTimesYSum);
    printf("%f, %f, %f, %f", xSum, xSquaredSum, ySum, xTimesYSum);
}

int main(void) {
    struct ChemicalElement periodicTable[PERIODIC_TABLE_SIZE];
    loadData(periodicTable);
    solveToFunction(periodicTable);
    return 0;
}
