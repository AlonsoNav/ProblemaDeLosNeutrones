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

void crearMatrizDeEcuaciones(double xSum, double xSquaredSum, double ySum, double xTimesYSum, double matrizEcuaciones[2][2]) {
    matrizEcuaciones[0][0] = xSum;
    matrizEcuaciones[0][1] = xSquaredSum;
    matrizEcuaciones[1][0] = ySum;
    matrizEcuaciones[1][1] = xTimesYSum;

    printf("\nMatriz Ecuaciones:\n");
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            printf("%.2f ", matrizEcuaciones[i][j]);
        }
        printf("\n");
    }
}

void calcularDeterminantes(double matrizEcuaciones[2][2], double *determinante, double *determinanteK, double *determinanteB) {
    double matrizAumentada[2][3] = {
            {matrizEcuaciones[0][1], matrizEcuaciones[0][0] , matrizEcuaciones[1][1]},
            {matrizEcuaciones[0][0], 100, matrizEcuaciones[1][0]}
    };

    // Imprimir la matriz aumentada (opcional, solo para verificar)
    printf("\nMatriz aumentada:\n");
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            printf("%.6f ", matrizAumentada[i][j]);
        }
        printf("\n");
    };

    // Calcular el determinante
    *determinante = matrizAumentada[0][0] * matrizAumentada[0][1] - matrizAumentada[1][0] * matrizAumentada[1][1];

    // Calcular los determinantes K y B
    *determinanteK = matrizAumentada[0][2] * matrizAumentada[0][1] - matrizAumentada[1][2] * matrizAumentada[1][1];
    *determinanteB = matrizAumentada[0][0] * matrizAumentada[0][2] - matrizAumentada[1][0] * matrizAumentada[1][2];
}

void resolverAyB(double *determinante, double *determinanteK, double *determinanteB, double *a, double *b) {
    double k = (*determinanteK)/(*determinante);
    *a = pow(2.71828, k);
    *b = (*determinanteB)/(*determinante);
}


void solveToFunction(struct ChemicalElement periodicTable[], double a, double b){
    double xSum = 0;
    double xSquaredSum = 0;
    double ySum = 0;
    double xTimesYSum = 0;
    double matrizEcuaciones[2][2];
    double determinante, determinanteK, determinanteB;
    calculateSum(periodicTable, &xSum, &xSquaredSum, &ySum, &xTimesYSum);
    printf("%f, %f, %f, %f", xSum, xSquaredSum, ySum, xTimesYSum);
    crearMatrizDeEcuaciones( xSum, xSquaredSum, ySum, xTimesYSum, matrizEcuaciones);
    calcularDeterminantes(matrizEcuaciones, &determinante, &determinanteK, &determinanteB);
    printf("%f, %f, %f", determinante, determinanteK, determinanteB);
    resolverAyB(&determinante, &determinanteK, &determinanteB, &a, &b);
    printf("\n\nA = %f, B = %f", a, b);
}

int main(void) {
    struct ChemicalElement periodicTable[PERIODIC_TABLE_SIZE];
    double a;
    double b;
    loadData(periodicTable);
    solveToFunction(periodicTable, a, b);
    return 0;
}
