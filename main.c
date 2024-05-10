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

void crearMatrizDeEcuaciones(double sumaX, double sumaX2, double sumaY, double sumaXY, double matrizEcuaciones[2][2]) {
    matrizEcuaciones[0][0] = sumaX2;
    matrizEcuaciones[0][1] = sumaX;
    matrizEcuaciones[1][0] = sumaX;
    matrizEcuaciones[1][1] = 100;

    // Vector de términos independientes
    double vectorB[2] = {sumaXY, sumaY};

    printf("\nMatriz Ecuaciones:\n");
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            printf("%.2f ", matrizEcuaciones[i][j]);
        }
        printf("\n");
    }

    // Matriz aumentada
    double matrizAumentada[2][3];
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            matrizAumentada[i][j] = matrizEcuaciones[i][j];
        }
        matrizAumentada[i][2] = vectorB[i];
    }

    // Imprimir la matriz aumentada (opcional, solo para verificar)
    printf("\n Matriz aumentada:\n");
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            printf("%.2f ", matrizAumentada[i][j]);
        }
        printf("\n");
    }
}

//void calcularDeterminantes(double matrizEcuaciones[1][1], double *determinante, double *determinanteK, double *determinanteB) {
//    *determinante = matrizEcuaciones[0][0] * matrizEcuaciones[1][1] - matrizEcuaciones[0][1] * matrizEcuaciones[1][0];
//    *determinanteK = matrizEcuaciones[1][1];
//    *determinanteB = matrizEcuaciones[0][1];
//}


void solveToFunction(struct ChemicalElement periodicTable[], double a, double b){
    double xSum = 0;
    double xSquaredSum = 0;
    double ySum = 0;
    double xTimesYSum = 0;
    double matrizEcuaciones[2][2];
    //double determinante, determinanteK, determinanteB;
    calculateSum(periodicTable, &xSum, &xSquaredSum, &ySum, &xTimesYSum);
    printf("%f, %f, %f, %f", xSum, xSquaredSum, ySum, xTimesYSum);
    crearMatrizDeEcuaciones( xSum, xSquaredSum, ySum, xTimesYSum, matrizEcuaciones);
    //calcularDeterminantes(matrizEcuaciones, &determinante, &determinanteK, &determinanteB);
    //printf("%f, %f, %f", determinante, determinanteK, determinanteB);
}

int main(void) {
    struct ChemicalElement periodicTable[PERIODIC_TABLE_SIZE];
    double a;
    double b;
    loadData(periodicTable);
    solveToFunction(periodicTable, a, b);
    return 0;
}
