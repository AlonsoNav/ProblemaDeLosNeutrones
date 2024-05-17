#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cJSON.h"

const int PERIODIC_TABLE_SIZE = 100;

struct ChemicalElement
{
    char name[50];
    int number;
    int actualNeutrons;
    double functionNeutrons;
};

void loadFile(char **data)
{
    FILE *file = fopen("..\\PeriodicTableJSON.json", "r");
    if (!file)
    {
        fprintf(stderr, "No se pudo abrir el archivo.\n");
        return;
    }

    // Read the content into data string
    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);
    *data = (char *)malloc(length + 1);
    fread(*data, 1, length, file);
    (*data)[length] = '\0';
    fclose(file);
}

void loadData(struct ChemicalElement periodicTable[])
{
    char *data;
    loadFile(&data);
    cJSON *json = cJSON_Parse(data);
    if (!json)
    {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr != NULL)
        {
            fprintf(stderr, "Error antes de: %s\n", error_ptr);
        }
        return;
    }

    cJSON *elements = cJSON_GetObjectItem(json, "elements");

    // Iterating the elements
    cJSON *element;
    int index = 0;
    cJSON_ArrayForEach(element, elements)
    {
        if (index >= PERIODIC_TABLE_SIZE)
            break;
        cJSON *name = cJSON_GetObjectItem(element, "name");
        cJSON *number = cJSON_GetObjectItem(element, "number");
        cJSON *atomic_mass = cJSON_GetObjectItem(element, "atomic_mass");

        // Assigning values to periodicTable
        snprintf(periodicTable[index].name, 50, "%s", name->valuestring);
        periodicTable[index].number = number->valueint;
        periodicTable[index].actualNeutrons = (int)round(atomic_mass->valuedouble - number->valueint);

        index++;
    }

    cJSON_Delete(json);
    free(data);
}

void calculateSum(struct ChemicalElement periodicTable[], double *xSum, double *xSquaredSum, double *ySum, double *xTimesYSum)
{
    for (int i = 0; i < PERIODIC_TABLE_SIZE; i++)
    {
        double logX = log(periodicTable[i].number);
        double logY = (periodicTable[i].actualNeutrons != 0) ? log(periodicTable[i].actualNeutrons) : 0;
        *xSum += logX;
        *xSquaredSum += pow(logX, 2);
        *ySum += logY;
        *xTimesYSum += logX * logY;
    }
}

void crearMatrizDeEcuaciones(double xSum, double xSquaredSum, double ySum, double xTimesYSum, double matrizEcuaciones[2][3])
{
    matrizEcuaciones[0][0] = PERIODIC_TABLE_SIZE;
    matrizEcuaciones[0][1] = xSum;
    matrizEcuaciones[0][2] = ySum;
    matrizEcuaciones[1][0] = xSum;
    matrizEcuaciones[1][1] = xSquaredSum;
    matrizEcuaciones[1][2] = xTimesYSum;
}

void calcularDeterminantes(double matrizEcuaciones[2][3], double *determinante, double *determinanteK, double *determinanteB)
{
    *determinante = matrizEcuaciones[0][0] * matrizEcuaciones[1][1] - matrizEcuaciones[0][1] * matrizEcuaciones[1][0];
    *determinanteK = matrizEcuaciones[0][2] * matrizEcuaciones[1][1] - matrizEcuaciones[0][1] * matrizEcuaciones[1][2];
    *determinanteB = matrizEcuaciones[0][0] * matrizEcuaciones[1][2] - matrizEcuaciones[0][2] * matrizEcuaciones[1][0];
}

void resolverAyB(double determinante, double determinanteK, double determinanteB, double *a, double *b)
{
    double k = determinanteK / determinante;
    *a = exp(k);
    *b = determinanteB / determinante;
}

void solveToFunction(struct ChemicalElement periodicTable[], double *a, double *b)
{
    double xSum, xSquaredSum, ySum, xTimesYSum, determinante, determinanteK, determinanteB;
    xSum = xSquaredSum = ySum = xTimesYSum = 0;
    double matrizEcuaciones[2][3];
    calculateSum(periodicTable, &xSum, &xSquaredSum, &ySum, &xTimesYSum);
    crearMatrizDeEcuaciones(xSum, xSquaredSum, ySum, xTimesYSum, matrizEcuaciones);
    calcularDeterminantes(matrizEcuaciones, &determinante, &determinanteK, &determinanteB);
    resolverAyB(determinante, determinanteK, determinanteB, a, b);
}

void calculateLinearFit(struct ChemicalElement periodicTable[], double *a, double *b)
{
    double sumZ = 0, sumZSquared = 0, sumN = 0, sumZTimesN = 0;
    for (int i = 0; i < PERIODIC_TABLE_SIZE; i++)
    {
        double Z = (double)periodicTable[i].number;
        double N = (double)periodicTable[i].actualNeutrons;
        sumZ += Z;
        sumZSquared += Z * Z;
        sumN += N;
        sumZTimesN += Z * N;
    }
    *a = (PERIODIC_TABLE_SIZE * sumZTimesN - sumZ * sumN) / (PERIODIC_TABLE_SIZE * sumZSquared - sumZ * sumZ);
    *b = (sumN - (*a) * sumZ) / PERIODIC_TABLE_SIZE;
}

void printPeriodicTable(struct ChemicalElement periodicTable[], double a, double b)
{
    printf("Tabla de Comparacion:\n");
    printf("-------------------------------------------------------------------------------------------------------------\n");
    printf("| NumE | N Real | N Calculado | redondeo al mas cercano | redondeo hacia cero | DifMasCercano | DifHaciaCero |\n");
    printf("-------------------------------------------------------------------------------------------------------------\n");

    for (int i = 0; i < PERIODIC_TABLE_SIZE; i++)
    {
        double Z = (double)periodicTable[i].number;
        double neutReal = (double)periodicTable[i].actualNeutrons;
        double neutPredicho = a * pow(Z,b);
        int predichoRounded = (int)round(neutPredicho);
        int predichoFloor = (int)floor(neutPredicho);
        int difRounded = abs(neutReal - predichoRounded);
        int difFloor = abs(neutReal - predichoFloor);

        printf("|%6d|%8d|%13.2f|%25d|%21d|%15d|%14d|\n", periodicTable[i].number, periodicTable[i].actualNeutrons, neutPredicho, predichoRounded, predichoFloor, difRounded, difFloor);
    }

    printf("--------------------------------------------------------------------------------------------------------------\n");
}

int main(void)
{
    struct ChemicalElement periodicTable[PERIODIC_TABLE_SIZE];
    double a, b;
    loadData(periodicTable);

    printf("Parametros a y b de la funcion de ajuste con determinantes:\n");
    solveToFunction(periodicTable, &a, &b);
    printPeriodicTable(periodicTable, a, b);

    return 0;
}
