#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <sstream>

using namespace std;

// Funcion para leer la matriz de puntuacion
vector<vector<int>> leerCSV(const string &matriz_archivo) {
    vector<vector<int>> matriz;
    ifstream file(matriz_archivo);

    if (!file) { // Caso de error
        cerr << "Error: No se pudo abrir el archivo " << matriz_archivo << endl; 
        exit(1);
    }

    string line;
    bool primeraFila = true; // Ignorar la primera fila

    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        vector<int> fila;

        if (primeraFila) { // Ignorar encabezados de la primera fila
            primeraFila = false;
            continue;
        }

        bool primeraColumna = true;
        while (getline(ss, value, ',')) {
            if (primeraColumna) { // Ignorar encabezados de la primera columna
                primeraColumna = false;
                continue;
            }

            try {
                fila.push_back(stoi(value)); // Convertir los valores a enteros
            } catch (const invalid_argument &e) {
                cerr << "Error: Valor no válido en el archivo CSV: \"" << value << "\"" << endl;
                exit(1);
            }
        }
        matriz.push_back(fila);
    }

    file.close();
    return matriz;
}

// Función para imprimir la matriz
void imprimirMatriz(const vector<vector<int>> &matriz) {
    for (const auto &fila : matriz) {
        for (int valor : fila) {
            cout << valor << " ";
        }
        cout << endl;
    }
}

// Función para inicializar la matriz con las penalizaciones
vector<vector<int>> matriz_inicial(vector<vector<int>> &matriz, int puntaje_penalidad, int filas, int columnas) {
    // Primera fila
    for (int j = 0; j < columnas; ++j) {
        matriz[0][j] = j * puntaje_penalidad;
    }

    // Primera columna
    for (int i = 0; i < filas; ++i) {
        matriz[i][0] = i * puntaje_penalidad;
    }
    return matriz;
}

// Funcion para rellenar la matriz de direcciones
vector<vector<int>> matriz_direcciones(vector<vector<int>> &matriz_direccion, int puntuacion, int fila, int columna, int caso_match, int caso_izquierda, int caso_arriba){
    // DIAGONAL = 0
    // IZQUIERDA = 1
    // ARRIBA = 2
    
    if (puntuacion == caso_match) {// Prioriza la diagonal
        matriz_direccion[fila][columna] = 0; 
        return matriz_direccion;
    } else {
        if (puntuacion == caso_arriba) {
            matriz_direccion[fila][columna] = 2; 
            return matriz_direccion;
        } else {
            matriz_direccion[fila][columna] = 1; 
            return matriz_direccion;
        }
    }
}

// Funcion para rellenar la matriz utilizando Needleman-Wunsch
vector<vector<int>> needleman_wunsch(vector<vector<int>> &matriz, const vector<vector<int>> &matriz_puntuacion, const vector<char> &secuencia_HORIZONTAL, const vector<char> &secuencia_VERTICAL, int puntaje_penalidad, const string arreglo_ADN[], int filas, int columnas) {
    // Matriz de direcciones vacia
    vector<vector<int>> matriz_direccion(filas, vector<int>(columnas, 0));
    
    // Aqui mapea los nucleotidos para comparar
    auto indice_nucleotido = [&arreglo_ADN](char nucleotido) {
        for (int i = 0; i < 4; ++i) {
            if (arreglo_ADN[i][0] == nucleotido) {
                return i;
            }
        }
        return -1;
    };

    // Iterar para rellenar matriz
    for (int i = 1; i < filas; ++i) { 
        for (int j = 1; j < columnas; ++j) {

            int indice_HORIZONTAL = indice_nucleotido(secuencia_HORIZONTAL[j - 1]);
            int indice_VERTICAL = indice_nucleotido(secuencia_VERTICAL[i - 1]);

            // Revisar los 3 casos
            int caso_match = matriz[i - 1][j - 1] + matriz_puntuacion[indice_VERTICAL][indice_HORIZONTAL]; // CASO: diagonal

            int caso_izquierda = matriz[i - 1][j] + puntaje_penalidad; // CASO: izquierda

            int caso_arriba = matriz[i][j - 1] + puntaje_penalidad; // CASO: arriba

            // Seleccionar el máximo de los tres valores calculados
            int puntuacion = max(caso_match, max(caso_arriba, caso_izquierda));
            matriz[i][j] = puntuacion;

            // Rellena matriz de direcciones
            matriz_direcciones(matriz_direccion,
                puntuacion,
                i, j,
                caso_match, caso_izquierda, caso_arriba
            );
        }
    }
    return matriz; // Retornar la matriz modificada
}




int main(int argc, char **argv) { //proyecto secuenciaH.txt secuenciaV.txt matriz.csv penitencia_puntos_valor
    //Variables
    const int puntaje_penalidad = stoi(argv[4]);
    string linea;
    string arreglo_ADN[4] = {"A", "C", "G", "T"};

    vector<char> secuencia_HORIZONTAL;
    vector<char> secuencia_VERTICAL; 

    vector<vector<int>> matriz_puntuacion = leerCSV(argv[3]);
    vector<vector<int>> secuencias_comparadas;


    // Revisar caracteres validos
    auto nucleotido_valido = [&arreglo_ADN](char c) {
        for (const string &nucleotido : arreglo_ADN) {
            if (nucleotido[0] == c) {
                return true;
            }
        }
        return false;
    };

    // Abrir el primer archivo (HORIZONTAL)
    ifstream archivo_sec1(argv[1]);
    if (!archivo_sec1) {
        cerr << "Error: No se pudo abrir el archivo " << argv[1] << endl;
        return 1;
    }

    // Leer datos del primer archivo
    while (getline(archivo_sec1, linea)) {
        if (linea[0] == '>') continue; // skip header
        for (char lineac : linea) {
            if (nucleotido_valido(lineac)) {
                secuencia_HORIZONTAL.push_back(lineac);
            }
        }
    }
    archivo_sec1.close();

    // Abrir el segundo archivo (VERTICAL)
    ifstream archivo_sec2(argv[2]);
    if (!archivo_sec2) {
        cerr << "Error: No se pudo abrir el archivo " << argv[2] << endl;
        return 1;
    }

    // Leer datos del segundo archivo
    while (getline(archivo_sec2, linea)){
        if (linea[0] == '>') continue; 
        for (char lineac : linea) {
            if (nucleotido_valido(lineac)) {
                secuencia_VERTICAL.push_back(lineac);
            }
        }
    }
    archivo_sec2.close();

    // Se crea matriz segun longitud de secuencias
    int filas = secuencia_HORIZONTAL.size() +1;
    int columnas = secuencia_VERTICAL.size() +1;

    vector<vector<int>> matriz(filas, vector<int>(columnas, 0));

    matriz_inicial(
        matriz,
        puntaje_penalidad,
        filas, columnas
    );

    // Aquí se va al proceso de needleman_wunsch
    needleman_wunsch(
        matriz, matriz_puntuacion,
        secuencia_HORIZONTAL, secuencia_VERTICAL,
        puntaje_penalidad,
        arreglo_ADN,
        filas, columnas
    );


    /*// Mostrar los contenidos leidos de ambos archivos
    cout << "Contenido del archivo 1:" << endl;
    for (char c : secuencia_HORIZONTAL) {
        cout << c;
    }/*

    cout << "\nContenido del archivo 2:" << endl;
    for (char c : secuencia_VERTICAL) {
        cout << c;
    }*/

    // Imprimir la matriz de Needleman
    cout << "Matriz de Needleman wunsch:" << endl;
    imprimirMatriz(matriz);

    return 0;
}