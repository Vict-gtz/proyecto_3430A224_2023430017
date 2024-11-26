#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <sstream>

using namespace std;


vector<vector<int>> leerCSV(const string &matriz_archivo) {
    vector<vector<int>> matriz;
    ifstream file(matriz_archivo);

    // Caso de error
    if (!file) {
        cerr << "Error: No se pudo abrir el archivo " << matriz_archivo << endl; 
        exit(1);
    }

    string line;
    bool esPrimeraFila = true; // Ignorar la primera fila

    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        vector<int> fila;

        if (esPrimeraFila) { // Ignorar encabezados de la primera fila
            esPrimeraFila = false;
            continue;
        }

        bool esPrimeraColumna = true;
        while (getline(ss, value, ',')) {
            if (esPrimeraColumna) { // Ignorar encabezados de la primera columna
                esPrimeraColumna = false;
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



int main(int argc, char **argv) { //proyecto secuenciaH.txt secuenciaV.txt matriz.csv penitencia_puntos_valor
    //Variables
    const int puntaje_penalidad = stoi(argv[4]);

    string linea;

    vector<char> secuencia_HORIZONTAL;
    vector<char> secuencia_VERTICAL; 
    string arreglo_ADN[4] = {"A", "C", "G", "T"};
    vector<vector<int>> matriz_puntuacion = leerCSV(argv[3]); //ATCG -> orden

    // Revisar caracteres validos
    auto is_valid_nucleotide = [](char c) {
        return c == 'A' || c == 'T' || c == 'C' || c == 'G';
    };

    // Abrir el primer archivo (HORIZONTAL)
    ifstream archivo_sec1(argv[1]);
    if (!archivo_sec1) {
        cerr << "Error: No se pudo abrir el archivo " << argv[1] << endl;
        return 1;
    }

    // Leer datos del primer archivo
    while (getline(archivo_sec1, linea)) {
        if (linea[0] == '>') continue; 
        for (char lineac : linea) {
            if (is_valid_nucleotide(lineac)) {
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
            if (is_valid_nucleotide(lineac)) {
                secuencia_VERTICAL.push_back(lineac);
            }
        }
    }
    archivo_sec2.close();

    // Mostrar los contenidos leídos de ambos archivos
    cout << "Contenido del archivo 1:" << endl;
    for (char c : secuencia_HORIZONTAL) {
        cout << c;
    }

    cout << "\nContenido del archivo 2:" << endl;
    for (char c : secuencia_VERTICAL) {
        cout << c;
    }

    // Imprimir la matriz
    cout << "Matriz leída desde el archivo CSV:" << endl;
    imprimirMatriz(matriz_puntuacion);

    return 0;
}