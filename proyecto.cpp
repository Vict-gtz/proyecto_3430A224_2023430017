#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <sstream>
#include <climits>
#include <utility>

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

// Funcion para crear un arreglo de camino de regreso segun matriz de direcciones
vector<int> arreglo_CaminoRegreso(const vector<vector<int>> &matriz_direccion) {
    // Se crea un arreglo para guardar el camino de regreso
    vector<int> camino;
    int fila = matriz_direccion.size() - 1;
    int columna = matriz_direccion[0].size() - 1;

    while (fila > 0 && columna > 0) {
        camino.insert(camino.begin(), matriz_direccion[fila][columna]);

        if (matriz_direccion[fila][columna] == 0) { // Se mueve en la diagonal
            --fila; 
            --columna; 
        } else if (matriz_direccion[fila][columna] == 1) { // Se mueve hacia arriba
            --fila; 
        } else { // Se mueve hacia la izquierda
            --columna; 
        }
    }

    return camino;
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

// Función modificada para realizar Needleman-Wunsch y rellenar la matriz de direcciones
pair<vector<vector<int>>, vector<vector<int>>> needleman_wunsch(const vector<vector<int>> &matriz_puntuacion, const vector<char> &secuencia_HORIZONTAL, const vector<char> &secuencia_VERTICAL, int puntaje_penalidad, const string arreglo_ADN[], int filas, int columnas) {
    // Matriz de puntuación
    vector<vector<int>> matriz(filas, vector<int>(columnas, 0));

    // Matriz de direcciones (0: diagonal, 1: arriba, 2: izquierda)
    vector<vector<int>> matriz_direcciones(filas, vector<int>(columnas, 2));

    // Función para mapear nucleótido a índice
    auto indice_nucleotido = [&arreglo_ADN](char nucleotido) {
        for (int i = 0; i < 4; ++i) {
            if (arreglo_ADN[i][0] == nucleotido) {
                return i;
            }
        }
        return -1; // Valor inválido
    };

    // Inicializar la primera fila y columna con las penalizaciones
    for (int i = 0; i < filas; ++i) {
        matriz[i][0] = i * puntaje_penalidad;
        matriz_direcciones[i][0] = 1; // Arriba
    }
    for (int j = 0; j < columnas; ++j) {
        matriz[0][j] = j * puntaje_penalidad;
        matriz_direcciones[0][j] = 2; // Izquierda
    }

    // Llenar las matrices de puntuaciones y direcciones
    for (int i = 1; i < filas; ++i) {
        for (int j = 1; j < columnas; ++j) {
            // Índices en la matriz de puntuación
            int indice_HORIZONTAL = indice_nucleotido(secuencia_HORIZONTAL[j - 1]);
            int indice_VERTICAL = indice_nucleotido(secuencia_VERTICAL[i - 1]);

            // Revisar los 3 casos
            int caso_match = matriz[i - 1][j - 1] + matriz_puntuacion[indice_VERTICAL][indice_HORIZONTAL];
            int caso_arriba = matriz[i - 1][j] + puntaje_penalidad;
            int caso_izquierda = matriz[i][j - 1] + puntaje_penalidad;

            // Determinar el valor máximo y su dirección
            if (caso_match >= caso_arriba && caso_match >= caso_izquierda) {
                matriz[i][j] = caso_match;
                matriz_direcciones[i][j] = 0; // Diagonal
            } else if (caso_arriba >= caso_izquierda) {
                matriz[i][j] = caso_arriba;
                matriz_direcciones[i][j] = 1; // Arriba
            } else {
                matriz[i][j] = caso_izquierda;
                matriz_direcciones[i][j] = 2; // Izquierda
            }
        }
    }

    return make_pair(matriz, matriz_direcciones);
}


// Funcion para cambiar las secuencias
pair<vector<char>, vector<char>>  cambio_secuencias(vector<int> arreglo_direcciones, vector<char> &secuencia_HORIZONTAL, vector<char> &secuencia_VERTICAL){
    int offset_horizontal = 0; // Desplazamiento acumulado para secuencia_HORIZONTAL
    int offset_vertical = 0;   // Desplazamiento acumulado para secuencia_VERTICAL
    
    for (int i = arreglo_direcciones.size(); i > 0; --i) { 
        if (arreglo_direcciones[i] == 1){ // CASO: arriba
            secuencia_HORIZONTAL.insert(secuencia_HORIZONTAL.begin() + i + offset_horizontal, '-');
            offset_horizontal++;

        } else if (arreglo_direcciones[i] == 2) { // CASO: izquierda
            secuencia_VERTICAL.insert(secuencia_VERTICAL.begin() + i + offset_vertical, '-');
            offset_vertical++;
        }
        
    }

    return make_pair(secuencia_HORIZONTAL, secuencia_VERTICAL);
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



    // Revisar caracteres validos (para leer secuencias)
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

    // Aqui la matriz se llena con los valores del gap
    matriz_inicial(
        matriz,
        puntaje_penalidad,
        filas, columnas
    );

    // Aqui se va al proceso de needleman_wunsch
    auto resultado = needleman_wunsch(
        matriz_puntuacion,
        secuencia_HORIZONTAL,
        secuencia_VERTICAL,
        puntaje_penalidad,
        arreglo_ADN,
        filas, columnas
    );

    vector<vector<int>> matriz_punt = resultado.first;           // Matriz de puntuación
    vector<vector<int>> matriz_direcciones = resultado.second; // Matriz de direcciones


    // Obtener el camino de regreso a partir de la matriz de direcciones en un arreglo
    vector<int> arreglo_direcciones = arreglo_CaminoRegreso(matriz_direcciones);


    // Imprimir matrices
    cout << "Matriz de puntuación:\nA  T  C  G" << endl;
    imprimirMatriz(matriz_puntuacion);

    cout << "\n\nMatriz de Needleman Wusnch:" << endl;
    imprimirMatriz(matriz_punt);

    cout << "\nMatriz de Direcciones (camino de vuelta):" << endl;
    imprimirMatriz(matriz_direcciones);

    
    // Imprimir el camino
    cout << "\nCamino de regreso:";
    for (int i = 0; i < arreglo_direcciones.size(); ++i) { 
        cout << arreglo_direcciones[i] << " ";
    }

    // Cambio de secuencias
    auto secuencias_cambiadas = cambio_secuencias(arreglo_direcciones, secuencia_HORIZONTAL, secuencia_VERTICAL);

    vector<char> nueva_HORIZONTAL = secuencias_cambiadas.first;
    vector<char> nueva_VERTICAL = secuencias_cambiadas.second;


    // Mostrar los contenidos leidos de ambos archivos
    cout << "\nSecuencia horizontal:";
    for (char c : secuencia_HORIZONTAL) {
        cout << c;
    }

    cout << "\nSecuencia vertical:  ";
    for (char c : secuencia_VERTICAL) {
        cout << c;
    }

    // Guardar las secuencias arregladas en un archivo
    ofstream archivo_salida("secuencias_arregladas.txt");
    if (!archivo_salida) {
        cerr << "\nError: No se pudo crear el archivo de salida 'secuencias_arregladas.txt'" << endl;
        return 1;
    }

    archivo_salida << "> Secuencia 1:\n";
    for (char c : nueva_HORIZONTAL) {
        archivo_salida << c;
    }
    archivo_salida << "\n> Secuencia 2:\n";
    for (char c : nueva_VERTICAL) {
        archivo_salida << c;
    }

    archivo_salida.close();
    cout << "\n\nLas secuencias ajustadas se han guardado en 'secuencias_arregladas.txt'.\n";

    return 0;
}