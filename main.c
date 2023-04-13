#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>


#define elementSize 16

/**стандартная свертка разреженного и не разреженного векторов F2
 * M - длина неразреженного вектора
 * arrSize - количество ненулевых элементов резреженного вектора
 * inda - вектор номеров позиций ненулевых координт разреженного вектора
 * b - неразреженный вектор F2
 * c - результирующий вектор (F2)*/
void
computationF2(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c);

/**стандартная свертка разреженного и не разреженного векторов в Z
 * M - длина неразреженного вектора
 * arrSize - количество ненулевых элементов резреженного вектора
 * inda - вектор номеров позиций ненулевых координт разреженного вектора
 * b - неразреженный вектор F2
 * c - результирующий вектор (Z)*/
void
computationZ(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c);


/**перевод последовательности unsigned short в строку*/
char *toBinary(unsigned short n, int len) {
    char *binary = (char *) malloc(sizeof(char) * len);
    int k = 0;
    for (unsigned i = (1 << (len - 1)); i > 0; i = i / 2) {
        binary[k++] = (n & i) ? '1' : '0';
    }
    binary[k] = '\0';
    return binary;
}

/**генерирует разреженный массив, как массив индексов ненулевых элементов, упорядоченных в порядке возрастания.
 * arr - результирующий вектор индексов
 * m - максимальный индекс в векторе
 * res - количество элементов в результирующем векторе
 * */
void generateSparseArray(unsigned short *arr, int m, unsigned short resSize) {
    unsigned short i = 0;
    unsigned short a = 0;
    while (i < resSize && a < m) {
        unsigned int c = rand();
        a = a + (c % ((m - a) / (resSize - i)));
        arr[i++] = a++;
    }
}


/**преоброзует обычный вектор в "компактный"*/
void createDenseArray(unsigned short *res, int resSize, const unsigned short *donor, int donorSize) {
    int i;
    for (i = 0; i < donorSize / elementSize; i++) {
        unsigned short temp = 0;
        for (int j = 0; j < elementSize; j++) {
            temp = temp << 1;
            temp += donor[i * elementSize + j];
        }
        res[i] = temp;
    }
    if (donorSize % elementSize != 0) {

        unsigned short temp = 0;
        for (int j = i * elementSize; j < donorSize; j++) {
            temp = temp << 1;
            temp += donor[j];
        }
        temp = temp << (elementSize - (donorSize % elementSize));
        res[resSize - 1] = temp;
    }

}

/**генерирует обычный массив (b)*/
void generateDenseArray(unsigned short *arr, int m) {
    unsigned short temp = 0;
    for (int i = 0; i < m; i++) {
        temp = rand() % 20;
        if (temp > 10) {
            arr[i] = 0;
        } else {
            arr[i] = 1;
        }
    }
}

void calculateSparsedAndUsual2(unsigned short M, int arrSize, int m, const unsigned short *inda, unsigned short *b2,
                               unsigned short *res) {
    unsigned short masks[elementSize];
    unsigned short negMasks[elementSize];
    unsigned short modG = M % elementSize;
    unsigned short modGNeg = (elementSize - modG) % elementSize;
    unsigned short temp = -1;


    for (int i = 0; i < elementSize; i++) {
        masks[i] = temp;
        negMasks[i] = ~temp;
        temp -= 1 << (i);
        //printf("%s %s\n", toBinary(masks[i], elementSize), toBinary(negMasks[i], elementSize));
    }
    if (inda[0] < elementSize) {
        if (inda[arrSize - 1] >= elementSize) {
            int i = 0;
            for (i = 0; inda[i] < elementSize; i++) {
                unsigned short mod = inda[i] % elementSize;
                unsigned short modNeg = (elementSize - mod) % elementSize;

                //printf("mod %d modNeg %d\n", mod, modNeg);


                res[0] = res[0] ^ (b2[0] & masks[mod]) >> mod;
                printf("start\n");
                for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");
                unsigned short it = 1;
                for (int j = 1; j < m - 1; j++) {
                    printf("%s %s %s %s\n", toBinary(((b2[j] & masks[mod]) >> mod), elementSize),
                           toBinary((b2[j] & masks[mod]), elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]) << modNeg, elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]), elementSize));
                    res[it] = res[it] ^ (((b2[j] & masks[mod]) >> mod) + ((b2[j - 1] & negMasks[mod]) << modNeg));
                    it++;
                }

                for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");

                // todo подумать что будет с окончаниями элементов
                unsigned short temp2;
                temp2 = ((b2[m - 1] >> mod) + ((b2[m - 2] & negMasks[mod]) << modNeg));
                res[it] = res[it] ^ temp2;

                printf("%s %s\n", toBinary((b2[m - 1] >> modGNeg), elementSize),
                       toBinary(((b2[m - 2] & negMasks[modGNeg]) << modG), elementSize));

                unsigned short temp3 = (b2[m - 1] >> modGNeg) + ((b2[m - 2] & negMasks[modGNeg]) << modG);
                printf("%s\n", toBinary(temp3, elementSize));
                unsigned short temp4 = (temp3 & negMasks[mod]) << modNeg;
                res[0] = res[0] ^ temp4;


            }




            /**for по элементам, которые больше elementSize*/
            for (; i < arrSize; i++) {
                unsigned short mod = inda[i] % elementSize;
                unsigned short modNeg = (elementSize - mod) % elementSize;
                unsigned short start = inda[i] / elementSize;

                //printf("mod %d modNeg %d\n", mod, modNeg);


                res[start] = res[start] ^ (b2[0] & masks[mod]) >> mod;
                /*printf("start\n");*/
                /*for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }*/
                /*printf("\n");*/
                unsigned short it = start + 1;
                int j = 1;
                for (; it < m; it++, j++) {
                    /*printf("%s %s %s %s\n", toBinary(((b2[j] & masks[mod]) >> mod), elementSize),
                           toBinary((b2[j] & masks[mod]), elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]) << modNeg, elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]), elementSize));*/
                    res[it] = res[it] ^ (((b2[j] & masks[mod]) >> mod) + ((b2[j - 1] & negMasks[mod]) << modNeg));
                }

                //todo выдаст ошибку при modG == 0
                unsigned short newMod = (mod + modG) % elementSize;
                unsigned short newModNeg = (elementSize - newMod) % elementSize;
                /*printf("newMod %d newModNeg %d\n", newMod, newModNeg);

                printf("%d\n", j);*/
                it = 0;
                for (; j < m - 1; j++, it++) {
                    /*printf("new values %s %s %s %s\n", toBinary(((b2[j] & masks[newModNeg]) >> newMod), elementSize),
                           toBinary((b2[j] & masks[newModNeg]), elementSize),
                           toBinary((b2[j - 1] & negMasks[newModNeg]) << newMod, elementSize),
                           toBinary((b2[j - 1] & negMasks[newModNeg]), elementSize));*/
                    res[it] =
                            res[it] ^
                            (((b2[j] & masks[newModNeg]) >> newMod) + ((b2[j - 1] & negMasks[newModNeg]) << newMod));
                }

                /*for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");*/

                // todo подумать что будет с окончаниями элементов


                unsigned short temp3 = (b2[m - 1] >> modGNeg) + ((b2[m - 2] & negMasks[modGNeg]) << modG);

                res[it] = res[it] ^ temp3;
            }


        } else {
            for (int i = 0; i < arrSize; i++) {
                unsigned short mod = inda[i] % elementSize;
                unsigned short modNeg = (elementSize - mod) % elementSize;

                //printf("mod %d modNeg %d\n", mod, modNeg);


                res[0] = res[0] ^ (b2[0] & masks[mod]) >> mod;
                printf("start\n");
                for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");
                unsigned short it = 1;
                for (int j = 1; j < m - 1; j++) {
                    printf("%s %s %s %s\n", toBinary(((b2[j] & masks[mod]) >> mod), elementSize),
                           toBinary((b2[j] & masks[mod]), elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]) << modNeg, elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]), elementSize));
                    res[it] = res[it] ^ (((b2[j] & masks[mod]) >> mod) + ((b2[j - 1] & negMasks[mod]) << modNeg));
                    it++;
                }

                for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");

                // todo подумать что будет с окончаниями элементов
                unsigned short temp2;
                temp2 = ((b2[m - 1] >> mod) + ((b2[m - 2] & negMasks[mod]) << modNeg));
                res[it] = res[it] ^ temp2;

                printf("%s %s\n", toBinary((b2[m - 1] >> modGNeg), elementSize),
                       toBinary(((b2[m - 2] & negMasks[modGNeg]) << modG), elementSize));

                unsigned short temp3 = (b2[m - 1] >> modGNeg) + ((b2[m - 2] & negMasks[modGNeg]) << modG);
                printf("%s\n", toBinary(temp3, elementSize));
                unsigned short temp4 = (temp3 & negMasks[mod]) << modNeg;
                res[0] = res[0] ^ temp4;
            }

        }
    } else {
        for (int i = 0; i < arrSize; i++) {
            unsigned short mod = inda[i] % elementSize;
            unsigned short modNeg = (elementSize - mod) % elementSize;
            unsigned short start = inda[i] / elementSize;

            //printf("mod %d modNeg %d\n", mod, modNeg);


            res[start] = res[start] ^ (b2[0] & masks[mod]) >> mod;
            /*printf("start\n");*/
            /*for (int z = 0; z < m; z++) {
                printf("%s ", toBinary(res[z], elementSize));
            }*/
            /*printf("\n");*/
            unsigned short it = start + 1;
            int j = 1;
            for (; it < m; it++, j++) {
                /*printf("%s %s %s %s\n", toBinary(((b2[j] & masks[mod]) >> mod), elementSize),
                       toBinary((b2[j] & masks[mod]), elementSize),
                       toBinary((b2[j - 1] & negMasks[mod]) << modNeg, elementSize),
                       toBinary((b2[j - 1] & negMasks[mod]), elementSize));*/
                res[it] = res[it] ^ (((b2[j] & masks[mod]) >> mod) + ((b2[j - 1] & negMasks[mod]) << modNeg));
            }

            //todo выдаст ошибку при modG == 0
            unsigned short newMod = (mod + modG) % elementSize;
            unsigned short newModNeg = (elementSize - newMod) % elementSize;
            /*printf("newMod %d newModNeg %d\n", newMod, newModNeg);

            printf("%d\n", j);*/
            it = 0;
            for (; j < m - 1; j++, it++) {
                /*printf("new values %s %s %s %s\n", toBinary(((b2[j] & masks[newModNeg]) >> newMod), elementSize),
                       toBinary((b2[j] & masks[newModNeg]), elementSize),
                       toBinary((b2[j - 1] & negMasks[newModNeg]) << newMod, elementSize),
                       toBinary((b2[j - 1] & negMasks[newModNeg]), elementSize));*/
                res[it] =
                        res[it] ^
                        (((b2[j] & masks[newModNeg]) >> newMod) + ((b2[j - 1] & negMasks[newModNeg]) << newMod));
            }

            /*for (int z = 0; z < m; z++) {
                printf("%s ", toBinary(res[z], elementSize));
            }
            printf("\n");*/

            // todo подумать что будет с окончаниями элементов


            unsigned short temp3 = (b2[m - 1] >> modGNeg) + ((b2[m - 2] & negMasks[modGNeg]) << modG);

            res[it] = res[it] ^ temp3;
        }
    }
    res[m - 1] = res[m - 1] & masks[modG];

    for (int z = 0; z < m; z++) {
        printf("%s ", toBinary(res[z], elementSize));
    }
    printf("\n");
}

/**разжимает вектор*/
void decompressArray(unsigned short M, unsigned short m, unsigned short *arrTo, unsigned short *arrFrom) {
    int it = M - 1;
    unsigned short temp = arrFrom[m - 1];
    unsigned short mod = M % elementSize;
    temp = temp >> (elementSize - (M % elementSize)) % elementSize;
    //обработка последнего элемента массива
    for (int i = 0; i < mod; i++) {
        arrTo[it--] = temp & 1;
        temp = temp >> 1;
    }

    //обработка основной части элементов
    for (int i = m - 2; i >= 0; i--) {
        temp = arrFrom[i];
        for (int j = 0; j < elementSize; j++) {
            arrTo[it--] = temp & 1;
            temp = temp >> 1;
        }
    }

}

int main2() {

    srand(time(0));

    /**|e1| + |e2| <= t
     * |e1| = |e2| - не обязательно
     *
     * Для примера
     * n = 4801
     * |h1| = |h2| = 45
     * |e1| + |e2| = 84
     * T = 45/2+4*/
    unsigned short hLength = 45;
    unsigned short eLength = 42;
    unsigned T = hLength / 2 + 5;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 4801;

    /**компактное хранение h1*/
    unsigned short h1Compact[hLength];

    /**компактное хранение h2*/
    unsigned short h2Compact[hLength];

    /**компактное хранение h1 в обратном порядке*/
    unsigned short h1TransCompact[hLength];

    /**компактное хранение h2 в обратном порядке*/
    unsigned short h2TransCompact[hLength];

    /**компактное хранение e1*/
    unsigned short e1Compact[eLength];

    /**компактное хранение e2*/
    unsigned short e2Compact[eLength];

    /**заполняем вектора данными*/
    generateSparseArray(h1Compact, n, hLength);
    generateSparseArray(h2Compact, n, hLength);
    generateSparseArray(e1Compact, n, eLength);
    generateSparseArray(e2Compact, n, eLength);

    /*for (int i = 0; i < hLength; i++) {
        h1TransCompact[hLength - 1 - i] = n - 1 - h1Compact[i];
        h2TransCompact[hLength - 1 - i] = n - 1 - h2Compact[i];
    }*/
    for (int i = 0; i < hLength; i++) {
        h1TransCompact[i] = n - h1Compact[i];
        h2TransCompact[i] = n - h2Compact[i];
    }

    unsigned short e1[n];
    unsigned short e2[n];
    for (int i = 0; i < n; i++)
        e1[i] = e2[i] = 0;
    for (int i = 0; i < eLength; i++)
        e1[e1Compact[i]] = 1;
    for (int i = 0; i < eLength; i++)
        e2[e2Compact[i]] = 1;

    unsigned short c1[n];
    unsigned short c2[n];
    unsigned short s[n];

    computationF2(n, hLength, h1Compact, e1, c1);
    computationF2(n, hLength, h2Compact, e2, c2);
    for (int i = 0; i < n; i++) {
        s[i] = c1[i] ^ c2[i];
    }

    /**посчитали S*/
    for (int i = 0; i < 150; i++)
        printf("%d", s[i]);


    printf("\n");

    unsigned short u[n], v[n];
    for (int i = 0; i < n; i++)
        u[i] = v[i] = 0;

    unsigned short sTemp[n];
    for (int i = 0; i < n; i++) {
        sTemp[i] = s[i];
    }

    for (int z = 0; z < num_it; z++) {
        unsigned short upc1[n], upc2[n];
        for (int i = 0; i < n; i++)
            upc1[i] = upc2[i] = 0;
        computationZ(n, hLength, h1TransCompact, sTemp, upc1);
        computationZ(n, hLength, h2TransCompact, sTemp, upc2);
        for (int j = 0; j < n; j++) {
            if (upc1[j] >= T)
                u[j] = u[j] ^ 1;
            if (upc2[j] >= T)
                v[j] = v[j] ^ 1;
        }

        for (int i = 0; i < n; i++) {
            sTemp[i] = c1[i] = c2[i] = 0;
        }
        computationF2(n, hLength, h1Compact, u, c1);
        computationF2(n, hLength, h2Compact, v, c2);
        for (int i = 0; i < n; i++) {
            sTemp[i] = s[i] ^ c1[i] ^ c2[i];
        }

        /**проверка s` на ноль*/
        flag = 0;
        for (int i = 0; i < n; i++) {
            if (sTemp[i] != 0) {
                flag = 1;
                break;
            }
        }
        if (!flag) {
            exitFlag = 1;
            break;
        }

    }

    if (exitFlag) {
        printf("e1\n");
        for (int i = 0; i < n; i++) {
            printf("%d", e1[i]);
        }
        printf("\n");
        printf("e2\n");
        for (int i = 0; i < n; i++) {
            printf("%d", e2[i]);
        }
        printf("\n");
        printf("u\n");
        for (int i = 0; i < n; i++) {
            printf("%d", u[i]);
        }
        printf("\n");
        printf("v\n");
        for (int i = 0; i < n; i++) {
            printf("%d", v[i]);
        }
        printf("\n");
    } else {
        printf("fail\n");
    }
    return 0;
}

int main() {
    /*unsigned short inda[1] = {0};
    unsigned short b[7] = {1, 0, 1, 0, 0, 0, 0};
    unsigned short c[7] = {0,0,0,0,0,0,0};

    */
    /*for(int j = 0;j<7;j++) {
        for(int i = 0;i<7;i++)
            c[i] = 0;
        inda[0] = j;
        computationF2(7, 1, inda, b, c);
        for (int i = 0; i < 7; ++i) {
            printf("%d ", c[i]);
        }
        printf("\n");
    }*//*

    calculateSparsedAndUsual2(7,1,1,inda,b,c);*/

    srand(time(0));


    unsigned short M = 53;
    unsigned short m = (M + elementSize - 1) / elementSize;
    unsigned short inda[6] = {1,2,3,4,5,6};
    unsigned short b[M];
    unsigned short c2[M];
    for (int i = 0; i < M; i++) {
        c2[i] = 0;
    }
    unsigned short b2[m], c[m];
    for (int i = 0; i < m; i++)
        b2[i] = c[i] = 0;

    generateDenseArray(b, M);
    createDenseArray(b2, m, b, M);
    for (int i = 0; i < M; i++) {
        printf("%d", b[i]);
    }
    printf("\n");
    for (int i = 0; i < m; i++) {
        printf("%s ", toBinary(b2[i], elementSize));
    }
    printf("\n");
    printf("\n");

    calculateSparsedAndUsual2(M, 6, m, inda, b2, c);
    computationF2(M, 6, inda, b, c2);
    for (int i = 0; i < M; i++) {
        if (i != 0 && (i % elementSize) == 0)
            printf(" ");
        printf("%d", c2[i]);
    }
    printf("\n");

    /*printf("%d",8 << -1);*/

    //main2();
    return 0;
}


void
computationF2(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < M; j++) {
            c[j] = (c[j] + b[(j - inda[i] + M) % M]) % 2;
        }
    }
}

void
computationZ(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < M; j++) {
            c[j] = (c[j] + b[(j - inda[i] + M) % M]);
        }
    }
}
