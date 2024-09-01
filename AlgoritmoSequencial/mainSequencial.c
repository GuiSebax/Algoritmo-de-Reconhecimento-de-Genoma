#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define A 0
#define T 1
#define G 2
#define C 3
#define X 4

#define maxSeq 1000
#define sair 11

char mapaBases[5] = {'A', 'T', 'G', 'C', '-'};

int seqMaior[maxSeq], seqMenor[maxSeq];
int alinhaGMaior[maxSeq], alinhaGMenor[maxSeq];
int matrizEscores[maxSeq + 1][maxSeq + 1];

int tamSeqMaior, tamSeqMenor, tamAlinha, penalGap;
int grauMuta, escoreDiag, escoreLin, escoreCol;
int matrizPesos[4][4] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

int indRef, nTrocas, linPMaior, colPMaior, PMaior;
int linUMaior, colUMaior, UMaior, blocoTamanho;

int leTamMaior(void) {
    printf("\nLeitura do Tamanho da Sequencia Maior:");
    do {
        printf("\nDigite 0 < valor < %d = ", maxSeq);
        scanf("%d", &tamSeqMaior);
    } while ((tamSeqMaior < 1) || (tamSeqMaior > maxSeq));
}

int leTamMenor(void) {
    printf("\nLeitura do Tamanho da Sequencia Menor:");
    do {
        printf("\nDigite 0 < valor <= %d = ", tamSeqMaior);
        scanf("%d", &tamSeqMenor);
    } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));
}

int lePenalidade(void) {
    int penal;
    printf("\nLeitura da Penalidade de Gap:");
    do {
        printf("\nDigite valor >= 0 = ");
        scanf("%d", &penal);
    } while (penal < 0);
    return penal;
}

void leMatrizPesos() {
    int i, j;
    printf("\nLeitura da Matriz de Pesos:\n");
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            printf("Digite valor %c x %c = ", mapaBases[i], mapaBases[j]);
            scanf("%d", &(matrizPesos[i][j]));
        }
        printf("\n");
    }
}

void mostraMatrizPesos(void) {
    int i, j;
    printf("\nMatriz de Pesos Atual:");
    printf("\n%4c%4c%4c%4c%4c\n", ' ', 'A', 'T', 'G', 'C');
    for (i = 0; i < 4; i++) {
        printf("%4c", mapaBases[i]);
        for (j = 0; j < 4; j++)
            printf("%4d", matrizPesos[i][j]);
        printf("\n");
    }
}

int leGrauMutacao(void) {
    int prob;
    printf("\nLeitura da Porcentagem Maxima de Mutacao Aleatoria:\n");
    do {
        printf("\nDigite 0 <= valor <= 100 = ");
        scanf("%d", &prob);
    } while ((prob < 0) || (prob > 100));
    return prob;
}

void leSequenciasArquivo(char *file1, char *file2) {
    FILE *fp1 = fopen(file1, "r");
    FILE *fp2 = fopen(file2, "r");
    if (fp1 == NULL || fp2 == NULL) {
        printf("Erro ao abrir os arquivos de sequências.\n");
        exit(1);
    }

    char base;
    tamSeqMaior = 0;
    tamSeqMenor = 0;

    while (fscanf(fp1, "%c", &base) != EOF && tamSeqMaior < maxSeq) {
        switch (base) {
            case 'A': seqMaior[tamSeqMaior++] = A; break;
            case 'T': seqMaior[tamSeqMaior++] = T; break;
            case 'G': seqMaior[tamSeqMaior++] = G; break;
            case 'C': seqMaior[tamSeqMaior++] = C; break;
        }
    }

    while (fscanf(fp2, "%c", &base) != EOF && tamSeqMenor < maxSeq) {
        switch (base) {
            case 'A': seqMenor[tamSeqMenor++] = A; break;
            case 'T': seqMenor[tamSeqMenor++] = T; break;
            case 'G': seqMenor[tamSeqMenor++] = G; break;
            case 'C': seqMenor[tamSeqMenor++] = C; break;
        }
    }

    fclose(fp1);
    fclose(fp2);
}

void leSequenciasTeclado(void) {
    int i, erro;
    char seqMaiorAux[maxSeq], seqMenorAux[maxSeq];

    printf("\nLeitura das Sequencias:\n");

    // Sequência Maior
    do {
        printf("\nPara a Sequencia Maior, Digite apenas caracteres 'A', 'T', 'G' e 'C'\n> ");
        fgets(seqMaiorAux, maxSeq, stdin);
        tamSeqMaior = strlen(seqMaiorAux) - 1;
    } while (tamSeqMaior < 1);

    for (i = 0; i < tamSeqMaior; i++) {
        switch (seqMaiorAux[i]) {
            case 'A': seqMaior[i] = A; break;
            case 'T': seqMaior[i] = T; break;
            case 'G': seqMaior[i] = G; break;
            case 'C': seqMaior[i] = C; break;
            default: erro = 1;
        }
    }

    // Sequência Menor
    do {
        printf("\nPara a Sequencia Menor, Digite apenas caracteres 'A', 'T', 'G' e 'C'\n> ");
        fgets(seqMenorAux, maxSeq, stdin);
        tamSeqMenor = strlen(seqMenorAux) - 1;
    } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));

    for (i = 0; i < tamSeqMenor; i++) {
        switch (seqMenorAux[i]) {
            case 'A': seqMenor[i] = A; break;
            case 'T': seqMenor[i] = T; break;
            case 'G': seqMenor[i] = G; break;
            case 'C': seqMenor[i] = C; break;
            default: erro = 1;
        }
    }
}

void geraSequenciasAleatorias(void) {
    int i, dif, probAux;
    char base;

    printf("\nGeracao Aleatoria das Sequencias:\n");

    // Sequência Maior
    for (i = 0; i < tamSeqMaior; i++) {
        base = rand() % 4;
        seqMaior[i] = base;
    }

    dif = tamSeqMaior - tamSeqMenor;
    indRef = 0;
    if (dif > 0)
        indRef = rand() % dif;

    for (i = 0; i < tamSeqMenor; i++)
        seqMenor[i] = seqMaior[indRef + i];

    i = 0;
    nTrocas = 0;
    while ((i < tamSeqMenor) && (nTrocas < ((grauMuta * tamSeqMenor) / 100))) {
        probAux = rand() % 100 + 1;
        if (probAux <= grauMuta) {
            seqMenor[i] = (seqMenor[i] + (rand() % 3) + 1) % 4;
            nTrocas++;
        }
        i++;
    }

    printf("\nSequencias Geradas: Dif = %d, IndRef = %d, NTrocas = %d\n", dif, indRef, nTrocas);
}

void geraMatrizEscoresMPI(int rank, int size, int blocoTamanho) {
    int lin, col, peso, inicio, fim;
    int numBlocos, blocoAtual, elementosNoBloco;

    inicio = (tamSeqMenor * rank) / size + 1;
    fim = (tamSeqMenor * (rank + 1)) / size;

    if (rank == 0) {
        // Inicializa linha e coluna de penalidades
        for (col = 0; col <= tamSeqMaior; col++)
            matrizEscores[0][col] = -1 * (col * penalGap);

        for (lin = 0; lin <= tamSeqMenor; lin++)
            matrizEscores[lin][0] = -1 * (lin * penalGap);
    }

    for (lin = inicio; lin <= fim; lin++) {
        for (col = 1; col <= tamSeqMaior; col++) {
            peso = matrizPesos[seqMenor[lin - 1]][seqMaior[col - 1]];
            escoreDiag = matrizEscores[lin - 1][col - 1] + peso;
            escoreLin = matrizEscores[lin][col - 1] - penalGap;
            escoreCol = matrizEscores[lin - 1][col] - penalGap;

            if (escoreDiag >= escoreLin && escoreDiag >= escoreCol)
                matrizEscores[lin][col] = escoreDiag;
            else if (escoreLin > escoreCol)
                matrizEscores[lin][col] = escoreLin;
            else
                matrizEscores[lin][col] = escoreCol;
        }

        // Envio da linha em blocos para o processo 0
        numBlocos = (tamSeqMaior + 1) / blocoTamanho;
        if ((tamSeqMaior + 1) % blocoTamanho != 0) numBlocos++;
        
        for (blocoAtual = 0; blocoAtual < numBlocos; blocoAtual++) {
            int inicioBloco = blocoAtual * blocoTamanho;
            elementosNoBloco = blocoTamanho;
            if (inicioBloco + elementosNoBloco > tamSeqMaior + 1) {
                elementosNoBloco = tamSeqMaior + 1 - inicioBloco;
            }

            if (rank != 0) {
                MPI_Send(&matrizEscores[lin][inicioBloco], elementosNoBloco, MPI_INT, 0, lin, MPI_COMM_WORLD);
            } else {
                if (rank != 0) {
                    MPI_Recv(&matrizEscores[lin][inicioBloco], elementosNoBloco, MPI_INT, rank, lin, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }

    if (rank == 0) {
        linPMaior = 1;
        colPMaior = 1;
        PMaior = matrizEscores[1][1];
        linUMaior = 1;
        colUMaior = 1;
        UMaior = matrizEscores[1][1];

        for (lin = 1; lin <= tamSeqMenor; lin++) {
            for (col = 1; col <= tamSeqMaior; col++) {
                if (PMaior < matrizEscores[lin][col]) {
                    linPMaior = lin;
                    colPMaior = col;
                    PMaior = matrizEscores[lin][col];
                }
                if (UMaior <= matrizEscores[lin][col]) {
                    linUMaior = lin;
                    colUMaior = col;
                    UMaior = matrizEscores[lin][col];
                }
            }
        }
        printf("\nMatriz de escores Gerada.");
        printf("\nPrimeiro Maior escore = %d na celula [%d,%d]", PMaior, linPMaior, colPMaior);
        printf("\nUltimo Maior escore = %d na celula [%d,%d]", UMaior, linUMaior, colUMaior);
    }
}


void mostraMatrizEscores(void) {
    int i, lin, col;

    printf("\nMatriz de escores Atual:\n");

    printf("%4c%4c", ' ', ' ');
    for (i = 0; i <= tamSeqMaior; i++)
        printf("%4d", i);
    printf("\n");

    printf("%4c%4c%4c", ' ', ' ', '-');
    for (i = 0; i < tamSeqMaior; i++)
        printf("%4c", mapaBases[seqMaior[i]]);
    printf("\n");

    printf("%4c%4c", '0', '-');
    for (col = 0; col <= tamSeqMaior; col++)
        printf("%4d", matrizEscores[0][col]);
    printf("\n");

    for (lin = 1; lin <= tamSeqMenor; lin++) {
        printf("%4d%4c", lin, mapaBases[seqMenor[lin - 1]]);
        for (col = 0; col <= tamSeqMaior; col++) {
            printf("%4d", matrizEscores[lin][col]);
        }
        printf("\n");
    }
}

void gravaMatrizEscoresEmArquivo(void) {
    FILE *fp = fopen("matriz_escores.txt", "w");
    if (fp == NULL) {
        printf("Erro ao abrir o arquivo de saída.\n");
        exit(1);
    }

    int lin, col;

    fprintf(fp, "\nMatriz de escores:\n");

    for (lin = 0; lin <= tamSeqMenor; lin++) {
        for (col = 0; col <= tamSeqMaior; col++) {
            fprintf(fp, "%4d", matrizEscores[lin][col]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Matriz de escores gravada em matriz_escores.txt\n");
}

void traceBack(int tipo) {
    int tbLin, tbCol, peso, pos, aux, i;

    if (tipo == 1) {
        printf("\nGeracao do Primeiro Maior Alinhamento Global:\n");
        tbLin = linPMaior;
        tbCol = colPMaior;
    } else {
        printf("\nGeracao do Ultimo Maior Alinhamento Global:\n");
        tbLin = linUMaior;
        tbCol = colUMaior;
    }

    pos = 0;
    do {
        peso = matrizPesos[seqMenor[tbLin - 1]][seqMaior[tbCol - 1]];
        escoreDiag = matrizEscores[tbLin - 1][tbCol - 1] + peso;
        escoreLin = matrizEscores[tbLin][tbCol - 1] - penalGap;
        escoreCol = matrizEscores[tbLin - 1][tbCol] - penalGap;

        if ((escoreDiag >= escoreLin) && (escoreDiag >= escoreCol)) {
            if (seqMenor[tbLin - 1] != seqMaior[tbCol - 1]) {
                alinhaGMenor[pos] = X;
                alinhaGMaior[pos] = seqMaior[tbCol - 1];
                tbCol--;
                pos++;
                alinhaGMenor[pos] = seqMenor[tbLin - 1];
                alinhaGMaior[pos] = X;
                tbLin--;
                pos++;
            } else {
                alinhaGMenor[pos] = seqMenor[tbLin - 1];
                tbLin--;
                alinhaGMaior[pos] = seqMaior[tbCol - 1];
                tbCol--;
                pos++;
            }
        } else if (escoreLin >= escoreCol) {
            alinhaGMenor[pos] = X;
            alinhaGMaior[pos] = seqMaior[tbCol - 1];
            tbCol--;
            pos++;
        } else {
            alinhaGMenor[pos] = seqMenor[tbLin - 1];
            alinhaGMaior[pos] = X;
            tbLin--;
            pos++;
        }
    } while ((tbLin != 0) && (tbCol != 0));

    while (tbLin > 0) {
        alinhaGMenor[pos] = seqMenor[tbLin - 1];
        alinhaGMaior[pos] = X;
        tbLin--;
        pos++;
    }

    while (tbCol > 0) {
        alinhaGMenor[pos] = X;
        alinhaGMaior[pos] = seqMaior[tbCol - 1];
        tbCol--;
        pos++;
    }

    tamAlinha = pos;

    for (i = 0; i < (tamAlinha / 2); i++) {
        aux = alinhaGMenor[i];
        alinhaGMenor[i] = alinhaGMenor[tamAlinha - i - 1];
        alinhaGMenor[tamAlinha - i - 1] = aux;

        aux = alinhaGMaior[i];
        alinhaGMaior[i] = alinhaGMaior[tamAlinha - i - 1];
        alinhaGMaior[tamAlinha - i - 1] = aux;
    }

    printf("\nAlinhamento Global Gerado.");
}

int menuOpcao(void) {
    int op;
    char enter;

    do {
        printf("\nMenu de Opcao:");
        printf("\n<01> Ler Matriz de Pesos");
        printf("\n<02> Mostrar Matriz de Pesos");
        printf("\n<03> Ler Penalidade de Gap");
        printf("\n<04> Mostrar Penalidade");
        printf("\n<05> Definir Sequencias Genomicas");
        printf("\n<06> Mostrar Sequencias");
        printf("\n<07> Definir Tamanho do Bloco");
        printf("\n<08> Gerar Matriz de Escores");
        printf("\n<09> Mostrar Matriz de Escores");
        printf("\n<10> Gerar Alinhamento Global");
        printf("\n<11> Mostrar Alinhamento Global");
        printf("\n<12> Sair");
        printf("\nDigite a opcao => ");
        scanf("%d", &op);
        scanf("%c", &enter);  // Remove o enter
    } while ((op < 1) || (op > sair));

    return (op);
}

void trataOpcao(int op, int rank, int size) {
    int resp;
    char enter;
    char file1[100], file2[100];

    switch (op) {
        case 1:
            leMatrizPesos();
            break;
        case 2:
            mostraMatrizPesos();
            break;
        case 3:
            penalGap = lePenalidade();
            break;
        case 4:
            printf("\nPenalidade = %d", penalGap);
            break;
        case 5:
            printf("\nDeseja Definicao: <1>MANUAL, <2>ALEATORIA, <3>ARQUIVO? = ");
            scanf("%d", &resp);
            scanf("%c", &enter); /* remove o enter */
            if (resp == 1) {
                leSequenciasTeclado();
            } else if (resp == 2) {
                leTamMaior();
                leTamMenor();
                grauMuta = leGrauMutacao();
                geraSequenciasAleatorias();
            } else {
                printf("\nDigite o nome do arquivo para a sequencia 1: ");
                scanf("%s", file1);
                printf("\nDigite o nome do arquivo para a sequencia 2: ");
                scanf("%s", file2);
                leSequenciasArquivo(file1, file2);
            }
            break;
        case 6:
            mostraSequencias();
            break;
        case 7:
            printf("\nDigite o tamanho do bloco para transmissão MPI: ");
            scanf("%d", &blocoTamanho);
            break;
        case 8:
            geraMatrizEscoresMPI(rank, size, blocoTamanho);
            if (rank == 0) gravaMatrizEscoresEmArquivo();
            break;
        case 9:
            if (rank == 0) mostraMatrizEscores();
            break;
        case 10:
            if (rank == 0) {
                printf("\nDeseja: <1> Primeiro Maior ou <2> Ultimo Maior? = ");
                scanf("%d", &resp);
                scanf("%c", &enter); /* remove o enter */
                traceBack(resp);
            }
            break;
        case 11:
            if (rank == 0) mostraAlinhamentoGlobal();
            break;
    }
}

void main(int argc, char *argv[]) {
    int opcao, rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL));

    do {
        if (rank == 0) {
            printf("\n\nPrograma Needleman-Wunsch com MPI\n");
        }
        opcao = menuOpcao();
        MPI_Bcast(&opcao, 1, MPI_INT, 0, MPI_COMM_WORLD);
        trataOpcao(opcao, rank, size);
    } while (opcao != sair);

    MPI_Finalize();
}
