#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define A 0 // A
#define T 1
#define G 2
#define C 3
#define sair 14

#define maxSeq 10000

char baseMapa[5] = {'A', 'T', 'G', 'C', '-'};

char seqMaior[maxSeq], seqMenor[maxSeq], alinhaMaior[maxSeq][maxSeq], alinhaMenor[maxSeq][maxSeq];

int matrizScores[maxSeq + 1][maxSeq + 1];

int tamSeqMaior = 6, tamSeqMenor = 6, tamAlinha[maxSeq], penalGap = 0, grauMuta = 0, diagScore, linScore, colScore, k = 1;

int matrizPesos[4][4] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

int alinhamentoScores[maxSeq];

typedef struct
{
    int tbLin;
    int tbCol;
    int pos;
    int rank;
} traceback_data_t;

// Declaração de funções
int leTamMaior(void);
int leTamMenor(void);
int lePenalidade(void);
void leMatrizPesos();
void mostraMatrizPesos(void);
int leGrauMutacao(void);
void leSequenciasDeArquivo();
void leSequencias();
void geraSequencias();
void mostraSequencias(void);
void mostraAlinhamentoGlobal(void);
void salvaMatrizScores(void);
void mostraMatrizScores();
void preencheMatrizScores(int rank, int size);
void traceBack(int rank, int size, int k);
int leNumeroDeAlinhamentos(void);
int leNumeroDeThreads(void);
int menuOpcao(void);
void trataOpcao(int op, int rank, int size);

// Implementação das funções
int leTamMaior(void)
{
    printf("\nLeitura do Tamanho da Sequencia Maior:");
    do
    {
        printf("\nDigite 0 < valor < %d = ", maxSeq);
        scanf("%d", &tamSeqMaior);
    } while ((tamSeqMaior < 1) || (tamSeqMaior > maxSeq));
    return tamSeqMaior;
}

int leTamMenor(void)
{
    printf("\nLeitura do Tamanho da Sequencia Menor:");
    do
    {
        printf("\nDigite 0 < valor <= %d = ", tamSeqMaior);
        scanf("%d", &tamSeqMenor);
    } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));
    return tamSeqMenor;
}

int lePenalidade(void)
{
    int penal;

    printf("\nLeitura da Penalidade de Gap:");
    do
    {
        printf("\nDigite valor >= 0 = ");
        scanf("%d", &penal);
    } while (penal < 0);

    return penal;
}

void leMatrizPesos()
{
    int i, j;

    printf("\nLeitura da Matriz de Pesos:\n");
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            printf("Digite valor %c x %c = ", baseMapa[i], baseMapa[j]);
            scanf("%d", &(matrizPesos[i][j]));
        }
        printf("\n");
    }
}

void mostraMatrizPesos(void)
{
    int i, j;

    printf("\nMatriz de Pesos Atual:");
    printf("\n%4c%4c%4c%4c%4c\n", ' ', 'A', 'T', 'G', 'C');
    for (i = 0; i < 4; i++)
    {
        printf("%4c", baseMapa[i]);
        for (j = 0; j < 4; j++)
            printf("%4d", matrizPesos[i][j]);
        printf("\n");
    }
}

int leGrauMutacao(void)
{
    int prob;

    printf("\nLeitura da Porcentagem Maxima de Mutacao Aleatoria:\n");
    do
    {
        printf("\nDigite 0 <= valor <= 100 = ");
        scanf("%d", &prob);
    } while ((prob < 0) || (prob > 100));

    return prob;
}

void leSequenciasDeArquivo()
{
    FILE *fileMaior, *fileMenor;
    char nomeArquivoMaior[100], nomeArquivoMenor[100];

    printf("\nDigite o nome do arquivo para a Sequencia Maior: ");
    scanf("%s", nomeArquivoMaior);
    fileMaior = fopen(nomeArquivoMaior, "r");
    if (fileMaior == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", nomeArquivoMaior);
        exit(1);
    }

    printf("\nDigite o nome do arquivo para a Sequencia Menor: ");
    scanf("%s", nomeArquivoMenor);
    fileMenor = fopen(nomeArquivoMenor, "r");
    if (fileMenor == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", nomeArquivoMenor);
        exit(1);
    }

    int oldTamMaior = tamSeqMaior;
    int oldTamMenor = tamSeqMenor;
    char oldMaior[maxSeq];
    char oldMenor[maxSeq];
    memcpy(oldMaior, seqMaior, maxSeq);
    memcpy(oldMenor, seqMenor, maxSeq);

    tamSeqMaior = fread(seqMaior, sizeof(char), maxSeq, fileMaior);
    tamSeqMenor = fread(seqMenor, sizeof(char), maxSeq, fileMenor);

    for (int i = 0; i < tamSeqMaior; i++)
    {
        if (seqMaior[i] == 'A')
        {
            seqMaior[i] = 0;
        }
        else if (seqMaior[i] == 'T')
        {
            seqMaior[i] = 1;
        }
        else if (seqMaior[i] == 'G')
        {
            seqMaior[i] = 2;
        }
        else if (seqMaior[i] == 'C')
        {
            seqMaior[i] = 3;
        }
        else
        {
            tamSeqMaior = oldTamMaior;
            memcpy(seqMaior, oldMaior, maxSeq);
            printf("Sequencia maior invalida, valor original restaurado.\n");
        }
    }

    for (int i = 0; i < tamSeqMenor; i++)
    {
        if (seqMenor[i] == 'A')
        {
            seqMenor[i] = 0;
        }
        else if (seqMenor[i] == 'T')
        {
            seqMenor[i] = 1;
        }
        else if (seqMenor[i] == 'G')
        {
            seqMenor[i] = 2;
        }
        else if (seqMenor[i] == 'C')
        {
            seqMenor[i] = 3;
        }
        else
        {
            tamSeqMenor = oldTamMenor;
            memcpy(seqMenor, oldMenor, maxSeq);
            printf("Sequencia menor invalida, valor original restaurado. \n");
        }
    }

    fclose(fileMaior);
    fclose(fileMenor);

    if (tamSeqMenor > tamSeqMaior)
    {
        printf("Erro: Sequencia menor é maior que a sequencia maior.\n");
        exit(1);
    }
}

void leSequencias()
{
    int i, erro, opcao;

    printf("\nSelecione o modo de entrada das sequencias:\n");
    printf("1. Manual\n");
    printf("2. Arquivo\n");
    printf("3. Aleatoria\n");
    scanf("%d", &opcao);

    if (opcao == 1)
    {
        printf("\nLeitura das Sequencias:\n");
        do
        {
            printf("\nPara a Sequencia Maior,");
            printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
            do
            {
                printf("\n> ");
                fgets(seqMaior, maxSeq, stdin);
                tamSeqMaior = strlen(seqMaior) - 1;
            } while (tamSeqMaior < 1);
            printf("\ntamSeqMaior = %d\n", tamSeqMaior);
            i = 0;
            erro = 0;
            do
            {
                switch (seqMaior[i])
                {
                case 'A':
                    seqMaior[i] = (char)A;
                    break;
                case 'T':
                    seqMaior[i] = (char)T;
                    break;
                case 'G':
                    seqMaior[i] = (char)G;
                    break;
                case 'C':
                    seqMaior[i] = (char)C;
                    break;
                default:
                    erro = 1;
                }
                i++;
            } while ((erro == 0) && (i < tamSeqMaior));
        } while (erro == 1);

        do
        {
            printf("\nPara a Sequencia Menor, ");
            printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
            do
            {
                printf("\n> ");
                fgets(seqMenor, maxSeq, stdin);
                tamSeqMenor = strlen(seqMenor) - 1;
            } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));
            printf("\ntamSeqMenor = %d\n", tamSeqMenor);

            i = 0;
            erro = 0;
            do
            {
                switch (seqMenor[i])
                {
                case 'A':
                    seqMenor[i] = (char)A;
                    break;
                case 'T':
                    seqMenor[i] = (char)T;
                    break;
                case 'G':
                    seqMenor[i] = (char)G;
                    break;
                case 'C':
                    seqMenor[i] = (char)C;
                    break;
                default:
                    erro = 1;
                }
                i++;
            } while ((erro == 0) && (i < tamSeqMenor));
        } while (erro == 1);
    }
    else if (opcao == 2)
    {
        leSequenciasDeArquivo();
    }
    else if (opcao == 3)
    {
        leTamMaior();
        leTamMenor();
        grauMuta = leGrauMutacao();
        geraSequencias();
    }
    else
    {
        printf("Opcao invalida!\n");
        exit(1);
    }
}

void geraSequencias()
{
    int i, dif, probAux, ind, nTrocas;
    char base;

    srand(time(NULL));

    printf("\nGeracao Aleatoria das Sequencias:\n");

    for (i = 0; i < tamSeqMaior; i++)
    {
        base = (char)(rand() % 4);
        seqMaior[i] = base;
    }

    dif = tamSeqMaior - tamSeqMenor;
    ind = 0;
    if (dif > 0)
        ind = rand() % dif;

    for (i = 0; i < tamSeqMenor; i++)
        seqMenor[i] = seqMaior[ind + i];

    i = 0;
    nTrocas = 0;
    while ((i < tamSeqMenor) && (nTrocas < grauMuta))
    {
        probAux = rand() % 100;

        if (probAux < grauMuta)
        {
            seqMenor[i] = (seqMenor[i] + (rand() % 3) + 1) % 4;
            nTrocas++;
        }
    }

    printf("\nSequencias Geradas, Dif = %d Ind = %d\n", dif, ind);
}

void mostraSequencias(void)
{
    int i;

    printf("\nSequencias Atuais:\n");
    printf("\nSequencia Maior, Tam = %d\n", tamSeqMaior);
    printf("%c", baseMapa[(int)seqMaior[0]]);
    for (i = 1; i < tamSeqMaior; i++)
        printf("%c", baseMapa[(int)seqMaior[i]]);
    printf("\n");

    printf("\nSequencia Menor, Tam = %d\n", tamSeqMenor);
    printf("%c", baseMapa[(int)seqMenor[0]]);
    for (i = 1; i < tamSeqMenor; i++)
        printf("%c", baseMapa[(int)seqMenor[i]]);
    printf("\n");
}

void mostraAlinhamentoGlobal(void)
{
    int i, j;
    int bestProcess = 0;
    int bestScore = alinhamentoScores[0];

    printf("\nAlinhamentos Atuais - Mostrando %d Alinhamentos:\n", k);

    for (j = 0; j < k; j++)
    {
        printf("\nProcesso %d -> Alinhamento %d (Score: %d):\n", j, j + 1, alinhamentoScores[j]);
        printf("%c", baseMapa[(int)alinhaMaior[j][0]]);
        for (i = 1; i < tamAlinha[j]; i++)
            printf("%c", baseMapa[(int)alinhaMaior[j][i]]);
        printf("\n");

        printf("%c", baseMapa[(int)alinhaMenor[j][0]]);
        for (i = 1; i < tamAlinha[j]; i++)
            printf("%c", baseMapa[(int)alinhaMenor[j][i]]);
        printf("\n");

        if (alinhamentoScores[j] > bestScore)
        {
            bestScore = alinhamentoScores[j];
            bestProcess = j;
        }
    }

    printf("\nMelhor Alinhamento (Processo %d - Score: %d):\n", bestProcess, bestScore);
    printf("%c", baseMapa[(int)alinhaMaior[bestProcess][0]]);
    for (i = 1; i < tamAlinha[bestProcess]; i++)
        printf("%c", baseMapa[(int)alinhaMaior[bestProcess][i]]);
    printf("\n");

    printf("%c", baseMapa[(int)alinhaMenor[bestProcess][0]]);
    for (i = 1; i < tamAlinha[bestProcess]; i++)
        printf("%c", baseMapa[(int)alinhaMenor[bestProcess][i]]);
    printf("\n");
}

void salvaMatrizScores(void)
{
    FILE *file;
    char nomeArquivo[100];
    int i, lin, col;

    printf("\nDigite o nome do arquivo para salvar a Matriz de Scores: ");
    scanf("%s", nomeArquivo);
    file = fopen(nomeArquivo, "w");
    if (file == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", nomeArquivo);
        exit(1);
    }

    fprintf(file, "Matriz de Scores:\n");
    fprintf(file, "%4c%4c", ' ', ' ');
    for (i = 0; i <= tamSeqMaior; i++)
        fprintf(file, "%4d", i);
    fprintf(file, "\n");

    fprintf(file, "%4c%4c%4c", ' ', ' ', '-');
    for (i = 0; i < tamSeqMaior; i++)
        fprintf(file, "%4c", baseMapa[(int)(seqMaior[i])]);
    fprintf(file, "\n");

    fprintf(file, "%4c%4c", '0', '-');
    for (col = 0; col <= tamSeqMaior; col++)
        fprintf(file, "%4d", matrizScores[0][col]);
    fprintf(file, "\n");

    for (lin = 1; lin <= tamSeqMenor; lin++)
    {
        fprintf(file, "%4d%4c", lin, baseMapa[(int)(seqMenor[lin - 1])]);
        for (col = 0; col <= tamSeqMaior; col++)
        {
            fprintf(file, "%4d", matrizScores[lin][col]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void mostraMatrizScores()
{
    int i, lin, col;

    printf("\nMatriz de Scores Atual:\n");

    printf("%4c%4c", ' ', ' ');
    for (i = 0; i <= tamSeqMaior; i++)
        printf("%4d", i);
    printf("\n");

    printf("%4c%4c%4c", ' ', ' ', '-');
    for (i = 0; i < tamSeqMaior; i++)
        printf("%4c", baseMapa[(int)(seqMaior[i])]);
    printf("\n");

    printf("%4c%4c", '0', '-');
    for (col = 0; col <= tamSeqMaior; col++)
        printf("%4d", matrizScores[0][col]);
    printf("\n");

    for (lin = 1; lin <= tamSeqMenor; lin++)
    {
        printf("%4d%4c", lin, baseMapa[(int)(seqMenor[lin - 1])]);
        for (col = 0; col <= tamSeqMaior; col++)
        {
            printf("%4d", matrizScores[lin][col]);
        }
        printf("\n");
    }
}

void preencheMatrizScores(int rank, int size)
{
    int lin, col, peso;
    int blockSize = 10; // Defina um tamanho de bloco apropriado
    int linesPerProcess = (tamSeqMenor + size - 1) / size; // Divisão equidistante

    for (lin = rank; lin <= tamSeqMenor; lin += size)
    {
        for (col = 1; col <= tamSeqMaior; col++)
        {
            peso = matrizPesos[(int)(seqMenor[lin - 1])][(int)(seqMaior[col - 1])];
            diagScore = matrizScores[lin - 1][col - 1] + peso;
            linScore = matrizScores[lin][col - 1] - penalGap;
            colScore = matrizScores[lin - 1][col] - penalGap;

            if ((diagScore >= linScore) && (diagScore >= colScore))
                matrizScores[lin][col] = diagScore;
            else if (linScore > colScore)
                matrizScores[lin][col] = linScore;
            else
                matrizScores[lin][col] = colScore;
        }

        if (rank != 0)
        {
            for (int i = 0; i <= tamSeqMaior; i += blockSize)
            {
                MPI_Send(&matrizScores[lin][i], blockSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    if (rank == 0)
    {
        for (lin = 1; lin <= tamSeqMenor; lin++)
        {
            if (lin % size != 0)
            {
                for (int i = 0; i <= tamSeqMaior; i += blockSize)
                {
                    MPI_Recv(&matrizScores[lin][i], blockSize, MPI_INT, lin % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }
}


void traceBack(int rank, int size, int k)
{
    int tbLin = tamSeqMenor;
    int tbCol = tamSeqMaior;
    int pos = 0;
    int bestScore = 0;

    while (tbLin > 0 && tbCol > 0)
    {
        int peso = matrizPesos[(int)(seqMenor[tbLin - 1])][(int)(seqMaior[tbCol - 1])];
        diagScore = matrizScores[tbLin - 1][tbCol - 1] + peso;
        linScore = matrizScores[tbLin][tbCol - 1] - penalGap;
        colScore = matrizScores[tbLin - 1][tbCol] - penalGap;

        if ((diagScore >= linScore) && (diagScore >= colScore))
        {
            alinhaMenor[rank][pos] = seqMenor[tbLin - 1];
            alinhaMaior[rank][pos] = seqMaior[tbCol - 1];
            tbLin--;
            tbCol--;
            pos++;
            bestScore += peso;
        }
        else if (linScore > colScore)
        {
            alinhaMenor[rank][pos] = (char)4;
            alinhaMaior[rank][pos] = seqMaior[tbCol - 1];
            tbCol--;
            pos++;
            bestScore -= penalGap;
        }
        else
        {
            alinhaMenor[rank][pos] = seqMenor[tbLin - 1];
            alinhaMaior[rank][pos] = (char)4;
            tbLin--;
            pos++;
            bestScore -= penalGap;
        }
    }

    while (tbLin > 0)
    {
        alinhaMenor[rank][pos] = seqMenor[tbLin - 1];
        alinhaMaior[rank][pos] = (char)4;
        tbLin--;
        pos++;
        bestScore -= penalGap;
    }

    while (tbCol > 0)
    {
        alinhaMenor[rank][pos] = (char)4;
        alinhaMaior[rank][pos] = seqMaior[tbCol - 1];
        tbCol--;
        pos++;
        bestScore -= penalGap;
    }

    tamAlinha[rank] = pos;
    alinhamentoScores[rank] = bestScore;

    if (rank != 0)
    {
        MPI_Send(&alinhamentoScores[rank], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&tamAlinha[rank], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&alinhaMaior[rank][0], tamAlinha[rank], MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&alinhaMenor[rank][0], tamAlinha[rank], MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        for (int i = 1; i < k; i++)
        {
            MPI_Recv(&alinhamentoScores[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&tamAlinha[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&alinhaMaior[i][0], tamAlinha[i], MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&alinhaMenor[i][0], tamAlinha[i], MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}


int leNumeroDeAlinhamentos(void)
{
    int k;
    printf("\nLeitura do Numero de Alinhamentos a Mostrar (k):\n");
    do
    {
        printf("\nDigite valor > 0: ");
        scanf("%d", &k);
    } while (k <= 0);
    return k;
}

int menuOpcao(void)
{
    int op;
    char enter;

    do
    {
        printf("\nMenu de Opcao:");
        printf("\n<01> Definir o Numero de Threads");
        printf("\n<02> Ler Matriz de Pesos");
        printf("\n<03> Mostrar Matriz de Pesos");
        printf("\n<04> Ler Penalidade de Gap");
        printf("\n<05> Mostrar Penalidade");
        printf("\n<06> Definir Sequencias Genomicas");
        printf("\n<07> Mostrar Sequencias");
        printf("\n<08> Gerar Matriz de Scores");
        printf("\n<09> Mostrar Matriz de Scores");
        printf("\n<10> Salvar Matriz de Scores");
        printf("\n<11> Gerar Alinhamento Global");
        printf("\n<12> Mostrar Alinhamento Global");
        printf("\n<13> Sair");
        printf("\nDigite a opcao => ");
        scanf("%d", &op);
        scanf("%c", &enter);
    } while ((op < 1) || (op > 13));

    return (op);
}

void trataOpcao(int op, int rank, int size)
{
    int resp;
    char enter;

    switch (op)
    {
    case 1:
        //numThreads = leNumeroDeThreads();
        break;
    case 2:
        leMatrizPesos();
        break;
    case 3:
        mostraMatrizPesos();
        break;
    case 4:
        penalGap = lePenalidade();
        break;
    case 5:
        printf("\nPenalidade = %d", penalGap);
        break;
    case 6:
        printf("\nDeseja Definicao: <1>MANUAL, <2>ARQUIVO ou <3>ALEATORIA? = ");
        scanf("%d", &resp);
        scanf("%c", &enter); /* remove o enter */
        if (resp == 1)
        {
            leSequencias();
        }
        else if (resp == 2)
        {
            leSequenciasDeArquivo();
        }
        else if (resp == 3)
        {
            leTamMaior();
            leTamMenor();
            grauMuta = leGrauMutacao();
            geraSequencias();
        }
        break;
    case 7:
        mostraSequencias();
        break;
    case 8:
        preencheMatrizScores(rank, size);
        break;
    case 9:
        if (rank == 0)
            mostraMatrizScores();
        break;

    case 10:
        if (rank == 0)
            salvaMatrizScores();
        break;

    case 11:
        k = leNumeroDeAlinhamentos();
        traceBack(rank, size, k);
        break;

    case 12:
        if (rank == 0)
            mostraAlinhamentoGlobal();
        break;

    case 13:
        exit(0);
    }
}

int main(int argc, char **argv)
{
    int opcao;

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    do
    {
        if (rank == 0)
        {
            printf("\n\nPrograma Needleman-Wunsch MPI\n");
            opcao = menuOpcao();

            // Enviar a opção para todos os processos
            for (int i = 1; i < size; i++)
            {
                MPI_Send(&opcao, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        else
        {
            // Receber a opção do processo 0
            MPI_Recv(&opcao, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        trataOpcao(opcao, rank, size);
    } while (opcao != 14);

    MPI_Finalize();
    return 0;
}
