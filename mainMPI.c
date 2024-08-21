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
#define max_processos 5
#define max_bloc_size 1000

char baseMapa[5] = {'A', 'T', 'G', 'C', '-'};

char seqMaior[maxSeq], seqMenor[maxSeq], alinhaMaior[maxSeq], alinhaMenor[maxSeq];

int matrizScores[maxSeq + 1][maxSeq + 1];

int tamSeqMaior = 4, tamSeqMenor = 3, tamAlinha, penalGap = 0, grauMuta = 0, diagScore, linScore, colScore, k = 1;
int matrizPesos[4][4] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

int alinhamentoScores[maxSeq];

int rank, np, blocoTamanho;

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
void geraMatrizScoresMPI(int np, int blocoTamanho);
void traceBackMPI();
void mostraAlinhamentoGlobal(void);
void salvaMatrizScores(void);
void mostraMatrizScores();
int menuOpcao(void);
void trataOpcao(int op);

int max(int a, int b)
{
    return (a > b) ? a : b;
}

int min(int a, int b)
{
    return (a < b) ? a : b;
}

// Implementação das funções
int leTamMaior(void)
{
    printf("\nLeitura do Tamanho da Sequencia Maior:");
    do
    {
        printf("\nDigite 0 < valor < %d = ", maxSeq);
        fflush(stdout);
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
        fflush(stdout);
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
        fflush(stdout);
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
            fflush(stdout);
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
        fflush(stdout);
        scanf("%d", &prob);
    } while ((prob < 0) || (prob > 100));

    return prob;
}

void leSequenciasDeArquivo()
{
    FILE *fileMaior, *fileMenor;
    char nomeArquivoMaior[100], nomeArquivoMenor[100];

    printf("\nDigite o nome do arquivo para a Sequencia Maior: ");
    fflush(stdout);
    scanf("%s", nomeArquivoMaior);
    fileMaior = fopen(nomeArquivoMaior, "r");
    if (fileMaior == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", nomeArquivoMaior);
        exit(1);
    }

    printf("\nDigite o nome do arquivo para a Sequencia Menor: ");
    fflush(stdout);
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

    for (int i = 0; i < tamSeqMenor; i++)
    {
        int probAux = rand() % 100;
        if (probAux < grauMuta)
        {
            seqMenor[i] = (seqMenor[i] + (rand() % 3) + 1) % 4;
        }
    }

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
    printf("Digite sua opcao ==> ");
    fflush(stdout);
    scanf("%d", &opcao);

    // Solicitar o grau de mutação
    grauMuta = leGrauMutacao();

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

        // Aplicar mutações com base no grau de mutação
        for (i = 0; i < tamSeqMenor; i++)
        {
            int probAux = rand() % 100;
            if (probAux < grauMuta)
            {
                seqMenor[i] = (seqMenor[i] + (rand() % 3) + 1) % 4;
            }
        }
    }
    else if (opcao == 2)
    {
        leSequenciasDeArquivo();

        // Aplicar mutações com base no grau de mutação
        for (i = 0; i < tamSeqMenor; i++)
        {
            int probAux = rand() % 100;
            if (probAux < grauMuta)
            {
                seqMenor[i] = (seqMenor[i] + (rand() % 3) + 1) % 4;
            }
        }
    }
    else if (opcao == 3)
    {
        leTamMaior();
        leTamMenor();
        geraSequencias(); // Aqui a mutação já está aplicada dentro da função
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
    int i;

    printf("\nAlinhamento Obtido - Tamanho = %d:\n", tamAlinha);

    printf("%c", baseMapa[alinhaMaior[0]]);
    for (i = 1; i < tamAlinha; i++)
        printf("%c", baseMapa[alinhaMaior[i]]);
    printf("\n");

    printf("%c", baseMapa[alinhaMenor[0]]);
    for (i = 1; i < tamAlinha; i++)
        printf("%c", baseMapa[alinhaMenor[i]]);
    printf("\n");
}

void salvaMatrizScores(void)
{
    FILE *file;
    char nomeArquivo[100];
    int i, lin, col;

    printf("\nDigite o nome do arquivo para salvar a Matriz de Scores: ");
    fflush(stdout);
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

void geraMatrizScoresMPI(int np, int blocoTamanho) {
    int i, j;
    int linhasPorProcesso, linhasRestantes;
    int startLinha, endLinha;
    int numLinhas = tamSeqMenor + 1; // Total de linhas da matriz
    MPI_Status status;

    // Calcular número de linhas por processo
    linhasPorProcesso = numLinhas / np;
    linhasRestantes = numLinhas % np;

    // Distribuição das linhas entre os processos
    startLinha = rank * linhasPorProcesso;
    endLinha = startLinha + linhasPorProcesso - 1;

    if (rank < linhasRestantes) {
        startLinha += rank;
        endLinha += rank + 1;
    } else {
        startLinha += linhasRestantes;
        endLinha += linhasRestantes;
    }

    // Ajustar o range para o último processo
    if (rank == np - 1) {
        endLinha = numLinhas - 1;
    }

    // Inicializar a matriz local de scores para o processo atual
    int localMatrizScores[tamSeqMenor + 1][tamSeqMaior + 1];

    // Preencher a primeira linha da matriz
    if (rank == 0) {
        for (j = 0; j <= tamSeqMaior; j++) {
            matrizScores[0][j] = j * penalGap;
        }
        // Enviar a primeira linha para todos os outros processos e o tamanho do bloco
        for (i = 1; i < np; i++) {
            MPI_Send(&blocoTamanho, sizeof(int), MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(matrizScores[0], tamSeqMaior + 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    } else {
        // Receber a primeira linha da matriz
        MPI_Recv(matrizScores[0], tamSeqMaior + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    // Preencher a matriz de scores local para o processo atual
    for (i = startLinha; i <= endLinha; i++) {
        for (j = 1; j <= tamSeqMaior; j++) {
            if (i == 0) {
                localMatrizScores[i][j] = j * penalGap;
            } else {
                localMatrizScores[i][j] = max(localMatrizScores[i - 1][j] + penalGap,
                                              localMatrizScores[i][j - 1] + penalGap);

                if (seqMaior[j - 1] == seqMenor[i - 1]) {
                    diagScore = matrizPesos[seqMenor[i - 1]][seqMaior[j - 1]];
                } else {
                    diagScore = 0;
                }
                localMatrizScores[i][j] = max(localMatrizScores[i][j],
                                              localMatrizScores[i - 1][j - 1] + diagScore);
            }
        }
    }

    // Enviar as linhas calculadas para o processo 0
    if (rank != 0) {
        for (i = startLinha; i <= endLinha; i++) {
            MPI_Send(localMatrizScores[i], tamSeqMaior + 1, MPI_INT, 0, i, MPI_COMM_WORLD);
        }
    }

    // Receber as linhas calculadas de todos os processos
    if (rank == 0) {
        for (i = 1; i < np; i++) {
            int startLinhaRecebida = i * linhasPorProcesso;
            int endLinhaRecebida = startLinhaRecebida + linhasPorProcesso - 1;

            if (i < linhasRestantes) {
                startLinhaRecebida += i;
                endLinhaRecebida += i + 1;
            } else {
                startLinhaRecebida += linhasRestantes;
                endLinhaRecebida += linhasRestantes;
            }

            if (i == np - 1) {
                endLinhaRecebida = numLinhas - 1;
            }

            for (j = startLinhaRecebida; j <= endLinhaRecebida; j++) {
                MPI_Recv(matrizScores[j], tamSeqMaior + 1, MPI_INT, i, j, MPI_COMM_WORLD, &status);
            }
        }
    }
}


// Traceback
void traceBackMPI()
{
    int lin = tamSeqMenor;
    int col = tamSeqMaior;
    tamAlinha = 0;

    while (lin > 0 || col > 0)
    {
        if (lin > 0 && col > 0 && matrizScores[lin][col] == matrizScores[lin - 1][col - 1] + matrizPesos[seqMenor[lin - 1]][seqMaior[col - 1]])
        {
            alinhaMaior[tamAlinha] = seqMaior[col - 1];
            alinhaMenor[tamAlinha] = seqMenor[lin - 1];
            lin--;
            col--;
        }
        else if (lin > 0 && matrizScores[lin][col] == matrizScores[lin - 1][col] + penalGap)
        {
            alinhaMaior[tamAlinha] = 4;
            alinhaMenor[tamAlinha] = seqMenor[lin - 1];
            lin--;
        }
        else
        {
            alinhaMaior[tamAlinha] = seqMaior[col - 1];
            alinhaMenor[tamAlinha] = 4;
            col--;
        }
        tamAlinha++;
    }

    for (int i = 0; i < tamAlinha / 2; i++)
    {
        int temp = alinhaMaior[i];
        alinhaMaior[i] = alinhaMaior[tamAlinha - 1 - i];
        alinhaMaior[tamAlinha - 1 - i] = temp;

        temp = alinhaMenor[i];
        alinhaMenor[i] = alinhaMenor[tamAlinha - 1 - i];
        alinhaMenor[tamAlinha - 1 - i] = temp;
    }
}
int menuOpcao(void)
{
    int op;
    char enter;

    do
    {
        printf("\nMenu de Opcao:");
        printf("\n<01> Definir o tamanho do Bloco");
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
        fflush(stdout);
        scanf("%d", &op);
        fflush(stdout);
        scanf("%c", &enter);
    } while ((op < 1) || (op > 13));

    return (op);
}

void trataOpcao(int op)
{
    int resp;
    char enter;

    switch (op)
    {
    case 1:
        printf("\nDigite o tamanho do bloco: ");
        fflush(stdout);
        scanf("%d", &blocoTamanho);
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
        leSequencias();
        break;
    case 7:
        mostraSequencias();
        break;
    case 8:
        printf("Número de processos: %d, Tamanho do bloco: %d\n", np, blocoTamanho); // Verifique os valores
        geraMatrizScoresMPI(np, blocoTamanho);
        break;
    case 9:
        mostraMatrizScores();
        break;

    case 10:
        salvaMatrizScores();
        break;

    case 11:
        traceBackMPI();
        break;

    case 12:
        mostraAlinhamentoGlobal();
        break;

    case 13:
        exit(0);
    }
}

int main(int argc, char *argv[])
{
    int opcao;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Status status;

    if (rank == 0)
    {
        printf("Numero total de processos: %d\n", np);
        do
        {
            printf("\n\nPrograma Needleman-Wunsch Paralelo com MPI\n");
            opcao = menuOpcao();
            trataOpcao(opcao);
        } while (opcao != 13);
    }else {
        int tamRecv;
        MPI_Recv(&tamRecv, sizeof(int), MPI_INT, 0, rank, MPI_COMM_WORLD, &status);
        geraMatrizScoresMPI(rank, tamRecv);
    }
   
    MPI_Finalize();
    return 0;
}
