# Programação Concorrente - Análise de Genoma

Este trabalho consiste de 2 implementações com base no algoritmo de Needleman-Wunsch para análise de Genôma.

1 - Versão com Threads (necessário que o usuário tenha a biblioteca de Threads instalada no ambiente utilizado): 
        Para essa versão o algoritmo de Genôma foi paralelizado utilizando k Threads, a paralelização foi feita da seguinte forma:
        - O paralelismo ocorre em 2 etapas: 1) a construção/preenchimento da matriz de scores deve ser paralelizada entre k threads, de forma mais igualitária possível. Tal distribuição foi feita por linhas de forma que cada thread receba um conjunto de linhas parar gerar. Cada thread deve trabalhar com um conjunto de linhas equidistantes, não consecutivas, para promover a paralelização. Se necessário, deve ser utilizado mecanismos de sincronização; 2) a reconstrução das sequências alinhadas foi paralelizada em até k threads. Tal paralelização foi feita quando houver mais que uma possibilidade de traceback, ainda que com scores diferentes, com isso, o programa deverá mostrar até k pares de alinhamentos.
        - A paralelização ocorreu no preenchimento da matriz de scores e na reconstrução dos até k pares de alinhamentos possíveis. A versão paralela deve possuir as mesmas funcionalidades e características da versão sequencial também disponível.

Para executar o algoritmo basta utilizar os seguintes comandos:
```
$cc -lpthread main.c -o main
$./main
```

2 - Versão com MPI (necessário que o usuário tenha a biblioteca de MPI instalada no ambiente utilizado):
        Para essa versão o algoritmo de Genôma foi paralelizado utilizando MPI, a paralelização foi feita da seguinte forma:
        - O paralelismo foi realizado em apenas uma epata da aplicação: a construção/preenchimento da matriz de scores foi paralelizada entre np processos, de forma mais igualitária possível. Tal distribuição foi feita por linhas de forma que cada processo deve construir um conjunto de linhas. Cada processo deve trabalhar com um conjunto de linhas equidistantes, não consecutivass, para promover a paralelização. O traceback será único, identificando apenas um único alinhamento, executado pelo processo 0.
        - Para permitir o paralelismo na construção da matrriz de scores, na medida em que um processo vai gerando os scores de uma linha, ele vai transmitindo esta linha, bloco por bloco, para o próximo processo que precisa dela e para o processo 0 também. Esta transmissão(por mensagem MPI) foi feita por blocos, para não congestionar muito o subsistema de troca de mensagens. ATENÇÃO: Os processos não devem esperar a finalização de uma linha para somente depois poder transmiti-la, senão a execução será sequencializada.
        - O tamanho do bloco é definido pelo usuário. Como o processo 0 vai recebendo os blocos de todos os demais processos, ao final, ele terá a matriz de scores montada por completo e poderá assim executar o traceback. Nesse modelo de paralelização, o processo 0 não precisa trabalhar na paralelização, apenas no gerencimaneto da execução MPI.

Para executar o algortmo basta utilizar os seguintes comandos:
```
$mpicc mainMPI.c -o mainMPI
$mpirun --oversubscribe -np <numero de processos> ./mainMPI
```
O oversubscribe é para garantir que ele não passe dos processos disponíveis no seu computador.
