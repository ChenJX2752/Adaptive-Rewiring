#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// 数据结构定义
#define MAX_NODES 1000
#define Data_w "E:\\Adaptive\\Code\\dataput_new\\smt_6.txt"

typedef struct EdgeNode
{
    int node;
    struct EdgeNode *next;
} EdgeNode;

typedef struct
{
    EdgeNode **adjacencyList;
    int numEdges;
} Network;

typedef struct SimplexNode
{
    int node1;
    int node2;
    int node3;
    struct SimplexNode *next;
} SimplexNode;

typedef struct
{
    SimplexNode *simplices;
    int numSimplices;
} SimplexList;

typedef struct
{
    double *infectionDensity;
    double *averageDegree;
    double *averageSimplexCount;
} SimulationResult;

// 定义状态记录表格结构
typedef struct
{
    int checked; // 是否已经检查过
} InfectionPairRecord;

// 初始化感染记录表格
void initializeInfectionPairRecord(InfectionPairRecord *infectedPairs, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            infectedPairs[i * N + j].checked = 0;
        }
    }
}

// 边和单纯形的存储
void addEdge(Network *network, int u, int v)
{
    EdgeNode *newNode = (EdgeNode *)malloc(sizeof(EdgeNode));
    newNode->node = v;
    newNode->next = network->adjacencyList[u];
    network->adjacencyList[u] = newNode;
}

void addSimplex(SimplexList *simplexList, int u, int v, int w)
{
    SimplexNode *newSimplex = (SimplexNode *)malloc(sizeof(SimplexNode));
    newSimplex->node1 = u;
    newSimplex->node2 = v;
    newSimplex->node3 = w;
    newSimplex->next = simplexList->simplices;
    simplexList->simplices = newSimplex;
    simplexList->numSimplices++;
}

void readNetwork(const char *edgesFile, const char *simplicesFile, Network *network, SimplexList *simplexList)
{
    // 清空当前的网络和单纯形列表
    for (int i = 0; i < MAX_NODES; i++)
    {
        EdgeNode *temp = network->adjacencyList[i];
        while (temp != NULL)
        {
            EdgeNode *toFree = temp;
            temp = temp->next;
            free(toFree);
        }
        network->adjacencyList[i] = NULL;
    }
    network->numEdges = 0;

    SimplexNode *tempSimplex = simplexList->simplices;
    while (tempSimplex != NULL)
    {
        SimplexNode *toFree = tempSimplex;
        tempSimplex = tempSimplex->next;
        free(toFree);
    }
    simplexList->simplices = NULL;
    simplexList->numSimplices = 0;

    // Read simplices first
    FILE *file = fopen(simplicesFile, "r");
    if (!file)
    {
        fprintf(stderr, "Error opening simplices file\n");
        exit(1);
    }
    int u, v, w;
    while (fscanf(file, "%d %d %d", &u, &v, &w) != EOF)
    {
        addSimplex(simplexList, u, v, w);
    }
    fclose(file);
    // printf("%d\n", simplexList->numSimplices);
    //  Read edges
    file = fopen(edgesFile, "r");
    if (!file)
    {
        fprintf(stderr, "Error opening edges file\n");
        exit(1);
    }

    while (fscanf(file, "%d %d", &u, &v) != EOF)
    {
        // int inSimplex = 0;
        // SimplexNode *simplex = simplexList->simplices;
        // while (simplex != NULL)
        // {
        //     if ((simplex->node1 == u && simplex->node2 == v) || (simplex->node1 == v && simplex->node2 == u) ||
        //         (simplex->node1 == u && simplex->node3 == v) || (simplex->node1 == v && simplex->node3 == u) ||
        //         (simplex->node2 == u && simplex->node3 == v) || (simplex->node2 == v && simplex->node3 == u))
        //     {
        //         inSimplex = 1;
        //         break;
        //     }
        //     simplex = simplex->next;
        // }

        // if (!inSimplex)
        // {
        //     addEdge(network, u, v);
        //     addEdge(network, v, u);
        //     network->numEdges++;
        // }
        addEdge(network, u, v);
        addEdge(network, v, u);
        network->numEdges++;
    }
    fclose(file);
}

// 计算节点度分布
void calculateNodeDegreeDistribution(int *state, EdgeNode **adjacencyList, int N)
{
    int *degreeCount = (int *)calloc(N, sizeof(int));
    int maxDegree = 0;

    for (int i = 0; i < N; i++)
    {
        int degree = 0;
        EdgeNode *temp = adjacencyList[i];
        while (temp != NULL)
        {
            degree++;
            temp = temp->next;
        }
        degreeCount[degree]++;
        if (degree > maxDegree)
        {
            maxDegree = degree;
        }
    }

    int susceptibleCount = 0;
    int infectedCount = 0;
    for (int i = 0; i < N; i++)
    {
        if (state[i] == 0)
        {
            susceptibleCount++;
        }
        else if (state[i] == 1)
        {
            infectedCount++;
        }
    }

    printf("Degree distribution:\n");
    for (int i = 0; i <= maxDegree; i++)
    {
        if (degreeCount[i] > 0)
        {
            printf("Degree %d: %d nodes (S: %.2f%%, I: %.2f%%)\n", i, degreeCount[i],
                   100.0 * degreeCount[i] / susceptibleCount,
                   100.0 * degreeCount[i] / infectedCount);
        }
    }

    free(degreeCount);
}

// 辅助函数
int isConnected(EdgeNode **adjacencyList, int u, int v)
{
    EdgeNode *temp = adjacencyList[u];
    while (temp != NULL)
    {
        if (temp->node == v)
        {
            return 1;
        }
        temp = temp->next;
    }
    return 0;
}

// 辅助函数：检查三个节点是否在一个单纯形中
int isInSimplex(SimplexNode *simplexList, int u, int v, int w)
{
    SimplexNode *temp = simplexList;

    // 遍历单纯形链表，检查是否有包含 u, v, w 的单纯形
    while (temp != NULL)
    {
        // 检查单纯形的三个节点是否为 u, v, w （无序匹配）
        if ((temp->node1 == u || temp->node2 == u || temp->node3 == u) &&
            (temp->node1 == v || temp->node2 == v || temp->node3 == v) &&
            (temp->node1 == w || temp->node2 == w || temp->node3 == w))
        {
            return 1; // 找到包含 u, v, w 的单纯形，返回 1
        }
        temp = temp->next;
    }

    return 0; // 如果没有找到包含 u, v, w 的单纯形，返回 0
}

EdgeNode *removeEdge(EdgeNode *head, int node)
{
    EdgeNode *temp = head;
    EdgeNode *prev = NULL;

    while (temp != NULL)
    {
        if (temp->node == node)
        {
            if (prev == NULL)
            {
                EdgeNode *newHead = temp->next;
                free(temp);
                return newHead;
            }
            else
            {
                prev->next = temp->next;
                free(temp);
                return head;
            }
        }
        prev = temp;
        temp = temp->next;
    }
    return head;
}

SimplexNode *removeSimplex(SimplexNode *head, int u, int v, int w)
{
    SimplexNode *temp = head;
    SimplexNode *prev = NULL;

    while (temp != NULL)
    {
        if ((temp->node1 == u && temp->node2 == v && temp->node3 == w) || (temp->node1 == v && temp->node2 == u && temp->node3 == w) ||
            (temp->node1 == v && temp->node2 == w && temp->node3 == u) || (temp->node1 == u && temp->node2 == w && temp->node3 == v) ||
            (temp->node1 == w && temp->node2 == u && temp->node3 == v) || (temp->node1 == w && temp->node2 == v && temp->node3 == u))
        {
            if (prev == NULL)
            {
                SimplexNode *newHead = temp->next;
                free(temp);
                return newHead;
            }
            else
            {
                prev->next = temp->next;
                free(temp);
                return head;
            }
        }
        prev = temp;
        temp = temp->next;
    }
    return head;
}

int *getAvailableNodes(int *state, EdgeNode **adjacencyList, int i, int N, int requiredCount)
{
    int *availableNodes = (int *)malloc(N * sizeof(int));
    int count = 0;

    for (int n = 0; n < N; n++)
    {
        if (state[n] == 0 && !isConnected(adjacencyList, i, n))
        {
            availableNodes[count++] = n;
            if (count == requiredCount)
            {
                break;
            }
        }
    }

    if (count < requiredCount)
    {
        free(availableNodes);
        return NULL;
    }

    return availableNodes;
}

double highPrecisionRand()
{
    // return (double)rand() / ((double)RAND_MAX + 1);
    return rand() / (double)RAND_MAX;
}

void initializeInfection(int *state, int N, int initialInfected)
{
    for (int i = 0; i < N; i++)
    {
        state[i] = 0; // Initially, all nodes are susceptible
    }
    for (int i = 0; i < initialInfected; i++)
    {
        int randNode;
        do
        {
            randNode = rand() % N;
        } while (state[randNode] == 1); // Ensure no duplicates
        state[randNode] = 1; // Infect the node
    }
}

// Handle edge rewiring based on node states
void handleEdgeRewiring(SimplexList *simplexList, Network *network, int *state, double rewireProb)
{
    // Temporary storage for edges to be rewired
    typedef struct
    {
        int i;
        int neighbor;
    } EdgeToRewire;

    EdgeToRewire edgesToRewire[50000];
    int edgeCount = 0;

    // First loop: collect edges to be rewired
    for (int i = 0; i < MAX_NODES; i++)
    {
        if (state[i] == 0)
        {
            EdgeNode *temp = network->adjacencyList[i];
            while (temp != NULL)
            {
                int neighbor = temp->node;
                // int index_s = isInSimplex(simplexList->simplices, i, neighbor);
                if (state[neighbor] == 1 && highPrecisionRand() < rewireProb)
                {
                    // Store the edge for rewiring
                    edgesToRewire[edgeCount].i = i;
                    edgesToRewire[edgeCount].neighbor = neighbor;
                    edgeCount++;
                }
                temp = temp->next;
            }
        }
    }

    // Second loop: remove and rewire edges
    for (int k = 0; k < edgeCount; k++)
    {
        int i = edgesToRewire[k].i;
        int neighbor = edgesToRewire[k].neighbor;

        // Remove edge between i and neighbor
        network->adjacencyList[i] = removeEdge(network->adjacencyList[i], neighbor);
        network->adjacencyList[neighbor] = removeEdge(network->adjacencyList[neighbor], i);
        network->numEdges--;

        int newNeighbor_1;
        for (size_t y = 0; y < 30000; y++)
        {
            newNeighbor_1 = rand() % MAX_NODES;

            if (newNeighbor_1 != i && state[newNeighbor_1] != 1)
            {

                if (!isConnected(network->adjacencyList, i, newNeighbor_1))
                {
                    addEdge(network, i, newNeighbor_1);
                    addEdge(network, newNeighbor_1, i);
                    network->numEdges++;
                    break;
                }
            }
        }
    }
}

// Handle simplex rewiring based on node states
void handleSimplexRewiring(SimplexList *simplexList, Network *network, int *state)
{

    SimplexNode *current = simplexList->simplices;
    SimplexNode *prev = NULL;
    int simplex_deleted = 0;

    // 遍历现有的单纯形链表
    while (current != NULL)
    {
        int node1 = current->node1;
        int node2 = current->node2;
        int node3 = current->node3;

        // 检查这三个节点是否仍然两两相连
        if (!isConnected(network->adjacencyList, node1, node2) ||
            !isConnected(network->adjacencyList, node1, node3) ||
            !isConnected(network->adjacencyList, node2, node3))
        {
            // 如果不再连通，删除该单纯形
            if (prev == NULL)
            {
                simplexList->simplices = current->next;
            }
            else
            {
                prev->next = current->next;
            }
            SimplexNode *temp = current;
            current = current->next;
            free(temp); // 释放内存
            simplexList->numSimplices--;
            simplex_deleted++;
            continue; // 跳过剩下的逻辑，直接处理下一个单纯形
        }

        prev = current;
        current = current->next;
    }
    // printf("%d %d &&\n", simplex_deleted + simplexList->numSimplices, simplexList->numSimplices);
    int ss = 0;
    // 增加新的单纯形

    int u, v, w;
    // 查找三个两两相连且不在同一个单纯形中的节点
    while (1)
    {
        u = rand() % MAX_NODES;
        v = rand() % MAX_NODES;
        w = rand() % MAX_NODES;
        if (u != v && u != w && v != w)
        {
            if (isConnected(network->adjacencyList, u, v) &&
                isConnected(network->adjacencyList, u, w) &&
                isConnected(network->adjacencyList, v, w) &&
                !isInSimplex(simplexList->simplices, u, v, w))
            {
                if (simplex_deleted != 0)
                {
                    addSimplex(simplexList, u, v, w);
                    ss++;
                }

                if (ss >= simplex_deleted && simplexList->numSimplices == 1461)
                {
                    break;
                }
            }
        }
    }
}

// SIS模型
double runSISModel(Network *network, SimplexList *simplexList, int N, double beta, double gamma, double beta2, double w1, double w2, int steps, int initialInfected, SimulationResult *result)
{
    int *state = (int *)malloc(N * sizeof(int));
    int *new_state = (int *)malloc(N * sizeof(int));

    result->infectionDensity = (double *)calloc(steps, sizeof(double));
    result->averageDegree = (double *)calloc(steps, sizeof(double));
    result->averageSimplexCount = (double *)calloc(steps, sizeof(double));

    for (int step = 0; step < steps; step++)
    {
        result->infectionDensity[step] = 0;
        result->averageDegree[step] = 0;
        result->averageSimplexCount[step] = 0;
    }
    // printEdges(network);
    //  初始化感染
    initializeInfection(state, N, initialInfected);

    for (int i = 0; i < N; i++)
    {
        new_state[i] = state[i];
    }

    // 感染记录表格
    InfectionPairRecord *infectedPairs = (InfectionPairRecord *)malloc(N * N * sizeof(InfectionPairRecord));
    initializeInfectionPairRecord(infectedPairs, N);

    for (int step = 0; step < steps; step++)
    {

        int infectedCount = 0;
        int totalDegree = 0;
        int totalSimplices = 0;

        int q1[MAX_NODES];
        int q2[MAX_NODES];
        for (size_t i = 0; i < MAX_NODES; i++)
        {
            q1[i] = 0;
            q2[i] = 0;
        }

        // 更新感染状态
        for (int i = 0; i < N; i++)
        {
            if (state[i] == 1)
            { // 当前感染
                // 恢复
                if (highPrecisionRand() < gamma)
                {
                    new_state[i] = 0; // 恢复
                }
            }
            else
            { // 当前易感

                // 简单边传播
                EdgeNode *temp = network->adjacencyList[i];

                while (temp != NULL)
                {
                    if (state[temp->node] == 1)
                    {

                        q1[i]++;
                    }
                    temp = temp->next;
                }
            }
        }

        // Check for simplex-based infection
        SimplexNode *simplex = simplexList->simplices;
        while (simplex != NULL)
        {
            int node1 = simplex->node1;
            int node2 = simplex->node2;
            int node3 = simplex->node3;

            if (state[node1] == 0)
            { // Node 1 is susceptible
                if (state[node2] == 1 && state[node3] == 1)
                {
                    q2[node1]++;
                }
            }

            if (state[node2] == 0)
            { // Node 1 is susceptible
                if (state[node1] == 1 && state[node3] == 1)
                {
                    q2[node2]++;
                }
            }

            if (state[node3] == 0)
            { // Node 1 is susceptible
                if (state[node2] == 1 && state[node1] == 1)
                {
                    q2[node3]++;
                }
            }

            simplex = simplex->next;
        }
        for (size_t i = 0; i < MAX_NODES; i++)
        {
            //double B = (q1[i] * beta + q2[i] * beta2);
            double B = 1- pow(1-beta,q1[i])*pow(1-beta2,q2[i]);
            if (highPrecisionRand() < B && state[i] == 0)
            {
                new_state[i] = 1;
            }
        }

        handleEdgeRewiring(simplexList, network, state, w1);
        handleSimplexRewiring(simplexList, network, state);

        for (int i = 0; i < N; i++)
        {
            if (new_state[i] == 1)
            {
                infectedCount++;
            }
            EdgeNode *temp = network->adjacencyList[i];
            int degree = 0;
            while (temp != NULL)
            {
                degree++;
                temp = temp->next;
            }
            totalDegree += degree;

            SimplexNode *simplex = simplexList->simplices;
            int simplexCount = 0;
            while (simplex != NULL)
            {
                if (simplex->node1 == i || simplex->node2 == i || simplex->node3 == i)
                {
                    simplexCount++;
                }
                simplex = simplex->next;
            }
            totalSimplices += simplexCount;
            state[i] = new_state[i];
        }

        result->infectionDensity[step] = (double)infectedCount / N;
        result->averageDegree[step] = (double)totalDegree / N;
        result->averageSimplexCount[step] = (double)totalSimplices / N;
        printf("%lf %lf ||", result->averageDegree[step], result->averageSimplexCount[step]);
        printf("%d %lf\n", step, result->infectionDensity[step]);

        // 重置感染记录表格
        initializeInfectionPairRecord(infectedPairs, N);
        if (infectedCount == 0)
        {
            break;
        }
    }

    double sumi=0;
    for (int i = 0; i < steps; i++)
    {
        if (i>=4000)
        {
            sumi+=result->infectionDensity[i];
        }
        
    }

    FILE *fp;
    fp = fopen(Data_w, "w");
    for (size_t i = 0; i < 5000; i++)
    {

        if (fp)
        {
            fprintf(fp, "%lf\n", result->infectionDensity[i]);
        }
    }
    fclose(fp);

    free(state);
    free(new_state);
    free(infectedPairs);

    return sumi/1000;
}

// 平均化多个实验结果
void averageSimulationResults(SimulationResult *results, int numSimulations, int steps, SimulationResult *averagedResult)
{
    averagedResult->infectionDensity = (double *)calloc(steps, sizeof(double));
    averagedResult->averageDegree = (double *)calloc(steps, sizeof(double));
    averagedResult->averageSimplexCount = (double *)calloc(steps, sizeof(double));

    for (int i = 0; i < steps; i++)
    {
        for (int j = 0; j < numSimulations; j++)
        {
            averagedResult->infectionDensity[i] += results[j].infectionDensity[i];
            averagedResult->averageDegree[i] += results[j].averageDegree[i];
            averagedResult->averageSimplexCount[i] += results[j].averageSimplexCount[i];
        }
        averagedResult->infectionDensity[i] /= numSimulations;
        averagedResult->averageDegree[i] /= numSimulations;
        averagedResult->averageSimplexCount[i] /= numSimulations;
    }
}

double averageSimulationResults_Max(SimulationResult *results, int numSimulations, int steps, SimulationResult *averagedResult)
{
    double max = 0;
    averagedResult->infectionDensity = (double *)calloc(steps, sizeof(double));
    averagedResult->averageDegree = (double *)calloc(steps, sizeof(double));
    averagedResult->averageSimplexCount = (double *)calloc(steps, sizeof(double));

    for (int i = 4450; i < steps; i++)
    {
        for (int j = 0; j < numSimulations; j++)
        {
            if (max < results[j].infectionDensity[i])
            {
                max = results[j].infectionDensity[i];
            }
        }
    }

    return max;
}

double averageSimulationResults_Min(SimulationResult *results, int numSimulations, int steps, SimulationResult *averagedResult)
{
    double min = 1000;
    averagedResult->infectionDensity = (double *)calloc(steps, sizeof(double));
    averagedResult->averageDegree = (double *)calloc(steps, sizeof(double));
    averagedResult->averageSimplexCount = (double *)calloc(steps, sizeof(double));

    for (int i = 4450; i < steps; i++)
    {
        for (int j = 0; j < numSimulations; j++)
        {
            if (min > results[j].infectionDensity[i])
            {
                min = results[j].infectionDensity[i];
            }
        }
    }
    return min;
}

int main()
{
    const int N = MAX_NODES; // 假设最大节点数为 2000
    const int numSimulations = 1;
    const int steps = 5000;
    // const double beta = 0.015;
    const double gamma = 0.1;
    const double beta2 = 2.0 * gamma / 4.38;
    // const double beta2 = 1.0;
    const double w1 = 0.8;
    const double w2 = 0.0;
    const int initialInfected = 600;

    srand(time(NULL)); // 设置随机种子

    double adj_res[31];
    for (int z = 4; z < 31; z++)
    {
        SimulationResult *results = (SimulationResult *)malloc(numSimulations * sizeof(SimulationResult));
        double beta = z * 1.0 * gamma / 20.688;
        double max=0,min=1000;
        for (int i = 0; i < numSimulations; i++)
        {
            Network network = {(EdgeNode **)calloc(N, sizeof(EdgeNode *)), 0};
            SimplexList simplexList = {NULL, 0};
            double data;
            // 读取网络和单纯形数据
            readNetwork("edges2.txt", "simplices2.txt", &network, &simplexList);
            data=runSISModel(&network, &simplexList, N, beta, gamma, beta2, w1, w2, steps, initialInfected, &results[i]);
            if (data>max)
            {
                max=data;
            }

            if (data<min)
            {
                min=data;
            }
            
            //printf("%lf ",data);
            // 释放网络和单纯形数据
            for (int j = 0; j < N; j++)
            {
                EdgeNode *temp = network.adjacencyList[j];
                while (temp != NULL)
                {
                    EdgeNode *toFree = temp;
                    temp = temp->next;
                    free(toFree);
                }
            }
            free(network.adjacencyList);

            SimplexNode *tempSimplex = simplexList.simplices;
            while (tempSimplex != NULL)
            {
                SimplexNode *toFree = tempSimplex;
                tempSimplex = tempSimplex->next;
                free(toFree);
            }
        }

        SimulationResult averagedResult;
        averageSimulationResults(results, numSimulations, steps, &averagedResult);
        // 打印结果
        // for (int i = 0; i < steps; i++)
        // {
        //     printf("Step %d: Infection Density = %lf, Average Degree = %.2f, Average Simplex Count = %.2f\n",
        //            i, averagedResult.infectionDensity[i], averagedResult.averageDegree[i], averagedResult.averageSimplexCount[i]);
        //     // printf("%lf\n",averagedResult.infectionDensity[i]);
        // }
        adj_res[z] = 0;

        for (size_t c = 0; c < 5000; c++)
        {
            if (c >= 4000)
            {
                adj_res[z] += averagedResult.infectionDensity[c];
            }
        }
        //double min = averageSimulationResults_Min(results, numSimulations, steps, &averagedResult), max = averageSimulationResults_Max(results, numSimulations, steps, &averagedResult);
        printf("%lf %lf %lf %lf\n", z * 1.0, adj_res[z] / 1000,max,min);

        // 释放内存

        free(results);
        free(averagedResult.infectionDensity);
        free(averagedResult.averageDegree);
        free(averagedResult.averageSimplexCount);
        break;
    }

    system("pause");
    return 0;
}
