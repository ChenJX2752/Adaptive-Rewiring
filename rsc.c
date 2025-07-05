#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//结构体，存储1-单纯形（边）
typedef struct EdgeNode {
    int node1;
    int node2;
    struct EdgeNode* next;
} EdgeNode;

typedef struct {
    EdgeNode* head;
} EdgeList;

//添加和检查边
void addEdge(EdgeList* edgeList, int node1, int node2) {
    EdgeNode* newNode = (EdgeNode*)malloc(sizeof(EdgeNode));
    newNode->node1 = node1;
    newNode->node2 = node2;
    newNode->next = edgeList->head;
    edgeList->head = newNode;
}

int edgeExists(EdgeList* edgeList, int node1, int node2) {
    EdgeNode* current = edgeList->head;
    while (current != NULL) {
        if ((current->node1 == node1 && current->node2 == node2) ||
            (current->node1 == node2 && current->node2 == node1)) {
            return 1;
        }
        current = current->next;
    }
    return 0;
}

void freeEdgeList(EdgeList* edgeList) {
    EdgeNode* current = edgeList->head;
    EdgeNode* next;
    while (current != NULL) {
        next = current->next;
        free(current);
        current = next;
    }
    edgeList->head = NULL;
}
//随机数，和p2进行判断，p2可能是0.0000000几，一个((double)rand() * RAND_MAX精度不够
double highPrecisionRand() {
    //return ((double) rand() / RAND_MAX + (double) rand() / RAND_MAX ) / 2.0;
    return ((double)rand() * RAND_MAX + rand()) / (RAND_MAX * RAND_MAX);
}
//生成单纯复形，N节点数，K1和K2分别是边和2-单纯形平均度
void generateNetwork(int N, double k1, double k2, const char* edgesFile, const char* simplicesFile) {
    double p1 = (k1 - 2 * k2) / ((N - 1) - 2 * k2);
    double p2 = 2 * k2 / ((N - 1) * (N - 2));//p1和p2的公式见文章附录
    
    FILE *edges_fp = fopen(edgesFile, "w");
    FILE *simplices_fp = fopen(simplicesFile, "w");

    if (edges_fp == NULL || simplices_fp == NULL) {
        printf("Error opening file!\n");
        return;
    }

    // Seed the random number generator
    srand(time(NULL));

    EdgeList edgeList = {NULL};

    // Generate edges，遍历全部两个节点，用随机数判断是否小于p1，是则生成边，并存储和记录
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (((double) rand() / RAND_MAX) < p1) {
                addEdge(&edgeList, i, j);
                fprintf(edges_fp, "%d %d\n", i, j);
            }
        }
    }

    // Generate simplices，遍历全部三个节点，用随机数判断是否小于p2，是则生成2-单纯形，并判断是否需要补充内部的边
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            for (int k = j + 1; k < N; ++k) {
                if (highPrecisionRand() < p2) {
                    fprintf(simplices_fp, "%d %d %d\n", i, j, k);

                    // Add edges if they do not already exist
                    if (!edgeExists(&edgeList, i, j)) {
                        addEdge(&edgeList, i, j);
                        fprintf(edges_fp, "%d %d\n", i, j);
                    }
                    if (!edgeExists(&edgeList, i, k)) {
                        addEdge(&edgeList, i, k);
                        fprintf(edges_fp, "%d %d\n", i, k);
                    }
                    if (!edgeExists(&edgeList, j, k)) {
                        addEdge(&edgeList, j, k);
                        fprintf(edges_fp, "%d %d\n", j, k);
                    }
                }
            }
        }
    }

    fclose(edges_fp);
    fclose(simplices_fp);

    freeEdgeList(&edgeList);
}

int main() {
    int N = 1000;  // 节点数
    double k1 = 20.0;  // 平均度 k1
    double k2 = 4.0;  // 平均度 k2

    generateNetwork(N, k1, k2, "edges2.txt", "simplices2.txt");

    return 0;
}





