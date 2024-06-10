#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>
#include <queue>
#include <functional>

typedef std::vector<int> matrix_row;
typedef std::vector<matrix_row> matrix;

matrix read_matrix_from_file(const std::string &filename, int invalid_value = 999)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Erro ao abrir o arquivo: " << filename << std::endl;
        exit(1);
    }

    int num_nodes;
    file >> num_nodes;

    matrix adj(num_nodes, matrix_row(num_nodes, 0));

    for (int i = 0; i < num_nodes; ++i)
    {
        for (int j = 0; j < num_nodes; ++j)
        {
            file >> adj[i][j];
            if (adj[i][j] == invalid_value)
            {
                adj[i][j] = 0;
            }
        }
    }

    file.close();

    return adj;
}





matrix adjacency_to_incidence(const matrix &adj)
{
    int cols = adj.size();
    int rows = adj[0].size();
    int edge = 0;
    matrix incidence;

    for (int col = 0; col < cols; ++col)
    {
        for (int row = 0; row <= col; ++row)
        {
            if (adj[col][row] != 0)
            {
                incidence.push_back(matrix_row(cols, 0));
                incidence[edge][row] = incidence[edge][col] = 1;
                ++edge;
            }
        }
    }

    return incidence;
}

matrix incidence_to_adjacency(const matrix &inc)
{
    int edges = inc.size();
    assert(edges > 0);

    int vertices = inc[0].size();
    assert(vertices > 0);

    matrix adjacency(vertices, matrix_row(vertices, 0));

    for (int edge = 0; edge < edges; ++edge)
    {
        int a = -1, b = -1, vertex = 0;
        for (; vertex < vertices && a == -1; ++vertex)
        {
            if (inc[edge][vertex]) a = vertex;
        }
        for (; vertex < vertices && b == -1; ++vertex)
        {
            if (inc[edge][vertex]) b = vertex;
        }
        if (b == -1) b = a;
        adjacency[a][b] = adjacency[b][a] = 1;
    }

    return adjacency;
}

matrix adjacency_to_adjacencyList(const matrix &adj)
{
    matrix adjList(adj.size());

    for (int i = 0; i < adj.size(); i++)
    {
        adjList[i] = matrix_row();

        for (int j = 0; j < adj[i].size(); j++)
        {
            if (adj[i][j] > 0)
            {
                adjList[i].push_back(j);
            }
            else
            {
                adjList[i].push_back(-1);
            }
        }
    }

    return adjList;
}

void print_matrix(const matrix &m)
{
    int cols = m.size();
    if (cols == 0)
        return;
    int rows = m[0].size();
    if (rows == 0)
        return;

    for (int c = 0; c < cols; ++c)
    {
        for (int r = 0; r < rows; ++r)
        {
            std::cout << m[c][r] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_matrix_to_file(const matrix &m, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Erro ao abrir o arquivo: " << filename << std::endl;
        exit(1);
    }

    int cols = m.size();
    if (cols == 0)
        return;
    int rows = m[0].size();
    if (rows == 0)
        return;

    for (int c = 0; c < cols; ++c)
    {
        for (int r = 0; r < rows; ++r)
        {
            file << m[c][r] << " ";
        }
        file << std::endl;
    }

    file.close();
}

void dfs(const matrix &graph, std::vector<bool> &visited, int at)
{
    if (visited[at])
        return;

    visited[at] = true;
    std::cout << "No visitado: " << at << std::endl;

    const matrix_row &neighbors = graph[at];
    for (int next : neighbors)
    {
        if (next != -1)
        {
            dfs(graph, visited, next);
        }
    }
}

void bfs(const matrix &graph, std::vector<bool> &visited, int start_node)
{
    std::queue<int> q;

    q.push(start_node);
    visited[start_node] = true;

    while (!q.empty())
    {
        int current_node = q.front();
        q.pop();

        std::cout << "No visitado: " << current_node << std::endl;

        const matrix_row &neighbors = graph[current_node];
        for (int neighbor : neighbors)
        {
            if (neighbor != -1 && !visited[neighbor])
            {
                q.push(neighbor);
                visited[neighbor] = true;
            }
        }
    }
}

struct Edge
{
    int from, to, cost;

    Edge(int f, int t, int c) : from(f), to(t), cost(c) {}

    bool operator<(const Edge &other) const
    {
        return cost > other.cost;
    }
};

void addEdges(int nodeIndex, const matrix &adj, std::priority_queue<Edge> &pq, std::vector<bool> &visited)
{
    visited[nodeIndex] = true;

    for (int neighbor = 0; neighbor < adj[nodeIndex].size(); ++neighbor)
    {
        if (adj[nodeIndex][neighbor] > 0 && !visited[neighbor])
        {
            pq.push(Edge(nodeIndex, neighbor, adj[nodeIndex][neighbor]));
        }
    }
}

std::pair<int, std::vector<Edge>> lazyPrims(const matrix &adj, int start_node = 0)
{
    int n = adj.size();
    int m = n - 1;

    int edgeCount = 0, mstCost = 0;
    std::vector<Edge> mstEdges(m, Edge(-1, -1, -1));

    std::priority_queue<Edge> pq;
    std::vector<bool> visited(n, false);

    addEdges(start_node, adj, pq, visited);

    while (!pq.empty() && edgeCount != m)
    {
        Edge edge = pq.top();
        pq.pop();

        int nodeIndex = edge.to;

        if (visited[nodeIndex])
        {
            continue;
        }

        mstEdges[edgeCount++] = edge;
        mstCost += edge.cost;

        addEdges(nodeIndex, adj, pq, visited);
    }

    if (edgeCount != m)
    {
        return std::make_pair(-1, std::vector<Edge>());
    }

    return std::make_pair(mstCost, mstEdges);
}

struct WeightedEdge
{
    int to, cost;

    WeightedEdge(int t, int c) : to(t), cost(c) {}
};

std::vector<int> dijkstra(const matrix &graph, int start_node)
{
    int n = graph.size();
    std::vector<int> dist(n, 999);
    std::vector<bool> visited(n, false);

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    dist[start_node] = 0;
    pq.push({start_node, 0});

    while (!pq.empty())
    {
        int current_node = pq.top().first;
        int current_dist = pq.top().second;
        pq.pop();

        if (visited[current_node])
            continue;

        visited[current_node] = true;

        const matrix_row &neighbors = graph[current_node];
        for (int i = 0; i < neighbors.size(); ++i)
        {
            int neighbor = neighbors[i];
            if (neighbor != -1)
            {
                int new_dist = current_dist + graph[current_node][i];

                if (new_dist < dist[neighbor])
                {
                    dist[neighbor] = new_dist;
                    pq.push({neighbor, new_dist});
                }
            }
        }
    }

    return dist;
}

int main()
{
    std::string input_filename = "C:/Users/Renan/Documents/Trabalho Teoria dos Grafos/Matriz Adjacencia.txt";
    matrix adj = read_matrix_from_file(input_filename);

    std::cout << "Matriz de Adjacencia" << std::endl;
    print_matrix(adj);
    print_matrix_to_file(adj, "C:/Users/Renan/Documents/Trabalho Teoria dos Grafos/output_adjacency_matrix.txt");

    std::cout << "Matriz de Incidencia" << std::endl;
    matrix incidence = adjacency_to_incidence(adj);
    print_matrix(incidence);
    print_matrix_to_file(incidence, "C:/Users/Renan/Documents/Trabalho Teoria dos Grafos/output_incidence_matrix.txt");

    matrix adjacency = incidence_to_adjacency(incidence);
    print_matrix(adjacency);
    print_matrix_to_file(adjacency, "C:/Users/Renan/Documents/Trabalho Teoria dos Grafos/output_adjacency_matrix_from_incidence.txt");

    std::cout << "Lista de Adjacencia" << std::endl;
    matrix adjList = adjacency_to_adjacencyList(adj);
    print_matrix(adjList);

    std::vector<bool> visited(adjList.size(), false);
    int start_node = 0;
    std::cout << "Busca em Profundidade a partir do no " << start_node << ":" << std::endl;
    dfs(adjList, visited, start_node);

    std::cout << std::endl;

    std::vector<bool> bfs_visited(adjList.size(), false);
    int start_node_bfs = 0;
    std::cout << "Busca em Largura a partir do no " << start_node_bfs << ":" << std::endl;
    bfs(adjList, bfs_visited, start_node_bfs);

    std::cout << std::endl;

    int start_node_prim = 0;
    auto primResult = lazyPrims(adj, start_node_prim);

    if (primResult.first != -1)
    {
        std::cout << "Custo da Arvore Geradora Minima de Prim: " << primResult.first << std::endl;
        std::cout << "Arestas na AGM:" << std::endl;
        for (const auto &edge : primResult.second)
        {
            std::cout << "Aresta: " << edge.from << " - " << edge.to << ", Custo: " << edge.cost << std::endl;
        }
    }
    else
    {
        std::cout << "Nao existe Arvore Geradora Minima!" << std::endl;
    }

    int start_node_dijkstra = 0;
    std::vector<int> dijkstra_result = dijkstra(adj, start_node_dijkstra);

    std::cout << "Distancias Minimas de Dijkstra a partir do no " << start_node_dijkstra << ":\n";
    for (int i = 0; i < dijkstra_result.size(); ++i)
    {
        std::cout << "No " << i << ": " << dijkstra_result[i] << std::endl;
    }

    return 0;
}
